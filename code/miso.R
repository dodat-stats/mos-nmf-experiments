## Soft-motif mixture-of-submanifolds NMF.
##
## This is the random-gamma version of the MiSo model.  Each submanifold slot
## has a posterior distribution over factor rows instead of a hard factor index.
## The local xi variables are streamed over feature blocks and are not stored.

if (!exists("poisson_susie")) {
  source("code/poisson-susie.R")
}
if (!exists("poisson_susie_nmf")) {
  source("code/joint-learn-susie-poi-F.R")
}

log_sum_exp <- function(x) {
  xmax = max(x)
  xmax + log(sum(exp(x - xmax)))
}

softmax_rows <- function(log_w) {
  row_max = apply(log_w, 1, max)
  w = exp(log_w - row_max)
  w / rowSums(w)
}

loading_scores_from_mf <- poisson_susie_nmf_loading_scores

gamma_shape_from_moments <- function(E_lambda, E_log_lambda, n_iter = 8,
                                     min_shape = 1e-3, max_shape = 1e4,
                                     eps = 1e-12) {
  delta = log(pmax(E_lambda, eps)) - E_log_lambda
  delta = pmax(delta, eps)
  shape = ifelse(delta > 0.5, 1 / (2 * delta), 1 / delta)
  shape = pmin(pmax(shape, min_shape), max_shape)

  for (j in seq_len(n_iter)) {
    f = log(shape) - digamma(shape) - delta
    fp = 1 / shape - trigamma(shape)
    shape = shape - f / fp
    shape = pmin(pmax(shape, min_shape), max_shape)
  }
  shape
}

init_motifs_from_loading_scores <- function(scores, S, D, init_seed = NULL,
                                            min_share = 0.10,
                                            allow_repeats = TRUE,
                                            eps = 1e-12) {
  if (!is.null(init_seed)) set.seed(init_seed)
  K = ncol(scores)
  scores_norm = scores / pmax(rowSums(scores), eps)
  km = kmeans(scores_norm, centers = S, nstart = 20, iter.max = 100)
  motifs = matrix(NA_integer_, nrow = S, ncol = D)

  for (s in seq_len(S)) {
    center = km$centers[s, ] / pmax(sum(km$centers[s, ]), eps)
    ord = order(center, decreasing = TRUE)
    active = ord[center[ord] >= min_share]
    if (length(active) == 0) active = ord[1]
    if (length(active) > D) active = active[seq_len(D)]
    if (allow_repeats) {
      motifs[s, ] = rep(active, length.out = D)
    } else {
      motifs[s, ] = ord[seq_len(D)]
    }
  }

  list(motifs = motifs, cluster = km$cluster, centers = km$centers)
}

init_soft_gamma_from_loading_scores <- function(scores, S, D, init_seed = NULL,
                                                min_share = 0.10,
                                                gamma_floor = 0.05,
                                                surplus_slots = c("repeat",
                                                                  "uniform"),
                                                eps = 1e-12) {
  surplus_slots = match.arg(surplus_slots)
  hard_init = init_motifs_from_loading_scores(
    scores = scores,
    S = S,
    D = D,
    init_seed = init_seed,
    min_share = min_share,
    allow_repeats = TRUE,
    eps = eps
  )

  K = ncol(scores)
  gamma_bar = array(gamma_floor / K, dim = c(S, D, K))
  for (s in seq_len(S)) {
    seen = rep(FALSE, K)
    for (d in seq_len(D)) {
      k = hard_init$motifs[s, d]
      is_surplus_repeat = seen[k]
      if (!(surplus_slots == "uniform" && is_surplus_repeat)) {
        gamma_bar[s, d, k] = 1 - gamma_floor + gamma_floor / K
      }
      seen[k] = TRUE
    }
  }

  list(gamma_bar = gamma_bar, hard_init = hard_init)
}

gamma_bar_from_motifs <- function(motifs, K, gamma_floor = 0.05) {
  S = nrow(motifs)
  D = ncol(motifs)
  gamma_bar = array(gamma_floor / K, dim = c(S, D, K))
  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      gamma_bar[s, d, motifs[s, d]] = 1 - gamma_floor + gamma_floor / K
    }
  }
  gamma_bar
}

normalize_gamma_bar <- function(gamma_bar, eps = 1e-12) {
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  K = dim(gamma_bar)[3]
  gamma_bar = pmax(gamma_bar, eps)
  denom = rowSums(gamma_bar, dims = 2)
  for (k in seq_len(K)) gamma_bar[, , k] = gamma_bar[, , k] / denom
  gamma_bar
}

miso_beta_sd <- function(F, gamma_bar, beta0, slot_scale = NULL) {
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  F_sum = rowSums(F)
  if (is.null(slot_scale)) slot_scale = matrix(1, S, D)
  beta_sd = beta0
  for (s in seq_len(S)) {
    beta_sd[s, ] = beta0[s, ] +
      slot_scale[s, ] * as.vector(gamma_bar[s, , ] %*% F_sum)
  }
  beta_sd
}

miso_accumulate_xi <- function(Y, F, gamma_bar, alpha, alpha0, beta0,
                                   slot_scale = NULL,
                                   omega = NULL, block_size = 100,
                                   compute_C = FALSE,
                                   compute_component_elbo = FALSE,
                                   eps = 1e-12) {
  N = nrow(Y)
  M = ncol(Y)
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  K = nrow(F)
  log_F = log(pmax(F, eps))
  F_sum = rowSums(F)
  if (is.null(slot_scale)) slot_scale = matrix(1, S, D)
  beta_sd = miso_beta_sd(F, gamma_bar, beta0, slot_scale)

  allocated_count = array(0, dim = c(N, S, D))
  C = if (compute_C) array(0, dim = c(S, D, M)) else NULL
  component_elbo = if (compute_component_elbo) matrix(0, N, S) else NULL
  feature_blocks = split(seq_len(M), ceiling(seq_len(M) / block_size))

  for (s in seq_len(S)) {
    alpha_s = matrix(alpha[, s, ], nrow = N, ncol = D)
    beta_s = beta_sd[s, ]
    log_slot_scale_s = log(pmax(slot_scale[s, ], eps))
    E_log_lambda = digamma(alpha_s) -
      matrix(log(beta_s), nrow = N, ncol = D, byrow = TRUE)
    E_lambda = alpha_s /
      matrix(beta_s, nrow = N, ncol = D, byrow = TRUE)
    gamma_s = matrix(gamma_bar[s, , ], nrow = D, ncol = K)

    if (compute_C) {
      weight_s = matrix(omega[, s], nrow = N, ncol = block_size)
    }

    for (block in feature_blocks) {
      B = length(block)
      Yb = Y[, block, drop = FALSE]
      E_log_F_gamma = gamma_s %*% log_F[, block, drop = FALSE]

      log_xi = array(0, dim = c(N, D, B))
      for (d in seq_len(D)) {
        log_xi[, d, ] =
          log_slot_scale_s[d] +
          E_log_lambda[, d] +
          matrix(E_log_F_gamma[d, ], nrow = N, ncol = B, byrow = TRUE)
      }

      log_xi_max = log_xi[, 1, ]
      if (D > 1) {
        for (d in 2:D) log_xi_max = pmax(log_xi_max, log_xi[, d, ])
      }

      xi = array(0, dim = c(N, D, B))
      for (d in seq_len(D)) xi[, d, ] = exp(log_xi[, d, ] - log_xi_max)
      xi_sum = xi[, 1, ]
      if (D > 1) {
        for (d in 2:D) xi_sum = xi_sum + xi[, d, ]
      }
      for (d in seq_len(D)) xi[, d, ] = xi[, d, ] / pmax(xi_sum, eps)

      if (compute_C) {
        if (ncol(weight_s) != B) {
          weight_block = matrix(omega[, s], nrow = N, ncol = B)
        } else {
          weight_block = weight_s
        }
      }

      for (d in seq_len(D)) {
        xi_y = xi[, d, ] * Yb
        allocated_count[, s, d] =
          allocated_count[, s, d] + rowSums(xi_y)

        if (compute_C) {
          C[s, d, block] = colSums(weight_block * xi_y)
        }

        if (compute_component_elbo) {
          component_elbo[, s] = component_elbo[, s] +
            E_log_lambda[, d] * rowSums(xi_y) +
            rowSums(xi_y * matrix(E_log_F_gamma[d, ],
                                  nrow = N, ncol = B, byrow = TRUE))

          xi_log_xi = xi[, d, ] * log(pmax(xi[, d, ], eps))
          component_elbo[, s] = component_elbo[, s] -
            rowSums(Yb * xi_log_xi)
        }
      }
    }

    if (compute_component_elbo) {
      expected_F_sum = slot_scale[s, ] * as.vector(gamma_s %*% F_sum)
      for (d in seq_len(D)) {
        component_elbo[, s] = component_elbo[, s] -
          E_lambda[, d] * expected_F_sum[d]

        a0 = alpha0[s, d]
        b0 = beta0[s, d]
        kl = (alpha_s[, d] - a0) * digamma(alpha_s[, d]) +
          (log(beta_s[d]) - log(b0)) * a0 -
          (lgamma(alpha_s[, d]) - lgamma(a0)) -
          (beta_s[d] - b0) * alpha_s[, d] / beta_s[d]
        component_elbo[, s] = component_elbo[, s] - kl
      }
    }
  }

  list(
    allocated_count = allocated_count,
    C = C,
    component_elbo = component_elbo,
    beta_sd = beta_sd
  )
}

fit_miso_components <- function(Y, F, gamma_bar, alpha0, beta0,
                                    slot_scale = NULL,
                                    alpha_init = NULL, n_inner = 5,
                                    block_size = 100, omega = NULL,
                                    compute_C = FALSE,
                                    compute_component_elbo = TRUE,
                                    eps = 1e-12) {
  N = nrow(Y)
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  if (is.null(slot_scale)) slot_scale = matrix(1, S, D)

  if (is.null(alpha_init)) {
    alpha = array(0, dim = c(N, S, D))
    for (s in seq_len(S)) {
      for (d in seq_len(D)) alpha[, s, d] = alpha0[s, d] + rowSums(Y) / D
    }
  } else {
    alpha = alpha_init
  }

  for (inner in seq_len(n_inner)) {
    counts = miso_accumulate_xi(
      Y = Y,
      F = F,
      gamma_bar = gamma_bar,
      alpha = alpha,
      alpha0 = alpha0,
      beta0 = beta0,
      slot_scale = slot_scale,
      block_size = block_size,
      compute_C = FALSE,
      compute_component_elbo = FALSE,
      eps = eps
    )$allocated_count

    for (s in seq_len(S)) {
      for (d in seq_len(D)) alpha[, s, d] = alpha0[s, d] + counts[, s, d]
    }
  }

  final_counts = miso_accumulate_xi(
    Y = Y,
    F = F,
    gamma_bar = gamma_bar,
    alpha = alpha,
    alpha0 = alpha0,
    beta0 = beta0,
    slot_scale = slot_scale,
    block_size = block_size,
    compute_C = FALSE,
    compute_component_elbo = FALSE,
    eps = eps
  )$allocated_count
  for (s in seq_len(S)) {
    for (d in seq_len(D)) alpha[, s, d] = alpha0[s, d] + final_counts[, s, d]
  }

  stats = miso_accumulate_xi(
    Y = Y,
    F = F,
    gamma_bar = gamma_bar,
    alpha = alpha,
    alpha0 = alpha0,
    beta0 = beta0,
    slot_scale = slot_scale,
    omega = omega,
    block_size = block_size,
    compute_C = compute_C,
    compute_component_elbo = compute_component_elbo,
    eps = eps
  )

  list(
    alpha = alpha,
    beta_sd = stats$beta_sd,
    allocated_count = stats$allocated_count,
    C = stats$C,
    component_elbo = stats$component_elbo
  )
}

update_miso_gamma <- function(F, gamma_bar, C, alpha, beta_sd, omega,
                                  slot_scale = NULL,
                                  rho_prior = NULL, gamma_step = 1,
                                  eps = 1e-12) {
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  K = nrow(F)
  log_F = log(pmax(F, eps))
  F_sum = rowSums(F)
  if (is.null(slot_scale)) slot_scale = matrix(1, S, D)

  if (is.null(rho_prior)) {
    rho_prior = array(1 / K, dim = c(S, D, K))
  }

  gamma_new = gamma_bar
  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      E_lambda = alpha[, s, d] / beta_sd[s, d]
      penalty = slot_scale[s, d] * sum(omega[, s] * E_lambda) * F_sum
      score = as.vector(log_F %*% C[s, d, ]) - penalty +
        log(pmax(rho_prior[s, d, ], eps))
      score = score - max(score)
      prob = exp(score)
      prob = prob / sum(prob)
      gamma_new[s, d, ] = (1 - gamma_step) * gamma_bar[s, d, ] +
        gamma_step * prob
    }
  }

  normalize_gamma_bar(gamma_new, eps)
}

update_miso_slot_scale <- function(F, gamma_bar, C, alpha, beta_sd,
                                       omega, slot_scale,
                                       slot_scale_step = 1,
                                       min_scale = 1e-3,
                                       max_scale = 1e3,
                                       normalize_within_submanifold = TRUE,
                                       eps = 1e-12) {
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  F_sum = rowSums(F)
  scale_new = slot_scale

  for (s in seq_len(S)) {
    gamma_s = matrix(gamma_bar[s, , ], nrow = D, ncol = nrow(F))
    expected_F_sum = as.vector(gamma_s %*% F_sum)
    for (d in seq_len(D)) {
      count_sd = sum(C[s, d, ])
      E_lambda = alpha[, s, d] / beta_sd[s, d]
      exposure_sd = sum(omega[, s] * E_lambda) * expected_F_sum[d]
      scale_cavi = count_sd / pmax(exposure_sd, eps)
      scale_new[s, d] = (1 - slot_scale_step) * slot_scale[s, d] +
        slot_scale_step * scale_cavi
    }

    if (normalize_within_submanifold) {
      scale_new[s, ] = scale_new[s, ] / pmax(mean(scale_new[s, ]), eps)
    }
  }

  pmin(pmax(scale_new, min_scale), max_scale)
}

update_miso_F <- function(F, gamma_bar, C,
                              F_step = 1,
                              F_pseudocount = .Machine$double.eps,
                              eps = 1e-12) {
  K = nrow(F)
  M = ncol(F)
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  A = matrix(F_pseudocount, K, M)

  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      A = A + gamma_bar[s, d, ] %o% C[s, d, ]
    }
  }

  F_cavi = normalize_rows(A, eps)
  normalize_rows((1 - F_step) * F + F_step * F_cavi, eps)
}

update_miso_priors <- function(alpha, beta_sd, omega, alpha0, beta0,
                                   min_shape = 1e-3, max_shape = 1e4,
                                   eps = 1e-12) {
  S = ncol(omega)
  D = ncol(alpha0)

  for (s in seq_len(S)) {
    w = omega[, s]
    w_sum = sum(w)
    if (w_sum <= eps) next

    for (d in seq_len(D)) {
      E_lambda = sum(w * alpha[, s, d] / beta_sd[s, d]) / w_sum
      E_log_lambda = sum(w * (digamma(alpha[, s, d]) -
                                log(beta_sd[s, d]))) / w_sum

      shape = gamma_shape_from_moments(
        E_lambda = E_lambda,
        E_log_lambda = E_log_lambda,
        min_shape = min_shape,
        max_shape = max_shape,
        eps = eps
      )
      alpha0[s, d] = shape
      beta0[s, d] = shape / pmax(E_lambda, eps)
    }
  }

  list(alpha0 = alpha0, beta0 = beta0)
}

miso_gamma_kl <- function(gamma_bar, rho_prior = NULL, eps = 1e-12) {
  K = dim(gamma_bar)[3]
  if (is.null(rho_prior)) {
    rho_prior = array(1 / K, dim = dim(gamma_bar))
  }

  sum(gamma_bar * (log(pmax(gamma_bar, eps)) -
                     log(pmax(rho_prior, eps))))
}

miso_fixed_gamma <- function(Y, F, gamma_bar, max_iters = 50,
                                n_inner = 5, update_prior = TRUE,
                                prior_shape = 1, prior_beta = 1,
                                pi_init = NULL, rho_prior = NULL,
                                tol = 1e-5, min_iters = 5,
                                patience = 3, update_F = FALSE,
                                update_gamma = TRUE,
                                update_slot_scale = FALSE,
                                slot_scale_init = NULL,
                                slot_scale_step_init = 0.2,
                                slot_scale_step_ramp = 10,
                                min_slot_scale = 1e-3,
                                max_slot_scale = 1e3,
                                gamma_step_init = 0.5,
                                gamma_step_ramp = 10,
                                F_step_init = 0.2,
                                F_step_ramp = 20,
                                F_pseudocount = .Machine$double.eps,
                                block_size = 100, eps = 1e-12) {
  N = nrow(Y)
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  K = nrow(F)

  F = normalize_rows(F, eps)
  gamma_bar = normalize_gamma_bar(gamma_bar, eps)
  if (is.null(rho_prior)) rho_prior = array(1 / K, dim = c(S, D, K))
  if (is.null(slot_scale_init)) {
    slot_scale = matrix(1, S, D)
  } else {
    slot_scale = slot_scale_init
  }

  alpha0 = matrix(prior_shape, nrow = S, ncol = D)
  beta0 = matrix(prior_beta, nrow = S, ncol = D)
  if (is.null(pi_init)) {
    pi = rep(1 / S, S)
  } else {
    pi = pi_init / sum(pi_init)
  }

  omega = matrix(1 / S, nrow = N, ncol = S)
  alpha = NULL
  elbo = rep(NA_real_, max_iters)
  gamma_history = vector("list", max_iters)
  small_improve_count = 0
  converged = FALSE
  n_iter = max_iters
  component_fit = NULL

  for (iter in seq_len(max_iters)) {
    gamma_step = min(1, gamma_step_init + (1 - gamma_step_init) *
                       (iter - 1) / max(gamma_step_ramp - 1, 1))
    F_step = min(1, F_step_init + (1 - F_step_init) *
                   (iter - 1) / max(F_step_ramp - 1, 1))
    slot_scale_step = min(1, slot_scale_step_init +
                            (1 - slot_scale_step_init) * (iter - 1) /
                            max(slot_scale_step_ramp - 1, 1))

    component_fit = fit_miso_components(
      Y = Y,
      F = F,
      gamma_bar = gamma_bar,
      alpha0 = alpha0,
      beta0 = beta0,
      slot_scale = slot_scale,
      alpha_init = alpha,
      n_inner = n_inner,
      block_size = block_size,
      omega = omega,
      compute_C = TRUE,
      compute_component_elbo = TRUE,
      eps = eps
    )
    alpha = component_fit$alpha

    log_resp = component_fit$component_elbo +
      matrix(log(pmax(pi, eps)), nrow = N, ncol = S, byrow = TRUE)
    omega = softmax_rows(log_resp)
    pi = pmax(colMeans(omega), eps)
    pi = pi / sum(pi)

    component_fit = fit_miso_components(
      Y = Y,
      F = F,
      gamma_bar = gamma_bar,
      alpha0 = alpha0,
      beta0 = beta0,
      slot_scale = slot_scale,
      alpha_init = alpha,
      n_inner = 1,
      block_size = block_size,
      omega = omega,
      compute_C = TRUE,
      compute_component_elbo = TRUE,
      eps = eps
    )
    alpha = component_fit$alpha

    if (update_prior) {
      prior_fit = update_miso_priors(
        alpha = alpha,
        beta_sd = component_fit$beta_sd,
        omega = omega,
        alpha0 = alpha0,
        beta0 = beta0,
        eps = eps
      )
      alpha0 = prior_fit$alpha0
      beta0 = prior_fit$beta0
    }

    if (update_gamma) {
      gamma_bar = update_miso_gamma(
        F = F,
        gamma_bar = gamma_bar,
        C = component_fit$C,
        alpha = alpha,
        beta_sd = component_fit$beta_sd,
        omega = omega,
        slot_scale = slot_scale,
        rho_prior = rho_prior,
        gamma_step = gamma_step,
        eps = eps
      )
    }

    if (update_slot_scale) {
      slot_scale = update_miso_slot_scale(
        F = F,
        gamma_bar = gamma_bar,
        C = component_fit$C,
        alpha = alpha,
        beta_sd = component_fit$beta_sd,
        omega = omega,
        slot_scale = slot_scale,
        slot_scale_step = slot_scale_step,
        min_scale = min_slot_scale,
        max_scale = max_slot_scale,
        eps = eps
      )
    }

    if (update_F) {
      F = update_miso_F(
        F = F,
        gamma_bar = gamma_bar,
        C = component_fit$C,
        F_step = F_step,
        F_pseudocount = F_pseudocount,
        eps = eps
      )
    }

    elbo[iter] = sum(apply(log_resp, 1, log_sum_exp)) -
      miso_gamma_kl(gamma_bar, rho_prior = rho_prior, eps = eps)
    gamma_history[[iter]] = gamma_bar

    if (iter > 1) {
      rel_improve = (elbo[iter] - elbo[iter - 1]) / (abs(elbo[iter - 1]) + 1)
      if (iter >= min_iters && rel_improve < tol) {
        small_improve_count = small_improve_count + 1
      } else {
        small_improve_count = 0
      }
      if (iter >= min_iters && small_improve_count >= patience) {
        converged = TRUE
        n_iter = iter
        break
      }
    }
  }

  component_fit = fit_miso_components(
    Y = Y,
    F = F,
    gamma_bar = gamma_bar,
    alpha0 = alpha0,
    beta0 = beta0,
    slot_scale = slot_scale,
    alpha_init = alpha,
    n_inner = n_inner,
    block_size = block_size,
    omega = omega,
    compute_C = TRUE,
    compute_component_elbo = TRUE,
    eps = eps
  )
  log_resp = component_fit$component_elbo +
    matrix(log(pmax(pi, eps)), nrow = N, ncol = S, byrow = TRUE)
  omega = softmax_rows(log_resp)
  pi = pmax(colMeans(omega), eps)
  pi = pi / sum(pi)

  list(
    omega = omega,
    z_hat = max.col(omega),
    pi = pi,
    gamma_bar = gamma_bar,
    alpha = component_fit$alpha,
    beta_sd = component_fit$beta_sd,
    slot_scale = slot_scale,
    alpha0 = alpha0,
    beta0 = beta0,
    component_elbo = component_fit$component_elbo,
    C = component_fit$C,
    elbo = elbo,
    converged = converged,
    n_iter = n_iter,
    gamma_history = gamma_history[seq_len(n_iter)],
    F = F
  )
}

miso <- function(Y, K, S, D, max_iters = 50, n_inner = 5,
                    mf_max_iters = 60, mf_nmf_iters = 100,
                    init_seed = NULL, update_prior = TRUE,
                    prior_shape = 1, prior_beta = 1,
                    tol = 1e-5, min_iters = 5, patience = 3,
                    update_F = TRUE, update_gamma = TRUE,
                    update_slot_scale = FALSE,
                    slot_scale_step_init = 0.2,
                    slot_scale_step_ramp = 10,
                    min_slot_scale = 1e-3,
                    max_slot_scale = 1e3,
                    gamma_init_floor = 0.05,
                    motif_min_share = 0.10,
                    surplus_slots = c("repeat", "uniform"),
                    gamma_step_init = 0.5,
                    gamma_step_ramp = 10,
                    F_step_init = 0.2,
                    F_step_ramp = 20,
                    F_pseudocount = .Machine$double.eps,
                    block_size = 100) {
  surplus_slots = match.arg(surplus_slots)
  mf_fit = poisson_susie_nmf(
    Y = Y,
    K = K,
    D = D,
    max_iters = mf_max_iters,
    update_prior = TRUE,
    prior_shape = prior_shape,
    prior_beta = prior_beta,
    init_seed = init_seed,
    init_F = "poisson_nmf",
    nmf_iters = mf_nmf_iters,
    init_gamma_from_nmf = TRUE,
    elbo_every = 5,
    tol = 5e-4,
    min_iters = 15,
    patience = 2
  )

  gamma_init = init_soft_gamma_from_loading_scores(
    scores = loading_scores_from_mf(mf_fit),
    S = S,
    D = D,
    init_seed = init_seed,
    min_share = motif_min_share,
    gamma_floor = gamma_init_floor,
    surplus_slots = surplus_slots
  )

  fit = miso_fixed_gamma(
    Y = Y,
    F = mf_fit$F,
    gamma_bar = gamma_init$gamma_bar,
    max_iters = max_iters,
    n_inner = n_inner,
    update_prior = update_prior,
    prior_shape = prior_shape,
    prior_beta = prior_beta,
    tol = tol,
    min_iters = min_iters,
    patience = patience,
    update_F = update_F,
    update_gamma = update_gamma,
    update_slot_scale = update_slot_scale,
    slot_scale_step_init = slot_scale_step_init,
    slot_scale_step_ramp = slot_scale_step_ramp,
    min_slot_scale = min_slot_scale,
    max_slot_scale = max_slot_scale,
    gamma_step_init = gamma_step_init,
    gamma_step_ramp = gamma_step_ramp,
    F_step_init = F_step_init,
    F_step_ramp = F_step_ramp,
    F_pseudocount = F_pseudocount,
    block_size = block_size
  )

  fit$mf_fit = mf_fit
  fit$gamma_init = gamma_init
  fit
}

miso_fit_beta_sd <- function(fit) {
  if (!is.null(fit$beta_sd)) return(fit$beta_sd)
  if (!is.null(fit$lambda_sd)) return(fit$lambda_sd)
  stop("The fit contains neither a beta_sd nor a legacy lambda_sd field")
}

miso_fit_beta0 <- function(fit) {
  if (!is.null(fit$beta0)) return(fit$beta0)
  if (!is.null(fit$lambda0)) return(fit$lambda0)
  stop("The fit contains neither a beta0 nor a legacy lambda0 field")
}

motif_factor_scores <- function(fit) {
  S = dim(fit$gamma_bar)[1]
  D = dim(fit$gamma_bar)[2]
  K = dim(fit$gamma_bar)[3]
  scores = matrix(0, S, K)
  lambda_mean = matrix(0, S, D)
  beta_sd = miso_fit_beta_sd(fit)
  if (is.null(fit$slot_scale)) {
    slot_scale = matrix(1, S, D)
  } else {
    slot_scale = fit$slot_scale
  }

  for (s in seq_len(S)) {
    w = fit$omega[, s]
    w_sum = pmax(sum(w), 1e-12)
    for (d in seq_len(D)) {
      lambda_mean[s, d] =
        sum(w * fit$alpha[, s, d] / beta_sd[s, d]) / w_sum
      scores[s, ] = scores[s, ] +
        slot_scale[s, d] * lambda_mean[s, d] * fit$gamma_bar[s, d, ]
    }
  }

  list(scores = scores, lambda_mean = lambda_mean)
}

motif_gamma_entropy <- function(gamma_bar) {
  entropy = -gamma_bar * log(gamma_bar)
  entropy[is.nan(entropy)] = 0
  rowSums(entropy, dims = 2) / log(dim(gamma_bar)[3])
}

smallest_credible_set <- function(probability, coverage = 0.9,
                                  eps = 1e-12) {
  if (coverage <= 0 || coverage > 1) {
    stop("coverage must be strictly positive and at most one")
  }
  probability = pmax(probability, 0)
  probability = probability / pmax(sum(probability), eps)
  order_probability = order(probability, decreasing = TRUE)
  set_size = which(cumsum(probability[order_probability]) >= coverage)[1]
  sort(order_probability[seq_len(set_size)])
}

credible_set_purity <- function(credible_set, factor_correlation) {
  if (length(credible_set) <= 1) return(1)
  correlation = abs(factor_correlation[credible_set, credible_set,
                                       drop = FALSE])
  pairwise = correlation[lower.tri(correlation)]
  pairwise[is.na(pairwise)] = 0
  min(pairwise)
}

motif_credible_sets <- function(fit, coverage = 0.9,
                                min_abs_correlation = 0.5,
                                deduplicate = TRUE, eps = 1e-12) {
  if (min_abs_correlation < 0 || min_abs_correlation > 1) {
    stop("min_abs_correlation must be between zero and one")
  }

  S = dim(fit$gamma_bar)[1]
  D = dim(fit$gamma_bar)[2]
  factor_correlation = stats::cor(t(fit$F))
  dimensions = do.call(rbind, lapply(seq_len(S), function(s) {
    do.call(rbind, lapply(seq_len(D), function(d) {
      probability = fit$gamma_bar[s, d, ]
      credible_set = smallest_credible_set(probability, coverage, eps)
      data.frame(
        motif = s,
        dimension = d,
        credible_set = paste(credible_set, collapse = ","),
        credible_set_size = length(credible_set),
        posterior_coverage = sum(probability[credible_set]),
        purity = credible_set_purity(credible_set, factor_correlation),
        top_factor = which.max(probability),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(dimensions) = NULL
  dimensions$retained = dimensions$purity >= min_abs_correlation
  dimensions$selected = dimensions$retained

  if (deduplicate) {
    dimensions$selected = FALSE
    for (s in seq_len(S)) {
      retained_index = which(dimensions$motif == s & dimensions$retained)
      if (length(retained_index) == 0) next
      unique_index = retained_index[
        !duplicated(dimensions$credible_set[retained_index])
      ]
      dimensions$selected[unique_index] = TRUE
    }
  }

  selection = do.call(rbind, lapply(seq_len(S), function(s) {
    selected = dimensions$motif == s & dimensions$selected
    data.frame(
      motif = s,
      D_hat = sum(selected),
      selected_credible_sets =
        paste(dimensions$credible_set[selected], collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  rownames(selection) = NULL

  list(
    dimensions = dimensions,
    selection = selection,
    factor_correlation = factor_correlation,
    coverage = coverage,
    min_abs_correlation = min_abs_correlation
  )
}

motif_slot_activity <- function(fit, eps = 1e-12) {
  scores = motif_factor_scores(fit)
  lambda_mean = scores$lambda_mean
  S = nrow(lambda_mean)
  D = ncol(lambda_mean)
  if (is.null(fit$slot_scale)) {
    slot_scale = matrix(1, S, D)
  } else {
    slot_scale = fit$slot_scale
  }
  effective_lambda = lambda_mean * slot_scale
  prior_mean = fit$alpha0 / pmax(miso_fit_beta0(fit), eps)
  gamma_entropy = motif_gamma_entropy(fit$gamma_bar)
  top_gamma = apply(fit$gamma_bar, c(1, 2), max)

  out = do.call(rbind, lapply(seq_len(S), function(s) {
    lambda_total = sum(effective_lambda[s, ])
    prior_total = sum(prior_mean[s, ])
    data.frame(
      submanifold = s,
      slot = seq_len(D),
      lambda_mean = lambda_mean[s, ],
      slot_scale = slot_scale[s, ],
      effective_lambda = effective_lambda[s, ],
      lambda_fraction = effective_lambda[s, ] / pmax(lambda_total, eps),
      prior_mean = prior_mean[s, ],
      prior_fraction = prior_mean[s, ] / pmax(prior_total, eps),
      gamma_entropy = gamma_entropy[s, ],
      top_gamma = top_gamma[s, ]
    )
  }))
  rownames(out) = NULL
  out
}

select_motif_dimensions <- function(fit, prior_fraction_threshold = 0.05,
                                    lambda_fraction_threshold = 0.05,
                                    entropy_threshold = 0.80,
                                    eps = 1e-12) {
  slot_tab = motif_slot_activity(fit, eps = eps)
  slot_tab$active_ard =
    slot_tab$prior_fraction >= prior_fraction_threshold
  slot_tab$active_lambda =
    slot_tab$lambda_fraction >= lambda_fraction_threshold
  slot_tab$inactive_diffuse =
    slot_tab$lambda_fraction < lambda_fraction_threshold &
    slot_tab$gamma_entropy >= entropy_threshold
  slot_tab$active_ard_or_lambda =
    slot_tab$active_ard | slot_tab$active_lambda
  slot_tab$active_diffuse_heuristic = !slot_tab$inactive_diffuse

  selection = do.call(rbind, lapply(split(slot_tab, slot_tab$submanifold),
                                    function(tab) {
    data.frame(
      submanifold = unique(tab$submanifold),
      D_hat_ard = sum(tab$active_ard),
      D_hat_lambda = sum(tab$active_lambda),
      D_hat_ard_or_lambda = sum(tab$active_ard_or_lambda),
      D_hat_diffuse_heuristic = sum(tab$active_diffuse_heuristic),
      max_prior_fraction = max(tab$prior_fraction),
      min_active_prior_fraction =
        if (any(tab$active_ard)) min(tab$prior_fraction[tab$active_ard])
        else NA_real_
    )
  }))
  rownames(selection) = NULL

  list(slots = slot_tab, selection = selection)
}
