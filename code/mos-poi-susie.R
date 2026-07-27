## Mixture-of-submanifolds Poisson factor model.
##
## This is an initial, deliberately concrete version of the MoS layer.  It
## assumes each submanifold has D hard factor slots, and fits observation-level
## mixture responsibilities conditional on a shared nonnegative factor matrix F.
## The hard slots are a useful first bridge from MF-Poisson-SuSiE to the full
## model where each submanifold slot has its own posterior over factors.

if (!exists("normalize_rows")) {
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

loading_scores_from_mf <- function(fit) {
  E_beta = fit$alpha / fit$lambda
  apply(fit$gamma_bar * E_beta, c(1, 3), sum)
}

posterior_inclusion_scores_from_mf <- function(fit) {
  1 - apply(1 - fit$gamma_bar, c(1, 3), prod)
}

gamma_shape_from_moments <- function(E_beta, E_log_beta, n_iter = 8,
                                     min_shape = 1e-3, max_shape = 1e4,
                                     eps = 1e-12) {
  delta = log(pmax(E_beta, eps)) - E_log_beta
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

cluster_accuracy <- function(z_hat, z_true) {
  z_hat = as.integer(as.factor(z_hat))
  z_true = as.integer(as.factor(z_true))
  S_hat = max(z_hat)
  S_true = max(z_true)
  tab = table(z_hat, z_true)

  if (!requireNamespace("clue", quietly = TRUE)) {
    ## Greedy fallback for small diagnostics when clue is not installed.
    used_hat = rep(FALSE, S_hat)
    used_true = rep(FALSE, S_true)
    correct = 0
    for (j in seq_len(min(S_hat, S_true))) {
      score = tab
      score[used_hat, ] = -1
      score[, used_true] = -1
      idx = arrayInd(which.max(score), dim(score))
      used_hat[idx[1]] = TRUE
      used_true[idx[2]] = TRUE
      correct = correct + tab[idx[1], idx[2]]
    }
  } else {
    assignment = clue::solve_LSAP(max(tab) - tab)
    correct = sum(tab[cbind(seq_len(length(assignment)), assignment)])
  }

  as.numeric(correct) / length(z_true)
}

init_motifs_from_loading_scores <- function(scores, S, D, init_seed = NULL,
                                            min_share = 0.10,
                                            allow_repeats = TRUE,
                                            eps = 1e-12) {
  if (!is.null(init_seed)) set.seed(init_seed)

  K = ncol(scores)
  scores_norm = scores / pmax(rowSums(scores), eps)
  km = kmeans(scores_norm, centers = S, nstart = 20)

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

simulate_mos_data <- function(N = 240, M = 300, K = 6,
                              motifs = rbind(c(1, 2), c(1, 3),
                                             c(4, 5), c(4, 6)),
                              shape = 15, rate = 1, factor_shape = 0.15,
                              factor_rate = 0.03, init_seed = NULL) {
  if (!is.null(init_seed)) set.seed(init_seed)

  S = nrow(motifs)
  D = ncol(motifs)
  z = rep(seq_len(S), length.out = N)
  z = sample(z)

  F = matrix(rgamma(K * M, shape = factor_shape, rate = factor_rate),
             nrow = K, ncol = M)
  F = normalize_rows(F)

  L = matrix(0, nrow = N, ncol = K)
  beta = matrix(0, nrow = N, ncol = D)
  for (i in seq_len(N)) {
    support = motifs[z[i], ]
    beta[i, ] = rgamma(D, shape = shape, rate = rate)
    L[i, support] = L[i, support] + beta[i, ]
  }

  Y = matrix(rpois(N * M, L %*% F), nrow = N, ncol = M)

  list(Y = Y, F = F, L = L, z = z, motifs = motifs, beta = beta)
}

fit_hard_motif_components <- function(Y, F, motifs, alpha0, lambda0,
                                      n_inner = 8, eps = 1e-12) {
  N = nrow(Y)
  S = nrow(motifs)
  D = ncol(motifs)
  M = ncol(Y)
  F_sum = rowSums(F)

  components = vector("list", S)
  component_elbo = matrix(0, nrow = N, ncol = S)

  for (s in seq_len(S)) {
    F_s = F[motifs[s, ], , drop = FALSE]
    log_F_s = log(pmax(F_s, eps))
    F_sum_s = F_sum[motifs[s, ]]

    alpha = matrix(alpha0[s, ], nrow = N, ncol = D, byrow = TRUE)
    lambda = matrix(lambda0[s, ] + F_sum_s, nrow = N, ncol = D,
                    byrow = TRUE)
    xi = array(1 / D, dim = c(N, D, M))

    for (inner in seq_len(n_inner)) {
      E_log_beta = digamma(alpha) - log(lambda)
      log_xi = array(0, dim = c(N, D, M))
      for (d in seq_len(D)) {
        log_xi[, d, ] =
          E_log_beta[, d] + matrix(log_F_s[d, ], nrow = N, ncol = M,
                                   byrow = TRUE)
      }

      log_xi_max = log_xi[, 1, ]
      if (D > 1) {
        for (d in 2:D) log_xi_max = pmax(log_xi_max, log_xi[, d, ])
      }
      for (d in seq_len(D)) xi[, d, ] = exp(log_xi[, d, ] - log_xi_max)
      xi_sum = xi[, 1, ]
      if (D > 1) {
        for (d in 2:D) xi_sum = xi_sum + xi[, d, ]
      }
      for (d in seq_len(D)) xi[, d, ] = xi[, d, ] / pmax(xi_sum, eps)

      for (d in seq_len(D)) {
        alpha[, d] = rowSums(xi[, d, ] * Y) + alpha0[s, d]
        lambda[, d] = lambda0[s, d] + F_sum_s[d]
      }
    }

    E_log_beta = digamma(alpha) - log(lambda)
    E_beta = alpha / lambda

    row_elbo = rep(0, N)
    for (d in seq_len(D)) {
      xi_y = xi[, d, ] * Y
      row_elbo = row_elbo + E_log_beta[, d] * rowSums(xi_y)
      row_elbo = row_elbo + rowSums(xi_y *
                                      matrix(log_F_s[d, ], nrow = N,
                                             ncol = M, byrow = TRUE))
      xi_log_xi = xi[, d, ] * log(pmax(xi[, d, ], eps))
      row_elbo = row_elbo - rowSums(Y * xi_log_xi)
      row_elbo = row_elbo - E_beta[, d] * F_sum_s[d]

      kl = (alpha[, d] - alpha0[s, d]) * digamma(alpha[, d]) +
        (log(lambda[, d]) - log(lambda0[s, d])) * alpha0[s, d] -
        (lgamma(alpha[, d]) - lgamma(alpha0[s, d])) -
        (lambda[, d] - lambda0[s, d]) * alpha[, d] / lambda[, d]
      row_elbo = row_elbo - kl
    }

    components[[s]] = list(alpha = alpha, lambda = lambda, xi = xi)
    component_elbo[, s] = row_elbo
  }

  list(components = components, component_elbo = component_elbo)
}

update_mos_priors <- function(components, omega, alpha0, lambda0,
                              min_shape = 1e-3, max_shape = 1e4,
                              eps = 1e-12) {
  S = length(components)
  D = ncol(alpha0)

  for (s in seq_len(S)) {
    w = omega[, s]
    w_sum = sum(w)
    if (w_sum <= eps) next

    for (d in seq_len(D)) {
      alpha = components[[s]]$alpha[, d]
      lambda = components[[s]]$lambda[, d]
      E_beta = sum(w * alpha / lambda) / w_sum
      E_log_beta = sum(w * (digamma(alpha) - log(lambda))) / w_sum

      shape = gamma_shape_from_moments(
        E_beta,
        E_log_beta,
        min_shape = min_shape,
        max_shape = max_shape,
        eps = eps
      )

      alpha0[s, d] = shape
      lambda0[s, d] = shape / pmax(E_beta, eps)
    }
  }

  list(alpha0 = alpha0, lambda0 = lambda0)
}

update_mos_F <- function(Y, F, motifs, components, omega,
                         F_pseudocount = .Machine$double.eps,
                         eps = 1e-12) {
  K = nrow(F)
  M = ncol(F)
  S = nrow(motifs)
  D = ncol(motifs)
  N = nrow(Y)

  A = matrix(F_pseudocount, nrow = K, ncol = M)
  for (s in seq_len(S)) {
    w = matrix(omega[, s], nrow = N, ncol = M)
    for (d in seq_len(D)) {
      k = motifs[s, d]
      A[k, ] = A[k, ] + colSums(w * components[[s]]$xi[, d, ] * Y)
    }
  }

  normalize_rows(A, eps)
}

update_mos_motifs <- function(Y, F, motifs, components, omega, alpha0,
                              lambda0, min_rel_improve = 0,
                              min_share = 0.10, n_loading_iters = 50,
                              eps = 1e-12) {
  S = nrow(motifs)
  D = ncol(motifs)
  K = nrow(F)
  Y_s = t(omega) %*% Y
  total_s = pmax(rowSums(Y_s), eps)

  H = matrix(total_s / K, nrow = S, ncol = K)
  for (s in seq_len(S)) {
    motif_tab = tabulate(motifs[s, ], nbins = K)
    H[s, ] = H[s, ] + total_s[s] * motif_tab / D
  }

  for (iter in seq_len(n_loading_iters)) {
    rate = H %*% F + eps
    H = H * ((Y_s / rate) %*% t(F)) /
      matrix(rowSums(F), nrow = S, ncol = K, byrow = TRUE)
    H = pmax(H, eps)
  }

  share = H / pmax(rowSums(H), eps)
  new_motifs = motifs

  for (s in seq_len(S)) {
    ord = order(share[s, ], decreasing = TRUE)
    active = ord[share[s, ord] >= min_share]
    if (length(active) == 0) active = ord[1]
    if (length(active) > D) active = active[seq_len(D)]
    new_motifs[s, ] = rep(active, length.out = D)
  }

  list(motifs = new_motifs, score = share)
}

mos_poi_susie_fixed_motifs <- function(Y, F, motifs, max_iters = 50,
                                       n_inner = 8, update_prior = TRUE,
                                       prior_shape = 1, prior_rate = 1,
                                       pi_init = NULL, tol = 1e-5,
                                       min_iters = 5, patience = 3,
                                       update_F = FALSE,
                                       update_motifs = FALSE,
                                       motif_update_every = 1,
                                       motif_min_rel_improve = 0,
                                       motif_min_share = 0.10,
                                       motif_loading_iters = 50,
                                       F_step_init = 1,
                                       F_step_ramp = 1,
                                       F_pseudocount = .Machine$double.eps,
                                       eps = 1e-12) {
  N = nrow(Y)
  S = nrow(motifs)
  D = ncol(motifs)

  F = normalize_rows(F, eps)
  motif_update_every = max(1, as.integer(motif_update_every))
  if (is.null(pi_init)) {
    pi = rep(1 / S, S)
  } else {
    pi = pi_init / sum(pi_init)
  }

  alpha0 = matrix(prior_shape, nrow = S, ncol = D)
  lambda0 = matrix(prior_rate, nrow = S, ncol = D)

  elbo = rep(NA_real_, max_iters)
  small_improve_count = 0
  converged = FALSE
  n_iter = max_iters
  component_fit = NULL
  omega = matrix(1 / S, nrow = N, ncol = S)
  motif_history = vector("list", max_iters)

  for (iter in seq_len(max_iters)) {
    rho_F = min(1, F_step_init + (1 - F_step_init) * (iter - 1) /
                  max(F_step_ramp - 1, 1))

    component_fit = fit_hard_motif_components(
      Y = Y,
      F = F,
      motifs = motifs,
      alpha0 = alpha0,
      lambda0 = lambda0,
      n_inner = n_inner,
      eps = eps
    )

    log_resp = component_fit$component_elbo +
      matrix(log(pmax(pi, eps)), nrow = N, ncol = S, byrow = TRUE)
    omega = softmax_rows(log_resp)
    pi = pmax(colMeans(omega), eps)
    pi = pi / sum(pi)

    if (update_prior) {
      prior_fit = update_mos_priors(component_fit$components, omega,
                                    alpha0, lambda0, eps = eps)
      alpha0 = prior_fit$alpha0
      lambda0 = prior_fit$lambda0
    }

    elbo[iter] = sum(apply(log_resp, 1, log_sum_exp))

    if (update_F) {
      F_cavi = update_mos_F(
        Y = Y,
        F = F,
        motifs = motifs,
        components = component_fit$components,
        omega = omega,
        F_pseudocount = F_pseudocount,
        eps = eps
      )
      F = normalize_rows((1 - rho_F) * F + rho_F * F_cavi, eps)
    }

    if (update_motifs && iter %% motif_update_every == 0) {
      motif_fit = update_mos_motifs(
        Y = Y,
        F = F,
        motifs = motifs,
        components = component_fit$components,
        omega = omega,
        alpha0 = alpha0,
        lambda0 = lambda0,
        min_rel_improve = motif_min_rel_improve,
        min_share = motif_min_share,
        n_loading_iters = motif_loading_iters,
        eps = eps
      )
      motifs = motif_fit$motifs
    }
    motif_history[[iter]] = motifs

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

  component_fit = fit_hard_motif_components(
    Y = Y,
    F = F,
    motifs = motifs,
    alpha0 = alpha0,
    lambda0 = lambda0,
    n_inner = n_inner,
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
    motifs = motifs,
    alpha0 = alpha0,
    lambda0 = lambda0,
    components = component_fit$components,
    component_elbo = component_fit$component_elbo,
    elbo = elbo,
    converged = converged,
    n_iter = n_iter,
    motif_history = motif_history[seq_len(n_iter)],
    F = F
  )
}

mos_poi_susie <- function(Y, K, S, D, max_iters = 50, n_inner = 8,
                          mf_max_iters = 60, mf_nmf_iters = 100,
                          init_seed = NULL, update_prior = TRUE,
                          tol = 1e-5, min_iters = 5, patience = 3,
                          update_F = TRUE, update_motifs = TRUE,
                          motif_update_every = 1,
                          motif_min_rel_improve = 0,
                          motif_min_share = 0.10,
                          motif_loading_iters = 50,
                          F_step_init = 0.2,
                          F_step_ramp = 20,
                          F_pseudocount = .Machine$double.eps) {
  mf_fit = mf_poi_susie(
    Y = Y,
    K = K,
    D = D,
    max_iters = mf_max_iters,
    update_prior = TRUE,
    init_seed = init_seed,
    init_F = "poisson_nmf",
    nmf_iters = mf_nmf_iters,
    init_gamma_from_nmf = TRUE,
    elbo_every = 5,
    tol = 5e-4,
    min_iters = 15,
    patience = 2
  )

  scores = loading_scores_from_mf(mf_fit)
  motif_init = init_motifs_from_loading_scores(scores, S = S, D = D,
                                               init_seed = init_seed)

  mos_fit = mos_poi_susie_fixed_motifs(
    Y = Y,
    F = mf_fit$F,
    motifs = motif_init$motifs,
    max_iters = max_iters,
    n_inner = n_inner,
    update_prior = update_prior,
    tol = tol,
    min_iters = min_iters,
    patience = patience,
    update_F = update_F,
    update_motifs = update_motifs,
    motif_update_every = motif_update_every,
    motif_min_rel_improve = motif_min_rel_improve,
    motif_min_share = motif_min_share,
    motif_loading_iters = motif_loading_iters,
    F_step_init = F_step_init,
    F_step_ramp = F_step_ramp,
    F_pseudocount = F_pseudocount
  )

  mos_fit$mf_fit = mf_fit
  mos_fit$motif_init = motif_init
  mos_fit
}
