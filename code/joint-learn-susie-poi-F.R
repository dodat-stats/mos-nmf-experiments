## ELBO for matrix factorization Poisson-SuSiE
# Y: (N, M)
# F: (K, M)
# xi: (N, D, M)
# gamma_bar: (N, D, K)
# alpha: (N, D, K)
# lambda: (N, D, K)
# alpha0: (N, D)
# lambda0: (N, D)
mf_ELBO <- function(Y, F, xi, gamma_bar, alpha, lambda, alpha0, lambda0) {
  N = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  K = dim(gamma_bar)[3]
  log_F = log(F)
  F_sum = rowSums(F)  # (K,)

  E_log_beta_cond = digamma(alpha) - log(lambda)  # (N, D, K)
  E_log_beta = rowSums(gamma_bar * E_log_beta_cond, dims = 2)  # (N, D)

  term1 = 0
  term2 = 0
  term4 = 0
  for (dd in 1:D) {
    xi_d_Y = xi[, dd, ] * Y  # (N, M)

    # Σ_{i,m} ξ_{imd} y_{im} E[log β_{id}]
    term1 = term1 + sum(E_log_beta[, dd] * rowSums(xi_d_Y))

    # Σ_{i,m} ξ_{imd} y_{im} Σ_k γ̄_{idk} log(F_{km})
    E_log_F_d = gamma_bar[, dd, ] %*% log_F  # (N, M)
    term2 = term2 + sum(xi_d_Y * E_log_F_d)

    # -Σ_{i,m} y_{im} ξ_{imd} log(ξ_{imd})
    xl = xi[, dd, ] * log(xi[, dd, ])
    xl[is.nan(xl)] = 0
    term4 = term4 - sum(xl * Y)
  }

  # -Σ_{i,d,k} γ̄_{idk} (α_{idk}/λ_{idk}) F_sum_k
  gbe = matrix(gamma_bar * alpha / lambda, nrow = N * D, ncol = K)
  term3 = -sum(gbe %*% F_sum)

  # KL discrete: Σ_{i,d,k} γ̄_{idk} (log γ̄_{idk} - log(1/K))
  gb_log_gb = gamma_bar * log(gamma_bar)
  gb_log_gb[is.nan(gb_log_gb)] = 0
  KL_disc = sum(gb_log_gb) + N * D * log(K)

  # KL continuous — loop over k to avoid broadcasting (N,D) -> (N,D,K)
  KL_cont = 0
  for (k in 1:K) {
    KL_cont = KL_cont + sum(gamma_bar[, , k] * (
      (alpha[, , k] - alpha0) * digamma(alpha[, , k])
      + (log(lambda[, , k]) - log(lambda0)) * alpha0
      - (lgamma(alpha[, , k]) - lgamma(alpha0))
      - (lambda[, , k] - lambda0) * alpha[, , k] / lambda[, , k]
    ))
  }

  ELBO = term1 + term2 + term3 + term4 - KL_disc - KL_cont
  return(ELBO)
}

## vectorized inverse digamma via Newton's method (Minka 2000)
inv_digamma <- function(target, n_iter = 5) {
  x = ifelse(target >= -2.22, exp(target) + 0.5, -1 / (target - digamma(1)))
  for (j in 1:n_iter) {
    x = x - (digamma(x) - target) / trigamma(x)
  }
  return(x)
}

normalize_rows <- function(X, eps = .Machine$double.eps) {
  X = pmax(X, eps)
  X / rowSums(X)
}

poisson_nmf_init <- function(Y, K, max_iters = 200, init_seed = NULL,
                             eps = 1e-12) {
  if (!is.null(init_seed)) set.seed(init_seed)

  N = nrow(Y)
  M = ncol(Y)
  total_count = pmax(rowSums(Y), eps)

  F = matrix(rgamma(K * M, shape = 1, rate = 1), nrow = K, ncol = M)
  F = normalize_rows(F, eps)

  L = matrix(rgamma(N * K, shape = 1, rate = 1), nrow = N, ncol = K)
  L = L / rowSums(L)
  L = L * total_count

  for (iter in seq_len(max_iters)) {
    rate = L %*% F + eps
    L = L * ((Y / rate) %*% t(F)) /
      matrix(rowSums(F), nrow = N, ncol = K, byrow = TRUE)
    L = pmax(L, eps)

    rate = L %*% F + eps
    F = F * (t(L) %*% (Y / rate)) /
      matrix(colSums(L), nrow = K, ncol = M)
    F = pmax(F, eps)

    F_scale = rowSums(F)
    F = F / F_scale
    L = sweep(L, 2, F_scale, "*")
  }

  list(L = L, F = F)
}

init_gamma_from_L <- function(L, D, gamma_floor = 1e-4) {
  N = nrow(L)
  K = ncol(L)
  gamma_bar = array(gamma_floor / K, dim = c(N, D, K))

  for (i in seq_len(N)) {
    top_k = order(L[i, ], decreasing = TRUE)
    for (d in seq_len(D)) {
      k = top_k[((d - 1) %% K) + 1]
      gamma_bar[i, d, k] = 1 - gamma_floor + gamma_floor / K
    }
  }

  gamma_bar
}

update_mf_xi <- function(F, gamma_bar, alpha, lambda) {
  N = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  M = ncol(F)
  K = dim(gamma_bar)[3]
  log_F = log(F)

  E_log_beta_cond = digamma(alpha) - log(lambda)  # (N, D, K)
  E_log_beta = rowSums(gamma_bar * E_log_beta_cond, dims = 2)  # (N, D)

  ## R stores arrays with the first index varying fastest, so this reshape maps
  ## row (n, d) to the same order as as.vector(E_log_beta).
  log_xi_mat = matrix(gamma_bar, nrow = N * D, ncol = K) %*% log_F
  log_xi_mat = log_xi_mat + as.vector(E_log_beta)
  log_xi = array(log_xi_mat, dim = c(N, D, M))

  log_xi_max = log_xi[, 1, ]
  if (D > 1) {
    for (d in 2:D) log_xi_max = pmax(log_xi_max, log_xi[, d, ])
  }

  xi = array(0, dim = c(N, D, M))
  for (d in seq_len(D)) xi[, d, ] = exp(log_xi[, d, ] - log_xi_max)
  xi_sum = xi[, 1, ]
  if (D > 1) {
    for (d in 2:D) xi_sum = xi_sum + xi[, d, ]
  }
  for (d in seq_len(D)) xi[, d, ] = xi[, d, ] / xi_sum

  xi
}

## Internal matrix factorization Poisson-SuSiE worker. Prefer the public
## mf_poi_susie(Y, K, D, ...) and mf_poi_susie_fixed_F(Y, F, D, ...) wrappers.
.mf_poi_susie_fit <- function(Y, F, alpha0, lambda0, max_iters = 100,
                              update_prior = TRUE, update_F = TRUE,
                              break_symmetry = FALSE, tol_dist_sim = 1e-4,
                              init_seed = NULL, gamma_step_init = 1,
                              gamma_step_ramp = 50,
                              init_F = c("given", "poisson_nmf"),
                              nmf_iters = 200, init_gamma_from_nmf = FALSE,
                              gamma_init_floor = 1e-4,
                              F_step_init = 1, F_step_ramp = 50,
                              F_pseudocount = .Machine$double.eps,
                              elbo_every = 1, tol = NULL, min_iters = 10,
                              patience = 3) {
  N = nrow(Y)
  M = ncol(Y)
  K = nrow(F)
  D = ncol(alpha0)

  init_F = match.arg(init_F)
  elbo_every = max(1, as.integer(elbo_every))
  min_iters = max(1, as.integer(min_iters))
  patience = max(1, as.integer(patience))

  ## initialization
  if (!is.null(init_seed)) set.seed(init_seed)
  nmf_fit = NULL
  if (init_F == "poisson_nmf") {
    nmf_fit = poisson_nmf_init(Y, K, max_iters = nmf_iters)
    F = nmf_fit$F
  }

  F = normalize_rows(F)
  log_F = log(F)
  F_sum = rowSums(F)  # (K,)

  if (!is.null(nmf_fit) && init_gamma_from_nmf) {
    gamma_bar = init_gamma_from_L(nmf_fit$L, D, gamma_init_floor)
  } else {
    gamma_bar = array(rgamma(N * D * K, 1, 1), dim = c(N, D, K))
    gb_sums = rowSums(gamma_bar, dims = 2)  # (N, D)
    for (k in 1:K) gamma_bar[, , k] = gamma_bar[, , k] / gb_sums
  }

  alpha = array(rgamma(N * D * K, 1, 1), dim = c(N, D, K))
  lambda = array(0, dim = c(N, D, K))
  for (k in 1:K) lambda[, , k] = lambda0 + F_sum[k]

  xi = array(1 / D, dim = c(N, D, M))
  elbo_hist = rep(NA_real_, max_iters)
  last_elbo = NA_real_
  small_improve_count = 0
  converged = FALSE
  n_iter = max_iters

  for (iter in 1:max_iters) {
    rho = min(1, gamma_step_init + (1 - gamma_step_init) * (iter - 1) /
                max(gamma_step_ramp - 1, 1))
    rho_F = min(1, F_step_init + (1 - F_step_init) * (iter - 1) /
                  max(F_step_ramp - 1, 1))

    if (break_symmetry && D > 1) {
      unif = rep(1 / K, K)
      unif_mat_N = matrix(unif, N, K, byrow = TRUE)
      is_active = matrix(TRUE, N, D)
      for (dd in 1:D) {
        is_active[, dd] = sqrt(rowSums((gamma_bar[, dd, ] - unif_mat_N)^2)) > tol_dist_sim
      }
      repeat {
        any_merged = FALSE
        for (d1 in 1:(D - 1)) {
          for (d2 in (d1 + 1):D) {
            pdist = sqrt(rowSums((gamma_bar[, d1, ] - gamma_bar[, d2, ])^2))
            to_merge = which(pdist < tol_dist_sim & is_active[, d1] & is_active[, d2])
            if (length(to_merge) > 0) {
              any_merged = TRUE
              gamma_bar[to_merge, d1, ] = (gamma_bar[to_merge, d1, ] + gamma_bar[to_merge, d2, ]) / 2
              gamma_bar[to_merge, d2, ] = unif_mat_N[to_merge, , drop = FALSE]
              is_active[to_merge, d2] = FALSE
            }
          }
        }
        if (!any_merged) break
      }
    }

    for (d in 1:D) {
      ## Gauss-Seidel-style xi refresh using all effects' latest parameters.
      xi = update_mf_xi(F, gamma_bar, alpha, lambda)

      ## update effect d for all observations
      xi_d_Y = xi[, d, ] * Y  # (N, M)
      xi_d_Y_sum = rowSums(xi_d_Y)  # (N,)

      alpha_id = xi_d_Y_sum + alpha0[, d]  # (N,)
      for (k in 1:K) alpha[, d, k] = alpha_id
      for (k in 1:K) lambda[, d, k] = lambda0[, d] + F_sum[k]

      # logprob: (N, K)
      logprob = xi_d_Y %*% t(log_F) - alpha_id * log(outer(lambda0[, d], F_sum, "+"))
      logprob = logprob - apply(logprob, 1, max)
      gamma_bar_cavi = exp(logprob)
      gamma_bar_cavi = gamma_bar_cavi / rowSums(gamma_bar_cavi)
      gamma_bar[, d, ] = (1 - rho) * gamma_bar[, d, ] + rho * gamma_bar_cavi

      if (update_prior) {
        E_beta_bar = rowSums(gamma_bar[, d, ] * alpha[, d, ] / lambda[, d, ])  # (N,)
        E_log_beta_bar = rowSums(gamma_bar[, d, ] *
                                   (digamma(alpha[, d, ]) - log(lambda[, d, ])))  # (N,)
        target = E_log_beta_bar + log(lambda0[, d])  # (N,)
        alpha0[, d] = inv_digamma(target)
        lambda0[, d] = alpha0[, d] / E_beta_bar
      }
    }

    ## Refresh xi after the last gamma/alpha/lambda update before using it for
    ## the shared-F M-step or objective evaluation.
    xi = update_mf_xi(F, gamma_bar, alpha, lambda)

    if (update_F) {
      ## F_{km} ∝ A_{km} = Σ_{i,d} γ̄_{idk} ξ_{imd} y_{im}
      A = matrix(0, K, M)
      for (dd in 1:D) {
        A = A + t(gamma_bar[, dd, ]) %*% (xi[, dd, ] * Y)  # (K,N) %*% (N,M) = (K,M)
      }
      A = A + F_pseudocount
      F_cavi = normalize_rows(A)
      F = normalize_rows((1 - rho_F) * F + rho_F * F_cavi)
      log_F = log(F)
      F_sum = rowSums(F)  # rep(1, K)

      xi = update_mf_xi(F, gamma_bar, alpha, lambda)
    }

    if (iter == 1 || iter == max_iters || iter %% elbo_every == 0) {
      elbo_hist[iter] = mf_ELBO(Y, F, xi, gamma_bar, alpha, lambda, alpha0, lambda0)
      if (!is.null(tol) && !is.na(last_elbo)) {
        rel_improve = (elbo_hist[iter] - last_elbo) / (abs(last_elbo) + 1)
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
      last_elbo = elbo_hist[iter]
    }
  }

  res = list(gamma_bar = gamma_bar, alpha = alpha, lambda = lambda,
             alpha0 = alpha0, lambda0 = lambda0, xi = xi, F = F,
             elbo = elbo_hist, elbo_iter = which(!is.na(elbo_hist)),
             converged = converged, n_iter = n_iter)
  return(res)
}

mf_poi_susie <- function(Y, K, D, max_iters = 100, update_prior = TRUE,
                         break_symmetry = FALSE, tol_dist_sim = 1e-4,
                         init_seed = NULL, prior_shape = 1, prior_rate = 1,
                         init_F = c("poisson_nmf", "random"),
                         nmf_iters = 200, init_gamma_from_nmf = TRUE,
                         gamma_init_floor = 1e-4,
                         gamma_step_init = 0.1, gamma_step_ramp = 50,
                         F_step_init = 0.05, F_step_ramp = 75,
                         F_pseudocount = .Machine$double.eps,
                         elbo_every = 1, tol = NULL, min_iters = 10,
                         patience = 3) {
  N = nrow(Y)
  M = ncol(Y)
  init_F = match.arg(init_F)

  alpha0 = matrix(prior_shape, N, D)
  lambda0 = matrix(prior_rate, N, D)
  F = matrix(1 / M, nrow = K, ncol = M)
  fit_init_F = "poisson_nmf"

  if (init_F == "random") {
    if (!is.null(init_seed)) set.seed(init_seed)
    F = matrix(rgamma(K * M, shape = 1, rate = 1), nrow = K, ncol = M)
    F = normalize_rows(F)
    fit_init_F = "given"
    init_gamma_from_nmf = FALSE
  }

  .mf_poi_susie_fit(
    Y = Y,
    F = F,
    alpha0 = alpha0,
    lambda0 = lambda0,
    max_iters = max_iters,
    update_prior = update_prior,
    update_F = TRUE,
    break_symmetry = break_symmetry,
    tol_dist_sim = tol_dist_sim,
    init_seed = init_seed,
    gamma_step_init = gamma_step_init,
    gamma_step_ramp = gamma_step_ramp,
    init_F = fit_init_F,
    nmf_iters = nmf_iters,
    init_gamma_from_nmf = init_gamma_from_nmf,
    gamma_init_floor = gamma_init_floor,
    F_step_init = F_step_init,
    F_step_ramp = F_step_ramp,
    F_pseudocount = F_pseudocount,
    elbo_every = elbo_every,
    tol = tol,
    min_iters = min_iters,
    patience = patience
  )
}

mf_poi_susie_fixed_F <- function(Y, F, D, max_iters = 100,
                                 update_prior = TRUE,
                                 break_symmetry = FALSE,
                                 tol_dist_sim = 1e-4,
                                 init_seed = NULL, prior_shape = 1,
                                 prior_rate = 1,
                                 gamma_step_init = 0.1,
                                 gamma_step_ramp = 50,
                                 elbo_every = 1, tol = NULL,
                                 min_iters = 10, patience = 3) {
  N = nrow(Y)
  alpha0 = matrix(prior_shape, N, D)
  lambda0 = matrix(prior_rate, N, D)

  .mf_poi_susie_fit(
    Y = Y,
    F = F,
    alpha0 = alpha0,
    lambda0 = lambda0,
    max_iters = max_iters,
    update_prior = update_prior,
    update_F = FALSE,
    break_symmetry = break_symmetry,
    tol_dist_sim = tol_dist_sim,
    init_seed = init_seed,
    gamma_step_init = gamma_step_init,
    gamma_step_ramp = gamma_step_ramp,
    init_F = "given",
    init_gamma_from_nmf = FALSE,
    elbo_every = elbo_every,
    tol = tol,
    min_iters = min_iters,
    patience = patience
  )
}
