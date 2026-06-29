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

## matrix factorization Poisson-SuSiE
# Y: (N, M) count matrix
# F: (K, M) shared factor matrix
# alpha0: (N, D) prior shape
# lambda0: (N, D) prior rate
mf_poi_susie <- function(Y, F, alpha0, lambda0, max_iters = 100, update_prior = TRUE,
                         update_F = FALSE, break_symmetry = FALSE,
                         fine_init_gamma = FALSE,
                         tol_dist_sim = 1e-4, init_seed = NULL,
                         step_init = 1, step_ramp = 50) {
  N = nrow(Y)
  M = ncol(Y)
  K = nrow(F)
  D = ncol(alpha0)

  if (update_F) F = F / rowSums(F)
  log_F = log(F)
  F_sum = rowSums(F)  # (K,)

  ## initialization
  if (!is.null(init_seed)) set.seed(init_seed)
  if (fine_init_gamma) {
    gamma_bar = array(0, dim = c(N, D, K))
    Ft = t(F)  # (M, K)
    for (i in 1:N) {
      fit_i = susieR::susie(Ft, Y[i, ], L = D, intercept = FALSE,
                            estimate_residual_variance = FALSE,
                            residual_variance = max(var(Y[i, ]), 1))
      pip = pmax(colSums(fit_i$alpha), 0)  # (K,) combined evidence per factor
      gb_i = matrix(0, D, K)
      for (dd in 1:D) {
        k_best = which.max(pip)
        gb_i[dd, k_best] = 1
        pip[k_best] = -Inf
      }
      gamma_bar[i, , ] = gb_i
    }
  } else {
    gamma_bar = array(rgamma(N * D * K, 1, 1), dim = c(N, D, K))
    gb_sums = rowSums(gamma_bar, dims = 2)  # (N, D)
    for (k in 1:K) gamma_bar[, , k] = gamma_bar[, , k] / gb_sums
  }

  alpha = array(rgamma(N * D * K, 1, 1), dim = c(N, D, K))
  lambda = array(0, dim = c(N, D, K))
  for (k in 1:K) lambda[, , k] = lambda0 + F_sum[k]

  xi = array(1 / D, dim = c(N, D, M))
  elbo_hist = numeric(max_iters)

  for (iter in 1:max_iters) {
    rho = min(1, step_init + (1 - step_init) * (iter - 1) / max(step_ramp - 1, 1))

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
      ## update xi for all observations
      E_log_beta_cond = digamma(alpha) - log(lambda)  # (N, D, K)
      E_log_beta = rowSums(gamma_bar * E_log_beta_cond, dims = 2)  # (N, D)

      for (dd in 1:D) {
        xi[, dd, ] = gamma_bar[, dd, ] %*% log_F + E_log_beta[, dd]
      }
      log_xi_max = xi[, 1, ]
      for (dd in 2:D) log_xi_max = pmax(log_xi_max, xi[, dd, ])
      for (dd in 1:D) xi[, dd, ] = xi[, dd, ] - log_xi_max
      xi = exp(xi)
      xi_sum = xi[, 1, ]
      if (D > 1) for (dd in 2:D) xi_sum = xi_sum + xi[, dd, ]
      for (dd in 1:D) xi[, dd, ] = xi[, dd, ] / xi_sum

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

    if (update_F) {
      ## F_{km} ∝ A_{km} = Σ_{i,d} γ̄_{idk} ξ_{imd} y_{im}
      A = matrix(0, K, M)
      for (dd in 1:D) {
        A = A + t(gamma_bar[, dd, ]) %*% (xi[, dd, ] * Y)  # (K,N) %*% (N,M) = (K,M)
      }
      A = A + .Machine$double.eps
      F = A / rowSums(A)
      log_F = log(F)
      F_sum = rowSums(F)  # rep(1, K)
    }

    elbo_hist[iter] = mf_ELBO(Y, F, xi, gamma_bar, alpha, lambda, alpha0, lambda0)
  }

  res = list(gamma_bar = gamma_bar, alpha = alpha, lambda = lambda,
             alpha0 = alpha0, lambda0 = lambda0, xi = xi, F = F, elbo = elbo_hist)
  return(res)
}
