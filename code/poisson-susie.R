## poisson-susie-ELBO-fn (eq 12, 14)
# y: (M,)
# F: (K, M)
# xi: (D, M)
# gamma_bar: (D, K)
# alpha: (D, K)
# lambda: (D, K)
# alpha0: (D,)
# lambda0: (D,)
one_data_ELBO <- function(y, F, xi, gamma_bar, alpha, lambda, alpha0, lambda0){
  D = nrow(gamma_bar)
  K = ncol(gamma_bar)

  E_log_beta_given_gamma = digamma(alpha) - log(lambda)  # (D, K)
  E_log_beta = rowSums(gamma_bar * E_log_beta_given_gamma)  # (D,)

  # Σ_{d,m} ξ_md y_m E[log β_d]
  term_E_llh_1 = sum(E_log_beta * (xi %*% y))

  # Σ_{d,m} ξ_md y_m E[log F_{γ_d m}]
  E_log_F = gamma_bar %*% log(F)  # (D, M)
  term_E_llh_2 = sum(sweep(xi * E_log_F, 2, y, "*"))

  # -Σ_{d,m} E[β_d F_{γ_d m}]
  E_beta = alpha / lambda  # (D, K)
  term_E_llh_3 = -sum((gamma_bar * E_beta) %*% F)

  # -Σ_{m,d} y_m ξ_md log(ξ_md)
  xi_log_xi = xi * log(xi)
  xi_log_xi[is.nan(xi_log_xi)] = 0
  term_E_llh_4 = -sum(sweep(xi_log_xi, 2, y, "*"))

  # KL(q_d || p_d) — eq 14
  gb_log_gb = gamma_bar * log(gamma_bar)
  gb_log_gb[is.nan(gb_log_gb)] = 0
  term_KL_discrete = sum(gb_log_gb) + D * log(K)

  alpha0_mat = matrix(alpha0, nrow = D, ncol = K)
  lambda0_mat = matrix(lambda0, nrow = D, ncol = K)
  term_KL_cont = sum(gamma_bar * (
    (alpha - alpha0_mat) * digamma(alpha)
    + (log(lambda) - log(lambda0_mat)) * alpha0_mat
    - (lgamma(alpha) - lgamma(alpha0_mat))
    - (lambda - lambda0_mat) * alpha / lambda
  ))

  ELBO = term_E_llh_1 + term_E_llh_2 + term_E_llh_3 + term_E_llh_4
         - term_KL_discrete - term_KL_cont
  return(ELBO)
}

## poisson-susie
poi_susie <- function(y, F, alpha0, lambda0, max_iters=100, update_prior=TRUE,
                      break_symmetry=FALSE, tol_dist_sim=1e-4, init_seed=NULL){
  K = dim(F)[1]
  M = dim(F)[2]
  D = length(alpha0)


  ## initialization
  if (!is.null(init_seed)) set.seed(init_seed)
  gamma_bar = matrix(rgamma(D * K, shape = 1, rate = 1), nrow = D, ncol = K)
  gamma_bar = gamma_bar / rowSums(gamma_bar) # (D, K)
  alpha = matrix(rgamma(D * K, shape = 1, rate = 1), nrow = D, ncol = K) # (D, K)
  lambda = matrix(rgamma(D * K, shape = 1, rate = 1), nrow = D, ncol = K) # (D, K)
  elbo_hist = numeric(max_iters)
  for (iter in c(1:max_iters)){

    if (break_symmetry && D > 1) {
      unif = rep(1 / K, K)
      dist_to_unif = sqrt(rowSums((gamma_bar - rep(1, D) %o% unif)^2))
      active = which(dist_to_unif > tol_dist_sim)
      if (length(active) >= 2) {
        dists = as.matrix(dist(gamma_bar[active, , drop = FALSE], method = "euclidean"))
        diag(dists) = Inf
        while (min(dists) < tol_dist_sim) {
          idx = which(dists == min(dists), arr.ind = TRUE)[1, ]
          d1 = active[idx[1]]; d2 = active[idx[2]]
          gamma_bar[d1, ] = (gamma_bar[d1, ] + gamma_bar[d2, ]) / 2
          gamma_bar[d2, ] = unif
          active = active[active != d2]
          if (length(active) < 2) break
          dists = as.matrix(dist(gamma_bar[active, , drop = FALSE], method = "euclidean"))
          diag(dists) = Inf
        }
      }
    }

    for (d in c(1:D)){
      ## update xi using all effects' latest params
      E_log_beta_given_gamma = digamma(alpha) - log(lambda)  # (D, K)
      E_log_beta = rowSums(gamma_bar * E_log_beta_given_gamma)  # (D,)
      E_log_F_gamma_d = gamma_bar %*% log(F) # (D, M)

      log_xi = sweep(E_log_F_gamma_d, 1, E_log_beta, "+") # (D, M)
      log_xi = sweep(log_xi, 2, apply(log_xi, 2, max), "-")
      xi = exp(log_xi)
      xi = sweep(xi, 2, colSums(xi), "/")

      ## update effect d
      alpha[d, ] = rep(sum(xi[d, ] * y) + alpha0[d], K)
      lambda[d, ] = rowSums(F) + lambda0[d]

      log_gamma_bar_d = (log(F) %*% (xi[d, ] * y)
                            - (sum(xi[d, ] * y) + alpha0[d]) * log(rowSums(F) + lambda0[d]))
      log_gamma_bar_d = log_gamma_bar_d - max(log_gamma_bar_d)
      gamma_bar_d = exp(log_gamma_bar_d)
      gamma_bar[d, ] = gamma_bar_d / sum(gamma_bar_d)

      if (update_prior) {
        ## update prior parameters (VEB)
        E_beta_bar = sum(gamma_bar[d, ] * alpha[d, ] / lambda[d, ])
        E_log_beta_bar = sum(gamma_bar[d, ] * (digamma(alpha[d, ]) - log(lambda[d, ])))
        target = E_log_beta_bar + log(lambda0[d])
        alpha0[d] = uniroot(function(a) digamma(a) - target,
                             lower = 1e-8, upper = 1e8)$root
        lambda0[d] = alpha0[d] / E_beta_bar
      }

    }
    elbo_hist[iter] = one_data_ELBO(y, F, xi, gamma_bar, alpha, lambda, alpha0, lambda0)
  }

  res = list(gamma_bar = gamma_bar,
             alpha = alpha,
             lambda = lambda,
             alpha0 = alpha0,
             lambda0 = lambda0,
             xi = xi,
             elbo = elbo_hist)
  return(res)
}
