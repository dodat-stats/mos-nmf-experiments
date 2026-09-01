## Poisson SuSiE for one observation.
##
## Notation follows the manuscript: lambda is the nonnegative latent loading,
## while beta is the rate of its Gamma distribution.

poisson_susie_elbo <- function(y, F, xi, gamma_bar, alpha, beta,
                               alpha0, beta0) {
  D = nrow(gamma_bar)
  K = ncol(gamma_bar)

  E_log_lambda_given_gamma = digamma(alpha) - log(beta)
  E_log_lambda = rowSums(gamma_bar * E_log_lambda_given_gamma)

  term_E_llh_1 = sum(E_log_lambda * (xi %*% y))

  E_log_F = gamma_bar %*% log(F)
  term_E_llh_2 = sum(sweep(xi * E_log_F, 2, y, "*"))

  E_lambda = alpha / beta
  term_E_llh_3 = -sum((gamma_bar * E_lambda) %*% F)

  xi_log_xi = xi * log(xi)
  xi_log_xi[is.nan(xi_log_xi)] = 0
  term_E_llh_4 = -sum(sweep(xi_log_xi, 2, y, "*"))

  gb_log_gb = gamma_bar * log(gamma_bar)
  gb_log_gb[is.nan(gb_log_gb)] = 0
  term_KL_discrete = sum(gb_log_gb) + D * log(K)

  alpha0_mat = matrix(alpha0, nrow = D, ncol = K)
  beta0_mat = matrix(beta0, nrow = D, ncol = K)
  term_KL_cont = sum(gamma_bar * (
    (alpha - alpha0_mat) * digamma(alpha) +
      (log(beta) - log(beta0_mat)) * alpha0_mat -
      (lgamma(alpha) - lgamma(alpha0_mat)) -
      (beta - beta0_mat) * alpha / beta
  ))

  term_E_llh_1 + term_E_llh_2 + term_E_llh_3 + term_E_llh_4 -
    term_KL_discrete - term_KL_cont
}

#' Fit single-observation Poisson SuSiE
#'
#' `beta0` is the Gamma prior rate. The returned `beta` and `beta0` fields
#' always denote Gamma rates; the latent loading itself is lambda.
poisson_susie <- function(y, F, alpha0, beta0, max_iters = 100,
                          update_prior = TRUE, break_symmetry = FALSE,
                          tol_dist_sim = 1e-4, init_seed = NULL) {
  K = nrow(F)
  M = ncol(F)
  D = length(alpha0)

  if (length(beta0) != D) stop("alpha0 and beta0 must have the same length")
  if (length(y) != M) stop("length(y) must equal ncol(F)")
  if (any(F <= 0)) stop("F must be strictly positive")

  if (!is.null(init_seed)) set.seed(init_seed)
  gamma_bar = matrix(rgamma(D * K, shape = 1, rate = 1),
                     nrow = D, ncol = K)
  gamma_bar = gamma_bar / rowSums(gamma_bar)
  alpha = matrix(rgamma(D * K, shape = 1, rate = 1),
                 nrow = D, ncol = K)
  beta = matrix(rgamma(D * K, shape = 1, rate = 1),
                nrow = D, ncol = K)
  xi = matrix(1 / D, nrow = D, ncol = M)
  elbo_hist = numeric(max_iters)

  for (iter in seq_len(max_iters)) {
    if (break_symmetry && D > 1) {
      unif = rep(1 / K, K)
      dist_to_unif = sqrt(rowSums((gamma_bar - rep(1, D) %o% unif)^2))
      active = which(dist_to_unif > tol_dist_sim)
      if (length(active) >= 2) {
        dists = as.matrix(dist(gamma_bar[active, , drop = FALSE]))
        diag(dists) = Inf
        while (min(dists) < tol_dist_sim) {
          idx = which(dists == min(dists), arr.ind = TRUE)[1, ]
          d1 = active[idx[1]]
          d2 = active[idx[2]]
          gamma_bar[d1, ] = (gamma_bar[d1, ] + gamma_bar[d2, ]) / 2
          gamma_bar[d2, ] = unif
          active = active[active != d2]
          if (length(active) < 2) break
          dists = as.matrix(dist(gamma_bar[active, , drop = FALSE]))
          diag(dists) = Inf
        }
      }
    }

    for (d in seq_len(D)) {
      E_log_lambda_given_gamma = digamma(alpha) - log(beta)
      E_log_lambda = rowSums(gamma_bar * E_log_lambda_given_gamma)
      E_log_F_gamma_d = gamma_bar %*% log(F)

      log_xi = sweep(E_log_F_gamma_d, 1, E_log_lambda, "+")
      log_xi = sweep(log_xi, 2, apply(log_xi, 2, max), "-")
      xi = exp(log_xi)
      xi = sweep(xi, 2, colSums(xi), "/")

      allocated_count = sum(xi[d, ] * y)
      alpha[d, ] = rep(allocated_count + alpha0[d], K)
      beta[d, ] = rowSums(F) + beta0[d]

      log_gamma_bar_d = log(F) %*% (xi[d, ] * y) -
        (allocated_count + alpha0[d]) * log(rowSums(F) + beta0[d])
      log_gamma_bar_d = log_gamma_bar_d - max(log_gamma_bar_d)
      gamma_bar_d = exp(log_gamma_bar_d)
      gamma_bar[d, ] = gamma_bar_d / sum(gamma_bar_d)

      if (update_prior) {
        E_lambda_bar = sum(gamma_bar[d, ] * alpha[d, ] / beta[d, ])
        E_log_lambda_bar = sum(
          gamma_bar[d, ] * (digamma(alpha[d, ]) - log(beta[d, ]))
        )
        target = E_log_lambda_bar + log(beta0[d])
        alpha0[d] = uniroot(
          function(a) digamma(a) - target,
          lower = 1e-8,
          upper = 1e8
        )$root
        beta0[d] = alpha0[d] / E_lambda_bar
      }
    }

    elbo_hist[iter] = poisson_susie_elbo(
      y, F, xi, gamma_bar, alpha, beta, alpha0, beta0
    )
  }

  list(
    gamma_bar = gamma_bar,
    alpha = alpha,
    beta = beta,
    alpha0 = alpha0,
    beta0 = beta0,
    xi = xi,
    elbo = elbo_hist
  )
}

## Compatibility wrapper for analyses written before the public API was frozen.
poi_susie <- function(y, F, alpha0, beta0 = NULL, lambda0 = NULL, ...) {
  if (is.null(beta0)) beta0 = lambda0
  if (is.null(beta0)) stop("beta0 must be supplied")
  poisson_susie(y = y, F = F, alpha0 = alpha0, beta0 = beta0, ...)
}

one_data_ELBO <- poisson_susie_elbo
