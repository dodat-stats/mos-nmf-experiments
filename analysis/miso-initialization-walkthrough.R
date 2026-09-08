## MiSo initialization walkthrough ------------------------------------------
##
## This script is deliberately self-contained. It uses only base R and keeps
## all important intermediate objects in the global environment so that they
## are easy to inspect in RStudio.
##
## Suggested workflow:
##   1. Change the settings below.
##   2. Run one section at a time with Cmd/Ctrl + Shift + Enter.
##   3. Inspect objects such as Y, nmf_fit, ps_fit,
##      miso_initializations, and miso_fits with View() or str().


## 1. User settings ---------------------------------------------------------

## These are fitted dimensions. The simulated truth remains K = S = D = 1.
K_fit = 2L
S_fit = 1L
D_fit = 2L

## Initialization to leave in the convenient object `miso_fit`.
## Choices are "distinct", "threshold_recycle", and "aligned_ps".
initialization_to_inspect = "aligned_ps"

## TRUE fits the same MiSo model from every compatible initialization.
## FALSE fits only `initialization_to_inspect`.
run_all_initializations = TRUE

## Simulation and fitting controls. These defaults are intentionally small.
N = 500L
M = 1000L
simulation_seed = 1L
nmf_seed = 2L
initialization_seed = 3L
nmf_iterations = 300L
poisson_susie_iterations = 60L
miso_iterations = 80L
miso_inner_iterations = 5L

## Initialization controls.
minimum_factor_share = 0.10
gamma_floor = 0.05

## Keeping F fixed makes the effect of gamma initialization easier to study.
## Set this to TRUE if you also want MiSo to update the factor dictionary.
update_F_in_miso = FALSE

## Optional empirical-Bayes updates for the Gamma loading priors.
update_poisson_susie_priors = FALSE
update_miso_priors = TRUE

stopifnot(K_fit >= 1L, S_fit >= 1L, D_fit >= 1L)
stopifnot(K_fit <= N, S_fit <= N)
allowed_initializations = c("distinct", "threshold_recycle", "aligned_ps")
stopifnot(initialization_to_inspect %in% allowed_initializations)


## 2. Small numerical helpers ----------------------------------------------

softmax = function(x) {
  x = x - max(x)
  answer = exp(x)
  answer / sum(answer)
}

softmax_rows = function(log_weights) {
  log_weights = as.matrix(log_weights)
  answer = matrix(0, nrow(log_weights), ncol(log_weights))
  for (i in seq_len(nrow(log_weights))) {
    answer[i, ] = softmax(log_weights[i, ])
  }
  answer
}

log_sum_exp = function(x) {
  largest = max(x)
  largest + log(sum(exp(x - largest)))
}

normalize_rows = function(x, epsilon = 1e-12) {
  x = pmax(x, epsilon)
  x / pmax(rowSums(x), epsilon)
}

gamma_shape_from_moments = function(expected_lambda, expected_log_lambda,
                                    n_steps = 8L,
                                    minimum_shape = 1e-3,
                                    maximum_shape = 1e4) {
  ## Solve log(a) - digamma(a) = log(E[lambda]) - E[log(lambda)].
  delta = log(pmax(expected_lambda, 1e-12)) - expected_log_lambda
  delta = pmax(delta, 1e-12)
  shape = if (delta > 0.5) 1 / (2 * delta) else 1 / delta
  shape = min(max(shape, minimum_shape), maximum_shape)

  for (step in seq_len(n_steps)) {
    value = log(shape) - digamma(shape) - delta
    derivative = 1 / shape - trigamma(shape)
    shape = shape - value / derivative
    shape = min(max(shape, minimum_shape), maximum_shape)
  }
  shape
}

poisson_deviance_per_entry = function(Y, fitted_mean) {
  fitted_mean = pmax(fitted_mean, 1e-12)
  contribution = fitted_mean - Y
  positive = Y > 0
  contribution[positive] = contribution[positive] +
    Y[positive] * log(Y[positive] / fitted_mean[positive])
  2 * sum(contribution) / length(Y)
}


## 3. Generate Y from the simplest model: true K = S = D = 1 ------------

set.seed(simulation_seed)

## A single smooth, strictly positive, row-normalized true factor.
feature_position = seq(-3, 3, length.out = M)
F_true = matrix(dnorm(feature_position, mean = 0, sd = 0.9) + 0.01,
                nrow = 1L)
F_true = normalize_rows(F_true)

## With one factor and one loading, E[sum_m Y_im | L_i] = L_i.
L_true = matrix(rgamma(N, shape = 6, rate = 0.10), ncol = 1L)
true_mean = L_true %*% F_true
Y = matrix(rpois(N * M, lambda = as.vector(true_mean)), N, M)

cat("Simulated Y with dimensions", N, "x", M, "\n")
cat("True dimensions: K = S = D = 1\n")
cat("Fitted dimensions: K =", K_fit,
    ", S =", S_fit, ", D =", D_fit, "\n\n")


## 4. Fit vanilla Poisson/KL NMF -------------------------------------------

fit_poisson_nmf = function(Y, K, n_iterations = 300L, seed = 1L,
                           epsilon = 1e-10) {
  ## Multiplicative updates minimize generalized KL divergence, equivalently
  ## maximizing the Poisson likelihood up to constants.
  set.seed(seed)
  N = nrow(Y)
  M = ncol(Y)

  L = matrix(rexp(N * K, rate = 1), N, K)
  F = matrix(rexp(K * M, rate = 1), K, M)
  F = normalize_rows(F, epsilon)
  L = L * matrix(rowSums(Y) / pmax(rowSums(L), epsilon), N, K)
  objective = numeric(n_iterations)

  for (iteration in seq_len(n_iterations)) {
    fitted_mean = L %*% F + epsilon
    L = L * ((Y / fitted_mean) %*% t(F)) /
      matrix(rowSums(F), N, K, byrow = TRUE)
    L = pmax(L, epsilon)

    fitted_mean = L %*% F + epsilon
    F = F * (t(L) %*% (Y / fitted_mean)) /
      matrix(colSums(L), K, M)
    F = pmax(F, epsilon)

    ## Remove the NMF scale ambiguity by making every row of F sum to one.
    factor_scale = rowSums(F)
    F = F / factor_scale
    L = sweep(L, 2, factor_scale, "*")

    fitted_mean = L %*% F + epsilon
    objective[iteration] = sum(fitted_mean - Y * log(fitted_mean))
  }

  list(
    L = L,
    F = F,
    fitted_mean = L %*% F,
    objective = objective
  )
}

nmf_fit = fit_poisson_nmf(
  Y = Y,
  K = K_fit,
  n_iterations = nmf_iterations,
  seed = nmf_seed
)

cat("Vanilla NMF deviance per matrix entry:",
    signif(poisson_deviance_per_entry(Y, nmf_fit$fitted_mean), 4), "\n\n")


#### I want to understand relationship between fit and true
sum(F_true)
sum(nmf_fit$F[1, ])
sum(nmf_fit$F[2, ])
sum((nmf_fit$F[1, ] - F_true)^2)
sum((nmf_fit$F[2, ] - F_true)^2)

(nmf_fit$L[, 1] + nmf_fit$L[, 2] - L_true) / L_true


## 5. Row-wise Poisson SuSiE with the fitted F frozen -----------------------
poisson_susie_one_row = function(y, F, D, n_iterations = 60L,
                                 update_prior = FALSE,
                                 initial_factor_order = NULL,
                                 epsilon = 1e-12) {
  ## Variational objects for this observation:
  ##   gamma_bar[d,k] = q(the factor selected by slot d is k)
  ##   Gamma(alpha[d,k], beta[d,k]) is q(lambda_d | selected factor k)
  ##   xi[d,m] allocates observed count y[m] to slot d.
  K = nrow(F)
  M = ncol(F)
  factor_exposure = rowSums(F)
  log_F = log(pmax(F, epsilon))

  alpha0 = rep(1, D)
  beta0 = rep(1 / pmax(sum(y) / D, 1), D)
  gamma_bar = matrix(epsilon, D, K)

  if (is.null(initial_factor_order)) {
    initial_factor_order = seq_len(K)
  }
  initial_factor_order = rep(initial_factor_order, length.out = D)
  for (d in seq_len(D)) gamma_bar[d, initial_factor_order[d]] = 1
  gamma_bar = normalize_rows(gamma_bar)

  alpha = matrix(alpha0 + sum(y) / D, D, K)
  beta = matrix(beta0, D, K) +
    matrix(factor_exposure, D, K, byrow = TRUE)
  xi = matrix(1 / D, D, M)

  for (iteration in seq_len(n_iterations)) {
    ## Coordinate ascent over the D single-effect slots.
    for (d in seq_len(D)) {
      expected_log_lambda = rowSums(
        gamma_bar * (digamma(alpha) - log(beta))
      )
      expected_log_factor = gamma_bar %*% log_F

      log_xi = expected_log_factor +
        matrix(expected_log_lambda, D, M)
      log_xi = sweep(log_xi, 2, apply(log_xi, 2, max), "-")
      xi = exp(log_xi)
      xi = sweep(xi, 2, colSums(xi), "/")

      allocated_count = sum(xi[d, ] * y)
      alpha[d, ] = alpha0[d] + allocated_count
      beta[d, ] = beta0[d] + factor_exposure

      log_probability = as.vector(log_F %*% (xi[d, ] * y)) -
        (alpha0[d] + allocated_count) *
        log(beta0[d] + factor_exposure)
      gamma_bar[d, ] = softmax(log_probability)

      if (update_prior) {
        expected_lambda = sum(gamma_bar[d, ] * alpha[d, ] / beta[d, ])
        expected_log_lambda = sum(
          gamma_bar[d, ] * (digamma(alpha[d, ]) - log(beta[d, ]))
        )
        alpha0[d] = gamma_shape_from_moments(
          expected_lambda, expected_log_lambda
        )
        beta0[d] = alpha0[d] / pmax(expected_lambda, epsilon)
      }
    }
  }

  list(
    gamma_bar = gamma_bar,
    alpha = alpha,
    beta = beta,
    alpha0 = alpha0,
    beta0 = beta0,
    xi = xi
  )
}

fit_rowwise_poisson_susie = function(Y, F, D, n_iterations = 60L,
                                     update_prior = FALSE) {
  N = nrow(Y)
  K = nrow(F)
  M = ncol(F)

  gamma_bar = array(0, c(N, D, K))
  alpha = array(0, c(N, D, K))
  beta = array(0, c(N, D, K))
  alpha0 = matrix(0, N, D)
  beta0 = matrix(0, N, D)
  xi = array(0, c(N, D, M))

  ## NMF loadings provide only a symmetry-breaking order. F stays fixed.
  nmf_loading_order = matrix(0L, N, K)
  for (i in seq_len(N)) {
    nmf_loading_order[i, ] = order(nmf_fit$L[i, ], decreasing = TRUE)
  }

  for (i in seq_len(N)) {
    row_fit = poisson_susie_one_row(
      y = Y[i, ],
      F = F,
      D = D,
      n_iterations = n_iterations,
      update_prior = update_prior,
      initial_factor_order = nmf_loading_order[i, ]
    )
    gamma_bar[i, , ] = row_fit$gamma_bar
    alpha[i, , ] = row_fit$alpha
    beta[i, , ] = row_fit$beta
    alpha0[i, ] = row_fit$alpha0
    beta0[i, ] = row_fit$beta0
    xi[i, , ] = row_fit$xi
  }

  ## Expected loading assigned to each fitted NMF factor.
  expected_factor_loading = matrix(0, N, K)
  for (i in seq_len(N)) {
    for (d in seq_len(D)) {
      expected_factor_loading[i, ] = expected_factor_loading[i, ] +
        gamma_bar[i, d, ] * alpha[i, d, ] / beta[i, d, ]
    }
  }

  list(
    gamma_bar = gamma_bar,
    alpha = alpha,
    beta = beta,
    alpha0 = alpha0,
    beta0 = beta0,
    xi = xi,
    expected_factor_loading = expected_factor_loading
  )
}

ps_fit = fit_rowwise_poisson_susie(
  Y = Y,
  F = nmf_fit$F,
  D = D_fit,
  n_iterations = poisson_susie_iterations,
  update_prior = update_poisson_susie_priors
)

cat("Finished", N, "row-wise Poisson SuSiE fits.\n")
cat("Inspect ps_fit$gamma_bar[i,,] to see observation i's slot posteriors.\n\n")


## 6. Cluster observations before constructing MiSo initializations --------

cluster_loading_scores = function(scores, S, seed = 1L,
                                   epsilon = 1e-12) {
  normalized_scores = scores / pmax(rowSums(scores), epsilon)

  if (S == 1L) {
    cluster = rep(1L, nrow(scores))
  } else {
    ## If the normalized rows are identical (for example K = 1), composition
    ## cannot define clusters. In that case split by total expected loading.
    number_distinct = nrow(unique(round(normalized_scores, 10)))
    if (number_distinct < S) {
      rank_by_size = rank(rowSums(scores), ties.method = "first")
      cluster = pmin(S, ceiling(S * rank_by_size / nrow(scores)))
    } else {
      set.seed(seed)
      cluster = kmeans(
        normalized_scores, centers = S, nstart = 20, iter.max = 100
      )$cluster
    }
  }

  centers = matrix(0, S, ncol(scores))
  for (s in seq_len(S)) {
    centers[s, ] = colMeans(normalized_scores[cluster == s, , drop = FALSE])
  }
  centers = normalize_rows(centers)

  list(cluster = as.integer(cluster), centers = centers,
       normalized_scores = normalized_scores)
}

preliminary_clusters = cluster_loading_scores(
  scores = ps_fit$expected_factor_loading,
  S = S_fit,
  seed = initialization_seed
)


## 7. Three MiSo initializations -------------------------------------------

gamma_from_factor_indices = function(factor_index, K, floor = 0.05) {
  S = nrow(factor_index)
  D = ncol(factor_index)
  gamma_bar = array(floor / K, c(S, D, K))
  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      gamma_bar[s, d, factor_index[s, d]] =
        1 - floor + floor / K
    }
  }
  gamma_bar
}

initialize_distinct = function(cluster_fit, D, floor = 0.05) {
  K = ncol(cluster_fit$centers)
  S = nrow(cluster_fit$centers)
  if (D > K) stop("Distinct initialization requires D_fit <= K_fit.")

  factor_index = matrix(0L, S, D)
  for (s in seq_len(S)) {
    factor_index[s, ] = order(
      cluster_fit$centers[s, ], decreasing = TRUE
    )[seq_len(D)]
  }

  list(
    gamma_bar = gamma_from_factor_indices(factor_index, K, floor),
    factor_index = factor_index,
    cluster = cluster_fit$cluster,
    description = "Top D distinct factors in each cluster center"
  )
}

initialize_threshold_recycle = function(cluster_fit, D,
                                         minimum_share = 0.10,
                                         floor = 0.05) {
  K = ncol(cluster_fit$centers)
  S = nrow(cluster_fit$centers)
  factor_index = matrix(0L, S, D)

  for (s in seq_len(S)) {
    center = cluster_fit$centers[s, ]
    ordered_factors = order(center, decreasing = TRUE)
    active = ordered_factors[center[ordered_factors] >= minimum_share]
    if (length(active) == 0L) active = ordered_factors[1L]
    if (length(active) > D) active = active[seq_len(D)]
    factor_index[s, ] = rep(active, length.out = D)
  }

  list(
    gamma_bar = gamma_from_factor_indices(factor_index, K, floor),
    factor_index = factor_index,
    cluster = cluster_fit$cluster,
    description = "Threshold active factors, then recycle them across slots"
  )
}

all_permutations = function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(j) {
    cbind(x[j], all_permutations(x[-j]))
  }))
}

best_one_to_one_assignment = function(score) {
  ## Rows are target MiSo slots; columns are row-wise Poisson-SuSiE slots.
  D = nrow(score)
  if (D == 1L) return(1L)

  if (D <= 8L) {
    candidates = all_permutations(seq_len(D))
    objective = apply(candidates, 1, function(candidate) {
      sum(score[cbind(seq_len(D), candidate)])
    })
    return(as.integer(candidates[which.max(objective), ]))
  }

  warning("Using greedy slot alignment because D > 8.")
  assignment = rep(NA_integer_, D)
  available = seq_len(D)
  target_order = order(apply(score, 1, max), decreasing = TRUE)
  for (target_slot in target_order) {
    chosen = available[which.max(score[target_slot, available])]
    assignment[target_slot] = chosen
    available = setdiff(available, chosen)
  }
  assignment
}

initialize_aligned_ps = function(ps_fit, distinct_initialization,
                                  floor = 0.05,
                                  epsilon = 1e-12) {
  ## The row-wise Poisson-SuSiE slots are exchangeable. For observation i,
  ## align them to the distinct seed factors of i's preliminary cluster.
  ## Then average their full gamma_bar[i,d,k] distributions instead of
  ## discarding this posterior information.
  N = dim(ps_fit$gamma_bar)[1]
  D = dim(ps_fit$gamma_bar)[2]
  K = dim(ps_fit$gamma_bar)[3]
  S = nrow(distinct_initialization$factor_index)

  expected_lambda = ps_fit$alpha / pmax(ps_fit$beta, epsilon)
  slot_loading = matrix(0, N, D)
  for (i in seq_len(N)) {
    for (d in seq_len(D)) {
      slot_loading[i, d] = sum(
        ps_fit$gamma_bar[i, d, ] * expected_lambda[i, d, ]
      )
    }
  }
  slot_fraction = slot_loading / pmax(rowSums(slot_loading), epsilon)

  posterior_sum = array(0, c(S, D, K))
  posterior_weight = matrix(0, S, D)
  slot_alignment = matrix(0L, N, D)

  for (i in seq_len(N)) {
    s = distinct_initialization$cluster[i]
    compatibility = matrix(0, D, D)

    for (target_slot in seq_len(D)) {
      seed_factor = distinct_initialization$factor_index[s, target_slot]
      for (local_slot in seq_len(D)) {
        compatibility[target_slot, local_slot] =
          ps_fit$gamma_bar[i, local_slot, seed_factor] *
          expected_lambda[i, local_slot, seed_factor]
      }
    }

    assignment = best_one_to_one_assignment(compatibility)
    slot_alignment[i, ] = assignment

    for (target_slot in seq_len(D)) {
      local_slot = assignment[target_slot]
      weight = slot_fraction[i, local_slot]
      posterior_sum[s, target_slot, ] =
        posterior_sum[s, target_slot, ] +
        weight * ps_fit$gamma_bar[i, local_slot, ]
      posterior_weight[s, target_slot] =
        posterior_weight[s, target_slot] + weight
    }
  }

  posterior_average = array(0, c(S, D, K))
  gamma_bar = array(0, c(S, D, K))
  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      if (posterior_weight[s, d] > epsilon) {
        posterior_average[s, d, ] =
          posterior_sum[s, d, ] / posterior_weight[s, d]
      } else {
        seed_factor = distinct_initialization$factor_index[s, d]
        posterior_average[s, d, seed_factor] = 1
      }
      posterior_average[s, d, ] =
        posterior_average[s, d, ] /
        sum(posterior_average[s, d, ])
      gamma_bar[s, d, ] =
        (1 - floor) * posterior_average[s, d, ] + floor / K
    }
  }

  list(
    gamma_bar = gamma_bar,
    factor_index = distinct_initialization$factor_index,
    cluster = distinct_initialization$cluster,
    posterior_average = posterior_average,
    posterior_weight = posterior_weight,
    slot_alignment = slot_alignment,
    description = paste(
      "Align row-wise Poisson-SuSiE slots, then average their",
      "factor-selection posteriors"
    )
  )
}

miso_initializations = list()

if (D_fit <= K_fit) {
  distinct_initialization = initialize_distinct(
    preliminary_clusters, D_fit, gamma_floor
  )
  miso_initializations$distinct = distinct_initialization
  miso_initializations$aligned_ps = initialize_aligned_ps(
    ps_fit, distinct_initialization, gamma_floor
  )
} else {
  warning("Skipping distinct and aligned_ps because D_fit > K_fit.")
}

miso_initializations$threshold_recycle = initialize_threshold_recycle(
  preliminary_clusters,
  D = D_fit,
  minimum_share = minimum_factor_share,
  floor = gamma_floor
)

if (!initialization_to_inspect %in% names(miso_initializations)) {
  stop(initialization_to_inspect,
       " is unavailable for these fitted dimensions. Use D_fit <= K_fit, ",
       "or choose threshold_recycle.")
}

cat("Available initializations:",
    paste(names(miso_initializations), collapse = ", "), "\n\n")


## 8. A readable, full MiSo variational algorithm --------------------------

gamma_kl = function(gamma_bar, epsilon = 1e-12) {
  ## KL[q(gamma_sd) || Uniform(1,...,K)].
  K = dim(gamma_bar)[3]
  sum(gamma_bar * (log(pmax(gamma_bar, epsilon)) + log(K)))
}

gamma_distribution_kl = function(alpha, beta, alpha0, beta0) {
  ## KL[Gamma(alpha,beta) || Gamma(alpha0,beta0)], using rate parameters.
  (alpha - alpha0) * digamma(alpha) +
    (log(beta) - log(beta0)) * alpha0 -
    (lgamma(alpha) - lgamma(alpha0)) -
    (beta - beta0) * alpha / beta
}

update_miso_local_variables = function(Y, F, gamma_bar, alpha0, beta0,
                                        alpha = NULL, n_inner = 5L,
                                        epsilon = 1e-12) {
  ## For every observation i and candidate submanifold s, update
  ## q(lambda_isd) and xi_isdm while holding global parameters fixed.
  N = nrow(Y)
  M = ncol(Y)
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  log_F = log(pmax(F, epsilon))
  factor_exposure = rowSums(F)

  beta = matrix(0, S, D)
  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      beta[s, d] = beta0[s, d] +
        sum(gamma_bar[s, d, ] * factor_exposure)
    }
  }

  if (is.null(alpha)) {
    alpha = array(0, c(N, S, D))
    for (s in seq_len(S)) {
      for (d in seq_len(D)) {
        alpha[, s, d] = alpha0[s, d] + rowSums(Y) / D
      }
    }
  }

  xi = array(1 / D, c(N, S, D, M))

  for (inner in seq_len(n_inner)) {
    for (s in seq_len(S)) {
      expected_log_factor =
        matrix(gamma_bar[s, , ], D, nrow(F)) %*% log_F

      for (i in seq_len(N)) {
        expected_log_lambda =
          digamma(alpha[i, s, ]) - log(beta[s, ])
        log_xi = expected_log_factor +
          matrix(expected_log_lambda, D, M)
        log_xi = sweep(log_xi, 2, apply(log_xi, 2, max), "-")
        xi_i = exp(log_xi)
        xi_i = sweep(xi_i, 2, colSums(xi_i), "/")
        xi[i, s, , ] = xi_i

        for (d in seq_len(D)) {
          alpha[i, s, d] = alpha0[s, d] + sum(xi_i[d, ] * Y[i, ])
        }
      }
    }
  }

  ## The component-specific variational lower bound drives q(z_i = s).
  component_lower_bound = matrix(0, N, S)
  for (s in seq_len(S)) {
    expected_log_factor =
      matrix(gamma_bar[s, , ], D, nrow(F)) %*% log_F
    expected_factor_exposure = numeric(D)
    for (d in seq_len(D)) {
      expected_factor_exposure[d] =
        sum(gamma_bar[s, d, ] * factor_exposure)
    }

    for (i in seq_len(N)) {
      for (d in seq_len(D)) {
        xi_id = xi[i, s, d, ]
        allocated_counts = xi_id * Y[i, ]
        expected_log_lambda = digamma(alpha[i, s, d]) - log(beta[s, d])
        expected_lambda = alpha[i, s, d] / beta[s, d]

        component_lower_bound[i, s] = component_lower_bound[i, s] +
          sum(allocated_counts *
                (expected_log_lambda + expected_log_factor[d, ] -
                   log(pmax(xi_id, epsilon)))) -
          expected_lambda * expected_factor_exposure[d] -
          gamma_distribution_kl(
            alpha[i, s, d], beta[s, d], alpha0[s, d], beta0[s, d]
          )
      }
    }
  }

  list(alpha = alpha, beta = beta, xi = xi,
       component_lower_bound = component_lower_bound)
}

accumulate_miso_counts = function(Y, xi, omega) {
  ## C[s,d,m] is the responsibility-weighted count allocated to slot d.
  N = nrow(Y)
  M = ncol(Y)
  S = ncol(omega)
  D = dim(xi)[3]
  C = array(0, c(S, D, M))

  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      for (i in seq_len(N)) {
        C[s, d, ] = C[s, d, ] +
          omega[i, s] * xi[i, s, d, ] * Y[i, ]
      }
    }
  }
  C
}

update_miso_gamma = function(F, gamma_bar, C, alpha, beta, omega,
                              step_size = 0.5, epsilon = 1e-12) {
  ## Coordinate update for q(gamma_sd), with a uniform factor prior.
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  K = nrow(F)
  log_F = log(pmax(F, epsilon))
  factor_exposure = rowSums(F)
  updated = gamma_bar

  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      expected_lambda_sum = sum(
        omega[, s] * alpha[, s, d] / beta[s, d]
      )
      log_probability = as.vector(log_F %*% C[s, d, ]) -
        expected_lambda_sum * factor_exposure - log(K)
      cavi_probability = softmax(log_probability)
      updated[s, d, ] =
        (1 - step_size) * gamma_bar[s, d, ] +
        step_size * cavi_probability
      updated[s, d, ] = updated[s, d, ] / sum(updated[s, d, ])
    }
  }
  updated
}

update_miso_loading_priors = function(alpha, beta, omega, alpha0, beta0,
                                       epsilon = 1e-12) {
  ## Responsibility-weighted empirical-Bayes moment updates.
  S = ncol(omega)
  D = ncol(alpha0)

  for (s in seq_len(S)) {
    total_weight = sum(omega[, s])
    if (total_weight <= epsilon) next

    for (d in seq_len(D)) {
      expected_lambda = sum(
        omega[, s] * alpha[, s, d] / beta[s, d]
      ) / total_weight
      expected_log_lambda = sum(
        omega[, s] *
          (digamma(alpha[, s, d]) - log(beta[s, d]))
      ) / total_weight
      alpha0[s, d] = gamma_shape_from_moments(
        expected_lambda, expected_log_lambda
      )
      beta0[s, d] = alpha0[s, d] / pmax(expected_lambda, epsilon)
    }
  }
  list(alpha0 = alpha0, beta0 = beta0)
}

update_miso_F = function(F, gamma_bar, C, step_size = 0.25,
                          epsilon = 1e-12) {
  ## Optional dictionary update. The expected allocated counts are averaged
  ## over each slot's factor-selection posterior, then rows are normalized.
  K = nrow(F)
  M = ncol(F)
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  expected_counts = matrix(epsilon, K, M)

  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      expected_counts = expected_counts + gamma_bar[s, d, ] %o% C[s, d, ]
    }
  }
  cavi_F = normalize_rows(expected_counts)
  normalize_rows((1 - step_size) * F + step_size * cavi_F)
}

fit_miso_simple = function(Y, F, initialization,
                            n_iterations = 80L, n_inner = 5L,
                            update_F = FALSE, update_priors = TRUE,
                            gamma_step_size = 0.5,
                            epsilon = 1e-12) {
  ## Model:
  ##   z_i ~ Categorical(pi)
  ##   gamma_sd ~ Categorical(1/K)
  ##   lambda_isd ~ Gamma(alpha0_sd, beta0_sd)
  ##   Y_im | z_i=s ~ Poisson(sum_d lambda_isd F[gamma_sd,m]).
  ##
  ## The variational algorithm alternates local loading/allocation updates,
  ## submanifold responsibilities, and global gamma/pi/prior updates.
  N = nrow(Y)
  S = dim(initialization$gamma_bar)[1]
  D = dim(initialization$gamma_bar)[2]
  gamma_bar = initialization$gamma_bar

  omega = matrix((1 - 0.90) / S, N, S)
  for (i in seq_len(N)) {
    omega[i, initialization$cluster[i]] = 0.90 + (1 - 0.90) / S
  }
  omega = omega / rowSums(omega)
  pi = colMeans(omega)

  cluster_mean_count = numeric(S)
  for (s in seq_len(S)) {
    cluster_mean_count[s] = mean(rowSums(Y)[initialization$cluster == s])
  }
  alpha0 = matrix(2, S, D)
  beta0 = matrix(0, S, D)
  for (s in seq_len(S)) {
    beta0[s, ] = 2 / pmax(cluster_mean_count[s] / D, 1)
  }

  alpha = NULL
  objective = numeric(n_iterations)
  gamma_history = vector("list", n_iterations)

  for (iteration in seq_len(n_iterations)) {
    local = update_miso_local_variables(
      Y, F, gamma_bar, alpha0, beta0,
      alpha = alpha, n_inner = n_inner, epsilon = epsilon
    )
    alpha = local$alpha

    log_responsibility = sweep(
      local$component_lower_bound, 2, log(pmax(pi, epsilon)), "+"
    )
    omega = softmax_rows(log_responsibility)
    pi = colMeans(omega)

    C = accumulate_miso_counts(Y, local$xi, omega)
    gamma_bar = update_miso_gamma(
      F, gamma_bar, C, alpha, local$beta, omega,
      step_size = gamma_step_size, epsilon = epsilon
    )

    if (update_priors) {
      prior_update = update_miso_loading_priors(
        alpha, local$beta, omega, alpha0, beta0, epsilon
      )
      alpha0 = prior_update$alpha0
      beta0 = prior_update$beta0
    }

    if (update_F) {
      F = update_miso_F(F, gamma_bar, C, epsilon = epsilon)
    }

    objective[iteration] = sum(apply(log_responsibility, 1, log_sum_exp)) -
      gamma_kl(gamma_bar, epsilon)
    gamma_history[[iteration]] = gamma_bar

    if (iteration > 5L) {
      relative_change = abs(objective[iteration] - objective[iteration - 1L]) /
        (1 + abs(objective[iteration - 1L]))
      if (relative_change < 1e-7) {
        objective = objective[seq_len(iteration)]
        gamma_history = gamma_history[seq_len(iteration)]
        break
      }
    }
  }

  ## Recompute local variables so every returned object agrees with final
  ## F, gamma_bar, alpha0, and beta0.
  local = update_miso_local_variables(
    Y, F, gamma_bar, alpha0, beta0,
    alpha = alpha, n_inner = n_inner, epsilon = epsilon
  )
  log_responsibility = sweep(
    local$component_lower_bound, 2, log(pmax(pi, epsilon)), "+"
  )
  omega = softmax_rows(log_responsibility)
  pi = colMeans(omega)

  expected_factor_loading = matrix(0, N, nrow(F))
  for (i in seq_len(N)) {
    for (s in seq_len(S)) {
      for (d in seq_len(D)) {
        expected_factor_loading[i, ] = expected_factor_loading[i, ] +
          omega[i, s] * local$alpha[i, s, d] / local$beta[s, d] *
          gamma_bar[s, d, ]
      }
    }
  }
  fitted_mean = expected_factor_loading %*% F

  list(
    F = F,
    gamma_bar = gamma_bar,
    gamma_bar_initial = initialization$gamma_bar,
    omega = omega,
    pi = pi,
    alpha = local$alpha,
    beta = local$beta,
    alpha0 = alpha0,
    beta0 = beta0,
    xi = local$xi,
    expected_factor_loading = expected_factor_loading,
    fitted_mean = fitted_mean,
    objective = objective,
    gamma_history = gamma_history,
    initialization = initialization
  )
}


## 9. Fit MiSo from the chosen initialization(s) ----------------------------

initializations_to_fit = if (run_all_initializations) {
  names(miso_initializations)
} else {
  initialization_to_inspect
}

miso_fits = list()
for (initialization_name in initializations_to_fit) {
  cat("Fitting MiSo from", initialization_name, "initialization ...\n")
  miso_fits[[initialization_name]] = fit_miso_simple(
    Y = Y,
    F = nmf_fit$F,
    initialization = miso_initializations[[initialization_name]],
    n_iterations = miso_iterations,
    n_inner = miso_inner_iterations,
    update_F = update_F_in_miso,
    update_priors = update_miso_priors
  )
}

## This alias makes the selected fit convenient to inspect interactively.
miso_fit = miso_fits[[initialization_to_inspect]]


## 10. Compact comparisons and optional plots ------------------------------

summarize_gamma = function(gamma_bar) {
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  answer = data.frame(
    submanifold = integer(S * D),
    slot = integer(S * D),
    top_factor = integer(S * D),
    top_probability = numeric(S * D),
    entropy = numeric(S * D)
  )
  row = 0L
  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      row = row + 1L
      probability = gamma_bar[s, d, ]
      answer[row, ] = c(
        s, d, which.max(probability), max(probability),
        -sum(probability * log(pmax(probability, 1e-12)))
      )
    }
  }
  answer
}

initialization_summaries = lapply(
  miso_initializations,
  function(x) summarize_gamma(x$gamma_bar)
)

fit_comparison = do.call(rbind, lapply(names(miso_fits), function(name) {
  fit = miso_fits[[name]]
  initial_top = summarize_gamma(fit$gamma_bar_initial)$top_factor
  final_top = summarize_gamma(fit$gamma_bar)$top_factor
  data.frame(
    initialization = name,
    iterations = length(fit$objective),
    final_objective = tail(fit$objective, 1),
    deviance_per_entry = poisson_deviance_per_entry(Y, fit$fitted_mean),
    mean_assignment_probability = mean(apply(fit$omega, 1, max)),
    initial_distinct_factors = length(unique(initial_top)),
    final_distinct_factors = length(unique(final_top))
  )
}))
rownames(fit_comparison) = NULL

cat("\nMiSo fit comparison:\n")
print(fit_comparison, digits = 4, row.names = FALSE)
cat("\nSelected fit mixture weights:\n")
print(miso_fit$pi, digits = 4)
cat("\nSelected fit initial gamma summary:\n")
print(summarize_gamma(miso_fit$gamma_bar_initial), digits = 4,
      row.names = FALSE)
cat("\nSelected fit final gamma summary:\n")
print(summarize_gamma(miso_fit$gamma_bar), digits = 4,
      row.names = FALSE)

## These plots appear when the script is sourced interactively in RStudio.
if (interactive()) {
  local({
    old_graphics_parameters = par(no.readonly = TRUE)
    on.exit(par(old_graphics_parameters), add = TRUE)
    par(mfrow = c(1, 2))

    plot(nmf_fit$objective, type = "l", lwd = 2,
         xlab = "NMF iteration", ylab = "Poisson objective",
         main = "Vanilla NMF")

    objective_range = range(unlist(lapply(miso_fits, function(x) x$objective)))
    if (diff(objective_range) == 0) objective_range = objective_range + c(-1, 1)
    plot(NA, xlim = c(1, max(vapply(
      miso_fits, function(x) length(x$objective), integer(1)
    ))),
    ylim = objective_range, xlab = "MiSo iteration",
    ylab = "Variational objective", main = "Initialization comparison")
    colors = seq_along(miso_fits)
    for (j in seq_along(miso_fits)) {
      lines(miso_fits[[j]]$objective, col = colors[j], lwd = 2)
    }
    legend("bottomright", legend = names(miso_fits), col = colors,
           lty = 1, lwd = 2, bty = "n")
  })
}

cat("\nUseful RStudio objects:\n")
cat("  View(Y)\n")
cat("  View(nmf_fit$F)\n")
cat("  View(ps_fit$gamma_bar[1, , ])\n")
cat("  initialization_summaries\n")
cat("  fit_comparison\n")
cat("  View(miso_fit$gamma_bar)\n")
cat("  View(miso_fit$omega)\n")

miso_fit$alpha0
miso_fit$beta0

## Gamma prior means and standard deviations.
miso_fit$alpha0 / miso_fit$beta0
sqrt(miso_fit$alpha0) / miso_fit$beta0

for (s in seq_len(S_fit)) {
  cat("\nSubmanifold", s, "\n")
  print(round(miso_fit$gamma_bar[s, , ], 4))
}

gamma_summary = summarize_gamma(miso_fit$gamma_bar)
gamma_summary$normalized_entropy =
  gamma_summary$entropy / log(K_fit)
gamma_summary$effective_number_of_factors =
  exp(gamma_summary$entropy)

print(gamma_summary, digits = 4, row.names = FALSE)

print(miso_fit$pi, digits = 4)
