## This analysis shows that learning K is actually harder than I think.
## MiSo does not automatically learns K very well.
## The posterior is certain for K = 2 instead of being diffused.
## K true = 1 here.

## MiSo single-initialization walkthrough -----------------------------------
##
## This script is deliberately self-contained. It uses only base R and keeps
## all important intermediate objects in the global environment so that they
## are easy to inspect in RStudio. It constructs and fits exactly one MiSo
## initialization so the algorithm can be studied without repeated fits.
##
## Suggested workflow:
##   1. Change the settings below.
##   2. Run one section at a time with Cmd/Ctrl + Shift + Enter.
##   3. Inspect objects such as Y, nmf_fit, ps_fit, miso_initialization,
##      and miso_fit with View() or str().


## 1. User settings ---------------------------------------------------------

## These are fitted dimensions. The simulated truth remains K = S = D = 1.
K_fit = 2L
S_fit = 1L
D_fit = 2L

## Exactly one initialization is constructed and fitted.
## Choices are "distinct", "threshold_recycle", and "aligned_ps".
initialization_to_use = "distinct"

## Simulation and fitting controls. These defaults are intentionally small.
N = 1000L
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
update_F_in_miso = TRUE

## Optional empirical-Bayes updates for the Gamma loading priors.
update_poisson_susie_priors = FALSE
update_miso_priors = TRUE

stopifnot(K_fit >= 1L, S_fit >= 1L, D_fit >= 1L)
stopifnot(K_fit <= N, S_fit <= N)
allowed_initializations = c("distinct", "threshold_recycle", "aligned_ps")
stopifnot(initialization_to_use %in% allowed_initializations)


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

## Spike-and-slab parameters for the single true factor. Change these values
## here to investigate different levels and forms of sparsity.
##
## `F_true_zero_probability` controls the expected fraction of exact zeros.
## `F_true_gamma_shape` controls heterogeneity among nonzero entries: values
## below one produce many small values and a few large values, whereas larger
## values make the nonzero entries more similar.
## `F_true_gamma_rate` controls their scale before normalization. Because
## F_true is subsequently normalized to sum to one, a common rate has no
## effect on its relative shape.
F_true_zero_probability = 0.50
F_true_gamma_shape = 2
F_true_gamma_rate = 1

stopifnot(
  F_true_zero_probability >= 0,
  F_true_zero_probability < 1,
  F_true_gamma_shape > 0,
  F_true_gamma_rate > 0
)

## Draw the spike indicators. In the extremely unlikely event that every
## feature is zero, activate one randomly chosen feature so normalization is
## well-defined.
F_true_is_nonzero =
  rbinom(M, size = 1, prob = 1 - F_true_zero_probability) == 1
if (!any(F_true_is_nonzero)) {
  F_true_is_nonzero[sample.int(M, size = 1L)] = TRUE
}

F_true_unnormalized = numeric(M)
F_true_unnormalized[F_true_is_nonzero] = rgamma(
  sum(F_true_is_nonzero),
  shape = F_true_gamma_shape,
  rate = F_true_gamma_rate
)

## Normalize to unit L1 length while preserving the point masses at zero.
F_true = matrix(
  F_true_unnormalized / sum(F_true_unnormalized),
  nrow = 1L
)

## With one factor and one loading, E[sum_m Y_im | L_i] = L_i.
L_true = matrix(rgamma(N, shape = 6, rate = 0.10), ncol = 1L)
true_mean = L_true %*% F_true
Y = matrix(rpois(N * M, lambda = as.vector(true_mean)), N, M)

cat("Simulated Y with dimensions", N, "x", M, "\n")
cat("True dimensions: K = S = D = 1\n")
cat("True factor nonzero entries:", sum(F_true_is_nonzero), "of", M, "\n")
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


## The fitted loading total is identifiable because every fitted factor sums
## to one. Individual columns of L are not identifiable when K is overfitted.
nmf_loading_check = data.frame(
  true_loading = L_true[, 1],
  observed_row_total = rowSums(Y),
  fitted_loading_total = rowSums(nmf_fit$L)
)

## Compare each fitted factor with the true factor using a relative L2 distance.
nmf_factor_relative_L2 = vapply(seq_len(K_fit), function(k) {
  sqrt(
    sum((nmf_fit$F[k, ] - F_true[1, ])^2) /
      sum(F_true[1, ]^2)
  )
}, numeric(1))

## When K >= 2, inspect whether the first two fitted factors are genuinely
## similar to one another.
if (K_fit >= 2L) {
  nmf_factor_1_2_cosine =
    sum(nmf_fit$F[1, ] * nmf_fit$F[2, ]) /
    sqrt(sum(nmf_fit$F[1, ]^2) * sum(nmf_fit$F[2, ]^2))
}


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


## 7. Initialization choices ------------------------------------------------
## The functions are kept together for study, but only `initialization_to_use`
## is constructed below and only that initialization is fitted.

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

## Construct only the initialization selected at the top of the script.
if (initialization_to_use %in% c("distinct", "aligned_ps") &&
    D_fit > K_fit) {
  stop(initialization_to_use,
       " requires D_fit <= K_fit. Choose threshold_recycle when D_fit > K_fit.")
}

if (initialization_to_use == "distinct") {
  miso_initialization = initialize_distinct(
    preliminary_clusters,
    D = D_fit,
    floor = gamma_floor
  )
} else if (initialization_to_use == "threshold_recycle") {
  miso_initialization = initialize_threshold_recycle(
    preliminary_clusters,
    D = D_fit,
    minimum_share = minimum_factor_share,
    floor = gamma_floor
  )
} else if (initialization_to_use == "aligned_ps") {
  ## Alignment needs distinct factors only as slot labels. Only the resulting
  ## aligned Poisson-SuSiE initialization will be fitted.
  distinct_seed = initialize_distinct(
    preliminary_clusters,
    D = D_fit,
    floor = gamma_floor
  )
  miso_initialization = initialize_aligned_ps(
    ps_fit,
    distinct_initialization = distinct_seed,
    floor = gamma_floor
  )
}

cat("Using", initialization_to_use, "initialization:\n")
cat(" ", miso_initialization$description, "\n\n")


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
  final_elbo = sum(apply(log_responsibility, 1, log_sum_exp)) -
    gamma_kl(gamma_bar, epsilon)

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
    final_elbo = final_elbo,
    gamma_history = gamma_history,
    initialization = initialization
  )
}


## 9. Fit MiSo once ---------------------------------------------------------

cat("Fitting MiSo from", initialization_to_use, "initialization ...\n")

miso_fit = fit_miso_simple(
  Y = Y,
  F = nmf_fit$F,
  initialization = miso_initialization,
  n_iterations = miso_iterations,
  n_inner = miso_inner_iterations,
  update_F = update_F_in_miso,
  update_priors = update_miso_priors
)


## 10. Fit a matched K = 1, S = 1, D = 1 reference -------------------------

cat("Fitting the K = 1, S = 1, D = 1 reference ...\n")

nmf_k1_fit = fit_poisson_nmf(
  Y = Y,
  K = 1L,
  n_iterations = nmf_iterations,
  seed = nmf_seed
)

miso_k1_d1_initialization = list(
  gamma_bar = array(1, dim = c(1L, 1L, 1L)),
  factor_index = matrix(1L, nrow = 1L, ncol = 1L),
  cluster = rep(1L, N),
  description = "The only slot selects the only rank-one NMF factor"
)

miso_k1_d1_fit = fit_miso_simple(
  Y = Y,
  F = nmf_k1_fit$F,
  initialization = miso_k1_d1_initialization,
  n_iterations = miso_iterations,
  n_inner = miso_inner_iterations,
  update_F = update_F_in_miso,
  update_priors = update_miso_priors
)

elbo_comparison = data.frame(
  model = c(
    "MiSo K=1, S=1, D=1 from rank-one NMF",
    paste0(
      "MiSo K=", K_fit, ", S=", S_fit, ", D=", D_fit,
      " from rank-", K_fit, " NMF (", initialization_to_use, ")"
    )
  ),
  final_elbo = c(miso_k1_d1_fit$final_elbo, miso_fit$final_elbo),
  deviance_per_entry = c(
    poisson_deviance_per_entry(Y, miso_k1_d1_fit$fitted_mean),
    poisson_deviance_per_entry(Y, miso_fit$fitted_mean)
  ),
  iterations = c(
    length(miso_k1_d1_fit$objective),
    length(miso_fit$objective)
  )
)

elbo_difference_k2_minus_k1 =
  miso_fit$final_elbo - miso_k1_d1_fit$final_elbo


## 11. Focused diagnostics and optional plots -------------------------------

summarize_gamma = function(gamma_bar) {
  S = dim(gamma_bar)[1]
  D = dim(gamma_bar)[2]
  K = dim(gamma_bar)[3]
  answer = data.frame(
    submanifold = integer(S * D),
    slot = integer(S * D),
    top_factor = integer(S * D),
    top_probability = numeric(S * D),
    normalized_entropy = numeric(S * D),
    effective_number_of_factors = numeric(S * D)
  )

  row = 0L
  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      row = row + 1L
      probability = gamma_bar[s, d, ]
      entropy = -sum(probability * log(pmax(probability, 1e-12)))
      normalized_entropy = if (K == 1L) 0 else entropy / log(K)
      answer[row, ] = c(
        s,
        d,
        which.max(probability),
        max(probability),
        normalized_entropy,
        exp(entropy)
      )
    }
  }
  answer
}

initial_gamma_summary = summarize_gamma(miso_fit$gamma_bar_initial)
final_gamma_summary = summarize_gamma(miso_fit$gamma_bar)

## alpha0 and beta0 are learned Gamma-prior parameters. beta0 is a rate.
learned_prior_mean = miso_fit$alpha0 / miso_fit$beta0
learned_prior_sd = sqrt(miso_fit$alpha0) / miso_fit$beta0

slot_summary = do.call(rbind, lapply(seq_len(S_fit), function(s) {
  weights = miso_fit$omega[, s]
  weights = weights / sum(weights)

  do.call(rbind, lapply(seq_len(D_fit), function(d) {
    probability = miso_fit$gamma_bar[s, d, ]
    posterior_loading =
      miso_fit$alpha[, s, d] / miso_fit$beta[s, d]

    data.frame(
      submanifold = s,
      slot = d,
      prior_alpha = miso_fit$alpha0[s, d],
      prior_beta = miso_fit$beta0[s, d],
      prior_mean = learned_prior_mean[s, d],
      posterior_mean = sum(weights * posterior_loading),
      top_factor = which.max(probability),
      top_probability = max(probability),
      normalized_entropy =
        if (K_fit == 1L) 0 else
          -sum(probability * log(pmax(probability, 1e-12))) / log(K_fit)
    )
  }))
}))

cat("\nInitialization:", initialization_to_use, "\n")
cat("MiSo iterations:", length(miso_fit$objective), "\n")
cat("Poisson deviance per entry:",
    signif(poisson_deviance_per_entry(Y, miso_fit$fitted_mean), 4), "\n")

cat("\nELBO comparison with the explicit K = 1, S = 1, D = 1 fit:\n")
print(elbo_comparison, digits = 8, row.names = FALSE)
cat("K =", K_fit, ", D =", D_fit, "minus K = 1, D = 1 ELBO:",
    signif(elbo_difference_k2_minus_k1, 10), "\n")
if (elbo_difference_k2_minus_k1 > 0) {
  cat("The fitted K =", K_fit, ", D =", D_fit,
      "solution has the larger ELBO.\n")
} else if (elbo_difference_k2_minus_k1 < 0) {
  cat("The K = 1, D = 1 solution has the larger ELBO.\n")
} else {
  cat("The two solutions have equal ELBOs to numerical precision.\n")
}

cat("\nLearned submanifold proportions pi:\n")
print(
  setNames(miso_fit$pi, paste0("submanifold_", seq_len(S_fit))),
  digits = 4
)

cat("\nLearned Gamma prior shapes alpha0:\n")
print(miso_fit$alpha0, digits = 4)

cat("\nLearned Gamma prior rates beta0:\n")
print(miso_fit$beta0, digits = 4)

cat("\nLearned Gamma prior means alpha0 / beta0:\n")
print(learned_prior_mean, digits = 4)

cat("\nInitial gamma summary:\n")
print(initial_gamma_summary, digits = 4, row.names = FALSE)

cat("\nFinal gamma summary:\n")
print(final_gamma_summary, digits = 4, row.names = FALSE)

cat("\nSlot activity summary:\n")
print(slot_summary, digits = 4, row.names = FALSE)

for (s in seq_len(S_fit)) {
  cat("\nFull gamma_bar for submanifold", s, "\n")
  print(miso_fit$gamma_bar[s, , ], digits = 6)
}

if (interactive()) {
  old_graphics = par(no.readonly = TRUE)
  on.exit(par(old_graphics), add = TRUE)
  par(mfrow = c(1, 2))

  plot(
    nmf_fit$objective,
    type = "l",
    lwd = 2,
    xlab = "Iteration",
    ylab = "Objective",
    main = "Poisson NMF"
  )
  plot(
    miso_fit$objective,
    type = "l",
    lwd = 2,
    xlab = "Iteration",
    ylab = "Objective",
    main = paste("MiSo:", initialization_to_use)
  )
}

# cat("\nUseful objects to inspect in RStudio:\n")
# cat("  View(nmf_loading_check)\n")
# cat("  nmf_factor_relative_L2\n")
# cat("  nmf_factor_1_2_cosine  # available when K_fit >= 2\n")
# cat("  View(ps_fit$gamma_bar[1, , ])\n")
# cat("  View(miso_initialization$gamma_bar)\n")
# cat("  View(miso_fit$gamma_bar)\n")
# cat("  View(miso_fit$omega)\n")
# cat("  elbo_comparison\n")
# cat("  miso_k1_d1_fit\n")
# cat("  learned_prior_mean\n")
# cat("  slot_summary\n")
