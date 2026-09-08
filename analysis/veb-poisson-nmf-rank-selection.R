## VEB Poisson NMF: can it remove an overfitted factor? ----------------------
##
## This script is deliberately self-contained and uses only base R. It studies
## the smallest rank-selection experiment:
##
##   true K = 1, fitted K = 2.
##
## The fitted loading columns have independent factor-specific Gamma priors,
##
##   L_ik ~ Gamma(alpha0_k, beta0_k),
##
## where beta0_k is a rate. Mean-field variational inference estimates q(L),
## and variational empirical Bayes learns alpha0_k and beta0_k. There is no
## hard pruning inside the algorithm. We ask whether one beta0_k grows and its
## prior mean alpha0_k / beta0_k approaches zero.


## 1. Settings ---------------------------------------------------------------

K_true = 1L
K_fit = 2L

## The original experiment used N = 300, M = 200, and a mean row total of 60.
## Here both matrix dimensions are doubled. Because each row of F sums to one,
## we also double the mean row total so that increasing M does not halve the
## expected count per nonzero feature.
N = 600L
M = 400L
simulation_seed = 1L
nmf_seed = 2L

## Parameters of the true spike-and-slab factor.
F_true_zero_probability = 0.80
F_true_gamma_shape = 2
F_true_gamma_rate = 1

## Parameters of the true loading distribution. Its mean is 120.
L_true_gamma_shape = 6
L_true_gamma_mean = 120
L_true_gamma_rate = L_true_gamma_shape / L_true_gamma_mean

nmf_iterations = 300L
veb_iterations = 300L
veb_inner_iterations = 2L

## Update F during VEB to study the full VEB-NMF algorithm. Set this to FALSE
## to isolate prior learning while holding the vanilla-NMF dictionary fixed.
update_F_in_veb = TRUE

## Keep the ordinary vanilla-NMF initialization as a reference, but do not use
## it to initialize the random starts.
include_vanilla_nmf_reference_start = TRUE

## Independent random VEB starts. A concentration below one generates a useful
## mixture of balanced and naturally asymmetric global loading shares without
## explicitly choosing which factor should be removed. Set it to 1 for a
## uniform draw on the simplex; such starts spend much less time near a
## rank-deficient boundary.
number_random_starts = 8L
random_start_seed = 100L
random_global_share_concentration = 0.10
random_row_share_concentration = 30
random_factor_pseudocount = 1e-3

## Use "random_data_partition" for the full experiment. The alternative is a
## theoretical control with F_1 = F_2 = F_true; normally combine that control
## with update_F_in_veb = FALSE.
random_dictionary_mode = "random_data_partition"

## A factor is called effectively active only for reporting after convergence.
## This threshold is not used in any update.
minimum_relevance_share = 1e-3
minimum_allocated_count_share = 1e-3

## Numerical controls.
epsilon = 1e-12
minimum_gamma_shape = 1e-4
maximum_gamma_shape = 1e6
maximum_gamma_rate = 1e12
stop_early = FALSE
minimum_veb_iterations = 100L
relative_elbo_tolerance = 1e-8

stopifnot(K_true == 1L, K_fit == 2L)
stopifnot(N >= 1L, M >= 1L)
stopifnot(F_true_zero_probability >= 0, F_true_zero_probability < 1)
stopifnot(F_true_gamma_shape > 0, F_true_gamma_rate > 0)
stopifnot(L_true_gamma_shape > 0, L_true_gamma_rate > 0)
stopifnot(number_random_starts >= 1L)
stopifnot(random_global_share_concentration > 0)
stopifnot(random_row_share_concentration > 0)
stopifnot(random_factor_pseudocount >= 0)
stopifnot(random_dictionary_mode %in%
            c("random_data_partition", "exact_duplicate_true"))


## 2. Numerical helpers ------------------------------------------------------

normalize_rows_preserve_zeros = function(x, epsilon = 1e-12) {
  x = as.matrix(x)
  answer = x
  for (row in seq_len(nrow(x))) {
    total = sum(x[row, ])
    if (total <= epsilon) stop("Cannot normalize a row with zero total mass.")
    answer[row, ] = x[row, ] / total
  }
  answer
}

softmax_rows = function(log_weights) {
  log_weights = as.matrix(log_weights)
  largest = apply(log_weights, 1, max)
  weights = exp(sweep(log_weights, 1, largest, "-"))
  weights / rowSums(weights)
}

poisson_deviance_per_entry = function(Y, fitted_mean) {
  fitted_mean = pmax(fitted_mean, 1e-12)
  contribution = fitted_mean - Y
  positive = Y > 0
  contribution[positive] = contribution[positive] +
    Y[positive] * log(Y[positive] / fitted_mean[positive])
  2 * sum(contribution) / length(Y)
}

gamma_kl = function(alpha, beta, alpha0, beta0) {
  ## KL[Gamma(alpha,beta) || Gamma(alpha0,beta0)], using rate parameters.
  (alpha - alpha0) * digamma(alpha) +
    alpha0 * (log(beta) - log(beta0)) -
    lgamma(alpha) + lgamma(alpha0) -
    alpha + alpha * beta0 / beta
}

solve_gamma_shape = function(expected_value, expected_log_value,
                              initial_shape = 1,
                              minimum_shape = 1e-4,
                              maximum_shape = 1e6,
                              n_steps = 20L) {
  ## Solve
  ##   log(alpha0) - digamma(alpha0)
  ##     = log(E[lambda]) - E[log(lambda)].
  difference = log(pmax(expected_value, epsilon)) - expected_log_value
  difference = pmax(difference, 1e-12)

  shape = min(max(initial_shape, minimum_shape), maximum_shape)
  for (step in seq_len(n_steps)) {
    value = log(shape) - digamma(shape) - difference
    derivative = 1 / shape - trigamma(shape)
    proposal = shape - value / derivative
    if (!is.finite(proposal) || proposal <= 0) proposal = shape / 2
    shape = min(max(proposal, minimum_shape), maximum_shape)
  }
  shape
}


## 3. Simulate a true rank-one Gamma-Poisson factorization ------------------

set.seed(simulation_seed)

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
F_true = matrix(
  F_true_unnormalized / sum(F_true_unnormalized),
  nrow = K_true
)

L_true = matrix(
  rgamma(N, shape = L_true_gamma_shape, rate = L_true_gamma_rate),
  nrow = N,
  ncol = K_true
)
true_mean = L_true %*% F_true
Y = matrix(rpois(N * M, as.vector(true_mean)), N, M)

cat("True K:", K_true, "  Fitted K:", K_fit, "\n")
cat("Y dimensions:", N, "x", M, "\n")
cat("True factor nonzero entries:", sum(F_true_is_nonzero), "of", M, "\n")
cat("Mean observed row total:", signif(mean(rowSums(Y)), 5), "\n\n")


## 4. Vanilla Poisson NMF initialization ------------------------------------

fit_poisson_nmf = function(Y, K, n_iterations = 300L, seed = 1L,
                           epsilon = 1e-12) {
  set.seed(seed)
  N = nrow(Y)
  M = ncol(Y)

  L = matrix(rexp(N * K), N, K)
  F = matrix(rexp(K * M), K, M)
  F = normalize_rows_preserve_zeros(F, epsilon)
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
  Y,
  K = K_fit,
  n_iterations = nmf_iterations,
  seed = nmf_seed,
  epsilon = epsilon
)

cat("Vanilla-NMF deviance per entry:",
    signif(poisson_deviance_per_entry(Y, nmf_fit$fitted_mean), 5), "\n")
cat("Cosine similarity between its two fitted factors:",
    signif(
      sum(nmf_fit$F[1, ] * nmf_fit$F[2, ]) /
        sqrt(sum(nmf_fit$F[1, ]^2) * sum(nmf_fit$F[2, ]^2)),
      5
    ), "\n\n")


## 5. VEB Poisson NMF --------------------------------------------------------

make_random_veb_start = function(Y, K, seed,
                                  global_share_concentration = 0.1,
                                  row_share_concentration = 30,
                                  factor_pseudocount = 1e-3,
                                  dictionary_mode = "random_data_partition",
                                  F_true = NULL,
                                  epsilon = 1e-12) {
  set.seed(seed)
  N = nrow(Y)
  M = ncol(Y)

  ## First draw a global point on the K-simplex. With concentration below one,
  ## some starts arise naturally near a simplex boundary.
  global_share = rgamma(K, shape = global_share_concentration, rate = 1)
  global_share = global_share / sum(global_share)

  ## Then draw observation-specific loading shares around that global point.
  component_shape = pmax(row_share_concentration * global_share, 1e-6)
  unnormalized_row_share = matrix(
    rgamma(N * K, shape = rep(component_shape, each = N), rate = 1),
    N,
    K
  )
  unnormalized_row_share = pmax(unnormalized_row_share, epsilon)
  row_share = unnormalized_row_share / rowSums(unnormalized_row_share)
  loading_start = pmax(row_share * rowSums(Y), epsilon)

  if (dictionary_mode == "exact_duplicate_true") {
    if (is.null(F_true)) stop("F_true is required for the duplicate control.")
    F_initial = matrix(rep(F_true[1, ], each = K), K, M)
  } else {
    ## Randomly weighted empirical profiles give positive, data-informed factor
    ## rows without running NMF. The pseudocount protects empty features.
    F_initial = matrix(0, K, M)
    for (k in seq_len(K)) {
      weighted_count = colSums(
        Y * matrix(row_share[, k], N, M)
      ) + factor_pseudocount
      F_initial[k, ] = weighted_count / sum(weighted_count)
    }
  }

  list(
    F = F_initial,
    L = loading_start,
    global_share = global_share,
    realized_loading_share = colSums(loading_start) / sum(loading_start),
    seed = seed
  )
}

initialize_gamma_prior = function(loading_start) {
  K = ncol(loading_start)
  alpha0 = numeric(K)
  beta0 = numeric(K)

  for (k in seq_len(K)) {
    expected_value = mean(loading_start[, k])
    expected_log_value = mean(log(pmax(loading_start[, k], epsilon)))
    alpha0[k] = solve_gamma_shape(
      expected_value,
      expected_log_value,
      initial_shape = 1,
      minimum_shape = minimum_gamma_shape,
      maximum_shape = maximum_gamma_shape
    )
    beta0[k] = alpha0[k] / pmax(expected_value, epsilon)
  }

  list(alpha0 = alpha0, beta0 = beta0)
}

update_allocations = function(Y, F, alpha, beta, epsilon = 1e-12) {
  N = nrow(Y)
  M = ncol(Y)
  K = nrow(F)
  expected_log_loading = digamma(alpha) - log(beta)
  log_F = log(pmax(F, epsilon))
  xi = array(0, c(N, M, K))

  ## For each feature, normalize its K component scores for every observation.
  for (m in seq_len(M)) {
    log_score = expected_log_loading +
      matrix(log_F[, m], N, K, byrow = TRUE)
    xi[, m, ] = softmax_rows(log_score)
  }
  xi
}

update_loading_posterior = function(Y, F, xi, alpha0, beta0) {
  N = nrow(Y)
  K = nrow(F)
  allocated_count = matrix(0, N, K)

  for (k in seq_len(K)) {
    allocated_count[, k] = rowSums(Y * xi[, , k])
  }

  alpha = sweep(allocated_count, 2, alpha0, "+")
  beta = matrix(0, N, K)
  for (k in seq_len(K)) {
    beta[, k] = beta0[k] + sum(F[k, ])
  }

  list(alpha = alpha, beta = beta, allocated_count = allocated_count)
}

update_factor_dictionary = function(Y, xi, F_previous,
                                     epsilon = 1e-12) {
  K = dim(xi)[3]
  F = F_previous
  allocated_feature_count = matrix(0, K, ncol(Y))

  for (k in seq_len(K)) {
    allocated_feature_count[k, ] = colSums(Y * xi[, , k])
    if (sum(allocated_feature_count[k, ]) > epsilon) {
      F[k, ] = allocated_feature_count[k, ] /
        sum(allocated_feature_count[k, ])
    }
  }

  list(F = F, allocated_feature_count = allocated_feature_count)
}

update_gamma_priors = function(alpha, beta, alpha0_previous) {
  K = ncol(alpha)
  alpha0 = numeric(K)
  beta0 = numeric(K)

  for (k in seq_len(K)) {
    expected_value = mean(alpha[, k] / beta[, k])
    expected_log_value = mean(digamma(alpha[, k]) - log(beta[, k]))
    alpha0[k] = solve_gamma_shape(
      expected_value,
      expected_log_value,
      initial_shape = alpha0_previous[k],
      minimum_shape = minimum_gamma_shape,
      maximum_shape = maximum_gamma_shape
    )
    beta0[k] = min(
      alpha0[k] / pmax(expected_value, epsilon),
      maximum_gamma_rate
    )
  }

  list(alpha0 = alpha0, beta0 = beta0)
}

calculate_elbo = function(Y, F, xi, alpha, beta, alpha0, beta0,
                           epsilon = 1e-12) {
  N = nrow(Y)
  M = ncol(Y)
  K = nrow(F)
  expected_log_loading = digamma(alpha) - log(beta)
  expected_loading = alpha / beta
  answer = -sum(lgamma(Y + 1))

  for (k in seq_len(K)) {
    log_contribution =
      matrix(expected_log_loading[, k], N, M) +
      matrix(log(pmax(F[k, ], epsilon)), N, M, byrow = TRUE) -
      log(pmax(xi[, , k], epsilon))
    answer = answer + sum(Y * xi[, , k] * log_contribution)
    answer = answer - sum(expected_loading[, k]) * sum(F[k, ])
    answer = answer - sum(gamma_kl(
      alpha[, k], beta[, k], alpha0[k], beta0[k]
    ))
  }
  answer
}

fit_veb_poisson_nmf = function(Y, F_initial, loading_start,
                                n_iterations = 400L,
                                n_inner = 3L,
                                update_F = TRUE,
                                stop_early = FALSE,
                                minimum_iterations = 100L,
                                tolerance = 1e-8,
                                epsilon = 1e-12) {
  N = nrow(Y)
  K = nrow(F_initial)
  F = F_initial

  prior = initialize_gamma_prior(loading_start)
  alpha0 = prior$alpha0
  beta0 = prior$beta0

  ## Initialize q(L) around the positive vanilla-NMF loading start.
  alpha = matrix(alpha0, N, K, byrow = TRUE) + loading_start
  beta = matrix(beta0 + 1, N, K, byrow = TRUE)

  elbo = rep(NA_real_, n_iterations)
  alpha0_history = matrix(NA_real_, n_iterations, K)
  beta0_history = matrix(NA_real_, n_iterations, K)
  relevance_history = matrix(NA_real_, n_iterations, K)
  allocated_share_history = matrix(NA_real_, n_iterations, K)

  for (iteration in seq_len(n_iterations)) {
    for (inner in seq_len(n_inner)) {
      xi = update_allocations(Y, F, alpha, beta, epsilon)
      posterior = update_loading_posterior(Y, F, xi, alpha0, beta0)
      alpha = posterior$alpha
      beta = posterior$beta
    }

    if (update_F) {
      factor_update = update_factor_dictionary(Y, xi, F, epsilon)
      F = factor_update$F
    }

    prior = update_gamma_priors(alpha, beta, alpha0)
    alpha0 = prior$alpha0
    beta0 = prior$beta0

    elbo[iteration] = calculate_elbo(
      Y, F, xi, alpha, beta, alpha0, beta0, epsilon
    )
    alpha0_history[iteration, ] = alpha0
    beta0_history[iteration, ] = beta0
    relevance_history[iteration, ] = alpha0 / beta0
    allocated_share_history[iteration, ] =
      colSums(posterior$allocated_count) / sum(Y)

    if (stop_early && iteration >= minimum_iterations && iteration > 1L) {
      relative_change = abs(elbo[iteration] - elbo[iteration - 1L]) /
        (1 + abs(elbo[iteration - 1L]))
      if (relative_change < tolerance) {
        keep = seq_len(iteration)
        elbo = elbo[keep]
        alpha0_history = alpha0_history[keep, , drop = FALSE]
        beta0_history = beta0_history[keep, , drop = FALSE]
        relevance_history = relevance_history[keep, , drop = FALSE]
        allocated_share_history =
          allocated_share_history[keep, , drop = FALSE]
        break
      }
    }
  }

  ## Recompute local quantities so all returned objects use the final globals.
  for (inner in seq_len(n_inner)) {
    xi = update_allocations(Y, F, alpha, beta, epsilon)
    posterior = update_loading_posterior(Y, F, xi, alpha0, beta0)
    alpha = posterior$alpha
    beta = posterior$beta
  }
  fitted_mean = (alpha / beta) %*% F
  final_elbo = calculate_elbo(
    Y, F, xi, alpha, beta, alpha0, beta0, epsilon
  )
  final_allocated_share = colSums(posterior$allocated_count) / sum(Y)
  final_relevance = alpha0 / beta0
  relevance_share = final_relevance / sum(final_relevance)
  active = relevance_share >= minimum_relevance_share &
    final_allocated_share >= minimum_allocated_count_share

  list(
    F = F,
    alpha = alpha,
    beta = beta,
    alpha0 = alpha0,
    beta0 = beta0,
    relevance = final_relevance,
    relevance_share = relevance_share,
    allocated_share = final_allocated_share,
    active = active,
    effective_K = sum(active),
    xi = xi,
    fitted_mean = fitted_mean,
    elbo = elbo,
    final_elbo = final_elbo,
    alpha0_history = alpha0_history,
    beta0_history = beta0_history,
    relevance_history = relevance_history,
    allocated_share_history = allocated_share_history
  )
}


## 6. Fit the explicit K = 1 reference model --------------------------------

explicit_k1_start = list(
  F = matrix(colSums(Y) / sum(Y), nrow = 1L),
  L = matrix(rowSums(Y), ncol = 1L)
)

cat("Running the explicit K = 1 reference model ...\n")
explicit_k1_fit = fit_veb_poisson_nmf(
  Y,
  F_initial = explicit_k1_start$F,
  loading_start = explicit_k1_start$L,
  n_iterations = veb_iterations,
  n_inner = veb_inner_iterations,
  update_F = update_F_in_veb,
  stop_early = stop_early,
  minimum_iterations = minimum_veb_iterations,
  tolerance = relative_elbo_tolerance,
  epsilon = epsilon
)


## 7. Run K = 2 VEB from independent random starts --------------------------

veb_starts = list()

if (include_vanilla_nmf_reference_start) {
  veb_starts$vanilla_nmf_reference = list(
    F = nmf_fit$F,
    L = nmf_fit$L,
    global_share = colSums(nmf_fit$L) / sum(nmf_fit$L),
    realized_loading_share = colSums(nmf_fit$L) / sum(nmf_fit$L),
    seed = nmf_seed
  )
}

for (start_index in seq_len(number_random_starts)) {
  start_name = sprintf("random_%02d", start_index)
  veb_starts[[start_name]] = make_random_veb_start(
    Y,
    K = K_fit,
    seed = random_start_seed + start_index - 1L,
    global_share_concentration = random_global_share_concentration,
    row_share_concentration = random_row_share_concentration,
    factor_pseudocount = random_factor_pseudocount,
    dictionary_mode = random_dictionary_mode,
    F_true = F_true,
    epsilon = epsilon
  )
}

veb_fits = list()
for (start_name in names(veb_starts)) {
  cat("Running VEB from", start_name, "...\n")
  start = veb_starts[[start_name]]
  fit = fit_veb_poisson_nmf(
    Y,
    F_initial = start$F,
    loading_start = start$L,
    n_iterations = veb_iterations,
    n_inner = veb_inner_iterations,
    update_F = update_F_in_veb,
    stop_early = stop_early,
    minimum_iterations = minimum_veb_iterations,
    tolerance = relative_elbo_tolerance,
    epsilon = epsilon
  )
  fit$initial_global_share = start$global_share
  fit$initial_realized_loading_share = start$realized_loading_share
  fit$initialization_seed = start$seed
  veb_fits[[start_name]] = fit
}


## 8. Compare the stationary solutions --------------------------------------

fit_summary = do.call(rbind, lapply(names(veb_fits), function(start_name) {
  fit = veb_fits[[start_name]]
  data.frame(
    initialization = start_name,
    seed = fit$initialization_seed,
    initial_share_1 = fit$initial_realized_loading_share[1],
    initial_share_2 = fit$initial_realized_loading_share[2],
    iterations = length(fit$elbo),
    final_elbo = fit$final_elbo,
    deviance_per_entry = poisson_deviance_per_entry(Y, fit$fitted_mean),
    effective_K = fit$effective_K,
    alpha0_1 = fit$alpha0[1],
    alpha0_2 = fit$alpha0[2],
    beta0_1 = fit$beta0[1],
    beta0_2 = fit$beta0[2],
    prior_mean_1 = fit$relevance[1],
    prior_mean_2 = fit$relevance[2],
    prior_cv_1 = 1 / sqrt(fit$alpha0[1]),
    prior_cv_2 = 1 / sqrt(fit$alpha0[2]),
    allocated_share_1 = fit$allocated_share[1],
    allocated_share_2 = fit$allocated_share[2],
    factor_cosine = sum(fit$F[1, ] * fit$F[2, ]) /
      sqrt(sum(fit$F[1, ]^2) * sum(fit$F[2, ]^2))
  )
}))
rownames(fit_summary) = NULL

best_initialization = fit_summary$initialization[which.max(fit_summary$final_elbo)]
veb_fit = veb_fits[[best_initialization]]
random_start_rows = grepl("^random_", fit_summary$initialization)
number_random_starts_selecting_one = sum(
  fit_summary$effective_K[random_start_rows] == 1L
)
best_random_row = which(
  random_start_rows &
    fit_summary$final_elbo == max(fit_summary$final_elbo[random_start_rows])
)[1]
rank_one_rows = which(fit_summary$effective_K == 1L)
rank_two_rows = which(fit_summary$effective_K == 2L)

cat("\nVEB solutions:\n")
print(fit_summary, digits = 6, row.names = FALSE)

cat("\nRandom starts selecting effective K = 1:",
    number_random_starts_selecting_one, "of", number_random_starts, "\n")
cat("Highest-ELBO random start:",
    fit_summary$initialization[best_random_row], "\n")
cat("Effective K of highest-ELBO random start:",
    fit_summary$effective_K[best_random_row], "\n")

cat("\nExplicit K = 1 ELBO:",
    signif(explicit_k1_fit$final_elbo, 8), "\n")
cat("Explicit K = 1 deviance per entry:",
    signif(poisson_deviance_per_entry(Y, explicit_k1_fit$fitted_mean), 8),
    "\n")

if (length(rank_one_rows) > 0L && length(rank_two_rows) > 0L) {
  best_rank_one_elbo = max(fit_summary$final_elbo[rank_one_rows])
  best_rank_two_elbo = max(fit_summary$final_elbo[rank_two_rows])
  cat("Best effective-K=1 ELBO:", signif(best_rank_one_elbo, 8), "\n")
  cat("Best effective-K=2 ELBO:", signif(best_rank_two_elbo, 8), "\n")
  cat("Best K=2 minus best K=1 ELBO:",
      signif(best_rank_two_elbo - best_rank_one_elbo, 8), "\n")
  cat("ELBO difference per observed count:",
      signif((best_rank_two_elbo - best_rank_one_elbo) / sum(Y), 8), "\n")
  cat("Best K=2 minus explicit K=1 ELBO:",
      signif(best_rank_two_elbo - explicit_k1_fit$final_elbo, 8), "\n")
  cat("Best pruned K=2 minus explicit K=1 ELBO:",
      signif(best_rank_one_elbo - explicit_k1_fit$final_elbo, 8), "\n")
}

if (fit_summary$effective_K[best_random_row] == 1L) {
  cat("Conclusion: random multistart plus maximum ELBO selected K = 1.\n")
} else if (number_random_starts_selecting_one > 0L) {
  cat(paste0(
    "Conclusion: random starts found the K = 1 boundary, but maximum ",
    "training ELBO still selected K = 2.\n"
  ))
} else {
  cat("Conclusion: these random starts did not reach the K = 1 boundary.\n")
}

cat("\nHighest-ELBO initialization:", best_initialization, "\n")
cat("Learned alpha0:\n")
print(veb_fit$alpha0, digits = 7)
cat("Learned beta0:\n")
print(veb_fit$beta0, digits = 7)
cat("Learned prior means alpha0 / beta0:\n")
print(veb_fit$relevance, digits = 7)
cat("Learned prior coefficients of variation 1 / sqrt(alpha0):\n")
print(1 / sqrt(veb_fit$alpha0), digits = 7)
cat("Relative prior-mean shares:\n")
print(veb_fit$relevance_share, digits = 7)
cat("Allocated count shares:\n")
print(veb_fit$allocated_share, digits = 7)
cat("Effective K (reporting threshold only):", veb_fit$effective_K, "\n")

if (interactive()) {
  old_graphics = par(no.readonly = TRUE)
  on.exit(par(old_graphics), add = TRUE)
  par(mfrow = c(1, 3))

  matplot(
    veb_fit$beta0_history,
    type = "l",
    lty = 1,
    log = "y",
    xlab = "VEB iteration",
    ylab = "Learned prior rate beta0",
    main = "Prior rates"
  )
  legend("topleft", paste0("factor ", seq_len(K_fit)),
         col = seq_len(K_fit), lty = 1, bty = "n")

  matplot(
    veb_fit$relevance_history,
    type = "l",
    lty = 1,
    log = "y",
    xlab = "VEB iteration",
    ylab = "alpha0 / beta0",
    main = "Prior means"
  )
  legend("topright", paste0("factor ", seq_len(K_fit)),
         col = seq_len(K_fit), lty = 1, bty = "n")

  plot(
    veb_fit$elbo,
    type = "l",
    lwd = 2,
    xlab = "VEB iteration",
    ylab = "ELBO",
    main = paste("ELBO:", best_initialization)
  )
}

cat("\nUseful RStudio objects:\n")
cat("  fit_summary\n")
cat("  explicit_k1_fit\n")
cat("  veb_fits\n")
cat("  veb_fit$beta0_history\n")
cat("  veb_fit$relevance_history\n")
cat("  veb_fit$allocated_share_history\n")
cat("  View(veb_fit$F)\n")
cat("  plot(veb_fit$beta0_history[, 1], type = 'l', log = 'y')\n")
cat("  lines(veb_fit$beta0_history[, 2], col = 2)\n")
