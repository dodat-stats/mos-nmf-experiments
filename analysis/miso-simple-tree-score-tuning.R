## Tune the PIP and posterior-expected-loading score thresholds in a simple
## tree simulation with K_true = 3, S_true = 2, and D_true = 2.

source("code/miso-benchmark-utils.R")

N = 480L
M = 500L
K_TRUE = 3L
D_TRUE = 2L
TRUE_SUPPORTS = list(c(1L, 2L), c(1L, 3L))
LOADING_SHAPE = c(18, 18)
FACTOR_SHAPE = 0.1
FACTOR_RATE = 0.01
SIMULATION_SEED = 52L
FITTING_SEED = 3L

TRUNCATE_LEVELS = as.numeric(strsplit(
  Sys.getenv("MISO_TRUNCATE_LEVELS", unset = "0.50,0.60,0.70,0.80,0.90,0.95"),
  split = ",",
  fixed = TRUE
)[[1]])
MISO_MAX_ITERS = as.integer(Sys.getenv("MISO_MAX_ITERS", unset = "100"))
MISO_N_INNER = 3L
PS_MAX_ITERS = 30L
NMF_MAX_ITERS = 60L

cache_dir = file.path("analysis", "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

maximum_rectangular_assignment = function(score) {
  n = max(nrow(score), ncol(score))
  padded = matrix(0, n, n)
  padded[seq_len(nrow(score)), seq_len(ncol(score))] = score
  if (requireNamespace("clue", quietly = TRUE)) {
    return(as.integer(clue::solve_LSAP(padded, maximum = TRUE)))
  }

  permutations = all_permutations(seq_len(n))
  objective = apply(permutations, 1, function(candidate) {
    sum(padded[cbind(seq_len(n), candidate)])
  })
  as.integer(permutations[which.max(objective), ])
}

rectangular_cluster_accuracy = function(responsibility, true_group) {
  true_group = as.integer(as.factor(true_group))
  S_fit = ncol(responsibility)
  S_true = max(true_group)
  score = matrix(0, S_fit, S_true)
  for (g in seq_len(S_true)) {
    score[, g] = colSums(
      responsibility[true_group == g, , drop = FALSE]
    )
  }
  assignment = maximum_rectangular_assignment(score)
  valid = which(assignment[seq_len(S_fit)] <= S_true)
  sum(score[cbind(valid, assignment[valid])]) / nrow(responsibility)
}

preliminary_cache = file.path(
  cache_dir,
  sprintf(
    "simple-tree-preliminary-N%d-M%d-K%d-seed%d.rds",
    N, M, K_TRUE, SIMULATION_SEED
  )
)

if (file.exists(preliminary_cache)) {
  preliminary = readRDS(preliminary_cache)
} else {
  dat = simulate_support_scenario(
    N = N,
    M = M,
    K = K_TRUE,
    group_supports = TRUE_SUPPORTS,
    shape_vec = LOADING_SHAPE,
    factor_shape = FACTOR_SHAPE,
    factor_rate = FACTOR_RATE,
    seed = SIMULATION_SEED
  )
  nmf = poisson_nmf_init(
    dat$Y,
    K = K_TRUE,
    max_iters = NMF_MAX_ITERS,
    init_seed = FITTING_SEED
  )
  ps = poisson_susie_nmf_fixed_F(
    Y = dat$Y,
    F = nmf$F,
    D = D_TRUE,
    max_iters = PS_MAX_ITERS,
    update_prior = TRUE,
    init_seed = FITTING_SEED,
    init_L = nmf$L,
    gamma_init_floor = 1e-4,
    elbo_every = 5,
    tol = 5e-4,
    min_iters = 15,
    patience = 2
  )
  preliminary = list(dat = dat, nmf = nmf, ps = ps)
  saveRDS(preliminary, preliminary_cache)
}

fit_candidate = function(score_name, truncate_level) {
  initialization = init_pip_set_miso(
    ps_fit = preliminary$ps,
    pip_truncate_level = truncate_level,
    support_score = score_name,
    gamma_floor = 0.05
  )
  score_tag = if (score_name == "pip") "pip" else "R"
  fit_cache = file.path(
    cache_dir,
    sprintf(
      "simple-tree-%s-tau%03d-seed%d-iter%d.rds",
      score_tag,
      round(100 * truncate_level),
      FITTING_SEED,
      MISO_MAX_ITERS
    )
  )

  if (file.exists(fit_cache)) {
    saved = readRDS(fit_cache)
    return(saved)
  }

  message(
    "Fitting score = ", score_name,
    ", truncate level = ", format(truncate_level, nsmall = 2),
    ", inferred S = ", length(initialization$cluster_size)
  )
  fit = miso_fixed_gamma(
    Y = preliminary$dat$Y,
    F = preliminary$nmf$F,
    gamma_bar = initialization$gamma_bar,
    max_iters = MISO_MAX_ITERS,
    n_inner = MISO_N_INNER,
    update_prior = TRUE,
    alpha0_init = initialization$alpha0,
    beta0_init = initialization$beta0,
    pi_init = initialization$pi_init,
    omega_init = initialization$omega_init,
    tol = 1e-5,
    min_iters = 5,
    patience = 2,
    update_F = TRUE,
    update_gamma = TRUE,
    update_slot_scale = FALSE,
    gamma_step_init = 0.5,
    gamma_step_ramp = 8,
    F_step_init = 0.2,
    F_step_ramp = 12,
    block_size = 100
  )
  saved = list(fit = fit, initialization = initialization)
  saveRDS(saved, fit_cache)
  saved
}

summarize_candidate = function(score_name, truncate_level, candidate) {
  fit = candidate$fit
  initialization = candidate$initialization
  dat = preliminary$dat
  factor_match = match_factor_rows(fit$F, dat$F0)
  factor_cosine_by_true = factor_match$table$cosine[
    order(factor_match$table$matched_true_factor)
  ]
  observation_scores = miso_observation_factor_scores(fit)
  fitted_mean = observation_scores %*% fit$F
  true_mean = dat$L %*% dat$F0
  hard_responsibility = matrix(0, nrow(fit$omega), ncol(fit$omega))
  hard_responsibility[
    cbind(seq_len(nrow(fit$omega)), max.col(fit$omega))
  ] = 1

  data.frame(
    score = score_name,
    truncate_level = truncate_level,
    inferred_S = ncol(fit$omega),
    effective_S_pi_gt_0.01 = sum(fit$pi > 0.01),
    converged = fit$converged,
    iterations = fit$n_iter,
    final_elbo = tail(fit$elbo[is.finite(fit$elbo)], 1),
    deviance_per_entry = poisson_deviance(dat$Y, fitted_mean) / length(dat$Y),
    true_mean_deviance_per_entry = poisson_deviance(
      true_mean, fitted_mean
    ) / length(true_mean),
    relative_true_mean_error = sqrt(
      sum((fitted_mean - true_mean)^2) / sum(true_mean^2)
    ),
    mean_absolute_true_mean_error = mean(abs(fitted_mean - true_mean)),
    soft_motif_accuracy = rectangular_cluster_accuracy(fit$omega, dat$grp),
    hard_motif_accuracy = rectangular_cluster_accuracy(
      hard_responsibility, dat$grp
    ),
    exact_support_accuracy = support_accuracy(
      observation_scores, factor_match$learned_to_true, dat
    ),
    support_jaccard = support_jaccard(
      observation_scores, factor_match$learned_to_true, dat
    ),
    cosine_true_F1 = factor_cosine_by_true[1],
    cosine_true_F2 = factor_cosine_by_true[2],
    cosine_true_F3 = factor_cosine_by_true[3],
    mean_factor_cosine = mean(factor_match$table$cosine),
    minimum_factor_cosine = min(factor_match$table$cosine),
    learned_to_true_factor_permutation = paste(
      factor_match$learned_to_true, collapse = ","
    ),
    final_pi_original_order = paste(
      sprintf("%.6f", fit$pi), collapse = ","
    ),
    final_pi_decreasing = paste(
      sprintf("%.6f", sort(fit$pi, decreasing = TRUE)), collapse = ","
    ),
    minimum_cluster_size = min(initialization$cluster_size),
    median_cluster_size = median(initialization$cluster_size),
    maximum_cluster_size = max(initialization$cluster_size)
  )
}

candidate_grid = expand.grid(
  score = c("pip", "expected_loading"),
  truncate_level = TRUNCATE_LEVELS,
  stringsAsFactors = FALSE
)

candidates = lapply(seq_len(nrow(candidate_grid)), function(j) {
  fit_candidate(
    score_name = candidate_grid$score[j],
    truncate_level = candidate_grid$truncate_level[j]
  )
})

results = do.call(rbind, lapply(seq_len(nrow(candidate_grid)), function(j) {
  summarize_candidate(
    score_name = candidate_grid$score[j],
    truncate_level = candidate_grid$truncate_level[j],
    candidate = candidates[[j]]
  )
}))
row.names(results) = NULL

true_mean = preliminary$dat$L %*% preliminary$dat$F0
nmf_mean = preliminary$nmf$L %*% preliminary$nmf$F
benchmarks = data.frame(
  model = c("true Poisson mean", "rank-3 vanilla NMF"),
  deviance_per_entry = c(
    poisson_deviance(preliminary$dat$Y, true_mean) / length(preliminary$dat$Y),
    poisson_deviance(preliminary$dat$Y, nmf_mean) / length(preliminary$dat$Y)
  ),
  true_mean_deviance_per_entry = c(
    0,
    poisson_deviance(true_mean, nmf_mean) / length(true_mean)
  ),
  relative_true_mean_error = c(
    0,
    sqrt(sum((nmf_mean - true_mean)^2) / sum(true_mean^2))
  )
)

write.csv(
  results,
  file.path(cache_dir, "simple-tree-score-tuning-results.csv"),
  row.names = FALSE
)
write.csv(
  benchmarks,
  file.path(cache_dir, "simple-tree-reconstruction-benchmarks.csv"),
  row.names = FALSE
)

print(benchmarks, digits = 5, row.names = FALSE)
print(results[order(-results$final_elbo), ], digits = 5, row.names = FALSE)
