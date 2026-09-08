## Repeat the simple-tree score-threshold experiment with
## q(pi) = Dirichlet(delta) and delta_prior = (1, ..., 1).

source("code/miso-benchmark-utils.R")
source("code/miso-dirichlet-pi.R")

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
DIRICHLET_PRIOR = 1

TRUNCATE_LEVELS = as.numeric(strsplit(
  Sys.getenv("MISO_TRUNCATE_LEVELS", unset = "0.50,0.60,0.70,0.80,0.90,0.95"),
  split = ",",
  fixed = TRUE
)[[1]])
MISO_MAX_ITERS = as.integer(Sys.getenv("MISO_MAX_ITERS", unset = "100"))
MISO_N_INNER = 3L
MIXTURE_MAX_ITERS = 20L
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
      "simple-tree-dirichlet1-%s-tau%03d-seed%d-iter%d.rds",
      score_tag,
      round(100 * truncate_level),
      FITTING_SEED,
      MISO_MAX_ITERS
    )
  )

  if (file.exists(fit_cache)) return(readRDS(fit_cache))

  message(
    "Fitting Dirichlet-pi MiSo: score = ", score_name,
    ", truncate level = ", format(truncate_level, nsmall = 2),
    ", initial S* = ", length(initialization$cluster_size)
  )
  fit = miso_fixed_gamma_dirichlet_pi(
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
    dirichlet_prior = DIRICHLET_PRIOR,
    mixture_max_iters = MIXTURE_MAX_ITERS,
    mixture_tol = 1e-10,
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

  data.frame(
    score = score_name,
    truncate_level = truncate_level,
    initial_S_star = ncol(fit$omega),
    allocation_components_gt_0.01 = sum(fit$allocation_mass > 0.01),
    posterior_mean_pi_gt_0.01 = sum(fit$pi > 0.01),
    converged = fit$converged,
    iterations = fit$n_iter,
    final_elbo = tail(fit$elbo[is.finite(fit$elbo)], 1),
    true_mean_deviance_per_entry = poisson_deviance(
      true_mean, fitted_mean
    ) / length(true_mean),
    relative_true_mean_error = sqrt(
      sum((fitted_mean - true_mean)^2) / sum(true_mean^2)
    ),
    soft_motif_accuracy = rectangular_cluster_accuracy(fit$omega, dat$grp),
    cosine_true_F1 = factor_cosine_by_true[1],
    cosine_true_F2 = factor_cosine_by_true[2],
    cosine_true_F3 = factor_cosine_by_true[3],
    mean_factor_cosine = mean(factor_match$table$cosine),
    allocation_mass_decreasing = paste(
      sprintf("%.8f", sort(fit$allocation_mass, decreasing = TRUE)),
      collapse = ","
    ),
    posterior_mean_pi_decreasing = paste(
      sprintf("%.8f", sort(fit$pi, decreasing = TRUE)),
      collapse = ","
    ),
    posterior_expected_count_decreasing = paste(
      sprintf(
        "%.4f",
        sort(fit$dirichlet_posterior - fit$dirichlet_prior,
             decreasing = TRUE)
      ),
      collapse = ","
    )
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

prior_means = do.call(rbind, lapply(seq_len(nrow(candidate_grid)), function(j) {
  fit = candidates[[j]]$fit
  initialization = candidates[[j]]$initialization
  cluster_order = order(fit$allocation_mass, decreasing = TRUE)

  do.call(rbind, lapply(seq_along(cluster_order), function(cluster_rank) {
    s = cluster_order[cluster_rank]
    do.call(rbind, lapply(seq_len(ncol(fit$alpha0)), function(d) {
      data.frame(
        score = candidate_grid$score[j],
        truncate_level = candidate_grid$truncate_level[j],
        initial_S_star = length(cluster_order),
        cluster_mass_rank = cluster_rank,
        original_cluster = s,
        initial_cluster_size = initialization$cluster_size[s],
        allocation_mass = fit$allocation_mass[s],
        posterior_mean_pi = fit$pi[s],
        dimension = d,
        alpha0 = fit$alpha0[s, d],
        beta0 = fit$beta0[s, d],
        prior_mean_alpha0_over_beta0 = fit$alpha0[s, d] / fit$beta0[s, d],
        top_factor = which.max(fit$gamma_bar[s, d, ]),
        top_factor_probability = max(fit$gamma_bar[s, d, ])
      )
    }))
  }))
}))
row.names(prior_means) = NULL

write.csv(
  results,
  file.path(cache_dir, "simple-tree-dirichlet1-score-tuning-results.csv"),
  row.names = FALSE
)
write.csv(
  prior_means,
  file.path(cache_dir, "simple-tree-dirichlet1-prior-means.csv"),
  row.names = FALSE
)

point_results_file = file.path(
  cache_dir, "simple-tree-score-tuning-results.csv"
)
if (file.exists(point_results_file)) {
  point_results = read.csv(point_results_file, stringsAsFactors = FALSE)
  comparison = merge(
    point_results,
    results,
    by = c("score", "truncate_level"),
    suffixes = c("_point_pi", "_dirichlet_pi")
  )
  comparison = comparison[, c(
    "score",
    "truncate_level",
    "initial_S_star",
    "effective_S_pi_gt_0.01",
    "allocation_components_gt_0.01",
    "posterior_mean_pi_gt_0.01",
    "soft_motif_accuracy_point_pi",
    "soft_motif_accuracy_dirichlet_pi",
    "final_pi_decreasing",
    "allocation_mass_decreasing",
    "posterior_mean_pi_decreasing",
    "true_mean_deviance_per_entry_point_pi",
    "true_mean_deviance_per_entry_dirichlet_pi",
    "relative_true_mean_error_point_pi",
    "relative_true_mean_error_dirichlet_pi"
  )]
  names(comparison)[names(comparison) == "effective_S_pi_gt_0.01"] =
    "point_pi_components_gt_0.01"
  names(comparison)[names(comparison) == "final_pi_decreasing"] =
    "point_pi_decreasing"
  write.csv(
    comparison,
    file.path(cache_dir, "simple-tree-point-vs-dirichlet1-pi.csv"),
    row.names = FALSE
  )
}

print(
  results[order(results$score, results$truncate_level), ],
  digits = 6,
  row.names = FALSE
)
