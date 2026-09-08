## Simple-tree experiment with genuinely broader two-dimensional loadings.
## Active loadings are Gamma(shape = 2, rate = 1/9), so their mean remains 18.

source("code/miso-benchmark-utils.R")
source("code/miso-dirichlet-pi.R")

N = 480L
M = 500L
K_TRUE = 3L
D_TRUE = 2L
TRUE_SUPPORTS = list(c(1L, 2L), c(1L, 3L))
LOADING_SHAPE = c(2, 2)
LOADING_RATE = 1 / 9
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
    "simple-tree-gamma2-rate1over9-preliminary-N%d-M%d-K%d-seed%d.rds",
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
    loading_rate = LOADING_RATE,
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

fit_candidate = function(score_name, truncate_level,
                          pi_method = c("point", "dirichlet1")) {
  pi_method = match.arg(pi_method)
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
      "simple-tree-gamma2-rate1over9-%s-%s-tau%03d-seed%d-iter%d.rds",
      pi_method,
      score_tag,
      round(100 * truncate_level),
      FITTING_SEED,
      MISO_MAX_ITERS
    )
  )
  if (file.exists(fit_cache)) return(readRDS(fit_cache))

  message(
    "Fitting ", pi_method,
    ": score = ", score_name,
    ", truncate level = ", format(truncate_level, nsmall = 2),
    ", initial S* = ", length(initialization$cluster_size)
  )
  common_arguments = list(
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
  if (pi_method == "point") {
    fit = do.call(miso_fixed_gamma, common_arguments)
  } else {
    fit = do.call(
      miso_fixed_gamma_dirichlet_pi,
      c(
        common_arguments,
        list(
          dirichlet_prior = DIRICHLET_PRIOR,
          mixture_max_iters = MIXTURE_MAX_ITERS,
          mixture_tol = 1e-10
        )
      )
    )
  }
  saved = list(fit = fit, initialization = initialization)
  saveRDS(saved, fit_cache)
  saved
}

aggregate_prior_mean_by_true_group = function(fit, dat, learned_to_true) {
  S = nrow(fit$alpha0)
  K = nrow(fit$F)
  D = ncol(fit$alpha0)
  motif_mean_learned = matrix(0, S, K)
  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      motif_mean_learned[s, ] = motif_mean_learned[s, ] +
        (fit$alpha0[s, d] / fit$beta0[s, d]) * fit$gamma_bar[s, d, ]
    }
  }
  motif_mean_true = matrix(0, S, K)
  for (learned_k in seq_len(K)) {
    motif_mean_true[, learned_to_true[learned_k]] =
      motif_mean_learned[, learned_k]
  }

  group_weight = sapply(seq_along(TRUE_SUPPORTS), function(g) {
    colMeans(fit$omega[dat$grp == g, , drop = FALSE])
  })
  t(group_weight) %*% motif_mean_true
}

summarize_candidate = function(score_name, truncate_level, pi_method,
                                candidate) {
  fit = candidate$fit
  dat = preliminary$dat
  allocation_mass = if (pi_method == "point") fit$pi else fit$allocation_mass
  factor_match = match_factor_rows(fit$F, dat$F0)
  factor_cosine_by_true = factor_match$table$cosine[
    order(factor_match$table$matched_true_factor)
  ]
  observation_scores = miso_observation_factor_scores(fit)
  fitted_mean = observation_scores %*% fit$F
  true_mean = dat$L %*% dat$F0
  aggregate_prior_mean = aggregate_prior_mean_by_true_group(
    fit, dat, factor_match$learned_to_true
  )
  target_prior_mean = rbind(c(18, 18, 0), c(18, 0, 18))
  prior_mean_relative_error = sqrt(
    sum((aggregate_prior_mean - target_prior_mean)^2) /
      sum(target_prior_mean^2)
  )
  active_cluster = allocation_mass > 0.01
  cluster_prior_mean = fit$alpha0 / fit$beta0
  smaller_to_larger_prior_mean = apply(
    cluster_prior_mean, 1, function(x) min(x) / max(x)
  )

  data.frame(
    pi_method = pi_method,
    score = score_name,
    truncate_level = truncate_level,
    initial_S_star = ncol(fit$omega),
    active_S_mass_gt_0.01 = sum(active_cluster),
    active_clusters_with_two_prior_means_gt_1.8 = sum(
      active_cluster & rowSums(cluster_prior_mean > 1.8) == D_TRUE
    ),
    allocation_weighted_smaller_to_larger_prior_mean = sum(
      allocation_mass * smaller_to_larger_prior_mean
    ),
    converged = fit$converged,
    iterations = fit$n_iter,
    final_elbo = tail(fit$elbo[is.finite(fit$elbo)], 1),
    data_deviance_per_entry = poisson_deviance(
      dat$Y, fitted_mean
    ) / length(dat$Y),
    true_mean_deviance_per_entry = poisson_deviance(
      true_mean, fitted_mean
    ) / length(true_mean),
    relative_true_mean_error = sqrt(
      sum((fitted_mean - true_mean)^2) / sum(true_mean^2)
    ),
    soft_motif_accuracy = rectangular_cluster_accuracy(fit$omega, dat$grp),
    exact_support_accuracy = support_accuracy(
      observation_scores, factor_match$learned_to_true, dat
    ),
    cosine_true_F1 = factor_cosine_by_true[1],
    cosine_true_F2 = factor_cosine_by_true[2],
    cosine_true_F3 = factor_cosine_by_true[3],
    mean_factor_cosine = mean(factor_match$table$cosine),
    aggregate_prior_mean_relative_error = prior_mean_relative_error,
    group1_prior_means_true_factor_order = paste(
      sprintf("%.5f", aggregate_prior_mean[1, ]), collapse = ","
    ),
    group2_prior_means_true_factor_order = paste(
      sprintf("%.5f", aggregate_prior_mean[2, ]), collapse = ","
    ),
    allocation_mass_decreasing = paste(
      sprintf("%.8f", sort(allocation_mass, decreasing = TRUE)),
      collapse = ","
    ),
    cluster_prior_mean_pairs_by_mass = paste(
      apply(
        cluster_prior_mean[order(allocation_mass, decreasing = TRUE), ,
                           drop = FALSE],
        1,
        function(x) sprintf("(%.4f;%.4f)", x[1], x[2])
      ),
      collapse = ","
    )
  )
}

candidate_grid = expand.grid(
  pi_method = c("point", "dirichlet1"),
  score = c("pip", "expected_loading"),
  truncate_level = TRUNCATE_LEVELS,
  stringsAsFactors = FALSE
)

candidates = lapply(seq_len(nrow(candidate_grid)), function(j) {
  fit_candidate(
    score_name = candidate_grid$score[j],
    truncate_level = candidate_grid$truncate_level[j],
    pi_method = candidate_grid$pi_method[j]
  )
})

results = do.call(rbind, lapply(seq_len(nrow(candidate_grid)), function(j) {
  summarize_candidate(
    score_name = candidate_grid$score[j],
    truncate_level = candidate_grid$truncate_level[j],
    pi_method = candidate_grid$pi_method[j],
    candidate = candidates[[j]]
  )
}))
row.names(results) = NULL

true_composition = unlist(lapply(seq_along(TRUE_SUPPORTS), function(g) {
  support = TRUE_SUPPORTS[[g]]
  rows = preliminary$dat$grp == g
  preliminary$dat$L[rows, support[1]] /
    rowSums(preliminary$dat$L[rows, support, drop = FALSE])
}))
true_loading_singular_value_ratio = sapply(
  seq_along(TRUE_SUPPORTS), function(g) {
    rows = preliminary$dat$grp == g
    singular_values = svd(
      preliminary$dat$L[rows, TRUE_SUPPORTS[[g]], drop = FALSE],
      nu = 0,
      nv = 0
    )$d
    singular_values[2] / singular_values[1]
  }
)
true_mean_singular_value_ratio = sapply(
  seq_along(TRUE_SUPPORTS), function(g) {
    rows = preliminary$dat$grp == g
    noiseless_mean = preliminary$dat$L[
      rows, TRUE_SUPPORTS[[g]], drop = FALSE
    ] %*% preliminary$dat$F0[TRUE_SUPPORTS[[g]], , drop = FALSE]
    singular_values = svd(noiseless_mean, nu = 0, nv = 0)$d
    singular_values[2] / singular_values[1]
  }
)
true_loading_geometry = data.frame(
  empirical_active_loading_mean = mean(preliminary$dat$L[preliminary$dat$L > 0]),
  empirical_total_loading_mean = mean(rowSums(preliminary$dat$L)),
  composition_q05 = unname(quantile(true_composition, 0.05)),
  composition_median = median(true_composition),
  composition_q95 = unname(quantile(true_composition, 0.95)),
  group1_loading_singular_value_ratio = true_loading_singular_value_ratio[1],
  group2_loading_singular_value_ratio = true_loading_singular_value_ratio[2],
  group1_mean_singular_value_ratio = true_mean_singular_value_ratio[1],
  group2_mean_singular_value_ratio = true_mean_singular_value_ratio[2]
)

true_mean = preliminary$dat$L %*% preliminary$dat$F0
nmf_mean = preliminary$nmf$L %*% preliminary$nmf$F
nmf_factor_match = match_factor_rows(preliminary$nmf$F, preliminary$dat$F0)
baseline = data.frame(
  model = "rank-3 vanilla NMF",
  true_mean_deviance_per_entry = poisson_deviance(
    true_mean, nmf_mean
  ) / length(true_mean),
  relative_true_mean_error = sqrt(
    sum((nmf_mean - true_mean)^2) / sum(true_mean^2)
  ),
  mean_factor_cosine = mean(nmf_factor_match$table$cosine),
  minimum_factor_cosine = min(nmf_factor_match$table$cosine)
)

write.csv(
  results,
  file.path(cache_dir, "simple-tree-gamma2-rate1over9-results.csv"),
  row.names = FALSE
)
write.csv(
  true_loading_geometry,
  file.path(cache_dir, "simple-tree-gamma2-rate1over9-geometry.csv"),
  row.names = FALSE
)
write.csv(
  baseline,
  file.path(cache_dir, "simple-tree-gamma2-rate1over9-baseline.csv"),
  row.names = FALSE
)

print(true_loading_geometry, digits = 6, row.names = FALSE)
print(baseline, digits = 6, row.names = FALSE)
print(
  results[order(results$pi_method, results$score, results$truncate_level), ],
  digits = 6,
  row.names = FALSE
)
