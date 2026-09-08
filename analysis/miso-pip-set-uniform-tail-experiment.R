## Known-K experiment for score-set MiSo initializations.
##
## This script fits the anchor and tree simulations at their true K and D.
## S is not supplied: it is the number of distinct truncated row-wise
## Poisson-SuSiE score sets. MISO_SUPPORT_SCORE selects either ordinary PIP or
## posterior expected loading R_ik. Factors inside a set initialize
## concentrated motif slots; unsupported tail slots are uniform over factors.

source("code/miso-benchmark-utils.R")

PIP_TRUNCATE_LEVEL = as.numeric(
  Sys.getenv("MISO_PIP_TRUNCATE_LEVEL", unset = "0.95")
)
SUPPORT_SCORE = match.arg(
  Sys.getenv("MISO_SUPPORT_SCORE", unset = "pip"),
  c("pip", "expected_loading")
)
CACHE_PREFIX = if (SUPPORT_SCORE == "pip") "pip-set" else "loading-set"
INITIALIZATION_LABEL = if (SUPPORT_SCORE == "pip") {
  "truncated PIP sets + uniform tail"
} else {
  "truncated expected-loading sets + uniform tail"
}
MISO_MAX_ITERS = as.integer(
  Sys.getenv("MISO_MAX_ITERS", unset = "30")
)
MISO_N_INNER = 3L
PS_MAX_ITERS = 30L
NMF_MAX_ITERS = 60L

cache_dir = file.path("analysis", "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

anchor_supports = list(
  1, 2, 3, 4, 5,
  c(1, 2, 3), c(1, 2, 4), c(2, 3, 5)
)
tree_supports = list(
  c(1, 2, 3), c(1, 2, 4), c(1, 5, 6),
  c(1, 5, 7), c(1, 2), c(1, 5)
)

scenario_configuration = list(
  anchor = list(
    supports = anchor_supports,
    N = 480L,
    M = 500L,
    K = 5L,
    D = 3L,
    shape = c(rep(10, 5), rep(20, 3)),
    simulation_seed = 42L,
    fitting_seed = 1L
  ),
  tree = list(
    supports = tree_supports,
    N = 480L,
    M = 500L,
    K = 7L,
    D = 3L,
    shape = rep(18, length(tree_supports)),
    simulation_seed = 47L,
    fitting_seed = 2L
  )
)

get_preliminary_fit = function(scenario, configuration) {
  cache_file = file.path(
    cache_dir,
    sprintf(
      "pip-set-%s-preliminary-K%d-seed%d.rds",
      scenario, configuration$K, configuration$fitting_seed
    )
  )
  if (file.exists(cache_file)) return(readRDS(cache_file))

  dat = simulate_support_scenario(
    N = configuration$N,
    M = configuration$M,
    K = configuration$K,
    group_supports = configuration$supports,
    shape_vec = configuration$shape,
    seed = configuration$simulation_seed
  )
  nmf = poisson_nmf_init(
    dat$Y,
    K = configuration$K,
    max_iters = NMF_MAX_ITERS,
    init_seed = configuration$fitting_seed
  )
  ps = poisson_susie_nmf_fixed_F(
    Y = dat$Y,
    F = nmf$F,
    D = configuration$D,
    max_iters = PS_MAX_ITERS,
    update_prior = TRUE,
    init_seed = configuration$fitting_seed,
    init_L = nmf$L,
    gamma_init_floor = 1e-4,
    elbo_every = 5,
    tol = 5e-4,
    min_iters = 15,
    patience = 2
  )
  out = list(dat = dat, nmf = nmf, ps = ps)
  saveRDS(out, cache_file)
  out
}

fit_score_set_miso = function(scenario, configuration) {
  preliminary = get_preliminary_fit(scenario, configuration)
  initialization = init_pip_set_miso(
    preliminary$ps,
    pip_truncate_level = PIP_TRUNCATE_LEVEL,
    support_score = SUPPORT_SCORE,
    gamma_floor = 0.05
  )
  fit_cache = file.path(
    cache_dir,
    sprintf(
      "%s-%s-K%d-D%d-tau%03d-seed%d-iter%d.rds",
      CACHE_PREFIX,
      scenario,
      configuration$K,
      configuration$D,
      round(100 * PIP_TRUNCATE_LEVEL),
      configuration$fitting_seed,
      MISO_MAX_ITERS
    )
  )

  if (file.exists(fit_cache)) {
    fit = readRDS(fit_cache)
  } else {
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
    fit$gamma_init = initialization
    fit$ps_fit = preliminary$ps
    fit$nmf_fit = preliminary$nmf
    fit$mf_fit = preliminary$ps
    fit$motif_initialization = initialization$motif_initialization
    saveRDS(fit, fit_cache)
  }

  list(
    dat = preliminary$dat,
    initialization = initialization,
    fit = fit,
    supports = configuration$supports
  )
}

maximum_rectangular_assignment = function(score) {
  n = max(nrow(score), ncol(score))
  padded = matrix(0, n, n)
  padded[seq_len(nrow(score)), seq_len(ncol(score))] = score
  if (requireNamespace("clue", quietly = TRUE)) {
    return(as.integer(clue::solve_LSAP(padded, maximum = TRUE)))
  }

  assignment = rep(NA_integer_, n)
  available_rows = seq_len(n)
  available_columns = seq_len(n)
  while (length(available_rows) > 0) {
    candidate = which(
      padded[available_rows, available_columns, drop = FALSE] ==
        max(padded[available_rows, available_columns, drop = FALSE]),
      arr.ind = TRUE
    )[1, ]
    row = available_rows[candidate[1]]
    column = available_columns[candidate[2]]
    assignment[row] = column
    available_rows = setdiff(available_rows, row)
    available_columns = setdiff(available_columns, column)
  }
  assignment
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

summarize_fit = function(scenario, experiment) {
  dat = experiment$dat
  fit = experiment$fit
  initialization = experiment$initialization
  factor_match = match_factor_rows(fit$F, dat$F0)
  observation_scores = miso_observation_factor_scores(fit)
  fitted_mean = observation_scores %*% fit$F

  hard_responsibility = matrix(0, nrow(fit$omega), ncol(fit$omega))
  hard_responsibility[cbind(seq_len(nrow(fit$omega)), max.col(fit$omega))] = 1
  initial_prior_mean = fit$alpha0_initial / fit$beta0_initial
  final_prior_mean = fit$alpha0 / fit$beta0
  elbo_value = tail(fit$elbo[!is.na(fit$elbo)], 1)

  data.frame(
    scenario = scenario,
    initialization = INITIALIZATION_LABEL,
    support_score = SUPPORT_SCORE,
    pip_truncate_level = PIP_TRUNCATE_LEVEL,
    true_K = nrow(dat$F0),
    fitted_K = nrow(fit$F),
    true_S = length(experiment$supports),
    inferred_S = ncol(fit$omega),
    true_D = max(lengths(experiment$supports)),
    fitted_D = dim(fit$gamma_bar)[2],
    clusters_with_uniform_tail = sum(initialization$support_size <
                                       dim(fit$gamma_bar)[2]),
    minimum_cluster_size = min(initialization$cluster_size),
    median_cluster_size = median(initialization$cluster_size),
    maximum_cluster_size = max(initialization$cluster_size),
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
    mean_factor_cosine = mean(factor_match$table$cosine),
    minimum_factor_cosine = min(factor_match$table$cosine),
    mean_assignment_certainty = mean(apply(fit$omega, 1, max)),
    effective_motifs_pi_gt_0.01 = sum(fit$pi > 0.01),
    final_elbo = elbo_value,
    deviance_per_entry = poisson_deviance(dat$Y, fitted_mean) / length(dat$Y),
    iterations = fit$n_iter,
    initial_prior_mean_min = min(initial_prior_mean),
    initial_prior_mean_median = median(initial_prior_mean),
    initial_prior_mean_max = max(initial_prior_mean),
    final_prior_mean_min = min(final_prior_mean),
    final_prior_mean_median = median(final_prior_mean),
    final_prior_mean_max = max(final_prior_mean),
    mean_initial_gamma_entropy = mean(
      motif_gamma_entropy(initialization$gamma_bar)
    ),
    mean_final_gamma_entropy = mean(motif_gamma_entropy(fit$gamma_bar))
  )
}

experiments = lapply(names(scenario_configuration), function(scenario) {
  message("Fitting ", scenario, " scenario")
  fit_score_set_miso(scenario, scenario_configuration[[scenario]])
})
names(experiments) = names(scenario_configuration)

performance = do.call(rbind, lapply(names(experiments), function(scenario) {
  summarize_fit(scenario, experiments[[scenario]])
}))
row.names(performance) = NULL

cluster_details = do.call(rbind, lapply(names(experiments), function(scenario) {
  initialization = experiments[[scenario]]$initialization
  dat = experiments[[scenario]]$dat
  learned_to_true = match_factor_rows(
    experiments[[scenario]]$fit$nmf_fit$F, dat$F0
  )$learned_to_true
  group_count = table(initialization$cluster, dat$grp)
  do.call(rbind, lapply(seq_along(initialization$cluster_size), function(s) {
    count_s = group_count[s, ]
    represented_groups = which(count_s > 0)
    data.frame(
      scenario = scenario,
      cluster = s,
      cluster_size = initialization$cluster_size[s],
      support_size = initialization$support_size[s],
      learned_factor_support = paste(
        initialization$support_sets[[s]], collapse = ","
      ),
      matched_true_factor_support = paste(
        sort(learned_to_true[initialization$support_sets[[s]]]),
        collapse = ","
      ),
      number_uniform_tail_slots = scenario_configuration[[scenario]]$D -
        initialization$support_size[s],
      dominant_true_group = paste0("G", which.max(count_s)),
      cluster_purity = max(count_s) / sum(count_s),
      true_group_composition = paste(
        paste0("G", represented_groups, "=", count_s[represented_groups]),
        collapse = ";"
      )
    )
  }))
}))
row.names(cluster_details) = NULL

write.csv(
  performance,
  file.path(
    cache_dir,
    sprintf("%s-known-K-performance-tau%03d.csv", CACHE_PREFIX,
            round(100 * PIP_TRUNCATE_LEVEL))
  ),
  row.names = FALSE
)
write.csv(
  cluster_details,
  file.path(
    cache_dir,
    sprintf("%s-known-K-clusters-tau%03d.csv", CACHE_PREFIX,
            round(100 * PIP_TRUNCATE_LEVEL))
  ),
  row.names = FALSE
)

print(performance, digits = 4)
print(with(cluster_details, table(scenario, support_size)))
