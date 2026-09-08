## Does the small redundant cluster persist after doubling both N and M?
##
## Simple tree:
##   K = 3, D = 2, true supports {1, 2} and {1, 3}.
##   Active loadings are Gamma(shape = 2, rate = 1/9), with mean 18.
##
## This script compares the PIP-tau-0.70 S* fit against an S = 2 fit obtained
## by merging its smallest cluster into the most similar larger cluster.  The
## merge only constructs an initialization; all MiSo parameters are then refit.

source("code/miso-benchmark-utils.R")
source("code/miso-dirichlet-pi.R")

N = 960L
M = 1000L
K_TRUE = 3L
D_TRUE = 2L
TRUE_SUPPORTS = list(c(1L, 2L), c(1L, 3L))
LOADING_SHAPE = c(2, 2)
LOADING_RATE = 1 / 9
FACTOR_SHAPE = 0.1
FACTOR_RATE = 0.01
SIMULATION_SEED = 52L
FITTING_SEED = 3L
PIP_TRUNCATE_LEVEL = 0.70
DIRICHLET_PRIOR = 1

MISO_MAX_ITERS = as.integer(Sys.getenv("MISO_MAX_ITERS", unset = "100"))
MISO_N_INNER = 3L
MIXTURE_MAX_ITERS = 20L
PS_MAX_ITERS = 30L
NMF_MAX_ITERS = 60L

cache_dir = file.path("analysis", "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
cache_tag = sprintf("simple-tree-gamma2-N%d-M%d-pip070", N, M)

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

preliminary_cache = file.path(cache_dir, paste0(cache_tag, "-preliminary.rds"))
if (file.exists(preliminary_cache)) {
  preliminary = readRDS(preliminary_cache)
} else {
  message("Simulating N = ", N, ", M = ", M, " and fitting NMF/Poisson-SuSiE")
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

initialization = init_pip_set_miso(
  ps_fit = preliminary$ps,
  pip_truncate_level = PIP_TRUNCATE_LEVEL,
  support_score = "pip",
  gamma_floor = 0.05
)
message(
  "PIP tau = ", PIP_TRUNCATE_LEVEL,
  " produces initial S* = ", length(initialization$cluster_size),
  "; cluster sizes = ", paste(sort(initialization$cluster_size,
                                    decreasing = TRUE), collapse = ", ")
)

fit_miso_from_initialization = function(initialization, pi_method, stage) {
  cache_file = file.path(
    cache_dir,
    sprintf("%s-%s-%s-iter%d.rds", cache_tag, pi_method, stage,
            MISO_MAX_ITERS)
  )
  if (file.exists(cache_file)) return(readRDS(cache_file))

  message("Full MiSo refit: ", pi_method, ", ", stage)
  common_arguments = list(
    Y = preliminary$dat$Y,
    F = initialization$F,
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
      c(common_arguments, list(
        dirichlet_prior = DIRICHLET_PRIOR,
        mixture_max_iters = MIXTURE_MAX_ITERS,
        mixture_tol = 1e-10
      ))
    )
  }
  saveRDS(fit, cache_file)
  fit
}

base_initialization = list(
  F = preliminary$nmf$F,
  gamma_bar = initialization$gamma_bar,
  alpha0 = initialization$alpha0,
  beta0 = initialization$beta0,
  pi_init = initialization$pi_init,
  omega_init = initialization$omega_init
)
base_fits = list(
  point = fit_miso_from_initialization(base_initialization, "point", "unmerged"),
  dirichlet1 = fit_miso_from_initialization(
    base_initialization, "dirichlet1", "unmerged"
  )
)

merge_miso_clusters = function(fit, target_cluster, source_cluster,
                               eps = 1e-12) {
  S = ncol(fit$omega)
  D = ncol(fit$alpha0)
  if (S < 2 || target_cluster == source_cluster ||
      !all(c(target_cluster, source_cluster) %in% seq_len(S))) {
    stop("target_cluster and source_cluster must be distinct valid clusters.")
  }

  target_mean = fit$alpha0[target_cluster, ] / fit$beta0[target_cluster, ]
  source_mean = fit$alpha0[source_cluster, ] / fit$beta0[source_cluster, ]
  target_signature = target_mean * fit$gamma_bar[target_cluster, , ]
  source_signature = source_mean * fit$gamma_bar[source_cluster, , ]
  slot_permutations = all_permutations(seq_len(D))
  alignment_score = apply(slot_permutations, 1, function(permutation) {
    sum(target_signature * source_signature[permutation, , drop = FALSE])
  })
  source_slot_order = slot_permutations[which.max(alignment_score), ]

  retained_cluster = setdiff(seq_len(S), source_cluster)
  target_new = which(retained_cluster == target_cluster)
  gamma_bar = fit$gamma_bar[retained_cluster, , , drop = FALSE]
  alpha0 = fit$alpha0[retained_cluster, , drop = FALSE]
  beta0 = fit$beta0[retained_cluster, , drop = FALSE]
  omega = fit$omega[, retained_cluster, drop = FALSE]
  omega[, target_new] = omega[, target_new] + fit$omega[, source_cluster]
  pi_init = colMeans(omega)

  target_mass = sum(fit$omega[, target_cluster])
  source_mass = sum(fit$omega[, source_cluster])
  cluster_weight = c(target_mass, source_mass) /
    (target_mass + source_mass)

  for (d in seq_len(D)) {
    source_d = source_slot_order[d]
    component_mean = c(
      fit$alpha0[target_cluster, d] / fit$beta0[target_cluster, d],
      fit$alpha0[source_cluster, source_d] /
        fit$beta0[source_cluster, source_d]
    )
    component_variance = c(
      fit$alpha0[target_cluster, d] / fit$beta0[target_cluster, d]^2,
      fit$alpha0[source_cluster, source_d] /
        fit$beta0[source_cluster, source_d]^2
    )
    merged_mean = sum(cluster_weight * component_mean)
    merged_second_moment = sum(
      cluster_weight * (component_variance + component_mean^2)
    )
    merged_variance = pmax(merged_second_moment - merged_mean^2, eps)
    alpha0[target_new, d] = merged_mean^2 / merged_variance
    beta0[target_new, d] = merged_mean / merged_variance

    expected_factor_loading =
      cluster_weight[1] * component_mean[1] *
        fit$gamma_bar[target_cluster, d, ] +
      cluster_weight[2] * component_mean[2] *
        fit$gamma_bar[source_cluster, source_d, ]
    if (sum(expected_factor_loading) <= eps) {
      gamma_bar[target_new, d, ] =
        cluster_weight[1] * fit$gamma_bar[target_cluster, d, ] +
        cluster_weight[2] * fit$gamma_bar[source_cluster, source_d, ]
    } else {
      gamma_bar[target_new, d, ] =
        expected_factor_loading / sum(expected_factor_loading)
    }
  }

  list(
    F = fit$F,
    gamma_bar = normalize_gamma_bar(gamma_bar),
    alpha0 = alpha0,
    beta0 = beta0,
    omega_init = omega,
    pi_init = pi_init,
    source_cluster = source_cluster,
    target_cluster = target_cluster,
    source_slot_order = source_slot_order,
    target_weight = cluster_weight[1],
    source_weight = cluster_weight[2]
  )
}

cluster_factor_signature = function(fit, cluster) {
  colSums(
    (fit$alpha0[cluster, ] / fit$beta0[cluster, ]) *
      fit$gamma_bar[cluster, , , drop = FALSE][1, , ]
  )
}

choose_merge = function(fit, pi_method) {
  mass = if (pi_method == "point") fit$pi else fit$allocation_mass
  source_cluster = which.min(mass)
  candidates = setdiff(seq_along(mass), source_cluster)
  source_signature = cluster_factor_signature(fit, source_cluster)
  cosine = sapply(candidates, function(target_cluster) {
    target_signature = cluster_factor_signature(fit, target_cluster)
    sum(source_signature * target_signature) /
      sqrt(sum(source_signature^2) * sum(target_signature^2))
  })
  target_cluster = candidates[which.max(cosine)]
  list(
    source_cluster = source_cluster,
    target_cluster = target_cluster,
    source_mass = mass[source_cluster],
    target_mass = mass[target_cluster],
    signature_cosine = cosine,
    candidates = candidates
  )
}

merge_choices = lapply(names(base_fits), function(pi_method) {
  fit = base_fits[[pi_method]]
  if (ncol(fit$omega) != 3L) {
    stop("This focused experiment expects PIP tau = 0.70 to produce S* = 3.")
  }
  choose_merge(fit, pi_method)
})
names(merge_choices) = names(base_fits)

merged_initializations = lapply(names(base_fits), function(pi_method) {
  choice = merge_choices[[pi_method]]
  merge_miso_clusters(
    base_fits[[pi_method]],
    target_cluster = choice$target_cluster,
    source_cluster = choice$source_cluster
  )
})
names(merged_initializations) = names(base_fits)

merged_fits = list(
  point = fit_miso_from_initialization(
    merged_initializations$point, "point", "merged-S2"
  ),
  dirichlet1 = fit_miso_from_initialization(
    merged_initializations$dirichlet1, "dirichlet1", "merged-S2"
  )
)

aggregate_prior_mean_by_true_group = function(fit, learned_to_true) {
  S = nrow(fit$alpha0)
  motif_mean_learned = matrix(0, S, K_TRUE)
  for (s in seq_len(S)) {
    for (d in seq_len(D_TRUE)) {
      motif_mean_learned[s, ] = motif_mean_learned[s, ] +
        (fit$alpha0[s, d] / fit$beta0[s, d]) * fit$gamma_bar[s, d, ]
    }
  }
  motif_mean_true = matrix(0, S, K_TRUE)
  for (learned_k in seq_len(K_TRUE)) {
    motif_mean_true[, learned_to_true[learned_k]] =
      motif_mean_learned[, learned_k]
  }
  group_weight = sapply(seq_along(TRUE_SUPPORTS), function(g) {
    colMeans(fit$omega[preliminary$dat$grp == g, , drop = FALSE])
  })
  t(group_weight) %*% motif_mean_true
}

summarize_fit = function(fit, pi_method, stage) {
  dat = preliminary$dat
  mass = if (pi_method == "point") fit$pi else fit$allocation_mass
  factor_match = match_factor_rows(fit$F, dat$F0)
  factor_cosine_by_true = factor_match$table$cosine[
    order(factor_match$table$matched_true_factor)
  ]
  observation_scores = miso_observation_factor_scores(fit)
  fitted_mean = observation_scores %*% fit$F
  true_mean = dat$L %*% dat$F0
  aggregate_prior_mean = aggregate_prior_mean_by_true_group(
    fit, factor_match$learned_to_true
  )
  target_prior_mean = rbind(c(18, 18, 0), c(18, 0, 18))

  data.frame(
    pi_method = pi_method,
    stage = stage,
    S = ncol(fit$omega),
    converged = fit$converged,
    iterations = fit$n_iter,
    final_elbo = tail(fit$elbo[is.finite(fit$elbo)], 1),
    data_deviance_per_entry = poisson_deviance(dat$Y, fitted_mean) /
      length(dat$Y),
    true_mean_deviance_per_entry = poisson_deviance(true_mean, fitted_mean) /
      length(true_mean),
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
    aggregate_prior_mean_relative_error = sqrt(
      sum((aggregate_prior_mean - target_prior_mean)^2) /
        sum(target_prior_mean^2)
    ),
    allocation_mass_decreasing = paste(
      sprintf("%.8f", sort(mass, decreasing = TRUE)), collapse = ","
    )
  )
}

results = do.call(rbind, lapply(names(base_fits), function(pi_method) {
  rbind(
    summarize_fit(base_fits[[pi_method]], pi_method, "unmerged"),
    summarize_fit(merged_fits[[pi_method]], pi_method, "merged_and_refit")
  )
}))
row.names(results) = NULL
results$elbo_change_from_unmerged = NA_real_
for (pi_method in names(base_fits)) {
  rows = results$pi_method == pi_method
  results$elbo_change_from_unmerged[rows] =
    results$final_elbo[rows] -
    results$final_elbo[rows & results$stage == "unmerged"]
}

merge_table = do.call(rbind, lapply(names(merge_choices), function(pi_method) {
  choice = merge_choices[[pi_method]]
  data.frame(
    pi_method = pi_method,
    source_cluster = choice$source_cluster,
    source_mass = choice$source_mass,
    target_cluster = choice$target_cluster,
    target_mass = choice$target_mass,
    candidate_clusters = paste(choice$candidates, collapse = ","),
    signature_cosines = paste(sprintf("%.6f", choice$signature_cosine),
                              collapse = ",")
  )
}))
row.names(merge_table) = NULL

cluster_parameter_table = function(fit, pi_method, stage) {
  mass = if (pi_method == "point") fit$pi else fit$allocation_mass
  factor_match = match_factor_rows(fit$F, preliminary$dat$F0)
  cluster_order = order(mass, decreasing = TRUE)
  do.call(rbind, lapply(seq_along(cluster_order), function(cluster_rank) {
    s = cluster_order[cluster_rank]
    do.call(rbind, lapply(seq_len(D_TRUE), function(d) {
      gamma_true_order = numeric(K_TRUE)
      for (learned_k in seq_len(K_TRUE)) {
        gamma_true_order[factor_match$learned_to_true[learned_k]] =
          fit$gamma_bar[s, d, learned_k]
      }
      data.frame(
        pi_method = pi_method,
        stage = stage,
        cluster_mass_rank = cluster_rank,
        internal_cluster = s,
        mass = mass[s],
        dimension = d,
        alpha0 = fit$alpha0[s, d],
        beta0 = fit$beta0[s, d],
        prior_mean = fit$alpha0[s, d] / fit$beta0[s, d],
        gamma_true_F1 = gamma_true_order[1],
        gamma_true_F2 = gamma_true_order[2],
        gamma_true_F3 = gamma_true_order[3]
      )
    }))
  }))
}

cluster_parameters = do.call(rbind, lapply(names(base_fits), function(pi_method) {
  rbind(
    cluster_parameter_table(base_fits[[pi_method]], pi_method, "unmerged"),
    cluster_parameter_table(
      merged_fits[[pi_method]], pi_method, "merged_and_refit"
    )
  )
}))
row.names(cluster_parameters) = NULL

true_composition = unlist(lapply(seq_along(TRUE_SUPPORTS), function(g) {
  support = TRUE_SUPPORTS[[g]]
  rows = preliminary$dat$grp == g
  preliminary$dat$L[rows, support[1]] /
    rowSums(preliminary$dat$L[rows, support, drop = FALSE])
}))
geometry = data.frame(
  N = N,
  M = M,
  empirical_active_loading_mean = mean(
    preliminary$dat$L[preliminary$dat$L > 0]
  ),
  composition_q05 = unname(quantile(true_composition, 0.05)),
  composition_median = median(true_composition),
  composition_q95 = unname(quantile(true_composition, 0.95)),
  nmf_mean_factor_cosine = mean(
    match_factor_rows(preliminary$nmf$F, preliminary$dat$F0)$table$cosine
  )
)

write.csv(results, file.path(cache_dir, paste0(cache_tag, "-results.csv")),
          row.names = FALSE)
write.csv(merge_table,
          file.path(cache_dir, paste0(cache_tag, "-merge-choice.csv")),
          row.names = FALSE)
write.csv(cluster_parameters,
          file.path(cache_dir, paste0(cache_tag, "-cluster-parameters.csv")),
          row.names = FALSE)
write.csv(geometry, file.path(cache_dir, paste0(cache_tag, "-geometry.csv")),
          row.names = FALSE)

print(geometry, digits = 7, row.names = FALSE)
print(merge_table, digits = 7, row.names = FALSE)
print(results, digits = 7, row.names = FALSE)
