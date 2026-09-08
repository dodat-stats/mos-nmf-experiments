## Merge the small third cluster from the Gamma(2, 1/9) PIP-0.70 fit into
## either dominant cluster, refit S = 2 MiSo, and compare ELBO and recovery.

invisible(capture.output(source(
  "analysis/miso-simple-tree-gamma2-experiment.R"
)))

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
    merged_variance = pmax(
      merged_second_moment - merged_mean^2,
      eps
    )
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
    omega = omega,
    pi_init = pi_init,
    retained_original_clusters = retained_cluster,
    target_new = target_new,
    source_slot_order = source_slot_order
  )
}

fit_merged_candidate = function(base_fit, pi_method, target_cluster,
                                 source_cluster) {
  initialization = merge_miso_clusters(
    fit = base_fit,
    target_cluster = target_cluster,
    source_cluster = source_cluster
  )
  cache_file = file.path(
    cache_dir,
    sprintf(
      paste0(
        "simple-tree-gamma2-rate1over9-%s-pip070-",
        "merge-source%d-target%d-iter%d.rds"
      ),
      pi_method,
      source_cluster,
      target_cluster,
      MISO_MAX_ITERS
    )
  )
  if (file.exists(cache_file)) return(readRDS(cache_file))

  message(
    "Refitting ", pi_method,
    " after merging original cluster ", source_cluster,
    " into original cluster ", target_cluster
  )
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
    omega_init = initialization$omega,
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
  saveRDS(saved, cache_file)
  saved
}

base_fits = lapply(c("point", "dirichlet1"), function(pi_method) {
  readRDS(file.path(
    cache_dir,
    sprintf(
      "simple-tree-gamma2-rate1over9-%s-pip-tau070-seed3-iter100.rds",
      pi_method
    )
  ))$fit
})
names(base_fits) = c("point", "dirichlet1")

base_cluster_order = lapply(names(base_fits), function(pi_method) {
  fit = base_fits[[pi_method]]
  mass = if (pi_method == "point") fit$pi else fit$allocation_mass
  order(mass, decreasing = TRUE)
})
names(base_cluster_order) = names(base_fits)

merged_fits = list()
for (pi_method in names(base_fits)) {
  source_cluster = base_cluster_order[[pi_method]][3]
  for (target_rank in 1:2) {
    target_cluster = base_cluster_order[[pi_method]][target_rank]
    key = paste(pi_method, target_rank, sep = "_")
    merged_fits[[key]] = fit_merged_candidate(
      base_fit = base_fits[[pi_method]],
      pi_method = pi_method,
      target_cluster = target_cluster,
      source_cluster = source_cluster
    )
  }
}

summarize_for_merge = function(fit, pi_method, model, merge_target_rank) {
  result = summarize_candidate(
    score_name = "pip",
    truncate_level = 0.70,
    pi_method = pi_method,
    candidate = list(fit = fit)
  )
  data.frame(
    model = model,
    merge_target_rank = merge_target_rank,
    result,
    stringsAsFactors = FALSE
  )
}

comparison = do.call(rbind, lapply(names(base_fits), function(pi_method) {
  base_summary = summarize_for_merge(
    fit = base_fits[[pi_method]],
    pi_method = pi_method,
    model = "unmerged S=3",
    merge_target_rank = NA_integer_
  )
  merge_summaries = do.call(rbind, lapply(1:2, function(target_rank) {
    summarize_for_merge(
      fit = merged_fits[[paste(pi_method, target_rank, sep = "_")]]$fit,
      pi_method = pi_method,
      model = if (target_rank == 2) {
        "merged S=2: small F2 cluster into F1+F2 motif"
      } else {
        "merged S=2: small F2 cluster into F1+F3 motif"
      },
      merge_target_rank = target_rank
    )
  }))
  rbind(base_summary, merge_summaries)
}))
row.names(comparison) = NULL

comparison$elbo_change_from_unmerged = NA_real_
for (pi_method in unique(comparison$pi_method)) {
  rows = comparison$pi_method == pi_method
  base_elbo = comparison$final_elbo[
    rows & comparison$model == "unmerged S=3"
  ]
  comparison$elbo_change_from_unmerged[rows] =
    comparison$final_elbo[rows] - base_elbo
}

write.csv(
  comparison,
  file.path(cache_dir, "simple-tree-gamma2-merge-small-cluster-results.csv"),
  row.names = FALSE
)

cluster_parameter_table = function(fit, mass, pi_method, stage) {
  factor_match = match_factor_rows(fit$F, preliminary$dat$F0)
  cluster_order = order(mass, decreasing = TRUE)
  do.call(rbind, lapply(seq_along(cluster_order), function(cluster_rank) {
    s = cluster_order[cluster_rank]
    do.call(rbind, lapply(seq_len(ncol(fit$alpha0)), function(d) {
      gamma_true_order = numeric(nrow(fit$F))
      for (learned_k in seq_len(nrow(fit$F))) {
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
        gamma_learned_F1 = fit$gamma_bar[s, d, 1],
        gamma_learned_F2 = fit$gamma_bar[s, d, 2],
        gamma_learned_F3 = fit$gamma_bar[s, d, 3],
        gamma_true_F1 = gamma_true_order[1],
        gamma_true_F2 = gamma_true_order[2],
        gamma_true_F3 = gamma_true_order[3]
      )
    }))
  }))
}

point_correct_merge = merged_fits$point_2
dirichlet_correct_merge = merged_fits$dirichlet1_2
merged_initialization_as_fit = list(
  F = point_correct_merge$initialization$F,
  gamma_bar = point_correct_merge$initialization$gamma_bar,
  alpha0 = point_correct_merge$initialization$alpha0,
  beta0 = point_correct_merge$initialization$beta0
)
cluster_parameters = rbind(
  cluster_parameter_table(
    fit = base_fits$point,
    mass = base_fits$point$pi,
    pi_method = "point",
    stage = "unmerged_S3_posterior"
  ),
  cluster_parameter_table(
    fit = merged_initialization_as_fit,
    mass = point_correct_merge$initialization$pi_init,
    pi_method = "point",
    stage = "merged_S2_initialization"
  ),
  cluster_parameter_table(
    fit = point_correct_merge$fit,
    mass = point_correct_merge$fit$pi,
    pi_method = "point",
    stage = "merged_S2_posterior"
  ),
  cluster_parameter_table(
    fit = base_fits$dirichlet1,
    mass = base_fits$dirichlet1$allocation_mass,
    pi_method = "dirichlet1",
    stage = "unmerged_S3_posterior"
  ),
  cluster_parameter_table(
    fit = dirichlet_correct_merge$fit,
    mass = dirichlet_correct_merge$fit$allocation_mass,
    pi_method = "dirichlet1",
    stage = "merged_S2_posterior"
  )
)
row.names(cluster_parameters) = NULL
write.csv(
  cluster_parameters,
  file.path(cache_dir, "simple-tree-gamma2-merge-cluster-parameters.csv"),
  row.names = FALSE
)

print(
  comparison[, c(
    "pi_method",
    "model",
    "initial_S_star",
    "final_elbo",
    "elbo_change_from_unmerged",
    "data_deviance_per_entry",
    "true_mean_deviance_per_entry",
    "relative_true_mean_error",
    "soft_motif_accuracy",
    "exact_support_accuracy",
    "mean_factor_cosine",
    "aggregate_prior_mean_relative_error",
    "allocation_mass_decreasing"
  )],
  digits = 7,
  row.names = FALSE
)
