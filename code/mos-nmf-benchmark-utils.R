## Benchmark utilities for MoS-NMF simulation studies.

if (!exists("mos_nmf")) {
  source("code/mos-nmf.R")
}

row_cosine <- function(A, B) {
  A_norm = sqrt(rowSums(A^2))
  B_norm = sqrt(rowSums(B^2))
  (A %*% t(B)) / (A_norm %o% B_norm)
}

all_permutations <- function(x) {
  if (length(x) == 1) return(matrix(x, nrow = 1))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], all_permutations(x[-i]))
  }))
}

match_factor_rows <- function(F_hat, F_true) {
  K = nrow(F_true)
  similarity = row_cosine(F_hat, F_true)
  perms = all_permutations(seq_len(K))
  scores = apply(perms, 1, function(p) {
    sum(similarity[cbind(seq_len(K), p)])
  })
  learned_to_true = perms[which.max(scores), ]

  list(
    table = data.frame(
      learned_factor = seq_len(K),
      matched_true_factor = learned_to_true,
      cosine = similarity[cbind(seq_len(K), learned_to_true)]
    ),
    learned_to_true = learned_to_true,
    similarity = similarity
  )
}

support_label <- function(support) paste(sort(unique(support)), collapse = ",")

simulate_support_scenario <- function(N, M, K, group_supports, shape_vec,
                                      factor_shape = 0.1,
                                      factor_rate = 0.01, seed = 1) {
  set.seed(seed)
  S = length(group_supports)
  n_per_group = N / S
  grp = rep(seq_len(S), each = n_per_group)

  F0 = matrix(rgamma(M * K, shape = factor_shape, rate = factor_rate),
              nrow = K, ncol = M)
  F0 = normalize_rows(F0)

  L = matrix(0, nrow = N, ncol = K)
  for (s in seq_len(S)) {
    idx = which(grp == s)
    support = group_supports[[s]]
    L[idx, support] = matrix(
      rgamma(length(idx) * length(support), shape = shape_vec[s], rate = 1),
      nrow = length(idx),
      ncol = length(support)
    )
  }

  Y = matrix(rpois(N * M, L %*% F0), nrow = N, ncol = M)

  list(
    Y = Y,
    F0 = F0,
    L = L,
    grp = grp,
    group_supports = group_supports,
    true_support_size = rowSums(L > 0),
    true_support_label = apply(L > 0, 1, function(x) {
      paste(which(x), collapse = ",")
    }),
    true_patterns = vapply(group_supports, support_label, character(1))
  )
}

normalize_scores <- function(scores, eps = 1e-12) {
  scores = pmax(scores, eps)
  scores / rowSums(scores)
}

kmeans_cluster_scores <- function(scores, S, seed = 1) {
  set.seed(seed)
  kmeans(normalize_scores(scores), centers = S, nstart = 20)$cluster
}

map_scores_to_true <- function(scores, learned_to_true) {
  K = length(learned_to_true)
  out = matrix(0, nrow = nrow(scores), ncol = K)
  for (k in seq_len(K)) {
    out[, learned_to_true[k]] = out[, learned_to_true[k]] + scores[, k]
  }
  out
}

support_labels_from_scores <- function(scores_true_space, support_sizes) {
  vapply(seq_len(nrow(scores_true_space)), function(i) {
    active = order(scores_true_space[i, ], decreasing = TRUE)[
      seq_len(support_sizes[i])
    ]
    paste(sort(active), collapse = ",")
  }, character(1))
}

support_accuracy <- function(scores, learned_to_true, dat) {
  scores_true = map_scores_to_true(scores, learned_to_true)
  mean(support_labels_from_scores(scores_true, dat$true_support_size) ==
         dat$true_support_label)
}

gamma_uncertainty_summary <- function(gamma_bar) {
  ent = motif_gamma_entropy(gamma_bar)
  top = apply(gamma_bar, c(1, 2), max)
  data.frame(
    mean_gamma_entropy = mean(ent),
    max_gamma_entropy = max(ent),
    mean_top_gamma = mean(top),
    min_top_gamma = min(top)
  )
}

empty_uncertainty_summary <- function() {
  data.frame(
    mean_gamma_entropy = NA_real_,
    max_gamma_entropy = NA_real_,
    mean_top_gamma = NA_real_,
    min_top_gamma = NA_real_
  )
}

fit_nmf_baseline <- function(dat, K, S, seed, nmf_iters = 60) {
  nmf_fit = poisson_nmf_init(dat$Y, K = K, max_iters = nmf_iters,
                             init_seed = seed)
  factor_match = match_factor_rows(nmf_fit$F, dat$F0)
  scores = nmf_fit$L
  z_hat = kmeans_cluster_scores(scores, S = S, seed = seed)

  cbind(
    data.frame(
      method = "Poisson NMF + kmeans",
      cluster_accuracy = cluster_accuracy(z_hat, dat$grp),
      support_accuracy = support_accuracy(scores, factor_match$learned_to_true,
                                          dat),
      mean_factor_cosine = mean(factor_match$table$cosine),
      min_factor_cosine = min(factor_match$table$cosine),
      mean_max_responsibility = NA_real_
    ),
    empty_uncertainty_summary()
  )
}

fit_mf_baseline <- function(dat, K, S, D, seed, max_iters = 25,
                            nmf_iters = 50) {
  mf_fit = mf_poi_susie(
    Y = dat$Y,
    K = K,
    D = D,
    max_iters = max_iters,
    update_prior = TRUE,
    init_seed = seed,
    init_F = "poisson_nmf",
    nmf_iters = nmf_iters,
    init_gamma_from_nmf = TRUE,
    elbo_every = 5,
    tol = 5e-4,
    min_iters = 12,
    patience = 2
  )

  factor_match = match_factor_rows(mf_fit$F, dat$F0)
  scores = loading_scores_from_mf(mf_fit)
  z_hat = kmeans_cluster_scores(scores, S = S, seed = seed)

  cbind(
    data.frame(
      method = "MF-Poisson-SuSiE + kmeans",
      cluster_accuracy = cluster_accuracy(z_hat, dat$grp),
      support_accuracy = support_accuracy(scores, factor_match$learned_to_true,
                                          dat),
      mean_factor_cosine = mean(factor_match$table$cosine),
      min_factor_cosine = min(factor_match$table$cosine),
      mean_max_responsibility = NA_real_
    ),
    empty_uncertainty_summary()
  )
}

fit_soft_mos_method <- function(dat, K, S, D, seed, max_iters = 12,
                                n_inner = 3, mf_max_iters = 25,
                                mf_nmf_iters = 50, block_size = 100,
                                surplus_slots = "repeat",
                                update_slot_scale = FALSE) {
  fit = mos_nmf(
    Y = dat$Y,
    K = K,
    S = S,
    D = D,
    max_iters = max_iters,
    n_inner = n_inner,
    mf_max_iters = mf_max_iters,
    mf_nmf_iters = mf_nmf_iters,
    init_seed = seed,
    update_prior = TRUE,
    update_F = TRUE,
    update_gamma = TRUE,
    update_slot_scale = update_slot_scale,
    gamma_init_floor = 0.05,
    surplus_slots = surplus_slots,
    gamma_step_init = 0.5,
    gamma_step_ramp = 8,
    F_step_init = 0.2,
    F_step_ramp = 12,
    tol = 1e-5,
    min_iters = 5,
    patience = 2,
    block_size = block_size
  )

  factor_match = match_factor_rows(fit$F, dat$F0)
  scores = fit$omega %*% motif_factor_scores(fit)$scores
  uncertainty = gamma_uncertainty_summary(fit$gamma_bar)

  cbind(
    data.frame(
      method = "Soft MoS-NMF",
      cluster_accuracy = cluster_accuracy(fit$z_hat, dat$grp),
      support_accuracy = support_accuracy(scores, factor_match$learned_to_true,
                                          dat),
      mean_factor_cosine = mean(factor_match$table$cosine),
      min_factor_cosine = min(factor_match$table$cosine),
      mean_max_responsibility = mean(apply(fit$omega, 1, max))
    ),
    uncertainty
  )
}

run_one_benchmark <- function(dat, K, S, D, seed,
                              nmf_iters = 60,
                              mf_max_iters = 25,
                              mf_nmf_iters = 50,
                              mos_max_iters = 12,
                              mos_n_inner = 3,
                              block_size = 100) {
  rbind(
    fit_nmf_baseline(dat, K = K, S = S, seed = seed,
                     nmf_iters = nmf_iters),
    fit_mf_baseline(dat, K = K, S = S, D = D, seed = seed,
                    max_iters = mf_max_iters, nmf_iters = mf_nmf_iters),
    fit_soft_mos_method(dat, K = K, S = S, D = D, seed = seed,
                        max_iters = mos_max_iters, n_inner = mos_n_inner,
                        mf_max_iters = mf_max_iters,
                        mf_nmf_iters = mf_nmf_iters,
                        block_size = block_size)
  )
}

summarize_benchmark_results <- function(results) {
  metrics = c("cluster_accuracy", "support_accuracy", "mean_factor_cosine",
              "min_factor_cosine", "mean_max_responsibility",
              "mean_gamma_entropy")
  mean_or_na = function(x) {
    if (all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
  }
  sd_or_na = function(x) {
    if (sum(!is.na(x)) < 2) return(NA_real_)
    sd(x, na.rm = TRUE)
  }

  do.call(rbind, lapply(split(results, list(results$scenario, results$method),
                              drop = TRUE), function(tab) {
    out = data.frame(
      scenario = unique(tab$scenario),
      method = unique(tab$method),
      n_seed = length(unique(tab$seed))
    )
    for (metric in metrics) {
      out[[paste0(metric, "_mean")]] = mean_or_na(tab[[metric]])
      out[[paste0(metric, "_sd")]] = sd_or_na(tab[[metric]])
    }
    out
  }))
}
