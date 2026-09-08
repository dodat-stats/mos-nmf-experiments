## Initialization methods for MiSo.
##
## This file contains both the loading-score initializations used by the
## original public miso() wrapper and the truncated-PIP initialization based
## on vanilla NMF followed by fixed-dictionary row-wise Poisson SuSiE.

init_motifs_from_loading_scores <- function(scores, S, D, init_seed = NULL,
                                            min_share = 0.10,
                                            allow_repeats = FALSE,
                                            eps = 1e-12) {
  if (!is.null(init_seed)) set.seed(init_seed)
  K = ncol(scores)
  if (!allow_repeats && D > K) {
    stop("D must be no greater than K for distinct initialization.")
  }
  scores_norm = scores / pmax(rowSums(scores), eps)
  km = kmeans(scores_norm, centers = S, nstart = 20, iter.max = 100)
  motifs = matrix(NA_integer_, nrow = S, ncol = D)

  for (s in seq_len(S)) {
    center = km$centers[s, ] / pmax(sum(km$centers[s, ]), eps)
    ord = order(center, decreasing = TRUE)
    active = ord[center[ord] >= min_share]
    if (length(active) == 0) active = ord[1]
    if (length(active) > D) active = active[seq_len(D)]
    if (allow_repeats) {
      motifs[s, ] = rep(active, length.out = D)
    } else {
      motifs[s, ] = ord[seq_len(D)]
    }
  }

  list(motifs = motifs, cluster = km$cluster, centers = km$centers)
}

gamma_bar_from_motifs <- function(motifs, K, gamma_floor = 0.05) {
  S = nrow(motifs)
  D = ncol(motifs)
  gamma_bar = array(gamma_floor / K, dim = c(S, D, K))
  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      gamma_bar[s, d, motifs[s, d]] =
        1 - gamma_floor + gamma_floor / K
    }
  }
  gamma_bar
}

.miso_all_permutations <- function(x) {
  if (length(x) == 1) return(matrix(x, nrow = 1))
  do.call(rbind, lapply(seq_along(x), function(j) {
    cbind(x[j], .miso_all_permutations(x[-j]))
  }))
}

.miso_greedy_slot_assignment <- function(score) {
  D = nrow(score)
  assignment = rep(NA_integer_, D)
  available = seq_len(D)
  row_order = order(apply(score, 1, max), decreasing = TRUE)
  for (d in row_order) {
    selected = available[which.max(score[d, available])]
    assignment[d] = selected
    available = setdiff(available, selected)
  }
  assignment
}

.miso_best_slot_assignment <- function(score, permutations = NULL) {
  if (nrow(score) != ncol(score)) {
    stop("The slot-alignment score matrix must be square.")
  }
  D = nrow(score)
  if (D == 1) return(1L)

  if (requireNamespace("clue", quietly = TRUE)) {
    return(as.integer(clue::solve_LSAP(score, maximum = TRUE)))
  }

  if (!is.null(permutations)) {
    objective = apply(permutations, 1, function(candidate) {
      sum(score[cbind(seq_len(D), candidate)])
    })
    return(as.integer(permutations[which.max(objective), ]))
  }

  .miso_greedy_slot_assignment(score)
}

.miso_supported_slot_assignment <- function(score, slot_loading,
                                             permutations = NULL) {
  ## Rows of score are supported cluster anchors and columns are row-wise
  ## Poisson-SuSiE slots. Match the supported anchors first. Any remaining
  ## local slots are assigned to uniform tail dimensions in decreasing order
  ## of posterior expected loading.
  n_supported = nrow(score)
  D = ncol(score)
  if (length(slot_loading) != D) {
    stop("slot_loading must contain one value for each Poisson-SuSiE slot.")
  }

  if (n_supported == 0) {
    return(as.integer(order(slot_loading, decreasing = TRUE)))
  }

  if (requireNamespace("clue", quietly = TRUE)) {
    supported_assignment = as.integer(
      clue::solve_LSAP(score, maximum = TRUE)
    )
  } else if (!is.null(permutations)) {
    objective = apply(permutations, 1, function(candidate) {
      sum(score[cbind(seq_len(n_supported),
                      candidate[seq_len(n_supported)])])
    })
    supported_assignment = as.integer(
      permutations[which.max(objective), seq_len(n_supported)]
    )
  } else {
    supported_assignment = rep(NA_integer_, n_supported)
    available = seq_len(D)
    row_order = order(apply(score, 1, max), decreasing = TRUE)
    for (d in row_order) {
      selected = available[which.max(score[d, available])]
      supported_assignment[d] = selected
      available = setdiff(available, selected)
    }
  }

  remaining = setdiff(seq_len(D), supported_assignment)
  remaining = remaining[order(slot_loading[remaining], decreasing = TRUE)]
  as.integer(c(supported_assignment, remaining))
}

.miso_gamma_shape_from_moments <- function(E_lambda, E_log_lambda,
                                           n_iter = 8,
                                           min_shape = 1e-3,
                                           max_shape = 1e4,
                                           eps = 1e-12) {
  delta = log(pmax(E_lambda, eps)) - E_log_lambda
  delta = pmax(delta, eps)
  shape = ifelse(delta > 0.5, 1 / (2 * delta), 1 / delta)
  shape = pmin(pmax(shape, min_shape), max_shape)

  for (j in seq_len(n_iter)) {
    objective = log(shape) - digamma(shape) - delta
    derivative = 1 / shape - trigamma(shape)
    shape = shape - objective / derivative
    shape = pmin(pmax(shape, min_shape), max_shape)
  }
  shape
}

## Construct the initialization described in the manuscript's "My attempt"
## paragraph. Observations are grouped by truncated row-wise score sets. The
## score can be either the ordinary PIP or the posterior expected loading
##
##   R_ik = sum_d gamma_bar_idk * alpha_idk / beta_idk.
##
## Within a group, supported factors become ordered near-point-mass anchors;
## dimensions beyond the truncated support are exactly uniform over factors.
## Gamma priors are initialized by matching cluster-averaged E(lambda) and
## E(log lambda) after aligning the exchangeable Poisson-SuSiE slots.
init_pip_set_miso <- function(ps_fit, pip_truncate_level = 0.95,
                              gamma_floor = 0.05,
                              support_score = c("pip", "expected_loading"),
                              min_shape = 1e-3, max_shape = 1e4,
                              eps = 1e-12) {
  support_score = match.arg(support_score)
  if (pip_truncate_level <= 0 || pip_truncate_level > 1) {
    stop("pip_truncate_level must be strictly positive and at most one.")
  }
  if (gamma_floor < 0 || gamma_floor >= 1) {
    stop("gamma_floor must be nonnegative and smaller than one.")
  }
  if (is.null(ps_fit$gamma_bar) || is.null(ps_fit$alpha)) {
    stop("Score-set initialization requires Poisson-SuSiE gamma_bar and alpha.")
  }
  beta = if (!is.null(ps_fit$beta)) ps_fit$beta else ps_fit$lambda
  if (is.null(beta)) {
    stop("Score-set initialization requires Poisson-SuSiE beta rates.")
  }

  N = dim(ps_fit$gamma_bar)[1]
  D = dim(ps_fit$gamma_bar)[2]
  K = dim(ps_fit$gamma_bar)[3]
  if (!all(dim(ps_fit$alpha) == c(N, D, K)) ||
      !all(dim(beta) == c(N, D, K))) {
    stop("Poisson-SuSiE gamma_bar, alpha, and beta must have the same dimensions.")
  }
  if (D > K) {
    stop("Score-set initialization currently requires D to be no greater than K.")
  }

  expected_lambda = ps_fit$alpha / pmax(beta, eps)
  expected_log_lambda = digamma(ps_fit$alpha) - log(pmax(beta, eps))
  pip = 1 - apply(1 - ps_fit$gamma_bar, c(1, 3), prod)
  R = apply(ps_fit$gamma_bar * expected_lambda, c(1, 3), sum)
  truncation_score = if (support_score == "pip") pip else R

  support_by_observation = vector("list", N)
  support_key = character(N)
  for (i in seq_len(N)) {
    factor_order = order(truncation_score[i, ], decreasing = TRUE)
    cumulative_share = cumsum(truncation_score[i, factor_order]) /
      pmax(sum(truncation_score[i, ]), eps)
    crossing = which(cumulative_share >= pip_truncate_level)
    truncated_size = if (length(crossing) == 0) D else crossing[1]
    truncated_size = min(truncated_size, D)
    support_by_observation[[i]] = sort(factor_order[seq_len(truncated_size)])
    support_key[i] = paste(support_by_observation[[i]], collapse = ",")
  }

  distinct_key = unique(support_key)
  cluster = match(support_key, distinct_key)
  S = length(distinct_key)
  cluster_members = lapply(seq_len(S), function(s) which(cluster == s))
  support_sets = lapply(cluster_members, function(index) {
    support_by_observation[[index[1]]]
  })
  support_size = lengths(support_sets)

  cluster_pip = matrix(0, S, K)
  cluster_R = matrix(0, S, K)
  cluster_score = matrix(0, S, K)
  anchor_factor = matrix(NA_integer_, S, D)
  gamma_bar = array(1 / K, dim = c(S, D, K))
  for (s in seq_len(S)) {
    cluster_pip[s, ] = colMeans(pip[cluster_members[[s]], , drop = FALSE])
    cluster_R[s, ] = colMeans(R[cluster_members[[s]], , drop = FALSE])
    cluster_score[s, ] = colMeans(
      truncation_score[cluster_members[[s]], , drop = FALSE]
    )
    supported = support_sets[[s]]
    supported = supported[
      order(cluster_score[s, supported], decreasing = TRUE)
    ]
    support_sets[[s]] = supported
    for (d in seq_along(supported)) {
      anchor_factor[s, d] = supported[d]
      gamma_bar[s, d, ] = gamma_floor / K
      gamma_bar[s, d, supported[d]] =
        1 - gamma_floor + gamma_floor / K
    }
  }

  slot_mean = rowSums(ps_fit$gamma_bar * expected_lambda, dims = 2)
  slot_mean_log = rowSums(
    ps_fit$gamma_bar * expected_log_lambda, dims = 2
  )

  slot_alignment = matrix(NA_integer_, N, D)
  cluster_mean = matrix(0, S, D)
  cluster_mean_log = matrix(0, S, D)
  permutations = NULL
  if (!requireNamespace("clue", quietly = TRUE) && D <= 8) {
    permutations = .miso_all_permutations(seq_len(D))
  } else if (!requireNamespace("clue", quietly = TRUE) && D > 8) {
    warning(
      "Package 'clue' is unavailable; using greedy supported-slot matching ",
      "because D > 8."
    )
  }

  for (i in seq_len(N)) {
    s = cluster[i]
    n_supported = support_size[s]
    compatibility = matrix(0, nrow = n_supported, ncol = D)
    for (d in seq_len(n_supported)) {
      anchor = anchor_factor[s, d]
      compatibility[d, ] =
        ps_fit$gamma_bar[i, , anchor] * expected_lambda[i, , anchor]
    }
    assignment = .miso_supported_slot_assignment(
      compatibility,
      slot_loading = slot_mean[i, ],
      permutations = permutations
    )
    slot_alignment[i, ] = assignment
    cluster_mean[s, ] = cluster_mean[s, ] + slot_mean[i, assignment]
    cluster_mean_log[s, ] =
      cluster_mean_log[s, ] + slot_mean_log[i, assignment]
  }

  cluster_size = lengths(cluster_members)
  cluster_mean = cluster_mean / cluster_size
  cluster_mean_log = cluster_mean_log / cluster_size
  alpha0 = .miso_gamma_shape_from_moments(
    E_lambda = cluster_mean,
    E_log_lambda = cluster_mean_log,
    min_shape = min_shape,
    max_shape = max_shape,
    eps = eps
  )
  beta0 = alpha0 / pmax(cluster_mean, eps)

  omega_init = matrix(0, N, S)
  omega_init[cbind(seq_len(N), cluster)] = 1
  pi_init = cluster_size / N

  list(
    gamma_bar = gamma_bar,
    alpha0 = alpha0,
    beta0 = beta0,
    omega_init = omega_init,
    pi_init = pi_init,
    pip = pip,
    R = R,
    truncation_score = truncation_score,
    support_score = support_score,
    cluster_pip = cluster_pip,
    cluster_R = cluster_R,
    cluster_score = cluster_score,
    cluster = as.integer(cluster),
    cluster_size = cluster_size,
    support_sets = support_sets,
    support_size = support_size,
    anchor_factor = anchor_factor,
    slot_alignment = slot_alignment,
    cluster_mean_lambda = cluster_mean,
    cluster_mean_log_lambda = cluster_mean_log,
    motif_initialization = if (support_score == "pip") {
      "pip_sets_uniform_tail"
    } else {
      "loading_sets_uniform_tail"
    },
    description = paste(
      if (support_score == "pip") "Truncated PIP-set clusters," else
        "Truncated expected-loading-set clusters,",
      "supported factor anchors, uniform tail",
      "dimensions, and moment-matched Gamma priors"
    )
  )
}

init_loading_set_miso <- function(ps_fit, loading_truncate_level = 0.95,
                                  gamma_floor = 0.05,
                                  min_shape = 1e-3, max_shape = 1e4,
                                  eps = 1e-12) {
  init_pip_set_miso(
    ps_fit = ps_fit,
    pip_truncate_level = loading_truncate_level,
    support_score = "expected_loading",
    gamma_floor = gamma_floor,
    min_shape = min_shape,
    max_shape = max_shape,
    eps = eps
  )
}

## Align each observation's exchangeable MF slots to the distinct seed factors
## of its preliminary cluster, then average the aligned MF factor-selection
## posteriors. Slot contributions are weighted by their share of that
## observation's expected total loading, so observations have equal total
## weight regardless of count depth.
init_aligned_mf_gamma <- function(mf_fit, hard_init, gamma_floor = 0.05,
                                  eps = 1e-12) {
  if (is.null(mf_fit$gamma_bar) || is.null(mf_fit$alpha)) {
    stop("aligned_mf initialization requires MF gamma_bar and alpha.")
  }
  beta = if (!is.null(mf_fit$beta)) mf_fit$beta else mf_fit$lambda
  if (is.null(beta)) {
    stop("aligned_mf initialization requires MF beta rates.")
  }

  N = dim(mf_fit$gamma_bar)[1]
  D = dim(mf_fit$gamma_bar)[2]
  K = dim(mf_fit$gamma_bar)[3]
  S = nrow(hard_init$motifs)
  if (!all(dim(mf_fit$alpha) == c(N, D, K)) ||
      !all(dim(beta) == c(N, D, K))) {
    stop("MF gamma_bar, alpha, and beta must have the same dimensions.")
  }
  if (ncol(hard_init$motifs) != D ||
      length(hard_init$cluster) != N) {
    stop("The MF fit and preliminary motif clustering have incompatible dimensions.")
  }
  if (any(hard_init$motifs < 1 | hard_init$motifs > K)) {
    stop("A distinct seed factor is outside the MF factor range.")
  }

  expected_lambda = mf_fit$alpha / pmax(beta, eps)
  slot_loading = rowSums(mf_fit$gamma_bar * expected_lambda, dims = 2)
  slot_fraction = slot_loading / pmax(rowSums(slot_loading), eps)

  aligned_average = array(0, dim = c(S, D, K))
  aligned_weight = matrix(0, nrow = S, ncol = D)
  slot_alignment = matrix(NA_integer_, nrow = N, ncol = D)

  permutations = NULL
  if (!requireNamespace("clue", quietly = TRUE) && D <= 8) {
    permutations = .miso_all_permutations(seq_len(D))
  } else if (!requireNamespace("clue", quietly = TRUE) && D > 8) {
    warning(
      "Package 'clue' is unavailable; using greedy aligned_mf slot matching ",
      "because D > 8."
    )
  }

  for (i in seq_len(N)) {
    s = hard_init$cluster[i]
    compatibility = matrix(0, nrow = D, ncol = D)
    for (d in seq_len(D)) {
      seed_factor = hard_init$motifs[s, d]
      compatibility[d, ] =
        mf_fit$gamma_bar[i, , seed_factor] *
        expected_lambda[i, , seed_factor]
    }
    assignment = .miso_best_slot_assignment(
      compatibility, permutations = permutations
    )
    slot_alignment[i, ] = assignment

    for (d in seq_len(D)) {
      local_slot = assignment[d]
      weight = slot_fraction[i, local_slot]
      aligned_average[s, d, ] = aligned_average[s, d, ] +
        weight * mf_fit$gamma_bar[i, local_slot, ]
      aligned_weight[s, d] = aligned_weight[s, d] + weight
    }
  }

  gamma_bar = array(0, dim = c(S, D, K))
  for (s in seq_len(S)) {
    for (d in seq_len(D)) {
      if (aligned_weight[s, d] > eps) {
        posterior_average = aligned_average[s, d, ] / aligned_weight[s, d]
        posterior_average = posterior_average / sum(posterior_average)
      } else {
        posterior_average = rep(0, K)
        posterior_average[hard_init$motifs[s, d]] = 1
      }
      aligned_average[s, d, ] = posterior_average
      gamma_bar[s, d, ] =
        (1 - gamma_floor) * posterior_average + gamma_floor / K
    }
  }

  list(
    gamma_bar = gamma_bar,
    aligned_average = aligned_average,
    aligned_weight = aligned_weight,
    slot_alignment = slot_alignment
  )
}

init_soft_gamma_from_loading_scores <- function(
    scores, S, D, init_seed = NULL, min_share = 0.10,
    gamma_floor = 0.05, surplus_slots = c("repeat", "uniform"),
    motif_initialization = c("distinct", "threshold", "aligned_mf"),
    eps = 1e-12, mf_fit = NULL) {
  surplus_slots = match.arg(surplus_slots)
  motif_initialization = match.arg(motif_initialization)
  uses_distinct_seeds = motif_initialization %in% c("distinct", "aligned_mf")
  hard_init = init_motifs_from_loading_scores(
    scores = scores,
    S = S,
    D = D,
    init_seed = init_seed,
    min_share = min_share,
    allow_repeats = !uses_distinct_seeds,
    eps = eps
  )

  if (motif_initialization == "aligned_mf") {
    aligned = init_aligned_mf_gamma(
      mf_fit = mf_fit,
      hard_init = hard_init,
      gamma_floor = gamma_floor,
      eps = eps
    )
    return(c(
      list(
        gamma_bar = aligned$gamma_bar,
        hard_init = hard_init,
        motif_initialization = motif_initialization
      ),
      aligned[c("aligned_average", "aligned_weight", "slot_alignment")]
    ))
  }

  K = ncol(scores)
  gamma_bar = array(gamma_floor / K, dim = c(S, D, K))
  for (s in seq_len(S)) {
    seen = rep(FALSE, K)
    for (d in seq_len(D)) {
      k = hard_init$motifs[s, d]
      is_surplus_repeat = seen[k]
      if (motif_initialization == "distinct" ||
          !(surplus_slots == "uniform" && is_surplus_repeat)) {
        gamma_bar[s, d, k] = 1 - gamma_floor + gamma_floor / K
      }
      seen[k] = TRUE
    }
  }

  list(
    gamma_bar = gamma_bar,
    hard_init = hard_init,
    motif_initialization = motif_initialization
  )
}
