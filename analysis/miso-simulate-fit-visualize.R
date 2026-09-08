## ================================================================
## Simulate, fit, and visualize MiSo in one RStudio-friendly script
## ================================================================
## Run this file from the root of the mos-nmf-experiments project.
## The statistical implementation is loaded from code/, but this file contains
## the complete experimental workflow: simulation, NMF, Poisson-SuSiE,
## MiSo initialization, MiSo fitting, and visualization.

source("code/miso-benchmark-utils.R")

## ----------------------------------------------------------------
## 1. Parameters to change
## ----------------------------------------------------------------

## These are the three main fitting choices requested.
K_FIT = 3L
D_FIT = 2L
INIT_PIP_TRUNCATE_LEVEL = 0.70

## Simulation size. The defaults reproduce the doubled simple-tree experiment.
N = 960L
M = 1000L

## Fixed true simple-tree model.
K_TRUE = 3L
S_TRUE = 2L
D_TRUE = 2L
TRUE_SUPPORTS = list(c(1L, 2L), c(1L, 3L))

## Active L_ik ~ Gamma(shape = 2, rate = 1/9), with mean 18.
TRUE_LOADING_SHAPE = 2
TRUE_LOADING_RATE = 1 / 9

## Sparse Gamma simulation for the factor rows, followed by row normalization.
TRUE_FACTOR_SHAPE = 0.1
TRUE_FACTOR_RATE = 0.01

SIMULATION_SEED = 52L
FITTING_SEED = 3L

NMF_MAX_ITERS = 60L
POISSON_SUSIE_MAX_ITERS = 30L
MISO_MAX_ITERS = 100L
MISO_N_INNER = 3L

OUTPUT_DIR = file.path("analysis", "figures")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  K_FIT >= 1,
  D_FIT >= 1,
  INIT_PIP_TRUNCATE_LEVEL > 0,
  INIT_PIP_TRUNCATE_LEVEL <= 1,
  N %% S_TRUE == 0
)

## ----------------------------------------------------------------
## 2. Simulate Y from the simple-tree model
## ----------------------------------------------------------------

set.seed(SIMULATION_SEED)

true_group = rep(seq_len(S_TRUE), each = N / S_TRUE)

F_true = matrix(
  rgamma(
    K_TRUE * M,
    shape = TRUE_FACTOR_SHAPE,
    rate = TRUE_FACTOR_RATE
  ),
  nrow = K_TRUE,
  ncol = M
)
F_true = F_true / rowSums(F_true)

L_true = matrix(0, nrow = N, ncol = K_TRUE)
for (s in seq_len(S_TRUE)) {
  rows = which(true_group == s)
  support = TRUE_SUPPORTS[[s]]
  L_true[rows, support] = matrix(
    rgamma(
      length(rows) * length(support),
      shape = TRUE_LOADING_SHAPE,
      rate = TRUE_LOADING_RATE
    ),
    nrow = length(rows),
    ncol = length(support)
  )
}

true_mean = L_true %*% F_true
Y = matrix(rpois(N * M, lambda = as.vector(true_mean)), nrow = N, ncol = M)

cat(
  "Simulated Y with dimension", N, "x", M,
  "from supports {1,2} and {1,3}.\n"
)

## ----------------------------------------------------------------
## 3. Fit vanilla Poisson NMF with the user-selected K_FIT
## ----------------------------------------------------------------

nmf_fit = poisson_nmf_init(
  Y = Y,
  K = K_FIT,
  max_iters = NMF_MAX_ITERS,
  init_seed = FITTING_SEED
)

## ----------------------------------------------------------------
## 4. Fit row-wise Poisson-SuSiE with F fixed at the NMF estimate
## ----------------------------------------------------------------

poisson_susie_fit = poisson_susie_nmf_fixed_F(
  Y = Y,
  F = nmf_fit$F,
  D = D_FIT,
  max_iters = POISSON_SUSIE_MAX_ITERS,
  update_prior = TRUE,
  init_seed = FITTING_SEED,
  init_L = nmf_fit$L,
  gamma_init_floor = 1e-4,
  elbo_every = 5,
  tol = 5e-4,
  min_iters = 15,
  patience = 2
)

## ----------------------------------------------------------------
## 5. Construct the PIP-set initialization
## ----------------------------------------------------------------

miso_initialization = init_pip_set_miso(
  ps_fit = poisson_susie_fit,
  pip_truncate_level = INIT_PIP_TRUNCATE_LEVEL,
  support_score = "pip",
  gamma_floor = 0.05
)

S_INITIAL = length(miso_initialization$cluster_size)
cat(
  "PIP truncation level", INIT_PIP_TRUNCATE_LEVEL,
  "produced S* =", S_INITIAL,
  "with initial cluster sizes",
  paste(sort(miso_initialization$cluster_size, decreasing = TRUE),
        collapse = ", "),
  "\n"
)

## ----------------------------------------------------------------
## 6. Fit MiSo, including updates of F, gamma_bar, Gamma priors,
##    responsibilities, and mixing proportions
## ----------------------------------------------------------------

miso_fit = miso_fixed_gamma(
  Y = Y,
  F = nmf_fit$F,
  gamma_bar = miso_initialization$gamma_bar,
  max_iters = MISO_MAX_ITERS,
  n_inner = MISO_N_INNER,
  update_prior = TRUE,
  alpha0_init = miso_initialization$alpha0,
  beta0_init = miso_initialization$beta0,
  pi_init = miso_initialization$pi_init,
  omega_init = miso_initialization$omega_init,
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

cat(
  "MiSo finished after", miso_fit$n_iter,
  "iterations; converged =", miso_fit$converged, "\n"
)

## ----------------------------------------------------------------
## 7. Quantities used in the plots
## ----------------------------------------------------------------

row_cosine_local = function(A, B) {
  (A %*% t(B)) /
    (sqrt(rowSums(A^2)) %o% sqrt(rowSums(B^2)))
}

## Find a one-to-one matching between learned and true factors for display.
## If K_FIT > K_TRUE, unmatched fitted factors are shown after F1, F2, F3.
match_factors_for_display = function(F_hat, F_true) {
  similarity = row_cosine_local(F_hat, F_true)
  K_hat = nrow(F_hat)
  K_truth = nrow(F_true)
  n = max(K_hat, K_truth)
  padded_similarity = matrix(0, nrow = n, ncol = n)
  padded_similarity[seq_len(K_hat), seq_len(K_truth)] = similarity

  learned_to_true = rep(NA_integer_, K_hat)
  if (requireNamespace("clue", quietly = TRUE)) {
    assignment = as.integer(clue::solve_LSAP(
      padded_similarity,
      maximum = TRUE
    ))
    learned_to_true[assignment[seq_len(K_hat)] <= K_truth] =
      assignment[seq_len(K_hat)][assignment[seq_len(K_hat)] <= K_truth]
  } else {
    ## Greedy fallback used only to order the visualization.
    remaining_learned = seq_len(K_hat)
    remaining_true = seq_len(K_truth)
    while (length(remaining_learned) > 0 && length(remaining_true) > 0) {
      candidate = which(
        similarity[remaining_learned, remaining_true, drop = FALSE] ==
          max(similarity[remaining_learned, remaining_true, drop = FALSE]),
        arr.ind = TRUE
      )[1, ]
      learned_k = remaining_learned[candidate[1]]
      true_k = remaining_true[candidate[2]]
      learned_to_true[learned_k] = true_k
      remaining_learned = setdiff(remaining_learned, learned_k)
      remaining_true = setdiff(remaining_true, true_k)
    }
  }

  matched_order = unlist(lapply(seq_len(K_truth), function(true_k) {
    which(learned_to_true == true_k)
  }))
  unmatched = setdiff(seq_len(K_hat), matched_order)
  display_order = c(matched_order, unmatched)

  labels_by_learned = paste0("extra ", seq_len(K_hat))
  matched = which(!is.na(learned_to_true))
  labels_by_learned[matched] = paste0("F", learned_to_true[matched])

  list(
    display_order = display_order,
    display_labels = labels_by_learned[display_order],
    learned_to_true = learned_to_true,
    similarity = similarity
  )
}

factor_display = match_factors_for_display(miso_fit$F, F_true)
factor_order = factor_display$display_order
factor_labels = factor_display$display_labels

## The first three colors are fixed so F1, F2, and F3 remain consistent.
true_factor_colors = c("#0072B2", "#E69F00", "#009E73")
extra_factor_colors = if (K_FIT > K_TRUE) {
  hcl.colors(K_FIT - K_TRUE, palette = "Dark 3")
} else {
  character(0)
}
colors_by_learned = hcl.colors(K_FIT, palette = "Dark 3")
for (learned_k in seq_len(K_FIT)) {
  true_k = factor_display$learned_to_true[learned_k]
  if (!is.na(true_k)) colors_by_learned[learned_k] = true_factor_colors[true_k]
}
if (length(extra_factor_colors) > 0) {
  colors_by_learned[is.na(factor_display$learned_to_true)] =
    extra_factor_colors
}
factor_colors = colors_by_learned[factor_order]

normalize_composition = function(scores, eps = 1e-12) {
  scores = pmax(scores, eps)
  scores / rowSums(scores)
}

nmf_composition = normalize_composition(nmf_fit$L[, factor_order, drop = FALSE])
miso_scores = miso_observation_factor_scores(miso_fit)
miso_composition = normalize_composition(
  miso_scores[, factor_order, drop = FALSE]
)

hard_cluster = max.col(miso_fit$omega)
cluster_mass = miso_fit$pi
cluster_order = order(cluster_mass, decreasing = TRUE)
cluster_rank = match(hard_cluster, cluster_order)

## S x K inclusion probabilities. Under the variational factorization,
## this is the probability that factor k is selected by at least one of the
## D_FIT submanifold dimensions.
factor_inclusion_probability = matrix(
  0,
  nrow = ncol(miso_fit$omega),
  ncol = K_FIT
)
for (s in seq_len(nrow(factor_inclusion_probability))) {
  for (k in seq_len(K_FIT)) {
    factor_inclusion_probability[s, k] =
      1 - prod(1 - miso_fit$gamma_bar[s, , k])
  }
}

factor_inclusion_display = factor_inclusion_probability[
  cluster_order,
  factor_order,
  drop = FALSE
]

cat("\nS x K factor-inclusion probabilities:\n")
printed_support = round(factor_inclusion_display, 4)
rownames(printed_support) = paste0("cluster_", seq_along(cluster_order))
colnames(printed_support) = factor_labels
print(printed_support)

## ----------------------------------------------------------------
## 8. Figure 1: ordinary NMF versus adjacent MiSo clusters
## ----------------------------------------------------------------

cluster_colors = hcl.colors(max(length(cluster_order), 3), "Dark 3")[
  seq_along(cluster_order)
]

plot_loading_comparison = function() {
  layout(matrix(seq_len(2), ncol = 1), heights = c(1, 1))
  compact_device = dev.size("in")[2] < 7 || dev.size("in")[1] < 8
  par(
    mar = if (compact_device) {
      c(2.2, 3.3, 2.3, 0.5)
    } else {
      c(3.5, 4.5, 4.2, 1.0)
    },
    oma = if (compact_device) {
      c(0.2, 0.2, 2.5, 0.2)
    } else {
      c(0.5, 0.5, 4.2, 0.5)
    },
    family = "sans"
  )

  set.seed(2026)
  arbitrary_order = sample(seq_len(N))
  barplot(
    t(nmf_composition[arbitrary_order, , drop = FALSE]),
    col = factor_colors,
    border = NA,
    space = 0,
    axes = FALSE,
    ylim = c(0, 1)
  )
  axis(2, at = c(0, 0.5, 1), las = 1, cex.axis = 0.85)
  box(col = "#777777", lwd = 0.8)
  title(
    main = sprintf("A. Standard NMF with fitted K = %d", K_FIT),
    adj = 0,
    font.main = 2,
    cex.main = 1.05
  )
  mtext("Factor composition", side = 2,
        line = if (compact_device) 2.1 else 2.8, cex = 0.85)
  mtext("Observations in arbitrary order", side = 1,
        line = if (compact_device) 1.0 else 1.6,
        cex = 0.85)
  legend(
    "topright",
    legend = factor_labels,
    fill = factor_colors,
    border = NA,
    ncol = min(K_FIT, 6),
    bty = "n",
    inset = c(0, -0.20),
    xpd = NA,
    cex = if (K_FIT <= 8) 0.85 else 0.65
  )

  clustered_rows = lapply(seq_along(cluster_order), function(rank) {
    internal_cluster = cluster_order[rank]
    rows = which(hard_cluster == internal_cluster)
    strongest_factor = which.max(factor_inclusion_display[rank, ])
    rows[order(miso_composition[rows, strongest_factor])]
  })
  clustered_order = unlist(clustered_rows)
  cluster_sizes = lengths(clustered_rows)
  cluster_starts = c(1, head(cumsum(cluster_sizes), -1) + 1)
  cluster_ends = cumsum(cluster_sizes)
  cluster_space = rep(0, length(clustered_order))
  if (length(cluster_starts) > 1) cluster_space[cluster_starts[-1]] = 8

  bar_midpoints = barplot(
    t(miso_composition[clustered_order, , drop = FALSE]),
    col = factor_colors,
    border = NA,
    space = cluster_space,
    axes = FALSE,
    ylim = c(0, 1.16)
  )
  axis(2, at = c(0, 0.5, 1), las = 1, cex.axis = 0.85)
  box(col = "#777777", lwd = 0.8)
  title(
    main = sprintf(
      "B. MiSo: S* = %d learned clusters placed side by side",
      length(cluster_order)
    ),
    adj = 0,
    font.main = 2,
    cex.main = 1.05
  )
  mtext("Factor composition", side = 2,
        line = if (compact_device) 2.1 else 2.8, cex = 0.85)
  mtext("Observations (block width equals sample count)", side = 1,
        line = if (compact_device) 1.0 else 1.6, cex = 0.85)

  for (rank in seq_along(cluster_order)) {
    block_indices = cluster_starts[rank]:cluster_ends[rank]
    block_midpoint = mean(bar_midpoints[block_indices])
    block_start = min(bar_midpoints[block_indices])
    block_end = max(bar_midpoints[block_indices])
    segments(
      block_start,
      1.035,
      block_end,
      1.035,
      col = cluster_colors[rank],
      lwd = 4,
      xpd = NA
    )
    relative_width = cluster_sizes[rank] / N
    cluster_label = if (relative_width >= 0.08) {
      sprintf("Cluster %d (n=%d)", rank, cluster_sizes[rank])
    } else {
      sprintf("C%d n=%d", rank, cluster_sizes[rank])
    }
    text(
      block_midpoint,
      1.09,
      cluster_label,
      cex = if (relative_width >= 0.08) 0.80 else 0.66,
      font = 2,
      col = cluster_colors[rank],
      xpd = NA
    )
  }

  mtext(
    sprintf(
      "Simple tree: fitted K=%d, D=%d, PIP truncation=%.2f",
      K_FIT,
      D_FIT,
      INIT_PIP_TRUNCATE_LEVEL
    ),
    side = 3,
    outer = TRUE,
    line = if (compact_device) 1.2 else 2.1,
    font = 2,
    cex = 1.2
  )
}

## ----------------------------------------------------------------
## 9. Figure 2: S x K supported-factor heatmap
## ----------------------------------------------------------------

plot_supported_factor_heatmap = function() {
  ## Reset the two-row layout left behind by plot_loading_comparison().
  ## Without this reset, RStudio tries to draw the heatmap in half a plot pane.
  layout(matrix(1))
  S_FIT = nrow(factor_inclusion_display)
  K_DISPLAY = ncol(factor_inclusion_display)
  heat_colors = colorRampPalette(c("#F7FBFF", "#6A51A3"))(101)
  compact_device = dev.size("in")[2] < 6 || dev.size("in")[1] < 7

  par(
    mar = if (compact_device) {
      c(4.0, 6.3, 3.8, 3.0)
    } else {
      c(6.5, 8.5, 5.5, 4.5)
    },
    family = "sans"
  )
  plot(
    NA,
    xlim = c(0.5, K_DISPLAY + 2.1),
    ylim = c(S_FIT + 0.5, 0.5),
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i",
    main = "Factors supported by each fitted submanifold",
    cex.main = 1.2,
    font.main = 2
  )

  for (s in seq_len(S_FIT)) {
    for (k in seq_len(K_DISPLAY)) {
      probability = factor_inclusion_display[s, k]
      rect(
        k - 0.5,
        s - 0.5,
        k + 0.5,
        s + 0.5,
        col = heat_colors[1 + round(100 * probability)],
        border = if (probability >= 0.5) "#222222" else "white",
        lwd = if (probability >= 0.5) 1.8 else 1
      )
      if (K_DISPLAY <= 12 && S_FIT <= 12) {
        text(
          k,
          s,
          sprintf("%.2f", probability),
          col = if (probability > 0.58) "white" else "#222222",
          cex = 0.85,
          font = if (probability >= 0.5) 2 else 1
        )
      }
    }
  }

  axis(1, at = seq_len(K_DISPLAY), labels = factor_labels, las = 2,
       tick = FALSE, line = -0.5, cex.axis = 0.9)
  cluster_axis_labels = sprintf(
    "Cluster %d  (pi=%.3f)",
    seq_len(S_FIT),
    cluster_mass[cluster_order]
  )
  axis(2, at = seq_len(S_FIT), labels = cluster_axis_labels, las = 1,
       tick = FALSE, cex.axis = 0.85)
  mtext("Fitted factor", side = 1,
        line = if (compact_device) 2.8 else 4.5, cex = 0.9)
  mtext(
    "Cell value: probability that factor k appears in at least one slot",
    side = 3,
    line = 1.0,
    cex = 0.82
  )

  ## Probability color scale.
  bar_left = K_DISPLAY + 0.75
  bar_right = K_DISPLAY + 1.05
  bar_top = 0.75
  bar_bottom = S_FIT + 0.25
  for (j in 0:99) {
    y_bottom = bar_bottom - j / 100 * (bar_bottom - bar_top)
    y_top = bar_bottom - (j + 1) / 100 * (bar_bottom - bar_top)
    rect(
      bar_left,
      y_bottom,
      bar_right,
      y_top,
      col = heat_colors[j + 1],
      border = NA
    )
  }
  rect(bar_left, bar_top, bar_right, bar_bottom, border = "#555555")
  text(bar_right + 0.12, bar_top, "1", pos = 4, cex = 0.8)
  text(bar_right + 0.12, (bar_top + bar_bottom) / 2, "0.5", pos = 4,
       cex = 0.8)
  text(bar_right + 0.12, bar_bottom, "0", pos = 4, cex = 0.8)
}

## Save each figure on its own graphics device. Nothing is drawn on the RStudio
## graphics device, so an invalid RStudio graphics state cannot affect output.
parameter_tag = sprintf(
  "K%d-D%d-tau%03d",
  K_FIT,
  D_FIT,
  round(100 * INIT_PIP_TRUNCATE_LEVEL)
)

save_png = function(filename, width, height, plot_function) {
  output_file = file.path(OUTPUT_DIR, filename)
  png(
    output_file,
    width = width,
    height = height,
    res = 150,
    bg = "white"
  )
  png_device = dev.cur()
  tryCatch(
    plot_function(),
    finally = {
      if (dev.cur() == png_device) dev.off()
    }
  )
  normalizePath(output_file)
}

loading_figure_file = save_png(
  filename = paste0(
    "miso-simulate-fit-loading-clusters-", parameter_tag, ".png"
  ),
  width = 1080,
  height = 1100,
  plot_function = plot_loading_comparison
)

support_figure_file = save_png(
  filename = paste0("miso-supported-factors-", parameter_tag, ".png"),
  width = max(900, 100 * K_FIT + 450),
  height = max(650, 105 * ncol(miso_fit$omega) + 300),
  plot_function = plot_supported_factor_heatmap
)

cat("\nSaved figures:\n")
cat(loading_figure_file, "\n")
cat(support_figure_file, "\n")
