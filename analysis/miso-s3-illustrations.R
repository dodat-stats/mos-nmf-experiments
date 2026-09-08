## Phone-friendly illustrations of what the S = 3 MiSo fit adds to NMF.
## Uses the doubled simple-tree Gamma(2, 1/9) experiment.

source("code/miso-benchmark-utils.R")

preliminary_file = paste0(
  "analysis/cache/simple-tree-gamma2-N960-M1000-pip070-preliminary.rds"
)
fit_file = paste0(
  "analysis/cache/simple-tree-gamma2-N960-M1000-pip070-",
  "point-unmerged-iter100.rds"
)
if (!file.exists(preliminary_file) || !file.exists(fit_file)) {
  stop(
    "Run analysis/miso-simple-tree-gamma2-large-merge-experiment.R first."
  )
}

preliminary = readRDS(preliminary_file)
fit = readRDS(fit_file)

output_dir = Sys.getenv(
  "MISO_ILLUSTRATION_DIR",
  unset = file.path("analysis", "figures")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

factor_colors = c("#0072B2", "#E69F00", "#009E73")
cluster_colors = c("#7A5195", "#008C95", "#D1495B")
factor_names = c(expression(F[1]), expression(F[2]), expression(F[3]))

align_columns_to_true_factors = function(scores, learned_to_true) {
  aligned = matrix(0, nrow(scores), ncol(scores))
  for (learned_k in seq_len(ncol(scores))) {
    aligned[, learned_to_true[learned_k]] = scores[, learned_k]
  }
  aligned
}

fit_factor_match = match_factor_rows(fit$F, preliminary$dat$F0)
nmf_factor_match = match_factor_rows(preliminary$nmf$F, preliminary$dat$F0)

miso_scores = align_columns_to_true_factors(
  miso_observation_factor_scores(fit),
  fit_factor_match$learned_to_true
)
nmf_scores = align_columns_to_true_factors(
  preliminary$nmf$L,
  nmf_factor_match$learned_to_true
)
miso_composition = normalize_scores(miso_scores)
nmf_composition = normalize_scores(nmf_scores)

hard_cluster = max.col(fit$omega)
cluster_mass = fit$pi
cluster_order = order(cluster_mass, decreasing = TRUE)
cluster_rank = match(hard_cluster, cluster_order)

cluster_motif_label = function(internal_cluster) {
  cluster_mean = colMeans(
    miso_composition[hard_cluster == internal_cluster, , drop = FALSE]
  )
  active = order(cluster_mean, decreasing = TRUE)[1:2]
  if (cluster_mean[active[2]] < 0.08) {
    return(sprintf("mostly F%d", active[1]))
  }
  sprintf("F%d + F%d", sort(active)[1], sort(active)[2])
}

open_png = function(filename, width, height) {
  png(
    file.path(output_dir, filename),
    width = width,
    height = height,
    res = 150,
    bg = "white"
  )
}

draw_membership_bars = function(composition, observation_order, title,
                                 show_legend = FALSE) {
  barplot(
    t(composition[observation_order, , drop = FALSE]),
    col = factor_colors,
    border = NA,
    space = 0,
    axes = FALSE,
    ylim = c(0, 1)
  )
  axis(2, at = c(0, 0.5, 1), las = 1, cex.axis = 0.85)
  box(col = "#777777", lwd = 0.8)
  title(main = title, adj = 0, font.main = 2, cex.main = 1.05)
  mtext("Factor composition", side = 2, line = 2.8, cex = 0.85)
  mtext("Observations", side = 1, line = 1.6, cex = 0.85)
  if (show_legend) {
    legend(
      "topright",
      legend = c("F1", "F2", "F3"),
      fill = factor_colors,
      border = NA,
      horiz = TRUE,
      bty = "n",
      inset = c(0, -0.20),
      xpd = NA,
      cex = 0.9
    )
  }
}

## Figure 1: the same observations as a single NMF stack and as adjacent MiSo
## clusters. The width of each MiSo block is proportional to its hard-assigned
## sample count, so a small cluster remains visibly small.
open_png("miso-s3-from-nmf-to-clusters.png", width = 1080, height = 1100)
layout(matrix(seq_len(2), ncol = 1), heights = c(1, 1))
par(
  mar = c(3.5, 4.5, 4.2, 1.0),
  oma = c(0.5, 0.5, 4.2, 0.5),
  family = "sans"
)
set.seed(2026)
arbitrary_order = sample(seq_len(nrow(nmf_composition)))
draw_membership_bars(
  nmf_composition,
  arbitrary_order,
  "A. Standard NMF: one stack in arbitrary sample order",
  show_legend = TRUE
)

clustered_rows = lapply(seq_along(cluster_order), function(rank) {
  internal_cluster = cluster_order[rank]
  rows = which(hard_cluster == internal_cluster)
  rows[order(miso_composition[rows, 1])]
})
clustered_order = unlist(clustered_rows)
cluster_sizes = lengths(clustered_rows)
cluster_starts = c(1, head(cumsum(cluster_sizes), -1) + 1)
cluster_ends = cumsum(cluster_sizes)
cluster_space = rep(0, length(clustered_order))
cluster_space[cluster_starts[-1]] = 8

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
  main = "B. MiSo: learned clusters placed side by side",
  adj = 0,
  font.main = 2,
  cex.main = 1.05
)
mtext("Factor composition", side = 2, line = 2.8, cex = 0.85)
mtext("Observations (block width equals sample count)", side = 1,
      line = 1.6, cex = 0.85)
for (rank in seq_along(cluster_order)) {
  block_midpoint = mean(
    bar_midpoints[cluster_starts[rank]:cluster_ends[rank]]
  )
  block_start = min(bar_midpoints[cluster_starts[rank]:cluster_ends[rank]])
  block_end = max(bar_midpoints[cluster_starts[rank]:cluster_ends[rank]])
  segments(
    block_start,
    1.035,
    block_end,
    1.035,
    col = cluster_colors[rank],
    lwd = 4,
    xpd = NA
  )
  cluster_label = if (rank < 3) {
    sprintf("Cluster %d (n=%d)", rank, cluster_sizes[rank])
  } else {
    sprintf("C3 (n=%d)", cluster_sizes[rank])
  }
  text(
    block_midpoint,
    1.09,
    cluster_label,
    cex = if (rank == 3) 0.68 else 0.82,
    font = 2,
    col = cluster_colors[rank],
    xpd = NA
  )
}
mtext(
  "MiSo turns an uninterpretable loading stack into a small motif atlas",
  side = 3,
  outer = TRUE,
  line = 2.1,
  font = 2,
  cex = 1.25
)
dev.off()

## Figure 2: each cluster's two submanifold dimensions as factor recipes.
open_png("miso-s3-cluster-recipes.png", width = 1080, height = 1450)
layout(matrix(seq_len(3), ncol = 1))
par(
  mar = c(4.8, 5.5, 4.2, 2.0),
  oma = c(0.5, 0.5, 4.5, 0.5),
  family = "sans"
)
for (rank in seq_along(cluster_order)) {
  s = cluster_order[rank]
  prior_mean = fit$alpha0[s, ] / fit$beta0[s, ]
  gamma_true_order = matrix(0, nrow = ncol(fit$alpha0), ncol = 3)
  for (learned_k in seq_len(3)) {
    gamma_true_order[, fit_factor_match$learned_to_true[learned_k]] =
      fit$gamma_bar[s, , learned_k]
  }
  expected_contribution = prior_mean * gamma_true_order
  mids = barplot(
    t(expected_contribution),
    horiz = TRUE,
    names.arg = c(expression(d[1]), expression(d[2])),
    col = factor_colors,
    border = NA,
    xlim = c(0, 41),
    axes = FALSE,
    las = 1,
    cex.names = 1.0
  )
  axis(1, at = seq(0, 40, 10), cex.axis = 0.85)
  box(col = "#777777", lwd = 0.8)
  abline(v = 18, col = "#555555", lty = 3, lwd = 1.2)
  text(18, max(mids) + 0.55, "true active mean = 18", cex = 0.75,
       pos = 4, offset = 0.15, xpd = NA, col = "#555555")
  for (d in seq_len(2)) {
    label_x = if (prior_mean[d] < 0.2) 1.0 else prior_mean[d] + 0.8
    label = if (prior_mean[d] < 0.2) {
      sprintf("inactive: mean %.3f, diffuse", prior_mean[d])
    } else {
      sprintf("mean %.1f", prior_mean[d])
    }
    text(label_x, mids[d], label, pos = 4, cex = 0.85, xpd = NA)
  }
  title(
    main = sprintf(
      "Cluster %d - %.1f%% mass - %s",
      rank,
      100 * cluster_mass[s],
      cluster_motif_label(s)
    ),
    adj = 0,
    font.main = 2,
    cex.main = 1.05
  )
  mtext("Prior mean x posterior factor probability", side = 1,
        line = 3.0, cex = 0.85)
  if (rank == 1) {
    legend(
      "topright",
      legend = c("F1", "F2", "F3"),
      fill = factor_colors,
      border = NA,
      horiz = TRUE,
      bty = "n",
      inset = c(0, -0.22),
      xpd = NA,
      cex = 0.9
    )
  }
}
mtext(
  "What each MiSo submanifold is made of",
  side = 3,
  outer = TRUE,
  line = 2.2,
  font = 2,
  cex = 1.25
)
dev.off()

## Figure 3: factor-composition simplex, colored by the inferred cluster.
open_png("miso-s3-composition-simplex.png", width = 1080, height = 1080)
par(
  mar = c(2.5, 2.5, 5.2, 2.5),
  family = "sans",
  xpd = NA
)
simplex_x = miso_composition[, 2] + 0.5 * miso_composition[, 3]
simplex_y = sqrt(3) / 2 * miso_composition[, 3]
plot(
  NA,
  xlim = c(-0.08, 1.08),
  ylim = c(-0.08, 0.98),
  axes = FALSE,
  xlab = "",
  ylab = "",
  asp = 1,
  main = "MiSo reveals two motifs - and a small redundant split",
  cex.main = 1.25,
  font.main = 2
)
polygon(
  c(0, 1, 0.5),
  c(0, 0, sqrt(3) / 2),
  border = "#444444",
  lwd = 2,
  col = "#FAFAFA"
)
segments(0, 0, 0.5, sqrt(3) / 2, col = "#DDDDDD", lwd = 7)
segments(0, 0, 1, 0, col = "#DDDDDD", lwd = 7)
text(0, -0.045, expression(F[1]), cex = 1.2, font = 2)
text(1, -0.045, expression(F[2]), cex = 1.2, font = 2)
text(0.5, sqrt(3) / 2 + 0.045, expression(F[3]), cex = 1.2, font = 2)

## Draw the small cluster last so that its localization is visible.
for (rank in seq_along(cluster_order)) {
  draw_rank = c(1, 2, 3)[rank]
  rows = which(cluster_rank == draw_rank)
  points(
    simplex_x[rows],
    simplex_y[rows],
    pch = 16,
    cex = if (draw_rank == 3) 0.75 else 0.55,
    col = adjustcolor(cluster_colors[draw_rank], alpha.f = 0.48)
  )
}
for (rank in seq_along(cluster_order)) {
  rows = cluster_rank == rank
  center = c(mean(simplex_x[rows]), mean(simplex_y[rows]))
  points(center[1], center[2], pch = 21, cex = 2.0, lwd = 2,
         bg = "white", col = cluster_colors[rank])
  text(center[1], center[2], labels = rank, cex = 0.9, font = 2,
       col = cluster_colors[rank])
}

small_rows = cluster_rank == 3
small_center = c(mean(simplex_x[small_rows]), mean(simplex_y[small_rows]))
text(
  small_center[1] + 0.20,
  small_center[2] + 0.05,
  "small F3-specialized split",
  cex = 0.9,
  col = cluster_colors[3]
)
arrows(
  small_center[1] + 0.17,
  small_center[2] + 0.035,
  small_center[1] + 0.025,
  small_center[2] + 0.005,
  length = 0.08,
  lwd = 1.3,
  col = cluster_colors[3]
)
legend(
  "topright",
  legend = sprintf(
    "Cluster %d (%.1f%%)",
    seq_along(cluster_order),
    100 * cluster_mass[cluster_order]
  ),
  pch = 16,
  col = cluster_colors,
  bty = "n",
  cex = 0.9
)
mtext("Posterior expected factor composition", side = 1, line = 0.8,
      cex = 0.9)
dev.off()

message("Wrote illustrations to ", normalizePath(output_dir))
