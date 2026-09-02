#!/usr/bin/env Rscript

# Phone-readable summary of the large sci-Space MiSo experiment.
#
# Outputs:
#   output/sci-space-miso-large/phone-report/01-goals-summary.png
#   output/sci-space-miso-large/phone-report/02-factor-graph.png
#   output/sci-space-miso-large/phone-report/03-factor-interpretation-a.png
#   output/sci-space-miso-large/phone-report/04-factor-interpretation-b.png
#   output/sci-space-miso-large/phone-report/factor_interpretation.csv
#   output/sci-space-miso-large/phone-report/factor_graph_edges.csv

options(stringsAsFactors = FALSE)

project_dir <- normalizePath(".", mustWork = TRUE)
result_dir <- file.path(project_dir, "output", "sci-space-miso-large")
output_dir <- file.path(result_dir, "phone-report")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fit_file <- file.path(
  result_dir,
  "fits-v2-n2000-m800-K18-S20-D4",
  "fit-seed-003.rds"
)
fit_object <- readRDS(fit_file)
factor_scores <- fit_object$miso$factor_scores
effective_dimension <- fit_object$miso$active_dimension
mixture_weight <- fit_object$miso$fit$pi

K <- ncol(factor_scores)
S <- nrow(factor_scores)
factor_ids <- paste0("F", seq_len(K))
submanifold_ids <- paste0("S", seq_len(S))
colnames(factor_scores) <- factor_ids
rownames(factor_scores) <- submanifold_ids

# The fitted effective dimension determines which factors belong to each
# submanifold. This avoids introducing a second arbitrary threshold.
active <- matrix(FALSE, nrow = S, ncol = K,
                 dimnames = list(submanifold_ids, factor_ids))
for (s in seq_len(S)) {
  active[s, order(factor_scores[s, ], decreasing = TRUE)[seq_len(effective_dimension[s])]] <- TRUE
}

# Node usage is the mixture-weighted factor share. An edge is added whenever
# two factors are active in the same fitted submanifold. Its weight is the
# submanifold mass times the geometric mean of the two within-submanifold
# shares. This makes common, balanced co-use visually strongest.
factor_usage <- colSums(sweep(factor_scores, 1, mixture_weight, `*`))
edge_rows <- list()
edge_index <- 0L
for (s in seq_len(S)) {
  factors <- which(active[s, ])
  if (length(factors) < 2L) next
  pairs <- combn(factors, 2L)
  for (j in seq_len(ncol(pairs))) {
    edge_index <- edge_index + 1L
    i1 <- pairs[1L, j]
    i2 <- pairs[2L, j]
    edge_rows[[edge_index]] <- data.frame(
      from = factor_ids[i1],
      to = factor_ids[i2],
      submanifold = submanifold_ids[s],
      mixture_weight = mixture_weight[s],
      co_use_weight = mixture_weight[s] * sqrt(factor_scores[s, i1] * factor_scores[s, i2])
    )
  }
}
edges <- do.call(rbind, edge_rows)

factor_notes <- data.frame(
  factor = factor_ids,
  family = c(
    "Mesenchymal / ECM", "Neural", "Mesenchymal / ECM", "Neural", "Neural", "Neural",
    "Shared developmental", "Muscle", "Neural", "Neural", "Muscle", "Mesenchymal / ECM",
    "Blood / liver", "Mesenchymal / ECM", "Neural", "Mesenchymal / ECM", "Blood / liver",
    "Mesenchymal / ECM"
  ),
  label = c(
    "Cartilage matrix", "Radial-glial regulation", "Vascular mesenchyme",
    "Neuronal adhesion", "Neuronal differentiation", "Neuroepithelial differentiation",
    "Broad immature/core", "Skeletal muscle", "Testis-associated developmental",
    "Choroid plexus secretory", "Cardiac muscle", "Mesenchymal cartilage ECM",
    "Fetal liver/hepatocyte", "Collagen-rich stroma", "Glial/radial-glial",
    "Broad mesenchymal stroma", "Erythroid", "Apoe/Ptgds stromal support"
  ),
  confidence = c(
    "high", "moderate", "moderate", "high", "high", "moderate", "moderate", "high", "low",
    "high", "high", "moderate", "high", "moderate", "high", "moderate", "high", "moderate"
  ),
  interpretation = c(
    "Collagen II/IX/XI and aggrecan form a clear chondrocyte cartilage-matrix program.",
    "Nfib, Nfia, Nrg1, Dcc and Tcf4 suggest neural progenitor or radial-glial regulation.",
    "Sparc, Col4a1, Rbms3 and Gpc6 support a vascular or mesenchymal extracellular-matrix program.",
    "Pcdh9, Lsamp, Dcc, Nrxn1, Nlgn1 and Ptprd form a neuronal adhesion and axon-wiring program.",
    "Tuba1a, Map1b and Dpysl3 with neuronal adhesion genes indicate differentiating neurons.",
    "Sox11, Tuba1a, Map1b, Meis2 and Nnat indicate an immature neuroepithelial program.",
    "H19, Ccnd2, Tubb5 and translation genes recur across many lineages; this is shared developmental signal, not a specific cell type.",
    "Ttn, Actc1, Neb, Acta1, Mylpf and Myh3 identify skeletal-muscle contraction.",
    "This factor is concentrated in the testis-associated submanifold, but its genes are broadly developmental; the label is provisional.",
    "Ttr, Trpm3, Htr2c and Igfbp2 identify a choroid-plexus secretory program.",
    "Ttn, Myh7, Actc1, Myl2, Tnnt2 and Myl3 identify cardiac-muscle contraction.",
    "Gpc6, Trps1, Fbn2, Vcan and collagens indicate connective or cartilage-associated mesenchymal matrix.",
    "Afp, Alb and Trf with globins identify fetal liver or hepatocyte signal mixed with erythroid output.",
    "Col1a1, Col1a2, Col3a1 and Sparc define collagen-rich stromal or vascular matrix.",
    "Fabp7, Slc1a3, Ptprz1, Vim, Nfib and Ptn identify glial or radial-glial identity.",
    "Col3a1, Col1a1/2 and Gpc6 form a broad connective-tissue or stromal program reused across mesodermal submanifolds.",
    "Hemoglobin genes with Slc4a1 form a clear erythroid program.",
    "Apoe, Ptgds and Selenop suggest a secretory support-cell program within a fibroblast-rich submanifold; identity remains mixed."
  )
)

factor_gene_file <- file.path(result_dir, "best_seed_factor_genes.csv")
factor_genes <- read.csv(factor_gene_file, check.names = FALSE)
factor_notes$top_genes <- vapply(factor_ids, function(f) {
  paste(head(factor_genes$gene[factor_genes$factor == f], 6L), collapse = ", ")
}, character(1))

factor_notes$usage <- as.numeric(factor_usage[factor_notes$factor])
factor_notes$submanifolds <- vapply(seq_len(K), function(k) {
  candidates <- which(active[, k])
  if (!length(candidates)) return("none")
  ranked <- candidates[order(mixture_weight[candidates] * factor_scores[candidates, k], decreasing = TRUE)]
  paste(submanifold_ids[head(ranked, 4L)], collapse = ", ")
}, character(1))

write.csv(factor_notes, file.path(output_dir, "factor_interpretation.csv"), row.names = FALSE)
write.csv(edges, file.path(output_dir, "factor_graph_edges.csv"), row.names = FALSE)

palette <- c(
  "Neural" = "#0072B2",
  "Mesenchymal / ECM" = "#009E73",
  "Muscle" = "#D55E00",
  "Blood / liver" = "#CC79A7",
  "Shared developmental" = "#666666"
)
status_palette <- c(
  "supported" = "#2E7D32",
  "mixed" = "#C47F00",
  "not supported yet" = "#B3261E"
)
ink <- "#18202A"
muted <- "#5E6875"
light_gray <- "#E6E9ED"
very_light <- "#F6F7F9"
blue <- "#1769AA"
orange <- "#D55E00"

open_phone_png <- function(filename) {
  png(
    file.path(output_dir, filename),
    width = 1800,
    height = 2600,
    res = 200,
    bg = "white",
    type = "quartz"
  )
  par(
    mar = rep(0, 4),
    oma = rep(0, 4),
    xaxs = "i",
    yaxs = "i",
    family = "Helvetica",
    fg = ink
  )
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = NA)
}

draw_wrapped <- function(x, y, label, width, cex = 1, col = ink, font = 1,
                         adj = c(0, 1), line_spacing = 1.1) {
  lines <- strwrap(label, width = width)
  text(x, y, paste(lines, collapse = "\n"), adj = adj, cex = cex,
       col = col, font = font, family = "Helvetica")
  invisible(length(lines) * line_spacing)
}

draw_header <- function(title, subtitle) {
  text(0.055, 0.965, title, adj = c(0, 1), cex = 2.0, font = 2, col = ink)
  draw_wrapped(0.055, 0.920, subtitle, width = 86, cex = 0.92, col = muted)
  segments(0.055, 0.892, 0.945, 0.892, col = light_gray, lwd = 2)
}

# -------------------------------------------------------------------------
# Page 1: goals and headline metrics
# -------------------------------------------------------------------------
open_phone_png("01-goals-summary.png")
draw_header(
  "MiSo on the sci-Space mouse embryo",
  "2,000 cells, 800 genes, K=18 selected by held-out NMF; MiSo fitted with S=20 and D=4 across three seeds."
)

metrics <- data.frame(
  name = c("Held-out deviance", "Lineage NMI", "Active factors per cell", "Spatial excess agreement"),
  direction = c("lower is better", "higher is better", "lower is better", "higher is better"),
  nmf = c(0.3451, 0.560, 4.17, 0.0436),
  miso = c(0.3587, 0.577, 2.13, 0.0446),
  maximum = c(0.40, 0.65, 5.0, 0.055),
  digits = c(4, 3, 2, 4)
)

metric_top <- 0.855
metric_height <- 0.090
for (i in seq_len(nrow(metrics))) {
  top <- metric_top - (i - 1) * metric_height
  text(0.055, top, metrics$name[i], adj = c(0, 1), cex = 1.05, font = 2, col = ink)
  text(0.945, top, metrics$direction[i], adj = c(1, 1), cex = 0.75, col = muted)
  bar_left <- 0.30
  bar_right <- 0.87
  bar_width <- bar_right - bar_left
  nmf_y <- top - 0.031
  miso_y <- top - 0.060
  text(bar_left - 0.018, nmf_y, "NMF", adj = c(1, 0.5), cex = 0.76, col = muted)
  text(bar_left - 0.018, miso_y, "MiSo", adj = c(1, 0.5), cex = 0.76, col = ink, font = 2)
  rect(bar_left, nmf_y - 0.008, bar_right, nmf_y + 0.008, col = very_light, border = NA)
  rect(bar_left, miso_y - 0.008, bar_right, miso_y + 0.008, col = very_light, border = NA)
  rect(bar_left, nmf_y - 0.008,
       bar_left + bar_width * metrics$nmf[i] / metrics$maximum[i], nmf_y + 0.008,
       col = "#9AA3AD", border = NA)
  rect(bar_left, miso_y - 0.008,
       bar_left + bar_width * metrics$miso[i] / metrics$maximum[i], miso_y + 0.008,
       col = blue, border = NA)
  text(0.945, nmf_y, formatC(metrics$nmf[i], format = "f", digits = metrics$digits[i]),
       adj = c(1, 0.5), cex = 0.77, col = muted)
  text(0.945, miso_y, formatC(metrics$miso[i], format = "f", digits = metrics$digits[i]),
       adj = c(1, 0.5), cex = 0.77, col = ink, font = 2)
}

goals <- read.csv(file.path(result_dir, "goal_checks.csv"), check.names = FALSE)
text(0.055, 0.482, "Goal evaluation", adj = c(0, 1), cex = 1.35, font = 2, col = ink)
goal_y <- seq(0.438, 0.105, length.out = nrow(goals))
short_evidence <- c(
  "2.13 active factors for MiSo versus 4.17 for NMF.",
  "Only a 3.9% deviance cost for the structural constraint.",
  "Lineage NMI is 0.577 for MiSo versus 0.560 for NMF.",
  "Positive spatial signal, but essentially tied with NMF.",
  "Assignment uncertainty is still driven by depth, not boundaries.",
  "Gene signatures are stable; hard partitions remain seed-sensitive.",
  "Submanifolds use 1-3 factors; the largest contains only 10.3% of cells."
)
for (i in seq_len(nrow(goals))) {
  y <- goal_y[i]
  points(0.070, y - 0.004, pch = 21, bg = status_palette[goals$status[i]],
         col = status_palette[goals$status[i]], cex = 1.55)
  text(0.095, y + 0.006, goals$goal[i], adj = c(0, 1), cex = 0.88, font = 2, col = ink)
  text(0.945, y + 0.006, goals$status[i], adj = c(1, 1), cex = 0.73,
       col = status_palette[goals$status[i]], font = 2)
  text(0.095, y - 0.019, short_evidence[i], adj = c(0, 1), cex = 0.75, col = muted)
}

rect(0.055, 0.025, 0.945, 0.075, col = "#EEF4FA", border = NA)
draw_wrapped(
  0.072, 0.063,
  "Bottom line: MiSo makes the representation substantially more local and readable while preserving prediction and lineage structure. Spatial recovery is neutral, and uncertainty calibration is the main unresolved goal.",
  width = 98, cex = 0.77, col = ink
)
dev.off()

# -------------------------------------------------------------------------
# Page 2: exact factor co-use graph
# -------------------------------------------------------------------------
open_phone_png("02-factor-graph.png")
draw_header(
  "Factor graph induced by the fitted submanifolds",
  "Two factors are connected when they are jointly active in a MiSo submanifold. Edge labels identify that submanifold; thicker edges represent more mixture mass and more balanced co-use."
)

legend_x <- c(0.065, 0.255, 0.455, 0.625, 0.815)
for (i in seq_along(palette)) {
  points(legend_x[i], 0.868, pch = 21, bg = palette[i], col = palette[i], cex = 1.4)
  text(legend_x[i] + 0.015, 0.868, names(palette)[i], adj = c(0, 0.5), cex = 0.63, col = ink)
}

# Hand-tuned, reproducible layout. It preserves the graph structure without
# importing graph packages and keeps labels readable on a phone.
positions <- data.frame(
  factor = factor_ids,
  x = c(0.53, 0.54, 0.19, 0.48, 0.28, 0.12, 0.37, 0.81, 0.36,
        0.22, 0.22, 0.70, 0.17, 0.67, 0.70, 0.62, 0.79, 0.86),
  y = c(0.285, 0.38, 0.50, 0.70, 0.77, 0.79, 0.55, 0.43, 0.34,
        0.305, 0.71, 0.28, 0.62, 0.62, 0.78, 0.47, 0.55, 0.67)
)
rownames(positions) <- positions$factor

edge_scale <- edges$co_use_weight / max(edges$co_use_weight)
for (i in seq_len(nrow(edges))) {
  x1 <- positions[edges$from[i], "x"]
  y1 <- positions[edges$from[i], "y"]
  x2 <- positions[edges$to[i], "x"]
  y2 <- positions[edges$to[i], "y"]
  segments(x1, y1, x2, y2, col = adjustcolor("#6E7781", alpha.f = 0.62),
           lwd = 1.6 + 7.5 * edge_scale[i])
  mx <- (x1 + x2) / 2
  my <- (y1 + y2) / 2
  dx <- x2 - x1
  dy <- y2 - y1
  norm <- sqrt(dx^2 + dy^2)
  offset <- 0.010
  text(mx - offset * dy / norm, my + offset * dx / norm, edges$submanifold[i],
       cex = 0.58, col = muted, font = 2)
}

usage_scaled <- sqrt(factor_usage / max(factor_usage))
for (i in seq_len(K)) {
  id <- factor_ids[i]
  note <- factor_notes[i, ]
  x <- positions[id, "x"]
  y <- positions[id, "y"]
  radius <- 0.017 + 0.017 * usage_scaled[i]
  symbols(x, y, circles = radius, inches = FALSE, add = TRUE,
          bg = palette[note$family], fg = "white", lwd = 2)
  text(x, y, id, cex = 0.74, font = 2, col = "white")
}

text(0.055, 0.225, "Working factor labels", adj = c(0, 1), cex = 1.08, font = 2, col = ink)
columns <- split(seq_len(K), rep(1:3, each = 6))
column_x <- c(0.055, 0.365, 0.675)
for (j in seq_along(columns)) {
  for (r in seq_along(columns[[j]])) {
    i <- columns[[j]][r]
    y <- 0.194 - (r - 1) * 0.0275
    points(column_x[j], y, pch = 21, bg = palette[factor_notes$family[i]],
           col = palette[factor_notes$family[i]], cex = 1.05)
    text(column_x[j] + 0.013, y,
         paste0(factor_notes$factor[i], "  ", factor_notes$label[i]),
         adj = c(0, 0.5), cex = 0.64, col = ink)
  }
}

rect(0.055, 0.012, 0.945, 0.048, col = "#FFF6E8", border = NA)
text(0.072, 0.030,
     "A tree is possible, but it would delete real loops from S3, S5, and S11. The co-use graph is the more faithful summary.",
     adj = c(0, 0.5), cex = 0.72, col = ink)
dev.off()

# -------------------------------------------------------------------------
# Pages 3 and 4: factor interpretation cards
# -------------------------------------------------------------------------
draw_factor_page <- function(indices, filename, page_title, subtitle) {
  open_phone_png(filename)
  draw_header(page_title, subtitle)
  row_top <- 0.860
  row_height <- 0.090
  for (j in seq_along(indices)) {
    i <- indices[j]
    note <- factor_notes[i, ]
    top <- row_top - (j - 1) * row_height
    bottom <- top - 0.078
    if (j %% 2 == 0) rect(0.055, bottom, 0.945, top + 0.004, col = very_light, border = NA)
    rect(0.055, bottom, 0.067, top + 0.004, col = palette[note$family], border = NA)
    text(0.082, top, note$factor, adj = c(0, 1), cex = 1.13, font = 2, col = ink)
    text(0.145, top, note$label, adj = c(0, 1), cex = 1.02, font = 2, col = ink)
    confidence_col <- if (note$confidence == "high") status_palette["supported"] else if (note$confidence == "low") status_palette["not supported yet"] else status_palette["mixed"]
    text(0.930, top, paste0(note$confidence, " confidence"), adj = c(1, 1),
         cex = 0.69, font = 2, col = confidence_col)
    text(0.082, top - 0.027,
         paste0("Top genes: ", note$top_genes), adj = c(0, 1), cex = 0.72, col = muted)
    text(0.930, top - 0.027,
         paste0("Used by: ", note$submanifolds), adj = c(1, 1), cex = 0.70, col = muted)
    draw_wrapped(0.082, top - 0.049, note$interpretation, width = 105, cex = 0.73, col = ink)
  }
  rect(0.055, 0.022, 0.945, 0.064, col = "#EEF4FA", border = NA)
  draw_wrapped(
    0.072, 0.054,
    "These are working names inferred from top-loading genes and the lineage composition of the fitted submanifolds. They should be validated with embryo-level replication, spatial localization, and marker-set enrichment.",
    width = 102, cex = 0.69, col = ink
  )
  dev.off()
}

draw_factor_page(
  1:9,
  "03-factor-interpretation-a.png",
  "Factor interpretation: F1-F9",
  "The labels emphasize biological programs, not mutually exclusive cell types. A factor can be reused by several submanifolds."
)
draw_factor_page(
  10:18,
  "04-factor-interpretation-b.png",
  "Factor interpretation: F10-F18",
  "Confidence is highest when both the marker genes and the external lineage enrichment agree."
)

cat("Created phone-report outputs in", output_dir, "\n")
