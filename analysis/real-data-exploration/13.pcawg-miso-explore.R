#!/usr/bin/env Rscript

## Simple PCAWG + MiSo exploration file.
##
## In RStudio: open the mos-nmf-experiments project, edit the settings below,
## and click Source. The script fits MiSo and vanilla Poisson NMF to the same
## raw SBS96 mutation counts, saves the fits, and leaves useful objects in the
## Global Environment for View() or plotting.

## ========================================================================
## 1. SETTINGS TO EDIT
## ========================================================================

K = 36                         # Number of shared factors
S = 30                         # Number of submanifolds (motifs)
D = 5                          # Maximum dimension of each submanifold
SEED = 1
REFIT = FALSE                  # TRUE ignores a cached fit and starts again

MIN_MUTATIONS = 500            # Remove very sparse mutation spectra
MAX_MUTATIONS = 50000          # Remove the hypermutated tail in this pilot
EXCLUDED_CANCER_TYPES = "Skin-Melanoma"

MISO_ITERS = 10                # Increase after a first exploratory run
MISO_INITIALIZATION = "distinct" # Top-D factors without replacement
NMF_ITERS = 200                # Multiplicative-update iterations for vanilla NMF
N_FACTORS_TO_SHOW = 12         # Number of factor spectra to plot
N_MOTIFS_TO_SHOW = 4           # Number of continuous motifs to plot
SECOND_FACTOR_SHARE = 0.05     # Minimum share for a two-factor motif plot
PCA_COMPONENTS = 10            # Hellinger-PCA dimensions to inspect

## ========================================================================
## 2. LOAD THE DATA AND THE MiSo CODE
## ========================================================================

if (!file.exists("code/miso-benchmark-utils.R")) {
  stop("Open the project in RStudio and run this file from the project root.")
}
source("code/miso-benchmark-utils.R")

DATA_DIR = "data/pcawg"
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

## Download the observed spectra and published signatures if they are absent.
## The published signatures are used only to label fitted factors afterward.
download_if_missing = function(filename, url) {
  path = file.path(DATA_DIR, filename)
  if (!file.exists(path)) {
    message("Downloading ", filename, " ...")
    utils::download.file(url, path, mode = "wb")
  }
  path
}

PCAWG7 = paste0(
  "https://raw.githubusercontent.com/steverozen/PCAWG7/",
  "v0.1.4-branch/"
)
spectra_file = download_if_missing(
  "WGS_PCAWG.96.csv",
  paste0(PCAWG7, "data-raw/spectra/WGS_PCAWG.96.csv")
)
signature_file = download_if_missing(
  "PCAWG_sigProfiler_SBS_signatures.csv",
  paste0(
    PCAWG7, "data-raw/sig.profiler.signatures/",
    "sigProfiler_SBS_signatures_2019_05_22.csv"
  )
)

## Rows in the downloaded file are mutation contexts; columns are tumours.
spectra = read.csv(spectra_file, check.names = FALSE)
sample_key_all = names(spectra)[-(1:2)]
cancer_type_all = sub("::.*$", "", sample_key_all)
sample_id_all = sub("^.*::", "", sample_key_all)

Y_all = t(as.matrix(spectra[, -(1:2), drop = FALSE]))
storage.mode(Y_all) = "double"
rownames(Y_all) = sample_id_all

mutation_class = spectra[[1]]
trinucleotide = spectra[[2]]
mutation_context = paste0(
  substr(trinucleotide, 1, 1), "[", mutation_class, "]",
  substr(trinucleotide, 3, 3)
)
colnames(Y_all) = mutation_context

## Filter by total mutation burden. Cancer labels do not enter model fitting.
mutation_burden_all = rowSums(Y_all)
keep = mutation_burden_all >= MIN_MUTATIONS &
  mutation_burden_all <= MAX_MUTATIONS &
  !(cancer_type_all %in% EXCLUDED_CANCER_TYPES)

Y = Y_all[keep, , drop = FALSE]
cancer_type = cancer_type_all[keep]
sample_id = sample_id_all[keep]
mutation_burden = mutation_burden_all[keep]
rm(Y_all)

## This command-line-only mode lets us test the file quickly. It is FALSE in
## RStudio and does not affect an ordinary analysis.
QUICK_TEST = "--quick" %in% commandArgs(trailingOnly = TRUE)
K_fit = if (QUICK_TEST) min(K, 8) else K
S_fit = if (QUICK_TEST) min(S, 6) else S
D_fit = if (QUICK_TEST) min(D, 3) else D
FIT_MISO_ITERS = if (QUICK_TEST) 2 else MISO_ITERS
FIT_NMF_ITERS = if (QUICK_TEST) 10 else NMF_ITERS
if (QUICK_TEST) {
  set.seed(SEED)
  chosen = sample(seq_len(nrow(Y)), min(300, nrow(Y)))
  Y = Y[chosen, , drop = FALSE]
  cancer_type = cancer_type[chosen]
  sample_id = sample_id[chosen]
  mutation_burden = mutation_burden[chosen]
}

## ========================================================================
## 3. FIT MiSo (OR LOAD THE CACHED FIT)
## ========================================================================

exclude_tag = if (length(EXCLUDED_CANCER_TYPES)) {
  gsub("[^A-Za-z0-9]+", "-", paste(EXCLUDED_CANCER_TYPES, collapse = "+"))
} else "none"
fit_tag = paste0(
  "K", K_fit, "-S", S_fit, "-D", D_fit, "-seed", SEED,
  "-iter", FIT_MISO_ITERS,
  "-init-", MISO_INITIALIZATION,
  "-n", nrow(Y), "-min", MIN_MUTATIONS, "-max", MAX_MUTATIONS,
  "-exclude-", exclude_tag, if (QUICK_TEST) "-quick" else ""
)
OUTPUT_DIR = file.path("output/pcawg-miso/explore", fit_tag)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
fit_file = file.path(OUTPUT_DIR, "miso-fit.rds")

if (file.exists(fit_file) && !REFIT) {
  message("Loading cached fit: ", fit_file)
  fit = readRDS(fit_file)
} else {
  message(
    "Fitting MiSo to ", nrow(Y), " tumours: K=", K_fit,
    ", S=", S_fit, ", D=", D_fit, " ..."
  )
  fit = miso(
    Y = Y,
    K = K_fit,
    S = S_fit,
    D = D_fit,
    init_seed = SEED,
    max_iters = FIT_MISO_ITERS,
    n_inner = if (QUICK_TEST) 1 else 2,
    mf_max_iters = if (QUICK_TEST) 4 else 20,
    mf_nmf_iters = if (QUICK_TEST) 20 else 80,
    update_prior = TRUE,
    update_F = TRUE,
    update_gamma = TRUE,
    update_slot_scale = FALSE,
    block_size = 48,
    motif_initialization = MISO_INITIALIZATION
  )
  saveRDS(fit, fit_file)
}

## Fit ordinary Poisson NMF to exactly the same matrix. This gives us the
## unconstrained part-based baseline that MiSo should improve upon or retain.
nmf_file = file.path(
  OUTPUT_DIR, paste0("nmf-fit-K", K_fit, "-iter", FIT_NMF_ITERS, ".rds")
)
if (file.exists(nmf_file) && !REFIT) {
  message("Loading cached NMF fit: ", nmf_file)
  nmf_fit = readRDS(nmf_file)
} else {
  message("Fitting vanilla Poisson NMF with K=", K_fit, " ...")
  nmf_fit = poisson_nmf_init(
    Y, K = K_fit, max_iters = FIT_NMF_ITERS, init_seed = SEED
  )
  saveRDS(nmf_fit, nmf_file)
}

## ========================================================================
## 4. CREATE SIMPLE POSTERIOR SUMMARIES
## ========================================================================

normalize_rows_simple = function(x, eps = 1e-12) {
  x = pmax(x, 0)
  x / pmax(rowSums(x), eps)
}

## Expected factor loading for every tumour (N x K). These are not constrained
## to sum to one. loading_share is only a display transformation.
loading = miso_observation_factor_scores(fit)
loading_share = normalize_rows_simple(loading)
colnames(loading) = colnames(loading_share) = paste0("F", seq_len(K_fit))
rownames(loading) = rownames(loading_share) = sample_id

## Vanilla NMF loadings for the same tumours. `nmf_loading` retains magnitude;
## `nmf_loading_share` is normalized only to make stacked bars comparable.
nmf_loading = nmf_fit$L
nmf_loading_share = normalize_rows_simple(nmf_loading)
colnames(nmf_loading) = colnames(nmf_loading_share) =
  paste0("NMF", seq_len(K_fit))
rownames(nmf_loading) = rownames(nmf_loading_share) = sample_id

## Expected factor strength in every submanifold (S x K). Row normalization
## makes the association heatmap easy to compare across motifs.
motif_factor_strength = motif_factor_scores(fit)$scores
motif_factor_share = normalize_rows_simple(motif_factor_strength)
rownames(motif_factor_strength) = rownames(motif_factor_share) =
  paste0("S", seq_len(S_fit))
colnames(motif_factor_strength) = colnames(motif_factor_share) =
  paste0("F", seq_len(K_fit))

factor_usage = colSums(loading) / sum(loading)
dominant_factor = max.col(loading_share)
dominant_submanifold = fit$z_hat
assignment_certainty = apply(fit$omega, 1, max)

## Match fitted factors to published PCAWG signatures for interpretation only.
reference = read.csv(signature_file, check.names = FALSE)
reference_context = paste0(
  substr(reference[[2]], 1, 1), "[", reference[[1]], "]",
  substr(reference[[2]], 3, 3)
)
reference_order = match(mutation_context, reference_context)
if (anyNA(reference_order)) stop("Could not align the SBS96 contexts.")
reference_F = t(as.matrix(reference[reference_order, -(1:2), drop = FALSE]))
reference_F = normalize_rows_simple(reference_F)

factor_similarity = row_cosine(fit$F, reference_F)
best_reference = max.col(factor_similarity)
matched_signature = rownames(reference_F)[best_reference]
matched_cosine = factor_similarity[cbind(seq_len(K_fit), best_reference)]
factor_label = paste0(
  "F", seq_len(K_fit), "\n", matched_signature,
  " (", sprintf("%.2f", matched_cosine), ")"
)
factor_summary = data.frame(
  factor = paste0("F", seq_len(K_fit)),
  matched_signature = matched_signature,
  cosine = matched_cosine,
  usage = factor_usage
)

nmf_factor_similarity = row_cosine(nmf_fit$F, reference_F)
nmf_best_reference = max.col(nmf_factor_similarity)
nmf_matched_signature = rownames(reference_F)[nmf_best_reference]
nmf_matched_cosine = nmf_factor_similarity[
  cbind(seq_len(K_fit), nmf_best_reference)
]
nmf_factor_label = paste0(
  "NMF", seq_len(K_fit), "\n", nmf_matched_signature,
  " (", sprintf("%.2f", nmf_matched_cosine), ")"
)
nmf_factor_usage = colSums(nmf_loading) / sum(nmf_loading)
nmf_factor_summary = data.frame(
  factor = paste0("NMF", seq_len(K_fit)),
  matched_signature = nmf_matched_signature,
  cosine = nmf_matched_cosine,
  usage = nmf_factor_usage
)

## Simple measures of "part-basedness." A non-clustering representation needs
## several factors to explain most tumours and has no single dominant loading.
active_factors_90 = function(x) {
  x = normalize_rows_simple(x)
  apply(x, 1, function(row) {
    which(cumsum(sort(row, decreasing = TRUE)) >= 0.9)[1]
  })
}
loading_entropy = function(x) {
  x = normalize_rows_simple(x)
  -rowSums(x * log(pmax(x, 1e-12))) / pmax(log(ncol(x)), 1e-12)
}
miso_active_factors_90 = active_factors_90(loading)
nmf_active_factors_90 = active_factors_90(nmf_loading)
miso_loading_entropy = loading_entropy(loading)
nmf_loading_entropy = loading_entropy(nmf_loading)

part_based_summary = data.frame(
  method = c("Vanilla Poisson NMF", "MiSo"),
  median_active_factors_90pct = c(
    median(nmf_active_factors_90), median(miso_active_factors_90)
  ),
  mean_active_factors_90pct = c(
    mean(nmf_active_factors_90), mean(miso_active_factors_90)
  ),
  median_largest_loading_share = c(
    median(apply(nmf_loading_share, 1, max)),
    median(apply(loading_share, 1, max))
  ),
  median_normalized_entropy = c(
    median(nmf_loading_entropy), median(miso_loading_entropy)
  )
)

## A compact description of every motif. `second_factor_share` is especially
## useful for seeing whether a fitted motif is genuinely multi-factor.
motif_ranked_factor = t(apply(
  motif_factor_share, 1, order, decreasing = TRUE
))
motif_effective_dimension = apply(motif_factor_share, 1, function(x) {
  which(cumsum(sort(x, decreasing = TRUE)) >= 0.9)[1]
})
dimension_diagnostics = select_motif_dimensions(fit)
dimension_index = match(
  seq_len(S_fit), dimension_diagnostics$selection$submanifold
)
motif_summary = data.frame(
  submanifold = paste0("S", seq_len(S_fit)),
  mixture_weight = fit$pi,
  assigned_tumours = tabulate(dominant_submanifold, nbins = S_fit),
  effective_factors_90pct = motif_effective_dimension,
  D_hat_ard_or_lambda =
    dimension_diagnostics$selection$D_hat_ard_or_lambda[dimension_index],
  first_factor = paste0("F", motif_ranked_factor[, 1]),
  first_signature = matched_signature[motif_ranked_factor[, 1]],
  first_factor_share = motif_factor_share[
    cbind(seq_len(S_fit), motif_ranked_factor[, 1])
  ],
  second_factor = paste0("F", motif_ranked_factor[, 2]),
  second_signature = matched_signature[motif_ranked_factor[, 2]],
  second_factor_share = motif_factor_share[
    cbind(seq_len(S_fit), motif_ranked_factor[, 2])
  ]
)

sample_summary = data.frame(
  sample_id = sample_id,
  cancer_type = cancer_type,
  mutation_burden = mutation_burden,
  submanifold = paste0("S", dominant_submanifold),
  assignment_certainty = assignment_certainty,
  dominant_factor = paste0("F", dominant_factor),
  miso_active_factors_90pct = miso_active_factors_90,
  nmf_active_factors_90pct = nmf_active_factors_90,
  miso_largest_loading_share = apply(loading_share, 1, max),
  nmf_largest_loading_share = apply(nmf_loading_share, 1, max)
)

## Cluster similar motifs, factors, and cancer types next to one another.
safe_cluster_order = function(x) {
  if (nrow(x) <= 1) return(seq_len(nrow(x)))
  stats::hclust(stats::dist(x), method = "average")$order
}
S_order = safe_cluster_order(motif_factor_share)
K_order = safe_cluster_order(t(motif_factor_share))

cancer_levels = unique(cancer_type)
cancer_mean_loading = do.call(rbind, lapply(cancer_levels, function(label) {
  colMeans(loading_share[cancer_type == label, , drop = FALSE])
}))
rownames(cancer_mean_loading) = cancer_levels
cancer_order = cancer_levels[safe_cluster_order(cancer_mean_loading)]

## Within each cancer type, put similar tumours next to one another.
sample_order = unlist(lapply(cancer_order, function(label) {
  index = which(cancer_type == label)
  local_factor = dominant_factor[index]
  local_weight = loading_share[cbind(index, local_factor)]
  index[order(dominant_submanifold[index], local_factor, -local_weight)]
}), use.names = FALSE)

## A separate within-cancer ordering reveals the structure of the NMF loading
## matrix instead of imposing MiSo's ordering on it.
nmf_dominant_factor = max.col(nmf_loading_share)
nmf_sample_order = unlist(lapply(cancer_order, function(label) {
  index = which(cancer_type == label)
  local_factor = nmf_dominant_factor[index]
  local_weight = nmf_loading_share[cbind(index, local_factor)]
  index[order(local_factor, -local_weight)]
}), use.names = FALSE)

## Cancer-type x submanifold responsibility map, used in one plot below.
cancer_motif = do.call(rbind, lapply(cancer_order, function(label) {
  colMeans(fit$omega[cancer_type == label, , drop = FALSE])
}))
rownames(cancer_motif) = cancer_order
colnames(cancer_motif) = paste0("S", seq_len(S_fit))

## PCA is exploratory only and does not change the matrix fitted by MiSo.
## Square-rooted row proportions give Hellinger geometry: this emphasizes the
## shape of each mutation spectrum instead of total mutation burden.
pca_n_components = min(PCA_COMPONENTS, ncol(Y), nrow(Y) - 1)
pca_input = sqrt(normalize_rows_simple(Y))
pca_fit = stats::prcomp(
  pca_input, center = TRUE, scale. = FALSE, rank. = pca_n_components
)
pca_scores = pca_fit$x[, seq_len(pca_n_components), drop = FALSE]
pca_variance = pca_fit$sdev^2 / sum(pca_fit$sdev^2)

## A raw-count PCA is retained only as a negative control: if its leading axis
## follows mutation burden, it should not be interpreted as biological geometry.
raw_pca_fit = stats::prcomp(Y, center = TRUE, scale. = FALSE, rank. = 2)
raw_pca_scores = raw_pca_fit$x[, 1:2, drop = FALSE]
raw_pca_variance = raw_pca_fit$sdev^2 / sum(raw_pca_fit$sdev^2)

pca_variance_summary = data.frame(
  PC = paste0("PC", seq_len(pca_n_components)),
  explained_variance = pca_variance[seq_len(pca_n_components)],
  cumulative_variance = cumsum(pca_variance)[seq_len(pca_n_components)]
)
pca_sample_summary = data.frame(
  sample_id = sample_id,
  cancer_type = cancer_type,
  mutation_burden = mutation_burden,
  submanifold = paste0("S", dominant_submanifold),
  pca_scores,
  check.names = FALSE
)

## Save small, human-readable tables. The R objects above remain available too.
utils::write.csv(factor_summary, file.path(OUTPUT_DIR, "factor-summary.csv"),
                 row.names = FALSE)
utils::write.csv(nmf_factor_summary,
                 file.path(OUTPUT_DIR, "nmf-factor-summary.csv"),
                 row.names = FALSE)
utils::write.csv(motif_summary, file.path(OUTPUT_DIR, "motif-summary.csv"),
                 row.names = FALSE)
utils::write.csv(dimension_diagnostics$slots,
                 file.path(OUTPUT_DIR, "motif-slot-summary.csv"),
                 row.names = FALSE)
utils::write.csv(sample_summary, file.path(OUTPUT_DIR, "sample-summary.csv"),
                 row.names = FALSE)
utils::write.csv(motif_factor_share,
                 file.path(OUTPUT_DIR, "motif-factor-association.csv"))
utils::write.csv(cancer_motif,
                 file.path(OUTPUT_DIR, "cancer-type-motif-association.csv"))
utils::write.csv(part_based_summary,
                 file.path(OUTPUT_DIR, "part-basedness-summary.csv"),
                 row.names = FALSE)
utils::write.csv(pca_variance_summary,
                 file.path(OUTPUT_DIR, "pca-explained-variance.csv"),
                 row.names = FALSE)
utils::write.csv(pca_sample_summary,
                 file.path(OUTPUT_DIR, "pca-sample-scores.csv"),
                 row.names = FALSE)

## ========================================================================
## 5. PLOTTING FUNCTIONS
## ========================================================================

factor_colors = grDevices::hcl.colors(K_fit, "Dynamic")
names(factor_colors) = paste0("F", seq_len(K_fit))
nmf_factor_colors = grDevices::hcl.colors(K_fit, "Dynamic")
names(nmf_factor_colors) = paste0("NMF", seq_len(K_fit))

## Save each plot as a PNG. In RStudio, draw it in the Plot pane as well.
save_and_show = function(filename, width, height, plot_function) {
  grDevices::png(
    file.path(OUTPUT_DIR, filename), width = width, height = height,
    res = 180, bg = "white"
  )
  plot_function()
  grDevices::dev.off()
  if (interactive()) plot_function()
}

## Plot 1: the requested S x K association map.
plot_motif_factor_map = function() {
  graphics::par(mar = c(10, 8, 4, 2))
  shown = motif_factor_share[S_order, K_order, drop = FALSE]
  graphics::image(
    seq_len(K_fit), seq_len(S_fit), t(shown),
    col = rev(grDevices::hcl.colors(100, "YlOrRd")),
    axes = FALSE, xlab = "Shared factor", ylab = "",
    main = "MiSo S x K motif-factor association (row-normalized)",
    useRaster = TRUE
  )
  graphics::axis(1, at = seq_len(K_fit), labels = factor_label[K_order],
                 las = 2, cex.axis = if (K_fit > 24) 0.42 else 0.55)
  graphics::axis(
    2, at = seq_len(S_fit),
    labels = paste0("S", S_order, "  pi=", sprintf("%.2f", fit$pi[S_order])),
    las = 1, cex.axis = 0.65
  )
  graphics::mtext("Darker cells indicate a larger within-motif factor share.",
                  side = 3, line = 0.5, cex = 0.72)
  graphics::box()
}

## Plot 2: a conventional loading barplot, grouped by the given cancer labels.
## Rows are normalized only for this plot; `loading` above remains unnormalized.
plot_loading_barplot = function() {
  graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(5, 0.9))
  graphics::par(mar = c(12, 4, 4, 1))
  midpoint = graphics::barplot(
    t(loading_share[sample_order, , drop = FALSE]),
    space = 0, border = NA, col = factor_colors, axes = FALSE,
    names.arg = rep("", length(sample_order)),
    xlab = "", ylab = "Loading share (display only)",
    main = "MiSo loadings arranged by cancer type"
  )
  graphics::axis(2, at = c(0, 0.5, 1), las = 1)
  ordered_label = factor(cancer_type[sample_order], levels = cancer_order)
  label_midpoint = tapply(midpoint, ordered_label, mean)
  graphics::axis(1, at = label_midpoint, labels = names(label_midpoint),
                 las = 2, cex.axis = 0.55)
  group_size = as.integer(table(ordered_label))
  group_end = cumsum(group_size)
  if (length(group_end) > 1) {
    boundary = (midpoint[group_end[-length(group_end)]] +
      midpoint[group_end[-length(group_end)] + 1]) / 2
    graphics::abline(v = boundary, col = "white", lwd = 1)
  }
  graphics::box()

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center", legend = names(factor_colors), fill = factor_colors,
    ncol = min(9, K_fit), bty = "n", cex = 0.62,
    title = "Shared factors"
  )
  graphics::layout(1)
}

## Plot 2b: ordinary NMF loadings, grouped by the same cancer labels. Comparing
## this with Plot 2 shows directly whether MiSo preserved part-based structure.
plot_nmf_loading_barplot = function() {
  graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(5, 0.9))
  graphics::par(mar = c(12, 4, 4, 1))
  midpoint = graphics::barplot(
    t(nmf_loading_share[nmf_sample_order, , drop = FALSE]),
    space = 0, border = NA, col = nmf_factor_colors, axes = FALSE,
    names.arg = rep("", length(nmf_sample_order)),
    xlab = "", ylab = "Loading share (display only)",
    main = "Vanilla Poisson NMF loadings arranged by cancer type"
  )
  graphics::axis(2, at = c(0, 0.5, 1), las = 1)
  ordered_label = factor(cancer_type[nmf_sample_order], levels = cancer_order)
  label_midpoint = tapply(midpoint, ordered_label, mean)
  graphics::axis(1, at = label_midpoint, labels = names(label_midpoint),
                 las = 2, cex.axis = 0.55)
  group_size = as.integer(table(ordered_label))
  group_end = cumsum(group_size)
  if (length(group_end) > 1) {
    boundary = (midpoint[group_end[-length(group_end)]] +
      midpoint[group_end[-length(group_end)] + 1]) / 2
    graphics::abline(v = boundary, col = "white", lwd = 1)
  }
  graphics::box()

  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center", legend = nmf_factor_label, fill = nmf_factor_colors,
    ncol = min(9, K_fit), bty = "n", cex = 0.54,
    title = "NMF factor -> closest published signature (cosine)"
  )
  graphics::layout(1)
}

## Plot 3: which cancer labels are associated with which submanifolds?
plot_cancer_motif_map = function() {
  graphics::par(mar = c(8, 16, 4, 2))
  shown = cancer_motif[, S_order, drop = FALSE]
  graphics::image(
    seq_len(S_fit), seq_len(nrow(shown)), t(shown),
    col = rev(grDevices::hcl.colors(100, "Blues 3")),
    axes = FALSE, xlab = "Submanifold", ylab = "",
    main = "Mean posterior submanifold responsibility by cancer type",
    useRaster = TRUE
  )
  graphics::axis(1, at = seq_len(S_fit), labels = paste0("S", S_order),
                 las = 2, cex.axis = 0.65)
  graphics::axis(2, at = seq_len(nrow(shown)), labels = rownames(shown),
                 las = 1, cex.axis = 0.58)
  graphics::box()
}

## Plot 4: magnitude and mixture diagnostics hidden by a normalized barplot.
plot_scale_and_assignment = function() {
  graphics::layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1.5, 1, 1))
  graphics::par(mar = c(5, 10, 4, 1))
  graphics::boxplot(
    log10(mutation_burden + 1) ~ factor(cancer_type, levels = cancer_order),
    horizontal = TRUE, las = 1, outline = FALSE, cex.axis = 0.55,
    xlab = "log10 observed SBS count", ylab = "",
    main = "Mutation burden"
  )
  graphics::par(mar = c(8, 4, 4, 1))
  graphics::barplot(
    fit$pi, names.arg = paste0("S", seq_len(S_fit)), las = 2,
    col = "#67A9CF", border = NA, ylab = "Mixture weight",
    main = "Submanifold prevalence", cex.names = 0.65
  )
  graphics::par(mar = c(5, 4, 4, 1))
  graphics::hist(
    assignment_certainty, breaks = 20, col = "#F4A582", border = "white",
    xlab = "Maximum posterior responsibility", main = "Assignment certainty"
  )
  graphics::layout(1)
}

## Plot 4b: quantitative check of the clustering-versus-parts concern.
plot_part_basedness = function() {
  method_colors = c("#67A9CF", "#F4A582")
  graphics::layout(matrix(seq_len(3), nrow = 1))
  graphics::par(mar = c(7, 4, 4, 1))
  graphics::boxplot(
    list(NMF = nmf_active_factors_90, MiSo = miso_active_factors_90),
    col = method_colors, border = "#555555", outline = FALSE,
    ylab = "Factors needed for 90% of loading", main = "Active parts"
  )
  graphics::boxplot(
    list(
      NMF = apply(nmf_loading_share, 1, max),
      MiSo = apply(loading_share, 1, max)
    ),
    col = method_colors, border = "#555555", outline = FALSE,
    ylim = c(0, 1), ylab = "Largest loading share",
    main = "Single-factor dominance"
  )
  graphics::boxplot(
    list(NMF = nmf_loading_entropy, MiSo = miso_loading_entropy),
    col = method_colors, border = "#555555", outline = FALSE,
    ylim = c(0, 1), ylab = "Normalized loading entropy",
    main = "Loading diversity"
  )
  graphics::layout(1)
}

## Plot 4c: PCs 1--10 of the Hellinger-transformed mutation spectra. Smooth
## arms or bridges suggest continuous geometry; isolated clouds suggest blobs.
plot_pca_geometry = function() {
  n_pairs = floor(pca_n_components / 2)
  pair_index = head(seq_len(n_pairs), 5)
  graphics::par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3.2, 1))
  point_color = grDevices::adjustcolor("#2B8CBE", alpha.f = 0.28)
  for (j in pair_index) {
    x_pc = 2 * j - 1
    y_pc = 2 * j
    graphics::plot(
      pca_scores[, x_pc], pca_scores[, y_pc],
      pch = 16, cex = 0.42, col = point_color,
      xlab = sprintf("PC%d (%.1f%%)", x_pc, 100 * pca_variance[x_pc]),
      ylab = sprintf("PC%d (%.1f%%)", y_pc, 100 * pca_variance[y_pc]),
      main = paste0("PC", x_pc, " versus PC", y_pc)
    )
  }
  graphics::barplot(
    100 * pca_variance[seq_len(pca_n_components)],
    names.arg = paste0("PC", seq_len(pca_n_components)), las = 2,
    col = "#74A9CF", border = NA,
    ylab = "Explained variance (%)", main = "Hellinger-PCA scree"
  )
}

## Plot 4d: the same PC pairs colored by the given cancer labels. Labels are
## not used to calculate the PCs or to fit MiSo.
plot_pca_by_cancer_type = function() {
  n_pairs = floor(pca_n_components / 2)
  pair_index = head(seq_len(n_pairs), 5)
  cancer_palette = setNames(
    grDevices::hcl.colors(length(cancer_order), "Dynamic"), cancer_order
  )
  point_color = grDevices::adjustcolor(
    cancer_palette[cancer_type], alpha.f = 0.55
  )
  graphics::par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3.2, 1))
  for (j in pair_index) {
    x_pc = 2 * j - 1
    y_pc = 2 * j
    graphics::plot(
      pca_scores[, x_pc], pca_scores[, y_pc],
      pch = 16, cex = 0.40, col = point_color,
      xlab = sprintf("PC%d (%.1f%%)", x_pc, 100 * pca_variance[x_pc]),
      ylab = sprintf("PC%d (%.1f%%)", y_pc, 100 * pca_variance[y_pc]),
      main = paste0("PC", x_pc, " versus PC", y_pc)
    )
  }
  graphics::plot.new()
  graphics::legend(
    "center", legend = cancer_order, col = cancer_palette[cancer_order],
    pch = 16, ncol = 3, bty = "n", cex = 0.58,
    title = "Cancer type (display only)"
  )
}

## Plot 4e: why we do not use PCA of raw counts for the geometry diagnostic.
plot_pca_burden_control = function() {
  log_burden = log10(mutation_burden)
  burden_palette = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE)
  color_index = cut(log_burden, breaks = 100, labels = FALSE,
                    include.lowest = TRUE)
  point_color = grDevices::adjustcolor(
    burden_palette[color_index], alpha.f = 0.55
  )
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 4, 1))
  raw_cor = abs(stats::cor(raw_pca_scores[, 1], log_burden))
  graphics::plot(
    raw_pca_scores[, 1], raw_pca_scores[, 2],
    pch = 16, cex = 0.42, col = point_color,
    xlab = sprintf("Raw PC1 (%.1f%%)", 100 * raw_pca_variance[1]),
    ylab = sprintf("Raw PC2 (%.1f%%)", 100 * raw_pca_variance[2]),
    main = sprintf("Raw counts: |cor(PC1, log burden)| = %.2f", raw_cor)
  )
  hellinger_cor = abs(stats::cor(pca_scores[, 1], log_burden))
  graphics::plot(
    pca_scores[, 1], pca_scores[, 2],
    pch = 16, cex = 0.42, col = point_color,
    xlab = sprintf("Hellinger PC1 (%.1f%%)", 100 * pca_variance[1]),
    ylab = sprintf("Hellinger PC2 (%.1f%%)", 100 * pca_variance[2]),
    main = sprintf(
      "Mutation proportions: |cor(PC1, log burden)| = %.2f",
      hellinger_cor
    )
  )
}

## Plot 5: spectra for the most-used learned factors.
plot_factor_spectra = function() {
  shown_factor = head(order(factor_usage, decreasing = TRUE),
                      min(N_FACTORS_TO_SHOW, K_fit))
  ncol_plot = 3
  nrow_plot = ceiling(length(shown_factor) / ncol_plot)
  graphics::par(mfrow = c(nrow_plot, ncol_plot), mar = c(2, 3, 3, 1))
  class_palette = setNames(
    c("#03BCEE", "#010101", "#E32926", "#CAC9C9", "#A1CE63", "#EBC6C4"),
    c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  )
  for (k in shown_factor) {
    spectrum_midpoint = graphics::barplot(
      fit$F[k, ], space = 0, border = NA,
      col = class_palette[mutation_class], axes = FALSE,
      names.arg = rep("", ncol(fit$F)),
      main = paste0(
        "F", k, " -> ", matched_signature[k],
        "  cos=", sprintf("%.2f", matched_cosine[k]),
        "  use=", sprintf("%.0f%%", 100 * factor_usage[k])
      ), cex.main = 0.75
    )
    graphics::axis(2, las = 1, cex.axis = 0.55)
    class_midpoint = tapply(
      spectrum_midpoint,
      factor(mutation_class, levels = names(class_palette)), mean
    )
    graphics::axis(1, at = class_midpoint, labels = names(class_midpoint),
                   tick = FALSE, cex.axis = 0.58)
    boundary = (spectrum_midpoint[c(16, 32, 48, 64, 80)] +
      spectrum_midpoint[c(17, 33, 49, 65, 81)]) / 2
    graphics::abline(v = boundary, col = "white", lwd = 0.8)
  }
  unused = nrow_plot * ncol_plot - length(shown_factor)
  if (unused > 0) for (j in seq_len(unused)) graphics::plot.new()
}

## Plot 6: continuous geometry in large motifs with two distinct factors.
## x retains loading magnitude; y is an angular display coordinate. The model
## itself does not constrain a tumour's loading vector to sum to one.
plot_continuous_motifs = function() {
  learned_cosine = row_cosine(fit$F, fit$F)
  motif_pair = lapply(seq_len(S_fit), function(s) {
    ordered = order(motif_factor_share[s, ], decreasing = TRUE)
    first = ordered[1]
    alternatives = ordered[
      ordered != first &
        motif_factor_share[s, ordered] >= SECOND_FACTOR_SHARE &
        learned_cosine[first, ordered] < 0.95
    ]
    if (!length(alternatives)) return(c(NA_integer_, NA_integer_))
    different_process = alternatives[
      matched_signature[alternatives] != matched_signature[first]
    ]
    second = if (length(different_process)) different_process[1]
    else alternatives[1]
    c(first, second)
  })
  valid = vapply(motif_pair, function(x) !anyNA(x), logical(1))
  assigned = tabulate(dominant_submanifold, nbins = S_fit)
  candidate = order(ifelse(valid & assigned >= 20, fit$pi, -Inf),
                    decreasing = TRUE)
  candidate = head(candidate[is.finite(fit$pi[candidate]) & valid[candidate] &
                               assigned[candidate] >= 20],
                   N_MOTIFS_TO_SHOW)

  if (!length(candidate)) {
    graphics::par(mar = c(2, 2, 2, 2))
    graphics::plot.new()
    graphics::text(
      0.5, 0.58,
      paste0(
        "No fitted motif has two distinct factors with at least ",
        sprintf("%.0f%%", 100 * SECOND_FACTOR_SHARE), " association each."
      ),
      cex = 1.25, font = 2
    )
    graphics::text(
      0.5, 0.44,
      paste0(
        "Inspect View(motif_summary), lower SECOND_FACTOR_SHARE, or revisit ",
        "factor/motif postprocessing."
      ),
      cex = 0.95
    )
    return(invisible(NULL))
  }

  ncol_plot = 2
  nrow_plot = ceiling(max(1, length(candidate)) / ncol_plot)
  graphics::par(mfrow = c(nrow_plot, ncol_plot), mar = c(5, 5, 4, 1))
  label_palette = setNames(
    grDevices::hcl.colors(length(unique(cancer_type)), "Dynamic"),
    unique(cancer_type)
  )
  for (s in candidate) {
    factors = motif_pair[[s]]
    index = dominant_submanifold == s
    common_label = names(head(sort(table(cancer_type[index]), decreasing = TRUE), 4))
    color_group = ifelse(cancer_type[index] %in% common_label,
                         cancer_type[index], "Other")
    local_palette = c(label_palette[common_label], Other = "#BDBDBD")

    first_loading = loading[index, factors[1]]
    second_loading = loading[index, factors[2]]
    amplitude = first_loading + second_loading
    angle = first_loading / pmax(amplitude, 1e-12)

    graphics::plot(
      log10(amplitude + 1), angle, pch = 19,
      col = grDevices::adjustcolor(local_palette[color_group], alpha.f = 0.65),
      xlab = paste0("log10 amplitude: F", factors[1], " + F", factors[2]),
      ylab = paste0("F", factors[1], "/(F", factors[1], "+F", factors[2], ")"),
      main = paste0(
        "S", s, ": n=", sum(index), ", pi=", sprintf("%.2f", fit$pi[s]),
        "\n", matched_signature[factors[1]], " <-> ",
        matched_signature[factors[2]]
      )
    )
    graphics::legend("topright", legend = names(local_palette),
                     col = local_palette, pch = 19, bty = "n", cex = 0.58)
  }
  unused = nrow_plot * ncol_plot - length(candidate)
  if (unused > 0) for (j in seq_len(unused)) graphics::plot.new()
}

## ========================================================================
## 6. MAKE THE PLOTS
## ========================================================================

save_and_show("01-motif-factor-map.png", 2200, 1700,
              plot_motif_factor_map)
save_and_show("02-loading-barplot-by-cancer-type.png", 3000, 1500,
              plot_loading_barplot)
save_and_show("03-cancer-type-motif-map.png", 1800, 2000,
              plot_cancer_motif_map)
save_and_show("04-scale-and-assignment.png", 2600, 1500,
              plot_scale_and_assignment)
save_and_show("05-factor-spectra.png", 2200, 1900,
              plot_factor_spectra)
save_and_show("06-continuous-motifs.png", 2200, 1900,
              plot_continuous_motifs)
save_and_show("07-nmf-loading-barplot-by-cancer-type.png", 3000, 1500,
              plot_nmf_loading_barplot)
save_and_show("08-nmf-vs-miso-part-basedness.png", 2400, 1200,
              plot_part_basedness)
save_and_show("09-pca-hellinger-pc1-pc10.png", 2400, 1700,
              plot_pca_geometry)
save_and_show("10-pca-hellinger-by-cancer-type.png", 2400, 1700,
              plot_pca_by_cancer_type)
save_and_show("11-pca-burden-control.png", 2200, 1100,
              plot_pca_burden_control)

cat(
  "\nFinished. Output:", OUTPUT_DIR,
  "\nUseful RStudio objects:",
  "\n  View(factor_summary)",
  "\n  View(nmf_factor_summary)",
  "\n  View(motif_summary)",
  "\n  View(sample_summary)",
  "\n  View(motif_factor_share)",
  "\n  View(cancer_motif)",
  "\n  View(part_based_summary)",
  "\n  View(pca_variance_summary)",
  "\n  View(pca_sample_summary)",
  "\nRe-run any plot, for example: plot_motif_factor_map()\n"
)
print(part_based_summary, row.names = FALSE)

if (!isTRUE(fit$converged)) {
  message(
    "Note: this exploratory fit did not meet the formal stopping rule after ",
    fit$n_iter, " iterations. Increase MISO_ITERS to refit longer."
  )
}
