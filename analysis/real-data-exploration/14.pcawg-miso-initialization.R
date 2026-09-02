#!/usr/bin/env Rscript

## Controlled PCAWG MiSo initialization experiment.
##
## This file keeps Y, K, S, D, and the MiSo updates fixed. It varies only how
## the S x D factor-selection distributions are initialized, except for the
## final variant, which deliberately bypasses Poisson-SuSiE-NMF and starts from
## the vanilla Poisson NMF factors and loadings.

## -------------------------------------------------------------------------
## 1. Settings
## -------------------------------------------------------------------------

K = 36
S = 30
D = 5
SEED = 1
FIT_ITERS = 15
REFIT = FALSE

MIN_MUTATIONS = 500
MAX_MUTATIONS = 50000
EXCLUDED_CANCER_TYPES = "Skin-Melanoma"

MAIN_FIT_DIR = paste0(
  "output/pcawg-miso/explore/",
  "K36-S30-D5-seed1-iter10-n2438-min500-max50000-",
  "exclude-Skin-Melanoma"
)
OUTPUT_DIR = file.path(
  "output/pcawg-miso/initialization-study",
  paste0("K", K, "-S", S, "-D", D, "-iter", FIT_ITERS)
)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists("code/miso-benchmark-utils.R")) {
  stop("Open the project in RStudio and run this file from the project root.")
}
source("code/miso-benchmark-utils.R")

## -------------------------------------------------------------------------
## 2. Load exactly the same raw count matrix used in the exploration file
## -------------------------------------------------------------------------

spectra = read.csv(
  "data/pcawg/WGS_PCAWG.96.csv", check.names = FALSE
)
sample_key = names(spectra)[-(1:2)]
cancer_type_all = sub("::.*$", "", sample_key)
sample_id_all = sub("^.*::", "", sample_key)

Y_all = t(as.matrix(spectra[, -(1:2), drop = FALSE]))
storage.mode(Y_all) = "double"
mutation_burden_all = rowSums(Y_all)
keep = mutation_burden_all >= MIN_MUTATIONS &
  mutation_burden_all <= MAX_MUTATIONS &
  !(cancer_type_all %in% EXCLUDED_CANCER_TYPES)

Y = Y_all[keep, , drop = FALSE]
cancer_type = cancer_type_all[keep]
sample_id = sample_id_all[keep]
mutation_burden = mutation_burden_all[keep]
rm(Y_all)

## Reuse the already fitted Poisson-SuSiE-NMF and vanilla NMF starts. This
## makes the experiment fast and isolates the MiSo initialization step.
main_fit = readRDS(file.path(MAIN_FIT_DIR, "miso-fit.rds"))
nmf_fit = readRDS(file.path(MAIN_FIT_DIR, "nmf-fit-K36-iter200.rds"))
if (!all(dim(main_fit$F) == c(K, ncol(Y)))) stop("MiSo dimensions differ.")
if (!all(dim(nmf_fit$F) == c(K, ncol(Y)))) stop("NMF dimensions differ.")

susie_scores = poisson_susie_nmf_loading_scores(main_fit$mf_fit)
nmf_scores = nmf_fit$L

## -------------------------------------------------------------------------
## 3. Define four initialization strategies
## -------------------------------------------------------------------------

## Legacy initializer: retain cluster-center factors above 10%, then repeat.
repeat_init = init_soft_gamma_from_loading_scores(
  susie_scores, S = S, D = D, init_seed = SEED,
  min_share = 0.10, gamma_floor = 0.05, surplus_slots = "repeat",
  motif_initialization = "threshold"
)$gamma_bar

## Top-D distinct factors from the same Poisson-SuSiE-NMF loading centers.
make_distinct_init = function(scores, gamma_floor = 0.05) {
  hard = init_motifs_from_loading_scores(
    scores, S = S, D = D, init_seed = SEED,
    min_share = 0, allow_repeats = FALSE
  )
  gamma_bar_from_motifs(hard$motifs, K = K, gamma_floor = gamma_floor)
}
susie_distinct_init = make_distinct_init(susie_scores)

## Independent soft random distributions break slot symmetry without using
## loading information. This is a sensitivity check, not a preferred default.
set.seed(SEED)
random_soft_init = array(
  rgamma(S * D * K, shape = 1, rate = 1), dim = c(S, D, K)
)
random_soft_init = normalize_gamma_bar(random_soft_init)

## Bypass Poisson-SuSiE-NMF for initialization: use vanilla NMF loadings to
## choose top-D distinct factors and use vanilla NMF factors as the starting F.
nmf_distinct_init = make_distinct_init(nmf_scores)

initializations = list(
  repeat_susie = list(
    label = "Repeated SuSiE factors (current)",
    F = main_fit$mf_fit$F, gamma = repeat_init
  ),
  distinct_susie = list(
    label = "Top-D distinct SuSiE factors",
    F = main_fit$mf_fit$F, gamma = susie_distinct_init
  ),
  random_soft = list(
    label = "Random soft factors",
    F = main_fit$mf_fit$F, gamma = random_soft_init
  ),
  distinct_nmf = list(
    label = "Top-D distinct vanilla-NMF factors",
    F = nmf_fit$F, gamma = nmf_distinct_init
  )
)

## -------------------------------------------------------------------------
## 4. Fit each variant, caching the results
## -------------------------------------------------------------------------

fits = vector("list", length(initializations))
names(fits) = names(initializations)
for (name in names(initializations)) {
  fit_file = file.path(OUTPUT_DIR, paste0("fit-", name, ".rds"))
  if (file.exists(fit_file) && !REFIT) {
    message("Loading ", name, " ...")
    fits[[name]] = readRDS(fit_file)
  } else {
    message("Fitting ", name, " ...")
    init = initializations[[name]]
    fits[[name]] = miso_fixed_gamma(
      Y = Y,
      F = init$F,
      gamma_bar = init$gamma,
      max_iters = FIT_ITERS,
      n_inner = 2,
      update_prior = TRUE,
      update_F = TRUE,
      update_gamma = TRUE,
      update_slot_scale = FALSE,
      gamma_step_init = 0.5,
      gamma_step_ramp = 10,
      F_step_init = 0.2,
      F_step_ramp = 20,
      block_size = 48
    )
    saveRDS(fits[[name]], fit_file)
  }
}

## -------------------------------------------------------------------------
## 5. Compare the resulting loading geometry and fit
## -------------------------------------------------------------------------

normalize_rows_simple = function(x, eps = 1e-12) {
  x / pmax(rowSums(x), eps)
}
active_factors_90 = function(x) {
  x = normalize_rows_simple(x)
  apply(x, 1, function(row) {
    which(cumsum(sort(row, decreasing = TRUE)) >= 0.9)[1]
  })
}
loading_entropy = function(x) {
  x = normalize_rows_simple(x)
  -rowSums(x * log(pmax(x, 1e-12))) / log(ncol(x))
}
unique_slot_winners = function(gamma_bar) {
  winner = apply(gamma_bar, c(1, 2), which.max)
  apply(winner, 1, function(row) length(unique(row)))
}

fit_summaries = lapply(names(fits), function(name) {
  fit = fits[[name]]
  loading = miso_observation_factor_scores(fit)
  loading_share = normalize_rows_simple(loading)
  fitted_mean = loading %*% fit$F
  final_elbo = tail(fit$elbo[is.finite(fit$elbo)], 1)
  init_unique = unique_slot_winners(initializations[[name]]$gamma)
  final_unique = unique_slot_winners(fit$gamma_bar)

  data.frame(
    initialization = name,
    label = initializations[[name]]$label,
    converged = fit$converged,
    iterations = fit$n_iter,
    final_elbo = final_elbo,
    poisson_deviance_per_entry =
      poisson_deviance(Y, fitted_mean) / length(Y),
    median_active_factors_90pct = median(active_factors_90(loading)),
    mean_active_factors_90pct = mean(active_factors_90(loading)),
    median_largest_loading_share = median(apply(loading_share, 1, max)),
    median_loading_entropy = median(loading_entropy(loading)),
    median_initial_distinct_slots = median(init_unique),
    median_final_distinct_slots = median(final_unique),
    motifs_with_two_or_more_final_factors = sum(final_unique >= 2),
    median_max_gamma = median(apply(fit$gamma_bar, c(1, 2), max)),
    median_max_responsibility = median(apply(fit$omega, 1, max))
  )
})
comparison = do.call(rbind, fit_summaries)
rownames(comparison) = NULL
utils::write.csv(
  comparison, file.path(OUTPUT_DIR, "initialization-comparison.csv"),
  row.names = FALSE
)

## Keep the full posterior loadings available in RStudio.
loadings = lapply(fits, miso_observation_factor_scores)
loading_shares = lapply(loadings, normalize_rows_simple)

## Use one fixed tumor ordering for every panel so differences are not created
## by rearranging samples. The order comes from vanilla NMF within cancer type.
nmf_share = normalize_rows_simple(nmf_scores)
nmf_dominant = max.col(nmf_share)
cancer_order = names(sort(table(cancer_type), decreasing = TRUE))
sample_order = unlist(lapply(cancer_order, function(label) {
  index = which(cancer_type == label)
  weight = nmf_share[cbind(index, nmf_dominant[index])]
  index[order(nmf_dominant[index], -weight)]
}), use.names = FALSE)

factor_colors = grDevices::hcl.colors(K, "Dynamic")

open_png = function(filename, width = 2600, height = 1900) {
  grDevices::png(
    file.path(OUTPUT_DIR, filename), width = width, height = height,
    res = 180, bg = "white"
  )
}

## Plot 1: posterior loading bars for all initializations on the same ordering.
open_png("01-loading-comparison.png", 2800, 2200)
graphics::par(mfrow = c(2, 2), mar = c(2.2, 4, 4, 1))
for (name in names(fits)) {
  row = comparison[comparison$initialization == name, ]
  graphics::barplot(
    t(loading_shares[[name]][sample_order, , drop = FALSE]),
    space = 0, border = NA, col = factor_colors, axes = FALSE,
    names.arg = rep("", length(sample_order)),
    ylab = "Loading share", xlab = "Tumors grouped by cancer type",
    main = paste0(
      initializations[[name]]$label,
      "\nmedian active factors = ", row$median_active_factors_90pct,
      "; multi-factor motifs = ",
      row$motifs_with_two_or_more_final_factors, "/", S
    )
  )
  graphics::axis(2, at = c(0, 0.5, 1), las = 1)
  graphics::box()
}
grDevices::dev.off()

## Plot 2: direct numerical comparison of part-basedness and fit.
open_png("02-initialization-diagnostics.png", 2500, 1500)
graphics::par(mfrow = c(2, 2), mar = c(9, 4.5, 4, 1))
short_label = c("repeat", "distinct\nSuSiE", "random\nsoft", "distinct\nNMF")
bar_color = c("#D73027", "#4575B4", "#74ADD1", "#66BD63")
graphics::barplot(
  comparison$median_active_factors_90pct, names.arg = short_label,
  col = bar_color, border = NA, ylab = "Median factors for 90% loading",
  main = "Part-based loading dimension"
)
graphics::barplot(
  comparison$motifs_with_two_or_more_final_factors, names.arg = short_label,
  col = bar_color, border = NA, ylim = c(0, S),
  ylab = paste0("Motifs with >=2 factor directions (of ", S, ")"),
  main = "Submanifold dimensionality"
)
graphics::barplot(
  comparison$median_largest_loading_share, names.arg = short_label,
  col = bar_color, border = NA, ylim = c(0, 1),
  ylab = "Median largest loading share",
  main = "Single-factor dominance"
)
graphics::barplot(
  comparison$poisson_deviance_per_entry, names.arg = short_label,
  col = bar_color, border = NA,
  ylab = "In-sample Poisson deviance per entry",
  main = "Reconstruction fit (lower is better)"
)
grDevices::dev.off()

## Plot 3: S x K motif-factor associations after each fit.
open_png("03-motif-factor-comparison.png", 2500, 1900)
graphics::par(mfrow = c(2, 2), mar = c(6, 5, 4, 1))
for (name in names(fits)) {
  motif_scores = motif_factor_scores(fits[[name]])$scores
  motif_share = normalize_rows_simple(motif_scores)
  graphics::image(
    seq_len(K), seq_len(S), t(motif_share),
    col = rev(grDevices::hcl.colors(100, "YlOrRd")), axes = FALSE,
    xlab = "Factor", ylab = "Submanifold",
    main = initializations[[name]]$label, useRaster = TRUE
  )
  graphics::axis(1, at = seq_len(K), labels = seq_len(K),
                 las = 2, cex.axis = 0.48)
  graphics::axis(2, at = seq_len(S), labels = seq_len(S),
                 las = 1, cex.axis = 0.48)
  graphics::box()
}
grDevices::dev.off()

print(comparison, row.names = FALSE)
cat(
  "\nFinished. Output:", OUTPUT_DIR,
  "\nUseful objects: View(comparison), loadings, loading_shares, fits\n"
)
