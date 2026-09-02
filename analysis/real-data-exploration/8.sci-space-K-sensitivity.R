## Quick sensitivity check for the global NMF rank K.
##
## Run from the project root:
##   Rscript --vanilla analysis/real-data-exploration/8.sci-space-K-sensitivity.R
##
## This reuses the exact cells, genes, and Poisson-thinning split from the
## larger sci-Space experiment. It fits only Poisson NMF, so it is much faster
## than repeating the complete MiSo experiment for every candidate K.

## -------------------------------------------------------------------------
## 1. Settings
## -------------------------------------------------------------------------

K_VALUES = c(14, 18, 22)
FIT_SEEDS = c(1, 2, 3)
S_FOR_KMEANS = 10
NMF_ITERS = 35
ACTIVE_FACTOR_THRESHOLD = 0.05
USED_FACTOR_THRESHOLD = 0.01
RELATIVE_DEVIANCE_TOLERANCE = 0.02
TRAIN_FRACTION = 0.80
DATA_SEED = 20260830

INPUT_FILE = paste0(
  "output/sci-space-miso-large/",
  "preprocessed-n2000-m800-seed20260830.rds"
)
OUTPUT_DIR = "output/sci-space-miso-large/K-sensitivity"

## -------------------------------------------------------------------------
## 2. Data and helpers
## -------------------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop("Run analysis/real-data-exploration/7.sci-space-miso-large.R first.")
}
source("code/miso-benchmark-utils.R")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

dat = readRDS(INPUT_FILE)
set.seed(DATA_SEED + 1)
Y_train = matrix(
  rbinom(length(dat$Y), size = as.vector(dat$Y), prob = TRAIN_FRACTION),
  nrow = nrow(dat$Y), ncol = ncol(dat$Y), dimnames = dimnames(dat$Y)
)
Y_test = dat$Y - Y_train
lineage = dat$meta$final_cluster_label

normalized_mutual_information = function(x, y, eps = 1e-12) {
  probability = table(x, y) / length(x)
  px = rowSums(probability)
  py = colSums(probability)
  expected = px %o% py
  positive = probability > 0
  mutual_information = sum(
    probability[positive] *
      log(probability[positive] / pmax(expected[positive], eps))
  )
  entropy_x = -sum(px[px > 0] * log(px[px > 0]))
  entropy_y = -sum(py[py > 0] * log(py[py > 0]))
  mutual_information / sqrt(pmax(entropy_x * entropy_y, eps))
}

evaluate_nmf = function(K, seed) {
  cache_file = file.path(
    OUTPUT_DIR, sprintf("nmf-K%02d-seed-%03d.rds", K, seed)
  )
  if (file.exists(cache_file)) {
    message("Loading K=", K, ", seed=", seed)
    fit = readRDS(cache_file)
  } else {
    message("Fitting K=", K, ", seed=", seed)
    fit = poisson_nmf_init(
      Y_train, K = K, max_iters = NMF_ITERS, init_seed = seed
    )
    saveRDS(fit, cache_file)
  }

  loading_proportion = normalize_scores(fit$L)
  set.seed(seed)
  cluster = kmeans(
    normalize_scores(fit$L),
    centers = S_FOR_KMEANS,
    nstart = 20,
    iter.max = 100
  )$cluster
  train_mean = fit$L %*% fit$F
  test_mean = pmax(
    (1 - TRAIN_FRACTION) / TRAIN_FRACTION * train_mean,
    1e-12
  )

  factor_similarity = row_cosine(fit$F, fit$F)
  diag(factor_similarity) = NA_real_
  upper_similarity = factor_similarity[upper.tri(factor_similarity)]
  factor_usage = colMeans(loading_proportion)

  data.frame(
    K = K,
    seed = seed,
    heldout_deviance_per_entry =
      poisson_deviance(Y_test, test_mean) / length(Y_test),
    lineage_nmi = normalized_mutual_information(cluster, lineage),
    mean_active_factors = mean(
      rowSums(loading_proportion >= ACTIVE_FACTOR_THRESHOLD)
    ),
    used_factors = sum(factor_usage >= USED_FACTOR_THRESHOLD),
    max_factor_cosine = max(upper_similarity, na.rm = TRUE),
    redundant_factor_pairs = sum(upper_similarity >= 0.90, na.rm = TRUE)
  )
}

## -------------------------------------------------------------------------
## 3. Fit and summarize
## -------------------------------------------------------------------------

jobs = expand.grid(K = K_VALUES, seed = FIT_SEEDS)
results = do.call(rbind, lapply(seq_len(nrow(jobs)), function(j) {
  evaluate_nmf(jobs$K[j], jobs$seed[j])
}))

summary = do.call(rbind, lapply(split(results, results$K), function(tab) {
  data.frame(
    K = unique(tab$K),
    heldout_deviance_mean = mean(tab$heldout_deviance_per_entry),
    heldout_deviance_sd = sd(tab$heldout_deviance_per_entry),
    lineage_nmi_mean = mean(tab$lineage_nmi),
    mean_active_factors = mean(tab$mean_active_factors),
    mean_used_factors = mean(tab$used_factors),
    mean_max_factor_cosine = mean(tab$max_factor_cosine),
    mean_redundant_factor_pairs = mean(tab$redundant_factor_pairs)
  )
}))
rownames(summary) = NULL

write.csv(results, file.path(OUTPUT_DIR, "metrics_by_seed.csv"), row.names = FALSE)
write.csv(summary, file.path(OUTPUT_DIR, "summary.csv"), row.names = FALSE)

improvement = -100 * diff(summary$heldout_deviance_mean) /
  head(summary$heldout_deviance_mean, -1)
nmi_gain = diff(summary$lineage_nmi_mean)

best_deviance = min(summary$heldout_deviance_mean)
eligible_K = summary$K[
  summary$heldout_deviance_mean <=
    best_deviance * (1 + RELATIVE_DEVIANCE_TOLERANCE)
]
selected_K = min(eligible_K)

selection = data.frame(
  selected_K = selected_K,
  best_tested_K = summary$K[which.min(summary$heldout_deviance_mean)],
  best_heldout_deviance = best_deviance,
  relative_deviance_tolerance = RELATIVE_DEVIANCE_TOLERANCE,
  rule = "smallest K within the relative tolerance of the best held-out deviance"
)
write.csv(
  selection, file.path(OUTPUT_DIR, "selected_K.csv"), row.names = FALSE
)
saveRDS(selected_K, file.path(OUTPUT_DIR, "selected_K.rds"))

recommendation = if (selected_K == 18) {
  paste0(
    "The held-out rule selects K=18. It is within ",
    100 * RELATIVE_DEVIANCE_TOLERANCE,
    "% of the best tested deviance and is simpler than K=22."
  )
} else {
  paste0("The held-out rule selects K=", selected_K, ".")
}

png(
  file.path(OUTPUT_DIR, "K_sensitivity.png"),
  width = 1500, height = 750, res = 150
)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
plot(
  summary$K, summary$heldout_deviance_mean, type = "b", pch = 16,
  xlab = "K", ylab = "Held-out deviance per entry",
  main = "Predictive fit (lower is better)"
)
plot(
  summary$K, summary$lineage_nmi_mean, type = "b", pch = 16,
  xlab = "K", ylab = "NMI", main = "Lineage agreement"
)
plot(
  summary$K, summary$mean_max_factor_cosine, type = "b", pch = 16,
  ylim = c(0, 1), xlab = "K", ylab = "Maximum factor cosine",
  main = "Factor redundancy"
)
dev.off()

cat("NMF rank sensitivity\n====================\n\n")
print(summary, row.names = FALSE)
cat(
  "\nHeld-out deviance improvement (%):\n",
  paste0(
    "K=", head(summary$K, -1), " to K=", tail(summary$K, -1),
    ": ", sprintf("%.2f", improvement), "%", collapse = "\n"
  ),
  "\n\nSelection rule: choose the smallest K within ",
  100 * RELATIVE_DEVIANCE_TOLERANCE,
  "% of the best held-out deviance.\n",
  recommendation, "\n",
  sep = ""
)
