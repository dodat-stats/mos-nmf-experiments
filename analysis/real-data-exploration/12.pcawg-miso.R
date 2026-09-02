#!/usr/bin/env Rscript

## PCAWG SBS96 pilot for MiSo.
##
## Run from the project root:
##   Rscript --vanilla analysis/real-data-exploration/12.pcawg-miso.R
##
## For a short end-to-end check:
##   Rscript --vanilla analysis/real-data-exploration/12.pcawg-miso.R --quick
##
## The model is fitted only to the observed mutation counts. Published PCAWG
## signatures, exposures, cancer types, and aetiologies are used afterward for
## interpretation and external validation.

options(stringsAsFactors = FALSE)

## -------------------------------------------------------------------------
## 1. Settings: these are the main quantities to edit
## -------------------------------------------------------------------------

QUICK = "--quick" %in% commandArgs(trailingOnly = TRUE)

DATA_SEED = 20260901
TRAIN_FRACTION = 0.80

## Exclude extremely sparse spectra and the most obvious high-burden regime.
## Skin melanoma is set aside because its very large UV burden was analyzed
## separately in the PCAWG signature paper.
MIN_MUTATIONS = 500
MAX_MUTATIONS = 50000
EXCLUDED_CANCER_TYPES = "Skin-Melanoma"

K_CANDIDATES = if (QUICK) c(12, 24, 36) else
  c(8, 12, 16, 20, 24, 28, 32, 36, 40, 48, 56, 64)
K_SELECTION_SEEDS = if (QUICK) 1 else c(1, 2)
K_NMF_ITERS = if (QUICK) 25 else 80
K_LOADING_ITERS = if (QUICK) 30 else 100

## The held-out curve is deliberately reported rather than treated as an
## automatic selector. It continues improving toward the saturated SBS96
## basis in these high-count data. K=36 is a literature-guided pilot value:
## PCAWG reported 31 and 35 signatures in its low-burden extraction.
K = if (QUICK) 24 else 36

S = if (QUICK) 12 else 30
D = 5
FIT_SEEDS = if (QUICK) 1 else c(1, 2)
MF_ITERS = if (QUICK) 4 else 20
MF_NMF_ITERS = if (QUICK) 25 else 80
MISO_ITERS = if (QUICK) 3 else 9
MISO_INNER_ITERS = if (QUICK) 1 else 2
MISO_INITIALIZATION = "distinct"
BLOCK_SIZE = 48

DATA_DIR = "data/pcawg"
OUTPUT_DIR = if (QUICK) "output/pcawg-miso/quick" else "output/pcawg-miso"
CACHE_VERSION = "sbs96-v1"

## -------------------------------------------------------------------------
## 2. Files and dependencies
## -------------------------------------------------------------------------

if (!file.exists("code/miso-benchmark-utils.R")) {
  stop("Run this script from the project root.")
}
source("code/miso-benchmark-utils.R")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

download_if_needed = function(filename, url, md5) {
  destination = file.path(DATA_DIR, filename)
  valid = file.exists(destination) &&
    identical(unname(tools::md5sum(destination)), md5)
  if (!valid) {
    message("Downloading ", filename, "...")
    utils::download.file(url, destination, mode = "wb", quiet = FALSE)
  }
  observed_md5 = unname(tools::md5sum(destination))
  if (!identical(observed_md5, md5)) {
    stop("Checksum failed for ", destination)
  }
  destination
}

github_base = paste0(
  "https://raw.githubusercontent.com/steverozen/PCAWG7/",
  "v0.1.4-branch/"
)

data_files = data.frame(
  filename = c(
    "WGS_PCAWG.96.csv",
    "PCAWG.sample.sheet.rda",
    "PCAWG_sigProfiler_SBS_signatures.csv",
    "PCAWG_sigProfiler_SBS_exposures.csv",
    "COSMIC-v3.2-SBS-aetiology.csv"
  ),
  url = c(
    paste0(github_base, "data-raw/spectra/WGS_PCAWG.96.csv"),
    paste0(github_base, "data/PCAWG.sample.sheet.rda"),
    paste0(
      github_base,
      "data-raw/sig.profiler.signatures/",
      "sigProfiler_SBS_signatures_2019_05_22.csv"
    ),
    paste0(
      github_base,
      "data-raw/sig.profiler.exposures/",
      "PCAWG_sigProfiler_SBS_signatures_in_samples.csv"
    ),
    paste0(
      github_base,
      "data-raw/package-data-related/source-files/",
      "COSMIC-v3.2-SBS-proposed-aetiology-short-form.csv"
    )
  ),
  md5 = c(
    "99e602b3b220e764a0ad7a6abd60fb86",
    "0ac7d9b7556b18619e5fda76e1ab6db4",
    "a306356d26e6429f191ab17a278ce123",
    "fda806b4092040fe78279158f95f131d",
    "f9799f67f619b789f6b33818a4e3152f"
  )
)

paths = mapply(
  download_if_needed,
  data_files$filename,
  data_files$url,
  data_files$md5,
  USE.NAMES = FALSE
)
names(paths) = data_files$filename

## -------------------------------------------------------------------------
## 3. Read and filter the observed SBS96 spectra
## -------------------------------------------------------------------------

spectra = read.csv(paths[["WGS_PCAWG.96.csv"]], check.names = FALSE)
if (nrow(spectra) != 96 || ncol(spectra) < 3) {
  stop("The PCAWG SBS96 file does not have the expected shape.")
}

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

mutation_burden_all = rowSums(Y_all)
keep = mutation_burden_all >= MIN_MUTATIONS &
  mutation_burden_all <= MAX_MUTATIONS &
  !(cancer_type_all %in% EXCLUDED_CANCER_TYPES)

Y = Y_all[keep, , drop = FALSE]
cancer_type = cancer_type_all[keep]
sample_id = sample_id_all[keep]
sample_key = sample_key_all[keep]
mutation_burden = mutation_burden_all[keep]
rm(Y_all)

data_summary = data.frame(
  source_tumors = length(sample_id_all),
  source_cancer_types = length(unique(cancer_type_all)),
  retained_tumors = nrow(Y),
  retained_cancer_types = length(unique(cancer_type)),
  mutation_categories = ncol(Y),
  minimum_burden = min(mutation_burden),
  median_burden = median(mutation_burden),
  maximum_burden = max(mutation_burden),
  excluded_low_burden = sum(mutation_burden_all < MIN_MUTATIONS),
  excluded_high_burden = sum(mutation_burden_all > MAX_MUTATIONS),
  excluded_named_type = sum(cancer_type_all %in% EXCLUDED_CANCER_TYPES)
)
utils::write.csv(
  data_summary, file.path(OUTPUT_DIR, "data-summary.csv"), row.names = FALSE
)

sample_table = data.frame(
  sample_id = sample_id,
  sample_key = sample_key,
  cancer_type = cancer_type,
  mutation_burden = mutation_burden
)
utils::write.csv(
  sample_table, file.path(OUTPUT_DIR, "retained-samples.csv"), row.names = FALSE
)

## -------------------------------------------------------------------------
## 4. Held-out K diagnostic on new tumors
## -------------------------------------------------------------------------

infer_nmf_loadings = function(Y_new, F, max_iters = 100, eps = 1e-12) {
  K_local = nrow(F)
  L = matrix(rowSums(Y_new) / K_local, nrow(Y_new), K_local)
  L = pmax(L, eps)
  F_sum = rowSums(F)
  for (iter in seq_len(max_iters)) {
    mu = L %*% F + eps
    L = L * ((Y_new / mu) %*% t(F)) /
      matrix(F_sum, nrow(Y_new), K_local, byrow = TRUE)
    L = pmax(L, eps)
  }
  L
}

stratified_validation_indices = function(labels, fraction, seed) {
  set.seed(seed)
  unlist(lapply(split(seq_along(labels), labels), function(index) {
    sample(index, max(1, round(fraction * length(index))))
  }), use.names = FALSE)
}

k_cache_dir = file.path(OUTPUT_DIR, paste0("K-selection-", CACHE_VERSION))
dir.create(k_cache_dir, recursive = TRUE, showWarnings = FALSE)
k_rows = list()
k_index = 0L

for (seed in K_SELECTION_SEEDS) {
  validation_index = stratified_validation_indices(
    cancer_type, fraction = 0.20, seed = 2000 + seed
  )
  training_index = setdiff(seq_len(nrow(Y)), validation_index)

  set.seed(3000 + seed)
  Y_validation_fit = matrix(
    rbinom(
      length(Y[validation_index, , drop = FALSE]),
      size = as.integer(Y[validation_index, , drop = FALSE]),
      prob = 0.50
    ),
    nrow = length(validation_index), ncol = ncol(Y)
  )
  Y_validation_test = Y[validation_index, , drop = FALSE] -
    Y_validation_fit

  for (K_candidate in K_CANDIDATES) {
    cache_file = file.path(
      k_cache_dir, paste0("K", K_candidate, "-seed", seed, ".rds")
    )
    if (file.exists(cache_file)) {
      result = readRDS(cache_file)
    } else {
      nmf_fit = poisson_nmf_init(
        Y[training_index, , drop = FALSE],
        K = K_candidate,
        max_iters = K_NMF_ITERS,
        init_seed = seed
      )
      validation_loading = infer_nmf_loadings(
        Y_validation_fit, nmf_fit$F, max_iters = K_LOADING_ITERS
      )
      validation_mean = validation_loading %*% nmf_fit$F
      result = list(
        K = K_candidate,
        seed = seed,
        mean_log_likelihood = mean(dpois(
          Y_validation_test,
          lambda = pmax(validation_mean, 1e-12), log = TRUE
        )),
        deviance_per_entry = poisson_deviance(
          Y_validation_test, validation_mean
        ) / length(Y_validation_test)
      )
      saveRDS(result, cache_file)
    }
    k_index = k_index + 1L
    k_rows[[k_index]] = data.frame(
      K = result$K,
      seed = result$seed,
      mean_log_likelihood = result$mean_log_likelihood,
      deviance_per_entry = result$deviance_per_entry
    )
  }
}

k_results = do.call(rbind, k_rows)
k_summary = do.call(rbind, lapply(split(k_results, k_results$K), function(x) {
  data.frame(
    K = unique(x$K),
    mean_log_likelihood = mean(x$mean_log_likelihood),
    se_log_likelihood = if (nrow(x) > 1) {
      sd(x$mean_log_likelihood) / sqrt(nrow(x))
    } else 0,
    mean_deviance = mean(x$deviance_per_entry)
  )
}))
rownames(k_summary) = NULL

utils::write.csv(
  k_results, file.path(OUTPUT_DIR, "K-selection-runs.csv"), row.names = FALSE
)
utils::write.csv(
  k_summary, file.path(OUTPUT_DIR, "K-selection-summary.csv"), row.names = FALSE
)

if (!(K %in% K_CANDIDATES)) stop("K must be included in K_CANDIDATES.")
message(
  "Using pilot K=", K,
  "; best held-out candidate in the displayed range is K=",
  k_summary$K[which.max(k_summary$mean_log_likelihood)], "."
)

## -------------------------------------------------------------------------
## 5. Fit NMF and MiSo on the same thinned count matrix
## -------------------------------------------------------------------------

set.seed(DATA_SEED)
Y_train = matrix(
  rbinom(length(Y), size = as.integer(Y), prob = TRAIN_FRACTION),
  nrow = nrow(Y), ncol = ncol(Y),
  dimnames = dimnames(Y)
)
Y_test = Y - Y_train
test_scale = (1 - TRAIN_FRACTION) / TRAIN_FRACTION

nmf_file = file.path(
  OUTPUT_DIR, paste0("nmf-", CACHE_VERSION, "-K", K, ".rds")
)
if (file.exists(nmf_file)) {
  nmf_fit = readRDS(nmf_file)
} else {
  message("Fitting the vanilla Poisson NMF baseline...")
  nmf_fit = poisson_nmf_init(
    Y_train, K = K, max_iters = MF_NMF_ITERS, init_seed = DATA_SEED
  )
  saveRDS(nmf_fit, nmf_file)
}

fit_dir = file.path(
  OUTPUT_DIR,
  paste0(
    "fits-", CACHE_VERSION, "-K", K, "-S", S, "-D", D,
    "-init-", MISO_INITIALIZATION
  )
)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

fits = vector("list", length(FIT_SEEDS))
for (j in seq_along(FIT_SEEDS)) {
  seed = FIT_SEEDS[j]
  fit_file = file.path(fit_dir, paste0("fit-seed-", seed, ".rds"))
  if (file.exists(fit_file)) {
    message("Loading cached MiSo seed ", seed, "...")
    fits[[j]] = readRDS(fit_file)
  } else {
    message("Fitting MiSo seed ", seed, "...")
    fits[[j]] = miso(
      Y = Y_train,
      K = K,
      S = S,
      D = D,
      max_iters = MISO_ITERS,
      n_inner = MISO_INNER_ITERS,
      mf_max_iters = MF_ITERS,
      mf_nmf_iters = MF_NMF_ITERS,
      init_seed = seed,
      update_prior = TRUE,
      update_F = TRUE,
      update_gamma = TRUE,
      update_slot_scale = FALSE,
      motif_initialization = MISO_INITIALIZATION,
      tol = 1e-5,
      min_iters = min(5, MISO_ITERS),
      patience = 2,
      block_size = BLOCK_SIZE
    )
    saveRDS(fits[[j]], fit_file)
  }
}

final_objective = vapply(fits, function(current_fit) {
  finite = current_fit$elbo[is.finite(current_fit$elbo)]
  if (length(finite)) tail(finite, 1) else -Inf
}, numeric(1))
best_index = which.max(final_objective)
fit = fits[[best_index]]
best_seed = FIT_SEEDS[best_index]
saveRDS(fit, file.path(OUTPUT_DIR, "miso-fit-best.rds"))

## -------------------------------------------------------------------------
## 6. Posterior summaries and external signature matching
## -------------------------------------------------------------------------

normalize_rows_simple = function(x, eps = 1e-12) {
  x = pmax(x, 0)
  x / pmax(rowSums(x), eps)
}

active_count_90 = function(x) {
  apply(normalize_rows_simple(x), 1, function(row) {
    if (sum(row) <= 1e-12) return(0L)
    as.integer(which(cumsum(sort(row, decreasing = TRUE)) >= 0.90)[1])
  })
}

normalized_mutual_information = function(x, y, eps = 1e-12) {
  tab = table(x, y)
  pxy = tab / sum(tab)
  px = rowSums(pxy)
  py = colSums(pxy)
  expected = px %o% py
  positive = pxy > 0
  mutual_information = sum(
    pxy[positive] * log(pxy[positive] / pmax(expected[positive], eps))
  )
  hx = -sum(px[px > 0] * log(px[px > 0]))
  hy = -sum(py[py > 0] * log(py[py > 0]))
  mutual_information / pmax(sqrt(hx * hy), eps)
}

group_r_squared = function(values, group) {
  total = sum((values - mean(values))^2)
  fitted = ave(values, group, FUN = mean)
  1 - sum((values - fitted)^2) / pmax(total, 1e-12)
}

nmf_scores = nmf_fit$L
miso_scores = miso_observation_factor_scores(fit)
motif_scores = motif_factor_scores(fit)$scores
motif_shares = normalize_rows_simple(motif_scores)
effective_dimension = active_count_90(motif_scores)
dimension_fit = select_motif_dimensions(fit)

reference = read.csv(
  paths[["PCAWG_sigProfiler_SBS_signatures.csv"]], check.names = FALSE
)
reference_context = paste0(
  substr(reference[[2]], 1, 1), "[", reference[[1]], "]",
  substr(reference[[2]], 3, 3)
)
reference_order = match(mutation_context, reference_context)
if (anyNA(reference_order)) stop("Could not align the reference signatures.")
reference_F = t(as.matrix(reference[reference_order, -(1:2), drop = FALSE]))
reference_F = normalize_rows_simple(reference_F)

similarity = row_cosine(fit$F, reference_F)
top_reference_index = max.col(similarity)
top_reference = rownames(reference_F)[top_reference_index]
top_cosine = similarity[cbind(seq_len(K), top_reference_index)]

nmf_similarity = row_cosine(nmf_fit$F, reference_F)
nmf_top_reference_index = max.col(nmf_similarity)
nmf_top_reference = rownames(reference_F)[nmf_top_reference_index]
nmf_top_cosine = nmf_similarity[
  cbind(seq_len(K), nmf_top_reference_index)
]

aetiology = read.csv(
  paths[["COSMIC-v3.2-SBS-aetiology.csv"]], check.names = FALSE
)
top_aetiology = aetiology$proposed.aetiology[
  match(top_reference, aetiology$name)
]
top_aetiology[is.na(top_aetiology) | !nzchar(top_aetiology)] =
  "No established aetiology"

published_exposure = read.csv(
  paths[["PCAWG_sigProfiler_SBS_exposures.csv"]], check.names = FALSE
)
published_order = match(sample_id, published_exposure[["Sample Names"]])
if (anyNA(published_order)) stop("Could not align the published exposures.")
published_exposure = published_exposure[published_order, , drop = FALSE]

exposure_correlation = vapply(seq_len(K), function(k) {
  signature = top_reference[k]
  if (!(signature %in% names(published_exposure))) return(NA_real_)
  suppressWarnings(cor(
    miso_scores[, k], published_exposure[[signature]], method = "spearman"
  ))
}, numeric(1))

factor_usage = colSums(miso_scores)
factor_usage = factor_usage / sum(factor_usage)
factor_type_mass = sapply(seq_len(K), function(k) {
  tapply(miso_scores[, k], cancer_type, sum)
})
if (is.null(dim(factor_type_mass))) {
  factor_type_mass = matrix(factor_type_mass, ncol = K)
}
factor_type_mass = apply(
  factor_type_mass, 2, function(x) x / pmax(sum(x), 1e-12)
)
if (is.null(dim(factor_type_mass))) factor_type_mass = matrix(factor_type_mass)
top_factor_type_index = max.col(t(factor_type_mass))
top_factor_type = rownames(factor_type_mass)[top_factor_type_index]
top_factor_type_share = factor_type_mass[cbind(top_factor_type_index, seq_len(K))]

factor_summary = data.frame(
  factor = paste0("F", seq_len(K)),
  matched_signature = top_reference,
  cosine_similarity = top_cosine,
  published_exposure_spearman = exposure_correlation,
  proposed_aetiology = top_aetiology,
  usage = factor_usage,
  top_cancer_type = top_factor_type,
  top_cancer_type_share = top_factor_type_share
)

## An over-fitted K can split one published signature over several learned
## factors. Aggregate only for external validation; the MiSo fit itself never
## sees these reference labels.
signature_levels = unique(top_reference[order(factor_usage, decreasing = TRUE)])
aggregated_signature_scores = sapply(signature_levels, function(signature) {
  rowSums(miso_scores[, top_reference == signature, drop = FALSE])
})
if (is.null(dim(aggregated_signature_scores))) {
  aggregated_signature_scores = matrix(
    aggregated_signature_scores, ncol = 1,
    dimnames = list(NULL, signature_levels)
  )
}
aggregated_exposure_correlation = vapply(signature_levels, function(signature) {
  if (!(signature %in% names(published_exposure))) return(NA_real_)
  suppressWarnings(cor(
    aggregated_signature_scores[, signature],
    published_exposure[[signature]], method = "spearman"
  ))
}, numeric(1))
aggregated_signature_summary = data.frame(
  matched_signature = signature_levels,
  learned_factors = vapply(signature_levels, function(signature) {
    paste0("F", which(top_reference == signature), collapse = ";")
  }, character(1)),
  number_of_learned_factors = vapply(
    signature_levels, function(signature) sum(top_reference == signature),
    integer(1)
  ),
  total_usage = vapply(signature_levels, function(signature) {
    sum(factor_usage[top_reference == signature])
  }, numeric(1)),
  usage_weighted_cosine_similarity = vapply(
    signature_levels, function(signature) {
      index = top_reference == signature
      sum(factor_usage[index] * top_cosine[index]) /
        pmax(sum(factor_usage[index]), 1e-12)
    }, numeric(1)
  ),
  published_exposure_spearman = aggregated_exposure_correlation,
  proposed_aetiology = aetiology$proposed.aetiology[
    match(signature_levels, aetiology$name)
  ]
)
aggregated_signature_summary$proposed_aetiology[
  is.na(aggregated_signature_summary$proposed_aetiology) |
    !nzchar(aggregated_signature_summary$proposed_aetiology)
] = "No established aetiology"
factor_summary$aggregated_exposure_spearman =
  aggregated_exposure_correlation[match(top_reference, signature_levels)]

## A reference-aggregated motif representation is used only in the external
## interpretation layer. The model was fitted without reference signatures.
motif_signature_shares = sapply(signature_levels, function(signature) {
  rowSums(motif_shares[, top_reference == signature, drop = FALSE])
})
if (is.null(dim(motif_signature_shares))) {
  motif_signature_shares = matrix(
    motif_signature_shares, ncol = 1,
    dimnames = list(NULL, signature_levels)
  )
}
R_signature = ncol(motif_signature_shares)
signature_usage = colSums(
  sweep(motif_signature_shares, 1, fit$pi, "*")
)
utils::write.csv(
  factor_summary, file.path(OUTPUT_DIR, "factor-summary.csv"), row.names = FALSE
)
utils::write.csv(
  aggregated_signature_summary,
  file.path(OUTPUT_DIR, "aggregated-signature-summary.csv"),
  row.names = FALSE
)

## Mixture-weighted cancer composition of each motif.
cancer_levels = names(sort(table(cancer_type), decreasing = TRUE))
motif_cancer_mass = matrix(
  0, nrow = S, ncol = length(cancer_levels),
  dimnames = list(paste0("S", seq_len(S)), cancer_levels)
)
for (g in seq_along(cancer_levels)) {
  motif_cancer_mass[, g] = colSums(
    fit$omega[cancer_type == cancer_levels[g], , drop = FALSE]
  )
}
motif_cancer_share = motif_cancer_mass /
  pmax(rowSums(motif_cancer_mass), 1e-12)

top_factors_text = vapply(seq_len(S), function(s) {
  order_s = order(motif_shares[s, ], decreasing = TRUE)
  active = order_s[seq_len(effective_dimension[s])]
  paste(
    paste0("F", active, "=", top_reference[active]),
    collapse = "; "
  )
}, character(1))

top_cancers_text = vapply(seq_len(S), function(s) {
  top = order(motif_cancer_share[s, ], decreasing = TRUE)[seq_len(3)]
  paste(
    paste0(
      cancer_levels[top], " ",
      sprintf("%.0f%%", 100 * motif_cancer_share[s, top])
    ),
    collapse = "; "
  )
}, character(1))

hard_assignment = fit$z_hat
motif_summary = data.frame(
  submanifold = paste0("S", seq_len(S)),
  mixture_weight = fit$pi,
  effective_dimension_90pct = effective_dimension,
  D_hat_ard_or_lambda = dimension_fit$selection$D_hat_ard_or_lambda,
  top_factors = top_factors_text,
  top_cancer_types = top_cancers_text,
  cancer_type_purity = apply(motif_cancer_share, 1, max),
  assigned_tumors = tabulate(hard_assignment, nbins = S),
  median_mutation_burden = vapply(seq_len(S), function(s) {
    index = hard_assignment == s
    if (any(index)) median(mutation_burden[index]) else NA_real_
  }, numeric(1))
)
utils::write.csv(
  motif_summary, file.path(OUTPUT_DIR, "motif-summary.csv"), row.names = FALSE
)
utils::write.csv(
  motif_shares, file.path(OUTPUT_DIR, "motif-factor-shares.csv")
)
utils::write.csv(
  motif_cancer_share, file.path(OUTPUT_DIR, "motif-cancer-shares.csv")
)

## -------------------------------------------------------------------------
## 7. Held-out evaluation and burden-confounding diagnostic
## -------------------------------------------------------------------------

nmf_test_mean = test_scale * (nmf_scores %*% nmf_fit$F)
miso_test_mean = test_scale * (miso_scores %*% fit$F)

evaluation = data.frame(
  method = c("Vanilla Poisson NMF", "MiSo"),
  heldout_mean_log_likelihood = c(
    mean(dpois(Y_test, pmax(nmf_test_mean, 1e-12), log = TRUE)),
    mean(dpois(Y_test, pmax(miso_test_mean, 1e-12), log = TRUE))
  ),
  heldout_deviance_per_entry = c(
    poisson_deviance(Y_test, nmf_test_mean) / length(Y_test),
    poisson_deviance(Y_test, miso_test_mean) / length(Y_test)
  ),
  active_factors_per_tumor_90pct = c(
    mean(active_count_90(nmf_scores)),
    mean(active_count_90(miso_scores))
  )
)
utils::write.csv(
  evaluation, file.path(OUTPUT_DIR, "goal-evaluation.csv"), row.names = FALSE
)

diagnostics = data.frame(
  best_seed = best_seed,
  mean_max_responsibility = mean(apply(fit$omega, 1, max)),
  median_max_responsibility = median(apply(fit$omega, 1, max)),
  cancer_type_NMI = normalized_mutual_information(
    hard_assignment, cancer_type
  ),
  burden_decile_NMI = normalized_mutual_information(
    hard_assignment,
    cut(
      mutation_burden,
      breaks = unique(quantile(mutation_burden, seq(0, 1, 0.1))),
      include.lowest = TRUE
    )
  ),
  log_burden_R2_by_motif = group_r_squared(
    log10(mutation_burden), hard_assignment
  ),
  median_factor_match_cosine = median(top_cosine),
  factors_with_cosine_at_least_0.85 = sum(top_cosine >= 0.85),
  unique_reference_matches = length(unique(top_reference)),
  nmf_median_factor_match_cosine = median(nmf_top_cosine),
  nmf_factors_with_cosine_at_least_0.85 = sum(nmf_top_cosine >= 0.85),
  nmf_unique_reference_matches = length(unique(nmf_top_reference)),
  median_published_exposure_spearman = median(
    exposure_correlation, na.rm = TRUE
  ),
  median_aggregated_exposure_spearman = median(
    aggregated_exposure_correlation, na.rm = TRUE
  ),
  mean_effective_motif_dimension = sum(fit$pi * effective_dimension)
)
utils::write.csv(
  diagnostics, file.path(OUTPUT_DIR, "diagnostics.csv"), row.names = FALSE
)

## -------------------------------------------------------------------------
## 8. Factor co-use graph
## -------------------------------------------------------------------------

co_use = matrix(0, K, K)
for (s in seq_len(S)) {
  co_use = co_use + fit$pi[s] *
    (motif_shares[s, ] %o% motif_shares[s, ])
}
diag(co_use) = 0
weighted_degree = rowSums(co_use)

upper_index = which(upper.tri(co_use), arr.ind = TRUE)
edge_table = data.frame(
  from = paste0("F", upper_index[, 1]),
  to = paste0("F", upper_index[, 2]),
  weight = co_use[upper_index]
)
edge_table = edge_table[order(edge_table$weight, decreasing = TRUE), ]
edge_table = head(edge_table, min(nrow(edge_table), 2 * K))

node_table = transform(
  factor_summary,
  weighted_degree = weighted_degree,
  graph_role = ifelse(
    weighted_degree >= quantile(weighted_degree, 0.75), "hub",
    ifelse(weighted_degree <= quantile(weighted_degree, 0.25), "leaf", "middle")
  )
)
utils::write.csv(
  node_table, file.path(OUTPUT_DIR, "factor-graph-nodes.csv"), row.names = FALSE
)
utils::write.csv(
  edge_table, file.path(OUTPUT_DIR, "factor-graph-edges.csv"), row.names = FALSE
)

## A reference-aggregated version is easier to read when K is over-specified.
## This is an interpretation layer, not part of model fitting.
signature_co_use = matrix(0, R_signature, R_signature)
for (s in seq_len(S)) {
  signature_co_use = signature_co_use + fit$pi[s] *
    (motif_signature_shares[s, ] %o% motif_signature_shares[s, ])
}
diag(signature_co_use) = 0
signature_degree = rowSums(signature_co_use)

signature_upper = which(upper.tri(signature_co_use), arr.ind = TRUE)
signature_edge_table = data.frame(
  from = signature_levels[signature_upper[, 1]],
  to = signature_levels[signature_upper[, 2]],
  weight = signature_co_use[signature_upper]
)
signature_edge_table = signature_edge_table[
  order(signature_edge_table$weight, decreasing = TRUE),
]
signature_edge_table = head(
  signature_edge_table,
  min(nrow(signature_edge_table), 2 * R_signature)
)
signature_node_table = data.frame(
  signature = signature_levels,
  usage = signature_usage,
  weighted_degree = signature_degree,
  published_exposure_spearman = aggregated_exposure_correlation,
  proposed_aetiology = aggregated_signature_summary$proposed_aetiology
)
utils::write.csv(
  signature_node_table,
  file.path(OUTPUT_DIR, "signature-graph-nodes.csv"), row.names = FALSE
)
utils::write.csv(
  signature_edge_table,
  file.path(OUTPUT_DIR, "signature-graph-edges.csv"), row.names = FALSE
)

## -------------------------------------------------------------------------
## 9. Phone-readable plots
## -------------------------------------------------------------------------

open_phone_png = function(filename, width = 1800, height = 2400) {
  grDevices::png(
    file.path(OUTPUT_DIR, filename), width = width, height = height,
    res = 200, bg = "white"
  )
}

palette_types = setNames(
  grDevices::hcl.colors(length(cancer_levels), "Dynamic"), cancer_levels
)

## Page 1: data and K diagnostic.
open_phone_png("01-data-and-K-diagnostic.png")
graphics::layout(matrix(1:3, ncol = 1), heights = c(0.85, 1.15, 1.15))
graphics::par(mar = c(1, 1, 2, 1))
graphics::plot.new()
graphics::text(
  0, 0.95, "PCAWG SBS96 pilot", adj = c(0, 1), cex = 2.0, font = 2
)
summary_lines = c(
  paste0(
    "Observed input: ", nrow(Y), " tumors × 96 mutation categories; ",
    length(unique(cancer_type)), " cancer types."
  ),
  paste0(
    "Filter: ", MIN_MUTATIONS, "–", format(MAX_MUTATIONS, big.mark = ","),
    " mutations; skin melanoma excluded from this pilot."
  ),
  paste0(
    "Pilot fit: K=", K, ", S=", S, ", D=", D,
    ". Labels and published signatures were not used for fitting."
  ),
  "The K curve is a diagnostic: likelihood alone favors an increasingly saturated mutation-category basis."
)
graphics::text(
  0, 0.73, paste(summary_lines, collapse = "\n"),
  adj = c(0, 1), cex = 1.05
)

graphics::par(mar = c(5, 5, 3, 1))
graphics::hist(
  log10(mutation_burden_all), breaks = 45, col = "#6BAED6",
  border = "white", xlab = "log10 total SBS mutations",
  main = "Mutation burden before filtering"
)
graphics::abline(
  v = log10(c(MIN_MUTATIONS, MAX_MUTATIONS)),
  col = "#B2182B", lty = 2, lwd = 2
)

graphics::par(mar = c(5, 5, 3, 1))
best_candidate = k_summary$K[which.max(k_summary$mean_log_likelihood)]
graphics::plot(
  k_summary$K, k_summary$mean_log_likelihood,
  type = "b", pch = 19, lwd = 2, col = "#2166AC",
  xlab = "Number of vanilla NMF factors K",
  ylab = "Held-out mean Poisson log likelihood",
  main = "Generalization to held-out tumors still favors large K"
)
if (any(k_summary$se_log_likelihood > 0)) {
  graphics::arrows(
    k_summary$K,
    k_summary$mean_log_likelihood - 2 * k_summary$se_log_likelihood,
    k_summary$K,
    k_summary$mean_log_likelihood + 2 * k_summary$se_log_likelihood,
    angle = 90, code = 3, length = 0.04, col = "#2166AC"
  )
}
graphics::abline(v = K, col = "#B2182B", lty = 2, lwd = 2)
graphics::legend(
  "bottomright",
  legend = c(paste0("pilot K = ", K), paste0("best displayed K = ", best_candidate)),
  col = c("#B2182B", "#2166AC"), lty = c(2, 1), lwd = 2, bty = "n"
)
grDevices::dev.off()

## Page 2: predictive fit and reference matches.
open_phone_png("02-fit-and-signature-matches.png")
graphics::layout(matrix(c(1, 2, 3), ncol = 1), heights = c(0.75, 0.85, 1.6))
graphics::par(mar = c(1, 1, 2, 1))
graphics::plot.new()
graphics::text(0, 0.95, "Fit quality and factor interpretation", adj = c(0, 1),
               cex = 1.9, font = 2)
graphics::text(
  0, 0.68,
  paste0(
    "MiSo uses ", sprintf("%.2f", evaluation$active_factors_per_tumor_90pct[2]),
    " factors per tumor (90% mass), versus ",
    sprintf("%.2f", evaluation$active_factors_per_tumor_90pct[1]),
    " for vanilla NMF.\n",
    "Weighted motif dimension: ",
    sprintf("%.2f", diagnostics$mean_effective_motif_dimension),
    "; median factor-to-PCAWG-signature cosine: ",
    sprintf("%.2f", diagnostics$median_factor_match_cosine),
    " (NMF ", sprintf("%.2f", diagnostics$nmf_median_factor_match_cosine),
    ").\nDistinct closest signature matches: ",
    diagnostics$unique_reference_matches, " for MiSo versus ",
    diagnostics$nmf_unique_reference_matches, " for NMF."
  ),
  adj = c(0, 1), cex = 1.05
)

graphics::par(mar = c(5, 5, 3, 1))
bar_values = evaluation$heldout_deviance_per_entry
graphics::barplot(
  bar_values, names.arg = c("NMF", "MiSo"),
  col = c("#969696", "#2166AC"), border = NA,
  ylab = "Held-out deviance per entry",
  main = "Prediction after identical 80/20 Poisson thinning",
  ylim = c(0, max(bar_values) * 1.18)
)
graphics::text(
  x = c(0.7, 1.9), y = bar_values,
  labels = sprintf("%.3f", bar_values), pos = 3, cex = 0.9
)

display_signature = head(
  order(aggregated_signature_summary$total_usage, decreasing = TRUE),
  min(15, nrow(aggregated_signature_summary))
)
display_signature = display_signature[
  order(aggregated_signature_summary$usage_weighted_cosine_similarity[
    display_signature
  ])
]
signature_labels = paste0(
  aggregated_signature_summary$matched_signature[display_signature],
  "  [", aggregated_signature_summary$number_of_learned_factors[
    display_signature
  ], " factors; use ",
  sprintf("%.0f%%", 100 * aggregated_signature_summary$total_usage[
    display_signature
  ]), "; ρ=",
  sprintf("%.2f", aggregated_signature_summary$published_exposure_spearman[
    display_signature
  ]), "]"
)
graphics::par(mar = c(5, 18, 3, 1))
graphics::barplot(
  aggregated_signature_summary$usage_weighted_cosine_similarity[
    display_signature
  ], names.arg = signature_labels,
  horiz = TRUE, las = 1, col = "#67A9CF", border = NA,
  cex.names = 0.78,
  xlim = c(0, 1),
  xlab = "Usage-weighted cosine similarity to published PCAWG signature",
  main = "Overfitted factors grouped by closest reference signature"
)
graphics::abline(v = 0.85, lty = 2, col = "#B2182B")
grDevices::dev.off()

## Page 3: motif atlas after grouping split factors by reference signature.
top_motifs = head(order(fit$pi, decreasing = TRUE), min(20, S))
top_signatures = head(
  order(signature_usage, decreasing = TRUE), min(20, R_signature)
)
atlas = motif_signature_shares[top_motifs, top_signatures, drop = FALSE]
atlas = atlas / pmax(apply(atlas, 1, max), 1e-12)

open_phone_png("03-motif-atlas.png")
graphics::par(mar = c(11, 7, 5, 2))
graphics::image(
  seq_len(ncol(atlas)), seq_len(nrow(atlas)),
  t(atlas[nrow(atlas):1, , drop = FALSE]),
  col = rev(grDevices::hcl.colors(100, "YlOrRd")), axes = FALSE,
  xlab = "", ylab = "", useRaster = TRUE,
  main = "Recurring MiSo motifs: relative mutational-process use"
)
graphics::axis(
  1, at = seq_len(ncol(atlas)),
  labels = signature_levels[top_signatures], las = 2, cex.axis = 0.82
)
graphics::axis(
  2, at = seq_len(nrow(atlas)),
  labels = paste0(
    "S", rev(top_motifs), "  π=",
    sprintf("%.2f", fit$pi[rev(top_motifs)]),
    "  D=", effective_dimension[rev(top_motifs)]
  ),
  las = 1, cex.axis = 0.78
)
graphics::box()
graphics::mtext(
  paste0(
    "Rows show the ", nrow(atlas), " largest motifs; columns show the ",
    ncol(atlas), " matched processes. Reference grouping is for display only."
  ),
  side = 3, line = 1, cex = 0.85
)
grDevices::dev.off()

## Page 4: co-use graph after reference matching.
graph_similarity = signature_co_use /
  pmax(
    sqrt(rowSums(signature_co_use^2) %o% rowSums(signature_co_use^2)),
    1e-12
  )
diag(graph_similarity) = 1
graph_order = stats::hclust(
  stats::as.dist(pmax(1 - graph_similarity, 0)), method = "average"
)$order
angles = seq(
  pi / 2, pi / 2 + 2 * pi, length.out = R_signature + 1
)[-1]
degree_scaled = (signature_degree - min(signature_degree)) /
  pmax(max(signature_degree) - min(signature_degree), 1e-12)
radii = 1 - 0.45 * degree_scaled
xy = matrix(0, R_signature, 2)
xy[graph_order, ] = cbind(
  radii[graph_order] * cos(angles),
  radii[graph_order] * sin(angles)
)

open_phone_png("04-signature-graph.png")
graphics::par(mar = c(2, 2, 5, 2))
graphics::plot(
  xy, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
  xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25),
  main = "Mutational-process relationships induced by MiSo motifs"
)
max_edge = max(signature_edge_table$weight)
for (j in seq_len(nrow(signature_edge_table))) {
  from = match(signature_edge_table$from[j], signature_levels)
  to = match(signature_edge_table$to[j], signature_levels)
  graphics::segments(
    xy[from, 1], xy[from, 2], xy[to, 1], xy[to, 2],
    col = grDevices::adjustcolor("#636363", alpha.f = 0.30),
    lwd = 0.5 + 7 * signature_edge_table$weight[j] / max_edge
  )
}
node_size = 1.2 + 4.0 * sqrt(signature_usage / max(signature_usage))
node_color = grDevices::hcl.colors(100, "Plasma")[
  1 + round(99 * degree_scaled)
]
graphics::points(
  xy, pch = 21, bg = node_color, col = "white", lwd = 1.2,
  cex = node_size
)
graphics::text(
  xy[, 1], xy[, 2], labels = signature_levels,
  cex = 0.68, font = 2
)
graphics::mtext(
  paste0(
    "Learned factors are grouped by their closest PCAWG signature for display ",
    "only; large/inward nodes are common/high-degree."
  ),
  side = 3, line = 1.2, cex = 0.85
)
grDevices::dev.off()

## Page 5: radial and angular coordinates inside four large two-ray motifs.
## The radial coordinate retains loading magnitude. The angular coordinate is
## a display transformation; the MiSo model never constrains L to sum to one.
factor_similarity = row_cosine(fit$F, fit$F)
motif_pairs = lapply(seq_len(S), function(s) {
  ordered = order(motif_shares[s, ], decreasing = TRUE)
  first = ordered[1]
  alternatives = ordered[
    ordered != first &
      motif_shares[s, ordered] >= 0.05 &
      factor_similarity[first, ordered] < 0.95
  ]
  if (!length(alternatives)) return(c(NA_integer_, NA_integer_))
  c(first, alternatives[1])
})
has_pair = vapply(motif_pairs, function(pair) !anyNA(pair), logical(1))
assigned_count = tabulate(hard_assignment, nbins = S)
candidate_motifs = order(
  ifelse(has_pair & assigned_count >= 20, fit$pi, -Inf),
  decreasing = TRUE
)
candidate_motifs = candidate_motifs[
  is.finite(ifelse(
    has_pair[candidate_motifs] & assigned_count[candidate_motifs] >= 20,
    fit$pi[candidate_motifs], -Inf
  ))
]
candidate_motifs = head(candidate_motifs, 4)

open_phone_png("05-continuous-submanifolds.png")
graphics::par(mfrow = c(2, 2), mar = c(5, 5, 4, 1))
for (s in candidate_motifs) {
  factors = motif_pairs[[s]]
  index = hard_assignment == s
  type_in_motif = sort(table(cancer_type[index]), decreasing = TRUE)
  main_types = names(head(type_in_motif, 4))
  color_group = ifelse(cancer_type[index] %in% main_types,
                       cancer_type[index], "Other")
  local_palette = c(palette_types[main_types], Other = "#BDBDBD")
  first_loading = miso_scores[index, factors[1]]
  second_loading = miso_scores[index, factors[2]]
  pair_amplitude = first_loading + second_loading
  x = log10(pair_amplitude + 1)
  y = first_loading / pmax(pair_amplitude, 1e-12)
  graphics::plot(
    x, y, pch = 19,
    col = grDevices::adjustcolor(local_palette[color_group], alpha.f = 0.65),
    xlab = paste0("log10 amplitude: F", factors[1], " + F", factors[2]),
    ylab = paste0(
      "angular coordinate: F", factors[1], "/(F", factors[1],
      "+F", factors[2], ")"
    ),
    main = paste0(
      "S", s, ": n=", sum(index), ", π=", sprintf("%.2f", fit$pi[s]),
      ", angular IQR=", sprintf("%.2f", diff(stats::quantile(y, c(.25, .75))))
    )
  )
  graphics::legend(
    "topright", legend = names(local_palette), col = local_palette,
    pch = 19, cex = 0.62, bty = "n"
  )
}
if (length(candidate_motifs) < 4) {
  for (unused in seq_len(4 - length(candidate_motifs))) graphics::plot.new()
}
grDevices::dev.off()

## Page 6: cancer composition of the largest motifs.
top_cancer_columns = head(seq_along(cancer_levels), min(20, length(cancer_levels)))
cancer_atlas = motif_cancer_share[top_motifs, top_cancer_columns, drop = FALSE]
open_phone_png("06-motif-cancer-map.png")
graphics::par(mar = c(12, 7, 5, 2))
graphics::image(
  seq_len(ncol(cancer_atlas)), seq_len(nrow(cancer_atlas)),
  t(cancer_atlas[nrow(cancer_atlas):1, , drop = FALSE]),
  col = rev(grDevices::hcl.colors(100, "Blues 3")), axes = FALSE,
  xlab = "", ylab = "", useRaster = TRUE,
  main = "Cancer-type composition of the largest motifs"
)
graphics::axis(
  1, at = seq_len(ncol(cancer_atlas)),
  labels = colnames(cancer_atlas), las = 2, cex.axis = 0.72
)
graphics::axis(
  2, at = seq_len(nrow(cancer_atlas)),
  labels = paste0("S", rev(top_motifs)), las = 1, cex.axis = 0.82
)
graphics::box()
graphics::mtext(
  paste0(
    "Cancer-type NMI = ", sprintf("%.2f", diagnostics$cancer_type_NMI),
    "; burden-decile NMI = ", sprintf("%.2f", diagnostics$burden_decile_NMI),
    "."
  ),
  side = 3, line = 1, cex = 0.85
)
grDevices::dev.off()

## A short machine-generated handoff beside the figures.
summary_text = c(
  "PCAWG SBS96 MiSo pilot",
  "",
  paste0("Data: ", nrow(Y), " retained tumors x 96 SBS categories."),
  paste0("Fit: K=", K, ", S=", S, ", D=", D, ", best seed=", best_seed, "."),
  paste0(
    "Held-out deviance per entry: NMF=",
    sprintf("%.4f", evaluation$heldout_deviance_per_entry[1]),
    ", MiSo=", sprintf("%.4f", evaluation$heldout_deviance_per_entry[2]), "."
  ),
  paste0(
    "Active factors per tumor (90% mass): NMF=",
    sprintf("%.2f", evaluation$active_factors_per_tumor_90pct[1]),
    ", MiSo=", sprintf("%.2f", evaluation$active_factors_per_tumor_90pct[2]),
    "."
  ),
  paste0(
    "Median learned-factor/PCAWG-signature cosine=",
    sprintf("%.3f", diagnostics$median_factor_match_cosine), "."
  ),
  paste0(
    "Cancer-type NMI=", sprintf("%.3f", diagnostics$cancer_type_NMI),
    "; burden-decile NMI=", sprintf("%.3f", diagnostics$burden_decile_NMI),
    "; log-burden R2 by motif=",
    sprintf("%.3f", diagnostics$log_burden_R2_by_motif), "."
  ),
  "",
  "Important: held-out likelihood does not give a finite automatic K here;",
  paste0("the near-saturated SBS basis continues to improve prediction. K=", K, " is a"),
  "transparent pilot choice, not a finalized model-selection result."
)
writeLines(summary_text, file.path(OUTPUT_DIR, "README.txt"))

message("Finished. Results are in ", OUTPUT_DIR)
