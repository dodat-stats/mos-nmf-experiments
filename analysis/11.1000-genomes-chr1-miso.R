#!/usr/bin/env Rscript

## Poisson MiSo analysis of 1000 Genomes phase-3 chromosome 1 genotypes.
##
## Run from the project root:
##   Rscript --vanilla analysis/11.1000-genomes-chr1-miso.R
##
## Fast end-to-end check:
##   Rscript --vanilla analysis/11.1000-genomes-chr1-miso.R --quick
##
## Population and super-population labels are never used for fitting. They are
## used only after fitting to interpret the learned factors and submanifolds.

options(stringsAsFactors = FALSE)

## -------------------------------------------------------------------------
## 1. Settings
## -------------------------------------------------------------------------

QUICK = "--quick" %in% commandArgs(trailingOnly = TRUE)

DATA_SEED = 20260831
MAX_PER_POPULATION = if (QUICK) 10 else 40
TRAIN_FRACTION = 0.80
K_CANDIDATES = if (QUICK) c(4, 6) else c(4, 6, 8, 10, 12)
K_SELECTION_SEEDS = if (QUICK) 1 else c(1, 2)
K_NMF_ITERS = if (QUICK) 8 else 100

S = if (QUICK) 8 else 18
D = if (QUICK) 3 else 4
FIT_SEEDS = if (QUICK) 1 else c(1, 2)
MF_ITERS = if (QUICK) 3 else 30
MF_NMF_ITERS = if (QUICK) 8 else 100
MISO_ITERS = if (QUICK) 2 else 10
MISO_INNER_ITERS = if (QUICK) 1 else 2
MISO_INITIALIZATION = "distinct"
BLOCK_SIZE = 128

DATA_DIR = "data/1000-genomes/phase3-chr1"
OUTPUT_DIR = "output/1000-genomes-chr1-miso"
CACHE_VERSION = "poisson-genotype-v3-converged-nmf"
DOSAGE_FILE = file.path(
  DATA_DIR, "chr1.phase3.maf05.thin250kb.dosage.tsv.gz"
)
PANEL_FILE = file.path(
  DATA_DIR, "integrated_call_samples_v3.20130502.ALL.panel"
)
RELATED_FILE = file.path(DATA_DIR, "20140625_related_individuals.txt")

SUPER_POP_ORDER = c("AFR", "AMR", "EAS", "EUR", "SAS")
SUPER_POP_COLORS = c(
  AFR = "#D55E00", AMR = "#CC79A7", EAS = "#009E73",
  EUR = "#0072B2", SAS = "#E69F00"
)

## -------------------------------------------------------------------------
## 2. Files and preprocessing
## -------------------------------------------------------------------------

if (!file.exists("code/miso-benchmark-utils.R")) {
  stop("Run this script from the project root.")
}
source("code/miso-benchmark-utils.R")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(DOSAGE_FILE)) {
  message("Preparing chromosome-1 dosage table...")
  status = system2("bash", "analysis/helpers/prepare-1000g-chr1.sh")
  if (status != 0 || !file.exists(DOSAGE_FILE)) {
    stop("Chromosome-1 preprocessing failed")
  }
}

message("Reading physically thinned chromosome-1 dosages...")
dosage_table = utils::read.delim(
  gzfile(DOSAGE_FILE), header = TRUE, check.names = FALSE,
  quote = "", comment.char = "", na.strings = "NA"
)
variant_info = dosage_table[, 1:6]
sample_ids = names(dosage_table)[-(1:6)]
Y_all = t(as.matrix(dosage_table[, -(1:6), drop = FALSE]))
storage.mode(Y_all) = "numeric"
rownames(Y_all) = sample_ids
colnames(Y_all) = variant_info$variant
rm(dosage_table)

panel = utils::read.delim(PANEL_FILE, check.names = FALSE)
related = utils::read.delim(RELATED_FILE, check.names = FALSE)
names(related)[1] = "sample"
panel = panel[panel$sample %in% sample_ids, ]
panel = panel[!panel$sample %in% related$sample, ]
panel = panel[match(rownames(Y_all), panel$sample), ]
eligible = which(!is.na(panel$sample))

## Use a balanced subset so the largest sampled populations do not dominate
## the unsupervised fit. Labels determine only how many samples are retained;
## they never enter the likelihood or the MiSo updates.
set.seed(DATA_SEED)
chosen = unlist(lapply(sort(unique(panel$pop[eligible])), function(population) {
  candidates = eligible[panel$pop[eligible] == population]
  sample(candidates, min(length(candidates), MAX_PER_POPULATION))
}))
chosen = sample(chosen)

Y = Y_all[chosen, , drop = FALSE]
meta = panel[chosen, , drop = FALSE]
rownames(meta) = meta$sample
rm(Y_all, panel)

missing_rate = mean(is.na(Y))
if (missing_rate > 0) {
  snp_mean = colMeans(Y, na.rm = TRUE)
  missing = which(is.na(Y), arr.ind = TRUE)
  Y[missing] = snp_mean[missing[, 2]]
  ## The current Poisson-thinning implementation requires integer counts.
  Y = round(Y)
}

meta$super_pop = factor(meta$super_pop, levels = SUPER_POP_ORDER)
population_order = unique(
  meta[order(meta$super_pop, meta$pop), c("pop", "super_pop")]
)$pop
meta$pop = factor(meta$pop, levels = population_order)

data_summary = data.frame(
  individuals = nrow(Y),
  populations = nlevels(meta$pop),
  super_populations = nlevels(meta$super_pop),
  chromosome = 1,
  retained_snps = ncol(Y),
  minimum_maf = min(pmin(variant_info$alt_frequency,
                         1 - variant_info$alt_frequency)),
  physical_bin_bp = 250000,
  missing_genotype_fraction = missing_rate,
  mean_alt_alleles_per_individual = mean(rowSums(Y))
)
utils::write.csv(data_summary, file.path(OUTPUT_DIR, "data-summary.csv"),
                 row.names = FALSE)
utils::write.csv(meta, file.path(OUTPUT_DIR, "sample-metadata.csv"),
                 row.names = FALSE)
utils::write.csv(variant_info, file.path(OUTPUT_DIR, "variant-metadata.csv"),
                 row.names = FALSE)

## PCA is used only as a familiar external summary of genetic similarity.
scaled_Y = scale(Y)
scaled_Y[, !is.finite(colSums(scaled_Y))] = 0
pca = stats::prcomp(scaled_Y, center = FALSE, scale. = FALSE, rank. = 2)
pca_scores = pca$x[, 1:2, drop = FALSE]

## -------------------------------------------------------------------------
## 3. Select K with Poisson thinning under vanilla NMF
## -------------------------------------------------------------------------

k_cache_dir = file.path(
  OUTPUT_DIR,
  sprintf("K-selection-%s-N%d-M%d", CACHE_VERSION, nrow(Y), ncol(Y))
)
dir.create(k_cache_dir, recursive = TRUE, showWarnings = FALSE)
k_rows = list()
k_index = 0L

for (seed in K_SELECTION_SEEDS) {
  set.seed(DATA_SEED + seed)
  Y_train = matrix(
    stats::rbinom(length(Y), size = as.integer(Y), prob = TRAIN_FRACTION),
    nrow = nrow(Y), ncol = ncol(Y)
  )
  Y_test = Y - Y_train

  for (K in K_CANDIDATES) {
    cache_file = file.path(k_cache_dir, sprintf("K%02d-seed%02d.rds", K, seed))
    if (file.exists(cache_file)) {
      result = readRDS(cache_file)
    } else {
      message("Selecting K: fitting K=", K, ", seed=", seed)
      nmf_fit = poisson_nmf_init(
        Y_train, K = K, max_iters = K_NMF_ITERS, init_seed = seed
      )
      test_mean = (1 - TRAIN_FRACTION) / TRAIN_FRACTION *
        (nmf_fit$L %*% nmf_fit$F)
      result = list(
        K = K,
        seed = seed,
        heldout_deviance_per_entry =
          poisson_deviance(Y_test, test_mean) / length(Y_test),
        heldout_mean_log_likelihood = mean(stats::dpois(
          Y_test, lambda = pmax(test_mean, 1e-12), log = TRUE
        ))
      )
      saveRDS(result, cache_file)
    }
    k_index = k_index + 1L
    k_rows[[k_index]] = as.data.frame(result)
  }
}

k_results = do.call(rbind, k_rows)
k_summary = do.call(rbind, lapply(split(k_results, k_results$K), function(x) {
  data.frame(
    K = unique(x$K),
    heldout_deviance_mean = mean(x$heldout_deviance_per_entry),
    heldout_deviance_se = if (nrow(x) > 1) {
      stats::sd(x$heldout_deviance_per_entry) / sqrt(nrow(x))
    } else 0,
    heldout_log_likelihood_mean = mean(x$heldout_mean_log_likelihood)
  )
}))
rownames(k_summary) = NULL
k_summary = k_summary[order(k_summary$K), ]

selected_K = k_summary$K[which.min(k_summary$heldout_deviance_mean)]
K_AT_BOUNDARY = selected_K == max(K_CANDIDATES)

utils::write.csv(k_results, file.path(OUTPUT_DIR, "K-selection-runs.csv"),
                 row.names = FALSE)
utils::write.csv(k_summary, file.path(OUTPUT_DIR, "K-selection-summary.csv"),
                 row.names = FALSE)
saveRDS(selected_K, file.path(OUTPUT_DIR, "selected-K.rds"))
message("Selected K = ", selected_K)

## -------------------------------------------------------------------------
## 4. Fit full-data Poisson NMF and MiSo
## -------------------------------------------------------------------------

nmf_cache = file.path(
  OUTPUT_DIR,
  sprintf("nmf-%s-N%d-M%d-K%d.rds", CACHE_VERSION, nrow(Y), ncol(Y), selected_K)
)
if (file.exists(nmf_cache)) {
  message("Loading cached full-data Poisson NMF...")
  nmf_fit = readRDS(nmf_cache)
} else {
  message("Fitting full-data Poisson NMF...")
  nmf_fit = poisson_nmf_init(
    Y, K = selected_K, max_iters = max(K_NMF_ITERS, MF_NMF_ITERS),
    init_seed = DATA_SEED
  )
  saveRDS(nmf_fit, nmf_cache)
}

fit_dir = file.path(
  OUTPUT_DIR,
  sprintf(
    "fits-%s-N%d-M%d-K%d-S%d-D%d-init-%s", CACHE_VERSION,
    nrow(Y), ncol(Y), selected_K, S, D, MISO_INITIALIZATION
  )
)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
fits = vector("list", length(FIT_SEEDS))
fit_deviance = numeric(length(FIT_SEEDS))

for (j in seq_along(FIT_SEEDS)) {
  seed = FIT_SEEDS[j]
  cache_file = file.path(fit_dir, sprintf("seed-%03d.rds", seed))
  if (file.exists(cache_file)) {
    message("Loading cached MiSo seed ", seed, "...")
    fits[[j]] = readRDS(cache_file)
  } else {
    message("Fitting MiSo seed ", seed, "...")
    fits[[j]] = miso(
      Y = Y,
      K = selected_K,
      S = S,
      D = min(D, selected_K),
      max_iters = MISO_ITERS,
      n_inner = MISO_INNER_ITERS,
      mf_max_iters = MF_ITERS,
      mf_nmf_iters = MF_NMF_ITERS,
      init_seed = seed,
      update_prior = TRUE,
      update_F = TRUE,
      update_gamma = TRUE,
      update_slot_scale = FALSE,
      gamma_init_floor = 0.05,
      motif_initialization = MISO_INITIALIZATION,
      gamma_step_init = 0.5,
      gamma_step_ramp = 8,
      F_step_init = 0.2,
      F_step_ramp = 12,
      tol = 1e-5,
      min_iters = min(5, MISO_ITERS),
      patience = 2,
      block_size = BLOCK_SIZE
    )
    saveRDS(fits[[j]], cache_file)
  }
  score = miso_observation_factor_scores(fits[[j]])
  fit_deviance[j] = poisson_deviance(Y, score %*% fits[[j]]$F) / length(Y)
}

best_fit_index = which.min(fit_deviance)
fit = fits[[best_fit_index]]
best_seed = FIT_SEEDS[best_fit_index]
saveRDS(fit, file.path(OUTPUT_DIR, "miso-fit-best.rds"))

## -------------------------------------------------------------------------
## 5. Summaries and interpretation
## -------------------------------------------------------------------------

nmf_scores = nmf_fit$L
miso_scores = miso_observation_factor_scores(fit)
nmf_share = normalize_scores(nmf_scores)
miso_share = normalize_scores(miso_scores)
motif_scores = motif_factor_scores(fit)$scores
motif_share = normalize_scores(motif_scores)

group_mean = function(x, group) {
  total = rowsum(x, group, reorder = FALSE)
  total / as.vector(table(factor(group, levels = rownames(total))))
}

pop_levels = levels(meta$pop)
pop_nmf = group_mean(nmf_scores, meta$pop)[pop_levels, , drop = FALSE]
pop_miso = group_mean(miso_scores, meta$pop)[pop_levels, , drop = FALSE]
pop_omega = group_mean(fit$omega, meta$pop)[pop_levels, , drop = FALSE]
pop_super = as.character(meta$super_pop[match(pop_levels, meta$pop)])

colnames(pop_nmf) = colnames(pop_miso) = paste0("F", seq_len(selected_K))
colnames(pop_omega) = paste0("S", seq_len(S))
utils::write.csv(pop_nmf, file.path(OUTPUT_DIR, "population-nmf-loadings.csv"))
utils::write.csv(pop_miso, file.path(OUTPUT_DIR, "population-miso-loadings.csv"))
utils::write.csv(pop_omega,
                 file.path(OUTPUT_DIR, "population-submanifold-responsibility.csv"))
utils::write.csv(motif_share,
                 file.path(OUTPUT_DIR, "submanifold-factor-share.csv"))

population_given_factor = matrix(0, nrow = length(pop_levels),
                                 ncol = selected_K)
for (p in seq_along(pop_levels)) {
  population_given_factor[p, ] = colSums(
    miso_scores[meta$pop == pop_levels[p], , drop = FALSE]
  )
}
population_given_factor = sweep(
  population_given_factor, 2, colSums(population_given_factor), "/"
)
population_given_factor[!is.finite(population_given_factor)] = 0
rownames(population_given_factor) = pop_levels
colnames(population_given_factor) = paste0("F", seq_len(selected_K))

population_given_motif = matrix(0, nrow = length(pop_levels), ncol = S)
for (p in seq_along(pop_levels)) {
  population_given_motif[p, ] = colSums(
    fit$omega[meta$pop == pop_levels[p], , drop = FALSE]
  )
}
population_given_motif = sweep(
  population_given_motif, 2, colSums(population_given_motif), "/"
)
population_given_motif[!is.finite(population_given_motif)] = 0
rownames(population_given_motif) = pop_levels
colnames(population_given_motif) = paste0("S", seq_len(S))

effective_dimension = apply(motif_share, 1, function(x) {
  which(cumsum(sort(x, decreasing = TRUE)) >= 0.90)[1]
})
ACTIVE_MOTIF_THRESHOLD = 0.01
active_motifs = which(fit$pi >= ACTIVE_MOTIF_THRESHOLD)
if (length(active_motifs) < 2) {
  active_motifs = order(fit$pi, decreasing = TRUE)[seq_len(min(2, S))]
}

## Compare population geometry without normalizing each loading row. Columns
## are standardized only to put learned factor coordinates on equal scales.
population_genotype = group_mean(scaled_Y, meta$pop)[pop_levels, , drop = FALSE]
genetic_distance = stats::dist(population_genotype)
standardize_coordinates = function(x) {
  x = scale(x)
  x[!is.finite(x)] = 0
  x
}
nmf_distance = stats::dist(standardize_coordinates(pop_nmf))
miso_distance = stats::dist(standardize_coordinates(pop_miso))

distance_correlations = data.frame(
  representation = c("Vanilla Poisson NMF", "MiSo"),
  spearman_with_genotype_distance = c(
    stats::cor(as.vector(genetic_distance), as.vector(nmf_distance),
               method = "spearman"),
    stats::cor(as.vector(genetic_distance), as.vector(miso_distance),
               method = "spearman")
  )
)

pair_distance = function(distance, first, second) {
  matrix_distance = as.matrix(distance)
  matrix_distance[first, second]
}

cdx_jpt = data.frame(
  representation = c("Genotype", "Vanilla Poisson NMF", "MiSo"),
  CDX_JPT_distance = c(
    pair_distance(genetic_distance, "CDX", "JPT"),
    pair_distance(nmf_distance, "CDX", "JPT"),
    pair_distance(miso_distance, "CDX", "JPT")
  )
)

nearest_populations = function(distance, population, n = 4) {
  d = as.matrix(distance)[population, ]
  d = sort(d[names(d) != population])
  paste(names(head(d, n)), collapse = ", ")
}

nmf_mean = nmf_scores %*% nmf_fit$F
miso_mean = miso_scores %*% fit$F
reference_center = attr(scaled_Y, "scaled:center")
reference_scale = attr(scaled_Y, "scaled:scale")
scale_like_genotypes = function(x) {
  x = sweep(x, 2, reference_center, "-")
  x = sweep(x, 2, pmax(reference_scale, 1e-12), "/")
  x[!is.finite(x)] = 0
  x
}
nmf_fitted_population_distance = stats::dist(
  group_mean(scale_like_genotypes(nmf_mean), meta$pop)[pop_levels, , drop = FALSE]
)
miso_fitted_population_distance = stats::dist(
  group_mean(scale_like_genotypes(miso_mean), meta$pop)[pop_levels, , drop = FALSE]
)
evaluation = data.frame(
  representation = c("Vanilla Poisson NMF", "MiSo"),
  in_sample_deviance_per_genotype = c(
    poisson_deviance(Y, nmf_mean) / length(Y),
    poisson_deviance(Y, miso_mean) / length(Y)
  ),
  population_distance_correlation =
    distance_correlations$spearman_with_genotype_distance,
  fitted_mean_distance_correlation = c(
    stats::cor(as.vector(genetic_distance),
               as.vector(nmf_fitted_population_distance), method = "spearman"),
    stats::cor(as.vector(genetic_distance),
               as.vector(miso_fitted_population_distance), method = "spearman")
  ),
  average_active_factors_per_individual_90pct = c(
    mean(apply(nmf_share, 1, function(x) {
      which(cumsum(sort(x, decreasing = TRUE)) >= 0.90)[1]
    })),
    mean(apply(miso_share, 1, function(x) {
      which(cumsum(sort(x, decreasing = TRUE)) >= 0.90)[1]
    }))
  ),
  average_submanifold_dimension_90pct = c(
    NA_real_, sum(fit$pi * effective_dimension)
  )
)

utils::write.csv(evaluation, file.path(OUTPUT_DIR, "goal-evaluation.csv"),
                 row.names = FALSE)
utils::write.csv(distance_correlations,
                 file.path(OUTPUT_DIR, "population-distance-correlations.csv"),
                 row.names = FALSE)
utils::write.csv(cdx_jpt, file.path(OUTPUT_DIR, "CDX-JPT-distance.csv"),
                 row.names = FALSE)
utils::write.csv(population_given_factor,
                 file.path(OUTPUT_DIR, "population-given-factor.csv"))
utils::write.csv(population_given_motif,
                 file.path(OUTPUT_DIR, "population-given-submanifold.csv"))

## -------------------------------------------------------------------------
## 6. Wasserstein submanifold tree and factor co-use graph
## -------------------------------------------------------------------------

row_cosine_similarity = function(x) {
  norm = sqrt(rowSums(x^2))
  similarity = (x %*% t(x)) / pmax(norm %o% norm, 1e-12)
  pmin(pmax(similarity, 0), 1)
}

sinkhorn_distance = function(p, q, cost, epsilon = 0.05,
                              max_iters = 300, tol = 1e-8) {
  p = pmax(p, 1e-12); p = p / sum(p)
  q = pmax(q, 1e-12); q = q / sum(q)
  kernel = pmax(exp(-cost / epsilon), 1e-300)
  u = rep(1, length(p)); v = rep(1, length(q))
  for (iteration in seq_len(max_iters)) {
    old_u = u
    u = p / pmax(as.vector(kernel %*% v), 1e-300)
    v = q / pmax(as.vector(t(kernel) %*% u), 1e-300)
    if (max(abs(u - old_u) / pmax(abs(old_u), 1)) < tol) break
  }
  sum(((u %o% v) * kernel) * cost)
}

factor_cost = 1 - row_cosine_similarity(fit$F)
submanifold_distance = matrix(0, S, S)
for (s1 in seq_len(S - 1)) {
  for (s2 in (s1 + 1):S) {
    submanifold_distance[s1, s2] = submanifold_distance[s2, s1] =
      0.5 * (
        sinkhorn_distance(motif_share[s1, ], motif_share[s2, ], factor_cost) +
        sinkhorn_distance(motif_share[s2, ], motif_share[s1, ], factor_cost)
      )
  }
}
rownames(submanifold_distance) = colnames(submanifold_distance) =
  paste0("S", seq_len(S))
submanifold_tree = stats::hclust(
  stats::as.dist(
    submanifold_distance[active_motifs, active_motifs, drop = FALSE]
  ),
  method = "average"
)
## The fitted dendrogram has a stable gap between 0.014 and 0.030, so 0.020
## gives a conservative post-fit merge without joining the larger branches.
MERGE_DISTANCE = 0.020
merge_assignment = stats::cutree(submanifold_tree, h = MERGE_DISTANCE)
merge_table = data.frame(
  motif = names(merge_assignment),
  merge_group = unname(merge_assignment),
  mixture_weight = fit$pi[active_motifs]
)
merge_table$group_weight = ave(
  merge_table$mixture_weight, merge_table$merge_group, FUN = sum
)
utils::write.csv(merge_table,
                 file.path(OUTPUT_DIR, "submanifold-merge-groups.csv"),
                 row.names = FALSE)
utils::write.csv(submanifold_distance,
                 file.path(OUTPUT_DIR, "submanifold-wasserstein-distance.csv"))

co_use = matrix(0, selected_K, selected_K)
for (s in seq_len(S)) {
  co_use = co_use + fit$pi[s] * (motif_share[s, ] %o% motif_share[s, ])
}
diag(co_use) = 0
weighted_degree = rowSums(co_use)
factor_usage = colSums(sweep(motif_share, 1, fit$pi, "*"))

upper_index = which(upper.tri(co_use) & co_use > 0, arr.ind = TRUE)
edge_table = data.frame(
  from = upper_index[, 1], to = upper_index[, 2], weight = co_use[upper_index]
)
edge_table = head(edge_table[order(edge_table$weight, decreasing = TRUE), ],
                  2 * selected_K)

factor_breadth = apply(population_given_factor, 2, function(x) {
  x = x[x > 0]
  if (length(x) < 2) return(0)
  -sum(x * log(x)) / log(length(pop_levels))
})
factor_node_table = data.frame(
  factor = paste0("F", seq_len(selected_K)),
  top_population = pop_levels[apply(population_given_factor, 2, which.max)],
  usage = factor_usage,
  weighted_degree = weighted_degree,
  population_breadth = factor_breadth
)
utils::write.csv(factor_node_table,
                 file.path(OUTPUT_DIR, "factor-graph-nodes.csv"),
                 row.names = FALSE)
utils::write.csv(edge_table, file.path(OUTPUT_DIR, "factor-graph-edges.csv"),
                 row.names = FALSE)

graph_similarity = co_use / pmax(
  sqrt(rowSums(co_use^2) %o% rowSums(co_use^2)), 1e-12
)
diag(graph_similarity) = 1
graph_order = stats::hclust(
  stats::as.dist(pmax(1 - graph_similarity, 0)), method = "average"
)$order
graph_angle = seq(pi / 2, pi / 2 + 2 * pi,
                  length.out = selected_K + 1)[-1]
degree_scaled = (weighted_degree - min(weighted_degree)) /
  pmax(max(weighted_degree) - min(weighted_degree), 1e-12)
graph_radius = 1 - 0.45 * degree_scaled
graph_xy = matrix(0, selected_K, 2)
graph_xy[graph_order, ] = cbind(
  graph_radius[graph_order] * cos(graph_angle),
  graph_radius[graph_order] * sin(graph_angle)
)

## -------------------------------------------------------------------------
## 7. Phone-readable figures
## -------------------------------------------------------------------------

factor_colors = grDevices::hcl.colors(selected_K, "Dark 3")
motif_colors = grDevices::hcl.colors(S, "Dynamic")

phone_png = function(filename, width, height, expression) {
  grDevices::png(file.path(OUTPUT_DIR, filename), width = width,
                 height = height, res = 150)
  on.exit(grDevices::dev.off())
  force(expression)
}

plot_k_selection = function() {
  graphics::par(mar = c(4.5, 5, 3.5, 1))
  graphics::plot(
    k_summary$K, k_summary$heldout_deviance_mean, type = "b", pch = 16,
    lwd = 2, xlab = "Number of factors K",
    ylab = "Held-out Poisson deviance per genotype",
    main = paste0("Held-out selection chooses K = ", selected_K)
  )
  if (any(k_summary$heldout_deviance_se > 0)) {
    graphics::arrows(
      k_summary$K,
      k_summary$heldout_deviance_mean - k_summary$heldout_deviance_se,
      k_summary$K,
      k_summary$heldout_deviance_mean + k_summary$heldout_deviance_se,
      angle = 90, code = 3, length = 0.04
    )
  }
  graphics::abline(v = selected_K, col = "firebrick", lty = 2)
}

plot_pca = function() {
  graphics::par(mar = c(4.5, 4.5, 3, 1))
  graphics::plot(
    pca_scores, pch = 16, cex = 0.65,
    col = SUPER_POP_COLORS[as.character(meta$super_pop)],
    xlab = "Genotype PC1", ylab = "Genotype PC2",
    main = "Population structure in the thinned chromosome-1 data"
  )
  centroid = aggregate(pca_scores, list(pop = meta$pop), mean)
  graphics::text(centroid[, 2], centroid[, 3], labels = centroid$pop,
                 cex = 0.65, font = 2)
  graphics::legend("topright", legend = SUPER_POP_ORDER,
                   col = SUPER_POP_COLORS, pch = 16, bty = "n")
}

loading_ymax = 1.03 * max(rowSums(nmf_scores), rowSums(miso_scores))

plot_loading_panels = function(scores, title) {
  graphics::par(mfrow = c(length(SUPER_POP_ORDER), 1),
                mar = c(3.2, 4.2, 2.1, 0.5), oma = c(0, 0, 3, 0))
  for (super_population in SUPER_POP_ORDER) {
    index = which(meta$super_pop == super_population)
    local_pop = droplevels(meta$pop[index])
    local_pc = pca_scores[index, 1]
    index = index[order(local_pop, local_pc)]
    local_pop = droplevels(meta$pop[index])
    graphics::barplot(
      t(scores[index, , drop = FALSE]), col = factor_colors, border = NA,
      space = 0, axes = FALSE, ylim = c(0, loading_ymax), xlab = ""
    )
    graphics::axis(2, las = 1, cex.axis = 0.7)
    runs = rle(as.character(local_pop))
    ends = cumsum(runs$lengths)
    starts = c(1, head(ends, -1) + 1)
    graphics::abline(v = head(ends, -1), col = "white", lwd = 0.7)
    graphics::axis(1, at = (starts + ends) / 2, labels = runs$values,
                   tick = FALSE, line = -0.5, cex.axis = 0.72)
    graphics::mtext(super_population, side = 3, line = 0.2, adj = 0,
                    font = 2, col = SUPER_POP_COLORS[super_population])
  }
  graphics::mtext(title, outer = TRUE, side = 3, line = 1, font = 2, cex = 1.15)
}

plot_population_loadings = function() {
  graphics::par(mar = c(7, 4.5, 3, 1))
  midpoint = graphics::barplot(
    t(pop_miso), col = factor_colors, border = "white", space = 0.2,
    names.arg = pop_levels, las = 2, cex.names = 0.75,
    ylab = "Mean unnormalized loading",
    main = "MiSo reveals shared genetic parts across populations"
  )
  boundaries = which(
    head(pop_super, -1) != tail(pop_super, -1)
  )
  if (length(boundaries)) {
    graphics::abline(v = (midpoint[boundaries] + midpoint[boundaries + 1]) / 2,
                     col = "grey60", lty = 2)
  }
  graphics::legend("topleft", legend = paste0("F", seq_len(selected_K)),
                   fill = factor_colors, ncol = ceiling(selected_K / 4),
                   cex = 0.7, bty = "n")
}

plot_population_motifs = function() {
  graphics::par(mar = c(6.5, 5, 3, 2))
  displayed_motifs = active_motifs
  graphics::image(
    seq_along(pop_levels), seq_along(displayed_motifs),
    pop_omega[, displayed_motifs, drop = FALSE],
    col = grDevices::hcl.colors(100, "YlOrRd", rev = TRUE), axes = FALSE,
    xlab = "", ylab = "", main = "Mean submanifold responsibility by population"
  )
  graphics::axis(1, at = seq_along(pop_levels), labels = pop_levels,
                 las = 2, cex.axis = 0.7)
  graphics::axis(2, at = seq_along(displayed_motifs),
                 labels = paste0("S", displayed_motifs),
                 las = 1, cex.axis = 0.75)
  graphics::box()
}

plot_submanifold_tree = function() {
  displayed_motifs = active_motifs
  top_pop = apply(
    population_given_motif[, displayed_motifs, drop = FALSE], 2, function(x) {
    paste(pop_levels[order(x, decreasing = TRUE)[1:2]], collapse = "/")
  })
  labels = paste0(
    "S", displayed_motifs, " ", top_pop, " D",
    effective_dimension[displayed_motifs]
  )
  tree = submanifold_tree
  tree$labels = labels
  graphics::par(mar = c(4, 2, 3, 18), xpd = NA)
  graphics::plot(stats::as.dendrogram(tree), horiz = TRUE,
                 main = "Wasserstein tree of genetic submanifolds",
                 xlab = "Submanifold distance")
}

plot_factor_graph = function() {
  graphics::par(mar = c(1, 1, 3, 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35), asp = 1)
  if (nrow(edge_table)) {
    edge_scale = edge_table$weight / max(edge_table$weight)
    for (e in seq_len(nrow(edge_table))) {
      from = edge_table$from[e]; to = edge_table$to[e]
      graphics::segments(
        graph_xy[from, 1], graph_xy[from, 2],
        graph_xy[to, 1], graph_xy[to, 2],
        lwd = 0.5 + 6 * edge_scale[e],
        col = grDevices::adjustcolor("grey25", 0.20 + 0.55 * edge_scale[e])
      )
    }
  }
  node_cex = 1.5 + 2.2 * factor_usage / max(factor_usage)
  graphics::points(graph_xy, pch = 21, bg = factor_colors,
                   col = "white", lwd = 1.5, cex = node_cex)
  graphics::text(graph_xy[, 1], graph_xy[, 2],
                 labels = paste0("F", seq_len(selected_K)),
                 font = 2, cex = 0.72)
  radial = graph_xy / pmax(sqrt(rowSums(graph_xy^2)), 1e-12)
  label_xy = graph_xy + 0.16 * radial
  graphics::text(label_xy[, 1], label_xy[, 2],
                 labels = factor_node_table$top_population,
                 cex = 0.55, col = "grey25")
  graphics::title("Factor co-use graph: hubs are shared across more motifs")
}

plot_distance_comparison = function() {
  graphics::par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  graphics::plot(
    as.vector(genetic_distance), as.vector(nmf_distance), pch = 16, cex = 0.65,
    col = grDevices::adjustcolor("#0072B2", 0.5),
    xlab = "Genotype distance", ylab = "NMF loading distance",
    main = sprintf("NMF: Spearman %.2f",
                   distance_correlations$spearman_with_genotype_distance[1])
  )
  graphics::plot(
    as.vector(genetic_distance), as.vector(miso_distance), pch = 16, cex = 0.65,
    col = grDevices::adjustcolor("#D55E00", 0.5),
    xlab = "Genotype distance", ylab = "MiSo loading distance",
    main = sprintf("MiSo: Spearman %.2f",
                   distance_correlations$spearman_with_genotype_distance[2])
  )
}

phone_png("01-K-selection.png", 1400, 1000, plot_k_selection())
phone_png("02-genotype-PCA.png", 1400, 1100, plot_pca())
phone_png("03-NMF-raw-loadings.png", 1500, 2600,
          plot_loading_panels(nmf_scores, "Vanilla Poisson NMF: raw loadings"))
phone_png("04-MiSo-raw-loadings.png", 1500, 2600,
          plot_loading_panels(miso_scores, "MiSo: raw, unnormalized loadings"))
phone_png("05-population-MiSo-loadings.png", 1600, 1100,
          plot_population_loadings())
phone_png("06-population-submanifolds.png", 1600, 1100,
          plot_population_motifs())
phone_png("07-submanifold-tree.png", 1600, 1300,
          plot_submanifold_tree())
phone_png("08-factor-graph.png", 1500, 1400, plot_factor_graph())
phone_png("09-distance-preservation.png", 1600, 800,
          plot_distance_comparison())

## -------------------------------------------------------------------------
## 8. Plain-language automatic interpretation
## -------------------------------------------------------------------------

interpretation_path = file.path(OUTPUT_DIR, "interpretation.txt")
connection = file(interpretation_path, open = "wt")
sink(connection)
cat("1000 Genomes chromosome 1: Poisson MiSo interpretation\n")
cat("=======================================================\n\n")
cat("This is a working-Poisson analysis of diploid alternate-allele counts.\n")
cat("The Poisson mean factorization is the target; the conditional variance and\n")
cat("support are misspecified relative to a binomial genotype model.\n\n")
cat("Markers were physically thinned to one common SNP per 250 kb bin. This\n")
cat("reduces local dependence but is not a replacement for formal LD pruning.\n\n")
cat("Fit configuration\n")
print(data_summary, row.names = FALSE)
cat("  Selected K:", selected_K, "\n")
cat("  K selected at candidate-grid boundary:", K_AT_BOUNDARY, "\n")
cat("  Fitted S:", S, "\n")
cat("  Maximum D:", min(D, selected_K), "\n")
cat("  Best MiSo seed:", best_seed, "\n\n")

cat("Goal evaluation\n")
print(evaluation, row.names = FALSE)
cat("\nNearest populations to CDX\n")
cat("  Genotypes:", nearest_populations(genetic_distance, "CDX"), "\n")
cat("  NMF loadings:", nearest_populations(nmf_distance, "CDX"), "\n")
cat("  MiSo loadings:", nearest_populations(miso_distance, "CDX"), "\n")
cat("\nNearest populations to JPT\n")
cat("  Genotypes:", nearest_populations(genetic_distance, "JPT"), "\n")
cat("  NMF loadings:", nearest_populations(nmf_distance, "JPT"), "\n")
cat("  MiSo loadings:", nearest_populations(miso_distance, "JPT"), "\n")

cat("\nFactor interpretations\n")
for (k in seq_len(selected_K)) {
  top = order(population_given_factor[, k], decreasing = TRUE)[1:3]
  cat(sprintf(
    "  F%d: %s; breadth=%.2f; graph degree=%.3f\n",
    k, paste(pop_levels[top], collapse = ", "), factor_breadth[k],
    weighted_degree[k]
  ))
}

cat("\nActive submanifold interpretations\n")
for (s in active_motifs) {
  top_pop = order(population_given_motif[, s], decreasing = TRUE)[1:3]
  top_factor = order(motif_share[s, ], decreasing = TRUE)[
    seq_len(min(effective_dimension[s], 4))
  ]
  cat(sprintf(
    "  S%d: populations %s; D90=%d; factors %s; mixture weight=%.3f\n",
    s, paste(pop_levels[top_pop], collapse = ", "), effective_dimension[s],
    paste0("F", top_factor, collapse = ", "), fit$pi[s]
  ))
}
inactive_motifs = setdiff(seq_len(S), active_motifs)
cat("  Truncated motifs (mixture weight <", ACTIVE_MOTIF_THRESHOLD, "): ",
    if (length(inactive_motifs)) paste0("S", inactive_motifs, collapse = ", ")
    else "none", "\n", sep = "")
cat("  At Wasserstein distance <=", MERGE_DISTANCE, ", ",
    length(active_motifs), " active motifs merge to ",
    length(unique(merge_assignment)), " groups.\n", sep = "")

hub_order = order(weighted_degree, decreasing = TRUE)
broad_order = order(factor_breadth, decreasing = TRUE)
cat("\nGraph interpretation\n")
cat("  Most connected factors:", paste0("F", head(hub_order, 4),
                                         collapse = ", "), "\n")
cat("  Broadest population sharing:", paste0("F", head(broad_order, 4),
                                              collapse = ", "), "\n")
if (K_AT_BOUNDARY) {
  cat("\nCaution: selected K lies at the largest candidate. Expand the K grid.\n")
}
sink()
close(connection)

message("Done. Results are in ", normalizePath(OUTPUT_DIR))
message("Start with 05-population-MiSo-loadings.png and interpretation.txt")
