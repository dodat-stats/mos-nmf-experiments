## Larger sci-Space experiment for evaluating the practical goals of MiSo.
##
## Run from the project root:
##   Rscript --vanilla analysis/7.sci-space-miso-large.R
##
## The script caches preprocessing and every fitted seed. If it is interrupted,
## run the same command again and it will continue from the cached results.

## -------------------------------------------------------------------------
## 1. Settings: edit these first
## -------------------------------------------------------------------------

DATA_SEED = 20260830
FIT_SEEDS = c(1, 2, 3)

N_CELLS = 2000
N_GENES = 800
TRAIN_FRACTION = 0.80

K = 18  # fallback if the held-out NMF selection file is unavailable
S = 20  # deliberately over-specified MiSo submanifolds
D = 4   # maximum dimensions per fitted submanifold

NMF_ITERS = 35
MF_ITERS = 12
MISO_ITERS = 8
MISO_INNER_ITERS = 2

ACTIVE_FACTOR_THRESHOLD = 0.05
N_NEIGHBORS = 10
N_SPATIAL_PERMUTATIONS = 25

DATA_DIR = "data/sci-space"
OUTPUT_DIR = "output/sci-space-miso-large"
CACHE_VERSION = "v2"
USE_HELDOUT_SELECTED_K = TRUE
K_SELECTION_FILE = file.path(
  OUTPUT_DIR, "K-sensitivity", "selected_K.rds"
)

## -------------------------------------------------------------------------
## 2. Files and packages
## -------------------------------------------------------------------------

if (!file.exists("code/miso-benchmark-utils.R")) {
  stop("Run this script from the mos-nmf-experiments project root.")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Install the Matrix package before running this script.")
}

source("code/miso-benchmark-utils.R")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (USE_HELDOUT_SELECTED_K && file.exists(K_SELECTION_FILE)) {
  K = readRDS(K_SELECTION_FILE)
  message("Using K=", K, " selected by held-out vanilla Poisson NMF.")
} else if (USE_HELDOUT_SELECTED_K) {
  message(
    "No held-out K selection file was found; using fallback K=", K,
    ". Run analysis/8.sci-space-K-sensitivity.R to create it."
  )
}
fit_dir = file.path(
  OUTPUT_DIR,
  paste0(
    "fits-", CACHE_VERSION,
    "-n", N_CELLS, "-m", N_GENES,
    "-K", K, "-S", S, "-D", D
  )
)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

count_file = file.path(DATA_DIR, "GSE166692_sciSpace_count_matrix.mtx.gz")
cell_file = file.path(DATA_DIR, "GSE166692_sciSpace_cell_metadata.tsv.gz")
gene_file = file.path(DATA_DIR, "GSE166692_sciSpace_gene_metadata.tsv.gz")

needed_files = c(count_file, cell_file, gene_file)
if (!all(file.exists(needed_files))) {
  stop(
    "Missing sci-Space files:\n",
    paste(needed_files[!file.exists(needed_files)], collapse = "\n")
  )
}

## -------------------------------------------------------------------------
## 3. Data helpers
## -------------------------------------------------------------------------

choose_balanced_cells = function(cell_meta, n_cells, seed) {
  set.seed(seed)

  label = cell_meta$final_cluster_label
  sample_id = cell_meta$sample
  umi = cell_meta$n.umi
  upper_umi = quantile(umi, 0.995, na.rm = TRUE)
  eligible = which(
    !is.na(label) & nzchar(label) &
      !is.na(sample_id) &
      !is.na(umi) & umi >= 500 & umi <= upper_umi
  )

  labels = sort(unique(label[eligible]))
  target = rep(floor(n_cells / length(labels)), length(labels))
  target[seq_len(n_cells - sum(target))] = target[seq_len(n_cells - sum(target))] + 1
  names(target) = labels

  selected = integer(0)
  for (current_label in labels) {
    candidates = eligible[label[eligible] == current_label]
    by_sample = split(candidates, sample_id[candidates])
    sample_target = rep(
      floor(target[current_label] / length(by_sample)),
      length(by_sample)
    )
    sample_target[seq_len(target[current_label] - sum(sample_target))] =
      sample_target[seq_len(target[current_label] - sum(sample_target))] + 1

    chosen = unlist(Map(function(index, number) {
      sample(index, min(length(index), number))
    }, by_sample, sample_target), use.names = FALSE)

    if (length(chosen) < target[current_label]) {
      remaining = setdiff(candidates, chosen)
      chosen = c(
        chosen,
        sample(
          remaining,
          min(length(remaining), target[current_label] - length(chosen))
        )
      )
    }
    selected = c(selected, chosen)
  }

  if (length(selected) < n_cells) {
    remaining = setdiff(eligible, selected)
    selected = c(selected, sample(remaining, n_cells - length(selected)))
  }
  sort(sample(selected, n_cells))
}

## The complete file contains 150 million nonzero entries. This function
## streams it and keeps only the requested cells, so the full matrix is never
## held in memory.
read_selected_mtx_columns = function(mtx_file, selected_columns, n_rows) {
  column_map_file = tempfile(fileext = ".txt")
  on.exit(unlink(column_map_file), add = TRUE)
  write.table(
    data.frame(
      original = selected_columns,
      selected = seq_along(selected_columns)
    ),
    column_map_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  awk_program = paste0(
    "NR==FNR {keep[$1]=$2; next} ",
    "FNR<=3 {next} ",
    "$2>max_col {exit} ",
    "($2 in keep) {print $1, keep[$2], $3}"
  )
  command = paste(
    "gzip -dc", shQuote(mtx_file), "|",
    "awk", paste0("-v max_col=", max(selected_columns)),
    shQuote(awk_program),
    shQuote(column_map_file), "-"
  )

  connection = pipe(command, open = "r")
  on.exit(close(connection), add = TRUE)
  entries = scan(
    connection,
    what = list(gene = integer(), cell = integer(), count = double()),
    quiet = TRUE
  )

  Matrix::sparseMatrix(
    i = entries$gene,
    j = entries$cell,
    x = entries$count,
    dims = c(n_rows, length(selected_columns)),
    giveCsparse = TRUE
  )
}

preprocess_file = file.path(
  OUTPUT_DIR,
  paste0(
    "preprocessed-n", N_CELLS,
    "-m", N_GENES,
    "-seed", DATA_SEED, ".rds"
  )
)

if (file.exists(preprocess_file)) {
  message("Loading cached preprocessed data...")
  dat = readRDS(preprocess_file)
} else {
  message("Reading metadata and selecting balanced cells from both samples...")
  cell_meta = read.delim(
    gzfile(cell_file), stringsAsFactors = FALSE, check.names = FALSE,
    quote = "", comment.char = ""
  )
  gene_meta = read.delim(
    gzfile(gene_file), stringsAsFactors = FALSE, check.names = FALSE,
    quote = "", comment.char = ""
  )

  needed_columns = c(
    "Cell", "sample", "n.umi", "final_cluster_label", "slide_id",
    "coords.x1", "coords.x2", "umap1", "umap2"
  )
  if (!all(needed_columns %in% names(cell_meta))) {
    stop("The cell metadata does not have the expected columns.")
  }

  selected_cells = choose_balanced_cells(cell_meta, N_CELLS, DATA_SEED)
  selected_meta = cell_meta[selected_cells, , drop = FALSE]
  rownames(selected_meta) = selected_meta$Cell

  message(
    "Streaming ", length(selected_cells),
    " cells from the full count matrix..."
  )
  counts = read_selected_mtx_columns(
    count_file, selected_cells, nrow(gene_meta)
  )
  colnames(counts) = selected_meta$Cell

  gene_symbol = gene_meta$gene_short_name
  missing_symbol = is.na(gene_symbol) | !nzchar(gene_symbol)
  gene_symbol[missing_symbol] = rownames(gene_meta)[missing_symbol]
  rownames(counts) = make.unique(gene_symbol)

  Y_all = Matrix::t(counts)
  rm(counts)

  message("Selecting ", N_GENES, " overdispersed genes...")
  gene_detection = Matrix::colSums(Y_all > 0)
  gene_mean = Matrix::colMeans(Y_all)
  Y_squared = Y_all
  Y_squared@x = Y_squared@x^2
  gene_variance = Matrix::colMeans(Y_squared) - gene_mean^2
  rm(Y_squared)

  gene_dispersion = gene_variance / pmax(gene_mean, 1e-8)
  technical_gene = grepl(
    "^(mt-|Rpl|Rps)", colnames(Y_all), ignore.case = TRUE
  )
  eligible_gene =
    gene_detection >= max(20, round(0.01 * nrow(Y_all))) &
    gene_detection <= round(0.80 * nrow(Y_all)) &
    gene_mean > 0 &
    !technical_gene

  gene_order = order(gene_dispersion, decreasing = TRUE, na.last = NA)
  gene_order = gene_order[eligible_gene[gene_order]]
  selected_genes = head(gene_order, N_GENES)
  if (length(selected_genes) < N_GENES) {
    warning("Only ", length(selected_genes), " genes passed the filters.")
  }

  Y = as.matrix(Y_all[, selected_genes, drop = FALSE])
  selected_gene_table = data.frame(
    gene = colnames(Y),
    mean = gene_mean[selected_genes],
    detection_rate = gene_detection[selected_genes] / nrow(Y),
    dispersion = gene_dispersion[selected_genes]
  )

  dat = list(
    Y = Y,
    meta = selected_meta,
    genes = selected_gene_table,
    original_cell_index = selected_cells
  )
  saveRDS(dat, preprocess_file)
}

write.csv(
  dat$meta, file.path(OUTPUT_DIR, "selected_cells.csv"), row.names = FALSE
)
write.csv(
  dat$genes, file.path(OUTPUT_DIR, "selected_genes.csv"), row.names = FALSE
)

message(
  "Analysis matrix: ", nrow(dat$Y), " cells x ", ncol(dat$Y),
  " genes; samples: ",
  paste(names(table(dat$meta$sample)), table(dat$meta$sample), collapse = ", ")
)

## -------------------------------------------------------------------------
## 4. Make a reproducible Poisson-thinning train/test split
## -------------------------------------------------------------------------

set.seed(DATA_SEED + 1)
Y_train = matrix(
  rbinom(length(dat$Y), size = as.vector(dat$Y), prob = TRAIN_FRACTION),
  nrow = nrow(dat$Y),
  ncol = ncol(dat$Y),
  dimnames = dimnames(dat$Y)
)
Y_test = dat$Y - Y_train

## -------------------------------------------------------------------------
## 5. Evaluation helpers
## -------------------------------------------------------------------------

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
  entropy_x = -sum(px[px > 0] * log(px[px > 0]))
  entropy_y = -sum(py[py > 0] * log(py[py > 0]))
  mutual_information / sqrt(pmax(entropy_x * entropy_y, eps))
}

partition_purity = function(cluster, label) {
  sum(apply(table(cluster, label), 1, max)) / length(label)
}

adjusted_rand_index = function(x, y) {
  tab = table(x, y)
  choose_two = function(value) value * (value - 1) / 2
  observed = sum(choose_two(tab))
  row_pairs = sum(choose_two(rowSums(tab)))
  column_pairs = sum(choose_two(colSums(tab)))
  total_pairs = choose_two(sum(tab))
  expected = row_pairs * column_pairs / pmax(total_pairs, 1)
  upper = (row_pairs + column_pairs) / 2
  if (abs(upper - expected) < 1e-12) return(1)
  (observed - expected) / (upper - expected)
}

knn_indices = function(coordinates, k, group = NULL) {
  coordinates = as.matrix(coordinates)
  n = nrow(coordinates)
  if (is.null(group)) group = rep("all", n)
  result = matrix(NA_integer_, nrow = n, ncol = k)

  for (index in split(seq_len(n), group)) {
    if (length(index) <= 1) next
    current_k = min(k, length(index) - 1)
    for (i in index) {
      distance = rowSums(
        (coordinates[index, , drop = FALSE] - coordinates[i, ])^2
      )
      distance[index == i] = Inf
      result[i, seq_len(current_k)] = index[head(order(distance), current_k)]
    }
  }
  result
}

neighbor_agreement = function(value, neighbor) {
  row_index = matrix(
    seq_len(nrow(neighbor)), nrow = nrow(neighbor), ncol = ncol(neighbor)
  )
  keep = !is.na(neighbor)
  mean(value[row_index[keep]] == value[neighbor[keep]])
}

local_label_disagreement = function(label, neighbor) {
  row_index = matrix(
    seq_len(nrow(neighbor)), nrow = nrow(neighbor), ncol = ncol(neighbor)
  )
  disagreement = matrix(NA_real_, nrow(neighbor), ncol(neighbor))
  keep = !is.na(neighbor)
  disagreement[keep] =
    label[row_index[keep]] != label[neighbor[keep]]
  rowMeans(disagreement, na.rm = TRUE)
}

permuted_neighbor_agreement = function(value, neighbor, group, n_perm, seed) {
  set.seed(seed)
  group_index = split(seq_along(value), group)
  mean(replicate(n_perm, {
    permuted = value
    for (index in group_index) permuted[index] = sample(permuted[index])
    neighbor_agreement(permuted, neighbor)
  }))
}

normalize_rows_simple = function(x, eps = 1e-12) {
  x / pmax(rowSums(x), eps)
}

cluster_gene_signatures = function(loadings, F, cluster, S) {
  loading_proportion = normalize_rows_simple(loadings)
  center = do.call(rbind, lapply(seq_len(S), function(s) {
    colMeans(loading_proportion[cluster == s, , drop = FALSE])
  }))
  normalize_rows_simple(center %*% F)
}

best_match_cosine = function(A, B, eps = 1e-12) {
  A = A / pmax(sqrt(rowSums(A^2)), eps)
  B = B / pmax(sqrt(rowSums(B^2)), eps)
  similarity = A %*% t(B)
  mean(c(apply(similarity, 1, max), apply(similarity, 2, max)))
}

umap_neighbor = knn_indices(
  dat$meta[, c("umap1", "umap2")], N_NEIGHBORS
)
spatial_neighbor = knn_indices(
  dat$meta[, c("coords.x1", "coords.x2")],
  N_NEIGHBORS,
  group = interaction(dat$meta$sample, dat$meta$slide_id, drop = TRUE)
)
lineage = dat$meta$final_cluster_label
lineage_disagreement = local_label_disagreement(lineage, umap_neighbor)
spatial_group = interaction(dat$meta$sample, dat$meta$slide_id, drop = TRUE)

heldout_metrics = function(train_mean) {
  test_mean = pmax(
    (1 - TRAIN_FRACTION) / TRAIN_FRACTION * train_mean,
    1e-12
  )
  data.frame(
    heldout_deviance_per_entry =
      poisson_deviance(Y_test, test_mean) / length(Y_test),
    heldout_mean_log_score = mean(
      dpois(Y_test, lambda = test_mean, log = TRUE)
    )
  )
}

partition_metrics = function(cluster, method, seed) {
  spatial_observed = neighbor_agreement(cluster, spatial_neighbor)
  spatial_null = permuted_neighbor_agreement(
    cluster, spatial_neighbor, spatial_group,
    N_SPATIAL_PERMUTATIONS, seed = 10000 + seed
  )
  data.frame(
    method = method,
    seed = seed,
    lineage_nmi = normalized_mutual_information(cluster, lineage),
    lineage_purity = partition_purity(cluster, lineage),
    umap_neighbor_agreement = neighbor_agreement(cluster, umap_neighbor),
    spatial_neighbor_agreement = spatial_observed,
    spatial_null_agreement = spatial_null,
    spatial_excess_agreement = spatial_observed - spatial_null
  )
}

## -------------------------------------------------------------------------
## 6. Fit NMF + k-means and MiSo for each seed
## -------------------------------------------------------------------------

run_one_seed = function(seed) {
  fit_file = file.path(fit_dir, sprintf("fit-seed-%03d.rds", seed))
  if (file.exists(fit_file)) {
    message("Loading cached seed ", seed, "...")
    return(readRDS(fit_file))
  }

  message("Seed ", seed, ": fitting Poisson NMF + k-means...")
  nmf_fit = poisson_nmf_init(
    Y_train, K = K, max_iters = NMF_ITERS, init_seed = seed
  )
  nmf_cluster = kmeans_cluster_scores(nmf_fit$L, S = S, seed = seed)
  nmf_proportion = normalize_rows_simple(nmf_fit$L)
  nmf_active = rowSums(nmf_proportion >= ACTIVE_FACTOR_THRESHOLD)
  nmf_partition = partition_metrics(
    nmf_cluster, "Poisson NMF + k-means", seed
  )
  nmf_prediction = heldout_metrics(nmf_fit$L %*% nmf_fit$F)
  nmf_metrics = cbind(
    nmf_partition,
    data.frame(
      mean_active_factors = mean(nmf_active),
      max_active_factors = max(nmf_active),
      compression_ratio = K / mean(nmf_active),
      mean_assignment_certainty = NA_real_,
      uncertainty_boundary_rho = NA_real_,
      uncertain_minus_certain_disagreement = NA_real_
    ),
    nmf_prediction
  )

  message("Seed ", seed, ": fitting MiSo...")
  miso_fit = miso(
    Y = Y_train,
    K = K,
    S = S,
    D = D,
    max_iters = MISO_ITERS,
    n_inner = MISO_INNER_ITERS,
    mf_max_iters = MF_ITERS,
    mf_nmf_iters = NMF_ITERS,
    init_seed = seed,
    update_prior = TRUE,
    update_F = TRUE,
    update_gamma = TRUE,
    surplus_slots = "uniform",
    block_size = 100
  )

  factor_scores = motif_factor_scores(miso_fit)$scores
  factor_scores = normalize_rows_simple(factor_scores)
  active_dimension = rowSums(
    factor_scores >= ACTIVE_FACTOR_THRESHOLD
  )
  certainty = apply(miso_fit$omega, 1, max)
  uncertainty = 1 - certainty
  uncertainty_rho = suppressWarnings(cor(
    uncertainty, lineage_disagreement,
    method = "spearman", use = "complete.obs"
  ))
  uncertainty_quartile = quantile(uncertainty, c(0.25, 0.75))
  disagreement_gap =
    mean(lineage_disagreement[uncertainty >= uncertainty_quartile[2]]) -
    mean(lineage_disagreement[uncertainty <= uncertainty_quartile[1]])

  miso_partition = partition_metrics(miso_fit$z_hat, "MiSo", seed)
  miso_train_mean = miso_observation_factor_scores(miso_fit) %*% miso_fit$F
  miso_prediction = heldout_metrics(miso_train_mean)
  miso_metrics = cbind(
    miso_partition,
    data.frame(
      mean_active_factors = sum(miso_fit$pi * active_dimension),
      max_active_factors = max(active_dimension),
      compression_ratio = K / sum(miso_fit$pi * active_dimension),
      mean_assignment_certainty = mean(certainty),
      uncertainty_boundary_rho = uncertainty_rho,
      uncertain_minus_certain_disagreement = disagreement_gap
    ),
    miso_prediction
  )

  result = list(
    seed = seed,
    metrics = rbind(nmf_metrics, miso_metrics),
    nmf = list(
      F = nmf_fit$F,
      loadings = nmf_fit$L,
      cluster = nmf_cluster,
      active_factors = nmf_active,
      signatures = cluster_gene_signatures(
        nmf_fit$L, nmf_fit$F, nmf_cluster, S
      )
    ),
    miso = list(
      fit = miso_fit,
      factor_scores = factor_scores,
      active_dimension = active_dimension,
      signatures = normalize_rows_simple(factor_scores %*% miso_fit$F)
    )
  )

  ## xi is a large temporary array from the initialization and is not needed
  ## for interpretation or evaluation after the fit.
  result$miso$fit$mf_fit$xi = NULL
  saveRDS(result, fit_file)
  result
}

fit_results = lapply(FIT_SEEDS, run_one_seed)
all_metrics = do.call(rbind, lapply(fit_results, `[[`, "metrics"))
rownames(all_metrics) = NULL
write.csv(
  all_metrics, file.path(OUTPUT_DIR, "metrics_by_seed.csv"), row.names = FALSE
)

## -------------------------------------------------------------------------
## 7. Multi-seed stability
## -------------------------------------------------------------------------

seed_pairs = combn(seq_along(fit_results), 2)
stability = do.call(rbind, lapply(seq_len(ncol(seed_pairs)), function(j) {
  a = seed_pairs[1, j]
  b = seed_pairs[2, j]
  rbind(
    data.frame(
      method = "Poisson NMF + k-means",
      seed_a = fit_results[[a]]$seed,
      seed_b = fit_results[[b]]$seed,
      partition_ari = adjusted_rand_index(
        fit_results[[a]]$nmf$cluster,
        fit_results[[b]]$nmf$cluster
      ),
      signature_cosine = best_match_cosine(
        fit_results[[a]]$nmf$signatures,
        fit_results[[b]]$nmf$signatures
      )
    ),
    data.frame(
      method = "MiSo",
      seed_a = fit_results[[a]]$seed,
      seed_b = fit_results[[b]]$seed,
      partition_ari = adjusted_rand_index(
        fit_results[[a]]$miso$fit$z_hat,
        fit_results[[b]]$miso$fit$z_hat
      ),
      signature_cosine = best_match_cosine(
        fit_results[[a]]$miso$signatures,
        fit_results[[b]]$miso$signatures
      )
    )
  )
}))
write.csv(
  stability, file.path(OUTPUT_DIR, "stability_by_seed_pair.csv"),
  row.names = FALSE
)

metric_names = setdiff(names(all_metrics), c("method", "seed"))
mean_or_na = function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}
sd_or_na = function(x) {
  if (sum(!is.na(x)) < 2) return(NA_real_)
  sd(x, na.rm = TRUE)
}
metric_summary = do.call(rbind, lapply(split(all_metrics, all_metrics$method),
                                       function(tab) {
  out = data.frame(method = unique(tab$method), n_seed = nrow(tab))
  for (metric in metric_names) {
    out[[paste0(metric, "_mean")]] = mean_or_na(tab[[metric]])
    out[[paste0(metric, "_sd")]] = sd_or_na(tab[[metric]])
  }
  out
}))
rownames(metric_summary) = NULL
write.csv(
  metric_summary, file.path(OUTPUT_DIR, "metric_summary.csv"), row.names = FALSE
)

stability_summary = do.call(rbind, lapply(split(stability, stability$method),
                                          function(tab) {
  data.frame(
    method = unique(tab$method),
    mean_partition_ari = mean(tab$partition_ari),
    mean_signature_cosine = mean(tab$signature_cosine)
  )
}))
rownames(stability_summary) = NULL
write.csv(
  stability_summary, file.path(OUTPUT_DIR, "stability_summary.csv"),
  row.names = FALSE
)

uncertainty_diagnostics = do.call(rbind, lapply(fit_results, function(result) {
  omega = result$miso$fit$omega
  uncertainty = 1 - apply(omega, 1, max)
  entropy = -rowSums(omega * log(pmax(omega, 1e-12)))
  motif_disagreement = local_label_disagreement(
    result$miso$fit$z_hat, umap_neighbor
  )
  data.frame(
    seed = result$seed,
    median_responsibility_entropy = median(entropy),
    uncertainty_99th_percentile = unname(quantile(uncertainty, 0.99)),
    uncertainty_log_umi_rho = cor(
      uncertainty, log1p(dat$meta$n.umi), method = "spearman"
    ),
    uncertainty_lineage_boundary_rho = cor(
      uncertainty, lineage_disagreement, method = "spearman"
    ),
    uncertainty_motif_boundary_rho = cor(
      uncertainty, motif_disagreement, method = "spearman"
    )
  )
}))
write.csv(
  uncertainty_diagnostics,
  file.path(OUTPUT_DIR, "uncertainty_diagnostics.csv"),
  row.names = FALSE
)

## -------------------------------------------------------------------------
## 8. Check the goals using explicit, deliberately modest thresholds
## -------------------------------------------------------------------------

mean_metric = function(method, metric) {
  mean(all_metrics[all_metrics$method == method, metric], na.rm = TRUE)
}
stability_metric = function(method, metric) {
  mean(stability[stability$method == method, metric], na.rm = TRUE)
}
status_from = function(supported, mixed = FALSE) {
  if (isTRUE(supported)) return("supported")
  if (isTRUE(mixed)) return("mixed")
  "not supported yet"
}

miso_active = mean_metric("MiSo", "mean_active_factors")
nmf_active = mean_metric("Poisson NMF + k-means", "mean_active_factors")
miso_deviance = mean_metric("MiSo", "heldout_deviance_per_entry")
nmf_deviance = mean_metric(
  "Poisson NMF + k-means", "heldout_deviance_per_entry"
)
miso_nmi = mean_metric("MiSo", "lineage_nmi")
nmf_nmi = mean_metric("Poisson NMF + k-means", "lineage_nmi")
miso_spatial = mean_metric("MiSo", "spatial_excess_agreement")
nmf_spatial = mean_metric(
  "Poisson NMF + k-means", "spatial_excess_agreement"
)
uncertainty_rho = mean_metric("MiSo", "uncertainty_boundary_rho")
uncertainty_gap = mean_metric(
  "MiSo", "uncertain_minus_certain_disagreement"
)
uncertainty_umi_rho = mean(uncertainty_diagnostics$uncertainty_log_umi_rho)
miso_ari = stability_metric("MiSo", "partition_ari")
miso_signature_stability = stability_metric("MiSo", "signature_cosine")
nmf_ari = stability_metric("Poisson NMF + k-means", "partition_ari")

goal_checks = rbind(
  data.frame(
    goal = "Compact local representation",
    status = status_from(
      miso_active <= D && miso_active <= 0.8 * nmf_active,
      miso_active <= D
    ),
    evidence = sprintf(
      "Mean active factors: MiSo %.2f versus NMF %.2f; fitted D=%d and K=%d.",
      miso_active, nmf_active, D, K
    )
  ),
  data.frame(
    goal = "Preserve held-out fit",
    status = status_from(
      miso_deviance <= 1.05 * nmf_deviance,
      miso_deviance <= 1.15 * nmf_deviance
    ),
    evidence = sprintf(
      "Held-out deviance per entry: MiSo %.4f versus NMF %.4f.",
      miso_deviance, nmf_deviance
    )
  ),
  data.frame(
    goal = "Recover external lineage structure",
    status = status_from(
      miso_nmi >= nmf_nmi - 0.02,
      miso_nmi >= nmf_nmi - 0.08
    ),
    evidence = sprintf(
      "NMI with held-out lineage labels: MiSo %.3f versus NMF %.3f.",
      miso_nmi, nmf_nmi
    )
  ),
  data.frame(
    goal = "Recover spatial organization",
    status = status_from(
      miso_spatial > 0.05 && miso_spatial >= nmf_spatial - 0.02,
      miso_spatial > 0
    ),
    evidence = sprintf(
      "Spatial neighbor agreement above permutation: MiSo %.3f versus NMF %.3f.",
      miso_spatial, nmf_spatial
    )
  ),
  data.frame(
    goal = "Meaningful assignment uncertainty",
    status = status_from(
      uncertainty_rho >= 0.10 && uncertainty_gap > 0.03,
      uncertainty_rho > 0 && uncertainty_gap > 0
    ),
    evidence = sprintf(
      paste0(
        "Spearman rho with local lineage disagreement %.3f; ",
        "uncertain-minus-certain gap %.3f; rho with log UMI %.3f."
      ),
      uncertainty_rho, uncertainty_gap, uncertainty_umi_rho
    )
  ),
  data.frame(
    goal = "Stable representation across seeds",
    status = status_from(
      miso_ari >= 0.50 && miso_signature_stability >= 0.80,
      miso_signature_stability >= 0.70
    ),
    evidence = sprintf(
      "MiSo pairwise ARI %.3f and signature cosine %.3f; NMF+k-means ARI %.3f.",
      miso_ari, miso_signature_stability, nmf_ari
    )
  )
)
write.csv(
  goal_checks, file.path(OUTPUT_DIR, "goal_checks.csv"), row.names = FALSE
)

## -------------------------------------------------------------------------
## 9. Interpret the best held-out MiSo fit
## -------------------------------------------------------------------------

miso_rows = all_metrics$method == "MiSo"
best_seed = all_metrics$seed[miso_rows][
  which.min(all_metrics$heldout_deviance_per_entry[miso_rows])
]
best = fit_results[[which(FIT_SEEDS == best_seed)]]
best_fit = best$miso$fit
best_factor_scores = best$miso$factor_scores

top_factor_genes = do.call(rbind, lapply(seq_len(K), function(k) {
  index = head(order(best_fit$F[k, ], decreasing = TRUE), 10)
  data.frame(
    factor = paste0("F", k),
    rank = seq_along(index),
    gene = colnames(Y_train)[index],
    weight = best_fit$F[k, index]
  )
}))

celltype_motif = do.call(rbind, lapply(sort(unique(lineage)), function(label) {
  colMeans(best_fit$omega[lineage == label, , drop = FALSE])
}))
rownames(celltype_motif) = sort(unique(lineage))
colnames(celltype_motif) = paste0("S", seq_len(S))

motif_summary = do.call(rbind, lapply(seq_len(S), function(s) {
  active = which(best_factor_scores[s, ] >= ACTIVE_FACTOR_THRESHOLD)
  if (length(active) == 0) active = which.max(best_factor_scores[s, ])
  active = active[order(best_factor_scores[s, active], decreasing = TRUE)]
  top_lineage = order(celltype_motif[, s], decreasing = TRUE)[seq_len(3)]
  signature = best$miso$signatures[s, ]
  top_gene = head(order(signature, decreasing = TRUE), 8)

  data.frame(
    motif = paste0("S", s),
    mixture_weight = best_fit$pi[s],
    assigned_cells = sum(best_fit$z_hat == s),
    effective_dimension = length(active),
    active_factors = paste0(
      "F", active, " (", sprintf("%.2f", best_factor_scores[s, active]), ")",
      collapse = "; "
    ),
    top_lineages = paste0(
      rownames(celltype_motif)[top_lineage], " (",
      sprintf("%.2f", celltype_motif[top_lineage, s]), ")",
      collapse = "; "
    ),
    top_signature_genes = paste(colnames(Y_train)[top_gene], collapse = "/")
  )
}))

write.csv(
  top_factor_genes, file.path(OUTPUT_DIR, "best_seed_factor_genes.csv"),
  row.names = FALSE
)
write.csv(
  motif_summary, file.path(OUTPUT_DIR, "best_seed_motif_summary.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(cell_type = rownames(celltype_motif), celltype_motif),
  file.path(OUTPUT_DIR, "best_seed_celltype_motif_scores.csv"),
  row.names = FALSE
)

largest_motif = motif_summary[which.max(motif_summary$mixture_weight), ]
goal_checks = rbind(
  goal_checks,
  data.frame(
    goal = "Readable recurring factor combinations",
    status = if (largest_motif$mixture_weight <= 0.30) "supported" else "mixed",
    evidence = sprintf(
      paste0(
        "Best fit has 1-%d active factors per submanifold and several coherent ",
        "gene programs; %s is a broad catch-all containing %.1f%% of cells."
      ),
      max(motif_summary$effective_dimension),
      largest_motif$motif,
      100 * largest_motif$mixture_weight
    )
  )
)
write.csv(
  goal_checks, file.path(OUTPUT_DIR, "goal_checks.csv"), row.names = FALSE
)

## -------------------------------------------------------------------------
## 10. Figures
## -------------------------------------------------------------------------

motif_colors = hcl.colors(S, "Dark 3")
nmf_cluster_colors = hcl.colors(S, "Set 3")
## Polychrome 36 is designed to keep many categorical colors distinguishable.
factor_colors = palette.colors(K, "Polychrome 36")

png(
  file.path(OUTPUT_DIR, "best_seed_umap.png"),
  width = 1800, height = 650, res = 150
)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
plot(
  dat$meta$umap1, dat$meta$umap2,
  col = nmf_cluster_colors[best$nmf$cluster], pch = 16, cex = 0.45,
  xlab = "UMAP 1", ylab = "UMAP 2", main = "NMF + k-means"
)
plot(
  dat$meta$umap1, dat$meta$umap2,
  col = motif_colors[best_fit$z_hat], pch = 16, cex = 0.45,
  xlab = "UMAP 1", ylab = "UMAP 2", main = "MiSo submanifold"
)
plot(
  dat$meta$umap1, dat$meta$umap2,
  col = gray(1 - apply(best_fit$omega, 1, max)), pch = 16, cex = 0.45,
  xlab = "UMAP 1", ylab = "UMAP 2",
  main = "MiSo certainty\n(darker is more certain)"
)
dev.off()

png(
  file.path(OUTPUT_DIR, "best_seed_atlas.png"),
  width = 1700, height = 850, res = 150
)
layout(matrix(c(1, 2), nrow = 1), widths = c(1.25, 1))
par(mar = c(4, 4, 4, 1))
cell_order = order(lineage, best$nmf$cluster)
barplot(
  t(normalize_rows_simple(best$nmf$loadings[cell_order, ])),
  space = 0, border = NA, col = factor_colors, axes = FALSE,
  names.arg = rep("", nrow(Y_train)),
  xlab = paste0(N_CELLS, " cells ordered by lineage"),
  main = paste0("NMF admixture: ", N_CELLS, " x ", K, " weights")
)
box()
par(mar = c(5, 5, 4, 3))
image(
  x = seq_len(K), y = seq_len(S),
  z = t(best_factor_scores),
  col = hcl.colors(50, "YlOrRd", rev = TRUE),
  axes = FALSE,
  xlab = "Shared factor", ylab = "MiSo submanifold",
  main = paste0("MiSo atlas: ", S, " x ", K, " weights")
)
axis(1, at = seq_len(K), labels = paste0("F", seq_len(K)), las = 2, cex.axis = 0.7)
axis(2, at = seq_len(S), labels = paste0("S", seq_len(S)), las = 2)
box()
dev.off()

## The usual all-cell admixture plot is unreadable at this scale. This version
## gives every annotated lineage its own panel and orders cells by their
## dominant NMF factor within the lineage.
lineage_levels = sort(unique(lineage))
panel_columns = 5
panel_rows = ceiling(length(lineage_levels) / panel_columns)
panel_id = seq_len(panel_rows * panel_columns)
panel_id[panel_id > length(lineage_levels)] = 0
legend_id = length(lineage_levels) + 1
layout_matrix = rbind(
  matrix(panel_id, nrow = panel_rows, byrow = TRUE),
  rep(legend_id, panel_columns)
)

png(
  file.path(OUTPUT_DIR, "nmf_admixture_by_lineage.png"),
  width = 1800, height = 1450, res = 150
)
layout(layout_matrix, heights = c(rep(1, panel_rows), 0.24))
loading_proportion = normalize_rows_simple(best$nmf$loadings)
for (j in seq_along(lineage_levels)) {
  current_lineage = lineage_levels[j]
  index = which(lineage == current_lineage)
  dominant_factor = max.col(loading_proportion[index, , drop = FALSE])
  dominant_weight = loading_proportion[cbind(index, dominant_factor)]
  index = index[order(dominant_factor, -dominant_weight)]

  par(mar = c(1, if (j %% panel_columns == 1) 3 else 1, 3, 1))
  barplot(
    t(loading_proportion[index, , drop = FALSE]),
    space = 0,
    border = NA,
    col = factor_colors,
    axes = FALSE,
    names.arg = rep("", length(index)),
    main = paste0(current_lineage, " (n=", length(index), ")"),
    cex.main = 0.70
  )
  if (j %% panel_columns == 1) {
    axis(2, at = c(0, 0.5, 1), labels = c("0", ".5", "1"), las = 1,
         cex.axis = 0.65)
  }
  box()
}
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = paste0("F", seq_len(K)),
  fill = factor_colors,
  ncol = 7,
  bty = "n",
  cex = 0.85,
  title = "Shared NMF factor"
)
dev.off()

png(
  file.path(OUTPUT_DIR, "best_seed_lineage_heatmap.png"),
  width = 1200, height = 1500, res = 150
)
heatmap(
  celltype_motif,
  Rowv = NA, Colv = NA, scale = "none",
  col = hcl.colors(50, "Blues 3", rev = TRUE),
  margins = c(6, 15),
  xlab = "MiSo submanifold", ylab = "Annotated lineage",
  main = paste0("Lineage enrichment, best seed ", best_seed)
)
dev.off()

png(
  file.path(OUTPUT_DIR, "goal_metrics.png"),
  width = 1700, height = 900, res = 150
)
par(mfrow = c(2, 3), mar = c(6, 4, 3, 1))
plot_metric = function(metric, main, lower_better = FALSE) {
  values = tapply(all_metrics[[metric]], all_metrics$method, mean, na.rm = TRUE)
  barplot(
    values,
    col = c("#2C7FB8", "gray65"),
    las = 2,
    ylab = metric,
    main = main,
    ylim = c(0, max(values) * 1.2)
  )
  if (lower_better) mtext("lower is better", side = 3, cex = 0.7)
}
plot_metric("mean_active_factors", "Local complexity", TRUE)
plot_metric("heldout_deviance_per_entry", "Held-out prediction", TRUE)
plot_metric("lineage_nmi", "Lineage agreement")
plot_metric("spatial_excess_agreement", "Spatial coherence")
barplot(
  tapply(stability$partition_ari, stability$method, mean),
  col = c("#2C7FB8", "gray65"), las = 2, ylab = "pairwise ARI",
  main = "Partition stability", ylim = c(0, 1)
)
barplot(
  tapply(stability$signature_cosine, stability$method, mean),
  col = c("#2C7FB8", "gray65"), las = 2, ylab = "best-match cosine",
  main = "Signature stability", ylim = c(0, 1)
)
dev.off()

## -------------------------------------------------------------------------
## 11. Plain-text report
## -------------------------------------------------------------------------

report = c(
  "MiSo larger sci-Space goal check",
  "=================================",
  "",
  sprintf(
    "Data: %d cells x %d genes across %d sample IDs; K=%d, S=%d, D=%d.",
    nrow(Y_train), ncol(Y_train), length(unique(dat$meta$sample)), K, S, D
  ),
  sprintf(
    "Fits: %d seeds; %.0f%% Poisson thinning for training and %.0f%% for testing.",
    length(FIT_SEEDS), 100 * TRAIN_FRACTION, 100 * (1 - TRAIN_FRACTION)
  ),
  paste0(
    "K was selected using held-out vanilla Poisson NMF before the MiSo fits; ",
    "S was then deliberately over-specified."
  ),
  "",
  "Goal checks",
  "-----------",
  capture.output(print(goal_checks, row.names = FALSE)),
  "",
  "Metric summary",
  "--------------",
  capture.output(print(metric_summary, row.names = FALSE)),
  "",
  "Stability summary",
  "-----------------",
  capture.output(print(stability_summary, row.names = FALSE)),
  "",
  "Uncertainty diagnostics",
  "-----------------------",
  capture.output(print(uncertainty_diagnostics, row.names = FALSE)),
  "",
  paste0("Best held-out MiSo seed: ", best_seed),
  capture.output(print(motif_summary, row.names = FALSE)),
  "",
  "The thresholds in goal_checks.csv are diagnostic, not hypothesis tests.",
  "Biological interpretation still requires inspecting genes, spatial maps,",
  "replication across embryos, and sensitivity to K, S, and D."
)
writeLines(report, file.path(OUTPUT_DIR, "experiment_report.txt"))

cat(paste(report, collapse = "\n"), "\n")
message("Finished. Results are in ", OUTPUT_DIR, ".")
