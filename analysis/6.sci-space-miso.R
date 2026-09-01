## MiSo pilot analysis of the sci-Space mouse embryo data.
##
## Run this file from the project root:
##   Rscript analysis/6.sci-space-miso.R
##
## The first run is intentionally small. It is meant to verify the full
## workflow and produce interpretable motif summaries before scaling up.

## -------------------------------------------------------------------------
## 1. Settings: these are the main lines to edit
## -------------------------------------------------------------------------

SEED = 20260830

N_CELLS = 800
N_GENES = 500
MAX_SOURCE_CELLS = 30000

K = 10  # number of shared NMF factors
S = 8   # number of MiSo motifs
D = 3   # maximum number of factors used by one motif

NMF_ITERS = 25
MF_ITERS = 12
MISO_ITERS = 8

DATA_DIR = "data/sci-space"
OUTPUT_DIR = "output/sci-space-miso"

## -------------------------------------------------------------------------
## 2. Files and packages
## -------------------------------------------------------------------------

if (!file.exists("code/miso.R")) {
  stop("Run this script from the mos-nmf-experiments project root.")
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Install the Matrix package before running this script.")
}

source("code/miso.R")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

count_file = file.path(
  DATA_DIR,
  "GSE166692_sciSpace_count_matrix.mtx.gz"
)
cell_file = file.path(
  DATA_DIR,
  "GSE166692_sciSpace_cell_metadata.tsv.gz"
)
gene_file = file.path(
  DATA_DIR,
  "GSE166692_sciSpace_gene_metadata.tsv.gz"
)

needed_files = c(count_file, cell_file, gene_file)
if (!all(file.exists(needed_files))) {
  stop(
    "Missing sci-Space files:\n",
    paste(needed_files[!file.exists(needed_files)], collapse = "\n")
  )
}

set.seed(SEED)

## -------------------------------------------------------------------------
## 3. Read metadata and choose a balanced pilot set of cells
## -------------------------------------------------------------------------

message("Reading cell and gene metadata...")

cell_meta = read.delim(
  gzfile(cell_file),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

gene_meta = read.delim(
  gzfile(gene_file),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

if (!all(c("Cell", "n.umi", "final_cluster_label") %in% names(cell_meta))) {
  stop("The cell metadata does not have the expected columns.")
}

choose_balanced_cells = function(cell_meta, n_cells, max_source_cells, seed) {
  set.seed(seed)

  source_index = seq_len(min(nrow(cell_meta), max_source_cells))
  label = cell_meta$final_cluster_label[source_index]
  umi = cell_meta$n.umi[source_index]
  upper_umi = quantile(umi, 0.995, na.rm = TRUE)

  eligible = source_index[
    !is.na(label) & nzchar(label) &
      !is.na(umi) & umi >= 500 & umi <= upper_umi
  ]

  eligible_by_label = split(
    eligible,
    cell_meta$final_cluster_label[eligible]
  )
  per_label = ceiling(n_cells / length(eligible_by_label))

  selected = unlist(lapply(eligible_by_label, function(index) {
    sample(index, min(length(index), per_label))
  }), use.names = FALSE)

  if (length(selected) > n_cells) {
    selected = sample(selected, n_cells)
  }

  if (length(selected) < n_cells) {
    remaining = setdiff(eligible, selected)
    selected = c(
      selected,
      sample(remaining, min(length(remaining), n_cells - length(selected)))
    )
  }

  sort(selected)
}

selected_cells = choose_balanced_cells(
  cell_meta = cell_meta,
  n_cells = N_CELLS,
  max_source_cells = MAX_SOURCE_CELLS,
  seed = SEED
)

pilot_meta = cell_meta[selected_cells, , drop = FALSE]
rownames(pilot_meta) = pilot_meta$Cell

message(
  "Selected ", length(selected_cells), " cells from ",
  length(unique(pilot_meta$final_cluster_label)), " annotated lineages."
)

## -------------------------------------------------------------------------
## 4. Stream only the selected cells from the large Matrix Market file
## -------------------------------------------------------------------------

## The full count matrix has 150 million nonzero entries and should not be
## loaded for a pilot analysis. This helper streams the compressed file and
## keeps only the selected columns. The Matrix Market file is sorted by cell.
read_selected_mtx_columns = function(mtx_file, selected_columns, n_rows) {
  column_map_file = tempfile(fileext = ".txt")
  on.exit(unlink(column_map_file), add = TRUE)

  column_map = data.frame(
    original = selected_columns,
    selected = seq_along(selected_columns)
  )
  write.table(
    column_map,
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

  message("Streaming selected cells from the count matrix...")
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

counts_gene_by_cell = read_selected_mtx_columns(
  mtx_file = count_file,
  selected_columns = selected_cells,
  n_rows = nrow(gene_meta)
)

colnames(counts_gene_by_cell) = pilot_meta$Cell

gene_symbol = gene_meta$gene_short_name
gene_symbol[is.na(gene_symbol) | !nzchar(gene_symbol)] =
  rownames(gene_meta)[is.na(gene_symbol) | !nzchar(gene_symbol)]
gene_symbol = make.unique(gene_symbol)
rownames(counts_gene_by_cell) = gene_symbol

Y_all = Matrix::t(counts_gene_by_cell)
rm(counts_gene_by_cell)

## -------------------------------------------------------------------------
## 5. Choose overdispersed genes and form a small dense matrix
## -------------------------------------------------------------------------

message("Selecting ", N_GENES, " overdispersed genes...")

gene_detection = Matrix::colSums(Y_all > 0)
gene_mean = Matrix::colMeans(Y_all)

Y_squared = Y_all
Y_squared@x = Y_squared@x^2
gene_variance = Matrix::colMeans(Y_squared) - gene_mean^2
rm(Y_squared)

gene_dispersion = gene_variance / pmax(gene_mean, 1e-8)
technical_gene = grepl("^(mt-|Rpl|Rps)", colnames(Y_all), ignore.case = TRUE)

eligible_gene =
  gene_detection >= max(10, round(0.01 * nrow(Y_all))) &
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
rm(Y_all)

pilot_genes = data.frame(
  gene = colnames(Y),
  mean = gene_mean[selected_genes],
  detection_rate = gene_detection[selected_genes] / nrow(Y),
  dispersion = gene_dispersion[selected_genes]
)

write.csv(
  pilot_meta,
  file.path(OUTPUT_DIR, "selected_cells.csv"),
  row.names = FALSE
)
write.csv(
  pilot_genes,
  file.path(OUTPUT_DIR, "selected_genes.csv"),
  row.names = FALSE
)

message("Pilot matrix: ", nrow(Y), " cells x ", ncol(Y), " genes.")
message("Median selected-gene count per cell: ", median(rowSums(Y)))

## -------------------------------------------------------------------------
## 6. Fit MiSo
## -------------------------------------------------------------------------

message(
  "Fitting MiSo with K=", K, ", S=", S, ", D=", D, ".",
  " This may take several minutes."
)

fit = miso(
  Y = Y,
  K = K,
  S = S,
  D = D,
  max_iters = MISO_ITERS,
  n_inner = 2,
  mf_max_iters = MF_ITERS,
  mf_nmf_iters = NMF_ITERS,
  init_seed = SEED,
  update_prior = TRUE,
  update_F = FALSE,
  update_gamma = TRUE,
  surplus_slots = "uniform",
  block_size = 100
)

saveRDS(fit, file.path(OUTPUT_DIR, "miso_fit.rds"))

## -------------------------------------------------------------------------
## 7. Summarize factors and motifs
## -------------------------------------------------------------------------

top_factor_genes = do.call(rbind, lapply(seq_len(K), function(k) {
  index = head(order(fit$F[k, ], decreasing = TRUE), 12)
  data.frame(
    factor = paste0("F", k),
    rank = seq_along(index),
    gene = colnames(Y)[index],
    weight = fit$F[k, index]
  )
}))

factor_scores = motif_factor_scores(fit)$scores
factor_scores = factor_scores / pmax(rowSums(factor_scores), 1e-12)
colnames(factor_scores) = paste0("F", seq_len(K))
rownames(factor_scores) = paste0("S", seq_len(S))

celltype = pilot_meta$final_cluster_label
celltype_levels = sort(unique(celltype))
celltype_motif = do.call(rbind, lapply(celltype_levels, function(label) {
  colMeans(fit$omega[celltype == label, , drop = FALSE])
}))
rownames(celltype_motif) = celltype_levels
colnames(celltype_motif) = paste0("S", seq_len(S))

summarize_motif = function(s) {
  ## A factor is called active when it explains at least 5% of the motif.
  ## This combines duplicated fitted dimensions that choose the same factor.
  factor_index = which(factor_scores[s, ] >= 0.05)
  if (length(factor_index) == 0) {
    factor_index = which.max(factor_scores[s, ])
  }
  factor_index = factor_index[
    order(factor_scores[s, factor_index], decreasing = TRUE)
  ]

  celltype_weight = sort(celltype_motif[, s], decreasing = TRUE)
  top_celltype = names(head(celltype_weight, 3))

  factor_labels = paste0(
    "F", factor_index, " (",
    formatC(factor_scores[s, factor_index], digits = 2, format = "f"),
    ")"
  )
  celltype_labels = paste0(
    top_celltype, " (",
    formatC(celltype_weight[top_celltype], digits = 2, format = "f"),
    ")"
  )
  representative_genes = vapply(factor_index, function(k) {
    paste(
      head(top_factor_genes$gene[top_factor_genes$factor == paste0("F", k)], 4),
      collapse = "/"
    )
  }, character(1))

  data.frame(
    motif = paste0("S", s),
    mixture_weight = fit$pi[s],
    assigned_cells = sum(fit$z_hat == s),
    mean_assignment_certainty = if (any(fit$z_hat == s)) {
      mean(apply(fit$omega[fit$z_hat == s, , drop = FALSE], 1, max))
    } else {
      NA_real_
    },
    effective_dimension = length(factor_index),
    top_cell_types = paste(celltype_labels, collapse = "; "),
    top_factors = paste(factor_labels, collapse = "; "),
    representative_genes = paste(representative_genes, collapse = "; ")
  )
}

motif_summary = do.call(rbind, lapply(seq_len(S), summarize_motif))
slot_summary = motif_slot_activity(fit)

write.csv(
  top_factor_genes,
  file.path(OUTPUT_DIR, "factor_top_genes.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(motif = rownames(factor_scores), factor_scores),
  file.path(OUTPUT_DIR, "motif_factor_scores.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(cell_type = rownames(celltype_motif), celltype_motif),
  file.path(OUTPUT_DIR, "celltype_motif_scores.csv"),
  row.names = FALSE
)
write.csv(
  motif_summary,
  file.path(OUTPUT_DIR, "motif_summary.csv"),
  row.names = FALSE
)
write.csv(
  slot_summary,
  file.path(OUTPUT_DIR, "motif_dimension_activity.csv"),
  row.names = FALSE
)

## -------------------------------------------------------------------------
## 8. Simple plots using base R
## -------------------------------------------------------------------------

motif_colors = grDevices::hcl.colors(S, palette = "Dark 3")

draw_umap_motifs = function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(
    pilot_meta$umap1,
    pilot_meta$umap2,
    col = motif_colors[fit$z_hat],
    pch = 16,
    cex = 0.6,
    xlab = "UMAP 1",
    ylab = "UMAP 2",
    main = "Dominant MiSo motif"
  )
  legend(
    "topright",
    legend = paste0("S", seq_len(S)),
    col = motif_colors,
    pch = 16,
    cex = 0.7,
    bty = "n"
  )
  plot(
    pilot_meta$umap1,
    pilot_meta$umap2,
    col = gray(1 - apply(fit$omega, 1, max)),
    pch = 16,
    cex = 0.6,
    xlab = "UMAP 1",
    ylab = "UMAP 2",
    main = "Assignment uncertainty\n(darker means more certain)"
  )
}

pdf(file.path(OUTPUT_DIR, "umap_motifs.pdf"), width = 11, height = 5)
draw_umap_motifs()
dev.off()
png(
  file.path(OUTPUT_DIR, "umap_motifs.png"),
  width = 1650, height = 750, res = 150
)
draw_umap_motifs()
dev.off()

draw_motif_factor_heatmap = function() {
  heatmap(
    factor_scores,
    Rowv = NA,
    Colv = NA,
    scale = "none",
    col = hcl.colors(40, "YlOrRd", rev = TRUE),
    margins = c(6, 6),
    xlab = "Shared NMF factor",
    ylab = "MiSo motif",
    main = "Factor composition of each motif"
  )
}

pdf(file.path(OUTPUT_DIR, "motif_factor_heatmap.pdf"), width = 8, height = 6)
draw_motif_factor_heatmap()
dev.off()
png(
  file.path(OUTPUT_DIR, "motif_factor_heatmap.png"),
  width = 1200, height = 900, res = 150
)
draw_motif_factor_heatmap()
dev.off()

draw_celltype_motif_heatmap = function() {
  heatmap(
    celltype_motif,
    Rowv = NA,
    Colv = NA,
    scale = "none",
    col = hcl.colors(40, "Blues 3", rev = TRUE),
    margins = c(6, 15),
    xlab = "MiSo motif",
    ylab = "Annotated lineage",
    main = "Average motif responsibility by lineage"
  )
}

pdf(file.path(OUTPUT_DIR, "celltype_motif_heatmap.pdf"), width = 8, height = 10)
draw_celltype_motif_heatmap()
dev.off()
png(
  file.path(OUTPUT_DIR, "celltype_motif_heatmap.png"),
  width = 1200, height = 1500, res = 150
)
draw_celltype_motif_heatmap()
dev.off()

## Plot the six slides with the most selected cells. These are pilot cells,
## not all cells on a slide, so use this plot to look for broad spatial
## coherence rather than fine tissue boundaries.
has_spatial_coordinates = all(
  c("slide_id", "coords.x1", "coords.x2") %in% names(pilot_meta)
)

if (has_spatial_coordinates) {
  slide_count = sort(table(pilot_meta$slide_id), decreasing = TRUE)
  spatial_slides = names(head(slide_count, 6))

  draw_spatial_motifs = function() {
    par(mfrow = c(2, 3), mar = c(3, 3, 3, 1))
    for (slide in spatial_slides) {
      index = which(pilot_meta$slide_id == slide)
      plot(
        pilot_meta$coords.x1[index],
        pilot_meta$coords.x2[index],
        col = motif_colors[fit$z_hat[index]],
        pch = 16,
        cex = 0.9,
        asp = 1,
        xlab = "x",
        ylab = "y",
        main = paste0(slide, " (n=", length(index), ")")
      )
      if (slide == spatial_slides[1]) {
        legend(
          "topright",
          legend = paste0("S", seq_len(S)),
          col = motif_colors,
          pch = 16,
          cex = 0.55,
          bty = "n"
        )
      }
    }
  }

  pdf(file.path(OUTPUT_DIR, "spatial_motifs.pdf"), width = 11, height = 7)
  draw_spatial_motifs()
  dev.off()
  png(
    file.path(OUTPUT_DIR, "spatial_motifs.png"),
    width = 1650, height = 1050, res = 150
  )
  draw_spatial_motifs()
  dev.off()
}

## -------------------------------------------------------------------------
## 9. Print and save a first interpretation
## -------------------------------------------------------------------------

interpretation = c(
  "MiSo sci-Space pilot interpretation",
  "===================================",
  "",
  paste0(
    "This is a pilot fit using ", nrow(Y), " cells, ", ncol(Y),
    " genes, K=", K, ", S=", S, ", and D=", D, "."
  ),
  "The factors were learned during MF-Poisson-SuSiE initialization and then",
  "held fixed while MiSo learned recurring factor combinations.",
  "An effective dimension counts factor scores of at least 0.05; duplicate",
  "overfitted dimensions assigned to the same factor are counted only once.",
  "",
  capture.output(print(motif_summary, row.names = FALSE)),
  "",
  "Interpret these motifs as hypotheses, not final biological labels.",
  "A convincing motif should have coherent top genes, concentrate in related",
  "cell types or spatial regions, and remain stable in a larger multi-seed fit."
)

writeLines(
  interpretation,
  file.path(OUTPUT_DIR, "interpretation.txt")
)

cat(paste(interpretation, collapse = "\n"), "\n")
message("Finished. Results are in ", OUTPUT_DIR, ".")
