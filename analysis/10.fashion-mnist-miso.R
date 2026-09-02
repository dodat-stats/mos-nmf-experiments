#!/usr/bin/env Rscript

## Fashion-MNIST experiment for MiSo.
##
## Run from the project root:
##   Rscript --vanilla analysis/10.fashion-mnist-miso.R
##
## For a fast end-to-end check:
##   Rscript --vanilla analysis/10.fashion-mnist-miso.R --quick
##
## The labels are never used to fit NMF or MiSo. They are used only afterward
## to interpret and evaluate the learned factors and submanifolds.

options(stringsAsFactors = FALSE)

## -------------------------------------------------------------------------
## 1. Settings
## -------------------------------------------------------------------------

QUICK = "--quick" %in% commandArgs(trailingOnly = TRUE)

DATA_SEED = 20260831
N_PER_CLASS = if (QUICK) 25 else 150
PIXEL_DIVISOR = 16
TRAIN_FRACTION = 0.80

K_CANDIDATES = if (QUICK) c(6, 10) else c(8, 12, 16, 20)
K_SELECTION_SEEDS = if (QUICK) 1 else c(1, 2)
K_NMF_ITERS = if (QUICK) 8 else 50

S = if (QUICK) 8 else 12
D = if (QUICK) 3 else 4
FIT_SEEDS = if (QUICK) 1 else c(1, 2)
MF_ITERS = if (QUICK) 3 else 25
MF_NMF_ITERS = if (QUICK) 8 else 50
MISO_ITERS = if (QUICK) 2 else 10
MISO_INNER_ITERS = if (QUICK) 1 else 2
MISO_INITIALIZATION = "distinct"
BLOCK_SIZE = 98

DATA_DIR = "data/fashion-mnist"
OUTPUT_DIR = "output/fashion-mnist-miso"
CACHE_VERSION = "beta-v1"

CLASS_NAMES = c(
  "T-shirt/top", "Trouser", "Pullover", "Dress", "Coat",
  "Sandal", "Shirt", "Sneaker", "Bag", "Ankle boot"
)

## -------------------------------------------------------------------------
## 2. Files and dependencies
## -------------------------------------------------------------------------

if (!file.exists("code/miso-benchmark-utils.R")) {
  stop("Run this script from the project root.")
}
source("code/miso-benchmark-utils.R")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

fashion_files = data.frame(
  filename = c(
    "train-images-idx3-ubyte.gz",
    "train-labels-idx1-ubyte.gz"
  ),
  url = c(
    paste0(
      "https://raw.githubusercontent.com/zalandoresearch/fashion-mnist/",
      "master/data/fashion/train-images-idx3-ubyte.gz"
    ),
    paste0(
      "https://raw.githubusercontent.com/zalandoresearch/fashion-mnist/",
      "master/data/fashion/train-labels-idx1-ubyte.gz"
    )
  ),
  md5 = c(
    "8d4fb7e6c68d591d4c3dfef9ec88bf0d",
    "25c81989df183df01b3e8a0aad5dffbe"
  )
)

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

data_paths = mapply(
  download_if_needed,
  fashion_files$filename,
  fashion_files$url,
  fashion_files$md5,
  USE.NAMES = FALSE
)

## -------------------------------------------------------------------------
## 3. Read IDX data and make a balanced analysis subset
## -------------------------------------------------------------------------

read_idx_images = function(path) {
  connection = gzfile(path, "rb")
  on.exit(close(connection))
  magic = readBin(connection, integer(), n = 1, size = 4, endian = "big")
  n_image = readBin(connection, integer(), n = 1, size = 4, endian = "big")
  n_row = readBin(connection, integer(), n = 1, size = 4, endian = "big")
  n_col = readBin(connection, integer(), n = 1, size = 4, endian = "big")
  if (magic != 2051L || n_row != 28L || n_col != 28L) {
    stop("Unexpected Fashion-MNIST image header")
  }
  pixels = readBin(
    connection, integer(), n = n_image * n_row * n_col,
    size = 1, signed = FALSE
  )
  matrix(pixels, nrow = n_image, ncol = n_row * n_col, byrow = TRUE)
}

read_idx_labels = function(path) {
  connection = gzfile(path, "rb")
  on.exit(close(connection))
  magic = readBin(connection, integer(), n = 1, size = 4, endian = "big")
  n_label = readBin(connection, integer(), n = 1, size = 4, endian = "big")
  if (magic != 2049L) stop("Unexpected Fashion-MNIST label header")
  readBin(connection, integer(), n = n_label, size = 1, signed = FALSE)
}

images = read_idx_images(data_paths[1])
labels = read_idx_labels(data_paths[2])
if (nrow(images) != length(labels)) stop("Image and label counts differ")

set.seed(DATA_SEED)
chosen = unlist(lapply(0:9, function(label) {
  candidates = which(labels == label)
  sample(candidates, N_PER_CLASS)
}))
chosen = sample(chosen)

## Pixel intensities are scaled to small integer pseudo-counts. We deliberately
## do not log-transform them: the Poisson objective operates on the original
## nonnegative scale, and the y! term is constant with respect to the fit.
Y = round(images[chosen, , drop = FALSE] / PIXEL_DIVISOR)
label = labels[chosen]
label_name = factor(CLASS_NAMES[label + 1], levels = CLASS_NAMES)
rm(images, labels)

data.frame(
  n_images = nrow(Y),
  n_pixels = ncol(Y),
  minimum = min(Y),
  maximum = max(Y),
  mean_total_intensity = mean(rowSums(Y))
) |>
  utils::write.csv(
    file.path(OUTPUT_DIR, "data-summary.csv"), row.names = FALSE
  )

## -------------------------------------------------------------------------
## 4. Select K by held-out Poisson likelihood under vanilla NMF
## -------------------------------------------------------------------------

k_cache_dir = file.path(OUTPUT_DIR, paste0("K-selection-", CACHE_VERSION))
dir.create(k_cache_dir, recursive = TRUE, showWarnings = FALSE)
k_rows = list()
k_index = 0L

for (K in K_CANDIDATES) {
  for (seed in K_SELECTION_SEEDS) {
    cache_file = file.path(k_cache_dir, paste0("K", K, "-seed", seed, ".rds"))
    if (file.exists(cache_file)) {
      result = readRDS(cache_file)
    } else {
      set.seed(1000L + seed)
      Y_train = matrix(
        rbinom(length(Y), size = as.integer(Y), prob = TRAIN_FRACTION),
        nrow = nrow(Y), ncol = ncol(Y)
      )
      Y_test = Y - Y_train
      nmf_fit = poisson_nmf_init(
        Y_train, K = K, max_iters = K_NMF_ITERS, init_seed = seed
      )
      test_mean = (1 - TRAIN_FRACTION) / TRAIN_FRACTION *
        (nmf_fit$L %*% nmf_fit$F)
      result = list(
        K = K,
        seed = seed,
        heldout_mean_log_likelihood = mean(
          dpois(Y_test, lambda = pmax(test_mean, 1e-12), log = TRUE)
        ),
        heldout_deviance_per_pixel =
          poisson_deviance(Y_test, test_mean) / length(Y_test)
      )
      saveRDS(result, cache_file)
    }
    k_index = k_index + 1L
    k_rows[[k_index]] = data.frame(
      K = result$K,
      seed = result$seed,
      heldout_mean_log_likelihood = result$heldout_mean_log_likelihood,
      heldout_deviance_per_pixel = result$heldout_deviance_per_pixel
    )
  }
}

k_results = do.call(rbind, k_rows)
k_summary = do.call(rbind, lapply(split(k_results, k_results$K), function(x) {
  data.frame(
    K = unique(x$K),
    mean_log_likelihood = mean(x$heldout_mean_log_likelihood),
    se_log_likelihood = if (nrow(x) > 1) {
      sd(x$heldout_mean_log_likelihood) / sqrt(nrow(x))
    } else {
      0
    },
    mean_deviance = mean(x$heldout_deviance_per_pixel)
  )
}))
rownames(k_summary) = NULL
selected_K = k_summary$K[which.max(k_summary$mean_log_likelihood)]
K_AT_BOUNDARY = selected_K == max(K_CANDIDATES)

utils::write.csv(k_results, file.path(OUTPUT_DIR, "K-selection-runs.csv"),
                 row.names = FALSE)
utils::write.csv(k_summary, file.path(OUTPUT_DIR, "K-selection-summary.csv"),
                 row.names = FALSE)
saveRDS(selected_K, file.path(OUTPUT_DIR, "selected-K.rds"))
message("Selected K = ", selected_K)
if (K_AT_BOUNDARY) {
  warning(
    "Held-out likelihood selected the largest candidate K. ",
    "Expand K_CANDIDATES before treating K as final."
  )
}

## -------------------------------------------------------------------------
## 5. Fit the vanilla NMF baseline and MiSo
## -------------------------------------------------------------------------

nmf_file = file.path(
  OUTPUT_DIR,
  paste0("nmf-full-", CACHE_VERSION, "-n", nrow(Y), "-K", selected_K, ".rds")
)
if (file.exists(nmf_file)) {
  nmf_full = readRDS(nmf_file)
} else {
  nmf_full = poisson_nmf_init(
    Y, K = selected_K, max_iters = max(K_NMF_ITERS, MF_NMF_ITERS),
    init_seed = DATA_SEED
  )
  saveRDS(nmf_full, nmf_file)
}

fit_dir = file.path(
  OUTPUT_DIR,
  paste0(
    "fits-", CACHE_VERSION, "-n", nrow(Y), "-K", selected_K,
    "-S", S, "-D", D, "-init-", MISO_INITIALIZATION
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
      Y = Y,
      K = selected_K,
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

final_objective = vapply(fits, function(fit) {
  finite = fit$elbo[is.finite(fit$elbo)]
  if (length(finite)) tail(finite, 1) else -Inf
}, numeric(1))
best_index = which.max(final_objective)
fit = fits[[best_index]]
best_seed = FIT_SEEDS[best_index]
saveRDS(fit, file.path(OUTPUT_DIR, "miso-fit-best.rds"))

## -------------------------------------------------------------------------
## 6. Posterior summaries and evaluation
## -------------------------------------------------------------------------

normalize_probability_rows = function(x, eps = 1e-12) {
  x = pmax(x, 0)
  x / pmax(rowSums(x), eps)
}

observation_factor_scores = miso_observation_factor_scores(fit)
observation_factor_share = normalize_probability_rows(observation_factor_scores)
nmf_factor_share = normalize_probability_rows(nmf_full$L)
submanifold_factor_scores = motif_factor_scores(fit)$scores
submanifold_factor_share = normalize_probability_rows(submanifold_factor_scores)

class_factor_mean = do.call(rbind, lapply(0:9, function(g) {
  colMeans(observation_factor_share[label == g, , drop = FALSE])
}))
rownames(class_factor_mean) = CLASS_NAMES
colnames(class_factor_mean) = paste0("F", seq_len(selected_K))

class_given_factor = sweep(class_factor_mean, 2, colSums(class_factor_mean), "/")
class_given_factor[!is.finite(class_given_factor)] = 0

class_submanifold_mean = do.call(rbind, lapply(0:9, function(g) {
  colMeans(fit$omega[label == g, , drop = FALSE])
}))
rownames(class_submanifold_mean) = CLASS_NAMES
colnames(class_submanifold_mean) = paste0("S", seq_len(S))

class_given_submanifold = matrix(0, nrow = 10, ncol = S)
for (g in 0:9) {
  class_given_submanifold[g + 1, ] = colSums(
    fit$omega[label == g, , drop = FALSE]
  )
}
class_given_submanifold = sweep(
  class_given_submanifold, 2, colSums(class_given_submanifold), "/"
)
class_given_submanifold[!is.finite(class_given_submanifold)] = 0
rownames(class_given_submanifold) = CLASS_NAMES
colnames(class_given_submanifold) = paste0("S", seq_len(S))

top_class_for_factor = apply(class_given_factor, 2, which.max)
top_class_for_submanifold = apply(class_given_submanifold, 2, which.max)

effective_dimension = apply(submanifold_factor_share, 1, function(x) {
  order_x = order(x, decreasing = TRUE)
  which(cumsum(x[order_x]) >= 0.90)[1]
})

nearest_centroid_accuracy = function(scores, labels, seed = 1) {
  set.seed(seed)
  train = unlist(lapply(0:9, function(g) {
    index = which(labels == g)
    sample(index, floor(0.7 * length(index)))
  }))
  test = setdiff(seq_along(labels), train)
  centroids = do.call(rbind, lapply(0:9, function(g) {
    colMeans(scores[train[labels[train] == g], , drop = FALSE])
  }))
  score_norm = rowSums(scores^2)
  center_norm = rowSums(centroids^2)
  distance = outer(score_norm[test], center_norm, "+") -
    2 * scores[test, , drop = FALSE] %*% t(centroids)
  prediction = max.col(-distance) - 1L
  mean(prediction == labels[test])
}

nmf_mean = nmf_full$L %*% nmf_full$F
miso_mean = observation_factor_scores %*% fit$F

evaluation = data.frame(
  representation = c("Vanilla Poisson NMF", "MiSo"),
  in_sample_deviance_per_pixel = c(
    poisson_deviance(Y, nmf_mean) / length(Y),
    poisson_deviance(Y, miso_mean) / length(Y)
  ),
  nearest_centroid_accuracy = c(
    nearest_centroid_accuracy(nmf_factor_share, label, DATA_SEED),
    nearest_centroid_accuracy(observation_factor_share, label, DATA_SEED)
  ),
  average_active_factors_per_image_90pct = c(
    mean(apply(nmf_factor_share, 1, function(x) {
      which(cumsum(sort(x, decreasing = TRUE)) >= 0.90)[1]
    })),
    mean(apply(observation_factor_share, 1, function(x) {
      which(cumsum(sort(x, decreasing = TRUE)) >= 0.90)[1]
    }))
  ),
  average_submanifold_dimension_90pct = c(NA_real_, mean(effective_dimension))
)

utils::write.csv(evaluation, file.path(OUTPUT_DIR, "goal-evaluation.csv"),
                 row.names = FALSE)
utils::write.csv(class_factor_mean,
                 file.path(OUTPUT_DIR, "class-factor-mean.csv"))
utils::write.csv(class_given_factor,
                 file.path(OUTPUT_DIR, "class-given-factor.csv"))
utils::write.csv(class_submanifold_mean,
                 file.path(OUTPUT_DIR, "class-submanifold-mean.csv"))
utils::write.csv(class_given_submanifold,
                 file.path(OUTPUT_DIR, "class-given-submanifold.csv"))
utils::write.csv(submanifold_factor_share,
                 file.path(OUTPUT_DIR, "submanifold-factor-share.csv"))

## -------------------------------------------------------------------------
## 7. Wasserstein distance between submanifolds
## -------------------------------------------------------------------------

row_cosine_similarity = function(x) {
  norms = sqrt(rowSums(x^2))
  similarity = (x %*% t(x)) / pmax(norms %o% norms, 1e-12)
  pmin(pmax(similarity, 0), 1)
}

factor_cost = 1 - row_cosine_similarity(fit$F)

sinkhorn_distance = function(p, q, cost, epsilon = 0.05,
                              max_iters = 300, tol = 1e-8) {
  p = pmax(p, 1e-12); p = p / sum(p)
  q = pmax(q, 1e-12); q = q / sum(q)
  kernel = exp(-cost / epsilon)
  kernel = pmax(kernel, 1e-300)
  u = rep(1, length(p))
  v = rep(1, length(q))
  for (iter in seq_len(max_iters)) {
    old_u = u
    u = p / pmax(as.vector(kernel %*% v), 1e-300)
    v = q / pmax(as.vector(t(kernel) %*% u), 1e-300)
    if (max(abs(u - old_u) / pmax(abs(old_u), 1)) < tol) break
  }
  plan = (u %o% v) * kernel
  sum(plan * cost)
}

submanifold_distance = matrix(0, S, S)
for (s1 in seq_len(S)) {
  if (s1 == S) next
  for (s2 in (s1 + 1):S) {
    forward = sinkhorn_distance(
      submanifold_factor_share[s1, ], submanifold_factor_share[s2, ],
      factor_cost
    )
    backward = sinkhorn_distance(
      submanifold_factor_share[s2, ], submanifold_factor_share[s1, ],
      factor_cost
    )
    submanifold_distance[s1, s2] = submanifold_distance[s2, s1] =
      (forward + backward) / 2
  }
}
rownames(submanifold_distance) = colnames(submanifold_distance) =
  paste0("S", seq_len(S))
submanifold_tree = hclust(as.dist(submanifold_distance), method = "average")
utils::write.csv(submanifold_distance,
                 file.path(OUTPUT_DIR, "submanifold-wasserstein-distance.csv"))

## -------------------------------------------------------------------------
## 8. Factor co-use graph
## -------------------------------------------------------------------------

co_use = matrix(0, selected_K, selected_K)
for (s in seq_len(S)) {
  co_use = co_use + fit$pi[s] *
    (submanifold_factor_share[s, ] %o% submanifold_factor_share[s, ])
}
diag(co_use) = 0
factor_usage = colSums(
  sweep(submanifold_factor_share, 1, fit$pi, "*")
)
weighted_degree = rowSums(co_use)

upper_index = which(upper.tri(co_use) & co_use > 0, arr.ind = TRUE)
edge_table = data.frame(
  from = upper_index[, 1],
  to = upper_index[, 2],
  weight = co_use[upper_index]
)
edge_table = edge_table[order(edge_table$weight, decreasing = TRUE), ]
edge_table = head(edge_table, min(nrow(edge_table), 2 * selected_K))

factor_node_table = data.frame(
  factor = seq_len(selected_K),
  top_class = CLASS_NAMES[top_class_for_factor],
  usage = factor_usage,
  weighted_degree = weighted_degree,
  role = ifelse(
    weighted_degree >= quantile(weighted_degree, 0.75), "hub",
    ifelse(weighted_degree <= quantile(weighted_degree, 0.25), "leaf", "middle")
  )
)
utils::write.csv(factor_node_table,
                 file.path(OUTPUT_DIR, "factor-graph-nodes.csv"),
                 row.names = FALSE)
utils::write.csv(edge_table, file.path(OUTPUT_DIR, "factor-graph-edges.csv"),
                 row.names = FALSE)

graph_similarity = co_use / pmax(sqrt(rowSums(co_use^2) %o% rowSums(co_use^2)),
                                  1e-12)
diag(graph_similarity) = 1
graph_distance = as.dist(pmax(1 - graph_similarity, 0))

## Use graph clustering to order a radial layout. Radius is inversely related
## to weighted degree, so broadly reused factors move inward and leaf-like
## factors remain on the perimeter. This is stable even when redundant or
## unused factors have exactly the same graph profile.
graph_order = hclust(graph_distance, method = "average")$order
graph_angle = seq(pi / 2, pi / 2 + 2 * pi, length.out = selected_K + 1)[
  -1
]
degree_scaled = (weighted_degree - min(weighted_degree)) /
  pmax(max(weighted_degree) - min(weighted_degree), 1e-12)
graph_radius = 1 - 0.45 * degree_scaled
graph_xy = matrix(0, selected_K, 2)
graph_xy[graph_order, ] = cbind(
  graph_radius[graph_order] * cos(graph_angle),
  graph_radius[graph_order] * sin(graph_angle)
)

## -------------------------------------------------------------------------
## 9. Plotting helpers
## -------------------------------------------------------------------------

phone_png = function(filename, width = 1400, height = 1600, expr) {
  grDevices::png(
    file.path(OUTPUT_DIR, filename), width = width, height = height,
    res = 180, bg = "white"
  )
  on.exit(grDevices::dev.off())
  force(expr)
}

draw_image = function(values, title = "", scale_each = TRUE) {
  image_matrix = matrix(values, nrow = 28, ncol = 28, byrow = TRUE)
  if (scale_each) image_matrix = image_matrix / pmax(max(image_matrix), 1e-12)
  graphics::image(
    t(image_matrix[28:1, , drop = FALSE]),
    col = gray.colors(256, start = 0, end = 1), axes = FALSE,
    useRaster = TRUE, main = title
  )
  graphics::box(col = "grey75")
}

draw_heatmap = function(x, main, x_labels, y_labels,
                         xlab = "", ylab = "", show_values = FALSE) {
  nr = nrow(x); nc = ncol(x)
  graphics::image(
    seq_len(nc), seq_len(nr), t(x[nr:1, , drop = FALSE]),
    col = hcl.colors(100, "YlOrRd"), axes = FALSE,
    xlab = xlab, ylab = ylab, main = main, useRaster = TRUE
  )
  graphics::axis(1, at = seq_len(nc), labels = x_labels, las = 2, cex.axis = 0.75)
  graphics::axis(2, at = seq_len(nr), labels = rev(y_labels), las = 2,
                 cex.axis = 0.75)
  graphics::box()
  if (show_values) {
    for (r in seq_len(nr)) {
      for (cc in seq_len(nc)) {
        graphics::text(cc, nr - r + 1, sprintf("%.2f", x[r, cc]), cex = 0.48)
      }
    }
  }
}

plot_k_selection = function() {
  graphics::par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))
  ylim = range(
    k_summary$mean_log_likelihood - 2 * k_summary$se_log_likelihood,
    k_summary$mean_log_likelihood + 2 * k_summary$se_log_likelihood
  )
  graphics::plot(
    k_summary$K, k_summary$mean_log_likelihood, type = "b", pch = 19,
    lwd = 2, xlab = "Number of NMF factors K",
    ylab = "Held-out mean Poisson log likelihood", ylim = ylim,
    main = paste0("Held-out selection chooses K = ", selected_K)
  )
  has_error_bar = k_summary$se_log_likelihood > 0
  if (any(has_error_bar)) {
    graphics::arrows(
      k_summary$K[has_error_bar],
      k_summary$mean_log_likelihood[has_error_bar] -
        2 * k_summary$se_log_likelihood[has_error_bar],
      k_summary$K[has_error_bar],
      k_summary$mean_log_likelihood[has_error_bar] +
        2 * k_summary$se_log_likelihood[has_error_bar],
      angle = 90, code = 3, length = 0.05
    )
  }
  graphics::abline(v = selected_K, lty = 2, col = "firebrick")
}

plot_factor_atlas = function() {
  n_col = 4
  n_row = ceiling(selected_K / n_col)
  graphics::par(mfrow = c(n_row, n_col), mar = c(0.5, 0.5, 2.3, 0.5))
  for (k in seq_len(selected_K)) {
    draw_image(
      fit$F[k, ],
      paste0("F", k, ": ", CLASS_NAMES[top_class_for_factor[k]])
    )
  }
  for (unused in seq_len(n_col * n_row - selected_K)) graphics::plot.new()
}

plot_submanifold_atlas = function() {
  n_col = 3
  n_row = ceiling(S / n_col)
  graphics::par(mfrow = c(n_row, n_col), mar = c(0.5, 0.5, 2.4, 0.5))
  for (s in seq_len(S)) {
    motif_image = as.vector(submanifold_factor_scores[s, ] %*% fit$F)
    draw_image(
      motif_image,
      paste0(
        "S", s, " (D90=", effective_dimension[s], "): ",
        CLASS_NAMES[top_class_for_submanifold[s]]
      )
    )
  }
  for (unused in seq_len(n_col * n_row - S)) graphics::plot.new()
}

plot_submanifold_usage = function() {
  graphics::layout(matrix(c(1, 2), nrow = 2), heights = c(1, 1.6))
  graphics::par(mar = c(7, 11, 3, 1), mgp = c(7, 1, 0))
  draw_heatmap(
    class_submanifold_mean,
    "Average submanifold responsibility within each clothing class",
    paste0("S", seq_len(S)), CLASS_NAMES,
    xlab = "Submanifold", ylab = "Fashion-MNIST class", show_values = TRUE
  )

  shown = unlist(lapply(0:9, function(g) head(which(label == g), 50)))
  shown = shown[order(label[shown])]
  omega_shown = fit$omega[shown, , drop = FALSE]
  graphics::par(mar = c(7, 8, 3, 1), mgp = c(5, 1, 0))
  graphics::image(
    seq_len(S), seq_len(nrow(omega_shown)), t(omega_shown[nrow(omega_shown):1, ]),
    col = hcl.colors(100, "Viridis"), axes = FALSE,
    xlab = "Submanifold", ylab = "Images grouped by class",
    main = "Individual responsibilities (up to 50 images per class)",
    useRaster = TRUE
  )
  graphics::axis(1, at = seq_len(S), labels = paste0("S", seq_len(S)),
                 las = 2, cex.axis = 0.75)
  class_mid = vapply(0:9, function(g) mean(which(label[shown] == g)), numeric(1))
  graphics::axis(2, at = nrow(omega_shown) - class_mid + 1,
                 labels = CLASS_NAMES, las = 2, cex.axis = 0.65)
  graphics::box()
}

representative_index = vapply(0:9, function(g) {
  index = which(label == g)
  class_mean = colMeans(observation_factor_share[index, , drop = FALSE])
  distance = rowSums(
    (observation_factor_share[index, , drop = FALSE] -
       matrix(class_mean, nrow = length(index), ncol = selected_K,
              byrow = TRUE))^2
  )
  index[which.min(distance)]
}, integer(1))

plot_reconstructions = function(classes) {
  graphics::par(mfrow = c(length(classes), 3), mar = c(0.2, 0.2, 1.8, 0.2),
                oma = c(0, 0, 4, 0))
  for (g in classes) {
    i = representative_index[g + 1]
    draw_image(Y[i, ], paste0(CLASS_NAMES[g + 1], ": observed"))
    draw_image(nmf_mean[i, ], "Vanilla NMF")
    draw_image(miso_mean[i, ], "MiSo")
  }
  graphics::mtext("Representative images and reconstructions", outer = TRUE,
                  line = 1, cex = 1.2)
}

plot_submanifold_tree = function() {
  old_labels = submanifold_tree$labels
  submanifold_tree$labels = paste0(
    "S", seq_len(S), "  ", CLASS_NAMES[top_class_for_submanifold]
  )
  graphics::par(mfrow = c(1, 1), mar = c(5, 3, 3, 13))
  graphics::plot(
    as.dendrogram(submanifold_tree), horiz = TRUE,
    main = "Submanifold tree from Sinkhorn/Wasserstein distance",
    xlab = "Transport distance between factor mixing measures"
  )
  submanifold_tree$labels = old_labels
}

plot_factor_graph = function() {
  graphics::par(mfrow = c(1, 1), mar = c(1, 1, 3, 1))
  x = graph_xy[, 1]; y = graph_xy[, 2]
  x_pad = 0.12 * diff(range(x)); y_pad = 0.12 * diff(range(y))
  graphics::plot(
    x, y, type = "n", axes = FALSE, xlab = "", ylab = "",
    xlim = range(x) + c(-x_pad, x_pad),
    ylim = range(y) + c(-y_pad, y_pad),
    main = "Factor graph: edges show co-use within fitted submanifolds"
  )
  edge_width = 1 + 7 * edge_table$weight / pmax(max(edge_table$weight), 1e-12)
  for (e in seq_len(nrow(edge_table))) {
    graphics::segments(
      x[edge_table$from[e]], y[edge_table$from[e]],
      x[edge_table$to[e]], y[edge_table$to[e]],
      lwd = edge_width[e], col = grDevices::adjustcolor("grey45", 0.5)
    )
  }
  class_colors = hcl.colors(10, "Dark 3")
  node_size = 1.8 + 3.2 * sqrt(factor_usage / max(factor_usage))
  graphics::points(
    x, y, pch = 21, cex = node_size,
    bg = class_colors[top_class_for_factor], col = "white", lwd = 1.5
  )
  graphics::text(x, y, labels = paste0("F", seq_len(selected_K)),
                 cex = if (selected_K > 15) 0.58 else 0.68, font = 2)
  graphics::text(
    1.10 * x, 1.10 * y,
    labels = CLASS_NAMES[top_class_for_factor],
    cex = if (selected_K > 15) 0.46 else 0.56
  )
}

## -------------------------------------------------------------------------
## 10. Write phone-readable PNG figures and a multi-page PDF
## -------------------------------------------------------------------------

phone_png("01-K-selection.png", 1400, 1100, plot_k_selection())
phone_png("02-factor-atlas.png", 1400, 400 * ceiling(selected_K / 4),
          plot_factor_atlas())
phone_png("03-submanifold-atlas.png", 1400, 440 * ceiling(S / 3),
          plot_submanifold_atlas())
phone_png("04-submanifold-usage.png", 1500, 1900, plot_submanifold_usage())
phone_png("05a-reconstructions-classes-0-4.png", 1200, 1900,
          plot_reconstructions(0:4))
phone_png("05b-reconstructions-classes-5-9.png", 1200, 1900,
          plot_reconstructions(5:9))
phone_png("06-submanifold-tree.png", 1500, 1300, plot_submanifold_tree())
phone_png("07-factor-graph.png", 1500, 1400, plot_factor_graph())

grDevices::pdf(
  file.path(OUTPUT_DIR, "fashion-mnist-miso-report.pdf"),
  width = 9, height = 10, onefile = TRUE
)
plot_k_selection()
plot_factor_atlas()
plot_submanifold_atlas()
plot_submanifold_usage()
plot_reconstructions(0:4)
plot_reconstructions(5:9)
plot_submanifold_tree()
plot_factor_graph()
grDevices::dev.off()

## -------------------------------------------------------------------------
## 11. Plain-language automatic interpretation
## -------------------------------------------------------------------------

interpretation_file = file.path(OUTPUT_DIR, "interpretation.txt")
connection = file(interpretation_file, open = "wt")
sink(connection)
cat("Fashion-MNIST MiSo interpretation\n")
cat("=================================\n\n")
cat("Fit configuration\n")
cat("  Images:", nrow(Y), "(", N_PER_CLASS, "per class)\n")
cat("  Selected K:", selected_K, "\n")
cat("  K selected at candidate-grid boundary:", K_AT_BOUNDARY, "\n")
cat("  Fitted S:", S, "\n")
cat("  Maximum D:", D, "\n")
cat("  Best seed:", best_seed, "\n\n")

cat("Goal evaluation\n")
print(evaluation, row.names = FALSE)
if (K_AT_BOUNDARY) {
  cat("\n  Caution: held-out likelihood was still increasing at the largest K.\n")
  cat("  Expand K_CANDIDATES before treating the selected K as final.\n")
}
cat("\nFactor interpretations (post-hoc class associations)\n")
for (k in seq_len(selected_K)) {
  top = order(class_given_factor[, k], decreasing = TRUE)[1:2]
  cat(sprintf(
    "  F%d: %s (%.2f), %s (%.2f); graph role=%s\n",
    k, CLASS_NAMES[top[1]], class_given_factor[top[1], k],
    CLASS_NAMES[top[2]], class_given_factor[top[2], k],
    factor_node_table$role[k]
  ))
}

cat("\nSubmanifold interpretations\n")
for (s in seq_len(S)) {
  top_class = order(class_given_submanifold[, s], decreasing = TRUE)[1:2]
  top_factor = order(submanifold_factor_share[s, ], decreasing = TRUE)[
    seq_len(min(effective_dimension[s], 4))
  ]
  cat(sprintf(
    "  S%d: classes %s (%.2f), %s (%.2f); D90=%d; factors %s\n",
    s,
    CLASS_NAMES[top_class[1]], class_given_submanifold[top_class[1], s],
    CLASS_NAMES[top_class[2]], class_given_submanifold[top_class[2], s],
    effective_dimension[s], paste0("F", top_factor, collapse = ", ")
  ))
}

hub_order = order(weighted_degree, decreasing = TRUE)
leaf_order = order(weighted_degree)
cat("\nGraph interpretation\n")
cat("  Most connected factors:",
    paste0("F", head(hub_order, 4), collapse = ", "), "\n")
cat("  Most leaf-like factors:",
    paste0("F", head(leaf_order, 4), collapse = ", "), "\n")
cat("  Inspect 02-factor-atlas.png to decide whether hubs encode broad shape\n")
cat("  or intensity features and leaves encode more class-specific parts.\n")

sink()
close(connection)

message("Done. Results are in ", normalizePath(OUTPUT_DIR))
message("Start with fashion-mnist-miso-report.pdf and interpretation.txt")
