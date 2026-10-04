#!/usr/bin/env Rscript

# Does getting the editing-rate distribution right change the TREE?
#
# The feature comparison showed the empirical rate model reproduces clone 84's
# homoplasy (0.279 against 0.265 observed) where the fitted gamma doubles it
# (0.551). Homoplasy is what breaks character-based reconstruction, so the
# question is whether that difference shows up in the trees themselves.
#
# The comparison is exactly paired. All three clone 84 configurations were run
# from the same seed and the lineage is drawn independently of the editing, so
# every configuration carries a byte-identical ground-truth tree (verified by
# checksum). A replicate therefore fixes the tree, the 250 sampled cells and the
# 5 integrations read out, and varies only the rate model. Nothing here is
# confounded by tree shape or by which cells were captured.
#
# Replicates resample cells and integrations from ONE simulation per rate model,
# so the spread is sampling variability at fixed truth, not variability over
# independent simulations. It answers "does this change the tree from this
# recording", not "how variable is the effect across recordings".
#
# Two references bound every metric:
#   perfect  neighbour joining on the TRUE patristic distances: the ceiling NJ
#            itself imposes, so the gap to a rate model is the recorder's cost
#   random   a random topology on the same tips: the floor
#
# Simulations are read UNMASKED. The feature comparison masked them with the
# real missingness pattern because every feature there is sensitive to capture;
# the published tree benchmarks do not mask, and masking would confound the rate
# model with the assay.
#
# Arguments:
#   --match-dir=<dir>     clone84_simmatch directory (required)
#   --output-prefix=<p>   destination prefix (required)
#   --integrations=<n>    recorders read out per replicate. Default 5.
#   --cells=<n>           cells per tree. Default 250.
#   --replicates=<n>      Default 10.
#   --models=<list>       Default empirical_rates,fitted_shape
#   --clone-counts=<list> cluster counts for the ARI. Default 5,10,20,25.
#   --seed=<n>            Default 1.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
match_dir <- value_after("--match-dir")
output_prefix <- value_after("--output-prefix")
if (is.null(match_dir) || is.null(output_prefix)) {
  stop("--match-dir and --output-prefix are required.", call. = FALSE)
}
match_dir <- normalizePath(match_dir, mustWork = TRUE)
integration_count <- as.integer(value_after("--integrations", "5"))
cell_count <- as.integer(value_after("--cells", "250"))
replicates <- as.integer(value_after("--replicates", "10"))
models <- strsplit(value_after("--models", "empirical_rates,fitted_shape"),
                   ",", fixed = TRUE)[[1L]]
seed <- as.integer(value_after("--seed", "1"))
clone_counts <- as.integer(strsplit(value_after("--clone-counts", "5,10,20,25"),
                                    ",", fixed = TRUE)[[1L]])

for (package in c("ape", "phangorn", "mclust", "data.table")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
source(file.path(dirname(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1])), "..", "tree_metrics.R"))

# ---- load once ---------------------------------------------------------------
tree_paths <- file.path(match_dir, models, "physicell_lineage_sampled.nwk")
if (!all(file.exists(tree_paths))) {
  stop("Missing lineage trees for: ",
       paste(models[!file.exists(tree_paths)], collapse = ", "), call. = FALSE)
}
checksums <- vapply(tree_paths, function(path) {
  as.character(tools::md5sum(path))
}, character(1))
if (length(unique(checksums)) != 1L) {
  stop("The rate models do not share a ground-truth tree; the paired design ",
       "this script assumes does not hold.", call. = FALSE)
}
cat(sprintf("[trees] one ground truth shared by %d rate models\n", length(models)))
# Division waypoints leave unary nodes, and a unary root is promoted to a tip by
# anything that unroots; collapsing is topology-preserving.
truth_full <- ape::collapse.singles(ape::read.tree(tree_paths[[1L]]))

matrices <- lapply(models, function(model) {
  path <- file.path(match_dir, model, "barcode_binary_score_matrix.csv.gz")
  score <- data.table::fread(path, showProgress = FALSE)
  labels <- score[[1L]]
  data.table::set(score, j = 1L, value = NULL)
  x <- as.matrix(score)
  storage.mode(x) <- "integer"
  x[x != 0L] <- 1L
  rownames(x) <- labels
  x
})
names(matrices) <- models
widths <- vapply(matrices, ncol, numeric(1))
if (length(unique(widths)) != 1L) stop("Rate models differ in width.", call. = FALSE)
integrations_of <- sub("_pos_.*$", "", colnames(matrices[[1L]]))
units <- unique(integrations_of)
cat(sprintf("[trees] %d cells x %d sites over %d integrations\n",
            nrow(matrices[[1L]]), ncol(matrices[[1L]]), length(units)))
if (integration_count > length(units)) {
  stop("Only ", length(units), " integrations available.", call. = FALSE)
}

# ---- one replicate -----------------------------------------------------------
# Informative characters only. The filter is on DISTINCT STATES, not on the
# edited count: a character every cell shares carries no split, and one that is
# constant is not a character at all.
informative <- function(x) {
  apply(x, 2L, function(column) length(unique(column)) > 1L)
}
build_tree <- function(x) {
  keep <- informative(x)
  if (sum(keep) < 2L) return(NULL)
  x <- x[, keep, drop = FALSE]
  distance <- stats::dist(x, method = "manhattan") / ncol(x)
  if (any(!is.finite(distance))) return(NULL)
  tryCatch(ape::nj(distance), error = function(e) NULL)
}

set.seed(seed)
rows <- list()
era_rows <- list()
for (replicate in seq_len(replicates)) {
  cells <- sort(sample(nrow(matrices[[1L]]), cell_count))
  chosen_units <- sort(sample(units, integration_count))
  columns <- which(integrations_of %in% chosen_units)
  labels <- rownames(matrices[[1L]])[cells]
  truth <- ape::keep.tip(truth_full, labels)
  tips <- truth$tip.label

  emit <- function(model, tree, characters = NA_integer_) {
    if (is.null(tree)) {
      cat(sprintf("[rep %02d] %s: no tree\n", replicate, model))
      return(invisible(NULL))
    }
    tree <- ape::keep.tip(tree, tips)
    scored <- score_pair(truth, tree, tips, clone_counts = clone_counts)
    era <- scored$era
    scored$era <- NULL
    rows[[length(rows) + 1L]] <<- cbind(
      data.frame(model = model, replicate = replicate,
                 integrations = integration_count, cells = cell_count,
                 characters = characters, stringsAsFactors = FALSE),
      as.data.frame(scored, stringsAsFactors = FALSE))
    if (!is.null(era)) {
      era$model <- model
      era$replicate <- replicate
      era_rows[[length(era_rows) + 1L]] <<- era
    }
    invisible(NULL)
  }

  for (model in models) {
    block <- matrices[[model]][cells, columns, drop = FALSE]
    emit(model, root_estimate(build_tree(block)), sum(informative(block)))
  }
  # References depend only on the truth tree, which both models share.
  emit("perfect", root_estimate(ape::nj(stats::as.dist(
    ape::cophenetic.phylo(truth)[tips, tips]))))
  emit("random", ape::rtree(length(tips), tip.label = tips))
  cat(sprintf("[rep %02d] %s\n", replicate,
              paste(chosen_units, collapse = " ")))
}

scores <- do.call(rbind, rows)
if (is.null(scores)) stop("Nothing scored.", call. = FALSE)
scores$topology_accuracy <- 1 - scores$unrooted_rf
utils::write.csv(scores, paste0(output_prefix, "_per_replicate.csv"),
                 row.names = FALSE)

metric_columns <- c("topology_accuracy", "unrooted_rf", "rooted_rf",
                    "clade_precision", "clade_recall", "triplet_agreement",
                    "ordering_agreement", "cophenetic_spearman",
                    sprintf("clone_ari_k%d", clone_counts))
summary_table <- do.call(rbind, lapply(split(scores, scores$model),
                                       function(part) {
  values <- lapply(metric_columns, function(column) mean(part[[column]],
                                                         na.rm = TRUE))
  spreads <- lapply(metric_columns, function(column) stats::sd(part[[column]],
                                                               na.rm = TRUE))
  names(values) <- metric_columns
  names(spreads) <- paste0(metric_columns, "_sd")
  cbind(data.frame(model = part$model[1], replicates = nrow(part),
                   stringsAsFactors = FALSE),
        as.data.frame(values), as.data.frame(spreads))
}))
utils::write.csv(summary_table, paste0(output_prefix, "_summary.csv"),
                 row.names = FALSE)

show <- function(columns, title) {
  cat(sprintf("\n== %s ==\n\n", title))
  part <- summary_table[, c("model", "replicates", columns)]
  print(part[order(-part[[columns[1]]]), ], row.names = FALSE, digits = 3)
}
cat(sprintf("\n%d replicates: %d cells, %d BASELINE recorders each\n",
            replicates, cell_count, integration_count))
show(c("topology_accuracy", "unrooted_rf"), "Topology accuracy (1 - normalised RF)")
show(c("triplet_agreement", "ordering_agreement"), "Rooted ancestry")
show(sprintf("clone_ari_k%d", clone_counts), "Clone assignment (ARI)")

cat("\n== Informative characters (of 1360 sites over 5 recorders) ==\n\n")
characters <- scores[!is.na(scores$characters), ]
print(aggregate(characters ~ model, data = characters,
                FUN = function(v) c(mean = mean(v), sd = stats::sd(v))),
      row.names = FALSE, digits = 4)

# Deep splits and recent splits fail for different reasons. A rate model that
# front-loads its editing resolves the top of the tree and runs out of variable
# characters near the tips; one that edits late does the reverse. A single RF
# number averages the two away.
era_scores <- do.call(rbind, era_rows)
if (!is.null(era_scores)) {
  era_summary <- do.call(rbind, lapply(
    split(era_scores, list(era_scores$model, era_scores$era), drop = TRUE),
    function(part) data.frame(model = part$model[1], era = part$era[1],
                              clades = nrow(part),
                              recall = mean(part$recovered),
                              stringsAsFactors = FALSE)))
  utils::write.csv(era_summary, paste0(output_prefix, "_by_era.csv"),
                   row.names = FALSE)
  cat("\n== Clade recall by when the split arose ==\n\n")
  wide <- stats::reshape(era_summary[, c("model", "era", "recall")],
                         idvar = "model", timevar = "era", direction = "wide")
  print(wide, row.names = FALSE, digits = 3)
}

# Paired tests: every replicate scored both models on the same tree, the same
# cells and the same integrations, so the difference is the estimate and a
# paired test is the right one.
if (length(models) == 2L) {
  cat("\n== Paired comparison, replicate by replicate ==\n\n")
  first <- scores[scores$model == models[1], ]
  second <- scores[scores$model == models[2], ]
  merged <- merge(first, second, by = "replicate", suffixes = c("_a", "_b"))
  for (metric in c("topology_accuracy", "triplet_agreement",
                   sprintf("clone_ari_k%d", clone_counts))) {
    a <- merged[[paste0(metric, "_a")]]
    b <- merged[[paste0(metric, "_b")]]
    test <- stats::wilcox.test(a, b, paired = TRUE, exact = FALSE)
    cat(sprintf("  %-20s %s %.3f vs %s %.3f  (diff %+.3f, wins %d/%d, p = %.4f)\n",
                metric, models[1], mean(a), models[2], mean(b), mean(a - b),
                sum(a > b), length(a), test$p.value))
  }
}
cat(sprintf("\nWrote:\n  %s\n  %s\n", paste0(output_prefix, "_summary.csv"),
            paste0(output_prefix, "_per_replicate.csv")))
