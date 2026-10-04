#!/usr/bin/env Rscript

# Lineage against spatial organisation in the simulated organoids.
#
# This is the organoid counterpart of the Visium spatial-lineage analysis in
# analysis/physicell_visium.R. That code could not be reused directly: it reads
# a PhysiCell run's node table and indexes MRCAs by binary lifting, whereas the
# organoid ground truth ships per-clone Newick trees and a flat cell table. The
# quantities are the same, computed from what the organoids actually provide.
#
# Two questions, at two scales:
#
#   between clones  are cells that sit near each other more likely to share a
#                   founder? This asks whether clones occupy coherent territory
#                   rather than being stirred through the tissue.
#   within a clone  among cells of one founder, does spatial separation track
#                   lineage separation? This is the finer question: a clone can
#                   hold territory while its internal branches are spatially
#                   scrambled.
#
# Both are reported against a permutation null that shuffles cell positions,
# which preserves the lineage structure and the spatial point pattern exactly
# and destroys only the correspondence between them. Without it a positive
# correlation is uninterpretable, because cells that divide recently are near
# each other for reasons of geometry alone.
#
# Lineage distance within a clone is the cophenetic distance on the truth tree,
# in the tree's own units (minutes of simulated time). Spatial distance is 3D
# Euclidean in microns.
#
# Arguments:
#   --ground-truth=<dir>  lineage_ground_truth directory (required)
#   --output-dir=<dir>    destination (required)
#   --organoids=<n>       organoids to sample. Default 30.
#   --pairs=<n>           between-clone pairs sampled per organoid. Default 200000.
#   --min-clone=<n>       smallest clone for the within-clone analysis. Default 8.
#   --permutations=<n>    position shuffles. Default 20.
#   --seed=<n>            Default 1.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
ground_truth <- value_after("--ground-truth")
output_dir <- value_after("--output-dir")
if (is.null(ground_truth) || is.null(output_dir)) {
  stop("--ground-truth and --output-dir are required.", call. = FALSE)
}
ground_truth <- normalizePath(ground_truth, mustWork = TRUE)
organoid_count <- as.integer(value_after("--organoids", "30"))
pair_count <- as.integer(value_after("--pairs", "200000"))
min_clone <- as.integer(value_after("--min-clone", "8"))
permutations <- as.integer(value_after("--permutations", "20"))
seed <- as.integer(value_after("--seed", "1"))

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The ape package is required.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

distance_breaks <- c(0, 10, 20, 30, 40, 60, 80, 120, Inf)
bin_labels <- c("0-10", "10-20", "20-30", "30-40", "40-60", "60-80",
                "80-120", ">120")

set.seed(seed)
available <- list.dirs(ground_truth, recursive = FALSE, full.names = FALSE)
available <- available[file.exists(file.path(ground_truth, available,
                                             "cells.csv"))]
if (!length(available)) stop("No organoids found.", call. = FALSE)
selected <- if (length(available) > organoid_count) {
  sort(sample(available, organoid_count))
} else available
cat(sprintf("[spatial] %d organoids\n", length(selected)))

between_rows <- list()
within_rows <- list()
clone_rows <- list()

for (organoid in selected) {
  source_dir <- file.path(ground_truth, organoid)
  cells <- utils::read.csv(file.path(source_dir, "cells.csv"),
                           stringsAsFactors = FALSE)
  clones <- utils::read.csv(file.path(source_dir, "clones.csv"),
                            stringsAsFactors = FALSE)
  if (!all(c("x", "y", "z", "founder", "cell_id") %in% names(cells))) next
  positions <- as.matrix(cells[, c("x", "y", "z")])

  # ---- between clones: does proximity predict a shared founder? ----
  pairs_wanted <- min(pair_count, choose(nrow(cells), 2))
  first <- sample.int(nrow(cells), pairs_wanted, replace = TRUE)
  second <- sample.int(nrow(cells), pairs_wanted, replace = TRUE)
  distinct <- first != second
  first <- first[distinct]
  second <- second[distinct]
  spatial <- sqrt(rowSums((positions[first, , drop = FALSE] -
                             positions[second, , drop = FALSE])^2))
  same_founder <- cells$founder[first] == cells$founder[second]
  bin <- cut(spatial, breaks = distance_breaks, labels = bin_labels,
             include.lowest = TRUE)

  # Null: shuffle which cell sits where. Lineage labels and the point cloud are
  # both untouched, so any surviving signal is geometry rather than ancestry.
  null_fractions <- vapply(seq_len(permutations), function(index) {
    shuffled <- sample.int(nrow(cells))
    mean(cells$founder[shuffled[first]] == cells$founder[shuffled[second]])
  }, numeric(1))

  for (level in bin_labels) {
    chosen <- which(bin == level)
    if (!length(chosen)) next
    between_rows[[length(between_rows) + 1L]] <- data.frame(
      organoid = organoid, distance_bin = level, pairs = length(chosen),
      mean_distance = mean(spatial[chosen]),
      same_founder_fraction = mean(same_founder[chosen]),
      null_same_founder_fraction = mean(null_fractions),
      stringsAsFactors = FALSE)
  }

  # ---- within a clone: does spatial separation track lineage separation? ----
  truth_trees <- ape::read.tree(file.path(source_dir, "truth_clones.nwk"))
  if (inherits(truth_trees, "phylo")) truth_trees <- list(truth_trees)
  scorable <- clones[clones$tree_line >= 0 &
                       clones$n_alive_leaves >= min_clone, , drop = FALSE]
  for (index in seq_len(nrow(scorable))) {
    founder <- scorable$founder[index]
    tree <- truth_trees[[scorable$tree_line[index] + 1L]]
    if (is.null(tree)) next
    # Division waypoints leave unary nodes; collapsing is topology-preserving.
    tree <- ape::collapse.singles(tree)
    member <- cells$cell_id %in% tree$tip.label & cells$founder == founder
    if (sum(member) < min_clone) next
    ids <- cells$cell_id[member]
    shared <- intersect(ids, tree$tip.label)
    if (length(shared) < min_clone) next
    tree <- ape::keep.tip(tree, shared)
    order_index <- match(shared, cells$cell_id)
    clone_positions <- positions[order_index, , drop = FALSE]
    lineage <- ape::cophenetic.phylo(tree)[shared, shared]
    spatial_matrix <- as.matrix(stats::dist(clone_positions))
    upper <- upper.tri(lineage)
    lineage_values <- lineage[upper]
    spatial_values <- spatial_matrix[upper]
    if (length(lineage_values) < 10L ||
        stats::sd(lineage_values) == 0) next
    observed <- suppressWarnings(stats::cor(spatial_values, lineage_values,
                                            method = "spearman"))
    null_values <- vapply(seq_len(permutations), function(iteration) {
      shuffled <- sample.int(length(shared))
      permuted <- as.matrix(stats::dist(clone_positions[shuffled, ,
                                                        drop = FALSE]))
      suppressWarnings(stats::cor(permuted[upper], lineage_values,
                                  method = "spearman"))
    }, numeric(1))
    clone_rows[[length(clone_rows) + 1L]] <- data.frame(
      organoid = organoid, founder = founder, cells = length(shared),
      spearman = observed, null_mean = mean(null_values, na.rm = TRUE),
      null_sd = stats::sd(null_values, na.rm = TRUE),
      stringsAsFactors = FALSE)
    clone_bin <- cut(spatial_values, breaks = distance_breaks,
                     labels = bin_labels, include.lowest = TRUE)
    for (level in bin_labels) {
      chosen <- which(clone_bin == level)
      if (!length(chosen)) next
      within_rows[[length(within_rows) + 1L]] <- data.frame(
        organoid = organoid, founder = founder, distance_bin = level,
        pairs = length(chosen), mean_distance = mean(spatial_values[chosen]),
        mean_lineage_distance = mean(lineage_values[chosen]),
        stringsAsFactors = FALSE)
    }
  }
  cat(sprintf("[spatial] %s done\n", organoid))
}

between <- do.call(rbind, between_rows)
within <- do.call(rbind, within_rows)
clone_level <- do.call(rbind, clone_rows)
if (is.null(between)) stop("Nothing scored.", call. = FALSE)

summarise_bins <- function(data, value_column) {
  do.call(rbind, lapply(split(data, data$distance_bin), function(part) {
    values <- part[[value_column]]
    standard_error <- stats::sd(values) / sqrt(nrow(part))
    data.frame(distance_bin = part$distance_bin[1], groups = nrow(part),
               pairs = sum(part$pairs), mean_distance = mean(part$mean_distance),
               value = mean(values),
               lower = mean(values) - 1.96 * standard_error,
               upper = mean(values) + 1.96 * standard_error,
               stringsAsFactors = FALSE)
  }))
}
between_summary <- summarise_bins(between, "same_founder_fraction")
between_summary$null_value <- vapply(
  split(between$null_same_founder_fraction, between$distance_bin), mean,
  numeric(1))[between_summary$distance_bin]
between_summary$distance_bin <- factor(between_summary$distance_bin,
                                       levels = bin_labels)
between_summary <- between_summary[order(between_summary$distance_bin), ]

utils::write.csv(between_summary, file.path(output_dir,
                                            "between_clone_by_distance.csv"),
                 row.names = FALSE)
if (!is.null(within)) {
  within_summary <- summarise_bins(within, "mean_lineage_distance")
  within_summary$distance_bin <- factor(within_summary$distance_bin,
                                        levels = bin_labels)
  within_summary <- within_summary[order(within_summary$distance_bin), ]
  utils::write.csv(within_summary,
                   file.path(output_dir, "within_clone_by_distance.csv"),
                   row.names = FALSE)
}
if (!is.null(clone_level)) {
  utils::write.csv(clone_level, file.path(output_dir, "clone_correlations.csv"),
                   row.names = FALSE)
}

cat("\n== Between clones: shared-founder fraction by spatial distance ==\n\n")
print(between_summary[, c("distance_bin", "pairs", "mean_distance", "value",
                          "lower", "upper", "null_value")],
      row.names = FALSE, digits = 3)
if (!is.null(within)) {
  cat("\n== Within clones: mean lineage distance by spatial distance ==\n\n")
  print(within_summary[, c("distance_bin", "groups", "pairs", "mean_distance",
                           "value", "lower", "upper")],
        row.names = FALSE, digits = 4)
}
if (!is.null(clone_level)) {
  observed_mean <- mean(clone_level$spearman, na.rm = TRUE)
  null_mean <- mean(clone_level$null_mean, na.rm = TRUE)
  test <- stats::t.test(clone_level$spearman, clone_level$null_mean,
                        paired = TRUE)
  cat(sprintf(
    "\n== Within-clone spatial/lineage correlation ==\n\n  clones %d   observed rho %.4f   permuted rho %.4f   paired p %.3g\n",
    nrow(clone_level), observed_mean, null_mean, test$p.value))
  cat(sprintf("  clones with rho > 0: %.1f%%\n",
              100 * mean(clone_level$spearman > 0, na.rm = TRUE)))
}
cat(sprintf("\nWrote outputs under %s\n", output_dir))
