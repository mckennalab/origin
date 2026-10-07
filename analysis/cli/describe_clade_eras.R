#!/usr/bin/env Rscript

# What the clade-recall eras actually are.
#
# The era panels bin true clades into earliest/early/late/latest, and those are
# QUARTILES OF SPLIT TIME computed per tree -- not fixed time thresholds. That
# has two consequences a reader needs stated:
#
#   the bins are not comparable across shapes in absolute time, because each
#   tree's quartiles come from its own depth distribution; and
#
#   early means big, because the trees are ultrametric, so a clade's root-to-node
#   time is when it arose and deep clades contain most of the cells. "Recovers
#   early splits" and "recovers large clades" are nearly the same claim.
#
# Balanced is the one shape that does not divide into four equal bins: its
# branch lengths take only a handful of distinct values, so many nodes tie at
# the same depth and cannot be separated at a quartile boundary.
#
# Usage:
#   Rscript describe_clade_eras.R <run-dir> [output-csv]

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments)) {
  stop("Usage: describe_clade_eras.R <run-dir> [output-csv]", call. = FALSE)
}
run_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
output_path <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(run_dir, "clade_era_definitions.csv")
}
if (!requireNamespace("ape", quietly = TRUE)) {
  stop("The ape package is required.", call. = FALSE)
}

era_levels <- c("earliest", "early", "late", "latest")
rows <- list()
shape_dirs <- list.dirs(file.path(run_dir, "simulation"), recursive = FALSE,
                        full.names = TRUE)
for (shape_dir in shape_dirs) {
  shape <- sub("^shape_", "", basename(shape_dir))
  for (seed_dir in list.dirs(shape_dir, recursive = FALSE, full.names = TRUE)) {
    seed <- basename(seed_dir)
    truth_path <- Sys.glob(file.path(seed_dir, "recorder_*", "parameter_*",
                                     "condition*", "ground_truth_tree.nwk"))
    if (!length(truth_path)) next
    # Unary division waypoints would otherwise be counted as nodes; collapsing
    # is topology-preserving.
    truth <- ape::collapse.singles(ape::read.tree(truth_path[[1L]]))
    parts <- ape::prop.part(truth)
    labels <- attr(parts, "labels")
    sizes <- lengths(parts)
    keep <- which(sizes >= 2L & sizes < length(labels))
    if (!length(keep)) next
    times <- ape::node.depth.edgelength(truth)[length(labels) + keep]
    breaks <- stats::quantile(times, probs = seq(0, 1, 0.25))
    era <- cut(times, breaks = breaks, include.lowest = TRUE,
               labels = era_levels)
    clade_sizes <- sizes[keep]
    for (level in era_levels) {
      members <- which(era == level)
      if (!length(members)) next
      rows[[length(rows) + 1L]] <- data.frame(
        shape = shape, seed = seed, era = level,
        clades = length(members),
        split_time_min = min(times[members]),
        split_time_max = max(times[members]),
        clade_size_min = min(clade_sizes[members]),
        clade_size_median = stats::median(clade_sizes[members]),
        clade_size_max = max(clade_sizes[members]),
        tips = length(labels),
        scorable_clades = length(keep),
        stringsAsFactors = FALSE)
    }
  }
}
table_out <- do.call(rbind, rows)
if (is.null(table_out)) stop("No truth trees found.", call. = FALSE)
table_out$era <- factor(table_out$era, levels = era_levels)
table_out <- table_out[order(table_out$shape, table_out$seed, table_out$era), ]
utils::write.csv(table_out, output_path, row.names = FALSE)

cat("\n== Clade-recall eras, averaged over seeds ==\n\n")
summary_out <- do.call(rbind, lapply(
  split(table_out, list(table_out$shape, table_out$era), drop = TRUE),
  function(part) data.frame(
    shape = part$shape[1], era = part$era[1], seeds = nrow(part),
    mean_clades = mean(part$clades),
    mean_split_time_min = mean(part$split_time_min),
    mean_split_time_max = mean(part$split_time_max),
    mean_clade_size_median = mean(part$clade_size_median),
    stringsAsFactors = FALSE)))
summary_out <- summary_out[order(summary_out$shape, summary_out$era), ]
print(summary_out, row.names = FALSE, digits = 3)
cat(sprintf("\nWrote %s (%d rows)\n", output_path, nrow(table_out)))
