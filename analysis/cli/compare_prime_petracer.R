#!/usr/bin/env Rscript

# Prime editing before and after the PEtracer replacement.
#
# The prime slot in the recorder parameter grid was a single-outcome prime
# editor: 6 targets per integration, one pegRNA each, so an edit was recorded as
# presence/absence. It now carries a PEtracer-style recorder -- 3 edit sites per
# integration, each installing one of 8 predefined marks -- at the same 5
# integrations as the rest of the panel.
#
# Two things changed together and the comparison cannot separate them:
#
#   the alphabet   an edited site now says WHICH mark it received, so two cells
#                  that edited a site independently are usually distinguishable
#                  rather than looking like shared ancestry
#   the geometry   3 sites x 8 marks against 6 binary sites, and a different
#                  rate ladder anchored on edit saturation
#
# Reporting them as one change is the honest framing: this is "the prime slot
# now models PEtracer", not "adding marks is worth X".
#
# Usage:
#   Rscript compare_prime_petracer.R <grid-dir> <backup-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2L) {
  stop("Usage: compare_prime_petracer.R <grid-dir> <backup-dir> [prefix]",
       call. = FALSE)
}
grid_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
backup_dir <- normalizePath(arguments[[2L]], mustWork = TRUE)
output_prefix <- if (length(arguments) >= 3L) arguments[[3L]] else {
  file.path(grid_dir, "prime_petracer_comparison")
}

read_accuracy <- function(root, file) {
  path <- file.path(root, file)
  if (!file.exists(path)) {
    stop("Missing ", path, call. = FALSE)
  }
  read.delim(path, stringsAsFactors = FALSE)
}
setting_column <- function(data) {
  grep("setting", names(data), value = TRUE)[1]
}
tier_levels <- c("very_low", "low", "mid", "high", "very_high")

new_overall <- read_accuracy(grid_dir, "accuracy_by_recorder_parameter.tsv")
old_overall <- read_accuracy(file.path(backup_dir, "summary_tables"),
                             "accuracy_by_recorder_parameter.tsv")
setting <- setting_column(new_overall)

tier_table <- function(data, label) {
  part <- data[data$system == "prime", , drop = FALSE]
  frame <- data.frame(
    tier = factor(part[[setting]], levels = tier_levels),
    mean_normalized_rf = part$mean_normalized_rf,
    stringsAsFactors = FALSE
  )
  frame <- frame[order(frame$tier), , drop = FALSE]
  stats::setNames(frame$mean_normalized_rf, as.character(frame$tier))
}
old_prime <- tier_table(old_overall)
new_prime <- tier_table(new_overall)
common <- intersect(names(old_prime), names(new_prime))

cat("== Prime slot: single-outcome vs PEtracer ==\n\n")
cat(sprintf("%-11s %10s %10s %10s\n", "tier", "old RF", "new RF", "change"))
for (tier in tier_levels) {
  old_value <- if (tier %in% names(old_prime)) old_prime[[tier]] else NA_real_
  new_value <- if (tier %in% names(new_prime)) new_prime[[tier]] else NA_real_
  cat(sprintf("%-11s %10.4f %10.4f %+10.4f\n", tier, old_value, new_value,
              new_value - old_value))
}
cat(sprintf("\nbest old: %s at %.4f\nbest new: %s at %.4f\n",
            names(which.min(old_prime)), min(old_prime),
            names(which.min(new_prime)), min(new_prime)))
# Accuracy is 1 - RF, so a drop in RF is the improvement.
cat(sprintf("improvement at the optimum: %.4f RF (%.1f%% of the old error)\n",
            min(old_prime) - min(new_prime),
            100 * (min(old_prime) - min(new_prime)) / min(old_prime)))

cat("\n== Against the rest of the panel (best tier per system, new grid) ==\n\n")
best_by_system <- do.call(rbind, lapply(
  split(new_overall, new_overall$system), function(part) {
    best <- part[which.min(part$mean_normalized_rf), , drop = FALSE]
    data.frame(system = best$system[1], tier = best[[setting]][1],
               mean_normalized_rf = best$mean_normalized_rf[1],
               stringsAsFactors = FALSE)
  }
))
best_by_system <- best_by_system[order(best_by_system$mean_normalized_rf), ]
print(best_by_system, row.names = FALSE, digits = 4)

# Per shape, because the ladder was anchored on saturation and the four
# topologies differ in depth; a single pooled number would hide that.
shape_file <- "accuracy_by_recorder_parameter_shape.tsv"
if (file.exists(file.path(grid_dir, shape_file))) {
  new_shape <- read_accuracy(grid_dir, shape_file)
  old_shape <- read_accuracy(file.path(backup_dir, "summary_tables"), shape_file)
  shape_setting <- setting_column(new_shape)
  cat("\n== Prime by tree shape: best tier and mean normalised RF ==\n\n")
  cat(sprintf("%-14s %-11s %8s  %-11s %8s\n",
              "shape", "old tier", "old RF", "new tier", "new RF"))
  for (shape in sort(unique(new_shape$shape))) {
    old_part <- old_shape[old_shape$system == "prime" &
                            old_shape$shape == shape, , drop = FALSE]
    new_part <- new_shape[new_shape$system == "prime" &
                            new_shape$shape == shape, , drop = FALSE]
    if (!nrow(old_part) || !nrow(new_part)) next
    old_best <- old_part[which.min(old_part$mean_normalized_rf), ]
    new_best <- new_part[which.min(new_part$mean_normalized_rf), ]
    cat(sprintf("%-14s %-11s %8.4f  %-11s %8.4f\n", shape,
                old_best[[shape_setting]], old_best$mean_normalized_rf,
                new_best[[shape_setting]], new_best$mean_normalized_rf))
  }
}

method_file <- "accuracy_by_recorder_parameter_method.tsv"
if (file.exists(file.path(grid_dir, method_file))) {
  new_method <- read_accuracy(grid_dir, method_file)
  old_method <- read_accuracy(file.path(backup_dir, "summary_tables"),
                              method_file)
  cat("\n== Prime by tree-building method, at each grid's best tier ==\n\n")
  cat(sprintf("%-12s %10s %10s %10s\n", "method", "old RF", "new RF", "change"))
  for (method in sort(unique(new_method$method))) {
    old_part <- old_method[old_method$system == "prime" &
                             old_method$method == method, , drop = FALSE]
    new_part <- new_method[new_method$system == "prime" &
                             new_method$method == method, , drop = FALSE]
    if (!nrow(old_part) || !nrow(new_part)) next
    old_value <- min(old_part$mean_normalized_rf)
    new_value <- min(new_part$mean_normalized_rf)
    cat(sprintf("%-12s %10.4f %10.4f %+10.4f\n", method, old_value, new_value,
                new_value - old_value))
  }
}

comparison <- data.frame(
  tier = tier_levels,
  old_mean_normalized_rf = as.numeric(old_prime[tier_levels]),
  new_mean_normalized_rf = as.numeric(new_prime[tier_levels]),
  stringsAsFactors = FALSE
)
comparison$change <- comparison$new_mean_normalized_rf -
  comparison$old_mean_normalized_rf
utils::write.csv(comparison, paste0(output_prefix, ".csv"), row.names = FALSE)
cat(sprintf("\nWrote: %s\n", paste0(output_prefix, ".csv")))
