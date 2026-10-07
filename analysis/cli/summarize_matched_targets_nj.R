#!/usr/bin/env Rscript

# NJ accuracy with recorder TARGETS matched across technologies.
#
# The recorders carry different numbers of targets per integration, so matching
# integrations compares constructs of very different capacity. Each recorder is
# read out at the integration count named in design.csv, which encodes whichever
# budget the run was built around -- matched targets or matched information.
#
# The comparison is capacity-matched but not cost-matched: 50 integrations is a
# far harder construct to build and recover than 2, and nothing here accounts
# for that. It also holds the target count fixed while letting integration count
# vary, so any effect of spreading the same targets over more independent
# integrations is folded into the recorder comparison rather than separated
# from it.
#
# Usage:
#   Rscript summarize_matched_targets_nj.R <run-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments)) {
  stop("Usage: summarize_matched_targets_nj.R <run-dir> [prefix]", call. = FALSE)
}
run_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(run_dir, "matched_targets_nj")
}

# The system-to-integration mapping is read from design.csv in the run
# directory rather than hardcoded, so the same assessment serves the
# target-matched and information-matched designs without a second copy that
# could drift out of step with the run it describes.
design_path <- file.path(run_dir, "design.csv")
if (!file.exists(design_path)) {
  stop("No design.csv in ", run_dir,
       "; it must name system, label, targets_per_integration and integrations.",
       call. = FALSE)
}
matched <- utils::read.csv(design_path, stringsAsFactors = FALSE)
required <- c("system", "label", "targets_per_integration", "integrations")
missing <- setdiff(required, names(matched))
if (length(missing)) {
  stop("design.csv is missing: ", paste(missing, collapse = ", "), call. = FALSE)
}
matched$targets <- matched$targets_per_integration * matched$integrations

read_one <- function(system, integrations) {
  path <- file.path(run_dir, "reconstruction",
                    sprintf("%s_k_%d", system, integrations),
                    "reconstruction_summary.csv.gz")
  if (!file.exists(path)) {
    warning("Missing ", path, call. = FALSE)
    return(NULL)
  }
  results <- read.delim(gzfile(path), sep = ",", stringsAsFactors = FALSE)
  results <- results[results$system == system & results$method == "nj", ,
                     drop = FALSE]
  if ("status" %in% names(results)) {
    failed <- sum(results$status != "success")
    if (failed) {
      warning(sprintf("%s: %d of %d NJ runs failed", system, failed,
                      nrow(results)), call. = FALSE)
    }
    results <- results[results$status == "success", , drop = FALSE]
  }
  results
}

rows <- lapply(seq_len(nrow(matched)), function(index) {
  results <- read_one(matched$system[index], matched$integrations[index])
  if (is.null(results) || !nrow(results)) return(NULL)
  value_column <- intersect(c("normalized_rf", "normalised_rf",
                              "mean_normalized_rf"), names(results))[1]
  if (is.na(value_column)) {
    stop("No normalised RF column in the reconstruction summary.", call. = FALSE)
  }
  values <- as.numeric(results[[value_column]])
  values <- values[is.finite(values)]
  standard_error <- stats::sd(values) / sqrt(length(values))
  character_column <- intersect(c("characters", "character_count"),
                                names(results))[1]
  cbind(
    matched[index, ],
    data.frame(
      replicates = length(values),
      mean_characters = if (!is.na(character_column)) {
        mean(as.numeric(results[[character_column]]))
      } else NA_real_,
      mean_normalized_rf = mean(values),
      sd_normalized_rf = stats::sd(values),
      lower = mean(values) - 1.96 * standard_error,
      upper = mean(values) + 1.96 * standard_error,
      stringsAsFactors = FALSE
    )
  )
})
summary_table <- do.call(rbind, Filter(Negate(is.null), rows))
if (is.null(summary_table)) stop("Nothing to summarise.", call. = FALSE)
summary_table <- summary_table[order(summary_table$mean_normalized_rf), ]

# The header describes the design that was actually run, read from design.csv,
# rather than naming a budget that a second design would silently contradict.
cat(sprintf("\n== NJ accuracy: %s ==\n\n",
            paste(sprintf("%s %dx%d", sub(" .*", "", matched$label),
                          matched$integrations, matched$targets_per_integration),
                  collapse = ", ")))
print(summary_table[, c("label", "targets_per_integration", "integrations",
                        "targets", "replicates", "mean_characters",
                        "mean_normalized_rf", "lower", "upper")],
      row.names = FALSE, digits = 4)

# Informative characters, not configured targets, are what tree building
# actually consumes; a recorder can carry 100 targets and still present far
# fewer usable characters.
if (!all(is.na(summary_table$mean_characters))) {
  cat("\nCharacters surviving the invariant filter, against the targets configured:\n")
  for (index in seq_len(nrow(summary_table))) {
    cat(sprintf("  %-20s %6.1f of %3d (%.0f%%)\n",
                summary_table$label[index],
                summary_table$mean_characters[index],
                summary_table$targets[index],
                100 * summary_table$mean_characters[index] /
                  summary_table$targets[index]))
  }
}

utils::write.csv(summary_table, paste0(output_prefix, "_summary.csv"),
                 row.names = FALSE)
cat(sprintf("\nWrote: %s\n", paste0(output_prefix, "_summary.csv")))
