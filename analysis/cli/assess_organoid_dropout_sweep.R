#!/usr/bin/env Rscript

# BASELINE recorder accuracy under integration and cell dropout.
#
# The baseline organoid assessment assumes perfect recovery: every integration
# observed in every cell, every live cell captured, no missing data. This sweep
# relaxes both assumptions over a grid and reports what accuracy survives.
#
# Dropout is applied after the replay, so the recordings written by
# assess_organoid_baseline_nj.R are reused rather than resimulated. Every
# condition therefore starts from the same underlying edits, and differences
# between conditions are attributable to the recovery model alone.
#
# The two kinds of loss are modelled differently because they are different:
#
#   integration dropout  each cell independently fails to recover each of its
#                        integrations with probability p. The affected
#                        characters become MISSING, not unedited -- recording
#                        "no edit seen" where nothing was observed would invent
#                        evidence of shared ancestry. Missing is passed to
#                        cliqueR as -1.
#   cell dropout         each live cell is independently not captured with
#                        probability p. The ground-truth tree is pruned to the
#                        captured cells before scoring, since comparing a
#                        subsampled reconstruction against the full truth tree
#                        would count absent tips as errors.
#
# Arguments:
#   --recordings=<dir>    directory of per-organoid replays (required)
#   --ground-truth=<dir>  lineage_ground_truth directory (required)
#   --output-dir=<dir>    destination (required)
#   --integration-dropout=<list>  default 0,0.05,0.1,0.25,0.5
#   --cell-dropout=<list>         default 0,0.1,0.25,0.5,0.8
#   --organoids=<n>       organoids to use. Default 50.
#   --min-leaves=<n>      smallest clone to score, after dropout. Default 4.
#   --workers=<n>         default 6.
#   --seed=<n>            default 1.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
recordings_dir <- value_after("--recordings")
ground_truth <- value_after("--ground-truth")
output_dir <- value_after("--output-dir")
if (is.null(recordings_dir) || is.null(ground_truth) || is.null(output_dir)) {
  stop("--recordings, --ground-truth and --output-dir are required.",
       call. = FALSE)
}
recordings_dir <- normalizePath(recordings_dir, mustWork = TRUE)
ground_truth <- normalizePath(ground_truth, mustWork = TRUE)
parse_list <- function(text) as.numeric(strsplit(text, ",", fixed = TRUE)[[1L]])
integration_dropout <- parse_list(value_after("--integration-dropout",
                                              "0,0.05,0.1,0.25,0.5"))
cell_dropout <- parse_list(value_after("--cell-dropout", "0,0.1,0.25,0.5,0.8"))
organoid_count <- as.integer(value_after("--organoids", "50"))
min_leaves <- as.integer(value_after("--min-leaves", "4"))
workers <- as.integer(value_after("--workers", "6"))
seed <- as.integer(value_after("--seed", "1"))

for (package in c("Matrix", "ape", "cliqueR", "parallel", "jsonlite")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

available <- list.dirs(recordings_dir, recursive = FALSE, full.names = FALSE)
available <- available[file.exists(file.path(
  recordings_dir, available, "barcode_binary_score_matrix_sparse.rds"
))]
if (!length(available)) stop("No recordings found.", call. = FALSE)
set.seed(seed)
selected <- if (length(available) > organoid_count) {
  sort(sample(available, organoid_count))
} else {
  available
}
grid <- expand.grid(integration_dropout = integration_dropout,
                    cell_dropout = cell_dropout,
                    KEEP.OUT.ATTRS = FALSE)
cat(sprintf("[sweep] %d organoids x %d conditions\n", length(selected),
            nrow(grid)))

#' Score every clone of one organoid under one dropout condition
#'
#' @param organoid_id Recording directory name.
#' @param score_matrix Full binary matrix, rows named c<ID>.
#' @param cells,clones,truth_trees Ground-truth tables and trees.
#' @param integration_rate,cell_rate Dropout probabilities.
#' @param condition_seed Seed for this organoid-condition pair.
#' @return Data frame with one row per scored clone.
score_condition <- function(organoid_id, score_matrix, cells, clones,
                            truth_trees, integration_rate, cell_rate,
                            condition_seed) {
  set.seed(condition_seed)
  working <- score_matrix

  if (integration_rate > 0) {
    # Each cell loses each integration independently. Whole integrations drop
    # together because that is how recovery fails in practice: an integration is
    # amplified and sequenced, or it is not.
    integration_of <- sub("_pos_.*$", "", colnames(working))
    integrations <- unique(integration_of)
    lost <- matrix(
      stats::runif(nrow(working) * length(integrations)) < integration_rate,
      nrow = nrow(working), ncol = length(integrations)
    )
    for (index in seq_along(integrations)) {
      rows <- which(lost[, index])
      if (length(rows)) {
        working[rows, integration_of == integrations[index]] <- -1L
      }
    }
  }

  retained_cells <- rownames(working)
  if (cell_rate > 0) {
    keep <- stats::runif(nrow(working)) >= cell_rate
    retained_cells <- rownames(working)[keep]
    working <- working[keep, , drop = FALSE]
  }
  # A cell that lost every integration carries no data at all; it is
  # unrecoverable in the same sense as an uncaptured cell.
  if (nrow(working)) {
    observed_per_cell <- rowSums(working != -1L)
    working <- working[observed_per_cell > 0, , drop = FALSE]
  }
  if (nrow(working) < min_leaves) return(NULL)

  scorable <- clones[clones$tree_line >= 0, , drop = FALSE]
  rows <- lapply(seq_len(nrow(scorable)), function(index) {
    founder <- scorable$founder[index]
    clone_cells <- cells$cell_id[cells$founder == founder]
    present <- rownames(working) %in% clone_cells
    if (sum(present) < min_leaves) return(NULL)
    clone_matrix <- working[present, , drop = FALSE]

    # Informative among OBSERVED entries only: a character seen in just one
    # cell says nothing about how the others relate.
    observed <- clone_matrix != -1L
    edited <- clone_matrix == 1L
    observed_count <- colSums(observed)
    edited_count <- colSums(edited)
    informative <- observed_count >= 2L & edited_count > 0L &
      edited_count < observed_count
    if (!any(informative)) return(NULL)
    clone_matrix <- clone_matrix[, informative, drop = FALSE]

    truth_tree <- truth_trees[[scorable$tree_line[index] + 1L]]
    estimate <- tryCatch({
      cliqueR::build_trees(clone_matrix, methods = "nj", missing = -1L)$nj
    }, error = function(condition) NULL)
    if (is.null(estimate)) return(NULL)

    # Prune the truth to the captured cells: absent tips are not errors.
    shared <- intersect(estimate$tip.label, truth_tree$tip.label)
    if (length(shared) < min_leaves) return(NULL)
    distance <- tryCatch(cliqueR::rf_distance(
      ape::keep.tip(estimate, shared), ape::keep.tip(truth_tree, shared),
      normalize = TRUE
    ), error = function(condition) {
      # Announce rather than swallow: a silent NA here once removed an entire
      # condition from the grid without leaving any trace in the output.
      cat(sprintf("[sweep] SCORING FAILED %s founder %s int=%g cell=%g: %s\n",
                  organoid_id, founder, integration_rate, cell_rate,
                  conditionMessage(condition)))
      NA_real_
    })
    if (!is.finite(distance)) return(NULL)

    data.frame(
      organoid_id = organoid_id, founder = founder,
      integration_dropout = integration_rate, cell_dropout = cell_rate,
      leaves = length(shared), full_leaves = scorable$n_alive_leaves[index],
      informative_characters = sum(informative),
      missing_fraction = mean(clone_matrix == -1L),
      normalized_rf = distance, stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

results <- parallel::mclapply(seq_along(selected), function(position) {
  organoid_id <- selected[position]
  source_dir <- file.path(ground_truth, organoid_id)
  score_matrix <- as.matrix(readRDS(file.path(
    recordings_dir, organoid_id, "barcode_binary_score_matrix_sparse.rds"
  )))
  rownames(score_matrix) <- sub("^cell_", "c", rownames(score_matrix))
  storage.mode(score_matrix) <- "integer"
  cells <- read.csv(file.path(source_dir, "cells.csv"), stringsAsFactors = FALSE)
  clones <- read.csv(file.path(source_dir, "clones.csv"),
                     stringsAsFactors = FALSE)
  # Collapse unary internal nodes, for the reason given in the Gillespie sweep:
  # rf_distance unroots, and unrooting a tree whose root has one child adds a
  # spurious tip. Topology-preserving, so scores are unchanged where the
  # comparison already succeeded.
  truth_trees <- lapply(ape::read.tree(file.path(source_dir, "truth_clones.nwk")),
                        ape::collapse.singles)

  per_condition <- lapply(seq_len(nrow(grid)), function(index) {
    score_condition(
      organoid_id, score_matrix, cells, clones, truth_trees,
      grid$integration_dropout[index], grid$cell_dropout[index],
      condition_seed = seed + 1000L * position + index
    )
  })
  cat(sprintf("[sweep] %s done\n", organoid_id))
  do.call(rbind, Filter(Negate(is.null), per_condition))
}, mc.cores = workers, mc.preschedule = FALSE)

clone_scores <- do.call(rbind, Filter(Negate(is.null), results))
if (is.null(clone_scores)) stop("No clone scored under any condition.",
                                call. = FALSE)

scores_path <- file.path(output_dir, "dropout_clone_scores.tsv.gz")
connection <- gzfile(scores_path, open = "wt")
write.table(clone_scores, connection, sep = "\t", row.names = FALSE,
            quote = FALSE)
close(connection)

summary_table <- do.call(rbind, lapply(
  split(clone_scores, list(clone_scores$integration_dropout,
                           clone_scores$cell_dropout), drop = TRUE),
  function(part) {
    standard_error <- stats::sd(part$normalized_rf) / sqrt(nrow(part))
    data.frame(
      integration_dropout = part$integration_dropout[1],
      cell_dropout = part$cell_dropout[1],
      clones = nrow(part), mean_leaves = mean(part$leaves),
      mean_informative_characters = mean(part$informative_characters),
      mean_normalized_rf = mean(part$normalized_rf),
      mean_accuracy = 1 - mean(part$normalized_rf),
      lower = mean(part$normalized_rf) - 1.96 * standard_error,
      upper = mean(part$normalized_rf) + 1.96 * standard_error,
      exact_fraction = mean(part$normalized_rf == 0),
      stringsAsFactors = FALSE
    )
  }
))
summary_table <- summary_table[order(summary_table$cell_dropout,
                                     summary_table$integration_dropout), ]
write.table(summary_table, file.path(output_dir, "dropout_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

jsonlite::write_json(list(
  recordings = recordings_dir, ground_truth = ground_truth,
  organoids = length(selected), integration_dropout = integration_dropout,
  cell_dropout = cell_dropout, min_leaves = min_leaves, seed = seed,
  note = paste("Dropout applied post-replay to shared recordings;",
               "missing encoded as -1; truth pruned to captured cells.")
), file.path(output_dir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)

cat(sprintf("\n%d clone scores across %d conditions\n", nrow(clone_scores),
            nrow(summary_table)))
print(summary_table[, c("integration_dropout", "cell_dropout", "clones",
                        "mean_leaves", "mean_normalized_rf", "exact_fraction")],
      row.names = FALSE, digits = 3)
cat(sprintf("\nWrote:\n  %s\n  %s\n", scores_path,
            file.path(output_dir, "dropout_summary.tsv")))
