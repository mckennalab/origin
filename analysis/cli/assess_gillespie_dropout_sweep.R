#!/usr/bin/env Rscript

# Recorder accuracy under integration and cell dropout, on the nested Gillespie
# integration series.
#
# The companion sweep (assess_organoid_dropout_sweep.R) does the same thing to
# the organoid recordings. Dropout semantics are identical between the two so
# the curves are comparable:
#
#   integration dropout  each cell independently fails to recover each of its
#                        integrations with probability p; affected characters
#                        become MISSING (-1), not unedited.
#   cell dropout         each sampled cell is independently not captured; the
#                        ground-truth tree is pruned to the captured cells
#                        before scoring.
#
# Source matrices are the level-100 conditions, i.e. the full recording, so
# dropout removes information from the same realization rather than from
# separately simulated smaller recorders. That distinction matters: the original
# series varies how many integrations EXIST and leaves the matrix complete,
# whereas dropout varies how many are RECOVERED per cell and makes it ragged.
#
# Mitochondrial recording has no integrations -- its characters are mt_<pos>_
# <allele> and its recovery axis is genome depth, which the original series
# already sweeps. It is therefore included only at integration dropout 0, and
# excluded from higher rates rather than having an arbitrary grouping imposed.
#
# Arguments:
#   --series=<dir>        neutral_best_nj_integration_series_nested_n250 root
#   --output-dir=<dir>    destination (required)
#   --integration-dropout=<list>  default 0,0.05,0.1,0.25,0.5
#   --cell-dropout=<list>         default 0,0.1,0.25,0.5,0.8
#   --level=<n>           source condition level. Default 100.
#   --min-leaves=<n>      default 4.
#   --workers=<n>         default 6.
#   --seed=<n>            default 1.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
series_dir <- value_after("--series")
output_dir <- value_after("--output-dir")
if (is.null(series_dir) || is.null(output_dir)) {
  stop("--series and --output-dir are required.", call. = FALSE)
}
series_dir <- normalizePath(series_dir, mustWork = TRUE)
parse_list <- function(text) as.numeric(strsplit(text, ",", fixed = TRUE)[[1L]])
integration_dropout <- parse_list(value_after("--integration-dropout",
                                              "0,0.05,0.1,0.25,0.5"))
cell_dropout <- parse_list(value_after("--cell-dropout", "0,0.1,0.25,0.5,0.8"))
level <- as.integer(value_after("--level", "100"))
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

condition_name <- sprintf("condition_level_%03d", level)
matrices <- list.files(
  file.path(series_dir, "simulation"),
  pattern = "^recording_character_matrix_sparse[.]rds$",
  recursive = TRUE, full.names = TRUE
)
matrices <- matrices[grepl(condition_name, matrices, fixed = TRUE)]
if (!length(matrices)) {
  # Level directories are zero-padded to two digits in some runs.
  condition_name <- sprintf("condition_level_%d", level)
  matrices <- list.files(
    file.path(series_dir, "simulation"),
    pattern = "^recording_character_matrix_sparse[.]rds$",
    recursive = TRUE, full.names = TRUE
  )
  matrices <- matrices[grepl(condition_name, matrices, fixed = TRUE)]
}
if (!length(matrices)) stop("No level-", level, " condition matrices found.",
                            call. = FALSE)

tasks <- data.frame(
  matrix_path = matrices,
  system = sub(".*/recorder_([^/]+)/.*", "\\1", matrices),
  seed_label = sub(".*/(seed_[0-9]+)/.*", "\\1", matrices),
  stringsAsFactors = FALSE
)
grid <- expand.grid(integration_dropout = integration_dropout,
                    cell_dropout = cell_dropout, KEEP.OUT.ATTRS = FALSE)
cat(sprintf("[sweep] %d source recordings (%s), %d conditions\n",
            nrow(tasks), paste(sort(unique(tasks$system)), collapse = ", "),
            nrow(grid)))

#' Score one recording under one dropout condition
#'
#' @param task One row of the task table.
#' @param score_matrix Character matrix, cells in rows.
#' @param truth_tree Ground-truth sampled tree.
#' @param integration_rate,cell_rate Dropout probabilities.
#' @param condition_seed Seed for reproducibility.
#' @return One-row data frame, or NULL when the condition leaves nothing to score.
score_condition <- function(task, score_matrix, truth_tree, integration_rate,
                            cell_rate, condition_seed) {
  integration_of <- sub("_pos_.*$", "", colnames(score_matrix))
  has_integrations <- any(grepl("^int_", colnames(score_matrix)))
  # Mitochondrial characters are not grouped into integrations; imposing a
  # grouping would invent a recovery model the recorder does not have.
  if (integration_rate > 0 && !has_integrations) return(NULL)

  set.seed(condition_seed)
  working <- score_matrix
  if (integration_rate > 0) {
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
  if (cell_rate > 0) {
    keep <- stats::runif(nrow(working)) >= cell_rate
    working <- working[keep, , drop = FALSE]
  }
  if (nrow(working)) {
    working <- working[rowSums(working != -1L) > 0, , drop = FALSE]
  }
  if (nrow(working) < min_leaves) return(NULL)

  observed_count <- colSums(working != -1L)
  edited_count <- colSums(working == 1L)
  informative <- observed_count >= 2L & edited_count > 0L &
    edited_count < observed_count
  if (!any(informative)) return(NULL)
  clone_matrix <- working[, informative, drop = FALSE]

  estimate <- tryCatch({
    cliqueR::build_trees(clone_matrix, methods = "nj", missing = -1L)$nj
  }, error = function(condition) NULL)
  if (is.null(estimate)) return(NULL)

  shared <- intersect(estimate$tip.label, truth_tree$tip.label)
  if (length(shared) < min_leaves) return(NULL)
  distance <- tryCatch(cliqueR::rf_distance(
    ape::keep.tip(estimate, shared), ape::keep.tip(truth_tree, shared),
    normalize = TRUE
  ), error = function(condition) {
    # Announce rather than swallow: a silent NA here once removed an entire
    # condition from the grid without leaving any trace in the output.
    cat(sprintf("[sweep] SCORING FAILED %s %s int=%g cell=%g: %s\n",
                task$system, task$seed_label, integration_rate, cell_rate,
                conditionMessage(condition)))
    NA_real_
  })
  if (!is.finite(distance)) return(NULL)

  data.frame(
    system = task$system, seed_label = task$seed_label,
    integration_dropout = integration_rate, cell_dropout = cell_rate,
    cells = length(shared), informative_characters = sum(informative),
    missing_fraction = mean(clone_matrix == -1L),
    normalized_rf = distance, stringsAsFactors = FALSE
  )
}

results <- parallel::mclapply(seq_len(nrow(tasks)), function(position) {
  task <- tasks[position, ]
  score_matrix <- as.matrix(readRDS(task$matrix_path))
  storage.mode(score_matrix) <- "integer"
  # Every division is kept as an internal node, so the tree is full of unary
  # waypoints and its root has a single child. rf_distance calls ape::unroot(),
  # which promotes a unary root into an extra tip and then rejects the pair on
  # mismatched labels. Collapsing singles preserves topology and summed branch
  # lengths while making the root bifurcating; scores are unchanged where the
  # comparison already worked.
  truth_tree <- ape::collapse.singles(ape::read.tree(file.path(
    dirname(task$matrix_path), "ground_truth_tree.nwk")))
  rows <- lapply(seq_len(nrow(grid)), function(index) {
    score_condition(task, score_matrix, truth_tree,
                    grid$integration_dropout[index], grid$cell_dropout[index],
                    condition_seed = seed + 1000L * position + index)
  })
  cat(sprintf("[sweep] %s %s done\n", task$system, task$seed_label))
  do.call(rbind, Filter(Negate(is.null), rows))
}, mc.cores = workers, mc.preschedule = FALSE)

scores <- do.call(rbind, Filter(Negate(is.null), results))
if (is.null(scores)) stop("Nothing scored.", call. = FALSE)

scores_path <- file.path(output_dir, "dropout_scores.tsv.gz")
connection <- gzfile(scores_path, open = "wt")
write.table(scores, connection, sep = "\t", row.names = FALSE, quote = FALSE)
close(connection)

summary_table <- do.call(rbind, lapply(
  split(scores, list(scores$system, scores$integration_dropout,
                     scores$cell_dropout), drop = TRUE),
  function(part) {
    standard_error <- if (nrow(part) > 1L) {
      stats::sd(part$normalized_rf) / sqrt(nrow(part))
    } else NA_real_
    data.frame(
      system = part$system[1],
      integration_dropout = part$integration_dropout[1],
      cell_dropout = part$cell_dropout[1], replicates = nrow(part),
      mean_cells = mean(part$cells),
      mean_informative_characters = mean(part$informative_characters),
      mean_normalized_rf = mean(part$normalized_rf),
      lower = mean(part$normalized_rf) - 1.96 * standard_error,
      upper = mean(part$normalized_rf) + 1.96 * standard_error,
      stringsAsFactors = FALSE
    )
  }
))
summary_table <- summary_table[order(summary_table$system,
                                     summary_table$cell_dropout,
                                     summary_table$integration_dropout), ]
write.table(summary_table, file.path(output_dir, "dropout_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

jsonlite::write_json(list(
  series = series_dir, level = level, integration_dropout = integration_dropout,
  cell_dropout = cell_dropout, min_leaves = min_leaves, seed = seed,
  note = paste("Dropout applied to level-", level, " recordings; missing as -1;",
               "truth pruned to captured cells; mitochondrial excluded from",
               "integration dropout (no integration grouping).", sep = "")
), file.path(output_dir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)

cat(sprintf("\n%d scores, %d summary rows\n", nrow(scores), nrow(summary_table)))
baseline_rows <- summary_table[summary_table$system == "baseline", ]
print(baseline_rows[, c("integration_dropout", "cell_dropout", "replicates",
                        "mean_cells", "mean_normalized_rf")],
      row.names = FALSE, digits = 3)
cat(sprintf("\nWrote:\n  %s\n  %s\n", scores_path,
            file.path(output_dir, "dropout_summary.tsv")))
