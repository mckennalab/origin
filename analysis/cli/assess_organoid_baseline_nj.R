#!/usr/bin/env Rscript

# BASELINE recorder accuracy on PhysiCell organoid ground-truth lineages.
#
# For each exported organoid this replays a 5-integration BASELINE base-editing
# recorder over the real division history, then reconstructs each founder clone
# independently with neighbour joining and scores it against that clone's
# ground-truth tree. The organoid lineage is known exactly, so the only error
# measured is what the recorder plus the reconstruction lose.
#
# Clones are scored separately rather than as one organoid tree because each
# founder is an independent day-zero cell: the clones share no ancestry after
# t = 0, so a joined tree would score the trivial between-clone splits along
# with the within-clone structure that is actually at issue.
#
# Arguments:
#   --ground-truth=<dir>  lineage_ground_truth directory (required)
#   --params=<path>       BASELINE parameter JSON (required)
#   --output-dir=<dir>    destination (required)
#   --integrations=<n>    barcode integrations per cell. Default 5.
#   --end-time=<x>        replay horizon in minutes. Default 10080 (7 days).
#   --min-leaves=<n>      smallest clone to score. Default 4; below four tips an
#                         unrooted topology has no internal split, so normalised
#                         RF is 0 by construction and carries no information.
#   --organoids=<list>    comma-separated organoid ids. Default all.
#   --limit=<n>           score only the first n organoids. Default all.
#   --workers=<n>         parallel organoids. Default 4.
#   --seed=<n>            recorder seed. Default 1.
#
# Writes clone_accuracy.tsv.gz (one row per scored clone), organoid_summary.tsv
# and settings.json.

arguments <- commandArgs(trailingOnly = TRUE)
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
ground_truth <- value_after("--ground-truth")
params_path <- value_after("--params")
output_dir <- value_after("--output-dir")
if (is.null(ground_truth) || is.null(params_path) || is.null(output_dir)) {
  stop("--ground-truth, --params and --output-dir are all required.",
       call. = FALSE)
}
ground_truth <- normalizePath(ground_truth, mustWork = TRUE)
params_path <- normalizePath(params_path, mustWork = TRUE)
integrations <- as.integer(value_after("--integrations", "5"))
end_time <- as.numeric(value_after("--end-time", "10080"))
min_leaves <- as.integer(value_after("--min-leaves", "4"))
requested <- value_after("--organoids")
limit <- as.integer(value_after("--limit", "0"))
workers <- as.integer(value_after("--workers", "4"))
seed <- as.integer(value_after("--seed", "1"))

for (package in c("Matrix", "ape", "cliqueR", "jsonlite", "parallel")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
work_dir <- file.path(output_dir, "recordings")
log_dir <- file.path(output_dir, "logs")
for (directory in c(work_dir, log_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

manifest <- read.csv(file.path(ground_truth, "manifest.csv"),
                     stringsAsFactors = FALSE)
if (!is.null(requested)) {
  wanted <- trimws(strsplit(requested, ",", fixed = TRUE)[[1L]])
  manifest <- manifest[manifest$organoid_id %in% wanted, , drop = FALSE]
}
if (limit > 0L) manifest <- utils::head(manifest, limit)
if (!nrow(manifest)) stop("No organoids selected.", call. = FALSE)
cat(sprintf("[organoid] %d organoids, %d integrations, horizon %g minutes\n",
            nrow(manifest), integrations, end_time))

#' Replay the recorder over one organoid and score each of its clones
#'
#' @param organoid_id Directory name under the ground-truth root.
#' @return Data frame with one row per scored clone, or NULL if the organoid
#'   could not be processed.
process_organoid <- function(organoid_id) {
  source_dir <- file.path(ground_truth, organoid_id)
  cells <- read.csv(file.path(source_dir, "cells.csv"), stringsAsFactors = FALSE)
  clones <- read.csv(file.path(source_dir, "clones.csv"),
                     stringsAsFactors = FALSE)
  truth_trees <- ape::read.tree(file.path(source_dir, "truth_clones.nwk"))

  # Founders that never divided are absent from the division log, so the
  # reconstruction cannot see them and rejects their live cells. Supplying the
  # founder set explicitly is what --founders exists for.
  founders_path <- file.path(work_dir, sprintf("%s_founders.csv", organoid_id))
  utils::write.csv(data.frame(ID = sort(unique(cells$founder))),
                   founders_path, row.names = FALSE)

  recording_dir <- file.path(work_dir, organoid_id)
  matrix_path <- file.path(recording_dir,
                           "barcode_binary_score_matrix_sparse.rds")
  if (!file.exists(matrix_path)) {
    output <- suppressWarnings(system2(
      "Rscript",
      c(file.path(repo_root, "origin", "inst", "scripts",
                  "simulate_physicell_lineage.R"),
        paste0("--lineage=", file.path(source_dir, "cell_divisions.csv")),
        paste0("--live-cells=", file.path(source_dir, "live_cells.csv")),
        paste0("--founders=", founders_path),
        paste0("--params=", params_path),
        paste0("--end-time=", format(end_time, scientific = FALSE)),
        paste0("--num-integrations=", integrations),
        "--modalities=barcode", paste0("--seed=", seed),
        "--progress=false",
        paste0("--output-dir=", recording_dir)),
      stdout = TRUE, stderr = TRUE
    ))
    writeLines(output, file.path(log_dir, paste0(organoid_id, ".log")))
    if (!file.exists(matrix_path)) return(NULL)
  }

  score_matrix <- readRDS(matrix_path)
  # Recorder rows are cell_<ID>; the ground truth uses c<ID> throughout.
  rownames(score_matrix) <- sub("^cell_", "c", rownames(score_matrix))

  scorable <- clones[clones$tree_line >= 0 &
                       clones$n_alive_leaves >= min_leaves, , drop = FALSE]
  if (!nrow(scorable)) return(NULL)

  rows <- lapply(seq_len(nrow(scorable)), function(index) {
    founder <- scorable$founder[index]
    clone_cells <- cells$cell_id[cells$founder == founder]
    present <- rownames(score_matrix) %in% clone_cells
    if (sum(present) < min_leaves) return(NULL)
    clone_matrix <- as.matrix(score_matrix[present, , drop = FALSE])
    # Characters constant across the clone carry no signal for it, even when
    # they are informative elsewhere in the organoid.
    edited <- colSums(clone_matrix != 0)
    informative <- edited > 0 & edited < nrow(clone_matrix)
    truth_tree <- truth_trees[[scorable$tree_line[index] + 1L]]

    if (!any(informative)) {
      # Every character constant: the recorder captured nothing for this clone.
      return(data.frame(
        organoid_id = organoid_id, founder = founder,
        leaves = sum(present), informative_characters = 0L,
        edited_fraction = mean(clone_matrix != 0),
        normalized_rf = NA_real_, status = "no_informative_characters",
        stringsAsFactors = FALSE
      ))
    }
    clone_matrix <- clone_matrix[, informative, drop = FALSE]
    estimate <- tryCatch({
      built <- cliqueR::build_trees(clone_matrix, methods = "nj")
      built$nj
    }, error = function(condition) NULL)
    if (is.null(estimate)) {
      return(data.frame(
        organoid_id = organoid_id, founder = founder,
        leaves = sum(present), informative_characters = sum(informative),
        edited_fraction = mean(clone_matrix != 0),
        normalized_rf = NA_real_, status = "nj_failed",
        stringsAsFactors = FALSE
      ))
    }
    shared <- intersect(estimate$tip.label, truth_tree$tip.label)
    if (length(shared) < min_leaves) {
      return(data.frame(
        organoid_id = organoid_id, founder = founder,
        leaves = sum(present), informative_characters = sum(informative),
        edited_fraction = mean(clone_matrix != 0),
        normalized_rf = NA_real_, status = "tip_mismatch",
        stringsAsFactors = FALSE
      ))
    }
    distance <- cliqueR::rf_distance(
      ape::keep.tip(estimate, shared), ape::keep.tip(truth_tree, shared),
      normalize = TRUE
    )
    data.frame(
      organoid_id = organoid_id, founder = founder, leaves = length(shared),
      informative_characters = sum(informative),
      edited_fraction = mean(clone_matrix != 0),
      normalized_rf = distance, status = "ok", stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

results <- parallel::mclapply(manifest$organoid_id, function(organoid_id) {
  outcome <- tryCatch(process_organoid(organoid_id), error = function(condition) {
    message(sprintf("[%s] %s", organoid_id, conditionMessage(condition)))
    NULL
  })
  cat(sprintf("[organoid] %s: %s clones scored\n", organoid_id,
              if (is.null(outcome)) 0L else nrow(outcome)))
  outcome
}, mc.cores = workers, mc.preschedule = FALSE)

clone_accuracy <- do.call(rbind, Filter(Negate(is.null), results))
if (is.null(clone_accuracy)) stop("No clone was scored.", call. = FALSE)
clone_accuracy <- merge(
  clone_accuracy,
  manifest[, c("organoid_id", "skew", "replicate", "n_alive", "max_generation")],
  by = "organoid_id", all.x = TRUE
)

accuracy_path <- file.path(output_dir, "clone_accuracy.tsv.gz")
connection <- gzfile(accuracy_path, open = "wt")
write.table(clone_accuracy, connection, sep = "\t", row.names = FALSE,
            quote = FALSE)
close(connection)

scored <- clone_accuracy[clone_accuracy$status == "ok", , drop = FALSE]
organoid_summary <- do.call(rbind, lapply(
  split(scored, scored$organoid_id), function(part) {
    data.frame(
      organoid_id = part$organoid_id[1], skew = part$skew[1],
      clones_scored = nrow(part), mean_leaves = mean(part$leaves),
      mean_normalized_rf = mean(part$normalized_rf),
      median_normalized_rf = stats::median(part$normalized_rf),
      perfect_clones = sum(part$normalized_rf == 0),
      stringsAsFactors = FALSE
    )
  }
))
write.table(organoid_summary, file.path(output_dir, "organoid_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

jsonlite::write_json(list(
  ground_truth = ground_truth, params = params_path,
  integrations = integrations, end_time = end_time, min_leaves = min_leaves,
  seed = seed, organoids = nrow(manifest),
  recorder = "BASELINE base editing, barcode modality only"
), file.path(output_dir, "settings.json"), pretty = TRUE, auto_unbox = TRUE)

cat(sprintf("\n%d clones scored across %d organoids (%d excluded)\n",
            nrow(scored), length(unique(scored$organoid_id)),
            nrow(clone_accuracy) - nrow(scored)))
if (nrow(clone_accuracy) > nrow(scored)) {
  print(table(clone_accuracy$status[clone_accuracy$status != "ok"]))
}
cat(sprintf("mean normalized RF %.4f, median %.4f, perfect %d (%.1f%%)\n",
            mean(scored$normalized_rf), stats::median(scored$normalized_rf),
            sum(scored$normalized_rf == 0),
            100 * mean(scored$normalized_rf == 0)))
cat(sprintf("\nWrote:\n  %s\n  %s\n", accuracy_path,
            file.path(output_dir, "organoid_summary.tsv")))
