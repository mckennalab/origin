#!/usr/bin/env Rscript

# Multi-metric assessment of the 5-integration recorder parameter grid.
#
# The grid holds integrations fixed at five and sweeps each recorder's editing
# rate over five tiers. That makes it the matched-INTEGRATION comparison, as
# opposed to the matched-target and matched-information runs, and the same
# metric set is applied here so the three designs can be read against each
# other.
#
# Every recorder-tier is scored against two references:
#
#   random    a random topology on the same tips: the floor
#   perfect   neighbour joining on the TRUE patristic distances: the ceiling the
#             reconstruction method itself imposes
#
# Both references depend only on the truth tree, and within a shape and seed
# every recorder is sampled from the same 250 cells and scored against the same
# truth tree. They are therefore computed once per seed and reused, rather than
# recomputed for all 25 recorder-tier combinations.
#
# Restricted to one tree shape (neutral by default): the shapes differ in depth
# and in how reconstructable they are, so pooling them would average a recorder
# that fails only on comb into one that fails everywhere.
#
# Usage:
#   Rscript assess_parameter_grid_trees.R <grid-dir> [shape] [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments)) {
  stop("Usage: assess_parameter_grid_trees.R <grid-dir> [shape] [prefix]",
       call. = FALSE)
}
grid_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
shape <- if (length(arguments) >= 2L) arguments[[2L]] else "neutral"
output_prefix <- if (length(arguments) >= 3L) arguments[[3L]] else {
  file.path(grid_dir, sprintf("parameter_grid_tree_metrics_%s", shape))
}
triplet_samples <- 20000L
ordering_samples <- 20000L
clone_counts <- c(5L, 10L, 25L)

for (package in c("ape", "phangorn", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
script_dir <- dirname(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(script_dir, "..", "tree_metrics.R"))

labels <- c(baseline = "BASELINE (cas12a)", wt_crispr = "WT-CRISPR (FLARE)",
            prime = "PEtracer (prime)", palincode = "PALINCODE",
            mitochondrial = "Mitochondrial")
tier_levels <- c("very_low", "low", "mid", "high", "very_high")

conditions <- list.dirs(file.path(grid_dir, "reconstruction", "conditions"),
                        recursive = FALSE, full.names = TRUE)
conditions <- conditions[startsWith(basename(conditions), paste0(shape, "_"))]
if (!length(conditions)) {
  stop("No ", shape, " conditions under ", grid_dir, call. = FALSE)
}
cat(sprintf("[assess] %d %s conditions\n", length(conditions), shape))

set.seed(1L)
reference_cache <- new.env(parent = emptyenv())
rows <- list()
era_rows <- list()

for (condition in conditions) {
  name <- basename(condition)
  estimate_path <- file.path(condition, "nj", "tree_nj.nwk")
  if (!file.exists(estimate_path)) next
  seed_label <- sub(".*_(seed_[0-9]+)_.*", "\\1", name)
  system <- sub(sprintf("^%s_seed_[0-9]+_(.*)_param_.*", shape), "\\1", name)
  tier <- sub(".*_param_(.*)_n_.*", "\\1", name)
  if (!system %in% names(labels)) next

  truth_path <- Sys.glob(file.path(
    grid_dir, "simulation", paste0("shape_", shape), seed_label,
    paste0("recorder_", system), paste0("parameter_", tier), "condition*",
    "ground_truth_tree.nwk"))
  if (!length(truth_path)) next
  # Unary division waypoints break anything that unroots; collapsing is
  # topology-preserving.
  truth <- ape::collapse.singles(ape::read.tree(truth_path[[1L]]))
  estimate <- ape::read.tree(estimate_path)
  tips <- intersect(truth$tip.label, estimate$tip.label)
  if (length(tips) < 10L) next
  truth <- ape::keep.tip(truth, tips)
  estimate <- ape::keep.tip(estimate, tips)

  emit <- function(variant, candidate, tier_name) {
    if (is.null(candidate)) return(invisible(NULL))
    scored <- score_pair(truth, candidate, tips, triplet_samples,
                         ordering_samples, clone_counts)
    era <- scored$era
    scored$era <- NULL
    rows[[length(rows) + 1L]] <<- cbind(
      data.frame(system = system, label = unname(labels[[system]]),
                 tier = tier_name, seed = seed_label, variant = variant,
                 stringsAsFactors = FALSE),
      as.data.frame(scored, stringsAsFactors = FALSE))
    if (!is.null(era)) {
      era$system <- system
      era$label <- unname(labels[[system]])
      era$tier <- tier_name
      era$variant <- variant
      era$seed <- seed_label
      era_rows[[length(era_rows) + 1L]] <<- era
    }
    invisible(NULL)
  }

  emit("observed", root_estimate(estimate), tier)

  # References depend only on the truth tree, which every recorder in this seed
  # shares, so they are scored once and carried under tier "reference".
  if (!exists(seed_label, envir = reference_cache, inherits = FALSE)) {
    perfect_tree <- ape::nj(stats::as.dist(
      ape::cophenetic.phylo(truth)[tips, tips]))
    emit("perfect", root_estimate(perfect_tree), "reference")
    emit("random", ape::rtree(length(tips), tip.label = tips), "reference")
    assign(seed_label, TRUE, envir = reference_cache)
  }
  cat(sprintf("[assess] %s %s %s\n", system, tier, seed_label))
}

scores <- do.call(rbind, rows)
if (is.null(scores)) stop("Nothing scored.", call. = FALSE)
era_scores <- do.call(rbind, era_rows)

metric_columns <- setdiff(names(scores),
                          c("system", "label", "tier", "seed", "variant"))
summary_table <- do.call(rbind, lapply(
  split(scores, list(scores$label, scores$tier, scores$variant), drop = TRUE),
  function(part) {
    values <- lapply(metric_columns, function(column) {
      mean(part[[column]], na.rm = TRUE)
    })
    names(values) <- metric_columns
    cbind(data.frame(label = part$label[1], tier = part$tier[1],
                     variant = part$variant[1], replicates = nrow(part),
                     stringsAsFactors = FALSE),
          as.data.frame(values))
  }
))
utils::write.csv(scores, paste0(output_prefix, "_per_condition.csv"),
                 row.names = FALSE)
utils::write.csv(summary_table, paste0(output_prefix, "_summary.csv"),
                 row.names = FALSE)

era_summary <- NULL
if (!is.null(era_scores)) {
  era_summary <- do.call(rbind, lapply(
    split(era_scores, list(era_scores$label, era_scores$tier,
                           era_scores$variant, era_scores$era), drop = TRUE),
    function(part) data.frame(
      label = part$label[1], tier = part$tier[1], variant = part$variant[1],
      era = part$era[1], clades = nrow(part), recall = mean(part$recovered),
      stringsAsFactors = FALSE)
  ))
  utils::write.csv(era_summary, paste0(output_prefix, "_by_era.csv"),
                   row.names = FALSE)
}

observed <- summary_table[summary_table$variant == "observed", , drop = FALSE]
best <- do.call(rbind, lapply(split(observed, observed$label), function(part) {
  part[which.min(part$unrooted_rf), , drop = FALSE]
}))
best <- best[order(best$unrooted_rf), ]

cat(sprintf("\n== %s trees, 5 integrations: each recorder at its best tier ==\n\n",
            shape))
print(best[, c("label", "tier", "unrooted_rf", "rooted_rf", "clade_recall",
               "triplet_agreement", "ordering_agreement",
               "cophenetic_spearman")],
      row.names = FALSE, digits = 3)
cat("\nclone assignment (ARI):\n\n")
print(best[, c("label", "tier", sprintf("clone_ari_k%d", clone_counts))],
      row.names = FALSE, digits = 3)

reference <- summary_table[summary_table$variant != "observed", , drop = FALSE]
reference <- do.call(rbind, lapply(split(reference, reference$variant),
                                   function(part) {
  values <- lapply(metric_columns, function(column) mean(part[[column]],
                                                         na.rm = TRUE))
  names(values) <- metric_columns
  cbind(data.frame(variant = part$variant[1], stringsAsFactors = FALSE),
        as.data.frame(values))
}))
cat("\nreferences (shared across recorders, one truth tree per seed):\n\n")
print(reference[, c("variant", "unrooted_rf", "clade_recall",
                    "triplet_agreement", "ordering_agreement",
                    "cophenetic_spearman", "clone_ari_k10")],
      row.names = FALSE, digits = 3)

if (!is.null(era_summary)) {
  cat("\n== Clade recall by when the split arose, at each best tier ==\n\n")
  keep <- merge(era_summary[era_summary$variant == "observed", ],
                best[, c("label", "tier")], by = c("label", "tier"))
  wide <- stats::reshape(keep[, c("label", "era", "recall")],
                         idvar = "label", timevar = "era", direction = "wide")
  print(wide, row.names = FALSE, digits = 3)
}
cat(sprintf("\nWrote:\n  %s\n  %s\n", paste0(output_prefix, "_summary.csv"),
            paste0(output_prefix, "_per_condition.csv")))
