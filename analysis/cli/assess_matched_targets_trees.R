#!/usr/bin/env Rscript

# Multi-metric assessment of the matched-target NJ reconstructions.
#
# Normalised RF answers one narrow question -- what fraction of unrooted splits
# match -- and treats every split as equally important. These recorders are
# compared on what a lineage tree is actually used for, so the assessment adds
# rooted structure, clade-level precision and recall, distance geometry,
# ancestry ordering, clone assignment, and whether the errors fall on early or
# late splits.
#
# Every metric is reported against two references, because a bare number is not
# interpretable on its own:
#
#   random    a random topology on the same tips. The floor. Any metric where a
#             recorder fails to beat this is carrying no usable signal.
#   perfect   neighbour joining run on the TRUE patristic distances. The
#             ceiling that NJ itself imposes given flawless input, so the gap
#             from `perfect` to a recorder is the recorder's cost, and the gap
#             from 0 to `perfect` is the reconstruction method's own.
#
# Rooting. NJ returns an unrooted tree while the truth is rooted at the founder,
# so estimates are midpoint-rooted. That is a choice, not a recovery of the true
# root, and it is applied identically to the observed, perfect and random trees
# so the comparison between them stays fair.
#
# Usage:
#   Rscript assess_matched_targets_trees.R <run-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments)) {
  stop("Usage: assess_matched_targets_trees.R <run-dir> [prefix]", call. = FALSE)
}
run_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(run_dir, "matched_targets_tree_metrics")
}
# A multi-shape run nests reconstructions one level deeper, under the shape.
# Single-shape runs keep the flat layout, so both are resolved rather than
# forcing the earlier runs to be reorganised.
shape <- if (length(arguments) >= 3L) arguments[[3L]] else "neutral"
triplet_samples <- 20000L
ordering_samples <- 20000L
clone_counts <- c(5L, 10L, 25L)
seed <- 1L

for (package in c("ape", "phangorn", "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
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

source(file.path(dirname(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1])), "..", "tree_metrics.R"))

# ---- run over every condition ----------------------------------------------

set.seed(seed)
rows <- list()
era_rows <- list()
for (index in seq_len(nrow(matched))) {
  system <- matched$system[index]
  integrations <- matched$integrations[index]
  nested <- file.path(run_dir, "reconstruction", shape,
                      sprintf("%s_k_%d", system, integrations), "conditions")
  flat <- file.path(run_dir, "reconstruction",
                    sprintf("%s_k_%d", system, integrations), "conditions")
  condition_root <- if (dir.exists(nested)) nested else flat
  if (!dir.exists(condition_root)) {
    warning("No reconstructions for ", system, call. = FALSE)
    next
  }
  conditions <- list.dirs(condition_root, recursive = FALSE, full.names = TRUE)
  for (condition in conditions) {
    estimate_path <- file.path(condition, "nj", "tree_nj.nwk")
    if (!file.exists(estimate_path)) next
    seed_label <- sub(".*_(seed_[0-9]+)_.*", "\\1", basename(condition))
    truth_path <- Sys.glob(file.path(
      run_dir, "simulation", paste0("shape_", shape), seed_label,
      paste0("recorder_", system), "parameter_*",
      sprintf("condition_level_%02d", integrations), "ground_truth_tree.nwk"))
    if (!length(truth_path)) {
      truth_path <- Sys.glob(file.path(
        run_dir, "simulation", paste0("shape_", shape), seed_label,
        paste0("recorder_", system), "parameter_*", "condition*",
        "ground_truth_tree.nwk"))
    }
    if (!length(truth_path)) {
      warning("No truth tree for ", basename(condition), call. = FALSE)
      next
    }
    # Division waypoints leave unary nodes, and a unary root breaks anything
    # that unroots; collapsing is topology-preserving.
    truth <- ape::collapse.singles(ape::read.tree(truth_path[[1L]]))
    estimate <- ape::read.tree(estimate_path)
    tips <- intersect(truth$tip.label, estimate$tip.label)
    if (length(tips) < 10L) next
    truth <- ape::keep.tip(truth, tips)
    estimate <- ape::keep.tip(estimate, tips)

    # The two references share the observed tree's tips and rooting treatment.
    perfect_tree <- ape::nj(stats::as.dist(
      ape::cophenetic.phylo(truth)[tips, tips]))
    random_tree <- ape::rtree(length(tips), tip.label = tips)

    variants <- list(
      observed = root_estimate(estimate),
      perfect = root_estimate(perfect_tree),
      random = random_tree
    )
    for (variant in names(variants)) {
      candidate <- variants[[variant]]
      if (is.null(candidate)) next
      scored <- score_pair(truth, candidate, tips)
      era <- scored$era
      scored$era <- NULL
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(system = system, label = matched$label[index],
                   integrations = integrations, targets = matched$targets[index],
                   seed = seed_label, variant = variant,
                   stringsAsFactors = FALSE),
        as.data.frame(scored, stringsAsFactors = FALSE)
      )
      if (!is.null(era)) {
        era$system <- system
        era$label <- matched$label[index]
        era$variant <- variant
        era$seed <- seed_label
        era_rows[[length(era_rows) + 1L]] <- era
      }
    }
    cat(sprintf("[assess] %s %s done\n", system, seed_label))
  }
}
scores <- do.call(rbind, rows)
if (is.null(scores)) stop("Nothing scored.", call. = FALSE)
era_scores <- do.call(rbind, era_rows)

metric_columns <- setdiff(names(scores), c("system", "label", "integrations",
                                           "targets", "seed", "variant"))
summary_table <- do.call(rbind, lapply(
  split(scores, list(scores$label, scores$variant), drop = TRUE),
  function(part) {
    values <- lapply(metric_columns, function(column) {
      mean(part[[column]], na.rm = TRUE)
    })
    names(values) <- metric_columns
    cbind(data.frame(label = part$label[1], variant = part$variant[1],
                     integrations = part$integrations[1],
                     targets = part$targets[1], replicates = nrow(part),
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
    split(era_scores, list(era_scores$label, era_scores$variant,
                           era_scores$era), drop = TRUE),
    function(part) data.frame(
      label = part$label[1], variant = part$variant[1], era = part$era[1],
      clades = nrow(part), recall = mean(part$recovered),
      stringsAsFactors = FALSE)
  ))
  utils::write.csv(era_summary, paste0(output_prefix, "_by_era.csv"),
                   row.names = FALSE)
}

show <- function(columns, title, digits = 3) {
  cat(sprintf("\n== %s ==\n\n", title))
  part <- summary_table[summary_table$variant == "observed", , drop = FALSE]
  part <- part[order(part$label), c("label", "targets", columns)]
  print(part, row.names = FALSE, digits = digits)
  cat("\nreferences:\n")
  reference <- summary_table[summary_table$variant != "observed", , drop = FALSE]
  reference <- reference[order(reference$variant, reference$label),
                         c("label", "variant", columns)]
  print(reference, row.names = FALSE, digits = digits)
}

show(c("unrooted_rf", "rooted_rf"), "Topology distance (lower is better)")
show(c("clade_precision", "clade_recall"), "Clade precision and recall")
show(c("triplet_agreement", "ordering_agreement", "cophenetic_spearman"),
     "Rooted ancestry and distance geometry (higher is better)")
show(sprintf("clone_ari_k%d", clone_counts), "Clone assignment (ARI)")
if (!is.null(era_summary)) {
  cat("\n== Clade recall by when the split arose ==\n\n")
  wide <- stats::reshape(
    era_summary[era_summary$variant == "observed", c("label", "era", "recall")],
    idvar = "label", timevar = "era", direction = "wide")
  print(wide, row.names = FALSE, digits = 3)
}
cat(sprintf("\nWrote:\n  %s\n  %s\n", paste0(output_prefix, "_summary.csv"),
            paste0(output_prefix, "_per_condition.csv")))
