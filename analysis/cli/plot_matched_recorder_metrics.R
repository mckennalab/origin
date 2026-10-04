#!/usr/bin/env Rscript

# Multi-metric comparison of lineage recorders at matched capacity.
#
# Seven metrics across four recorders does not read as a table, and the metrics
# do not share a scale or a direction, so the figure normalises nothing and
# instead groups panels by what they measure. Two reference lines appear on
# every panel:
#
#   random    a random topology on the same tips: the floor
#   perfect   neighbour joining on the TRUE patristic distances: the ceiling the
#             reconstruction method itself imposes
#
# Without those, a number like "triplet agreement 0.78" is uninterpretable --
# chance is 0.33, so it is well above the floor, but the ceiling is 1.00.
#
# Panel D is the one that carries the mechanism rather than the ranking: clade
# recall split by when the split arose. Binary recorders lose ancient splits
# because an early edit is indistinguishable from a later independent edit at
# the same site, and deep clades depend on exactly those shared early edits.
#
# Usage:
#   Rscript plot_matched_recorder_metrics.R <run-dir> [output-prefix] [subtitle] [input-stem] [exclude]
#
# exclude is a comma-separated list matched case-insensitively against recorder
# labels, for dropping a recorder that does not belong in a given comparison.
#
# panels selects which groups to draw, from topology, ancestry, clone and era,
# defaulting to all four. The figure height scales with the count so a shorter
# selection does not come out stretched.
#
# input-stem names the metric files to read, defaulting to the matched-capacity
# runs. The parameter-grid assessment sweeps an editing-rate tier as well, and
# when a `tier` column is present each recorder is shown at its best tier so the
# panel compares recorders rather than rate settings.

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments)) {
  stop("Usage: plot_matched_recorder_metrics.R <run-dir> [prefix] [subtitle]",
       call. = FALSE)
}
run_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(run_dir, "matched_recorder_metrics")
}
subtitle_text <- if (length(arguments) >= 3L) arguments[[3L]] else {
  "Neutral trees, 250 sampled cells, 5 seeds, neighbour joining"
}

if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) {
  cache_dir <- file.path(tempdir(), "xdg-cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = cache_dir)
}
for (package in c("ggplot2", "ragg")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
suppressPackageStartupMessages(library(ggplot2))

input_stem <- if (length(arguments) >= 4L) arguments[[4L]] else {
  "matched_targets_tree_metrics"
}
summary_table <- read.csv(
  file.path(run_dir, paste0(input_stem, "_summary.csv")),
  stringsAsFactors = FALSE)
per_condition <- read.csv(
  file.path(run_dir, paste0(input_stem, "_per_condition.csv")),
  stringsAsFactors = FALSE)
era_path <- file.path(run_dir, paste0(input_stem, "_by_era.csv"))
era_table <- if (file.exists(era_path)) {
  read.csv(era_path, stringsAsFactors = FALSE)
} else NULL

# Drop excluded recorders before anything is ranked or scaled, so the ordering,
# palette and shared axes are all computed on what is actually shown.
excluded <- if (length(arguments) >= 5L && nzchar(arguments[[5L]])) {
  trimws(strsplit(arguments[[5L]], ",", fixed = TRUE)[[1L]])
} else character()
if (length(excluded)) {
  drop_rows <- function(frame) {
    if (is.null(frame) || !"label" %in% names(frame)) return(frame)
    keep <- !Reduce(`|`, lapply(excluded, function(pattern) {
      grepl(pattern, frame$label, ignore.case = TRUE)
    }))
    frame[keep, , drop = FALSE]
  }
  summary_table <- drop_rows(summary_table)
  per_condition <- drop_rows(per_condition)
  era_table <- drop_rows(era_table)
  cat(sprintf("[plot] excluded: %s\n", paste(excluded, collapse = ", ")))
}

observed <- summary_table[summary_table$variant == "observed", , drop = FALSE]
# Recorder order and colour are fixed by identity, not by how a given run ranks
# them. Ranking per figure would repaint the same recorder a different colour in
# every comparison, which makes the target-matched, information-matched and
# 5-integration figures impossible to read against one another.
canonical_order <- c("PEtracer (prime)", "PALINCODE", "BASELINE (cas12a)",
                     "WT-CRISPR (FLARE)", "Mitochondrial")
canonical_palette <- stats::setNames(
  c("#B4436C", "#087E8B", "#2E7096", "#C98A4B", "#6A6A8A"), canonical_order)

# A rate-tier sweep carries several rows per recorder; show each at its best
# tier, so the figure compares recorders rather than rate settings.
best_tier <- NULL
if ("tier" %in% names(observed)) {
  observed <- do.call(rbind, lapply(split(observed, observed$label),
                                    function(part) {
    part[which.min(part$unrooted_rf), , drop = FALSE]
  }))
  best_tier <- stats::setNames(observed$tier, observed$label)
  per_condition <- per_condition[
    paste(per_condition$label, per_condition$tier) %in%
      paste(observed$label, observed$tier), , drop = FALSE]
}
recorder_order <- canonical_order[canonical_order %in% observed$label]
unknown_labels <- setdiff(observed$label, canonical_order)
if (length(unknown_labels)) {
  # Appending rather than dropping: a new recorder should appear, just without
  # a reserved colour, instead of silently vanishing from the figure.
  recorder_order <- c(recorder_order, sort(unknown_labels))
}
apply_order <- function(data) {
  data$label <- factor(data$label, levels = recorder_order)
  data
}
observed <- apply_order(observed)
per_condition <- apply_order(
  per_condition[per_condition$variant == "observed", , drop = FALSE])
palette <- canonical_palette[recorder_order]
names(palette) <- recorder_order
if (anyNA(palette)) {
  spare <- c("#5C6B73", "#8A6FA8", "#7A8C4E")
  missing <- which(is.na(palette))
  palette[missing] <- spare[seq_along(missing)]
}

reference_of <- function(metric, variant) {
  part <- summary_table[summary_table$variant == variant, , drop = FALSE]
  if (!nrow(part) || !metric %in% names(part)) return(NA_real_)
  mean(part[[metric]], na.rm = TRUE)
}

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.4, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  axis.text.x = ggplot2::element_text(size = 8),
  legend.position = "none"
)

#' One grouped panel of metrics sharing a direction and scale
metric_panel <- function(metrics, labels, title, subtitle, direction,
                         ylim = NULL) {
  long <- do.call(rbind, lapply(seq_along(metrics), function(index) {
    metric <- metrics[index]
    frame <- per_condition[, c("label", metric)]
    names(frame) <- c("label", "value")
    frame$metric <- factor(labels[index], levels = labels)
    frame
  }))
  means <- do.call(rbind, lapply(seq_along(metrics), function(index) {
    metric <- metrics[index]
    data.frame(label = observed$label, value = observed[[metric]],
               metric = factor(labels[index], levels = labels),
               stringsAsFactors = FALSE)
  }))
  references <- do.call(rbind, lapply(seq_along(metrics), function(index) {
    data.frame(
      metric = factor(rep(labels[index], 2L), levels = labels),
      variant = c("random", "perfect"),
      value = c(reference_of(metrics[index], "random"),
                reference_of(metrics[index], "perfect")),
      stringsAsFactors = FALSE)
  }))
  references <- references[is.finite(references$value), , drop = FALSE]
  # Points alone show the spread but not whether two bars actually differ;
  # a 95% interval on the mean is what separates a real gap from seed noise.
  intervals <- do.call(rbind, lapply(
    split(long, list(long$metric, long$label), drop = TRUE),
    function(part) {
      values <- part$value[is.finite(part$value)]
      if (length(values) < 2L) return(NULL)
      standard_error <- stats::sd(values) / sqrt(length(values))
      data.frame(label = part$label[1], metric = part$metric[1],
                 lower = mean(values) - 1.96 * standard_error,
                 upper = mean(values) + 1.96 * standard_error,
                 stringsAsFactors = FALSE)
    }))

  plot <- ggplot2::ggplot(means, ggplot2::aes(x = label, y = value,
                                              fill = label)) +
    ggplot2::geom_hline(
      data = references,
      ggplot2::aes(yintercept = value, linetype = variant),
      colour = "#6A6A6A", linewidth = 0.45) +
    ggplot2::geom_col(width = 0.68, alpha = 0.9) +
    ggplot2::geom_errorbar(data = intervals,
                           ggplot2::aes(x = label, ymin = lower, ymax = upper),
                           inherit.aes = FALSE, width = 0.18, linewidth = 0.4,
                           colour = "#2A2A2A") +
    ggplot2::geom_point(data = long, ggplot2::aes(x = label, y = value),
                        inherit.aes = FALSE, size = 0.9, alpha = 0.55,
                        colour = "#2A2A2A",
                        position = ggplot2::position_jitter(width = 0.12,
                                                            height = 0)) +
    ggplot2::facet_wrap(~metric, nrow = 1) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_linetype_manual(
      values = c(random = "dotted", perfect = "dashed"), guide = "none") +
    ggplot2::scale_x_discrete(labels = function(x) sub(" \\(.*", "", x)) +
    ggplot2::labs(title = title,
                  subtitle = paste0(subtitle, "  (", direction, ")"),
                  x = NULL, y = NULL) +
    panel_theme
  # coord_cartesian rather than scale limits: it clips the view without
  # dropping rows, so bars and intervals stay computed on every seed.
  if (!is.null(ylim)) plot <- plot + ggplot2::coord_cartesian(ylim = ylim)
  plot
}

topology_panel <- metric_panel(
  c("unrooted_rf", "rooted_rf", "clade_recall"),
  c("unrooted RF", "rooted RF", "clade recall"),
  "Topology recovery",
  "dashed = perfect recorder, dotted = random tree",
  "RF lower is better, recall higher")

ancestry_panel <- metric_panel(
  c("triplet_agreement", "ordering_agreement", "cophenetic_spearman"),
  c("rooted triplets", "ancestry ordering", "cophenetic rho"),
  "Rooted ancestry and distance geometry",
  "chance is 0.33 for triplets, 0.50 for ordering, 0.00 for rho",
  "higher is better")

clone_panel <- metric_panel(
  c("clone_ari_k5", "clone_ari_k10", "clone_ari_k25"),
  c("k = 5", "k = 10", "k = 25"),
  "Clone assignment (adjusted Rand index)",
  "both trees cut into the same number of groups, so this scores assignment",
  "higher is better",
  # Fixed 0-1 so the three granularities, and the three capacity designs, are
  # read on one scale. ARI's chance level is 0 and its maximum is 1, so a
  # slightly negative seed is worse-than-chance noise clipped from view, not
  # excluded from the mean.
  ylim = c(0, 1))

era_panel <- if (!is.null(era_table)) {
  era_observed <- era_table[era_table$variant == "observed", , drop = FALSE]
  if (!is.null(best_tier) && "tier" %in% names(era_observed)) {
    era_observed <- era_observed[
      paste(era_observed$label, era_observed$tier) %in%
        paste(names(best_tier), best_tier), , drop = FALSE]
  }
  era_observed <- era_observed[era_observed$label %in% recorder_order, ,
                               drop = FALSE]
  era_observed <- apply_order(era_observed)
  era_observed$era <- factor(era_observed$era,
                             levels = c("earliest", "early", "late", "latest"))
  ggplot2::ggplot(era_observed,
                  ggplot2::aes(x = era, y = recall, colour = label,
                               group = label)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Which splits are recovered: ancient or recent",
      subtitle = paste(
        "True clades binned by when they arose (the truth tree is ultrametric,",
        "so node time is literally split time).\nA steep slope means a recorder",
        "resolves recent structure but loses deep structure"),
      x = "When the split arose", y = "Clade recall"
    ) +
    panel_theme +
    ggplot2::theme(legend.position = "right",
                   legend.title = ggplot2::element_blank())
} else NULL

available <- list(topology = topology_panel, ancestry = ancestry_panel,
                  clone = clone_panel, era = era_panel)
requested <- if (length(arguments) >= 6L && nzchar(arguments[[6L]])) {
  trimws(strsplit(arguments[[6L]], ",", fixed = TRUE)[[1L]])
} else names(available)
unknown <- setdiff(requested, names(available))
if (length(unknown)) {
  stop("Unknown panel(s): ", paste(unknown, collapse = ", "),
       ". Choose from topology, ancestry, clone, era.", call. = FALSE)
}
panels <- Filter(Negate(is.null), available[requested])
if (!length(panels)) stop("No panels selected.", call. = FALSE)
cat(sprintf("[plot] panels: %s\n", paste(names(panels), collapse = ", ")))
combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(panels, ncol = 1) +
    patchwork::plot_annotation(
      title = "Lineage recorder comparison at matched capacity",
      subtitle = subtitle_text,
      caption = paste(
        "Bars are means over 5 seeds; points are individual seeds. The perfect-recorder reference is neighbour joining run on the true patristic",
        "\ndistances, and it recovers the tree exactly on every topology metric -- so NJ imposes no ceiling here and all observed error is the",
        "\nrecorder's, not the method's. Clade precision equals recall by construction: both trees are fully resolved and binary, so they contain the",
        "\nsame number of clades and the two ratios share a denominator."
      ),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
        plot.caption = ggplot2::element_text(size = 7.6, colour = "#4A4A4A",
                                             hjust = 0),
        plot.tag = ggplot2::element_text(face = "bold", size = 13)
      )
    )
} else {
  topology_panel
}

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
figure_height <- 2.4 + 3.4 * length(panels)
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 12, height = figure_height, units = "in", dpi = 200,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 12, height = figure_height, units = "in",
                bg = "#FFFFFF")
cat(sprintf("Wrote:\n  %s\n  %s\n", png_path, pdf_path))
