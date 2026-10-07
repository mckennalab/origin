#!/usr/bin/env Rscript

# Recorder comparison across tree shapes at capacity-controlled integration
# counts, neighbour joining.
#
# One figure rather than four, because the interesting result is that the
# recorders do NOT rank the same way on every topology, and a reader can only
# see that when the shapes sit side by side.
#
# Two things the panels are built to expose:
#
#   comb disagrees with itself   no recorder recovers comb clades (RF ~0.99,
#                                recall ~0.01) yet the distance geometry is
#                                captured well (rho 0.85-0.95). A single RF
#                                number would report total failure and hide that
#                                relative distances are largely right.
#   hierarchical does not separate  BASELINE has the lowest mean RF, but paired
#                                tests over the shared seeds put PEtracer and
#                                PALINCODE within noise of it (p = 0.45, 0.53).
#                                Ranking bare means there would report a winner
#                                the replicates do not support.
#
# Capacity is controlled, not matched: BASELINE carries 17-58% more bits than
# PEtracer in every shape, so PEtracer's wins come from a capacity deficit and
# BASELINE's lower hierarchical mean is not a significant win, and it comes
# with a capacity advantage besides.
#
# Usage:
#   Rscript plot_allshapes_recorder_metrics.R <run-dir> [output-prefix] [clone-k]
#
# clone-k selects which clone-assignment granularity to show (5, 10 or 25).
# Only one is drawn so the figure stays readable across four shapes; the other
# granularities are in the summary CSVs.

arguments <- commandArgs(trailingOnly = TRUE)
run_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(run_dir, "allshapes_recorder_metrics")
}
clone_k <- if (length(arguments) >= 3L) as.integer(arguments[[3L]]) else 5L
if (!clone_k %in% c(5L, 10L, 25L)) {
  stop("clone-k must be 5, 10 or 25.", call. = FALSE)
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

shapes <- c("balanced", "comb", "neutral", "hierarchical")
shape_labels <- c(balanced = "balanced (unit-length)", comb = "comb",
                  neutral = "neutral", hierarchical = "hierarchical")
load_shape <- function(shape, stem) {
  path <- file.path(run_dir, sprintf("metrics_%s_%s.csv", shape, stem))
  if (!file.exists(path)) return(NULL)
  frame <- read.csv(path, stringsAsFactors = FALSE)
  frame$shape <- shape
  frame
}
summary_table <- do.call(rbind, lapply(shapes, load_shape, stem = "summary"))
per_condition <- do.call(rbind, lapply(shapes, load_shape,
                                       stem = "per_condition"))
era_table <- do.call(rbind, lapply(shapes, load_shape, stem = "by_era"))
if (is.null(summary_table)) stop("No metric summaries found.", call. = FALSE)

recorder_order <- c("PEtracer (prime)", "PALINCODE", "BASELINE (cas12a)",
                    "WT-CRISPR (FLARE)")
palette <- stats::setNames(c("#B4436C", "#087E8B", "#2E7096", "#C98A4B"),
                           recorder_order)
prepare <- function(frame) {
  frame <- frame[frame$label %in% recorder_order, , drop = FALSE]
  frame$label <- factor(frame$label, levels = recorder_order)
  frame$shape <- factor(shape_labels[frame$shape],
                        levels = unname(shape_labels))
  frame
}
observed <- prepare(summary_table[summary_table$variant == "observed", ,
                                  drop = FALSE])
seeds <- prepare(per_condition[per_condition$variant == "observed", ,
                               drop = FALSE])
cat(sprintf("[plot] %d shape-recorder cells, %d seed replicates\n",
            nrow(observed), nrow(seeds)))

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.4, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  axis.text.x = ggplot2::element_blank(),
  axis.ticks.x = ggplot2::element_blank(),
  legend.position = "none"
)

# References are shared across recorders within a shape, so one line per shape.
reference_for <- function(metric, variant) {
  part <- prepare(summary_table[summary_table$variant == variant, ,
                                drop = FALSE])
  if (!nrow(part)) return(NULL)
  stats::aggregate(part[[metric]], by = list(shape = part$shape),
                   FUN = function(x) mean(x, na.rm = TRUE)) |>
    stats::setNames(c("shape", "value"))
}

shape_panel <- function(metric, title, subtitle) {
  frame <- observed[, c("shape", "label", metric)]
  names(frame) <- c("shape", "label", "value")
  # Bars alone hid that several differences sit inside seed-to-seed noise, so
  # every panel carries the five replicates and a 95% interval on the mean.
  points <- seeds[, c("shape", "label", metric)]
  names(points) <- c("shape", "label", "value")
  points <- points[is.finite(points$value), , drop = FALSE]
  intervals <- do.call(rbind, lapply(
    split(points, list(points$shape, points$label), drop = TRUE),
    function(part) {
      standard_error <- stats::sd(part$value) / sqrt(nrow(part))
      data.frame(shape = part$shape[1], label = part$label[1],
                 lower = mean(part$value) - 1.96 * standard_error,
                 upper = mean(part$value) + 1.96 * standard_error,
                 stringsAsFactors = FALSE)
    }))
  references <- do.call(rbind, lapply(c("random", "perfect"), function(v) {
    reference <- reference_for(metric, v)
    if (is.null(reference)) return(NULL)
    reference$variant <- v
    reference
  }))
  plot <- ggplot2::ggplot(frame, ggplot2::aes(x = label, y = value,
                                              fill = label))
  if (!is.null(references)) {
    plot <- plot + ggplot2::geom_hline(
      data = references,
      ggplot2::aes(yintercept = value, linetype = variant),
      colour = "#6A6A6A", linewidth = 0.4)
  }
  plot +
    ggplot2::geom_col(width = 0.72, alpha = 0.9) +
    ggplot2::geom_errorbar(data = intervals,
                           ggplot2::aes(x = label, ymin = lower, ymax = upper),
                           inherit.aes = FALSE, width = 0.2, linewidth = 0.4,
                           colour = "#2A2A2A") +
    ggplot2::geom_point(data = points, ggplot2::aes(x = label, y = value),
                        inherit.aes = FALSE, size = 0.8, alpha = 0.6,
                        colour = "#1A1A1A",
                        position = ggplot2::position_jitter(width = 0.13,
                                                            height = 0)) +
    ggplot2::facet_wrap(~shape, nrow = 1) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_linetype_manual(
      values = c(random = "dotted", perfect = "dashed"), guide = "none") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    panel_theme
}

panels <- list(
  shape_panel("unrooted_rf", "Topology error (normalised RF)",
              "lower is better; dashed = perfect recorder, dotted = random"),
  shape_panel("clade_recall", "Clade recall", "higher is better"),
  shape_panel("triplet_agreement", "Rooted triplet agreement",
              "higher is better; chance is 0.33"),
  shape_panel(sprintf("clone_ari_k%d", clone_k),
              sprintf("Clone assignment, k = %d (ARI)", clone_k),
              "higher is better"),
  shape_panel("cophenetic_spearman", "Cophenetic rank correlation",
              "higher is better; the perfect reference sits at 0.93, not 1.00")
)

if (!is.null(era_table)) {
  era_observed <- prepare(era_table[era_table$variant == "observed", ,
                                    drop = FALSE])
  era_observed$era <- factor(era_observed$era,
                             levels = c("earliest", "early", "late", "latest"))
  panels[[length(panels) + 1L]] <- ggplot2::ggplot(
    era_observed, ggplot2::aes(x = era, y = recall, colour = label,
                               group = label)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~shape, nrow = 1) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Which splits are recovered: ancient or recent",
      subtitle = "True clades binned by when they arose; a flat line means deep structure survives",
      x = NULL, y = "Clade recall") +
    panel_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7.5, angle = 30,
                                                       hjust = 1),
                   axis.ticks.x = ggplot2::element_line(),
                   legend.position = "right",
                   legend.title = ggplot2::element_blank())
}

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(panels, ncol = 1) +
    patchwork::plot_annotation(
      title = "Recorder comparison across tree shapes, capacity-controlled",
      subtitle = paste(
        "Neighbour joining, 250 cells, 5 seeds. BASELINE 2x50, WT-CRISPR 13x10,",
        "PEtracer 14x3 (9-state), PALINCODE 30x2"),
      caption = paste(
        "Capacity is controlled, not matched: BASELINE carries 17-58% more bits than PEtracer in every shape, so PEtracer's wins come from a",
        "\ncapacity deficit. Comb is the case where a single metric misleads -- no recorder recovers",
        "\nits clades (RF ~0.99, recall ~0.01) yet cophenetic correlation reaches 0.85-0.95, so relative distances are largely right even though the",
        "\ntopology is not. Error bars are 95% intervals over 5 seeds and points are the seeds themselves: on hierarchical BASELINE, PALINCODE and",
        "\nPEtracer are statistically indistinguishable (paired p = 0.45 and 0.53), and on comb all four are tied at failure, so only balanced and",
        "\nneutral separate the recorders -- there PEtracer beats every other recorder on paired tests (p = 0.001-0.046)."
      ),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
        plot.caption = ggplot2::element_text(size = 7.6, colour = "#4A4A4A",
                                             hjust = 0),
        plot.tag = ggplot2::element_text(face = "bold", size = 13)))
} else {
  panels[[1L]]
}

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 13, height = 16, units = "in", dpi = 190,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 13, height = 16, units = "in", bg = "#FFFFFF")
cat(sprintf("Wrote:\n  %s\n  %s\n", png_path, pdf_path))
