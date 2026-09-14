#!/usr/bin/env Rscript

# Recorder accuracy under integration and cell dropout, for both the organoid
# lineages and the Gillespie integration series.
#
# Two axes of loss are swept, and they are not symmetric:
#
#   integration dropout  removes characters from cells that are still present.
#                        This destroys signal without changing what is being
#                        reconstructed, so its effect on accuracy is real and
#                        directly interpretable.
#   cell dropout         removes tips. The tree being scored gets smaller, and
#                        normalised RF is a fraction of the internal splits of
#                        whatever tree remains. A smaller tree has fewer splits
#                        to get wrong, so RF falls for reasons that have nothing
#                        to do with recovery.
#
# Reading the cell-dropout axis as "dropout improves accuracy" is therefore a
# mistake, and panel C exists to say so quantitatively: it re-cuts the organoid
# clones by the size they ended up at, and asks whether cell dropout still moves
# accuracy once size is held fixed. If the curves lie on top of each other, the
# apparent improvement was size all along.
#
# Usage:
#   Rscript plot_dropout_sweeps.R <results-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
results_dir <- normalizePath(
  if (length(arguments)) arguments[[1L]] else ".", mustWork = TRUE
)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(results_dir, "dropout_sweeps")
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

organoid <- read.delim(file.path(results_dir, "organoid_dropout_summary.tsv"),
                       stringsAsFactors = FALSE)
gillespie <- read.delim(file.path(results_dir, "gillespie_dropout_summary.tsv"),
                        stringsAsFactors = FALSE)
clone_scores <- read.delim(
  gzfile(file.path(results_dir, "organoid_dropout_clone_scores.tsv.gz")),
  stringsAsFactors = FALSE
)
cat(sprintf("[plot] organoid %d conditions, gillespie %d conditions, %d clone scores\n",
            nrow(organoid), nrow(gillespie), nrow(clone_scores)))

dropout_levels <- sort(unique(organoid$cell_dropout))
dropout_palette <- stats::setNames(
  c("#1F3B4D", "#2E7096", "#5AA2B8", "#C98A4B", "#B4436C")[seq_along(dropout_levels)],
  format(dropout_levels)
)
as_percent <- function(x) sprintf("%g%%", 100 * x)

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.6, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  legend.key.height = grid::unit(0.8, "lines")
)

#' Dropout-grid panel: accuracy against integration dropout, one line per
#' cell-dropout level.
grid_panel <- function(data, title, subtitle) {
  data$cell_label <- factor(format(data$cell_dropout),
                            levels = format(dropout_levels))
  ggplot2::ggplot(data, ggplot2::aes(x = integration_dropout,
                                     y = mean_normalized_rf,
                                     colour = cell_label, fill = cell_label)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.9) +
    ggplot2::scale_colour_manual(name = "cell dropout", values = dropout_palette,
                                 labels = as_percent(dropout_levels)) +
    ggplot2::scale_fill_manual(name = "cell dropout", values = dropout_palette,
                               labels = as_percent(dropout_levels)) +
    ggplot2::scale_x_continuous(breaks = sort(unique(data$integration_dropout)),
                                labels = as_percent) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = "Integration dropout",
                  y = "Mean normalised RF (0 = exact)") +
    panel_theme
}

organoid_panel <- grid_panel(
  organoid,
  "Organoid lineages, BASELINE recorder",
  sprintf(paste("5 integrations, mean clone %.0f cells; %s-%s clones per condition",
                "(fewer survive at high dropout)"),
          organoid$mean_leaves[organoid$cell_dropout == 0 &
                                 organoid$integration_dropout == 0],
          format(min(organoid$clones), big.mark = ","),
          format(max(organoid$clones), big.mark = ","))
)

# A and B are not on the same scale, and the reason is the recorder size rather
# than the engine: the series condition used here carries 100 integrations
# against the organoid runs' 5, so it starts with far more signal per cell.
gillespie_baseline <- gillespie[gillespie$system == "baseline", , drop = FALSE]
gillespie_panel <- grid_panel(
  gillespie_baseline,
  "Gillespie series, BASELINE recorder",
  sprintf(paste("100 integrations, 250 sampled cells, %d replicates.",
                "Note the y-axis: 20x the recorder of panel A"),
          max(gillespie_baseline$replicates))
)

# Panel C: hold tree size fixed and ask whether cell dropout still matters.
size_breaks <- c(3, 5, 10, 20, 50, 100, Inf)
size_labels <- c("4-5", "6-10", "11-20", "21-50", "51-100", ">100")
matched <- clone_scores[clone_scores$integration_dropout == 0, , drop = FALSE]
matched$size_bin <- cut(matched$leaves, breaks = size_breaks, labels = size_labels)
matched <- matched[!is.na(matched$size_bin), , drop = FALSE]
matched_summary <- do.call(rbind, lapply(
  split(matched, list(matched$size_bin, matched$cell_dropout), drop = TRUE),
  function(part) {
    standard_error <- stats::sd(part$normalized_rf) / sqrt(nrow(part))
    data.frame(
      size_bin = as.character(part$size_bin[1]),
      cell_dropout = part$cell_dropout[1], clones = nrow(part),
      mean_normalized_rf = mean(part$normalized_rf),
      lower = mean(part$normalized_rf) - 1.96 * standard_error,
      upper = mean(part$normalized_rf) + 1.96 * standard_error,
      stringsAsFactors = FALSE
    )
  }
))
matched_summary$size_bin <- factor(matched_summary$size_bin, levels = size_labels)
matched_summary$cell_label <- factor(format(matched_summary$cell_dropout),
                                     levels = format(dropout_levels))

# How far apart the curves actually sit, quoted so the claim is not left to the
# eye. Split at 80%, because that level behaves differently from the rest: it is
# the only one where matched-size accuracy gets worse rather than staying flat.
curve_spread <- function(data) {
  values <- vapply(split(data, data$size_bin), function(part) {
    if (nrow(part) < 2L) return(NA_real_)
    diff(range(part$mean_normalized_rf))
  }, numeric(1))
  max(values, na.rm = TRUE)
}
moderate_spread <- curve_spread(
  matched_summary[matched_summary$cell_dropout <= 0.5, , drop = FALSE]
)
# Worst matched-size penalty at 80%, against full capture in the same size bin.
paired <- merge(
  matched_summary[matched_summary$cell_dropout == 0.8,
                  c("size_bin", "mean_normalized_rf")],
  matched_summary[matched_summary$cell_dropout == 0,
                  c("size_bin", "mean_normalized_rf")],
  by = "size_bin", suffixes = c("_high", "_full")
)
extreme_penalty <- max(paired$mean_normalized_rf_high -
                         paired$mean_normalized_rf_full)

matched_panel <- ggplot2::ggplot(
  matched_summary, ggplot2::aes(x = size_bin, y = mean_normalized_rf,
                                colour = cell_label, group = cell_label)
) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), width = 0.14,
                         linewidth = 0.4, alpha = 0.8) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.9) +
  ggplot2::scale_colour_manual(name = "cell dropout", values = dropout_palette,
                               labels = as_percent(dropout_levels)) +
  ggplot2::labs(
    title = "Cell dropout, with tree size held fixed",
    subtitle = sprintf(
      paste("Organoid clones at zero integration dropout, re-cut by surviving size.",
            "Up to 50%% the curves sit within %.3f RF;\nat 80%% accuracy is worse at matched size, by up to %.3f RF"),
      moderate_spread, extreme_penalty
    ),
    x = "Cells remaining in clone", y = "Mean normalised RF"
  ) +
  panel_theme

# Panel D: the same integration-dropout axis across recorder systems, with no
# cell dropout, so the comparison is between chemistries rather than sizes.
by_system <- gillespie[gillespie$cell_dropout == 0, , drop = FALSE]
system_palette <- c(baseline = "#B4436C", wt_crispr = "#1F3B4D",
                    prime = "#2E7096", palincode = "#087E8B",
                    mitochondrial = "#C98A4B")
system_panel <- ggplot2::ggplot(
  by_system[by_system$system != "mitochondrial", , drop = FALSE],
  ggplot2::aes(x = integration_dropout, y = mean_normalized_rf,
               colour = system, fill = system)
) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.13,
                       colour = NA) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.9) +
  # Mitochondrial characters are not grouped into integrations, so it has no
  # integration-dropout axis to sweep; it is shown as its single measured point.
  ggplot2::geom_point(
    data = by_system[by_system$system == "mitochondrial", , drop = FALSE],
    size = 2.6, shape = 17
  ) +
  ggplot2::scale_colour_manual(name = NULL, values = system_palette) +
  ggplot2::scale_fill_manual(name = NULL, values = system_palette) +
  ggplot2::scale_x_continuous(breaks = sort(unique(by_system$integration_dropout)),
                              labels = as_percent) +
  ggplot2::labs(
    title = "Integration dropout across recorder systems",
    subtitle = paste("Gillespie series at full cell capture, 100 integrations.",
                     "BASELINE leads throughout;\nmitochondrial (triangle) has no integrations to drop"),
    x = "Integration dropout", y = "Mean normalised RF"
  ) +
  panel_theme

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(organoid_panel, gillespie_panel, matched_panel,
                        system_panel, ncol = 2) +
    patchwork::plot_annotation(
      title = "Reconstruction accuracy under integration and cell dropout",
      subtitle = paste(
        "Integration dropout makes characters missing in cells that remain;",
        "cell dropout removes tips and the truth tree is pruned to match."
      ),
      caption = paste(
        "Falling RF along the cell-dropout axis in A and B is mostly a size effect, not a recovery effect: normalised RF is a fraction of the",
        "\ninternal splits of the surviving tree, and dropping cells leaves fewer splits to get wrong -- mean organoid clone size falls from 32.4",
        "\ncells to 10.0 across that axis. Panel C removes it by re-cutting on surviving size, and the improvement largely disappears; at 80%",
        "\ndropout accuracy is worse at matched size. Note that matched surviving size is not a matched clone: a 21-50 cell tree at 80% dropout",
        "\ncame from a clone of roughly 105-250 cells, so those bins are drawn from larger, older clones. Integration dropout, which removes",
        "\ncharacters without removing tips, is the axis that costs accuracy unambiguously."
      ),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
        plot.caption = ggplot2::element_text(size = 7.8, colour = "#4A4A4A",
                                             hjust = 0),
        plot.tag = ggplot2::element_text(face = "bold", size = 13)
      )
    )
} else {
  organoid_panel
}

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 13, height = 9.5, units = "in", dpi = 220,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 13, height = 9.5, units = "in", bg = "#FFFFFF")
utils::write.csv(matched_summary, paste0(output_prefix, "_size_matched.csv"),
                 row.names = FALSE)

cat("\nOrganoid, mean normalised RF:\n")
print(stats::xtabs(mean_normalized_rf ~ integration_dropout + cell_dropout,
                   data = organoid), digits = 3)
cat("\nGillespie BASELINE, mean normalised RF:\n")
print(stats::xtabs(mean_normalized_rf ~ integration_dropout + cell_dropout,
                   data = gillespie_baseline), digits = 3)
cat(sprintf("\nMean cells per tree, organoid: %.1f at 0%% cell dropout -> %.1f at 80%%\n",
            organoid$mean_leaves[organoid$cell_dropout == 0 &
                                   organoid$integration_dropout == 0],
            organoid$mean_leaves[organoid$cell_dropout == 0.8 &
                                   organoid$integration_dropout == 0]))
cat(sprintf("At fixed surviving size: curves within %.4f RF up to 50%% dropout; 80%% is worse by up to %.4f RF\n",
            moderate_spread, extreme_penalty))
cat(sprintf("\nWrote:\n  %s\n  %s\n", png_path, pdf_path))
