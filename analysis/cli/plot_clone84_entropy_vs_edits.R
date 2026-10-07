#!/usr/bin/env Rscript

# Per-integration entropy against edit burden, sized by how much of the barcode
# was actually captured.
#
# Entropy and edit count are not the same axis and can come apart: a barcode can
# carry many edits yet few distinct patterns if the same edits recur across
# cells, and a lightly edited barcode can still be highly diverse if its few
# edits fall in different places. Plotting them together shows which regime each
# integration is in.
#
# Point area is the number of non-missing cell-by-site entries for that barcode,
# because both axes are sensitive to capture: a poorly captured barcode has
# fewer observed edits AND fewer distinguishable patterns, so a low-left point
# may be an assay artefact rather than a property of the recorder.
#
# Usage:
#   Rscript plot_clone84_entropy_vs_edits.R <results-dir> [output-prefix] [cells]

arguments <- commandArgs(trailingOnly = TRUE)
results_dir <- normalizePath(
  if (length(arguments)) arguments[[1L]] else ".", mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(results_dir, "clone84_entropy_vs_edits")
}
cells <- if (length(arguments) >= 3L) as.integer(arguments[[3L]]) else 250L

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

display <- c(real = "clone 84 (real)", real_all = "clone 84 (all 34)",
             shipped_shape0.5 = "simulated, shape 0.5",
             fitted_shape = "simulated, shape fitted",
             empirical_rates = "simulated, empirical rates")
palette <- stats::setNames(c("#1F3B4D", "#1F3B4D", "#C98A4B", "#B4436C",
                             "#3E8E7E"), unname(display))
data <- read.csv(file.path(results_dir, "feature_by_integration.csv"),
                 stringsAsFactors = FALSE)
data$label <- factor(display[data$dataset], levels = unname(display))
data <- data[!is.na(data$label), , drop = FALSE]
# Captured entries: observed fraction of the cell-by-site block for the barcode.
per_row_cells <- if ("cells" %in% names(data)) data$cells else cells
data$captured <- data$observed_fraction * data$sites * per_row_cells
cat(sprintf("[plot] %d integrations x %d datasets; captured %0.0f-%0.0f entries\n",
            length(unique(data$integration_index)), length(unique(data$label)),
            min(data$captured), max(data$captured)))

# Labelled only where a point is far enough from the pack to be worth naming,
# so the panel does not become a list of integration numbers.
outliers <- data[data$dataset == "real" &
                   (data$normalised_entropy <
                      stats::quantile(data$normalised_entropy[
                        data$dataset == "real"], 0.15) |
                      data$edits_per_cell_mean <
                      stats::quantile(data$edits_per_cell_mean[
                        data$dataset == "real"], 0.15)), , drop = FALSE]

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
plot <- ggplot2::ggplot(
  data, ggplot2::aes(x = normalised_entropy, y = edits_per_cell_mean,
                     colour = label, size = captured)
) +
  ggplot2::geom_point(alpha = 0.75) +
  ggplot2::geom_text(
    data = outliers, inherit.aes = FALSE,
    mapping = ggplot2::aes(x = normalised_entropy, y = edits_per_cell_mean,
                           label = sprintf("int %d", integration_index)),
    hjust = -0.3, vjust = 0.4, size = 3, colour = "#2A2A2A") +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::scale_size_area(
    max_size = 11,
    labels = function(x) format(round(x), big.mark = ","),
    name = "captured entries") +
  ggplot2::labs(
    title = "Allele entropy against edit burden, per integration",
    subtitle = paste(
      "One point per integration. Point area is the number of non-missing",
      "cell-by-site entries for that barcode,\nsince both axes fall when capture falls"),
    x = "Normalised allele entropy", y = "Mean edits per cell",
    colour = NULL) +
  ggplot2::guides(colour = ggplot2::guide_legend(
    order = 1, override.aes = list(size = 4))) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12.5),
    plot.subtitle = ggplot2::element_text(size = 8.8, colour = "#4A4A4A"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "right")

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = plot, device = ragg::agg_png,
                width = 10, height = 7, units = "in", dpi = 220,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = plot, device = grDevices::cairo_pdf,
                width = 10, height = 7, units = "in", bg = "#FFFFFF")
cat(sprintf("Wrote:\n  %s\n  %s\n", png_path, pdf_path))
