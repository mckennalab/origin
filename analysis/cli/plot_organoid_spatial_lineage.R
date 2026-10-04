#!/usr/bin/env Rscript

# Lineage against spatial organisation in the simulated organoids.
#
# The organoid counterpart of the Visium correlogram. Two scales, because they
# can disagree: clones can hold coherent territory while their internal branches
# are spatially scrambled, and the reverse is possible too.
#
# Every panel carries a permutation null that shuffles which cell sits where.
# That preserves the lineage structure and the point cloud exactly and destroys
# only the correspondence between them, so it isolates ancestry from the fact
# that cells are packed in a finite volume. Without it, a falling shared-founder
# curve could be pure geometry.
#
# Usage:
#   Rscript plot_organoid_spatial_lineage.R <results-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
results_dir <- normalizePath(
  if (length(arguments)) arguments[[1L]] else ".", mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(results_dir, "organoid_spatial_lineage")
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

bin_labels <- c("0-10", "10-20", "20-30", "30-40", "40-60", "60-80",
                "80-120", ">120")
between <- read.csv(file.path(results_dir, "between_clone_by_distance.csv"),
                    stringsAsFactors = FALSE)
within <- read.csv(file.path(results_dir, "within_clone_by_distance.csv"),
                   stringsAsFactors = FALSE)
clones <- read.csv(file.path(results_dir, "clone_correlations.csv"),
                   stringsAsFactors = FALSE)
for (frame in c("between", "within")) {
  assign(frame, within(get(frame), {
    distance_bin <- factor(distance_bin, levels = bin_labels)
  }))
}
cat(sprintf("[plot] %d distance bins, %d clones\n", nrow(between),
            nrow(clones)))

accent <- "#B4436C"
null_colour <- "#8A8A8A"
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.6, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  legend.position = "none"
)

# Log scale: the shared-founder fraction falls by more than two orders of
# magnitude, and on a linear axis every bin past 40 microns would read as zero.
between_panel <- ggplot2::ggplot(
  between, ggplot2::aes(x = distance_bin, y = value, group = 1)
) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = null_value), colour = null_colour,
                      linetype = "dashed", linewidth = 0.5) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(lower, 1e-5), ymax = upper),
                       fill = accent, alpha = 0.2) +
  ggplot2::geom_line(colour = accent, linewidth = 0.9) +
  ggplot2::geom_point(colour = accent, size = 2.2) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Between clones: do neighbours share a founder?",
    subtitle = "Dashed line is the permutation null (positions shuffled). Log axis; ribbon is a 95% interval",
    x = "3D separation (microns)", y = "Fraction of pairs sharing a founder"
  ) +
  panel_theme

within_panel <- ggplot2::ggplot(
  within, ggplot2::aes(x = distance_bin, y = value, group = 1)
) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "#087E8B",
                       alpha = 0.2) +
  ggplot2::geom_line(colour = "#087E8B", linewidth = 0.9) +
  ggplot2::geom_point(colour = "#087E8B", size = 2.2) +
  ggplot2::labs(
    title = "Within a clone: does distance track ancestry?",
    subtitle = "Cophenetic distance on the truth tree, in minutes of simulated time",
    x = "3D separation (microns)", y = "Mean lineage distance (minutes)"
  ) +
  panel_theme

# Per-clone correlations against their own permuted values: one point per clone
# shows whether the pooled effect is general or driven by a few clones.
paired <- rbind(
  data.frame(source = "observed", value = clones$spearman,
             stringsAsFactors = FALSE),
  data.frame(source = "positions shuffled", value = clones$null_mean,
             stringsAsFactors = FALSE))
paired$source <- factor(paired$source,
                        levels = c("observed", "positions shuffled"))
clone_panel <- ggplot2::ggplot(paired, ggplot2::aes(x = source, y = value,
                                                    fill = source)) +
  ggplot2::geom_violin(colour = NA, alpha = 0.65, scale = "width") +
  ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
                        colour = "#2A2A2A", linewidth = 0.4) +
  ggplot2::geom_hline(yintercept = 0, colour = "#8A8A8A", linetype = "dotted",
                      linewidth = 0.4) +
  ggplot2::scale_fill_manual(values = c(observed = accent,
                                        `positions shuffled` = null_colour)) +
  ggplot2::labs(
    title = "Per-clone spatial/lineage correlation",
    subtitle = sprintf("%d clones; %.0f%% of clones have a positive correlation",
                       nrow(clones), 100 * mean(clones$spearman > 0,
                                                na.rm = TRUE)),
    x = NULL, y = "Spearman rho (spatial vs lineage distance)"
  ) +
  panel_theme

# Bigger clones span more tissue, so the correlation could be a size effect.
size_panel <- ggplot2::ggplot(clones, ggplot2::aes(x = cells, y = spearman)) +
  ggplot2::geom_hline(yintercept = 0, colour = "#8A8A8A", linetype = "dotted",
                      linewidth = 0.4) +
  ggplot2::geom_point(colour = accent, alpha = 0.55, size = 1.6) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                       colour = "#1F3B4D", linewidth = 0.8) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    title = "Is the correlation just clone size?",
    subtitle = sprintf(
      "Spearman(clone size, correlation) = %.2f",
      suppressWarnings(stats::cor(clones$cells, clones$spearman,
                                  method = "spearman", use = "complete.obs"))),
    x = "Cells in clone (log)", y = "Spearman rho"
  ) +
  panel_theme

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(between_panel, within_panel, clone_panel, size_panel,
                        ncol = 2) +
    patchwork::plot_annotation(
      title = "Lineage and spatial organisation in the simulated organoids",
      subtitle = "Ground-truth lineages and 3D cell positions; no recorder or reconstruction involved",
      caption = paste(
        "The null shuffles which cell occupies which position, leaving the lineage and the point cloud untouched, so it removes only the",
        "\ncorrespondence between them. This matters because cells packed in a finite volume are near each other for reasons that have nothing to do",
        "\nwith ancestry. Panel A is the clone-territory question and panel B the within-clone question; they are different claims and a recorder",
        "\ncould in principle recover one without the other."
      ),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
        plot.caption = ggplot2::element_text(size = 7.8, colour = "#4A4A4A",
                                             hjust = 0),
        plot.tag = ggplot2::element_text(face = "bold", size = 13)))
} else {
  between_panel
}

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 12, height = 9.5, units = "in", dpi = 220,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 12, height = 9.5, units = "in", bg = "#FFFFFF")
cat(sprintf("Wrote:\n  %s\n  %s\n", png_path, pdf_path))
