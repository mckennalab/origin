#!/usr/bin/env Rscript

# Does the simulation reproduce clone 84's BASELINE recording?
#
# Three distributions are compared: the real 10-barcode subset, the framework as
# shipped (per-target rates from Gamma with shape fixed at 0.5), and the
# framework with the dispersion fitted to the data.
#
# The mean editing rate is close in every configuration, so plotting means alone
# would suggest an agreement the distributions do not support. What separates the
# configurations is the SHAPE of the per-site distribution -- specifically how
# much mass sits at the two extremes, since those are the sites that carry
# lineage information (early-editing saturated sites) or none at all. The panels
# are therefore built around the full distribution and its tails.
#
# Usage:
#   Rscript plot_clone84_simulation_match.R <match-dir> <reference-dir> [prefix]

arguments <- commandArgs(trailingOnly = TRUE)
match_dir <- normalizePath(
  if (length(arguments)) arguments[[1L]] else ".", mustWork = TRUE)
reference_dir <- normalizePath(
  if (length(arguments) >= 2L) arguments[[2L]] else ".", mustWork = TRUE)
output_prefix <- if (length(arguments) >= 3L) arguments[[3L]] else {
  file.path(match_dir, "clone84_simulation_match")
}

if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) {
  cache_dir <- file.path(tempdir(), "xdg-cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = cache_dir)
}
for (package in c("ggplot2", "ragg", "data.table")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
suppressPackageStartupMessages(library(ggplot2))

rates <- read.csv(file.path(match_dir, "site_rates.csv"), stringsAsFactors = FALSE)
summary_table <- read.csv(file.path(match_dir, "match_summary.csv"),
                          stringsAsFactors = FALSE)
display <- c(reference = "clone 84 (real)",
             shipped_shape0.5 = "simulated, shape 0.5 (shipped)",
             fitted_shape = "simulated, shape fitted")
rates$label <- factor(display[rates$configuration], levels = unname(display))
rates <- rates[!is.na(rates$label), , drop = FALSE]
palette <- stats::setNames(c("#1F3B4D", "#C98A4B", "#B4436C"), unname(display))
cat(sprintf("[plot] %d site rates across %d configurations\n", nrow(rates),
            length(unique(rates$label))))

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.6, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  legend.position = "bottom", legend.title = ggplot2::element_blank(),
  legend.key.height = grid::unit(0.8, "lines")
)

# Histograms, one row per configuration: the bimodality is the thing to see, and
# overlaying three of them hides exactly the tails that matter.
histogram_panel <- ggplot2::ggplot(rates,
                                   ggplot2::aes(x = editing_rate, fill = label)) +
  ggplot2::geom_histogram(bins = 50, colour = NA, alpha = 0.9) +
  ggplot2::facet_wrap(~label, ncol = 1, scales = "free_y") +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = "Per-site editing rate",
    subtitle = "The real construct is U-shaped: most sites barely edit, a fifth saturate",
    x = "Cumulative editing rate", y = "Sites"
  ) +
  panel_theme +
  ggplot2::theme(legend.position = "none")

ecdf_panel <- ggplot2::ggplot(rates,
                              ggplot2::aes(x = editing_rate, colour = label)) +
  ggplot2::stat_ecdf(linewidth = 0.85) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = "Cumulative distribution",
    subtitle = "Vertical gaps are the mismatch; the shipped curve rises too early and saturates too late",
    x = "Cumulative editing rate", y = "Fraction of sites"
  ) +
  panel_theme

# Quantile-quantile against the real data: a configuration that reproduces the
# distribution lies on the diagonal at every quantile, not just on average.
probabilities <- seq(0.01, 0.99, by = 0.01)
observed <- stats::quantile(
  rates$editing_rate[rates$configuration == "reference"], probabilities)
quantile_frame <- do.call(rbind, lapply(
  setdiff(unique(rates$configuration), "reference"), function(configuration) {
    data.frame(
      configuration = configuration, label = display[[configuration]],
      observed = as.numeric(observed),
      simulated = as.numeric(stats::quantile(
        rates$editing_rate[rates$configuration == configuration],
        probabilities)),
      stringsAsFactors = FALSE
    )
  }))
quantile_frame$label <- factor(quantile_frame$label, levels = unname(display))
quantile_panel <- ggplot2::ggplot(
  quantile_frame, ggplot2::aes(x = observed, y = simulated, colour = label)
) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                       colour = "#8A8A8A", linewidth = 0.5) +
  ggplot2::geom_point(size = 1.5, alpha = 0.85) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = "Quantile agreement with the real data",
    subtitle = "Dashed line is exact agreement; points are the 1st to 99th percentiles",
    x = "Observed editing rate", y = "Simulated editing rate"
  ) +
  panel_theme

# Edit burden. Rates are read off the same matrices, so this is a check that
# matching the distribution also matches the absolute number of edits per cell.
burden <- list()
reference_cells <- read.csv(file.path(reference_dir, "reference_cells.csv"),
                            stringsAsFactors = FALSE)
burden[["reference"]] <- data.frame(
  configuration = "reference",
  edited_fraction = reference_cells$edited / reference_cells$observed,
  stringsAsFactors = FALSE)
for (configuration in setdiff(unique(rates$configuration), "reference")) {
  path <- file.path(match_dir, configuration,
                    "barcode_binary_score_matrix.csv.gz")
  if (!file.exists(path)) next
  score <- data.table::fread(path, showProgress = FALSE)
  if (!is.numeric(score[[1L]])) data.table::set(score, j = 1L, value = NULL)
  matrix_form <- as.matrix(score)
  burden[[configuration]] <- data.frame(
    configuration = configuration,
    edited_fraction = rowSums(matrix_form == 1L) / ncol(matrix_form),
    stringsAsFactors = FALSE)
}
burden_frame <- do.call(rbind, burden)
burden_frame$label <- factor(display[burden_frame$configuration],
                             levels = unname(display))
# Quoted in the caption, so computed rather than written down: these move
# whenever the calibration or the simulation seed changes.
burden_means <- vapply(split(burden_frame$edited_fraction, burden_frame$configuration),
                       mean, numeric(1))
names(burden_means) <- c(reference = "real", shipped_shape0.5 = "shipped",
                         fitted_shape = "fitted")[names(burden_means)]
burden_panel <- ggplot2::ggplot(
  burden_frame, ggplot2::aes(x = label, y = edited_fraction, fill = label)
) +
  ggplot2::geom_violin(colour = NA, alpha = 0.75, scale = "width") +
  ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
                        colour = "#2A2A2A", linewidth = 0.4) +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_x_discrete(labels = function(x) gsub(", ", ",\n", x)) +
  ggplot2::labs(
    title = "Edit burden per cell",
    subtitle = "Fraction of a cell's sites carrying an edit; the fitted shape is close, the shipped one undershoots",
    x = NULL, y = "Edited fraction of sites"
  ) +
  panel_theme +
  ggplot2::theme(legend.position = "none")

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(
    histogram_panel,
    patchwork::wrap_plots(ecdf_panel, quantile_panel, burden_panel, ncol = 1),
    ncol = 2, widths = c(1, 1.15)
  ) +
    patchwork::plot_annotation(
      title = "Reproducing the clone 84 BASELINE recording in simulation",
      subtitle = sprintf(
        "10 integrations x 272 sites, %d cells, %d cell divisions; per-target rates drawn from Gamma(shape, mean/shape)",
        summary_table$cells[summary_table$configuration == "reference"],
        summary_table$generation[summary_table$configuration == "reference"]
      ),
      caption = paste(
        "Mean editing rate is matched by construction in both simulated configurations, so agreement there is not evidence of anything. The test is",
        "\nthe distribution. With the shipped dispersion constant (shape = 0.5) the simulator puts too much mass in the middle and almost none at",
        "\nsaturation, missing the early-editing sites that resolve deep splits. Fitting the shape recovers the observed distribution closely.",
        sprintf("\nEdit burden per cell separates them far less than the distribution does (%s), so burden alone is a weak test.",
                paste(sprintf("%.0f%% %s", 100 * burden_means, names(burden_means)), collapse = ", "))
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
  histogram_panel
}

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 13, height = 10, units = "in", dpi = 220,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 13, height = 10, units = "in", bg = "#FFFFFF")

print(summary_table[, c("configuration", "shape", "site_rate_mean",
                        "site_rate_median", "fraction_saturated",
                        "fraction_low", "quantile_distance")],
      row.names = FALSE, digits = 4)
cat(sprintf("\nWrote:\n  %s\n  %s\n", png_path, pdf_path))
