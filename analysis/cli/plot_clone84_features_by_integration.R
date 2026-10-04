#!/usr/bin/env Rscript

# The clone 84 feature comparison, broken out by integration.
#
# The pooled comparison can hide a split result: one poorly captured barcode,
# or a handful of unusually fast sites, can carry a summary that then reads as a
# property of the recorder. Per-integration values separate "the model is wrong"
# from "a few barcodes are unusual".
#
# Real barcodes and simulated integrations correspond by POSITION, not identity:
# the i-th sorted barcode against int_i. They are not the same molecule, so a
# per-integration difference is only interpretable as a distribution over
# integrations, never as a claim about one barcode.
#
# Usage:
#   Rscript plot_clone84_features_by_integration.R <results-dir> [prefix]

arguments <- commandArgs(trailingOnly = TRUE)
results_dir <- normalizePath(
  if (length(arguments)) arguments[[1L]] else ".", mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(results_dir, "clone84_features_by_integration")
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

display <- c(real = "clone 84 (real)",
             shipped_shape0.5 = "simulated, shape 0.5",
             fitted_shape = "simulated, shape fitted",
             empirical_rates = "simulated, empirical rates")
palette <- stats::setNames(c("#1F3B4D", "#C98A4B", "#B4436C", "#3E8E7E"), unname(display))
data <- read.csv(file.path(results_dir, "feature_by_integration.csv"),
                 stringsAsFactors = FALSE)
data$label <- factor(display[data$dataset], levels = unname(display))
data <- data[!is.na(data$label), , drop = FALSE]
cat(sprintf("[plot] %d integrations x %d datasets\n",
            length(unique(data$integration_index)), length(unique(data$label))))

metrics <- c(site_rate_mean = "Mean site editing rate",
             edits_per_cell_mean = "Edits per cell",
             alleles = "Distinct alleles",
             normalised_entropy = "Normalised entropy",
             incompatible_fraction = "Homoplasy (incompatible pairs)",
             mean_abs_association = "Within-integration association")
long <- do.call(rbind, lapply(names(metrics), function(metric) {
  data.frame(integration_index = data$integration_index, label = data$label,
             metric = factor(metrics[[metric]], levels = unname(metrics)),
             value = data[[metric]], stringsAsFactors = FALSE)
}))

# Zero baseline. Every metric here is a count, a rate or a fraction with a real
# zero, and the comparisons drawn from the figure are proportional ones -- twice
# the homoplasy, half the association -- which only read honestly off an axis
# that starts at zero.
#
# Normalised entropy is the exception. It is bounded at 1 and all three datasets
# sit between 0.95 and 1.0, so a zero baseline collapses that panel into three
# indistinguishable points and hides a difference that is real. It stays zoomed,
# and the caption says so.
zero_metrics <- setdiff(unname(metrics), unname(metrics[["normalised_entropy"]]))
zero_anchor <- function(frame, value_column) {
  anchor <- frame[!duplicated(frame$metric), , drop = FALSE]
  anchor <- anchor[as.character(anchor$metric) %in% zero_metrics, , drop = FALSE]
  anchor[[value_column]] <- 0
  anchor
}

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
profile_panel <- ggplot2::ggplot(
  long, ggplot2::aes(x = factor(integration_index), y = value, colour = label,
                     group = label)) +
  ggplot2::geom_line(linewidth = 0.7, alpha = 0.9) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::geom_blank(data = zero_anchor(long, "value"),
                      ggplot2::aes(x = factor(integration_index), y = value),
                      inherit.aes = FALSE) +
  ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::labs(
    title = "Every feature, integration by integration",
    subtitle = paste(
      "Integrations are matched by position (i-th sorted barcode against int_i),",
      "not by identity"),
    x = "Integration", y = NULL) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.subtitle = ggplot2::element_text(size = 8.6, colour = "#4A4A4A"),
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 9),
    legend.position = "bottom", legend.title = ggplot2::element_blank())

# Level and spread together: mean across integrations, with whiskers at one
# standard deviation.
#
# Showing spread alone invited a reading that contradicted the panel above --
# real barcodes have the LOWEST homoplasy and, on a CV scale, the highest
# variability, which looks like a conflict until you notice the two panels were
# answering different questions. Plotting mean and SD in one place removes the
# ambiguity.
#
# SD rather than coefficient of variation: CV divides by the mean, and these
# datasets differ twofold in mean on several metrics, so a CV comparison reports
# the level difference a second time through the denominator. On homoplasy that
# turned a 1.2x difference in spread into a 2.5x difference in bar height.
spread <- do.call(rbind, lapply(names(metrics), function(metric) {
  do.call(rbind, lapply(split(data, data$label), function(part) {
    values <- part[[metric]]
    values <- values[is.finite(values)]
    if (!length(values)) return(NULL)
    data.frame(metric = factor(metrics[[metric]], levels = unname(metrics)),
               label = part$label[1], mean = mean(values),
               spread = stats::sd(values), stringsAsFactors = FALSE)
  }))
}))
spread_panel <- ggplot2::ggplot(spread, ggplot2::aes(x = label, y = mean,
                                                     colour = label)) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = mean - spread, ymax = mean + spread),
    width = 0.16, linewidth = 0.6) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::geom_blank(data = zero_anchor(spread, "mean"),
                      ggplot2::aes(x = label, y = mean), inherit.aes = FALSE) +
  ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 3) +
  ggplot2::scale_colour_manual(values = palette) +
  # Dropping everything after the comma made both simulated series read
  # "simulated", and this panel has no legend to tell them apart.
  ggplot2::scale_x_discrete(labels = function(x) {
    sub("simulated, ", "sim, ", sub(" \\(real\\)", "", x))
  }) +
  ggplot2::labs(
    title = "Level and spread across integrations",
    subtitle = "Point is the mean over integrations, whiskers are one standard deviation, in each metric's own units; axes start at zero",
    x = NULL, y = "Mean +/- SD across integrations") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.subtitle = ggplot2::element_text(size = 8.6, colour = "#4A4A4A"),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 7.5, angle = 20, hjust = 1),
    strip.text = ggplot2::element_text(size = 9),
    legend.position = "none")

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(profile_panel, spread_panel, ncol = 1,
                        heights = c(1.5, 1)) +
    patchwork::plot_annotation(
      title = "Clone 84 feature preservation, by integration",
      subtitle = sprintf(
        "%d integrations x %d sites, 250 cells; simulations masked with the real missingness pattern",
        length(unique(data$integration_index)), max(data$sites, na.rm = TRUE)),
      caption = paste(
        sprintf("Integrations correspond by position, not identity, so these are distributions over integrations rather than %d paired comparisons.",
                length(unique(data$integration_index))),
        "\nThe lower row shows level and spread together -- mean over integrations with one-SD whiskers -- because spread alone reads as though it",
        "\ncontradicts the panel above: real barcodes have the lowest homoplasy AND vary the most. SD is used rather than coefficient of variation,",
        "\nsince the datasets differ twofold in mean on several metrics and dividing by the mean would report that level difference a second time.",
        "\nAxes start at zero so the proportional comparisons read correctly. Normalised entropy is the one exception: it is bounded at 1 and every",
        "\ndataset sits above 0.95, so a zero baseline would hide a real difference. Read that panel, and only that panel, as a zoom."),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
        plot.caption = ggplot2::element_text(size = 7.6, colour = "#4A4A4A",
                                             hjust = 0),
        plot.tag = ggplot2::element_text(face = "bold", size = 13)))
} else {
  profile_panel
}
png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 12, height = 13, units = "in", dpi = 200,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 12, height = 13, units = "in", bg = "#FFFFFF")
cat(sprintf("Wrote:\n  %s\n  %s\n", png_path, pdf_path))
