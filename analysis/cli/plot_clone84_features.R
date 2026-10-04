#!/usr/bin/env Rscript

# Which features of clone 84 the fitted simulation preserves, and which it does
# not.
#
# The fit targeted two marginals: the per-site editing rate distribution and the
# edit burden per cell. Panels A and B are therefore a check that the fit still
# holds, not evidence about it. Panels C to E are the test: outcome diversity,
# homoplasy and inter-target dependence were never fitted, so the simulation is
# free to get them wrong, and where it does that is a statement about the
# generative model rather than about the calibration.
#
# All three datasets carry the real missingness pattern, so nothing here is
# explained by the simulations lacking dropout.
#
# Usage:
#   Rscript plot_clone84_features.R <results-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
results_dir <- normalizePath(
  if (length(arguments)) arguments[[1L]] else ".", mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(results_dir, "clone84_features")
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
relabel <- function(frame) {
  frame$label <- factor(display[frame$dataset], levels = unname(display))
  frame[!is.na(frame$label), , drop = FALSE]
}
sites <- relabel(read.csv(file.path(results_dir, "site_rates.csv"),
                          stringsAsFactors = FALSE))
cells <- relabel(read.csv(file.path(results_dir, "cell_edits.csv"),
                          stringsAsFactors = FALSE))
alleles <- relabel(read.csv(file.path(results_dir, "allele_diversity.csv"),
                            stringsAsFactors = FALSE))
pairs <- relabel(read.csv(file.path(results_dir, "pair_statistics.csv"),
                          stringsAsFactors = FALSE))
# The filtered whole-set run reports edits per integration; detect it from the
# summary rather than requiring the caller to remember which run this is.
summary_path <- file.path(results_dir, "feature_summary.csv")
per_integration_edits <- file.exists(summary_path) &&
  "edits_per_cell_per_integration" %in%
    names(read.csv(summary_path, nrows = 1, stringsAsFactors = FALSE))
cat(sprintf("[plot] %d datasets; edits reported %s\n",
            length(unique(sites$label)),
            ifelse(per_integration_edits, "per integration", "per cell")))

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.5, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  legend.position = "none")

site_panel <- ggplot2::ggplot(sites[!is.na(sites$site_rate), ],
                              ggplot2::aes(x = site_rate, colour = label)) +
  ggplot2::stat_ecdf(linewidth = 0.85) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(title = "A. Per-target editing frequency  [calibrated]",
                subtitle = "Cumulative distribution over sites",
                x = "Site editing rate", y = "Fraction of sites") +
  panel_theme + ggplot2::theme(legend.position = "bottom",
                               legend.title = ggplot2::element_blank())

cell_panel <- ggplot2::ggplot(cells, ggplot2::aes(x = edits, fill = label)) +
  ggplot2::geom_density(colour = NA, alpha = 0.55) +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::labs(title = "B. Edits per cell  [calibrated]",
                subtitle = if (per_integration_edits) {
                  "Per integration, so designs with different barcode counts compare"
                } else "Over observed sites only",
                x = if (per_integration_edits) "Edits per cell per integration"
                    else "Edits per cell",
                y = "Density") +
  panel_theme

allele_long <- rbind(
  data.frame(label = alleles$label, metric = "distinct alleles",
             value = alleles$alleles, stringsAsFactors = FALSE),
  data.frame(label = alleles$label, metric = "normalised entropy",
             value = alleles$normalised_entropy, stringsAsFactors = FALSE),
  data.frame(label = alleles$label, metric = "commonest allele",
             value = alleles$top_allele_fraction, stringsAsFactors = FALSE))
allele_long$metric <- factor(allele_long$metric,
                             levels = c("distinct alleles",
                                        "normalised entropy",
                                        "commonest allele"))
# Counts and fractions are read proportionally, so their panels are anchored at
# zero. Normalised entropy is not: it is bounded at 1 and every dataset sits
# above 0.95, so a zero baseline would render three identical boxes and hide a
# real difference. Panels D and E are bar charts and already start at zero.
allele_zero <- allele_long[!duplicated(allele_long$metric), , drop = FALSE]
allele_zero <- allele_zero[allele_zero$metric != "normalised entropy", ,
                           drop = FALSE]
allele_zero$value <- 0

allele_panel <- ggplot2::ggplot(allele_long,
                                ggplot2::aes(x = label, y = value,
                                             fill = label)) +
  ggplot2::geom_boxplot(width = 0.6, outlier.size = 0.7, alpha = 0.85) +
  ggplot2::geom_blank(data = allele_zero,
                      ggplot2::aes(x = label, y = value), inherit.aes = FALSE) +
  ggplot2::facet_wrap(~metric, scales = "free_y", nrow = 1) +
  ggplot2::scale_fill_manual(values = palette) +
  # Dropping everything after the comma left both simulated series reading
  # "simulated", with no legend on this panel to separate them.
  ggplot2::scale_x_discrete(labels = function(x) {
    sub("simulated, ", "sim, ", sub(" \\(real\\)", "", x))
  }) +
  ggplot2::labs(
    title = "C. Outcome diversity per integration  [not calibrated]",
    subtitle = "Distinct edit patterns across cells, and how evenly they are used; entropy is zoomed, the other two start at zero",
    x = NULL, y = NULL) +
  panel_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7.5, angle = 20,
                                                     hjust = 1))

homoplasy_panel <- ggplot2::ggplot(
  pairs, ggplot2::aes(x = stratum, y = incompatible_fraction, fill = label)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                    width = 0.7, alpha = 0.9) +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = "D. Homoplasy  [not calibrated]",
    subtitle = "Character pairs showing all four of 00/01/10/11, which no tree explains without a repeat",
    x = "Character pair stratum", y = "Incompatible pairs") +
  panel_theme + ggplot2::theme(legend.position = "bottom",
                               legend.title = ggplot2::element_blank())

dependence_panel <- ggplot2::ggplot(
  pairs, ggplot2::aes(x = stratum, y = mean_abs_association, fill = label)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                    width = 0.7, alpha = 0.9) +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::labs(
    title = "E. Inter-target dependence  [not calibrated]",
    subtitle = "Shared ancestry inflates both strata equally; excess within an integration is a cis effect",
    x = "Character pair stratum", y = "Mean |correlation|") +
  panel_theme

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(
    patchwork::wrap_plots(site_panel, cell_panel, ncol = 2),
    allele_panel,
    patchwork::wrap_plots(homoplasy_panel, dependence_panel, ncol = 2),
    ncol = 1, heights = c(1, 0.95, 1)) +
    patchwork::plot_annotation(
      title = "Which clone 84 features the BASELINE simulation preserves, by rate model",
      subtitle = sprintf(
        "%d integrations x %d comparable sites, 250 cells; simulations masked with the real missingness pattern",
        length(unique(alleles$integration[alleles$dataset == "real"])),
        max(alleles$sites, na.rm = TRUE)),
      caption = paste(
        "Three rate models: a gamma at the shipped dispersion of 0.5, a gamma with the dispersion fitted, and per-target rates sampled from the",
        "\nmeasured per-site rates of the full 9,349-cell recording. Panels A and B are the marginals each rate model was calibrated on, so",
        "\nagreement there is a check that the calibration holds, not evidence about the model. Panels C to E were never calibrated: the generative",
        "\nmodel is free to get them wrong, and a mismatch says something about the model rather than the fit. Homoplasy is measured without a",
        "\ntree, by four-gamete incompatibility, so it applies to the real data too."),
      tag_levels = NULL,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
        plot.caption = ggplot2::element_text(size = 7.6, colour = "#4A4A4A",
                                             hjust = 0)))
} else {
  site_panel
}
png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 13, height = 13, units = "in", dpi = 200,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 13, height = 13, units = "in", bg = "#FFFFFF")
cat(sprintf("Wrote:\n  %s\n  %s\n", png_path, pdf_path))
