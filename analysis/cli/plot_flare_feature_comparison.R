#!/usr/bin/env Rscript

# The FLARE/WT-Cas9 recorder against the GSM8791703 HL60 recording.
#
# Both sides are molecule-resolved: one row is one integration read out of one
# cell, eight target calls wide. Depth is pinned at 36 divisions (~44 days at a
# 29h HL60 doubling) and the editing rate was scaled until per-target
# saturation matched, because calibrating depth at the grid's rate would have
# needed 109 divisions and 7e10 cells.
#
# Panel A is therefore the calibrated quantity and proves nothing on its own --
# it is here so a failure elsewhere cannot be blamed on the marginal drifting.
# B to E were never calibrated.
#
# Usage:
#   Rscript plot_flare_feature_comparison.R <real-dir> <sim-dir> [prefix]

arguments <- commandArgs(trailingOnly = TRUE)
real_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
sim_dir <- normalizePath(arguments[[2L]], mustWork = TRUE)
output_prefix <- if (length(arguments) >= 3L) arguments[[3L]] else {
  file.path(sim_dir, "flare_feature_comparison")
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

palette <- c("HL60 recording (real)" = "#1F3B4D",
             "real, re-scored as the simulator encodes" = "#6B8CAE",
             "simulated, rate calibrated" = "#C05746")
# The comparison tables live with the run; the figure may be written elsewhere,
# so the two paths are resolved separately.
long <- read.csv(file.path(sim_dir, "flare_feature_comparison_long.csv"),
                 stringsAsFactors = FALSE)
long$dataset <- factor(c(real = "HL60 recording (real)",
                         real_collapsed = "real, re-scored as the simulator encodes",
                         simulated = "simulated, rate calibrated")[long$dataset],
                       levels = names(palette))
settings <- read.csv(file.path(sim_dir, "calibrated_settings.csv"),
                     stringsAsFactors = FALSE)
part <- function(name) long[long$feature == name, , drop = FALSE]

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.4, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  legend.position = "none")

panel_a <- ggplot2::ggplot(part("per-target editing rate"),
                           ggplot2::aes(x = factor(level), y = value,
                                        colour = dataset, group = dataset)) +
  ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2.6) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(title = "A. Per-target editing rate  [calibrated]",
                subtitle = "The real array is uniform; the simulation is ragged even after its mean is matched",
                x = "Target", y = "Fraction edited") +
  panel_theme + ggplot2::theme(legend.position = "bottom",
                               legend.title = ggplot2::element_blank())

panel_b <- ggplot2::ggplot(part("distinct calls per molecule"),
                           ggplot2::aes(x = factor(level), y = value,
                                        fill = dataset)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78),
                    width = 0.7, alpha = 0.9) +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::labs(title = "B. Independent characters per molecule",
                subtitle = "Mid blue is the same real molecules under the simulator's encoding: collapsing deletions costs a third of the characters",
                x = "Distinct calls across the 8 targets", y = "Share of molecules") +
  panel_theme

panel_c <- ggplot2::ggplot(part("largest shared span"),
                           ggplot2::aes(x = factor(level), y = value,
                                        fill = dataset)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78),
                    width = 0.7, alpha = 0.9) +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::labs(title = "C. Targets removed by one event",
                subtitle = "Compare the simulation against mid blue, not navy; span 0 is an untouched array, 11.4% real and none simulated",
                x = "Largest number of targets sharing one call", y = "Share of molecules") +
  panel_theme

alleles <- part("distinct alleles per target")
panel_d <- ggplot2::ggplot(alleles, ggplot2::aes(x = factor(level), y = value,
                                                 fill = dataset)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.78),
                    width = 0.7, alpha = 0.9) +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(title = "D. Distinct alleles per target  [log scale]",
                subtitle = "Every simulated deletion collapses to one state, so distinct events stop being distinguishable",
                x = "Target", y = "Distinct alleles") +
  panel_theme

pairs <- rbind(
  transform(part("pairwise phi"), measure = "association (phi)"),
  transform(part("P(same allele | both edited)"),
            measure = "P(same allele | both edited)"))
pairs <- pairs[is.finite(pairs$value), ]
pair_summary <- do.call(rbind, lapply(
  split(pairs, list(pairs$level, pairs$dataset, pairs$measure), drop = TRUE),
  function(p) data.frame(level = p$level[1], dataset = p$dataset[1],
                         measure = p$measure[1], value = mean(p$value))))
panel_e <- ggplot2::ggplot(pair_summary,
                           ggplot2::aes(x = factor(level), y = value,
                                        colour = dataset, group = dataset)) +
  ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 2.4) +
  ggplot2::facet_wrap(~measure, nrow = 1) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(
    title = "E. How the targets relate to each other, by distance along the array",
    subtitle = "Real targets edit together and share the event; simulated targets are near-independent",
    x = "Gap between targets", y = NULL) +
  panel_theme + ggplot2::theme(legend.position = "bottom",
                               legend.title = ggplot2::element_blank())

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(
    patchwork::wrap_plots(panel_a, panel_b, ncol = 2),
    patchwork::wrap_plots(panel_c, panel_d, ncol = 2),
    panel_e, ncol = 1, heights = c(1, 1, 0.95)) +
    patchwork::plot_annotation(
      title = "Does the FLARE recorder reproduce the HL60 Cas9 recording?",
      subtitle = sprintf(
        "13,500,552 real molecules against %s simulated, 8 targets each; depth fixed at 36 divisions and the rate scaled %.1fx to match saturation",
        format(settings$molecules, big.mark = ","), settings$rate_multiplier),
      caption = paste(
        "Depth could not be the calibrated knob: at the grid's rate the simulator needs 109 divisions to reach the observed saturation, which is 7e10 cells.",
        "\nThe rate had to be raised 6x the parameter grid's selected tier instead, which is itself a finding about that tier. Panel A is the calibrated",
        "\nquantity. Panels B to E were not calibrated and are where the model is tested: the simulation cannot produce an untouched array, collapses",
        "\nevery deletion to one indistinguishable state, and leaves the eight targets near-independent where the real array edits as one unit."),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
        plot.caption = ggplot2::element_text(size = 7.6, colour = "#4A4A4A",
                                             hjust = 0)))
} else {
  panel_e
}
ggplot2::ggsave(paste0(output_prefix, ".png"), plot = combined,
                device = ragg::agg_png, width = 13, height = 14, units = "in",
                dpi = 200, background = "#FFFFFF")
ggplot2::ggsave(paste0(output_prefix, ".pdf"), plot = combined,
                device = grDevices::cairo_pdf, width = 13, height = 14,
                units = "in", bg = "#FFFFFF")
cat(sprintf("Wrote:\n  %s.png\n  %s.pdf\n", output_prefix, output_prefix))
