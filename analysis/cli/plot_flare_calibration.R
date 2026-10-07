#!/usr/bin/env Rscript

# Calibrating the FLARE/WT-Cas9 recorder to the GSM8791703 HL60 recording.
#
# Reads left to right: what was measured and fitted (A), what the rate sweep did
# with it (B), and three features nobody fitted (C-E).
#
# Usage:
#   Rscript plot_flare_calibration.R <data-csv> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
data_path <- normalizePath(arguments[[1L]], mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(dirname(data_path), "flare_calibration")
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

REAL <- "#1F3B4D"; FIT <- "#C05746"; OLD <- "#9AA0A6"; ALT <- "#3E8E7E"
data <- read.csv(data_path, stringsAsFactors = FALSE)
part <- function(p) data[data$panel == p, , drop = FALSE]
theme_set(theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- theme(
  plot.title = element_text(face = "bold", size = 11.5),
  plot.subtitle = element_text(size = 8.4, colour = "#4A4A4A"),
  panel.grid.minor = element_blank(),
  legend.position = "bottom", legend.title = element_blank(),
  legend.margin = margin(t = -4))

a <- part("A")
a$series <- factor(a$series, levels = c("measured", "fitted mixture"))
panel_a <- ggplot(a, aes(x, y, colour = series)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.9) +
  scale_colour_manual(values = c("measured" = REAL, "fitted mixture" = FIT)) +
  scale_x_continuous(breaks = c(0, 26, 60, 100, 140, 180)) +
  expand_limits(y = 0) +
  labs(title = "A. Deletion length  [fitted]",
       subtitle = "Two classes, split at the 26bp target spacing: 33.5% repair locally, the rest resect",
       x = "Deletion length (bp)", y = "Share of events") + panel_theme

b <- part("B")
lines <- b[b$series %in% c("real, marginal", "real, excl. suppressed"), ]
points <- b[!b$series %in% c("real, marginal", "real, excl. suppressed"), ]
panel_b <- ggplot(points, aes(x, y, colour = series)) +
  geom_hline(data = lines, aes(yintercept = y, colour = series),
             linetype = "22", linewidth = 0.45) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 11)) +
  scale_colour_manual(values = c("marginal target" = OLD,
                                 "conditional target" = FIT,
                                 "real, marginal" = OLD,
                                 "real, excl. suppressed" = REAL)) +
  expand_limits(y = 0) +
  labs(title = "B. Rate calibration  [the target matters more than the knob]",
       subtitle = "The marginal 0.879 mixes in 11.4% of reads whose cells suppressed Cas9; the model cannot make those",
       x = "Editing rate, multiple of the parameter grid tier", y = "Per-target editing rate") +
  panel_theme

c_data <- part("C")
c_data <- c_data[!c_data$series %in% c("real", "simulated, marginal target"), ]
c_data$series <- factor(c_data$series, levels = c("real, excl. suppressed",
  "simulated, conditional target"))
panel_c <- ggplot(c_data, aes(factor(x), y, fill = series)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.66, alpha = 0.9) +
  scale_fill_manual(values = c("real, excl. suppressed" = REAL,
    "simulated, conditional target" = FIT)) +
  labs(title = "C. Edited targets per read  [not fitted]",
       subtitle = "The real array is all-or-nothing: 96% of touched reads lose all 8 targets, 0.4% sit in between",
       x = "Targets edited, of 8", y = "Share of reads") + panel_theme

d_data <- part("D")
d_data <- d_data[d_data$series != "simulated, marginal target", ]
d_data$series <- factor(d_data$series, levels = levels(c_data$series))
panel_d <- ggplot(d_data, aes(factor(x), y, fill = series)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.66, alpha = 0.9) +
  scale_fill_manual(values = c("real, excl. suppressed" = REAL,
    "simulated, conditional target" = FIT)) +
  labs(title = "D. Targets removed by one event  [not fitted]",
       subtitle = "Scored under the simulator's encoding on both sides, so the comparison is like for like",
       x = "Largest number of targets sharing one call", y = "Share of reads") +
  panel_theme

e_data <- part("E")
e_data$series <- factor(e_data$series, levels = c("real assay",
  "simulated, allele states", "simulated, old binary encoding"))
panel_e <- ggplot(e_data, aes(x, y, colour = series)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("real assay" = REAL,
    "simulated, allele states" = ALT,
    "simulated, old binary encoding" = OLD)) +
  labs(title = "E. Distinguishable outcomes, rarefied  [not fitted]",
       subtitle = "Both axes log. Unique counts grow with sampling depth, so the two sides are rarefied together",
       x = "Molecules sampled", y = "Unique whole-molecule outcomes") +
  panel_theme

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(
    patchwork::wrap_plots(panel_a, panel_b, ncol = 2),
    patchwork::wrap_plots(panel_c, panel_d, ncol = 2),
    panel_e, ncol = 1, heights = c(1, 1, 1)) +
    patchwork::plot_annotation(
      title = "Calibrating the FLARE recorder to the HL60 Cas9 recording",
      subtitle = "13.5M real molecules, 8 targets at 26bp spacing; deletion lengths fitted from the data, editing rate calibrated, everything else free",
      caption = paste(
        "Panels A and B are the calibration. A fits the deletion length mixture to 5.5M measured events; B sweeps the editing rate against both candidate",
        "\ntargets. The target was the choice that mattered: the 0.879 marginal mixes in the 11.4% of reads whose cells suppressed Cas9, which the model",
        "\ncannot produce, so hitting it forces the model to under-edit everything it does touch. Panels C to E show the run calibrated to the 0.992 rate",
        "\nconditional on being touched, and were not fitted. E needs the allele-state encoding: the binary matrix tree building used to read resolves an",
        "\norder of magnitude fewer outcomes than the construct actually produces."),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 9.5, colour = "#4A4A4A"),
                    plot.caption = element_text(size = 7.6, colour = "#4A4A4A",
                                                hjust = 0)))
} else panel_c
ggsave(paste0(output_prefix, ".png"), combined, device = ragg::agg_png,
       width = 13, height = 15, units = "in", dpi = 200, background = "#FFFFFF")
ggsave(paste0(output_prefix, ".pdf"), combined, device = grDevices::cairo_pdf,
       width = 13, height = 15, units = "in", bg = "#FFFFFF")
cat(sprintf("Wrote:\n  %s.png\n  %s.pdf\n", output_prefix, output_prefix))
