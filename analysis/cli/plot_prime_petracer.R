#!/usr/bin/env Rscript

# The prime slot before and after the PEtracer replacement.
#
# Two things changed together and the figure says so rather than implying one
# caused the other: the recorder gained an 8-mark alphabet per edit site (so an
# edit records WHICH mark, not just that one happened), and the geometry and
# rate ladder changed with it (3 sites x 8 marks against 6 binary sites).
#
# Lower normalised RF is better throughout, so improvement points downward.
# Panel B is the one to read carefully: the gain is large on three topologies
# and essentially absent on comb, and a pooled number would hide that.
#
# Usage:
#   Rscript plot_prime_petracer.R <grid-dir> <backup-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2L) {
  stop("Usage: plot_prime_petracer.R <grid-dir> <backup-dir> [prefix]",
       call. = FALSE)
}
grid_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
backup_dir <- normalizePath(file.path(arguments[[2L]], "summary_tables"),
                            mustWork = TRUE)
output_prefix <- if (length(arguments) >= 3L) arguments[[3L]] else {
  file.path(grid_dir, "prime_petracer_comparison")
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

tier_levels <- c("very_low", "low", "mid", "high", "very_high")
tier_labels <- c("VL", "L", "M", "H", "VH")
version_levels <- c("single-outcome prime", "PEtracer (8 marks)")
palette <- stats::setNames(c("#C98A4B", "#B4436C"), version_levels)

load_pair <- function(file) {
  new_data <- read.delim(file.path(grid_dir, file), stringsAsFactors = FALSE)
  old_data <- read.delim(file.path(backup_dir, file), stringsAsFactors = FALSE)
  new_data$version <- version_levels[2]
  old_data$version <- version_levels[1]
  combined <- rbind(new_data, old_data)
  combined <- combined[combined$system == "prime", , drop = FALSE]
  combined$version <- factor(combined$version, levels = version_levels)
  combined$tier <- factor(combined$parameter_setting, levels = tier_levels,
                          labels = tier_labels)
  combined
}

overall <- load_pair("accuracy_by_recorder_parameter.tsv")
by_shape <- load_pair("accuracy_by_recorder_parameter_shape.tsv")
by_method <- load_pair("accuracy_by_recorder_parameter_method.tsv")
cat(sprintf("[plot] %d overall rows, %d shape rows, %d method rows\n",
            nrow(overall), nrow(by_shape), nrow(by_method)))

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.6, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  legend.position = "bottom", legend.title = ggplot2::element_blank(),
  legend.key.height = grid::unit(0.8, "lines")
)

# The rate ladders are not the same numbers, so tiers are compared by RANK
# (VL..VH) and the actual rates are named in the subtitle.
has_interval <- all(c("ci95_lower_normalized_rf", "ci95_upper_normalized_rf")
                    %in% names(overall))
ladder_panel <- ggplot2::ggplot(
  overall, ggplot2::aes(x = tier, y = mean_normalized_rf, colour = version,
                        group = version)
)
if (has_interval) {
  ladder_panel <- ladder_panel +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci95_lower_normalized_rf,
                   ymax = ci95_upper_normalized_rf),
      width = 0.12, linewidth = 0.4, alpha = 0.85)
}
best_points <- do.call(rbind, lapply(split(overall, overall$version),
                                     function(part) {
  part[which.min(part$mean_normalized_rf), , drop = FALSE]
}))
ladder_panel <- ladder_panel +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::geom_point(data = best_points, shape = 21, size = 5, stroke = 1.1,
                      fill = NA, show.legend = FALSE) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::labs(
    title = "Editing-rate ladder, before and after",
    subtitle = paste(
      "Circled point is each version's best tier. Rates differ:",
      "old 0.01-0.16, new 0.0225-0.36 per cell cycle,\nso tiers align by rank rather than by value"
    ),
    x = "Editing-rate tier", y = "Mean normalised RF (lower is better)"
  ) +
  panel_theme

# Best tier per shape, as a dumbbell: the length of the segment is the gain.
shape_best <- do.call(rbind, lapply(
  split(by_shape, list(by_shape$shape, by_shape$version), drop = TRUE),
  function(part) part[which.min(part$mean_normalized_rf), , drop = FALSE]
))
shape_wide <- stats::reshape(
  shape_best[, c("shape", "version", "mean_normalized_rf")],
  idvar = "shape", timevar = "version", direction = "wide"
)
names(shape_wide) <- c("shape", "old", "new")
shape_wide <- shape_wide[order(shape_wide$new), ]
shape_wide$shape <- factor(shape_wide$shape, levels = shape_wide$shape)
shape_best$shape <- factor(shape_best$shape, levels = levels(shape_wide$shape))

shape_panel <- ggplot2::ggplot(shape_wide) +
  ggplot2::geom_segment(
    ggplot2::aes(x = old, xend = new, y = shape, yend = shape),
    colour = "#9A9A9A", linewidth = 0.9,
    arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed")
  ) +
  ggplot2::geom_point(data = shape_best,
                      ggplot2::aes(x = mean_normalized_rf, y = shape,
                                   colour = version), size = 3) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::labs(
    title = "Best tier by tree shape",
    subtitle = "Arrow runs from the old recorder to the new; comb barely moves",
    x = "Mean normalised RF at each version's best tier", y = NULL
  ) +
  panel_theme

method_best <- do.call(rbind, lapply(
  split(by_method, list(by_method$method, by_method$version), drop = TRUE),
  function(part) part[which.min(part$mean_normalized_rf), , drop = FALSE]
))
method_wide <- stats::reshape(
  method_best[, c("method", "version", "mean_normalized_rf")],
  idvar = "method", timevar = "version", direction = "wide"
)
names(method_wide) <- c("method", "old", "new")
method_wide <- method_wide[order(method_wide$new), ]
method_wide$method <- factor(method_wide$method, levels = method_wide$method)
method_best$method <- factor(method_best$method,
                             levels = levels(method_wide$method))

method_panel <- ggplot2::ggplot(method_wide) +
  ggplot2::geom_segment(
    ggplot2::aes(x = old, xend = new, y = method, yend = method),
    colour = "#9A9A9A", linewidth = 0.9,
    arrow = grid::arrow(length = grid::unit(0.10, "inches"), type = "closed")
  ) +
  ggplot2::geom_point(data = method_best,
                      ggplot2::aes(x = mean_normalized_rf, y = method,
                                   colour = version), size = 3) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::labs(
    title = "Best tier by tree-building method",
    subtitle = "Every method improves, which argues against a method-specific artefact",
    x = "Mean normalised RF at each version's best tier", y = NULL
  ) +
  panel_theme

# Where prime now sits in the panel. The other recorders are unchanged, so they
# are drawn once and prime twice.
new_all <- read.delim(file.path(grid_dir, "accuracy_by_recorder_parameter.tsv"),
                      stringsAsFactors = FALSE)
panel_best <- do.call(rbind, lapply(split(new_all, new_all$system),
                                    function(part) {
  best <- part[which.min(part$mean_normalized_rf), , drop = FALSE]
  data.frame(system = best$system[1], value = best$mean_normalized_rf[1],
             version = NA_character_, stringsAsFactors = FALSE)
}))
old_prime_best <- min(overall$mean_normalized_rf[
  overall$version == version_levels[1]])
panel_best$label <- ifelse(panel_best$system == "prime",
                           "prime (PEtracer)", panel_best$system)
panel_best <- rbind(panel_best, data.frame(
  system = "prime_old", value = old_prime_best, version = NA_character_,
  label = "prime (old)", stringsAsFactors = FALSE))
panel_best <- panel_best[order(panel_best$value), ]
panel_best$label <- factor(panel_best$label, levels = rev(panel_best$label))
panel_best$highlight <- ifelse(panel_best$system == "prime", "new",
                               ifelse(panel_best$system == "prime_old", "old",
                                      "other"))

context_panel <- ggplot2::ggplot(
  panel_best, ggplot2::aes(x = value, y = label, fill = highlight)
) +
  ggplot2::geom_col(width = 0.65, alpha = 0.9) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", value)),
                     hjust = -0.15, size = 3, colour = "#3A3A3A") +
  ggplot2::scale_fill_manual(
    values = c(new = "#B4436C", old = "#C98A4B", other = "#9FB3BF"),
    guide = "none") +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.14))) +
  ggplot2::labs(
    title = "Where prime sits in the panel",
    subtitle = "Best tier per recorder; every other recorder is unchanged by this work",
    x = "Mean normalised RF", y = NULL
  ) +
  panel_theme

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(ladder_panel, shape_panel, method_panel,
                        context_panel, ncol = 2) +
    patchwork::plot_annotation(
      title = "Prime slot: single-outcome prime editing replaced by PEtracer",
      subtitle = paste(
        "3 edit sites per integration with an 8-mark alphabet, 5 integrations,",
        "250 cells, 4 tree shapes x 5 seeds x 6 methods"
      ),
      caption = paste(
        "Two changes are confounded here and the figure does not separate them: the mark alphabet (an edit records WHICH of 8 marks was installed,",
        "\nso independent edits at a site are usually distinguishable from shared ancestry) and the geometry and rate ladder (3 sites x 8 marks against",
        "\n6 binary sites). Read this as 'the prime slot now models PEtracer', not as a measurement of what marks alone are worth -- that needs a",
        "\nmarks_per_target = 1 run at identical geometry. Comb trees barely improve, so the gain is not uniform across lineage topologies."
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
  ladder_panel
}

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 13, height = 10, units = "in", dpi = 220,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 13, height = 10, units = "in", bg = "#FFFFFF")
cat(sprintf("Wrote:\n  %s\n  %s\n", png_path, pdf_path))
