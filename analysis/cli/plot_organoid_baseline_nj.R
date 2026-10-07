#!/usr/bin/env Rscript

# BASELINE recorder accuracy on organoid ground-truth lineages, plotted.
#
# Every clone here has a known true tree, so normalised Robinson-Foulds is a
# direct measure of what the recorder and neighbour joining lose together.
# Lower is better; zero is an exact topology match.
#
# The panels answer, in order: how accurate is it overall, does accuracy depend
# on how many cells a clone has, does it depend on how much editing the clone
# actually accumulated, and is any of it affected by the mosaic skew dose.
#
# Usage:
#   Rscript plot_organoid_baseline_nj.R <results-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
results_dir <- normalizePath(
  if (length(arguments)) arguments[[1L]] else ".", mustWork = TRUE
)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(results_dir, "organoid_baseline_nj")
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

accuracy_path <- file.path(results_dir, "clone_accuracy.tsv.gz")
if (!file.exists(accuracy_path)) {
  stop("No clone_accuracy.tsv.gz in ", results_dir, call. = FALSE)
}
clones <- read.delim(gzfile(accuracy_path), stringsAsFactors = FALSE)
scored <- clones[clones$status == "ok" & is.finite(clones$normalized_rf), ,
                 drop = FALSE]
cat(sprintf("[plot] %d scored clones across %d organoids\n",
            nrow(scored), length(unique(scored$organoid_id))))

accent <- "#B4436C"
reference <- "#1F3B4D"
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.8, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  legend.position = "none"
)

overall_mean <- mean(scored$normalized_rf)
perfect_fraction <- mean(scored$normalized_rf == 0)

distribution_panel <- ggplot2::ggplot(
  scored, ggplot2::aes(x = normalized_rf)
) +
  ggplot2::geom_histogram(bins = 40, fill = accent, colour = NA, alpha = 0.85) +
  ggplot2::geom_vline(xintercept = overall_mean, colour = reference,
                      linetype = "dashed", linewidth = 0.6) +
  ggplot2::annotate("text", x = overall_mean, y = Inf,
                    label = sprintf("  mean %.3f", overall_mean),
                    hjust = 0, vjust = 1.6, size = 3.1, colour = reference) +
  ggplot2::labs(
    title = "Reconstruction accuracy per clone",
    subtitle = sprintf(
      "%d clones; %.0f%% recovered exactly (normalised RF = 0)",
      nrow(scored), 100 * perfect_fraction
    ),
    x = "Normalised Robinson-Foulds distance (0 = exact)", y = "Clones"
  ) +
  panel_theme

# Clone size is the dominant covariate: more tips means more splits to get
# right, and normalised RF is a fraction of those splits.
size_panel <- ggplot2::ggplot(
  scored, ggplot2::aes(x = leaves, y = normalized_rf)
) +
  ggplot2::geom_point(alpha = 0.18, size = 0.85, colour = accent) +
  ggplot2::geom_smooth(method = "loess", se = TRUE, colour = reference,
                       linewidth = 0.7, formula = y ~ x) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    title = "Accuracy against clone size",
    subtitle = "Each point is one clone; line is a loess fit",
    x = "Cells in clone (log)", y = "Normalised RF"
  ) +
  panel_theme

# How much signal the recorder actually wrote into this clone.
character_panel <- ggplot2::ggplot(
  scored, ggplot2::aes(x = informative_characters, y = normalized_rf)
) +
  ggplot2::geom_point(alpha = 0.18, size = 0.85, colour = "#087E8B") +
  ggplot2::geom_smooth(method = "loess", se = TRUE, colour = reference,
                       linewidth = 0.7, formula = y ~ x) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    title = "Accuracy against recorded signal",
    subtitle = "Characters edited in some but not all cells of the clone",
    x = "Informative characters (log)", y = "Normalised RF"
  ) +
  panel_theme

panels <- list(distribution_panel, size_panel, character_panel)

# The mosaic sweep varies a lineage-skew dose. If the recorder is unaffected by
# it, accuracy should be flat across skew; that is worth showing either way.
if ("skew" %in% names(scored) && sum(!is.na(scored$skew)) > 0) {
  skew_data <- scored[!is.na(scored$skew), , drop = FALSE]
  skew_summary <- do.call(rbind, lapply(
    split(skew_data, skew_data$skew), function(part) {
      n <- nrow(part)
      standard_error <- stats::sd(part$normalized_rf) / sqrt(n)
      data.frame(
        skew = part$skew[1], clones = n,
        mean_normalized_rf = mean(part$normalized_rf),
        lower = mean(part$normalized_rf) - 1.96 * standard_error,
        upper = mean(part$normalized_rf) + 1.96 * standard_error,
        stringsAsFactors = FALSE
      )
    }
  ))
  utils::write.csv(skew_summary, paste0(output_prefix, "_by_skew.csv"),
                   row.names = FALSE)
  panels[[4]] <- ggplot2::ggplot(
    skew_summary, ggplot2::aes(x = skew, y = mean_normalized_rf)
  ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = accent, alpha = 0.2) +
    ggplot2::geom_line(colour = accent, linewidth = 0.7) +
    ggplot2::geom_point(colour = accent, size = 1.8) +
    ggplot2::labs(
      title = "Accuracy against mosaic skew dose",
      subtitle = "Mean over clones with a 95% interval",
      x = "Lineage-skew dose", y = "Mean normalised RF"
    ) +
    panel_theme
}

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(panels, ncol = 2) +
    patchwork::plot_annotation(
      title = "BASELINE recorder accuracy on organoid ground-truth lineages",
      subtitle = sprintf(
        "%d clones from %d organoids, 5 integrations, neighbour joining, scored against the known clone tree",
        nrow(scored), length(unique(scored$organoid_id))
      ),
      caption = paste(
        "Each clone is reconstructed and scored independently: founders are",
        "separate day-zero cells, so a joined organoid tree would score the\n",
        "trivial between-clone splits alongside the within-clone structure at",
        "issue. Clones with fewer than four cells are excluded, since an\n",
        "unrooted topology on three tips has no internal split to get wrong."
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
  distribution_panel
}

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 12, height = 9, units = "in", dpi = 220,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 12, height = 9, units = "in", bg = "#FFFFFF")

size_bins <- cut(scored$leaves, breaks = c(3, 5, 10, 20, 50, 100, Inf),
                 labels = c("4-5", "6-10", "11-20", "21-50", "51-100", ">100"))
by_size <- do.call(rbind, lapply(split(scored, size_bins), function(part) {
  if (!nrow(part)) return(NULL)
  data.frame(size_range = as.character(size_bins[match(part$leaves[1],
                                                       part$leaves)])[1],
             clones = nrow(part),
             mean_normalized_rf = mean(part$normalized_rf),
             perfect = sum(part$normalized_rf == 0), stringsAsFactors = FALSE)
}))
by_size$size_range <- names(split(scored, size_bins))[
  match(by_size$clones, vapply(split(scored, size_bins), nrow, integer(1)))
]
utils::write.csv(by_size, paste0(output_prefix, "_by_size.csv"),
                 row.names = FALSE)

cat(sprintf("mean normalised RF %.4f  median %.4f  exact %.1f%%\n",
            overall_mean, stats::median(scored$normalized_rf),
            100 * perfect_fraction))
print(by_size, row.names = FALSE, digits = 3)
cat(sprintf("\nWrote:\n  %s\n  %s\n", png_path, pdf_path))
