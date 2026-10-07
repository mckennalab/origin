#!/usr/bin/env Rscript

# Reconstruction accuracy against clone size, for the organoid BASELINE runs.
#
# Accuracy here is 1 - normalised Robinson-Foulds, so higher is better and the
# scale matches how accuracy is reported elsewhere in the project.
#
# Two things make a plain scatter of 15,000 clones misleading, and both are
# handled rather than ignored:
#
#   overplotting   at this density a scatter shows only its own outline. The
#                  main panel bins the points and maps count to fill, so the
#                  mass of the distribution is visible instead of a silhouette.
#   discreteness   normalised RF on n tips can only take values k/(n-3), so
#                  small clones admit very few distinct values and produce
#                  diagonal banding that looks like structure but is arithmetic.
#                  The binned view absorbs it; the loess fit is drawn over the
#                  raw values regardless.
#
# The mean accuracy curve understates how sharply performance falls, because
# exact recovery collapses long before the mean does. The second panel reports
# that separately.
#
# Usage:
#   Rscript plot_organoid_accuracy_vs_size.R <results-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
results_dir <- normalizePath(
  if (length(arguments)) arguments[[1L]] else ".", mustWork = TRUE
)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(results_dir, "organoid_accuracy_vs_size")
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

clones <- read.delim(
  gzfile(file.path(results_dir, "clone_accuracy.tsv.gz")),
  stringsAsFactors = FALSE
)
scored <- clones[clones$status == "ok" & is.finite(clones$normalized_rf), ,
                 drop = FALSE]
scored$accuracy <- 1 - scored$normalized_rf
cat(sprintf("[plot] %d clones, sizes %d-%d\n", nrow(scored),
            min(scored$leaves), max(scored$leaves)))

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.8, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank()
)

scatter_panel <- ggplot2::ggplot(
  scored, ggplot2::aes(x = leaves, y = accuracy)
) +
  ggplot2::geom_bin2d(bins = 55) +
  ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
                       colour = "#B4436C", linewidth = 0.9) +
  ggplot2::scale_x_log10(
    breaks = c(4, 10, 30, 100, 300),
    labels = c("4", "10", "30", "100", "300")
  ) +
  ggplot2::scale_fill_gradient(
    name = "clones", low = "#DCE6EC", high = "#1F3B4D", trans = "log10"
  ) +
  ggplot2::labs(
    title = "Accuracy against clone size",
    subtitle = sprintf(
      "%d clones; fill is clone count per bin, line is a loess fit over the raw values",
      nrow(scored)
    ),
    x = "Cells in clone (log)",
    y = "Accuracy (1 - normalised RF)"
  ) +
  panel_theme

# Exact recovery, with a binomial interval. This is where the real degradation
# shows: the mean accuracy curve is almost flat over the same range.
size_breaks <- c(3, 5, 10, 20, 50, 100, Inf)
size_labels <- c("4-5", "6-10", "11-20", "21-50", "51-100", ">100")
scored$size_bin <- cut(scored$leaves, breaks = size_breaks, labels = size_labels)
by_size <- do.call(rbind, lapply(split(scored, scored$size_bin), function(part) {
  if (!nrow(part)) return(NULL)
  exact <- sum(part$normalized_rf == 0)
  interval <- stats::binom.test(exact, nrow(part))$conf.int
  data.frame(
    size_bin = part$size_bin[1], clones = nrow(part),
    mean_accuracy = mean(part$accuracy),
    exact_fraction = exact / nrow(part),
    lower = interval[1], upper = interval[2], stringsAsFactors = FALSE
  )
}))
by_size$size_bin <- factor(by_size$size_bin, levels = size_labels)

exact_panel <- ggplot2::ggplot(
  by_size, ggplot2::aes(x = size_bin, y = exact_fraction)
) +
  ggplot2::geom_col(fill = "#B4436C", alpha = 0.85, width = 0.68) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                         width = 0.18, colour = "#1F3B4D", linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("n=%d", clones)),
                     vjust = -0.6, size = 2.8, colour = "#4A4A4A") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                              expand = ggplot2::expansion(mult = c(0, 0.16))) +
  ggplot2::labs(
    title = "Exactly recovered clones",
    subtitle = "Mean accuracy barely moves over this range; exact recovery collapses",
    x = "Cells in clone", y = "Clones with normalised RF = 0"
  ) +
  panel_theme +
  ggplot2::theme(legend.position = "none")

# Informative characters track clone size almost perfectly, so an accuracy
# curve against character count is largely the size curve relabelled. Stating
# the confound is more useful than plotting it twice.
character_correlation <- stats::cor(scored$leaves,
                                    scored$informative_characters,
                                    method = "spearman")
confound_panel <- ggplot2::ggplot(
  scored, ggplot2::aes(x = leaves, y = informative_characters)
) +
  ggplot2::geom_bin2d(bins = 55) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  ggplot2::scale_fill_gradient(
    name = "clones", low = "#DCEBE8", high = "#0B4F49", trans = "log10"
  ) +
  ggplot2::labs(
    title = "Why signal and size cannot be separated here",
    subtitle = sprintf(
      "Spearman rho = %.2f: bigger clones accumulate more edits, so an accuracy curve against character count restates this one",
      character_correlation
    ),
    x = "Cells in clone (log)", y = "Informative characters (log)"
  ) +
  panel_theme

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(
    patchwork::wrap_plots(scatter_panel, exact_panel, ncol = 2,
                          widths = c(1.25, 1)),
    confound_panel, ncol = 1, heights = c(1, 0.85)
  ) +
    patchwork::plot_annotation(
      title = "BASELINE recorder: reconstruction accuracy against clone size",
      subtitle = sprintf(
        "%d clones from %d organoids, 5 integrations, neighbour joining against the known clone tree",
        nrow(scored), length(unique(scored$organoid_id))
      ),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
        plot.tag = ggplot2::element_text(face = "bold", size = 13)
      )
    )
} else {
  scatter_panel
}

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 12, height = 9.5, units = "in", dpi = 220,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 12, height = 9.5, units = "in", bg = "#FFFFFF")
utils::write.csv(by_size, paste0(output_prefix, "_by_size.csv"),
                 row.names = FALSE)

print(by_size, row.names = FALSE, digits = 3)
cat(sprintf("\nSpearman(clone size, informative characters) = %.3f\n",
            character_correlation))
cat(sprintf("Wrote:\n  %s\n  %s\n", png_path, pdf_path))
