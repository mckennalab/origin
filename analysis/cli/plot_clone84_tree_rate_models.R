#!/usr/bin/env Rscript

# Does calibrating the editing-rate distribution change the tree?
#
# Five independent simulations, each scored on 10 resamples of 250 cells and 5
# BASELINE recorders. The configurations share a byte-identical ground-truth
# tree within a seed, so a replicate varies only the rate model.
#
# Seed is the unit of replication and the figure is built around that. Within
# one seed the 10 replicates resample a single simulation, so a paired test
# there measures resampling noise: seeds 1 and 2 each produced a "significant"
# clone-assignment result, in opposite directions. Panel B therefore shows every
# seed's effect separately rather than pooling them into one interval, because
# the scatter between seeds is the result.
#
# Usage:
#   Rscript plot_clone84_tree_rate_models.R <results-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
results_dir <- normalizePath(
  if (length(arguments)) arguments[[1L]] else ".", mustWork = TRUE)
output_prefix <- if (length(arguments) >= 2L) arguments[[2L]] else {
  file.path(results_dir, "clone84_tree_rate_models")
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

display <- c(empirical_rates = "empirical rates",
             fitted_shape = "gamma, shape fitted",
             perfect = "perfect recorder",
             random = "random tree")
palette <- stats::setNames(c("#3E8E7E", "#B4436C", "#1F3B4D", "#9AA0A6"),
                           unname(display))

read_set <- function(suffix, tail) {
  path <- file.path(results_dir, sprintf("rate_model_trees%s_%s.csv", suffix, tail))
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE)
}
suffixes <- c("", paste0("_seed", 2:9))
replicates <- do.call(rbind, lapply(seq_along(suffixes), function(index) {
  frame <- read_set(suffixes[index], "per_replicate")
  if (is.null(frame)) return(NULL)
  frame$seed <- index
  frame
}))
if (is.null(replicates)) stop("No per-replicate results found.", call. = FALSE)
eras <- do.call(rbind, lapply(seq_along(suffixes), function(index) {
  frame <- read_set(suffixes[index], "by_era")
  if (is.null(frame)) return(NULL)
  frame$seed <- index
  frame
}))
seeds <- sort(unique(replicates$seed))
cat(sprintf("[plot] %d seeds x %d replicates\n", length(seeds),
            max(replicates$replicate)))

replicates$label <- factor(display[replicates$model], levels = unname(display))
by_seed <- function(metric) {
  do.call(rbind, lapply(split(replicates, list(replicates$seed,
                                               replicates$label), drop = TRUE),
                        function(part) data.frame(
    seed = part$seed[1], label = part$label[1],
    value = mean(part[[metric]], na.rm = TRUE), stringsAsFactors = FALSE)))
}

panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.4, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  legend.position = "none")
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))

# ---- A. topology accuracy, seed by seed --------------------------------------
topology <- by_seed("topology_accuracy")
references <- topology[topology$label %in% display[c("perfect", "random")], ]
models <- topology[topology$label %in% display[c("empirical_rates",
                                                 "fitted_shape")], ]
panel_a <- ggplot2::ggplot(models,
                           ggplot2::aes(x = factor(seed), y = value,
                                        colour = label)) +
  ggplot2::geom_hline(data = references,
                      ggplot2::aes(yintercept = value, colour = label),
                      linetype = "22", linewidth = 0.4, alpha = 0.55) +
  # Dodged: on three of five seeds the two models differ by less than the point
  # radius, and overplotting them would read as a single result.
  ggplot2::geom_line(ggplot2::aes(group = factor(seed)), colour = "#B0B0B0",
                     linewidth = 0.5,
                     position = ggplot2::position_dodge(width = 0.34)) +
  ggplot2::geom_point(size = 3,
                      position = ggplot2::position_dodge(width = 0.34)) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(
    title = "A. Topology accuracy, seed by seed",
    subtitle = "1 - normalised RF; dashed lines are the perfect and random references",
    x = "Simulation seed", y = "1 - nRF") +
  panel_theme + ggplot2::theme(legend.position = "bottom",
                               legend.title = ggplot2::element_blank())

# ---- B. effect size per metric, one point per seed ---------------------------
# The honest panel: a mean difference means little when the seeds straddle zero.
metrics <- c(topology_accuracy = "Topology accuracy",
             clade_recall = "Clade recall",
             triplet_agreement = "Triplet agreement",
             clone_ari_k25 = "Clone ARI, k=25",
             clone_ari_k20 = "Clone ARI, k=20",
             clone_ari_k10 = "Clone ARI, k=10",
             ordering_agreement = "Ordering agreement",
             cophenetic_spearman = "Cophenetic rho")
differences <- do.call(rbind, lapply(names(metrics), function(metric) {
  frame <- by_seed(metric)
  wide <- stats::reshape(frame, idvar = "seed", timevar = "label",
                         direction = "wide")
  data.frame(metric = metrics[[metric]], seed = wide$seed,
             difference = wide[[paste0("value.", display[["empirical_rates"]])]] -
               wide[[paste0("value.", display[["fitted_shape"]])]],
             stringsAsFactors = FALSE)
}))
means <- do.call(rbind, lapply(split(differences, differences$metric),
                               function(part) data.frame(
  metric = part$metric[1], difference = mean(part$difference),
  won = sum(part$difference > 0), total = nrow(part), stringsAsFactors = FALSE)))
order_by <- means$metric[order(means$difference)]
differences$metric <- factor(differences$metric, levels = order_by)
means$metric <- factor(means$metric, levels = order_by)
panel_b <- ggplot2::ggplot(differences,
                           ggplot2::aes(x = difference, y = metric)) +
  ggplot2::geom_vline(xintercept = 0, colour = "#4A4A4A", linewidth = 0.4) +
  ggplot2::geom_point(colour = "#9AA0A6", size = 2.4, alpha = 0.85) +
  ggplot2::geom_point(data = means, colour = "#3E8E7E", size = 4,
                      shape = 18) +
  ggplot2::geom_text(data = means, ggplot2::aes(
    x = max(differences$difference) * 1.18,
    label = sprintf("%d/%d", won, total)), size = 2.9, colour = "#4A4A4A") +
  ggplot2::labs(
    title = "B. Effect of calibrating the rate model, one point per seed",
    subtitle = "Grey: each seed. Green diamond: mean. Right: seeds favouring the empirical model",
    x = "Empirical minus fitted gamma", y = NULL) +
  panel_theme

# ---- C. clone assignment against k, with the ceiling -------------------------
clone_counts <- as.integer(sub(".*_k", "",
  grep("^clone_ari_k", names(replicates), value = TRUE)))
clone <- do.call(rbind, lapply(sort(clone_counts), function(k) {
  frame <- by_seed(sprintf("clone_ari_k%d", k))
  frame$k <- k
  frame
}))
clone_summary <- do.call(rbind, lapply(split(clone, list(clone$k, clone$label),
                                             drop = TRUE), function(part)
  data.frame(k = part$k[1], label = part$label[1], mean = mean(part$value),
             sd = stats::sd(part$value), stringsAsFactors = FALSE)))
clone_summary <- clone_summary[clone_summary$label != display[["random"]], ]
panel_c <- ggplot2::ggplot(clone_summary,
                           ggplot2::aes(x = factor(k), y = mean, colour = label,
                                        group = label)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd, ymax = mean + sd),
                         width = 0.12, linewidth = 0.5) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::scale_colour_manual(values = palette) +
  ggplot2::expand_limits(y = 0) +
  ggplot2::labs(
    title = "C. Clone assignment, and whether the cut is well posed",
    subtitle = "Mean +/- SD over seeds. Where the perfect recorder falls short of 1, the cut itself is ambiguous",
    x = "Clusters (k)", y = "Adjusted Rand index") +
  panel_theme + ggplot2::theme(legend.position = "bottom",
                               legend.title = ggplot2::element_blank())

# ---- D. what replicates: informative characters ------------------------------
characters <- do.call(rbind, lapply(split(
  replicates[!is.na(replicates$characters), ],
  list(replicates$seed[!is.na(replicates$characters)],
       replicates$label[!is.na(replicates$characters)]), drop = TRUE),
  function(part) data.frame(seed = part$seed[1], label = part$label[1],
                            value = mean(part$characters),
                            stringsAsFactors = FALSE)))
panel_d <- ggplot2::ggplot(characters,
                           ggplot2::aes(x = factor(seed), y = value,
                                        fill = label)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                    width = 0.68, alpha = 0.9) +
  ggplot2::scale_fill_manual(values = palette) +
  ggplot2::labs(
    title = "D. Informative characters, of 1360 sites",
    subtitle = "The one difference that replicates exactly: about 1.9x in every seed",
    x = "Simulation seed", y = "Characters with >1 state") +
  panel_theme + ggplot2::theme(legend.position = "bottom",
                               legend.title = ggplot2::element_blank())

# ---- E. clade recall by when the split arose ---------------------------------
panel_e <- NULL
if (!is.null(eras)) {
  eras$label <- factor(display[eras$model], levels = unname(display))
  eras <- eras[!is.na(eras$label) & eras$label != display[["random"]], ]
  eras$era <- factor(eras$era, levels = c("earliest", "early", "late", "latest"))
  era_summary <- do.call(rbind, lapply(split(eras, list(eras$era, eras$label),
                                             drop = TRUE), function(part)
    data.frame(era = part$era[1], label = part$label[1],
               mean = mean(part$recall), sd = stats::sd(part$recall),
               stringsAsFactors = FALSE)))
  panel_e <- ggplot2::ggplot(era_summary,
                             ggplot2::aes(x = era, y = mean, colour = label,
                                          group = label)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sd, ymax = mean + sd),
                           width = 0.12, linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::labs(
      title = "E. Clade recall by when the split arose",
      subtitle = "Mean +/- SD over seeds; deep splits are where the rate model matters most",
      x = NULL, y = "Fraction of true clades recovered") +
    panel_theme + ggplot2::theme(legend.position = "bottom",
                                 legend.title = ggplot2::element_blank())
}

panels <- Filter(Negate(is.null), list(panel_a, panel_b, panel_c, panel_d,
                                       panel_e))
combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(
    patchwork::wrap_plots(panel_a, panel_b, ncol = 2, widths = c(1, 1.15)),
    patchwork::wrap_plots(panel_c, panel_d, ncol = 2),
    panel_e, ncol = 1, heights = c(1, 1, 0.95)) +
    patchwork::plot_annotation(
      title = "Does calibrating the editing-rate model change the reconstructed tree?",
      subtitle = sprintf(
        "%d simulations, 10 resamples each of 250 cells and 5 BASELINE recorders; the rate models share a ground-truth tree within a seed",
        length(seeds)),
      caption = paste(
        "Seed is the unit of replication. Within one seed the replicates resample a single simulation, so a paired test there measures resampling",
        "\nnoise rather than the effect: seeds 1 and 2 each gave a significant clone-assignment result, in opposite directions. Panel B therefore plots",
        "\nevery seed rather than pooling. Topology accuracy and clade recall are the only differences that hold their sign across seeds, and they are",
        "\nthe same measurement twice. Clone assignment and distance geometry are unresolved at this sample size, and clone ARI at k=5 is not",
        "\ninterpretable at all: the perfect recorder reaches only 0.34 to 0.61 there, so the cut, not the recorder, is what the number describes."),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
        plot.caption = ggplot2::element_text(size = 7.6, colour = "#4A4A4A",
                                             hjust = 0)))
} else {
  panel_b
}
png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 13, height = 15, units = "in", dpi = 200,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 13, height = 15, units = "in", bg = "#FFFFFF")
cat(sprintf("Wrote:\n  %s\n  %s\n", png_path, pdf_path))
