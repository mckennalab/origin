#!/usr/bin/env Rscript

# Per-position editing rate for the clone 84 BASELINE recording.
#
# The input is a cell x position matrix, 34 barcodes each carrying 272 target
# sites, so 9,248 positions. Each entry is one of:
#
#   0   site observed, unedited
#   1   site observed, edited
#   ?   site not called in this cell
#   +   the whole barcode was not captured in this cell
#
# Editing rate is computed over OBSERVED entries only -- n1 / (n0 + n1). Both ?
# and + are excluded from the denominator rather than counted as unedited:
# scoring an uncalled site as "no edit" would understate the rate by whatever
# the dropout happens to be, and dropout is not uniform across positions.
#
# Because the denominator varies per position, the number of cells actually
# supporting each rate is reported alongside it. A position observed in 30 cells
# and one observed in 9,000 both yield a "rate", but only one of them is worth
# much, and the low-coverage ones sit disproportionately in the low-rate tail
# that the filter removes.
#
# Arguments:
#   --input=<csv>       cell x position matrix (required)
#   --output-dir=<dir>  destination for the table (required)
#   --figure-dir=<dir>  destination for plots. Default: --output-dir
#   --min-rate=<p>      drop positions below this editing rate. Default 0.0005
#                       (0.05%). Reported alongside a sensitivity table, since
#                       "too low" is a judgement the threshold should not hide.
#   --min-observed=<n>  positions observed in fewer than this many cells are
#                       reported but flagged. Default 0 (no flagging).
#   --generations=<n>   convert the cumulative rate to a per-division rate under
#                       an n-generation assumption. Default 30; 0 disables.
#
# Per-division rate. Base editing is irreversible, so a site is unedited at
# sampling only if it escaped every division: p = 1 - (1 - r)^n, hence
# r = 1 - (1 - p)^(1/n). The transform is monotone, so the median per-division
# rate equals the transform of the median cumulative rate -- but the MEAN is not
# preserved, and is dominated by saturated positions where r is barely
# identifiable at all. At p = 1 the estimate is r = 1 regardless of the true
# rate, because anything fast enough to saturate by generation n looks the same.
# Report the median, and read the mean only alongside the saturated count.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
input_path <- value_after("--input")
output_dir <- value_after("--output-dir")
if (is.null(input_path) || is.null(output_dir)) {
  stop("--input and --output-dir are required.", call. = FALSE)
}
input_path <- normalizePath(input_path, mustWork = TRUE)
min_rate <- as.numeric(value_after("--min-rate", "0.0005"))
min_observed <- as.integer(value_after("--min-observed", "0"))
generations <- as.numeric(value_after("--generations", "30"))

if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) {
  cache_dir <- file.path(tempdir(), "xdg-cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = cache_dir)
}
for (package in c("data.table", "ggplot2", "ragg")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
suppressPackageStartupMessages(library(ggplot2))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
figure_dir <- value_after("--figure-dir", output_dir)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- normalizePath(figure_dir, mustWork = TRUE)

# Everything is read as character: the matrix is not numeric, and letting fread
# type-guess would turn a column of "0"/"1" into integers while leaving a column
# containing "?" as character, so the tallies below would need two code paths.
cat("[rates] reading matrix\n")
matrix_table <- data.table::fread(input_path, header = TRUE,
                                  colClasses = "character",
                                  showProgress = FALSE)
cell_ids <- matrix_table[[1L]]
data.table::set(matrix_table, j = 1L, value = NULL)
cat(sprintf("[rates] %d cells x %d positions\n", length(cell_ids),
            ncol(matrix_table)))

# Tally one column at a time. A vectorised `matrix_table == "1"` would allocate
# a 9,349 x 9,248 logical matrix (~350 MB) per comparison; this stays flat.
tally <- vapply(matrix_table, function(column) {
  c(unedited = sum(column == "0"), edited = sum(column == "1"),
    uncalled = sum(column == "?"), not_captured = sum(column == "+"))
}, numeric(4))

positions <- data.frame(
  position = colnames(matrix_table),
  barcode = sub("-[0-9]+$", "", colnames(matrix_table)),
  site = as.integer(sub("^.*-", "", colnames(matrix_table))),
  unedited = tally["unedited", ], edited = tally["edited", ],
  uncalled = tally["uncalled", ], not_captured = tally["not_captured", ],
  stringsAsFactors = FALSE
)
positions$observed <- positions$unedited + positions$edited
# A position observed in no cell has no rate to report, as distinct from a rate
# of zero. Keeping these as NA stops them from inflating the low-rate tail.
positions$editing_rate <- ifelse(positions$observed > 0,
                                 positions$edited / positions$observed, NA_real_)
positions$observed_fraction <- positions$observed / length(cell_ids)
rownames(positions) <- NULL

unexpected <- sum(rowSums(tally)) - length(cell_ids) * ncol(matrix_table)
if (unexpected != 0) {
  stop(sprintf("%d entries were none of 0/1/?/+; the tallies would be wrong.",
               abs(unexpected)), call. = FALSE)
}

if (generations > 0) {
  positions$per_division_rate <- 1 - (1 - positions$editing_rate)^(1 / generations)
}
scored <- positions[!is.na(positions$editing_rate), , drop = FALSE]
scored$retained <- scored$editing_rate >= min_rate
retained <- scored[scored$retained, , drop = FALSE]
positions$retained <- !is.na(positions$editing_rate) &
  positions$editing_rate >= min_rate
utils::write.csv(positions, file.path(output_dir,
                                      "clone84_position_editing_rates.csv"),
                 row.names = FALSE)

#' Summarise a set of positions in one line.
describe <- function(data, label) {
  quantiles <- stats::quantile(data$editing_rate, c(0.25, 0.5, 0.75))
  cat(sprintf(
    "%-14s n=%5d  mean=%.4f  median=%.4f  IQR=%.4f-%.4f  min=%.5f  max=%.4f\n",
    label, nrow(data), mean(data$editing_rate), quantiles[2],
    quantiles[1], quantiles[3], min(data$editing_rate),
    max(data$editing_rate)))
}

cat("\n== Editing rate per position ==\n")
if (nrow(positions) > nrow(scored)) {
  cat(sprintf("%d position(s) observed in no cell; excluded from all summaries\n",
              nrow(positions) - nrow(scored)))
}
describe(scored, "pre-filter")
describe(retained, sprintf("post-filter"))
cat(sprintf("\nFilter: editing rate >= %g (%.3f%%) removed %d of %d positions (%.1f%%)\n",
            min_rate, 100 * min_rate, nrow(scored) - nrow(retained),
            nrow(scored), 100 * (1 - nrow(retained) / nrow(scored))))
cat(sprintf("Mean editing rate: %.4f before filtering, %.4f after\n",
            mean(scored$editing_rate), mean(retained$editing_rate)))

# "Too low" is a judgement call, so show what each candidate threshold costs
# rather than letting the default stand unexamined.
cat("\n== Threshold sensitivity ==\n")
for (candidate in c(0, 0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1)) {
  kept <- sum(scored$editing_rate >= candidate)
  cat(sprintf("  >= %-7g (%6.2f%%)  keeps %5d  drops %5d  mean rate %.4f\n",
              candidate, 100 * candidate, kept, nrow(scored) - kept,
              mean(scored$editing_rate[scored$editing_rate >= candidate])))
}

cat("\n== Coverage ==\n")
cat(sprintf("Observed entries per position: median %.0f of %d cells (%.1f%%)\n",
            stats::median(positions$observed), length(cell_ids),
            100 * stats::median(positions$observed_fraction)))
cat(sprintf("Uncalled '?' %.1f%% of all entries; not-captured '+' %.1f%%\n",
            100 * sum(positions$uncalled) / (length(cell_ids) * nrow(positions)),
            100 * sum(positions$not_captured) /
              (length(cell_ids) * nrow(positions))))
dropped <- scored[scored$editing_rate < min_rate, , drop = FALSE]
if (nrow(dropped)) {
  cat(sprintf("Dropped positions are observed in a median of %.0f cells, against %.0f for retained\n",
              stats::median(dropped$observed), stats::median(retained$observed)))
}
if (min_observed > 0L) {
  cat(sprintf("%d retained position(s) observed in fewer than %d cells\n",
              sum(retained$observed < min_observed), min_observed))
}

if (generations > 0) {
  cat(sprintf("\n== Per-division editing rate (%g generations) ==\n", generations))
  per_division_line <- function(data, label) {
    cat(sprintf("%-28s n=%5d  median %.6f (%.4f%%)  mean %.6f (%.4f%%)\n",
                label, nrow(data), stats::median(data$per_division_rate),
                100 * stats::median(data$per_division_rate),
                mean(data$per_division_rate),
                100 * mean(data$per_division_rate)))
  }
  per_division_line(scored, "all positions")
  per_division_line(retained, sprintf("post-filter (>=%.3f%%)", 100 * min_rate))
  # Above ~95% cumulative, a per-division rate is close to unidentifiable: many
  # different rates all saturate by generation n and land at the same p.
  identifiable <- retained[retained$editing_rate < 0.95, , drop = FALSE]
  per_division_line(identifiable, "post-filter, p < 95%")
  saturated <- scored[scored$editing_rate >= 0.99, , drop = FALSE]
  cat(sprintf("%d position(s) at p >= 99%% (r >= %.4f); %d at p = 1 where r is not identifiable\n",
              nrow(saturated), 1 - (1 - 0.99)^(1 / generations),
              sum(scored$editing_rate == 1)))
  cat(sprintf("These carry the mean: dropping them moves it from %.4f%% to %.4f%%\n",
              100 * mean(retained$per_division_rate),
              100 * mean(identifiable$per_division_rate)))
  cat("\nSensitivity to the generation count (median over p < 95%):\n")
  for (candidate in c(20, 25, 30, 35, 40)) {
    rate <- stats::median(1 - (1 - identifiable$editing_rate)^(1 / candidate))
    cat(sprintf("  %2d generations -> %.6f (%.4f%%)\n", candidate, rate,
                100 * rate))
  }
}

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11, base_family = "sans"))
panel_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(face = "bold", size = 11.5),
  plot.subtitle = ggplot2::element_text(size = 8.6, colour = "#4A4A4A"),
  panel.grid.minor = ggplot2::element_blank(),
  legend.position = "none"
)
accent <- "#B4436C"
reference <- "#1F3B4D"

#' Histogram of editing rate, on a shared x range so the two panels compare.
rate_panel <- function(data, title, subtitle, fill, show_threshold) {
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = editing_rate)) +
    ggplot2::geom_histogram(bins = 60, fill = fill, colour = NA, alpha = 0.88) +
    ggplot2::geom_vline(xintercept = mean(data$editing_rate), colour = reference,
                        linetype = "dashed", linewidth = 0.6) +
    ggplot2::annotate("text", x = mean(data$editing_rate), y = Inf,
                      label = sprintf("  mean %.3f", mean(data$editing_rate)),
                      hjust = 0, vjust = 1.7, size = 3.1, colour = reference) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    # coord_cartesian rather than scale limits: the two panels need a shared
    # range, but scale limits drop the bars that fall outside instead of
    # clipping the view, which silently removes counts from the histogram.
    ggplot2::coord_cartesian(xlim = range(scored$editing_rate)) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = "Editing rate (edited / observed)", y = "Positions") +
    panel_theme
  if (show_threshold && min_rate > 0) {
    plot <- plot + ggplot2::geom_vline(xintercept = min_rate, colour = "#087E8B",
                                       linewidth = 0.5)
  }
  plot
}

pre_panel <- rate_panel(
  scored, "Before filtering",
  sprintf("%d positions with at least one observation", nrow(scored)),
  accent, TRUE
)
post_panel <- rate_panel(
  retained, sprintf("After filtering at %.3f%%", 100 * min_rate),
  sprintf("%d positions retained, %d removed", nrow(retained),
          nrow(scored) - nrow(retained)),
  "#2E7096", FALSE
)

# The linear histograms put most of the mass in a few bins, so the tail the
# filter acts on is invisible there. A log axis is the only way to see what is
# actually being removed.
log_data <- scored[scored$editing_rate > 0, , drop = FALSE]
log_panel <- ggplot2::ggplot(log_data,
                             ggplot2::aes(x = editing_rate, fill = retained)) +
  ggplot2::geom_histogram(bins = 60, colour = NA, alpha = 0.88) +
  ggplot2::geom_vline(xintercept = min_rate, colour = "#087E8B",
                      linewidth = 0.5) +
  ggplot2::scale_x_log10(labels = function(x) sprintf("%g%%", 100 * x)) +
  ggplot2::scale_fill_manual(values = c(`FALSE` = "#C98A4B", `TRUE` = "#2E7096")) +
  ggplot2::labs(
    title = "The tail the filter acts on",
    subtitle = sprintf(
      "Log axis, non-zero rates only (%d position(s) at exactly zero are not shown); orange is removed",
      sum(scored$editing_rate == 0)
    ),
    x = "Editing rate (log)", y = "Positions"
  ) +
  panel_theme

combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
  patchwork::wrap_plots(
    patchwork::wrap_plots(pre_panel, post_panel, ncol = 2),
    log_panel, ncol = 1, heights = c(1, 0.85)
  ) +
    patchwork::plot_annotation(
      title = "Clone 84 BASELINE: editing rate per target position",
      subtitle = sprintf(
        "%d positions (%d barcodes x %d sites) across %d cells; rate is edited / observed, with '?' and '+' excluded from the denominator",
        nrow(positions), length(unique(positions$barcode)),
        length(unique(positions$site)), length(cell_ids)
      ),
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
        plot.tag = ggplot2::element_text(face = "bold", size = 13)
      )
    )
} else {
  pre_panel
}

png_path <- file.path(figure_dir, "clone84_editing_rate_distribution.png")
pdf_path <- file.path(figure_dir, "clone84_editing_rate_distribution.pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 12, height = 8.5, units = "in", dpi = 220,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 12, height = 8.5, units = "in", bg = "#FFFFFF")

cat(sprintf("\nWrote:\n  %s\n  %s\n  %s\n",
            file.path(output_dir, "clone84_position_editing_rates.csv"),
            png_path, pdf_path))
