#!/usr/bin/env Rscript

# Build a simulation-sized reference from the clone 84 BASELINE recording, and
# calibrate the simulator's mean per-target editing rate against it.
#
# The full recording is 9,349 cells x 34 barcodes x 272 sites. This subsets to a
# scale the simulator can match directly -- by default 10 barcodes and 250 cells
# -- and computes the two quantities the simulation has to reproduce:
#
#   edit burden    edits per cell, over observed sites only
#   heterogeneity  the distribution of per-site cumulative editing rate
#
# Calibration. The simulator draws a per-target per-division rate from
# Gamma(shape = 0.5, scale = m / 0.5) (draw_physicell_target_rates), so m is the
# MEAN per-division rate and the shape fixes the dispersion at CV = 1.41. Over n
# divisions a target's cumulative editing probability is 1 - (1 - r)^n. The
# value of m is therefore not free: it is whatever reproduces the observed
# cumulative distribution.
#
# m is not read off the empirical mean per-division rate, because that mean is
# corrupted by saturation -- targets edited in every cell invert to r = 1
# regardless of their true rate. Instead m is chosen to minimise the distance
# between the simulated and observed cumulative-rate distributions, which uses
# the whole distribution rather than a statistic the ceiling has distorted.
#
# Arguments:
#   --input=<csv>        clone 84 matrix (required)
#   --output-dir=<dir>   destination (required)
#   --integrations=<n>   barcodes to keep. Default 10.
#   --cells=<n>          cells to keep. Default 250.
#   --divisions=<n>      cell divisions in the time window. Default 30.
#   --seed=<n>           Default 1.

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
integration_count <- as.integer(value_after("--integrations", "10"))
cell_count <- as.integer(value_after("--cells", "250"))
divisions <- as.integer(value_after("--divisions", "30"))
seed <- as.integer(value_after("--seed", "1"))

for (package in c("data.table")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

cat("[reference] reading matrix\n")
full <- data.table::fread(input_path, header = TRUE, colClasses = "character",
                          showProgress = FALSE)
cell_ids <- full[[1L]]
data.table::set(full, j = 1L, value = NULL)
barcodes <- sub("-[0-9]+$", "", colnames(full))

set.seed(seed)
chosen_barcodes <- sort(sample(unique(barcodes), integration_count))
chosen_cells <- sort(sample(seq_along(cell_ids), min(cell_count, length(cell_ids))))
keep_columns <- which(barcodes %in% chosen_barcodes)
subset_table <- full[chosen_cells, keep_columns, with = FALSE]
subset_cells <- cell_ids[chosen_cells]
cat(sprintf("[reference] %d cells x %d positions (%d barcodes x %d sites)\n",
            nrow(subset_table), ncol(subset_table), length(chosen_barcodes),
            ncol(subset_table) / length(chosen_barcodes)))

# Per-site tallies, one column at a time for the reason given in the rate script.
tally <- vapply(subset_table, function(column) {
  c(unedited = sum(column == "0"), edited = sum(column == "1"),
    uncalled = sum(column == "?"), not_captured = sum(column == "+"))
}, numeric(4))
sites <- data.frame(
  position = colnames(subset_table),
  barcode = sub("-[0-9]+$", "", colnames(subset_table)),
  site = as.integer(sub("^.*-", "", colnames(subset_table))),
  unedited = tally["unedited", ], edited = tally["edited", ],
  uncalled = tally["uncalled", ], not_captured = tally["not_captured", ],
  stringsAsFactors = FALSE
)
sites$observed <- sites$unedited + sites$edited
sites$editing_rate <- ifelse(sites$observed > 0,
                             sites$edited / sites$observed, NA_real_)
sites$per_division_rate <- 1 - (1 - sites$editing_rate)^(1 / divisions)
rownames(sites) <- NULL

# Edit burden per cell, over observed sites only so that capture dropout does
# not read as a lower edit count.
per_cell <- vapply(seq_len(nrow(subset_table)), function(index) {
  row <- unlist(subset_table[index, ], use.names = FALSE)
  observed <- sum(row == "0") + sum(row == "1")
  edited <- sum(row == "1")
  c(observed = observed, edited = edited,
    fraction = if (observed > 0) edited / observed else NA_real_)
}, numeric(3))
cells <- data.frame(cell_id = subset_cells, observed = per_cell["observed", ],
                    edited = per_cell["edited", ],
                    edited_fraction = per_cell["fraction", ],
                    stringsAsFactors = FALSE)

scored <- sites[!is.na(sites$editing_rate), , drop = FALSE]
cat("\n== Reference subset ==\n")
cat(sprintf("sites            %d (%d with at least one observation)\n",
            nrow(sites), nrow(scored)))
cat(sprintf("editing rate     mean %.4f  median %.4f  IQR %.4f-%.4f\n",
            mean(scored$editing_rate), stats::median(scored$editing_rate),
            stats::quantile(scored$editing_rate, 0.25),
            stats::quantile(scored$editing_rate, 0.75)))
cat(sprintf("saturated >=99%%  %d (%.1f%%);  below 1%%  %d (%.1f%%)\n",
            sum(scored$editing_rate >= 0.99),
            100 * mean(scored$editing_rate >= 0.99),
            sum(scored$editing_rate < 0.01),
            100 * mean(scored$editing_rate < 0.01)))
cat(sprintf("edits per cell   mean %.1f  median %.0f  of %.0f observed sites\n",
            mean(cells$edited), stats::median(cells$edited),
            stats::median(cells$observed)))
cat(sprintf("capture          '+' %.1f%%  '?' %.1f%%\n",
            100 * sum(sites$not_captured) / (nrow(cells) * nrow(sites)),
            100 * sum(sites$uncalled) / (nrow(cells) * nrow(sites))))

# ---- Calibrate the simulator's mean per-target per-division rate -------------
#
# For a candidate m, the simulator's per-target rate is Gamma(0.5, m/0.5) and
# cumulative editing is 1 - (1 - r)^divisions. Compare that distribution against
# the observed one. Distance is the mean absolute difference over a dense grid
# of quantiles: it uses the whole distribution, and unlike a moment it is not
# dominated by the saturated tail.
observed_quantiles <- stats::quantile(scored$editing_rate,
                                      probs = seq(0.01, 0.99, by = 0.01))
simulate_cumulative <- function(mean_rate, draws = 200000L) {
  rate <- stats::rgamma(draws, shape = 0.5, scale = mean_rate / 0.5)
  rate <- pmin(rate, 1 - .Machine$double.eps)
  1 - (1 - rate)^divisions
}
distance_for <- function(mean_rate) {
  set.seed(seed)
  simulated <- simulate_cumulative(mean_rate)
  mean(abs(stats::quantile(simulated, probs = seq(0.01, 0.99, by = 0.01)) -
             observed_quantiles))
}
candidates <- exp(seq(log(0.0005), log(0.2), length.out = 60))
distances <- vapply(candidates, distance_for, numeric(1))
best <- candidates[which.min(distances)]
refined <- stats::optimize(distance_for,
                           interval = c(best / 2, best * 2))$minimum

set.seed(seed)
fitted <- simulate_cumulative(refined)
cat("\n== Calibration ==\n")
cat(sprintf("Gamma(shape=0.5, scale=m/0.5) over %d divisions\n", divisions))
cat(sprintf("fitted m = %.6f (%.4f%% per target per division)\n", refined,
            100 * refined))
cat(sprintf("quantile distance at fit: %.4f\n", distance_for(refined)))
cat("\n            observed   simulated\n")
for (probability in c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)) {
  cat(sprintf("  q%-5.2f   %8.4f   %9.4f\n", probability,
              stats::quantile(scored$editing_rate, probability),
              stats::quantile(fitted, probability)))
}
cat(sprintf("  mean     %8.4f   %9.4f\n", mean(scored$editing_rate),
            mean(fitted)))
cat(sprintf("  >=99%%    %8.4f   %9.4f\n", mean(scored$editing_rate >= 0.99),
            mean(fitted >= 0.99)))
cat(sprintf("  <1%%      %8.4f   %9.4f\n", mean(scored$editing_rate < 0.01),
            mean(fitted < 0.01)))

utils::write.csv(sites, file.path(output_dir, "reference_sites.csv"),
                 row.names = FALSE)
utils::write.csv(cells, file.path(output_dir, "reference_cells.csv"),
                 row.names = FALSE)
jsonlite_available <- requireNamespace("jsonlite", quietly = TRUE)
settings <- list(
  input = input_path, integrations = length(chosen_barcodes),
  barcodes = chosen_barcodes, cells = nrow(cells),
  sites_per_barcode = ncol(subset_table) / length(chosen_barcodes),
  divisions = divisions, seed = seed,
  fitted_mean_rate_per_target_per_division = refined
)
if (jsonlite_available) {
  jsonlite::write_json(settings, file.path(output_dir, "reference_settings.json"),
                       pretty = TRUE, auto_unbox = TRUE)
}
cat(sprintf("\nWrote:\n  %s\n  %s\n",
            file.path(output_dir, "reference_sites.csv"),
            file.path(output_dir, "reference_cells.csv")))
