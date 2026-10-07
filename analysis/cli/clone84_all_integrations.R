#!/usr/bin/env Rscript

# Per-integration features for every barcode in clone 84, not just the ten the
# simulation was matched to.
#
# The simulated construct carries ten integrations, so the matched comparison
# could only ever describe ten real barcodes. This widens the real side to all
# 34, which matters for any claim about how much barcodes DIFFER from one
# another: that spread was estimated from ten, and ten is few.
#
# The two sides are no longer geometry-matched, so this is a comparison of two
# distributions over integrations, not a paired one. Simulated rows are carried
# through unchanged from the matched run, where they were masked with the real
# missingness of the ten sampled barcodes; real rows here use each barcode's own
# missingness. That is the honest arrangement -- masking 34 real barcodes with a
# 10-barcode pattern would be meaningless -- but it does mean the simulated
# capture distribution reflects only those ten.
#
# Arguments:
#   --input=<csv>       clone 84 matrix (required)
#   --matched=<csv>     feature_by_integration.csv from the matched run, whose
#                       simulated rows are appended. Optional.
#   --output-dir=<dir>  destination (required)
#   --cells=<n>         cells to use. Default 250.
#   --pairs=<n>         character pairs sampled per integration. Default 2000.
#   --seed=<n>          Default 1.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
input_path <- value_after("--input")
matched_path <- value_after("--matched")
output_dir <- value_after("--output-dir")
if (is.null(input_path) || is.null(output_dir)) {
  stop("--input and --output-dir are required.", call. = FALSE)
}
cell_count <- as.integer(value_after("--cells", "250"))
pair_count <- as.integer(value_after("--pairs", "2000"))
seed <- as.integer(value_after("--seed", "1"))
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
script_dir <- dirname(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(script_dir, "..", "clone84_features.R"))

cat("[all] reading clone 84\n")
full <- data.table::fread(input_path, header = TRUE, colClasses = "character",
                          showProgress = FALSE)
cell_ids <- full[[1L]]
data.table::set(full, j = 1L, value = NULL)
# The same seed and draw as the matched run, so the cell set is identical and
# only the barcode set widens.
set.seed(seed)
chosen_cells <- sort(sample(seq_along(cell_ids),
                            min(cell_count, length(cell_ids))))
raw <- as.matrix(full[chosen_cells, ])
rm(full)
real <- matrix(NA_integer_, nrow(raw), ncol(raw), dimnames = dimnames(raw))
real[raw == "0"] <- 0L
real[raw == "1"] <- 1L
rm(raw)
integration_of <- sub("-[0-9]+$", "", colnames(real))
cat(sprintf("[all] %d cells x %d sites across %d integrations, %.1f%% observed\n",
            nrow(real), ncol(real), length(unique(integration_of)),
            100 * mean(!is.na(real))))

per_integration <- integration_features(real, integration_of, pair_count, seed)
per_integration$dataset <- "real_all"
per_integration$cells <- nrow(real)

combined <- per_integration
if (!is.null(matched_path) && file.exists(matched_path)) {
  matched <- utils::read.csv(matched_path, stringsAsFactors = FALSE)
  simulated <- matched[matched$dataset != "real", , drop = FALSE]
  if (nrow(simulated)) {
    if (!"cells" %in% names(simulated)) simulated$cells <- cell_count
    shared <- intersect(names(combined), names(simulated))
    combined <- rbind(combined[, shared, drop = FALSE],
                      simulated[, shared, drop = FALSE])
  }
}
utils::write.csv(combined, file.path(output_dir, "feature_by_integration.csv"),
                 row.names = FALSE)

report <- function(metric) {
  cat(sprintf("\n%s\n", metric))
  for (name in unique(combined$dataset)) {
    values <- combined[[metric]][combined$dataset == name]
    values <- values[is.finite(values)]
    if (!length(values)) next
    cat(sprintf("  %-18s n=%2d  mean %8.4f  sd %7.4f  CV %5.3f  range %.3f-%.3f\n",
                name, length(values), mean(values), stats::sd(values),
                stats::sd(values) / abs(mean(values)), min(values),
                max(values)))
  }
}
cat("\n== All real integrations against the simulated ten ==\n")
for (metric in c("site_rate_mean", "edits_per_cell_mean", "alleles",
                 "normalised_entropy", "incompatible_fraction",
                 "mean_abs_association")) {
  report(metric)
}
real_rows <- combined[combined$dataset == "real_all", , drop = FALSE]
low <- real_rows[order(real_rows$normalised_entropy), ][
  seq_len(min(5L, nrow(real_rows))), ]
cat("\nLowest-entropy real barcodes:\n\n")
print(low[, c("integration", "integration_index", "observed_fraction",
              "edits_per_cell_mean", "alleles", "normalised_entropy",
              "mean_abs_association")],
      row.names = FALSE, digits = 3)
cat(sprintf("\nWrote %s\n",
            file.path(output_dir, "feature_by_integration.csv")))
