#!/usr/bin/env Rscript

# Pooled feature comparison over every clone 84 barcode that clears an entropy
# threshold, against the ten-integration simulation.
#
# The matched run could only describe ten barcodes. This uses all 34, keeps
# those whose allele entropy clears the threshold, and compares them with the
# simulation. Two adjustments are forced by the sides no longer having the same
# number of integrations:
#
#   edits per cell is reported PER INTEGRATION. A cell carrying 24 barcodes
#   accumulates more edits than one carrying 10 for reasons that have nothing to
#   do with the recorder, so the raw count is not comparable across designs.
#
#   the simulation is masked with the missingness of the first ten qualifying
#   real barcodes, since a 24-barcode pattern cannot mask a 10-integration
#   matrix. Its capture therefore reflects those ten, not all 24.
#
# Selecting barcodes on entropy and then reporting entropy is circular: panel C
# agreement is guaranteed by the filter and is not evidence about the model. The
# features worth reading here are homoplasy and inter-target dependence, which
# the filter does not touch.
#
# Arguments:
#   --input=<csv>       clone 84 matrix (required)
#   --match-dir=<dir>   clone84_simmatch directory (required)
#   --integration-features=<csv>  per-integration table with real entropy,
#                       used to choose barcodes (required)
#   --output-dir=<dir>  destination (required)
#   --min-entropy=<p>   keep barcodes at or above this. Default 0.95.
#   --cells=<n>         Default 250.
#   --pairs=<n>         pairs per stratum. Default 20000.
#   --seed=<n>          Default 1.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
input_path <- value_after("--input")
match_dir <- value_after("--match-dir")
integration_path <- value_after("--integration-features")
output_dir <- value_after("--output-dir")
if (is.null(input_path) || is.null(match_dir) || is.null(integration_path) ||
    is.null(output_dir)) {
  stop("--input, --match-dir, --integration-features and --output-dir are required.",
       call. = FALSE)
}
min_entropy <- as.numeric(value_after("--min-entropy", "0.95"))
cell_count <- as.integer(value_after("--cells", "250"))
pair_count <- as.integer(value_after("--pairs", "20000"))
seed <- as.integer(value_after("--seed", "1"))
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
script_dir <- dirname(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(script_dir, "..", "clone84_features.R"))

per_integration <- utils::read.csv(integration_path, stringsAsFactors = FALSE)
real_rows <- per_integration[grepl("^real", per_integration$dataset), ,
                             drop = FALSE]
keep_barcodes <- real_rows$integration[
  real_rows$normalised_entropy >= min_entropy]
if (length(keep_barcodes) < 2L) {
  stop("Fewer than two barcodes clear the entropy threshold.", call. = FALSE)
}
cat(sprintf("[filtered] %d of %d barcodes clear entropy %.2f\n",
            length(keep_barcodes), nrow(real_rows), min_entropy))

cat("[filtered] reading clone 84\n")
full <- data.table::fread(input_path, header = TRUE, colClasses = "character",
                          showProgress = FALSE)
cell_ids <- full[[1L]]
data.table::set(full, j = 1L, value = NULL)
barcodes <- sub("-[0-9]+$", "", colnames(full))
set.seed(seed)
chosen_cells <- sort(sample(seq_along(cell_ids),
                            min(cell_count, length(cell_ids))))
keep_columns <- which(barcodes %in% keep_barcodes)
raw <- as.matrix(full[chosen_cells, keep_columns, with = FALSE])
rm(full)
real <- matrix(NA_integer_, nrow(raw), ncol(raw), dimnames = dimnames(raw))
real[raw == "0"] <- 0L
real[raw == "1"] <- 1L
rm(raw)
integration_of <- sub("-[0-9]+$", "", colnames(real))
real_units <- unique(integration_of)
cat(sprintf("[filtered] real: %d cells x %d sites across %d barcodes\n",
            nrow(real), ncol(real), length(real_units)))

# The simulation has ten integrations, so it is masked with the first ten
# qualifying barcodes rather than the whole filtered set.
mask_units <- real_units[seq_len(min(10L, length(real_units)))]
mask_columns <- which(integration_of %in% mask_units)
simulated_mask <- !is.na(real[, mask_columns, drop = FALSE])
load_simulated <- function(name) {
  path <- file.path(match_dir, name, "barcode_binary_score_matrix.csv.gz")
  if (!file.exists(path)) return(NULL)
  score <- data.table::fread(path, showProgress = FALSE)
  if (!is.numeric(score[[1L]])) data.table::set(score, j = 1L, value = NULL)
  matrix_form <- as.matrix(score)
  storage.mode(matrix_form) <- "integer"
  set.seed(seed)
  if (nrow(matrix_form) > nrow(real)) {
    matrix_form <- matrix_form[sort(sample(nrow(matrix_form), nrow(real))), ,
                               drop = FALSE]
  }
  if (ncol(matrix_form) != ncol(simulated_mask)) return(NULL)
  matrix_form[matrix_form != 0L] <- 1L
  matrix_form[!simulated_mask] <- NA_integer_
  matrix_form
}
datasets <- list(real = list(x = real, integration = integration_of))
for (name in c("shipped_shape0.5", "fitted_shape", "empirical_rates")) {
  simulated <- load_simulated(name)
  if (!is.null(simulated)) {
    datasets[[name]] <- list(x = simulated,
                             integration = sub("_pos_.*$", "",
                                               colnames(simulated)))
  }
}

rows <- list(); site_frames <- list(); cell_frames <- list()
allele_frames <- list(); pair_frames <- list()
for (name in names(datasets)) {
  x <- datasets[[name]]$x
  units <- datasets[[name]]$integration
  unit_count <- length(unique(units))
  rates <- site_rate(x)
  # Per integration, so a 24-barcode cell is comparable with a 10-barcode one.
  edits <- cell_edits(x) / unit_count
  alleles <- allele_diversity(x, units)
  pairs <- pair_statistics(x, units, pair_count, seed)
  site_frames[[name]] <- data.frame(dataset = name, site_rate = rates,
                                    stringsAsFactors = FALSE)
  cell_frames[[name]] <- data.frame(dataset = name, edits = edits,
                                    stringsAsFactors = FALSE)
  alleles$dataset <- name
  allele_frames[[name]] <- alleles
  if (!is.null(pairs)) {
    pairs$dataset <- name
    pair_frames[[name]] <- pairs
  }
  rows[[name]] <- data.frame(
    dataset = name, integrations = unit_count,
    site_rate_mean = mean(rates, na.rm = TRUE),
    site_rate_median = stats::median(rates, na.rm = TRUE),
    site_rate_saturated = mean(rates >= 0.99, na.rm = TRUE),
    edits_per_cell_per_integration = mean(edits),
    alleles_per_integration = mean(alleles$alleles),
    normalised_entropy = mean(alleles$normalised_entropy),
    top_allele_fraction = mean(alleles$top_allele_fraction),
    stringsAsFactors = FALSE)
  cat(sprintf("[filtered] %s done (%d integrations)\n", name, unit_count))
}
summary_table <- do.call(rbind, rows)
pair_table <- do.call(rbind, pair_frames)
utils::write.csv(summary_table, file.path(output_dir, "feature_summary.csv"),
                 row.names = FALSE)
utils::write.csv(do.call(rbind, site_frames),
                 file.path(output_dir, "site_rates.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, cell_frames),
                 file.path(output_dir, "cell_edits.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, allele_frames),
                 file.path(output_dir, "allele_diversity.csv"),
                 row.names = FALSE)
if (!is.null(pair_table)) {
  utils::write.csv(pair_table, file.path(output_dir, "pair_statistics.csv"),
                   row.names = FALSE)
}
cat("\n== Features 1-3 ==\n\n")
print(summary_table, row.names = FALSE, digits = 4)
if (!is.null(pair_table)) {
  cat("\n== Features 4-5 ==\n\n")
  print(pair_table[, c("dataset", "stratum", "pairs", "incompatible_fraction",
                       "mean_abs_association")], row.names = FALSE, digits = 4)
  cat("\nExcess association within an integration:\n")
  for (name in unique(pair_table$dataset)) {
    part <- pair_table[pair_table$dataset == name, ]
    cat(sprintf("  %-18s %+.4f\n", name,
                part$mean_abs_association[part$stratum == "within"] -
                  part$mean_abs_association[part$stratum == "between"]))
  }
}
cat(sprintf("\nWrote outputs under %s\n", output_dir))
