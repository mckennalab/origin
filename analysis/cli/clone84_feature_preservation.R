#!/usr/bin/env Rscript

# Does the fitted BASELINE simulation reproduce clone 84's structure, beyond
# its marginal editing rate?
#
# The earlier fit matched the per-site rate distribution and the edit burden.
# Those are marginal summaries: a simulation can match both while getting the
# joint structure wrong. This checks five features, three of which the fit was
# never asked to reproduce and could therefore fail:
#
#   per-target editing frequency  the marginal the fit targeted; included so a
#                                 failure elsewhere cannot be blamed on the fit
#                                 having drifted.
#   edits per cell                the other fitted marginal, as a full
#                                 distribution rather than a mean.
#   outcome diversity / entropy   BASELINE sites are binary, so diversity lives
#                                 at the integration level: how many distinct
#                                 edit patterns a barcode takes across cells,
#                                 and how evenly they are used.
#   homoplasy                     measured tree-free by four-gamete
#                                 incompatibility. Under perfect tree-like
#                                 inheritance with irreversible edits, no pair
#                                 of characters shows all four of 00/01/10/11;
#                                 pairs that do require a repeated or reverted
#                                 edit. This needs no truth tree, so it applies
#                                 to the real data as well as the simulations.
#   inter-target dependence       within-integration against between-integration
#                                 pairwise association. Shared ancestry inflates
#                                 both equally, so any excess within an
#                                 integration is a cis effect -- co-editing or
#                                 shared capture -- rather than lineage.
#
# All three datasets are masked with the REAL missingness pattern before any
# feature is computed. The simulations have no dropout of their own, and every
# feature here is sensitive to missingness, so comparing unmasked simulation
# against masked data would confound the recorder with the assay.
#
# Arguments:
#   --input=<csv>       clone 84 matrix (required)
#   --match-dir=<dir>   clone84_simmatch directory (required)
#   --reference=<dir>   clone84_reference directory (required)
#   --output-dir=<dir>  destination (required)
#   --cells=<n>         cells to compare. Default 250.
#   --pairs=<n>         character pairs sampled per stratum. Default 20000.
#   --exclude-integrations=<list>  integration POSITIONS to drop (1-based, in
#                       column order), applied to every dataset so the geometry
#                       stays matched. Dropping one only from the real data
#                       would leave the comparison uneven.
#   --seed=<n>          Default 1.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
input_path <- value_after("--input")
match_dir <- value_after("--match-dir")
reference_dir <- value_after("--reference")
output_dir <- value_after("--output-dir")
if (is.null(input_path) || is.null(match_dir) || is.null(reference_dir) ||
    is.null(output_dir)) {
  stop("--input, --match-dir, --reference and --output-dir are required.",
       call. = FALSE)
}
cell_count <- as.integer(value_after("--cells", "250"))
pair_count <- as.integer(value_after("--pairs", "20000"))
seed <- as.integer(value_after("--seed", "1"))
excluded_integrations <- {
  raw <- value_after("--exclude-integrations", "")
  if (nzchar(raw)) as.integer(strsplit(raw, ",", fixed = TRUE)[[1L]]) else
    integer()
}
for (package in c("data.table")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

settings <- jsonlite::read_json(file.path(reference_dir,
                                          "reference_settings.json"),
                                simplifyVector = TRUE)
chosen_barcodes <- settings$barcodes
divisions <- settings$divisions

# ---- real data, subset exactly as the reference was ----
cat("[features] reading clone 84\n")
full <- data.table::fread(input_path, header = TRUE, colClasses = "character",
                          showProgress = FALSE)
cell_ids <- full[[1L]]
data.table::set(full, j = 1L, value = NULL)
barcodes <- sub("-[0-9]+$", "", colnames(full))
set.seed(seed)
chosen_cells <- sort(sample(seq_along(cell_ids),
                            min(cell_count, length(cell_ids))))
keep_columns <- which(barcodes %in% chosen_barcodes)
real_raw <- as.matrix(full[chosen_cells, keep_columns, with = FALSE])
rm(full)
real <- matrix(NA_integer_, nrow(real_raw), ncol(real_raw),
               dimnames = dimnames(real_raw))
real[real_raw == "0"] <- 0L
real[real_raw == "1"] <- 1L
observed_mask <- !is.na(real)
integration_of <- sub("-[0-9]+$", "", colnames(real_raw))
if (length(excluded_integrations)) {
  units <- unique(integration_of)
  if (any(excluded_integrations < 1L | excluded_integrations > length(units))) {
    stop("--exclude-integrations must be positions in 1..", length(units),
         call. = FALSE)
  }
  dropped <- units[excluded_integrations]
  keep <- !(integration_of %in% dropped)
  real <- real[, keep, drop = FALSE]
  observed_mask <- observed_mask[, keep, drop = FALSE]
  integration_of <- integration_of[keep]
  # The same positions are dropped from the simulations further down, because
  # load_simulated() checks its width against `real`.
  cat(sprintf("[features] excluded integration(s) %s: %s\n",
              paste(excluded_integrations, collapse = ", "),
              paste(dropped, collapse = ", ")))
}
cat(sprintf("[features] real: %d cells x %d sites, %.1f%% observed\n",
            nrow(real), ncol(real), 100 * mean(observed_mask)))

# ---- simulated runs, masked with the real pattern ----
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
  if (length(excluded_integrations)) {
    simulated_units <- unique(sub("_pos_.*$", "", colnames(matrix_form)))
    simulated_drop <- simulated_units[excluded_integrations]
    matrix_form <- matrix_form[, !(sub("_pos_.*$", "",
                                       colnames(matrix_form)) %in%
                                     simulated_drop), drop = FALSE]
  }
  if (ncol(matrix_form) != ncol(real)) {
    stop(name, ": ", ncol(matrix_form), " sites against the real ", ncol(real),
         call. = FALSE)
  }
  # Any non-zero simulated state counts as edited; the real data is binary.
  matrix_form[matrix_form != 0L] <- 1L
  matrix_form[!observed_mask] <- NA_integer_
  dimnames(matrix_form) <- dimnames(real)
  matrix_form
}
datasets <- list(real = real)
for (name in c("shipped_shape0.5", "fitted_shape", "empirical_rates")) {
  simulated <- load_simulated(name)
  if (!is.null(simulated)) datasets[[name]] <- simulated
}
cat(sprintf("[features] datasets: %s\n", paste(names(datasets), collapse = ", ")))

source(file.path(dirname(sub("^--file=", "", grep("^--file=",
  commandArgs(trailingOnly = FALSE), value = TRUE)[1])), "..",
  "clone84_features.R"))

rows <- list()
integration_frames <- list()
site_frames <- list()
cell_frames <- list()
allele_frames <- list()
pair_frames <- list()
for (name in names(datasets)) {
  x <- datasets[[name]]
  rates <- site_rate(x)
  edits <- cell_edits(x)
  alleles <- allele_diversity(x, integration_of)
  pairs <- pair_statistics(x, integration_of, pair_count, seed)
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
    dataset = name,
    site_rate_mean = mean(rates, na.rm = TRUE),
    site_rate_median = stats::median(rates, na.rm = TRUE),
    site_rate_saturated = mean(rates >= 0.99, na.rm = TRUE),
    edits_per_cell_mean = mean(edits),
    edits_per_cell_sd = stats::sd(edits),
    comparable_cells = mean(alleles$cells),
    comparable_sites = mean(alleles$sites),
    alleles_per_integration = mean(alleles$alleles),
    normalised_entropy = mean(alleles$normalised_entropy),
    top_allele_fraction = mean(alleles$top_allele_fraction),
    stringsAsFactors = FALSE)
  per_integration <- integration_features(x, integration_of,
                                          max(2000L, pair_count %/% 10L), seed)
  per_integration$dataset <- name
  integration_frames[[name]] <- per_integration
  cat(sprintf("[features] %s done\n", name))
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

integration_table <- do.call(rbind, integration_frames)
utils::write.csv(integration_table,
                 file.path(output_dir, "feature_by_integration.csv"),
                 row.names = FALSE)

cat("\n== Features 1-3: marginals and outcome diversity ==\n\n")
print(summary_table, row.names = FALSE, digits = 4)
if (!is.null(pair_table)) {
  cat("\n== Features 4-5: homoplasy and inter-target dependence ==\n\n")
  print(pair_table[, c("dataset", "stratum", "pairs", "incompatible_fraction",
                       "mean_abs_association")],
        row.names = FALSE, digits = 4)
  cat("\nExcess association within an integration (within minus between):\n")
  for (name in unique(pair_table$dataset)) {
    part <- pair_table[pair_table$dataset == name, ]
    within_value <- part$mean_abs_association[part$stratum == "within"]
    between_value <- part$mean_abs_association[part$stratum == "between"]
    if (length(within_value) && length(between_value)) {
      cat(sprintf("  %-18s %+.4f\n", name, within_value - between_value))
    }
  }
}
cat("\n== Per integration, real against fitted ==\n\n")
wide <- integration_table[integration_table$dataset %in%
                            c("real", "fitted_shape"), , drop = FALSE]
for (metric in c("site_rate_mean", "edits_per_cell_mean", "alleles",
                 "incompatible_fraction", "mean_abs_association")) {
  values <- stats::reshape(
    wide[, c("integration_index", "dataset", metric)],
    idvar = "integration_index", timevar = "dataset", direction = "wide")
  names(values) <- c("integration", "real", "fitted")
  cat(sprintf("%s:\n", metric))
  cat(sprintf("  real   %s\n",
              paste(sprintf("%7.3f", values$real), collapse = "")))
  cat(sprintf("  fitted %s\n",
              paste(sprintf("%7.3f", values$fitted), collapse = "")))
}
cat(sprintf("\nWrote outputs under %s\n", output_dir))
