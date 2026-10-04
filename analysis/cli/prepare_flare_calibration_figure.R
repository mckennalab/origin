#!/usr/bin/env Rscript

# Assemble everything the FLARE calibration figure needs into one tidy file.
#
# The story runs inputs -> calibration -> validation: what was measured from the
# recording and fitted, what the rate sweep did with it, and whether the result
# reproduces features nobody fitted.
#
# Usage:
#   Rscript prepare_flare_calibration_figure.R <output-csv>

output <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(output)) stop("Usage: prepare_flare_calibration_figure.R <output-csv>",
                        call. = FALSE)
suppressPackageStartupMessages(library(data.table))
source("load_origin.R")
harness <- "/Users/aaronmck/Desktop/code/clique_2025_12_10/rust_cmd/test_harness/tree_building"
source(file.path(harness, "flare-wt-crispr-recorder.R"))
scratch <- "/private/tmp/claude-502/-Users-aaronmck-Desktop-code-remote-mito-clean/e7a58bc0-7417-4b36-b490-9f0e1a671f3d/scratchpad"
rows <- list()
add <- function(panel, series, x, y) {
  rows[[length(rows) + 1L]] <<- data.frame(panel = panel, series = series,
                                           x = x, y = y,
                                           stringsAsFactors = FALSE)
}

# ---- A. deletion length: what was measured, and the mixture fitted to it -----
lengths <- scan(file.path(scratch, "dellen.txt"), quiet = TRUE)
breaks <- c(0, 1, 2, 3, 5, 10, 20, 26, 40, 60, 80, 100, 120, 140, 160, 180, 200)
mids <- (head(breaks, -1) + breaks[-1]) / 2
add("A", "measured", mids,
    as.numeric(table(cut(lengths, breaks))) / sum(lengths <= 200))
set.seed(1)
draws <- replicate(200000, {
  if (stats::runif(1) < 0.335) {
    sum(ceiling(stats::rgamma(2, shape = 0.5931, rate = 0.15473))) + 1L
  } else {
    sum(ceiling(stats::rgamma(2, shape = 4.1475, rate = 0.06407))) + 1L
  }
})
add("A", "fitted mixture", mids,
    as.numeric(table(cut(draws, breaks))) / sum(draws <= 200))

# ---- B. the rate sweep, and which target it was aimed at --------------------
for (set in list(c("flare_local_0.335", "marginal target"),
                 c("flare_cond_0.335", "conditional target"))) {
  path <- file.path("analysis/results", set[1], "depth_calibration.csv")
  if (!file.exists(path)) next
  sweep <- read.csv(path)
  add("B", set[2], sweep$multiplier, sweep$target_rate)
}
add("B", "real, marginal", c(0.05, 20), c(0.8788, 0.8788))
add("B", "real, excl. suppressed", c(0.05, 20), c(0.9924, 0.9924))

# ---- C. edited targets per read ---------------------------------------------
real_profile <- c(0.1144, 0.0014, 0.0003, 0.0009, 0.0003, 0.0009, 0.0007,
                  0.0325, 0.8485)
add("C", "real", 0:8, real_profile)
add("C", "real, excl. suppressed", 0:8,
    c(0, real_profile[-1] / sum(real_profile[-1])))
molecules_of <- function(dir) {
  readRDS(sprintf("analysis/results/%s/calibrated_molecules.rds", dir))$molecules
}
for (set in list(c("flare_depth_18", "simulated, marginal target"),
                 c("flare_cond_0.335", "simulated, conditional target"))) {
  m <- molecules_of(set[1])
  h <- table(factor(rowSums(m != 0), levels = 0:8))
  add("C", set[2], 0:8, as.numeric(h) / sum(h))
}

# ---- D. targets removed by one event, matched encoding ----------------------
collapsed <- read.delim("analysis/results/flare_real/deletion_span_collapsed.tsv")
share <- collapsed$molecules[collapsed$largest_shared_span >= 1]
add("D", "real, excl. suppressed", 1:8, share / sum(share))
for (set in list(c("flare_depth_18", "simulated, marginal target"),
                 c("flare_cond_0.335", "simulated, conditional target"))) {
  m <- molecules_of(set[1])
  span <- apply(m, 1, function(r) {
    v <- r[r != 0]; if (!length(v)) 0L else max(table(v))
  })
  h <- table(factor(span[span >= 1], levels = 1:8))
  add("D", set[2], 1:8, as.numeric(h) / sum(h))
}

# ---- E. distinguishable outcomes, rarefied ----------------------------------
# Unique counts grow with sampling depth, so both sides are rarefied over the
# same molecule counts; an unrarefied comparison would say nothing.
real_tuples <- fread(file.path(scratch, "real_tuples.tsv"), header = FALSE,
                     sep = "\t", showProgress = FALSE)
alleles <- fread("analysis/results/flare_cond_0.335/calibrated/barcode_alleles.csv.gz",
                 showProgress = FALSE)
if (!is.numeric(alleles[[1L]])) set(alleles, j = 1L, value = NULL)
state <- as.matrix(alleles)
model <- list(wt_crispr = list(target_positions = 15L + 26L * (0:7),
                               barcode_length = 215L))
integrations <- length(unique(sub("_pos_.*$", "", colnames(state))))
characters <- as.matrix(flare_wt_crispr_logical_target_matrix(state, model,
                                                              integrations))
outcome <- function(mat) unlist(lapply(seq_len(integrations), function(i) {
  apply(mat[, (i - 1L) * 8L + seq_len(8L), drop = FALSE], 1, paste,
        collapse = "|")
}))
simulated <- outcome(characters)
binary <- outcome((characters != 0) * 1L)
set.seed(2)
n <- length(simulated)
index <- sample(nrow(real_tuples), min(n, nrow(real_tuples)))
grid <- c(250, 500, 1000, 2000, 5000, 10000, n)
grid <- grid[grid <= n]
add("E", "real assay", grid,
    vapply(grid, function(k) length(unique(real_tuples$V1[index[seq_len(k)]])),
           numeric(1)))
add("E", "simulated, allele states", grid,
    vapply(grid, function(k) length(unique(simulated[seq_len(k)])), numeric(1)))
add("E", "simulated, old binary encoding", grid,
    vapply(grid, function(k) length(unique(binary[seq_len(k)])), numeric(1)))

utils::write.csv(do.call(rbind, rows), output, row.names = FALSE)
cat(sprintf("Wrote %s (%d rows)\n", output, nrow(do.call(rbind, rows))))
