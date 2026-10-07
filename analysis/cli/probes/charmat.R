# Does the per-target encoding now keep what the profile knows?
source("load_origin.R")
H <- "/Users/aaronmck/Desktop/code/clique_2025_12_10/rust_cmd/test_harness/tree_building"
source(file.path(H, "flare-wt-crispr-recorder.R"))
suppressPackageStartupMessages(library(data.table))

a <- fread("analysis/results/flare_cond_0.335/calibrated/barcode_alleles.csv.gz", showProgress = FALSE)
if (!is.numeric(a[[1]])) set(a, j = 1L, value = NULL)
state <- as.matrix(a)
model <- list(wt_crispr = list(target_positions = 15L + 26L * (0:7),
                               barcode_length = 215L))
integrations <- length(unique(sub("_pos_.*$", "", colnames(state))))

t0 <- Sys.time()
chars <- flare_wt_crispr_logical_target_matrix(state, model, integrations)
cat(sprintf("built %d x %d character matrix in %.1fs\n", nrow(chars), ncol(chars),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

per_target <- apply(chars[, 1:8, drop = FALSE], 2, function(x) length(unique(x)))
cat(sprintf("\ndistinct states per target (integration 1): %s\n",
            paste(per_target, collapse = " ")))

# Whole-molecule outcomes, one integration at a time, as the assay reads them.
outcomes <- unlist(lapply(seq_len(integrations), function(i) {
  block <- chars[, (i - 1L) * 8L + seq_len(8L), drop = FALSE]
  apply(block, 1, paste, collapse = "|")
}))
cat(sprintf("unique whole-molecule outcomes: %d over %d molecules\n",
            length(unique(outcomes)), length(outcomes)))
cat("\nfor comparison, at 14,855 molecules:\n")
cat("  real assay, full resolution     411\n")
cat("  simulated, derived from profile 311\n")
cat("  simulated, OLD -1 encoding       28\n")
