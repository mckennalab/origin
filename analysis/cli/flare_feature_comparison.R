#!/usr/bin/env Rscript

# The FLARE recorder against the GSM8791703 HL60 recording, feature by feature.
#
# Both sides are molecule-resolved: one row is one integration read out of one
# cell, with eight target calls. The real side comes from a full-file scan of
# 13.5M PASS consensus molecules; the simulated side from a run at matched
# geometry whose DEPTH was calibrated to the observed per-target saturation,
# with the editing rate left at the parameter grid's selected tier.
#
# Because the rate was calibrated, per-target frequency agreeing proves
# nothing. The features that test the model are the ones describing how the
# eight targets relate to each other:
#
#   distinct calls per molecule  how many independent characters an integration
#                                actually yields, against a nominal eight
#   largest shared span          how many targets one event removes at once
#   allele spectrum              whether distinct events stay distinguishable
#   pairwise association         whether targets edit together, and whether
#                                co-edited targets carry the same event
#
# Usage:
#   Rscript flare_feature_comparison.R <real-dir> <sim-dir> [output-prefix]

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2L) {
  stop("Usage: flare_feature_comparison.R <real-dir> <sim-dir> [prefix]",
       call. = FALSE)
}
real_dir <- normalizePath(arguments[[1L]], mustWork = TRUE)
sim_dir <- normalizePath(arguments[[2L]], mustWork = TRUE)
output_prefix <- if (length(arguments) >= 3L) arguments[[3L]] else {
  file.path(sim_dir, "flare_feature_comparison")
}

simulated <- readRDS(file.path(sim_dir, "calibrated_molecules.rds"))
molecules <- simulated$molecules
targets <- ncol(molecules)
cat(sprintf("[compare] simulated %d molecules x %d targets at sim_length %.2f\n",
            nrow(molecules), targets, simulated$sim_length))

# ---- simulated features, computed exactly as the awk scan computed the real --
edited <- molecules != 0
sim_target_rate <- colMeans(edited)

row_stats <- t(apply(molecules, 1L, function(row) {
  distinct <- length(unique(row))
  values <- row[row != 0]
  span <- if (!length(values)) 0L else max(table(values))
  c(distinct = distinct, span = span, edited = sum(row != 0))
}))
sim_distinct <- table(factor(row_stats[, "distinct"], levels = seq_len(targets)))
sim_span <- table(factor(row_stats[, "span"], levels = 0:targets))

sim_alleles <- do.call(rbind, lapply(seq_len(targets), function(t) {
  counts <- table(molecules[, t])
  proportions <- counts / sum(counts)
  entropy <- -sum(proportions * log(proportions))
  # A target stuck in one state carries no diversity; log(1) would make that
  # 0/0 rather than the 0 it is.
  normalised <- if (length(counts) > 1L) entropy / log(length(counts)) else 0
  data.frame(target = t, distinct_alleles = length(counts), entropy = entropy,
             normalised_entropy = normalised,
             top_fraction = max(proportions),
             top_allele = names(counts)[which.max(counts)],
             stringsAsFactors = FALSE)
}))

sim_pairs <- do.call(rbind, lapply(seq_len(targets - 1L), function(a) {
  do.call(rbind, lapply(seq(a + 1L, targets), function(b) {
    ea <- edited[, a]; eb <- edited[, b]
    n11 <- sum(ea & eb); n10 <- sum(ea & !eb)
    n01 <- sum(!ea & eb); n00 <- sum(!ea & !eb)
    same <- sum(ea & eb & molecules[, a] == molecules[, b])
    data.frame(target_a = a, target_b = b, gap = b - a, n11 = n11, n10 = n10,
               n01 = n01, n00 = n00, same_call = same, stringsAsFactors = FALSE)
  }))
}))

# ---- real features -----------------------------------------------------------
real_rates <- read.delim(file.path(real_dir, "target_rates.tsv"))
real_distinct <- read.delim(file.path(real_dir, "distinct_calls.tsv"))
real_span <- read.delim(file.path(real_dir, "deletion_span.tsv"))
real_alleles <- read.delim(file.path(real_dir, "allele_spectrum.tsv"))
real_pairs <- read.delim(file.path(real_dir, "pairs.tsv"))
# The same real molecules re-scored under the simulator's encoding, where every
# deletion is one indistinguishable state. Without this the span and
# distinct-call panels compare different quantities: sharing "183D+64" means
# one event, sharing "-1" only means both targets are deleted.
collapsed_calls <- read.delim(file.path(real_dir, "distinct_calls_collapsed.tsv"))
collapsed_span <- read.delim(file.path(real_dir, "deletion_span_collapsed.tsv"))

phi_of <- function(frame) {
  for (column in c("n11", "n10", "n01", "n00", "same_call")) {
    frame[[column]] <- as.numeric(frame[[column]])
  }
  frame$phi <- (frame$n11 * frame$n00 - frame$n10 * frame$n01) /
    sqrt((frame$n11 + frame$n10) * (frame$n01 + frame$n00) *
           (frame$n11 + frame$n01) * (frame$n10 + frame$n00))
  frame$same_given_both <- frame$same_call / frame$n11
  frame
}
real_pairs <- phi_of(real_pairs)
sim_pairs <- phi_of(sim_pairs)
# phi needs both targets to vary. A target the simulation always edits makes the
# denominator zero, so those pairs are dropped and counted rather than silently
# becoming NaN in the mean.
usable <- function(frame) sum(is.finite(frame$phi))
cat(sprintf("[compare] phi defined on %d of %d real pairs, %d of %d simulated\n",
            usable(real_pairs), nrow(real_pairs), usable(sim_pairs),
            nrow(sim_pairs)))

share <- function(counts) as.numeric(counts) / sum(as.numeric(counts))
long <- rbind(
  data.frame(feature = "per-target editing rate", level = real_rates$target,
             dataset = "real", value = real_rates$rate),
  data.frame(feature = "per-target editing rate", level = seq_len(targets),
             dataset = "simulated", value = sim_target_rate),
  data.frame(feature = "distinct calls per molecule",
             level = real_distinct$distinct_calls, dataset = "real",
             value = share(real_distinct$molecules)),
  data.frame(feature = "distinct calls per molecule",
             level = collapsed_calls$distinct_calls, dataset = "real_collapsed",
             value = share(collapsed_calls$molecules)),
  data.frame(feature = "distinct calls per molecule", level = seq_len(targets),
             dataset = "simulated", value = share(sim_distinct)),
  data.frame(feature = "largest shared span",
             level = real_span$largest_shared_span, dataset = "real",
             value = share(real_span$molecules)),
  data.frame(feature = "largest shared span",
             level = collapsed_span$largest_shared_span,
             dataset = "real_collapsed",
             value = share(collapsed_span$molecules)),
  data.frame(feature = "largest shared span", level = 0:targets,
             dataset = "simulated", value = share(sim_span)),
  data.frame(feature = "distinct alleles per target",
             level = real_alleles$target, dataset = "real",
             value = real_alleles$distinct_alleles),
  data.frame(feature = "distinct alleles per target", level = seq_len(targets),
             dataset = "simulated", value = sim_alleles$distinct_alleles),
  data.frame(feature = "normalised entropy per target",
             level = real_alleles$target, dataset = "real",
             value = real_alleles$normalised_entropy),
  data.frame(feature = "normalised entropy per target", level = seq_len(targets),
             dataset = "simulated", value = sim_alleles$normalised_entropy),
  data.frame(feature = "pairwise phi", level = real_pairs$gap, dataset = "real",
             value = real_pairs$phi),
  data.frame(feature = "pairwise phi", level = sim_pairs$gap,
             dataset = "simulated", value = sim_pairs$phi),
  data.frame(feature = "P(same allele | both edited)", level = real_pairs$gap,
             dataset = "real", value = real_pairs$same_given_both),
  data.frame(feature = "P(same allele | both edited)", level = sim_pairs$gap,
             dataset = "simulated", value = sim_pairs$same_given_both)
)
utils::write.csv(long, paste0(output_prefix, "_long.csv"), row.names = FALSE)

weighted_mean <- function(levels, weights) {
  sum(as.numeric(levels) * share(weights))
}
headline <- data.frame(
  metric = c("mean per-target editing rate",
             "spread of per-target rate (max - min)",
             "mean distinct calls, real re-scored as the simulator encodes",
             "one event spans >= 7 targets, matched encoding",
             "mean distinct calls per molecule (of 8)",
             "molecules with no edited target",
             "molecules where one event spans >= 4 targets",
             "median distinct alleles per target",
             "mean normalised entropy per target",
             "mean pairwise phi between targets",
             "P(same allele | both edited)"),
  real = c(
    mean(real_rates$rate),
    diff(range(real_rates$rate)),
    weighted_mean(collapsed_calls$distinct_calls, collapsed_calls$molecules),
    sum(share(collapsed_span$molecules)[
      collapsed_span$largest_shared_span >= 7]),
    weighted_mean(real_distinct$distinct_calls, real_distinct$molecules),
    share(real_span$molecules)[real_span$largest_shared_span == 0],
    sum(share(real_span$molecules)[real_span$largest_shared_span >= 4]),
    stats::median(real_alleles$distinct_alleles),
    mean(real_alleles$normalised_entropy),
    mean(real_pairs$phi, na.rm = TRUE),
    sum(real_pairs$same_call) / sum(real_pairs$n11)),
  simulated = c(
    mean(sim_target_rate),
    diff(range(sim_target_rate)),
    weighted_mean(seq_len(targets), sim_distinct),
    sum(share(sim_span)[as.integer(names(sim_span)) >= 7]),
    weighted_mean(seq_len(targets), sim_distinct),
    share(sim_span)[1],
    sum(share(sim_span)[as.integer(names(sim_span)) >= 4]),
    stats::median(sim_alleles$distinct_alleles),
    mean(sim_alleles$normalised_entropy),
    mean(sim_pairs$phi, na.rm = TRUE),
    sum(sim_pairs$same_call) / sum(sim_pairs$n11)),
  stringsAsFactors = FALSE)
headline$ratio <- headline$simulated / headline$real
utils::write.csv(headline, paste0(output_prefix, "_headline.csv"),
                 row.names = FALSE)

cat("\n== FLARE: simulated against the HL60 recording ==\n\n")
print(headline, row.names = FALSE, digits = 4)
cat("\n== distinct calls per molecule (share) ==\n\n")
print(data.frame(calls = seq_len(targets),
                 real = share(real_distinct$molecules),
                 simulated = share(sim_distinct)), row.names = FALSE, digits = 3)
cat("\n== largest shared span (share) ==\n\n")
print(data.frame(span = 0:targets, real = share(real_span$molecules),
                 simulated = share(sim_span)), row.names = FALSE, digits = 3)
cat(sprintf("\nWrote:\n  %s\n  %s\n", paste0(output_prefix, "_headline.csv"),
            paste0(output_prefix, "_long.csv")))
