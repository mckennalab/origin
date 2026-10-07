# Per-integration feature definitions for the clone 84 comparison.
#
# Shared by clone84_feature_preservation.R (10 integrations, real against
# simulated) and clone84_all_integrations.R (all 34 real barcodes). Kept in one
# file because a hand-copied metric that drifts between call sites is a failure
# this project has already hit more than once.
#
# Every function takes the matrix as cells x sites with NA for unobserved, and
# an `integration` vector naming each column's barcode.

# ---- features ----
site_rate <- function(x) {
  edited <- colSums(x == 1L, na.rm = TRUE)
  observed <- colSums(!is.na(x))
  ifelse(observed > 0, edited / observed, NA_real_)
}
cell_edits <- function(x) rowSums(x == 1L, na.rm = TRUE)

#' Distinct edit patterns per integration, and how evenly they are used
#'
#' Treating a missing entry as its own symbol makes the allele a description of
#' the ASSAY rather than the recorder: because all three datasets carry the same
#' mask, the commonest "allele" becomes "barcode not captured" and every dataset
#' returns the same answer. Alleles are therefore defined only over cells whose
#' barcode was captured, and only at the sites observed in all of them, so two
#' cells are the same allele when their edits agree on comparable ground.
allele_diversity <- function(x, integration, min_observed = 0.5) {
  do.call(rbind, lapply(unique(integration), function(unit) {
    block <- x[, integration == unit, drop = FALSE]
    captured <- rowMeans(!is.na(block)) >= min_observed
    if (sum(captured) < 10L) return(NULL)
    block <- block[captured, , drop = FALSE]
    comparable <- which(colSums(is.na(block)) == 0L)
    if (length(comparable) < 5L) return(NULL)
    block <- block[, comparable, drop = FALSE]
    keys <- apply(block, 1L, paste, collapse = "")
    counts <- table(keys)
    proportions <- counts / sum(counts)
    entropy <- -sum(proportions * log(proportions))
    data.frame(integration = unit, cells = nrow(block),
               sites = length(comparable), alleles = length(counts),
               entropy = entropy,
               normalised_entropy = entropy / log(nrow(block)),
               top_allele_fraction = max(proportions),
               stringsAsFactors = FALSE)
  }))
}

#' Four-gamete incompatibility and pairwise association, by stratum
pair_statistics <- function(x, integration, pairs, rng_seed) {
  informative <- which(vapply(seq_len(ncol(x)), function(column) {
    values <- x[, column]
    values <- values[!is.na(values)]
    length(values) >= 20L && any(values == 1L) && any(values == 0L)
  }, logical(1)))
  if (length(informative) < 4L) return(NULL)
  set.seed(rng_seed)
  same <- outer(integration[informative], integration[informative], "==")
  draw <- function(want_same, wanted) {
    candidates <- which(upper.tri(same) & (same == want_same), arr.ind = TRUE)
    if (!nrow(candidates)) return(NULL)
    candidates[sample.int(nrow(candidates),
                          min(wanted, nrow(candidates))), , drop = FALSE]
  }
  strata <- list(within = draw(TRUE, pairs), between = draw(FALSE, pairs))
  do.call(rbind, lapply(names(strata), function(stratum) {
    index <- strata[[stratum]]
    if (is.null(index)) return(NULL)
    first <- informative[index[, 1L]]
    second <- informative[index[, 2L]]
    incompatible <- logical(length(first))
    association <- numeric(length(first))
    for (position in seq_along(first)) {
      a <- x[, first[position]]
      b <- x[, second[position]]
      usable <- !is.na(a) & !is.na(b)
      if (sum(usable) < 20L) {
        incompatible[position] <- NA
        association[position] <- NA_real_
        next
      }
      a <- a[usable]
      b <- b[usable]
      # All four gametes present means no tree can explain the pair without a
      # repeated or reverted edit.
      incompatible[position] <- any(a == 0L & b == 0L) && any(a == 0L & b == 1L) &&
        any(a == 1L & b == 0L) && any(a == 1L & b == 1L)
      association[position] <- if (stats::sd(a) == 0 || stats::sd(b) == 0) {
        NA_real_
      } else stats::cor(a, b)
    }
    data.frame(stratum = stratum, pairs = length(first),
               incompatible_fraction = mean(incompatible, na.rm = TRUE),
               mean_abs_association = mean(abs(association), na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
}

#' Every feature, restricted to one integration
#'
#' The pooled view can hide a split result: an integration whose barcode was
#' poorly captured, or whose sites happen to be fast, can carry the summary.
#' Per-integration values show whether a mismatch is a property of the recorder
#' or of a few barcodes. Real barcodes and simulated integrations correspond by
#' position -- the i-th sorted barcode against int_i -- not by identity.
integration_features <- function(x, integration, pairs, rng_seed) {
  units <- unique(integration)
  do.call(rbind, lapply(seq_along(units), function(index) {
    unit <- units[index]
    columns <- which(integration == unit)
    block <- x[, columns, drop = FALSE]
    rates <- site_rate(block)
    edits <- rowSums(block == 1L, na.rm = TRUE)
    observed <- rowSums(!is.na(block))
    alleles <- allele_diversity(block, rep(unit, ncol(block)))
    within <- pair_statistics(block, rep(unit, ncol(block)), pairs,
                              rng_seed + index)
    within_row <- if (!is.null(within)) {
      within[within$stratum == "within", , drop = FALSE]
    } else NULL
    data.frame(
      integration = unit, integration_index = index, sites = ncol(block),
      observed_fraction = mean(!is.na(block)),
      site_rate_mean = mean(rates, na.rm = TRUE),
      site_rate_median = stats::median(rates, na.rm = TRUE),
      site_rate_saturated = mean(rates >= 0.99, na.rm = TRUE),
      edits_per_cell_mean = mean(edits),
      edited_fraction = mean(edits / pmax(observed, 1L)),
      alleles = if (!is.null(alleles)) alleles$alleles[1] else NA_integer_,
      normalised_entropy = if (!is.null(alleles)) {
        alleles$normalised_entropy[1]
      } else NA_real_,
      top_allele_fraction = if (!is.null(alleles)) {
        alleles$top_allele_fraction[1]
      } else NA_real_,
      incompatible_fraction = if (!is.null(within_row) && nrow(within_row)) {
        within_row$incompatible_fraction[1]
      } else NA_real_,
      mean_abs_association = if (!is.null(within_row) && nrow(within_row)) {
        within_row$mean_abs_association[1]
      } else NA_real_,
      stringsAsFactors = FALSE)
  }))
}

