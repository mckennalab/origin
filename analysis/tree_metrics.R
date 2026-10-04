# Tree-comparison metrics shared by the recorder assessments.
#
# Normalised RF answers one narrow question -- what fraction of unrooted splits
# match -- and weights every split equally. These functions add rooted
# structure, clade-level precision and recall, distance geometry, ancestry
# ordering, clone assignment, and whether errors fall on early or late splits.
#
# Kept in one file because three separate assessments consume them, and a
# hand-copied metric that drifts between call sites is exactly the failure this
# project has already hit with duplicated rate labels.
#
# Sampling counts and clone granularities are arguments rather than globals, so
# a caller cannot silently change another caller's definition of a metric.

# ---- tree helpers -----------------------------------------------------------

#' Descendant tip sets of every non-trivial clade of a rooted tree
#'
#' Root and single-tip clades are dropped: both trees always share them, so
#' counting them would inflate every agreement score toward one.
clade_keys <- function(tree) {
  parts <- ape::prop.part(tree)
  labels <- attr(parts, "labels")
  sets <- lapply(parts, function(indices) sort(labels[indices]))
  sizes <- lengths(sets)
  keep <- sizes >= 2L & sizes < length(labels)
  vapply(sets[keep], paste, character(1), collapse = "\r")
}

#' Midpoint-root a tree, tolerating the negative branch lengths NJ produces
root_estimate <- function(tree) {
  working <- tree
  if (!is.null(working$edge.length)) {
    working$edge.length <- pmax(working$edge.length, 0)
  }
  rooted <- tryCatch(phangorn::midpoint(working), error = function(e) NULL)
  if (is.null(rooted)) {
    # Fall back to rooting on the most distant tip, which is still a defined
    # rooting rather than silently comparing a rooted tree with an unrooted one.
    distances <- ape::cophenetic.phylo(working)
    outgroup <- rownames(distances)[which.max(rowSums(distances))]
    rooted <- tryCatch(ape::root(working, outgroup = outgroup, resolve.root = TRUE),
                       error = function(e) NULL)
  }
  rooted
}

#' Topological depth of every node, in edges from the root
node_depths <- function(tree) {
  unit <- tree
  unit$edge.length <- rep.int(1, nrow(tree$edge))
  ape::node.depth.edgelength(unit)
}

#' Rooted triplet agreement
#'
#' For three tips the rooted topology names one pair as the cherry. Agreement is
#' the fraction of sampled triples where both trees name the same cherry, which
#' is the rooted counterpart of a split comparison and is sensitive to where the
#' root sits.
triplet_agreement <- function(truth, estimate, tips, samples) {
  if (length(tips) < 3L) return(NA_real_)
  truth_mrca <- ape::mrca(truth)
  estimate_mrca <- ape::mrca(estimate)
  truth_depth <- node_depths(truth)
  estimate_depth <- node_depths(estimate)
  draws <- matrix(
    tips[sample.int(length(tips), samples * 3L, replace = TRUE)],
    ncol = 3L
  )
  draws <- draws[draws[, 1L] != draws[, 2L] & draws[, 2L] != draws[, 3L] &
                   draws[, 1L] != draws[, 3L], , drop = FALSE]
  if (!nrow(draws)) return(NA_real_)
  cherry_of <- function(mrca_matrix, depth) {
    ab <- depth[mrca_matrix[cbind(draws[, 1L], draws[, 2L])]]
    ac <- depth[mrca_matrix[cbind(draws[, 1L], draws[, 3L])]]
    bc <- depth[mrca_matrix[cbind(draws[, 2L], draws[, 3L])]]
    winner <- max.col(cbind(ab, ac, bc), ties.method = "first")
    # A tie means the triple is unresolved in this tree; mark it so it cannot
    # be scored as agreement by accident.
    spread <- pmax(ab, ac, bc) - pmin(ab, ac, bc)
    ifelse(spread > 0, winner, NA_integer_)
  }
  truth_cherry <- cherry_of(truth_mrca, truth_depth)
  estimate_cherry <- cherry_of(estimate_mrca, estimate_depth)
  usable <- !is.na(truth_cherry) & !is.na(estimate_cherry)
  if (!any(usable)) return(NA_real_)
  mean(truth_cherry[usable] == estimate_cherry[usable])
}

#' Ancestor-descendant ordering
#'
#' For two disjoint pairs of tips, which pair coalesced more recently? This asks
#' whether the estimate orders divergence events correctly, using topological
#' MRCA depth rather than branch lengths, so it is not a restatement of the
#' cophenetic correlation below.
ordering_agreement <- function(truth, estimate, tips, samples) {
  if (length(tips) < 4L) return(NA_real_)
  truth_mrca <- ape::mrca(truth)
  estimate_mrca <- ape::mrca(estimate)
  truth_depth <- node_depths(truth)
  estimate_depth <- node_depths(estimate)
  draws <- matrix(
    tips[sample.int(length(tips), samples * 4L, replace = TRUE)],
    ncol = 4L
  )
  distinct <- apply(draws, 1L, function(row) length(unique(row)) == 4L)
  draws <- draws[distinct, , drop = FALSE]
  if (!nrow(draws)) return(NA_real_)
  order_of <- function(mrca_matrix, depth) {
    first <- depth[mrca_matrix[cbind(draws[, 1L], draws[, 2L])]]
    second <- depth[mrca_matrix[cbind(draws[, 3L], draws[, 4L])]]
    sign(first - second)
  }
  truth_order <- order_of(truth_mrca, truth_depth)
  estimate_order <- order_of(estimate_mrca, estimate_depth)
  # A zero means the two pairs coalesced at equal depth, i.e. that tree does not
  # order them at all. Scoring those as disagreements drags a random tree below
  # the 0.5 chance level, so -- as with the triplet metric -- only comparisons
  # both trees actually resolve are counted.
  usable <- truth_order != 0 & estimate_order != 0
  if (!any(usable)) return(NA_real_)
  mean(truth_order[usable] == estimate_order[usable])
}

#' Clone assignment accuracy
#'
#' Both trees are cut into the same number of groups and compared by adjusted
#' Rand index. Matching the group count scores assignment rather than rewarding
#' a tree for choosing the right number of clones.
clone_accuracy <- function(truth, estimate, tips, groups) {
  cut_tree <- function(tree) {
    distances <- stats::as.dist(ape::cophenetic.phylo(tree)[tips, tips])
    stats::cutree(stats::hclust(distances, method = "average"), k = groups)
  }
  truth_groups <- tryCatch(cut_tree(truth), error = function(e) NULL)
  estimate_groups <- tryCatch(cut_tree(estimate), error = function(e) NULL)
  if (is.null(truth_groups) || is.null(estimate_groups)) return(NA_real_)
  mclust::adjustedRandIndex(truth_groups, estimate_groups)
}

#' Recall of true clades split by how early they arose
#'
#' The truth tree is ultrametric, so a clade's root-to-node time is literally
#' when that split happened. Errors concentrated on late splits mean a recorder
#' resolves deep structure but not recent structure, which is a different
#' failure from uniform error and is invisible in a single RF number.
recall_by_era <- function(truth, estimate_keys) {
  parts <- ape::prop.part(truth)
  labels <- attr(parts, "labels")
  sets <- lapply(parts, function(indices) sort(labels[indices]))
  sizes <- lengths(sets)
  keep <- which(sizes >= 2L & sizes < length(labels))
  if (!length(keep)) return(NULL)
  keys <- vapply(sets[keep], paste, character(1), collapse = "\r")
  # prop.part returns one entry per internal node, in node order starting at
  # the root, so node id is the tip count plus the entry index.
  node_ids <- length(labels) + keep
  times <- ape::node.depth.edgelength(truth)[node_ids]
  era <- cut(times, breaks = stats::quantile(times, probs = seq(0, 1, 0.25)),
             include.lowest = TRUE,
             labels = c("earliest", "early", "late", "latest"))
  recovered <- keys %in% estimate_keys
  data.frame(
    era = factor(era, levels = c("earliest", "early", "late", "latest")),
    recovered = recovered, stringsAsFactors = FALSE
  )
}

#' Every metric for one truth/estimate pair
score_pair <- function(truth, estimate, tips, triplet_samples = 20000L,
                       ordering_samples = 20000L,
                       clone_counts = c(5L, 10L, 25L)) {
  truth_keys <- clade_keys(truth)
  estimate_keys <- clade_keys(estimate)
  shared <- length(intersect(truth_keys, estimate_keys))
  truth_distance <- ape::cophenetic.phylo(truth)[tips, tips]
  estimate_distance <- ape::cophenetic.phylo(estimate)[tips, tips]
  upper <- upper.tri(truth_distance)
  result <- list(
    unrooted_rf = tryCatch(
      phangorn::RF.dist(ape::unroot(truth), ape::unroot(estimate),
                        normalize = TRUE),
      error = function(e) NA_real_),
    rooted_rf = 1 - (2 * shared) / (length(truth_keys) + length(estimate_keys)),
    clade_precision = if (length(estimate_keys)) shared / length(estimate_keys)
      else NA_real_,
    clade_recall = if (length(truth_keys)) shared / length(truth_keys)
      else NA_real_,
    triplet_agreement = triplet_agreement(truth, estimate, tips,
                                          triplet_samples),
    ordering_agreement = ordering_agreement(truth, estimate, tips,
                                            ordering_samples),
    cophenetic_spearman = suppressWarnings(stats::cor(
      truth_distance[upper], estimate_distance[upper], method = "spearman"))
  )
  for (groups in clone_counts) {
    result[[sprintf("clone_ari_k%d", groups)]] <-
      clone_accuracy(truth, estimate, tips, groups)
  }
  result$era <- recall_by_era(truth, estimate_keys)
  result
}

