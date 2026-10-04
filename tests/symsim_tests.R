#!/usr/bin/env Rscript

# Tests for the SymSim lineage-profile integration in analysis/symsim_profiles.R.
#
# Kept separate from regression_tests.R because SymSim is a GitHub-only package
# that most checkouts will not have. The whole file skips cleanly when it is
# absent, so this can be run unconditionally:
#
#   Rscript tests/symsim_tests.R

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE),
                        value = TRUE)
script_path <- sub("^--file=", "", script_argument[[1L]])
repo_root <- normalizePath(file.path(dirname(script_path), ".."),
                           mustWork = TRUE)

if (!requireNamespace("SymSim", quietly = TRUE)) {
  cat("SymSim is not installed; skipping SymSim integration tests.\n")
  cat('Install with: remotes::install_github("YosefLab/SymSim")\n')
  quit(save = "no", status = 0)
}

failures <- 0L
expect_true <- function(condition, label) {
  if (!isTRUE(condition)) {
    failures <<- failures + 1L
    cat(sprintf("FAIL %s\n", label))
  } else {
    cat(sprintf("ok   %s\n", label))
  }
}
expect_error <- function(expression, pattern, label) {
  message_text <- tryCatch({
    force(expression)
    NA_character_
  }, error = function(condition) conditionMessage(condition))
  expect_true(!is.na(message_text) && grepl(pattern, message_text), label)
}

origin_include_analysis <- FALSE
source(file.path(repo_root, "load_origin.R"))
source(file.path(repo_root, "analysis", "symsim_profiles.R"))
suppressMessages(library(ape))

lineage_params <- list(
  num_init_cells = 1, sim_length = list(4), random_seed = 1,
  editing_induction = list(timepoint = 0, num_cells = NULL, frac_cells = 1),
  differentiation_induction = list(timepoint = 0, num_cells = NULL,
                                   frac_cells = 1),
  cell_type_dict = list(
    founder_cell_type = "p",
    cell_type_params = list(
      p = list(cell_cycle_length = 1, death_per_cell_cycle_prob = 0.02)
    ),
    uninduced_transition_matrix = list(list(1)),
    induced_transition_matrix = list(list(1))
  )
)
population <- simulate_gillespie_population(
  lineage_params, end_time = 4, seed = 1, max_cells = 100000L,
  show_progress = FALSE
)
alive <- population$nodes$physicell_id[
  population$nodes$alive_at_end %in% c(TRUE, "TRUE")
]
set.seed(11)
# 120 rather than a smaller number on purpose. The lineage-structure check
# below is a correlation over tip pairs, and at 45 tips it has so little power
# that it fails roughly one run in three on perfectly good signal (measured
# p-values of 0.030, 0.069 and 0.020 across three seeds). At 120 tips the
# correlation is 0.316-0.320 with p at the permutation floor across every seed
# tried, so the test detects a real effect instead of coin-flipping.
selected <- sample(alive, 120)

prepared <- prepare_symsim_tree(population$nodes, terminal_ids = selected)

expect_true(inherits(prepared$tree, "phylo"),
            "prepare_symsim_tree returns a phylo")
expect_true(ape::Ntip(prepared$tree) == length(selected),
            "one tip per requested terminal cell")
expect_true(all(prepared$tree$edge.length > 0),
            "no zero-length branches survive")
# Subsetting terminals leaves unary chains that SymSim cannot recurse through.
expect_true(prepared$collapsed_single_child_nodes > 0,
            "unary internal nodes were collapsed")
expect_true(prepared$tree$Nnode < ape::Ntip(prepared$tree),
            "internal nodes are fewer than tips after collapsing")
# The internal nodes must still exist: they are the hidden ancestors whose
# state descendants inherit.
expect_true(prepared$tree$Nnode > 1,
            "latent internal nodes are retained")

expect_error(
  prepare_symsim_tree(population$nodes, terminal_ids = head(alive, 2)),
  "At least three terminal cells",
  "too few terminal cells is rejected"
)

profiles <- suppressWarnings(simulate_symsim_profiles(
  prepared, ngenes = 300, oversample = 8, seed = 1
))

expect_true(nrow(profiles$true_counts) == 300,
            "gene count matches the request")
expect_true(ncol(profiles$true_counts) == ape::Ntip(prepared$tree),
            "one profile per terminal cell")
expect_true(all(colnames(profiles$true_counts) %in% prepared$tree$tip.label),
            "profile columns are named by terminal id")
expect_true(!anyDuplicated(colnames(profiles$true_counts)),
            "each terminal cell appears once")
expect_true(all(profiles$true_counts >= 0),
            "counts are non-negative")
expect_true(profiles$coverage$unassigned == 0,
            "every terminal cell received a profile")
expect_true(profiles$coverage$on_terminal_edge > 0,
            "some tips are assigned from their own terminal edge")
expect_true(
  profiles$coverage$on_terminal_edge + profiles$coverage$from_ancestor_edge +
    profiles$coverage$unassigned == profiles$coverage$terminal_nodes,
  "the coverage breakdown sums to the terminal count"
)

expect_error(
  simulate_symsim_profiles(prepared, ngenes = 50),
  "at least 100",
  "an unusably small gene count is rejected with a clear message"
)

# The point of using SymSim's tree mode at all: expression has to carry lineage
# structure. Cells closer in the tree should be closer in expression. Compared
# against a permuted null, because the raw correlation alone does not establish
# that the tree is what produced it.
exact <- profiles$assignment$terminal_id[
  profiles$assignment$assignment == "terminal_edge"
]
if (length(exact) >= 10L) {
  expression <- log1p(profiles$true_counts[, exact, drop = FALSE])
  subtree <- ape::keep.tip(profiles$tree, exact)
  tree_distance <- stats::cophenetic(subtree)[exact, exact]
  expression_distance <- as.matrix(stats::dist(t(expression)))
  upper <- upper.tri(tree_distance)
  observed <- stats::cor(tree_distance[upper], expression_distance[upper],
                         method = "spearman")
  set.seed(7)
  null_values <- replicate(100, {
    shuffled <- sample(exact)
    stats::cor(
      tree_distance[upper],
      as.matrix(stats::dist(t(log1p(
        profiles$true_counts[, shuffled, drop = FALSE]
      ))))[upper],
      method = "spearman"
    )
  })
  # A permutation p-value, not a comparison against max(null): the maximum of
  # 100 draws is an extreme order statistic, so requiring the observed value to
  # exceed it is a far stricter test than the nominal 0.05 and fails on
  # perfectly good signal.
  permutation_p <- (1 + sum(null_values >= observed)) / (1 + length(null_values))
  cat(sprintf(
    "     lineage-expression Spearman %.3f, permuted null %.3f (max %.3f), p = %.3f\n",
    observed, mean(null_values), max(null_values), permutation_p
  ))
  expect_true(permutation_p < 0.05,
              "expression distance tracks lineage distance beyond the permuted null")
} else {
  cat("     too few exact assignments to test lineage structure; skipped\n")
}

output_dir <- file.path(tempdir(), "symsim-tests")
unlink(output_dir, recursive = TRUE)
write_symsim_outputs(profiles, output_dir)
for (name in c("symsim_true_counts.csv.gz", "symsim_cell_assignment.csv.gz",
               "symsim_coverage.csv.gz", "symsim_lineage.nwk")) {
  expect_true(file.exists(file.path(output_dir, name)),
              paste("wrote", name))
}
written_tree <- ape::read.tree(file.path(output_dir, "symsim_lineage.nwk"))
expect_true(ape::Ntip(written_tree) == ape::Ntip(prepared$tree),
            "the written tree keeps every tip")
unlink(output_dir, recursive = TRUE)

if (failures > 0L) {
  cat(sprintf("\n%d SymSim test(s) failed.\n", failures))
  quit(save = "no", status = 1L)
}
cat("\nAll SymSim integration tests passed.\n")
