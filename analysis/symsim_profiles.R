# SymSim single-cell expression profiles over a simulated lineage.
#
# Generates one expression profile per terminal (sampled) cell of a lineage
# produced by the Gillespie or PhysiCell engines, while the internal nodes stay
# latent: they are never observed, but they carry the expression state that
# descendants inherit.
#
# SymSim (https://github.com/YosefLab/SymSim) is a GitHub-only package and is
# deliberately not a dependency of the origin package. Source this file and call
# the functions directly, the same way physicell_scdesign3.R is used.
#
# HOW SYMSIM IS USED, AND WHERE IT DOES NOT FIT
#
# SymSim's continuous mode walks a set of extrinsic variation factors (EVFs)
# down a tree as a Brownian motion, so cells sharing ancestry share expression
# state. That is exactly the model wanted here, and the lineage tree can be
# handed to it directly as `phyla`.
#
# What it does not do is place exactly one cell at each tip. SampleEdge()
# allocates cells to an edge in proportion to branch length and recurses with a
# shared budget, so some terminal edges receive no cell at all. Measured on a
# 12-tip tree: 10 of 12 terminal edges covered at 200 cells, 11 of 12 at 600.
# The edges missed were of ordinary length, so this is the sampler's budget and
# not a short-branch artifact -- oversampling improves coverage but never
# guarantees it.
#
# So assignment is explicit rather than assumed. Each terminal node takes the
# deepest cell sampled on its own terminal edge; if that edge received none, it
# takes the deepest cell from the nearest ancestral edge that did. Every profile
# is labelled with which of the two happened, and the coverage fraction is
# reported, so a caller can tell how much of the result is exact and how much is
# an ancestor standing in.
#
# Required: SymSim, ape. Optional: Matrix, for sparse output.

#' Check that SymSim is installed
#'
#' @return `TRUE` invisibly, or an error naming the install command.
require_symsim <- function() {
  if (!requireNamespace("SymSim", quietly = TRUE)) {
    stop(
      "The SymSim package is required and is not installed. Install it with:\n",
      '  remotes::install_github("YosefLab/SymSim")',
      call. = FALSE
    )
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("The ape package is required.", call. = FALSE)
  }
  # SymSim::SimulateTrueCounts() reads bundled parameter tables (match_params
  # and friends) by bare name rather than from its own namespace, so calling it
  # through SymSim:: fails with "object 'match_params' not found". The package
  # has to be attached, not merely loaded.
  if (!"package:SymSim" %in% search()) {
    suppressPackageStartupMessages(
      library("SymSim", character.only = TRUE)
    )
  }
  invisible(TRUE)
}

#' Build the SymSim input tree from a simulated lineage
#'
#' Converts an engine lineage table into the `phylo` object SymSim takes as
#' `phyla`. Only cells alive at the horizon become tips; every division point
#' between them stays an internal node and is never sampled.
#'
#' @param nodes Lineage node table, as returned in `population$nodes` by
#'   `simulate_gillespie_population()` or written to `lineage_nodes.csv`.
#' @param terminal_ids Optional PhysiCell ids to use as tips. Defaults to the
#'   cells flagged `alive_at_end`.
#' @param min_branch_length Replacement for zero or missing branch lengths.
#'   SymSim allocates cells in proportion to branch length and divides by the
#'   total, so a zero-length branch contributes nothing and can make the
#'   allocation degenerate.
#' @return A named list with `tree` (a `phylo`), `terminal_ids` in tip order,
#'   and `tip_edges`, the `parent_child` label SymSim uses for each tip's
#'   terminal edge.
#' @section Side effects: None.
prepare_symsim_tree <- function(nodes,
                                terminal_ids = NULL,
                                min_branch_length = 1e-6) {
  require_symsim()
  if (!all(c("physicell_id", "alive_at_end") %in% names(nodes))) {
    stop("nodes must carry physicell_id and alive_at_end columns.", call. = FALSE)
  }
  if (is.null(terminal_ids)) {
    terminal_ids <- nodes$physicell_id[
      nodes$alive_at_end %in% c(TRUE, "TRUE")
    ]
  }
  terminal_ids <- unique(as.character(terminal_ids))
  if (length(terminal_ids) < 3L) {
    stop("At least three terminal cells are needed to build a tree.",
         call. = FALSE)
  }
  if (!exists("physicell_lineage_to_newick", mode = "function")) {
    stop("physicell_lineage_to_newick() not found; load the origin package ",
         "or source load_origin.R first.", call. = FALSE)
  }

  newick <- physicell_lineage_to_newick(
    nodes, terminal_physicell_ids = terminal_ids
  )
  tree <- ape::read.tree(text = newick)
  if (is.null(tree)) stop("The lineage did not parse as a tree.", call. = FALSE)

  # Restricting to a subset of terminal cells leaves behind internal nodes with
  # a single child. SymSim's SampleSubtree() indexes children positionally and
  # fails on those with "missing value where TRUE/FALSE needed", so collapse
  # them. This preserves topology and summed branch lengths; only the unary
  # waypoints disappear, and they carry no branching information.
  singles_before <- tree$Nnode
  tree <- ape::collapse.singles(tree)
  collapsed_nodes <- singles_before - tree$Nnode

  # SymSim divides by the summed branch length, so non-positive branches have to
  # go before it sees the tree.
  invalid <- !is.finite(tree$edge.length) | tree$edge.length <= 0
  tree$edge.length[invalid] <- min_branch_length

  tip_indices <- tree$edge[, 2] <= ape::Ntip(tree)
  tip_edges <- stats::setNames(
    paste(tree$edge[tip_indices, 1], tree$edge[tip_indices, 2], sep = "_"),
    tree$tip.label[tree$edge[tip_indices, 2]]
  )
  list(
    tree = tree,
    terminal_ids = tree$tip.label,
    tip_edges = tip_edges,
    replaced_branches = sum(invalid),
    collapsed_single_child_nodes = collapsed_nodes
  )
}

#' Map each tip to one sampled SymSim cell
#'
#' Prefers a cell sampled on the tip's own terminal edge, taking the deepest one
#' because that is the point closest to the tip. Falls back to the nearest
#' ancestral edge when SymSim's budget left the terminal edge empty.
#'
#' @param cell_meta SymSim's `cell_meta`, carrying `pop` and `depth`.
#' @param prepared Output of `prepare_symsim_tree()`.
#' @return Data frame with one row per tip: `terminal_id`, the selected
#'   `cell_index`, and `assignment` of either `terminal_edge` or
#'   `ancestor_edge`. Tips with no usable cell get `NA` and
#'   `assignment = "unassigned"`.
assign_symsim_cells <- function(cell_meta, prepared) {
  tree <- prepared$tree
  n_tip <- ape::Ntip(tree)
  parents <- stats::setNames(tree$edge[, 1], tree$edge[, 2])
  edge_label <- function(child) {
    parent <- parents[[as.character(child)]]
    if (is.null(parent) || is.na(parent)) return(NA_character_)
    paste(parent, child, sep = "_")
  }
  by_pop <- split(seq_len(nrow(cell_meta)), cell_meta$pop)
  deepest <- vapply(by_pop, function(rows) {
    rows[which.max(cell_meta$depth[rows])]
  }, integer(1))

  # `deepest[[label]]` throws when the name is absent, and absent is the normal
  # case here: SymSim leaves some edges unsampled. Test membership first.
  has_cell <- function(label) !is.na(label) && label %in% names(deepest)
  rows <- lapply(seq_len(n_tip), function(tip) {
    label <- edge_label(tip)
    if (has_cell(label)) {
      return(data.frame(
        terminal_id = tree$tip.label[tip],
        cell_index = unname(deepest[[label]]),
        assignment = "terminal_edge", stringsAsFactors = FALSE
      ))
    }
    # Walk toward the root until an edge with a sampled cell turns up.
    current <- tip
    repeat {
      parent <- parents[[as.character(current)]]
      if (is.null(parent) || is.na(parent)) break
      label <- edge_label(parent)
      if (has_cell(label)) {
        return(data.frame(
          terminal_id = tree$tip.label[tip],
          cell_index = unname(deepest[[label]]),
          assignment = "ancestor_edge", stringsAsFactors = FALSE
        ))
      }
      current <- parent
    }
    data.frame(terminal_id = tree$tip.label[tip], cell_index = NA_integer_,
               assignment = "unassigned", stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Simulate SymSim expression profiles over a lineage
#'
#' Runs SymSim's continuous (tree) mode on the lineage, then reduces the sampled
#' cells to one profile per terminal node.
#'
#' @param prepared Output of `prepare_symsim_tree()`.
#' @param ngenes Number of genes. SymSim fails with an opaque error about array
#'   dimensions below roughly 100, so this is validated.
#' @param oversample Cells requested per tree edge. SymSim's coverage of
#'   terminal edges improves with this and is reported in the result.
#' @param nevf,n_de_evf,sigma,vary SymSim EVF settings, passed through.
#' @param seed Random seed.
#' @param observed Whether to also draw observed counts with
#'   `True2ObservedCounts()`, which adds capture efficiency and sequencing
#'   depth on top of the true counts.
#' @param protocol,alpha_mean,depth_mean Observed-count settings, used only when
#'   `observed` is `TRUE`.
#' @return Named list with `true_counts` (genes by terminal cells, columns named
#'   by terminal id), optional `observed_counts`, the `assignment` table, the
#'   full SymSim `simulation`, and a `coverage` summary.
#' @section Side effects: Calls `set.seed()` through SymSim.
simulate_symsim_profiles <- function(prepared,
                                     ngenes = 200L,
                                     oversample = 4L,
                                     nevf = 10L,
                                     n_de_evf = 6L,
                                     sigma = 0.4,
                                     vary = "s",
                                     seed = 1L,
                                     observed = FALSE,
                                     protocol = "UMI",
                                     alpha_mean = 0.1,
                                     depth_mean = 45000) {
  require_symsim()
  ngenes <- as.integer(ngenes)
  if (!is.finite(ngenes) || ngenes < 100L) {
    stop("ngenes must be at least 100; SymSim fails below roughly that size ",
         "with an opaque error about array dimensions.", call. = FALSE)
  }
  tree <- prepared$tree
  requested_cells <- max(
    as.integer(oversample) * nrow(tree$edge), 3L * ape::Ntip(tree)
  )

  simulation <- SymSim::SimulateTrueCounts(
    ncells_total = requested_cells,
    min_popsize = max(5L, as.integer(ape::Ntip(tree) / 10L)),
    ngenes = ngenes,
    evf_type = "continuous",
    nevf = nevf,
    n_de_evf = n_de_evf,
    vary = vary,
    Sigma = sigma,
    phyla = tree,
    randseed = seed
  )

  assignment <- assign_symsim_cells(simulation$cell_meta, prepared)
  usable <- !is.na(assignment$cell_index)
  if (!any(usable)) {
    stop("SymSim produced no cell that could be assigned to any terminal ",
         "node; raise oversample.", call. = FALSE)
  }

  true_counts <- simulation$counts[, assignment$cell_index[usable], drop = FALSE]
  colnames(true_counts) <- assignment$terminal_id[usable]
  rownames(true_counts) <- paste0("gene_", seq_len(nrow(true_counts)))

  observed_counts <- NULL
  if (isTRUE(observed)) {
    gene_lengths <- sample(SymSim::gene_len_pool, nrow(true_counts),
                           replace = FALSE)
    observed <- SymSim::True2ObservedCounts(
      true_counts = true_counts,
      meta_cell = simulation$cell_meta[assignment$cell_index[usable], ,
                                       drop = FALSE],
      protocol = protocol, alpha_mean = alpha_mean, alpha_sd = 0.002,
      gene_len = gene_lengths, depth_mean = depth_mean, depth_sd = 4500
    )
    observed_counts <- observed$counts
    dimnames(observed_counts) <- dimnames(true_counts)
  }

  coverage <- data.frame(
    terminal_nodes = nrow(assignment),
    on_terminal_edge = sum(assignment$assignment == "terminal_edge"),
    from_ancestor_edge = sum(assignment$assignment == "ancestor_edge"),
    unassigned = sum(assignment$assignment == "unassigned"),
    cells_simulated = ncol(simulation$counts),
    cells_requested = requested_cells,
    stringsAsFactors = FALSE
  )
  coverage$terminal_edge_fraction <-
    coverage$on_terminal_edge / coverage$terminal_nodes

  list(
    true_counts = true_counts,
    observed_counts = observed_counts,
    assignment = assignment,
    coverage = coverage,
    simulation = simulation,
    tree = tree
  )
}

#' Write SymSim profiles and their provenance
#'
#' @param profiles Output of `simulate_symsim_profiles()`.
#' @param output_dir Destination directory, created if absent.
#' @param compress Whether to gzip the CSV tables.
#' @return The output directory, invisibly.
#' @section Side effects: Writes `symsim_true_counts.csv`, optionally
#'   `symsim_observed_counts.csv`, plus `symsim_cell_assignment.csv`,
#'   `symsim_coverage.csv` and `symsim_lineage.nwk` under `output_dir`.
write_symsim_outputs <- function(profiles, output_dir, compress = TRUE) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_table <- function(value, name, row_names = FALSE) {
    path <- file.path(output_dir, name)
    if (compress) {
      path <- paste0(path, ".gz")
      connection <- gzfile(path, open = "wt")
      utils::write.csv(value, connection, row.names = row_names)
      close(connection)
    } else {
      utils::write.csv(value, path, row.names = row_names)
    }
    path
  }
  write_table(as.data.frame(profiles$true_counts), "symsim_true_counts.csv",
              row_names = TRUE)
  if (!is.null(profiles$observed_counts)) {
    write_table(as.data.frame(profiles$observed_counts),
                "symsim_observed_counts.csv", row_names = TRUE)
  }
  write_table(profiles$assignment, "symsim_cell_assignment.csv")
  write_table(profiles$coverage, "symsim_coverage.csv")
  ape::write.tree(profiles$tree, file.path(output_dir, "symsim_lineage.nwk"))
  invisible(output_dir)
}
