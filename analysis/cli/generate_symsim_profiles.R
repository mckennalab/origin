#!/usr/bin/env Rscript

# Generate SymSim single-cell expression profiles for the terminal cells of a
# simulated lineage, keeping the internal nodes latent.
#
# Reads a lineage node table written by either engine -- lineage_nodes.csv(.gz)
# from the Gillespie pipeline, or the equivalent from a PhysiCell replay -- and
# writes one expression profile per terminal cell.
#
# Arguments:
#   --lineage-nodes=<path>   lineage_nodes.csv or .csv.gz (required)
#   --output-dir=<path>      destination directory (required)
#   --terminal-cells=<path>  optional terminal_cells.csv(.gz); its physicell_id
#                            column selects the tips. Defaults to every cell
#                            flagged alive_at_end in the node table.
#   --max-terminals=<n>      randomly subsample to at most this many tips.
#                            SymSim's cost grows with the tree, so this is the
#                            knob for large populations. Default 200.
#   --genes=<n>              number of genes, at least 100. Default 500.
#   --oversample=<n>         cells requested per tree edge. Higher values raise
#                            the share of tips assigned from their own terminal
#                            edge rather than an ancestor. Default 6.
#   --nevf=<n>               SymSim EVFs. Default 10.
#   --n-de-evf=<n>           differential EVFs. Default 6.
#   --sigma=<x>              EVF random-walk standard deviation. Default 0.4.
#   --seed=<n>               random seed. Default 1.
#   --observed               also draw observed counts (capture efficiency and
#                            sequencing depth) as well as true counts.
#   --depth-mean=<x>         mean read depth for observed counts. Default 45000.
#
# Outputs, under --output-dir:
#   symsim_true_counts.csv.gz        genes by terminal cells
#   symsim_observed_counts.csv.gz    only with --observed
#   symsim_cell_assignment.csv.gz    which SymSim cell each tip took, and whether
#                                    it came from the tip's own edge
#   symsim_coverage.csv.gz           assignment counts and the exact fraction
#   symsim_lineage.nwk               the tree handed to SymSim
#
# Prints SYMSIM_OUTPUT_DIR=<path> on success.

arguments <- commandArgs(trailingOnly = TRUE)
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
# This script lives in analysis/cli, two levels below the repository root.
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
has_flag <- function(flag) any(arguments == flag)

usage <- function() {
  cat(readLines(script_path)[3:38], sep = "\n")
}
if (has_flag("--help") || has_flag("-h")) {
  usage()
  quit(save = "no", status = 0)
}

lineage_path <- value_after("--lineage-nodes")
output_dir <- value_after("--output-dir")
if (is.null(lineage_path) || is.null(output_dir)) {
  usage()
  stop("--lineage-nodes and --output-dir are both required.", call. = FALSE)
}
lineage_path <- normalizePath(lineage_path, mustWork = TRUE)
terminal_path <- value_after("--terminal-cells")
max_terminals <- as.integer(value_after("--max-terminals", "200"))
genes <- as.integer(value_after("--genes", "500"))
oversample <- as.integer(value_after("--oversample", "6"))
nevf <- as.integer(value_after("--nevf", "10"))
n_de_evf <- as.integer(value_after("--n-de-evf", "6"))
sigma <- as.numeric(value_after("--sigma", "0.4"))
seed <- as.integer(value_after("--seed", "1"))
depth_mean <- as.numeric(value_after("--depth-mean", "45000"))
observed <- has_flag("--observed")

if (any(!is.finite(c(max_terminals, genes, oversample, nevf, n_de_evf, seed))) ||
    max_terminals < 3L || oversample < 1L || nevf < 1L || n_de_evf < 0L) {
  stop("Terminal count, gene count, oversample, or EVF settings are invalid.",
       call. = FALSE)
}
if (genes < 100L) {
  stop("--genes must be at least 100; SymSim fails below roughly that size ",
       "with an opaque error about array dimensions.", call. = FALSE)
}

origin_include_analysis <- FALSE
source(file.path(repo_root, "load_origin.R"))
source(file.path(repo_root, "analysis", "symsim_profiles.R"))

read_maybe_gz <- function(path) {
  if (grepl("[.]gz$", path)) read.csv(gzfile(path)) else read.csv(path)
}
nodes <- read_maybe_gz(lineage_path)

terminal_ids <- NULL
if (!is.null(terminal_path)) {
  terminal_table <- read_maybe_gz(normalizePath(terminal_path, mustWork = TRUE))
  if (!"physicell_id" %in% names(terminal_table)) {
    stop("--terminal-cells table has no physicell_id column.", call. = FALSE)
  }
  terminal_ids <- as.character(terminal_table$physicell_id)
} else {
  terminal_ids <- as.character(
    nodes$physicell_id[nodes$alive_at_end %in% c(TRUE, "TRUE")]
  )
}
terminal_ids <- unique(terminal_ids)
cat(sprintf("[symsim] %d nodes, %d terminal cells\n",
            nrow(nodes), length(terminal_ids)))

if (length(terminal_ids) > max_terminals) {
  set.seed(seed)
  terminal_ids <- sample(terminal_ids, max_terminals)
  cat(sprintf("[symsim] subsampled to %d terminal cells\n", max_terminals))
}

prepared <- prepare_symsim_tree(nodes, terminal_ids = terminal_ids)
cat(sprintf(
  "[symsim] tree: %d tips, %d latent internal nodes (%d unary collapsed)\n",
  ape::Ntip(prepared$tree), prepared$tree$Nnode,
  prepared$collapsed_single_child_nodes
))

profiles <- simulate_symsim_profiles(
  prepared, ngenes = genes, oversample = oversample, nevf = nevf,
  n_de_evf = n_de_evf, sigma = sigma, seed = seed,
  observed = observed, depth_mean = depth_mean
)

coverage <- profiles$coverage
cat(sprintf(
  "[symsim] %d of %d tips took a cell from their own terminal edge (%.0f%%); %d fell back to an ancestor\n",
  coverage$on_terminal_edge, coverage$terminal_nodes,
  100 * coverage$terminal_edge_fraction, coverage$from_ancestor_edge
))
if (coverage$from_ancestor_edge > 0) {
  cat("[symsim] Ancestor fallbacks share expression with a parent branch ",
      "rather than being drawn at the tip. Raise --oversample to reduce them.\n",
      sep = "")
}
if (coverage$unassigned > 0) {
  warning(sprintf("%d tips could not be assigned any cell.",
                  coverage$unassigned), call. = FALSE)
}

if (!grepl("^/", output_dir)) output_dir <- file.path(repo_root, output_dir)
write_symsim_outputs(profiles, output_dir)
cat(sprintf("[symsim] counts: %d genes x %d cells\n",
            nrow(profiles$true_counts), ncol(profiles$true_counts)))
cat(sprintf("SYMSIM_OUTPUT_DIR=%s\n", normalizePath(output_dir)))
