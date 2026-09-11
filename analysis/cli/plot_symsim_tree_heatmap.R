#!/usr/bin/env Rscript

# Simulate a small lineage, give every terminal cell a SymSim expression
# profile, and draw the tree beside a terminal-cell correlation matrix whose
# rows and columns follow the tree's tip order.
#
# The figure is the check that the two halves agree: if expression carries
# lineage structure, the high-correlation blocks on the diagonal should line up
# with the clades drawn to their left. Nothing in the plot enforces that -- the
# heatmap is ordered by the tree, not by clustering the expression -- so any
# block structure that appears is a property of the simulation.
#
# Arguments:
#   --terminals=<n>    terminal cells to keep. Default 64.
#   --genes=<n>        genes to simulate, at least 100. Default 400.
#   --oversample=<n>   SymSim cells requested per tree edge. Default 100, which
#                      is high on purpose: see the note below.
#   --end-time=<x>     simulation horizon. Default 5.
#   --seed=<n>         seed for the lineage and the profiles. Default 4.
#   --output=<prefix>  output path prefix. Default analysis/figures/symsim_tree_heatmap
#   --method=<m>       correlation method, pearson or spearman. Default pearson.
#   --allow-duplicate-profiles  plot even if some tips share a profile.
#   --sampling=<mode>  clades (default) or uniform. See below.
#
# On --sampling: a pure-birth lineage sampled uniformly at one horizon is
# star-like. Measured on the default run, 55% of tip pairs sat at the maximum
# cophenetic distance, sibling pairs correlated at 0.693 against 0.680 for
# unrelated pairs, and the tree-expression association was only -0.147. There is
# almost nothing for the heatmap to show, because SymSim accumulates EVF
# variance along branches and tips that coalesce near the root are close to
# independent. Drawing the cells from a few deep clades instead gives the tips
# shared internal branches: the same measurement returns -0.725, with close
# pairs at 0.711 against 0.583 for distant ones. Both modes are honest; clades
# is the default because it is what makes the structure legible, and it matches
# how lineage-tracing experiments usually sample.
#
# On --oversample: SymSim does not guarantee a cell on every terminal edge, and
# a tip with no cell of its own borrows one from an ancestral edge. Two tips can
# then borrow the SAME cell and end up with identical profiles, which shows up
# as a spurious correlation of exactly 1. At 64 tips, oversample 8 left 9 tips
# borrowing and 4 profiles duplicated; oversample 40 gave all 64 tips their own
# cell and no duplicates. The script reports both counts and warns if any
# duplicate survives, because that would be an artifact of the assignment rather
# than of the biology.

arguments <- commandArgs(trailingOnly = TRUE)
script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

value_after <- function(prefix, default) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
terminals <- as.integer(value_after("--terminals", "64"))
genes <- as.integer(value_after("--genes", "400"))
oversample <- as.integer(value_after("--oversample", "100"))
end_time <- as.numeric(value_after("--end-time", "5"))
seed <- as.integer(value_after("--seed", "4"))
method <- tolower(value_after("--method", "pearson"))
sampling <- tolower(value_after("--sampling", "clades"))
if (!sampling %in% c("clades", "uniform")) {
  stop("--sampling must be clades or uniform.", call. = FALSE)
}
output_prefix <- value_after(
  "--output", file.path(repo_root, "analysis", "figures", "symsim_tree_heatmap")
)
if (!method %in% c("pearson", "spearman")) {
  stop("--method must be pearson or spearman.", call. = FALSE)
}

if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) {
  cache_dir <- file.path(tempdir(), "xdg-cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = cache_dir)
}
for (package in c("ggplot2", "ggtree", "aplot", "ragg", "ape")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required for this figure.", package),
         call. = FALSE)
  }
}

origin_include_analysis <- FALSE
source(file.path(repo_root, "load_origin.R"))
source(file.path(repo_root, "analysis", "symsim_profiles.R"))
suppressMessages(library(ggplot2))

# One cycling cell type, no differentiation: the only structure in the tree is
# the branching itself, so the heatmap cannot be explained by cell type.
params <- list(
  num_init_cells = 1, sim_length = list(end_time), random_seed = seed,
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
  params, end_time = end_time, seed = seed, max_cells = 100000L,
  show_progress = FALSE
)
alive <- population$nodes$physicell_id[
  population$nodes$alive_at_end %in% c(TRUE, "TRUE")
]
cat(sprintf("[tree] %d nodes, %d cells alive at the horizon\n",
            nrow(population$nodes), length(alive)))
if (length(alive) < terminals) {
  stop(sprintf("Only %d terminal cells available; asked for %d. Raise ",
               length(alive), terminals), "--end-time.", call. = FALSE)
}
#' Choose the terminal cells that become tips
#'
#' Uniform sampling spreads tips over the whole population, which for a
#' pure-birth lineage means most of them coalesce at the root. Clade sampling
#' takes equal groups from a few disjoint subtrees so the tips share deep
#' internal branches.
#'
#' @param nodes Lineage node table.
#' @param alive_ids Candidate terminal ids.
#' @param count Number of tips wanted.
#' @param mode Either "clades" or "uniform".
#' @param groups Number of clades to draw from.
#' @return Character vector of terminal ids.
choose_terminals <- function(nodes, alive_ids, count, mode, groups = 4L) {
  if (identical(mode, "uniform") || count < 2L * groups) {
    return(sample(alive_ids, count))
  }
  full <- prepare_symsim_tree(nodes, terminal_ids = alive_ids)
  tree <- full$tree
  per_clade <- count %/% groups
  internal <- (ape::Ntip(tree) + 1L):(ape::Ntip(tree) + tree$Nnode)
  sizes <- vapply(internal, function(node) {
    length(ape::extract.clade(tree, node)$tip.label)
  }, integer(1))
  # Cap clade size as well as flooring it. Without an upper bound a "clade" can
  # be most of the tree, the groups stop being distinct, and the structure the
  # mode exists to create disappears: allowing up to 80% of the population gave
  # an association of -0.217, while tight clades give about -0.7.
  candidates <- internal[sizes >= per_clade & sizes <= 3L * per_clade]
  if (!length(candidates)) {
    warning("No clade large enough; falling back to uniform sampling.",
            call. = FALSE)
    return(sample(alive_ids, count))
  }
  picked <- character()
  used <- character()
  for (node in sample(candidates)) {
    tips <- ape::extract.clade(tree, node)$tip.label
    if (length(intersect(tips, used))) next
    picked <- c(picked, sample(tips, per_clade))
    used <- c(used, tips)
    if (length(picked) >= count) break
  }
  if (length(picked) < count) {
    warning(sprintf(
      "Only %d tips from disjoint clades; topping up uniformly.", length(picked)
    ), call. = FALSE)
    remaining <- setdiff(alive_ids, sub("^[^0-9]*", "", picked))
    picked <- c(picked, sample(remaining, count - length(picked)))
  }
  # Tree tip labels carry a cell_ prefix; the lineage table keys on the bare id.
  sub("^[^0-9]*", "", picked[seq_len(count)])
}

set.seed(seed + 1L)
terminal_ids <- choose_terminals(population$nodes, alive, terminals, sampling)
cat(sprintf("[tree] sampling mode: %s\n", sampling))
prepared <- prepare_symsim_tree(population$nodes, terminal_ids = terminal_ids)
cat(sprintf("[tree] %d tips, %d latent internal nodes\n",
            ape::Ntip(prepared$tree), prepared$tree$Nnode))

profiles <- suppressWarnings(simulate_symsim_profiles(
  prepared, ngenes = genes, oversample = oversample, seed = seed
))
coverage <- profiles$coverage
allow_duplicates <- any(arguments == "--allow-duplicate-profiles")
duplicated_cells <- sum(duplicated(
  profiles$assignment$cell_index[!is.na(profiles$assignment$cell_index)]
))
cat(sprintf("[profiles] %d genes, %d cells; %d tips exact, %d from an ancestor\n",
            nrow(profiles$true_counts), ncol(profiles$true_counts),
            coverage$on_terminal_edge, coverage$from_ancestor_edge))
if (duplicated_cells > 0L && !allow_duplicates) {
  # This is not a cosmetic problem. Tips that borrow the same ancestor cell get
  # byte-identical profiles and appear as a solid block of r = 1 that looks
  # exactly like a real clade signal. Refuse rather than emit a figure whose
  # most eye-catching feature is an artifact.
  duplicate_ids <- profiles$assignment$terminal_id[
    duplicated(profiles$assignment$cell_index) |
      duplicated(profiles$assignment$cell_index, fromLast = TRUE)
  ]
  stop(sprintf(
    paste0(
      "%d terminal cells share a profile with another tip (%s). They would ",
      "render as a block of correlation 1 that is indistinguishable from a ",
      "real clade. Raise --oversample (currently %d) until every tip gets its ",
      "own cell, or pass --allow-duplicate-profiles to plot anyway."
    ),
    duplicated_cells, paste(utils::head(duplicate_ids, 6), collapse = ", "),
    oversample
  ), call. = FALSE)
}

expression <- log1p(profiles$true_counts)
correlation <- stats::cor(expression, method = method)

# Order strictly by the tree. aplot aligns on the shared discrete axis, so the
# tip labels have to be the factor levels and the level order has to be the
# order ggtree draws them in.
tree_plot <- ggtree::ggtree(profiles$tree, linewidth = 0.35) +
  ggtree::geom_tippoint(size = 0.9, colour = "#B4436C") +
  ggplot2::labs(title = sprintf("Lineage (%d terminal cells)",
                                ape::Ntip(profiles$tree))) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    plot.margin = ggplot2::margin(5, 2, 5, 5)
  )

# Take the tip order from the rendered tree and impose it on BOTH axes. aplot
# aligns only the shared axis (y); left to itself the x axis keeps ggplot's
# alphabetical order, and the two orders disagree. The visible symptom is the
# self-correlations scattering off the diagonal instead of forming it, which is
# also the giveaway that any apparent block structure would be meaningless.
tree_data <- tree_plot$data
tip_rows <- tree_data[tree_data$isTip, ]
tip_order <- tip_rows$label[order(tip_rows$y)]
if (!setequal(tip_order, rownames(correlation))) {
  stop("Tree tips and correlation matrix do not carry the same cells.",
       call. = FALSE)
}

long <- data.frame(
  row_cell = factor(rep(rownames(correlation), times = ncol(correlation)),
                    levels = tip_order),
  column_cell = factor(rep(colnames(correlation), each = nrow(correlation)),
                       levels = tip_order),
  correlation = as.vector(correlation),
  stringsAsFactors = FALSE
)
heatmap_plot <- ggplot2::ggplot(
  long, ggplot2::aes(x = column_cell, y = row_cell, fill = correlation)
) +
  ggplot2::geom_raster() +
  ggplot2::scale_fill_gradient2(
    name = sprintf("%s r", tools::toTitleCase(method)),
    low = "#3D67A8", mid = "#F7F7F5", high = "#B4436C",
    midpoint = stats::median(correlation[upper.tri(correlation)]),
    limits = range(correlation)
  ) +
  ggplot2::scale_x_discrete(expand = c(0, 0)) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::labs(
    title = "Terminal-cell expression correlation",
    subtitle = "Rows and columns follow the tree; the matrix is not clustered",
    x = NULL, y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 11),
    plot.subtitle = ggplot2::element_text(size = 8.5, colour = "#4A4A4A"),
    legend.position = "right",
    legend.key.width = ggplot2::unit(0.35, "cm")
  )

combined <- aplot::insert_left(heatmap_plot, tree_plot, width = 0.42)

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 11, height = 7.5, units = "in", dpi = 220,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 11, height = 7.5, units = "in", bg = "#FFFFFF")

# Does the block structure actually follow the tree? Compare correlation with
# cophenetic distance rather than leaving it to the eye.
tree_distance <- stats::cophenetic(profiles$tree)[
  rownames(correlation), colnames(correlation)
]
upper <- upper.tri(tree_distance)
association <- stats::cor(tree_distance[upper], correlation[upper],
                          method = "spearman")
counts_path <- paste0(output_prefix, "_correlation.csv.gz")
connection <- gzfile(counts_path, open = "wt")
utils::write.csv(correlation, connection)
close(connection)

cat(sprintf(
  "[check] Spearman(cophenetic distance, expression correlation) = %.3f (negative means closer relatives correlate more)\n",
  association
))
cat("\nWrote:\n", paste0("  ", c(png_path, pdf_path, counts_path),
                          collapse = "\n"), "\n")
