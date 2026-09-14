#!/usr/bin/env Rscript

# Standard Seurat clustering of SymSim profiles from a simulated lineage.
#
# Runs the conventional workflow -- normalise, variable features, scale, PCA,
# shared-nearest-neighbour graph, Louvain clustering, UMAP -- on the terminal
# cells of a simulated tree, and reports the cluster assignment. When a tree is
# supplied it also asks whether the clusters recover the lineage, which is the
# question the simulation exists to pose.
#
# Defaults are adjusted for a small simulated matrix rather than copied from a
# 10x tutorial. On 64 cells and a few hundred genes the usual settings are
# wrong in specific ways: PCA cannot return more components than cells, the
# 2000-variable-feature default silently takes every gene, and a resolution
# tuned for thousands of cells over-splits. Each is handled below and reported.
#
# Arguments:
#   --counts=<path>     genes-by-cells CSV(.gz), cells in columns (required)
#   --tree=<path>       Newick lineage; enables the cluster-vs-lineage check
#   --output=<prefix>   output prefix. Default beside the counts file.
#   --resolution=<x>    Louvain resolution. Default 0.5.
#   --dims=<n>          principal components to use. Default 15.
#   --variable=<n>      variable features to select. Default 200.
#   --seed=<n>          seed. Default 1.
#
# Outputs <prefix>_clusters.csv, <prefix>_summary.txt and <prefix>_umap.png.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
counts_path <- value_after("--counts")
if (is.null(counts_path)) {
  stop("--counts is required (genes by cells, cells in columns).", call. = FALSE)
}
counts_path <- normalizePath(counts_path, mustWork = TRUE)
tree_path <- value_after("--tree")
output_prefix <- value_after(
  "--output", sub("[.]csv([.]gz)?$", "", counts_path)
)
resolution <- as.numeric(value_after("--resolution", "0.5"))
requested_dims <- as.integer(value_after("--dims", "15"))
variable_features <- as.integer(value_after("--variable", "200"))
seed <- as.integer(value_after("--seed", "1"))

if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) {
  cache_dir <- file.path(tempdir(), "xdg-cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = cache_dir)
}
for (package in c("Seurat", "ggplot2", "ragg")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

counts <- as.matrix(read.csv(
  if (grepl("[.]gz$", counts_path)) gzfile(counts_path) else counts_path,
  row.names = 1, check.names = FALSE
))
cat(sprintf("[seurat] counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))
if (ncol(counts) < 10L) {
  stop("Too few cells to cluster meaningfully.", call. = FALSE)
}

# PCA cannot return more components than there are cells, and Seurat's npcs
# default of 50 exceeds a 64-cell matrix. Cap it, and cap the requested
# dimensions to what PCA can actually produce.
max_components <- min(nrow(counts), ncol(counts)) - 1L
npcs <- min(50L, max_components)
dims <- seq_len(min(requested_dims, npcs))
if (requested_dims > npcs) {
  cat(sprintf("[seurat] requested %d dimensions but only %d are available; using %d\n",
              requested_dims, npcs, length(dims)))
}
# FindVariableFeatures silently returns everything when nfeatures exceeds the
# gene count, which is not a selection at all. Say so rather than imply one.
if (variable_features >= nrow(counts)) {
  cat(sprintf("[seurat] %d variable features requested but only %d genes exist; using all\n",
              variable_features, nrow(counts)))
  variable_features <- nrow(counts)
}

object <- CreateSeuratObject(counts = counts, project = "symsim",
                             min.cells = 0, min.features = 0)
object <- NormalizeData(object, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
object <- FindVariableFeatures(object, selection.method = "vst",
                               nfeatures = variable_features, verbose = FALSE)
object <- ScaleData(object, features = rownames(object), verbose = FALSE)
object <- RunPCA(object, features = VariableFeatures(object), npcs = npcs,
                 seed.use = seed, verbose = FALSE)
object <- FindNeighbors(object, dims = dims, verbose = FALSE)
object <- FindClusters(object, resolution = resolution, random.seed = seed,
                       verbose = FALSE)
# n.neighbors must stay below the cell count for UMAP to run at all.
object <- RunUMAP(object, dims = dims, seed.use = seed, verbose = FALSE,
                  n.neighbors = min(30L, ncol(counts) - 1L))

clusters <- Idents(object)
cluster_table <- data.frame(
  cell = names(clusters), cluster = as.character(clusters),
  stringsAsFactors = FALSE
)
cat(sprintf("[seurat] %d clusters at resolution %.2f: sizes %s\n",
            length(unique(clusters)), resolution,
            paste(as.integer(table(clusters)), collapse = ", ")))

summary_lines <- c(
  sprintf("cells: %d", ncol(counts)),
  sprintf("genes: %d", nrow(counts)),
  sprintf("variable features: %d", length(VariableFeatures(object))),
  sprintf("principal components used: %d of %d computed", length(dims), npcs),
  sprintf("resolution: %.2f", resolution),
  sprintf("clusters: %d", length(unique(clusters))),
  sprintf("cluster sizes: %s", paste(as.integer(table(clusters)), collapse = ", "))
)

# Do the clusters recover the lineage? Compared against the tree's own deep
# split rather than asserted by eye.
if (!is.null(tree_path) && file.exists(tree_path) &&
    requireNamespace("ape", quietly = TRUE)) {
  tree <- ape::read.tree(tree_path)
  shared <- intersect(tree$tip.label, cluster_table$cell)
  if (length(shared) >= 10L) {
    distance <- stats::cophenetic(ape::keep.tip(tree, shared))[shared, shared]
    # Cut the tree into the same number of groups the clustering found, so the
    # comparison is like for like rather than rewarding a different granularity.
    lineage_groups <- stats::cutree(
      stats::hclust(stats::as.dist(distance), method = "average"),
      k = length(unique(clusters))
    )
    assigned <- cluster_table$cluster[match(shared, cluster_table$cell)]
    contingency <- table(lineage = lineage_groups, cluster = assigned)
    agreement <- if (requireNamespace("mclust", quietly = TRUE)) {
      mclust::adjustedRandIndex(lineage_groups, assigned)
    } else {
      NA_real_
    }
    cluster_table$lineage_group <- NA_integer_
    cluster_table$lineage_group[match(shared, cluster_table$cell)] <-
      unname(lineage_groups)
    cat("\n[lineage] clusters against lineage groups of matched granularity:\n")
    print(contingency)
    cat(sprintf("[lineage] adjusted Rand index = %.3f\n", agreement))
    summary_lines <- c(
      summary_lines,
      sprintf("adjusted Rand index vs lineage: %.3f", agreement),
      "contingency (rows lineage, columns cluster):",
      utils::capture.output(print(contingency))
    )
  }
}

utils::write.csv(cluster_table, paste0(output_prefix, "_clusters.csv"),
                 row.names = FALSE)
writeLines(summary_lines, paste0(output_prefix, "_summary.txt"))

umap_plot <- DimPlot(object, reduction = "umap", label = TRUE, pt.size = 2.2) +
  ggplot2::labs(
    title = sprintf("SymSim terminal cells, %d clusters (resolution %.2f)",
                    length(unique(clusters)), resolution),
    subtitle = sprintf("%d cells, %d genes, %d principal components",
                       ncol(counts), nrow(counts), length(dims))
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.subtitle = ggplot2::element_text(size = 9, colour = "#4A4A4A")
  )
ggplot2::ggsave(paste0(output_prefix, "_umap.png"), plot = umap_plot,
                device = ragg::agg_png, width = 7, height = 6, units = "in",
                dpi = 220, background = "#FFFFFF")

cat("\nWrote:\n", paste0("  ", c(
  paste0(output_prefix, "_clusters.csv"),
  paste0(output_prefix, "_summary.txt"),
  paste0(output_prefix, "_umap.png")
), collapse = "\n"), "\n")
