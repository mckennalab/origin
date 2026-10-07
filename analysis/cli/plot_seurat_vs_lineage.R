#!/usr/bin/env Rscript

# Seurat clusters against the lineage that generated them.
#
# Runs the standard workflow at several Louvain resolutions and draws, for each,
# the lineage with its tips coloured by cluster alongside the UMAP. Putting the
# tree beside the embedding is the point: the UMAP shows whether clusters are
# separable, the tree shows whether they correspond to clades, and those are
# different questions.
#
# Arguments:
#   --counts=<path>      genes-by-cells CSV(.gz) (required)
#   --tree=<path>        Newick lineage (required)
#   --resolutions=<list> comma-separated Louvain resolutions. Default 0.5,1.2,2
#   --dims=<n>           principal components. Default 15.
#   --variable=<n>       variable features. Default 200.
#   --seed=<n>           seed. Default 1.
#   --output=<prefix>    output prefix.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
counts_path <- value_after("--counts")
tree_path <- value_after("--tree")
if (is.null(counts_path) || is.null(tree_path)) {
  stop("--counts and --tree are both required.", call. = FALSE)
}
counts_path <- normalizePath(counts_path, mustWork = TRUE)
tree_path <- normalizePath(tree_path, mustWork = TRUE)
resolutions <- as.numeric(strsplit(
  value_after("--resolutions", "0.5,1.2,2"), ",", fixed = TRUE
)[[1L]])
requested_dims <- as.integer(value_after("--dims", "15"))
variable_features <- as.integer(value_after("--variable", "200"))
seed <- as.integer(value_after("--seed", "1"))
output_prefix <- value_after(
  "--output", file.path(dirname(counts_path), "seurat_vs_lineage")
)

if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) {
  cache_dir <- file.path(tempdir(), "xdg-cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(XDG_CACHE_HOME = cache_dir)
}
for (package in c("Seurat", "ggplot2", "ggtree", "patchwork", "ape", "ragg",
                  "mclust")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("The %s package is required.", package), call. = FALSE)
  }
}
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  # ggtree must be attached, not merely loaded: %<+% is an infix operator and
  # cannot be reached through ggtree::.
  library(ggtree)
})

counts <- as.matrix(read.csv(
  if (grepl("[.]gz$", counts_path)) gzfile(counts_path) else counts_path,
  row.names = 1, check.names = FALSE
))
tree <- ape::read.tree(tree_path)
cat(sprintf("[input] %d genes x %d cells; tree with %d tips\n",
            nrow(counts), ncol(counts), ape::Ntip(tree)))

max_components <- min(nrow(counts), ncol(counts)) - 1L
npcs <- min(50L, max_components)
dims <- seq_len(min(requested_dims, npcs))
if (variable_features >= nrow(counts)) variable_features <- nrow(counts)

object <- CreateSeuratObject(counts = counts, project = "symsim",
                             min.cells = 0, min.features = 0)
object <- NormalizeData(object, verbose = FALSE)
object <- FindVariableFeatures(object, selection.method = "vst",
                               nfeatures = variable_features, verbose = FALSE)
object <- ScaleData(object, features = rownames(object), verbose = FALSE)
object <- suppressWarnings(RunPCA(object, features = VariableFeatures(object),
                                  npcs = npcs, seed.use = seed, verbose = FALSE))
object <- FindNeighbors(object, dims = dims, verbose = FALSE)
object <- suppressWarnings(RunUMAP(
  object, dims = dims, seed.use = seed, verbose = FALSE,
  n.neighbors = min(30L, ncol(counts) - 1L)
))
embedding <- as.data.frame(Embeddings(object, "umap"))
names(embedding) <- c("umap_1", "umap_2")
embedding$cell <- rownames(embedding)

# A fixed palette so the same cluster id keeps its colour across panels.
cluster_palette <- c("#B4436C", "#087E8B", "#D29A18", "#3D67A8", "#708B36",
                     "#DC6B2F", "#8E6C9B", "#5A5A5A", "#B8B03A", "#2F6F5E")

tree_panels <- list()
umap_panels <- list()
summary_rows <- list()

for (index in seq_along(resolutions)) {
  resolution <- resolutions[[index]]
  object <- FindClusters(object, resolution = resolution, random.seed = seed,
                         verbose = FALSE)
  clusters <- Idents(object)
  n_clusters <- length(unique(clusters))

  # Cut the tree to the same number of groups the clustering found, so the
  # comparison rewards correct assignment rather than matching granularity.
  shared <- intersect(tree$tip.label, names(clusters))
  distance <- stats::cophenetic(ape::keep.tip(tree, shared))[shared, shared]
  lineage_groups <- stats::cutree(
    stats::hclust(stats::as.dist(distance), method = "average"), k = n_clusters
  )
  assigned <- as.character(clusters[shared])
  agreement <- mclust::adjustedRandIndex(lineage_groups, assigned)
  cat(sprintf("[resolution %.2f] %d clusters (sizes %s), ARI vs lineage %.3f\n",
              resolution, n_clusters,
              paste(as.integer(table(clusters)), collapse = ", "), agreement))
  summary_rows[[index]] <- data.frame(
    resolution = resolution, clusters = n_clusters,
    adjusted_rand_index = agreement, stringsAsFactors = FALSE
  )

  tip_data <- data.frame(
    label = shared, cluster = factor(assigned,
                                     levels = sort(unique(assigned))),
    stringsAsFactors = FALSE
  )
  values <- cluster_palette[seq_len(nlevels(tip_data$cluster))]
  names(values) <- levels(tip_data$cluster)

  tree_panels[[index]] <- (ggtree::ggtree(tree, linewidth = 0.3) %<+% tip_data) +
    ggtree::geom_tippoint(ggplot2::aes(colour = cluster), size = 1.9,
                          na.rm = TRUE) +
    ggplot2::scale_colour_manual(name = "cluster", values = values,
                                 na.value = "grey70") +
    ggplot2::labs(
      title = sprintf("resolution %.2f: %d clusters", resolution, n_clusters),
      subtitle = sprintf("ARI vs lineage %.3f", agreement)
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.subtitle = ggplot2::element_text(size = 8.5, colour = "#4A4A4A")
    )

  panel_data <- embedding
  panel_data$cluster <- factor(as.character(clusters[panel_data$cell]),
                               levels = levels(tip_data$cluster))
  umap_panels[[index]] <- ggplot2::ggplot(
    panel_data, ggplot2::aes(umap_1, umap_2, colour = cluster)
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(name = "cluster", values = values) +
    ggplot2::labs(title = sprintf("UMAP, resolution %.2f", resolution),
                  x = "UMAP 1", y = "UMAP 2") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 10)
    )
}

summary_table <- do.call(rbind, summary_rows)
utils::write.csv(summary_table, paste0(output_prefix, "_summary.csv"),
                 row.names = FALSE)

combined <- patchwork::wrap_plots(
  patchwork::wrap_plots(tree_panels, nrow = 1),
  patchwork::wrap_plots(umap_panels, nrow = 1),
  ncol = 1, heights = c(1.15, 1)
) +
  patchwork::plot_annotation(
    title = "Seurat clusters against the lineage that produced the cells",
    subtitle = sprintf(
      "%d terminal cells, %d genes, %d principal components. Top: lineage with tips coloured by cluster. Bottom: the same clusters on one shared UMAP.",
      ncol(counts), nrow(counts), length(dims)
    ),
    caption = paste(
      "Lineage groups for the ARI are obtained by cutting the tree to the same",
      "number of groups the clustering found, so agreement reflects assignment\n",
      "rather than matching granularity. The UMAP embedding is computed once",
      "and shared across panels; only the colouring changes."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 9.5, colour = "#4A4A4A"),
      plot.caption = ggplot2::element_text(size = 8, colour = "#4A4A4A",
                                           hjust = 0)
    )
  )

png_path <- paste0(output_prefix, ".png")
pdf_path <- paste0(output_prefix, ".pdf")
ggplot2::ggsave(png_path, plot = combined, device = ragg::agg_png,
                width = 13, height = 9, units = "in", dpi = 220,
                background = "#FFFFFF")
ggplot2::ggsave(pdf_path, plot = combined, device = grDevices::cairo_pdf,
                width = 13, height = 9, units = "in", bg = "#FFFFFF")

print(summary_table, row.names = FALSE, digits = 3)
cat("\nWrote:\n", paste0("  ", c(png_path, pdf_path,
                                 paste0(output_prefix, "_summary.csv")),
                          collapse = "\n"), "\n")
