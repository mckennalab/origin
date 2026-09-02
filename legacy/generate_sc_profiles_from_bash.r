#!/usr/bin/env Rscript
# generate_sc_profiles_from_bash.r
#
# Post-processing step that uses scDesign3 to generate single-cell expression
# profiles for every cell alive at every stopping point of a completed sim run
# (i.e. internal cells at all non-terminal stopping points + terminal cells at
# the final stopping point). Cell type comes straight from the simulator;
# pseudotime is derived from lineage-string depth and is optional.
#
# Required external R packages:
#   scDesign3, SingleCellExperiment, SummarizedExperiment, optparse, rjson, digest
# Optional:
#   Seurat       (only if --reference is a Seurat .rds)
#   zellkonverter (only if --reference is a .h5ad)
#
# Install (Bioconductor):
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("scDesign3", "SingleCellExperiment", "zellkonverter"))
#   install.packages(c("optparse", "rjson", "digest"))
#
# Invoked by bash_wrapper_all_combos.sh from the repository root as:
#   Rscript generate_sc_profiles_from_bash.r -I <run_id> -R <reference> \
#           --celltype_col <col> [--celltype_map <json>] [--use_pseudotime] \
#           [--ncores <n>] [--max_cells_per_timepoint <n>]
# The wrapper always passes -I, -R and --celltype_col, and appends each optional
# flag only when the matching wrapper option was supplied. -I and -R are
# required; -O, --cache_dir, --seed and --pseudotime_col are never passed by the
# wrapper and take the defaults declared in option_list below.
#
# Inputs (relative paths resolve against the working directory, which the
# wrapper sets to the repo root):
#   output/cell_populations/<run_id>/cell_population_*_time_<t>.rds
#   the reference dataset named by -R (.rds SCE/Seurat, or .h5ad)
#   the JSON named by --celltype_map, when supplied
#   <cache_dir>/fit_<key>.rds, reused when a compatible cached fit is present
#   scdesign3_helpers.R, sourced from this script's own directory
#
# Outputs, written under -O (default output/sc_profiles/<run_id>/):
#   sim_sce_all.rds         simulated counts + metadata for every cell
#   sim_sce_time_<t>.rds    one subset per stopping point
#   cell_metadata.csv       the covariate table used for the simulation
#   cell_count_summary.csv  cell counts by timepoint, cell type and terminal
#                           flag
# plus <cache_dir>/fit_<key>.rds (default output/scdesign3_fits/) whenever a new
# scDesign3 fit has to be computed.

script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- normalizePath(
  sub('^--file=', '', script_argument),
  mustWork = TRUE
)
# This script lives in legacy/, so the repository root is one level up. Its
# scDesign3 helpers moved to the analysis tier when the code was split.
repo_root <- normalizePath(file.path(dirname(script_path), '..'))
source(file.path(repo_root, 'analysis', 'scdesign3_helpers.R'))

suppressPackageStartupMessages({
  library(optparse)
  library(rjson)
})

option_list <- list(
  make_option(c('-I', '--run_id'), type = 'character', default = NULL,
              help = 'Simulator run id (the 13-digit name under output/cell_populations/)'),
  make_option(c('-R', '--reference'), type = 'character', default = NULL,
              help = 'Path to reference SCE / Seurat .rds or AnnData .h5ad'),
  make_option(c('-C', '--celltype_col'), type = 'character', default = 'cell_type',
              help = 'colData column in the reference holding cell-type labels [default: %default]'),
  make_option(c('--celltype_map'), type = 'character', default = NULL,
              help = 'Optional JSON mapping simulator cell types (ct1, ct2, ...) to reference labels'),
  make_option(c('-O', '--output_dir'), type = 'character', default = NULL,
              help = 'Output directory [default: output/sc_profiles/<run_id>/]'),
  make_option(c('--cache_dir'), type = 'character', default = 'output/scdesign3_fits',
              help = 'Where to cache scDesign3 fits keyed by reference content [default: %default]'),
  make_option(c('--ncores'), type = 'integer', default = 4,
              help = 'Cores for scDesign3 fitting / simulation [default: %default]'),
  make_option(c('--seed'), type = 'integer', default = 1,
              help = 'Random seed [default: %default]'),
  make_option(c('--use_pseudotime'), action = 'store_true', default = FALSE,
              help = 'Include pseudotime (lineage depth) as a smooth covariate; reference must also have a pseudotime column.'),
  make_option(c('--pseudotime_col'), type = 'character', default = 'pseudotime',
              help = 'Reference colData column for pseudotime when --use_pseudotime [default: %default]'),
  make_option(c('--max_cells_per_timepoint'), type = 'integer', default = NA,
              help = 'Optional cap on simulated cells per timepoint [default: no cap]')
)

opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

if (is.null(opt$run_id) || is.null(opt$reference)) {
  stop('--run_id and --reference are required.')
}

set.seed(opt$seed)

# ---- I/O paths ----------------------------------------------------------------
cell_pop_dir <- file.path('output', 'cell_populations', opt$run_id)
if (!dir.exists(cell_pop_dir)) {
  stop(sprintf('cell-population directory not found: %s', cell_pop_dir))
}

if (is.null(opt$output_dir)) {
  opt$output_dir <- file.path('output', 'sc_profiles', opt$run_id)
}
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$cache_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Reference loading --------------------------------------------------------
#' Load a reference dataset as a standardized SingleCellExperiment
#'
#' Accepts an AnnData `.h5ad` (via zellkonverter) or an `.rds` holding either a
#' `SingleCellExperiment` or a Seurat object, and normalizes it for scDesign3:
#' the first assay is aliased to `counts` if no `counts` assay exists,
#' `colData$cell_type` is set to the requested column as a factor, and
#' `colData$pseudotime` is set as numeric when pseudotime is requested.
#'
#' @param path Path to the reference `.rds` or `.h5ad`.
#' @param celltype_col Name of the `colData` column holding cell-type labels;
#'   must exist in the reference.
#' @param pseudotime_col Name of the `colData` column holding pseudotime; only
#'   required when `use_pseudotime` is TRUE.
#' @param use_pseudotime Logical; when TRUE the pseudotime column must be
#'   present and is copied to `colData$pseudotime`.
#' @return A `SingleCellExperiment` with a `counts` assay and standardized
#'   `cell_type` (and `pseudotime`) columns. Unsupported file extensions, object
#'   classes, missing columns, or an assay-free reference raise an error.
load_reference <- function(path, celltype_col, pseudotime_col, use_pseudotime) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SummarizedExperiment)
  })
  ext <- tolower(tools::file_ext(path))
  if (ext == 'h5ad') {
    if (!requireNamespace('zellkonverter', quietly = TRUE)) {
      stop('zellkonverter required for .h5ad input. BiocManager::install("zellkonverter")')
    }
    sce <- zellkonverter::readH5AD(path)
  } else if (ext == 'rds') {
    obj <- readRDS(path)
    if (is(obj, 'SingleCellExperiment')) {
      sce <- obj
    } else if (inherits(obj, 'Seurat')) {
      if (!requireNamespace('Seurat', quietly = TRUE)) {
        stop('Seurat required to ingest a Seurat .rds. install.packages("Seurat")')
      }
      sce <- Seurat::as.SingleCellExperiment(obj)
    } else {
      stop(sprintf('Unsupported reference R object class: %s',
                   paste(class(obj), collapse = ',')))
    }
  } else {
    stop(sprintf('Unrecognised reference extension: .%s (expected .rds or .h5ad)', ext))
  }

  cd <- as.data.frame(colData(sce))
  if (!(celltype_col %in% colnames(cd))) {
    stop(sprintf('Reference colData has no column "%s". Available: %s',
                 celltype_col, paste(colnames(cd), collapse = ', ')))
  }
  if (use_pseudotime && !(pseudotime_col %in% colnames(cd))) {
    stop(sprintf('--use_pseudotime requires reference colData column "%s"; not found.',
                 pseudotime_col))
  }

  if (!('counts' %in% assayNames(sce))) {
    if (length(assayNames(sce)) == 0) stop('Reference has no assays.')
    fallback <- assayNames(sce)[1]
    message(sprintf('Aliasing assay "%s" -> "counts" (no "counts" assay present).', fallback))
    assay(sce, 'counts') <- assay(sce, fallback)
  }
  colData(sce)$cell_type <- factor(as.character(colData(sce)[[celltype_col]]))
  if (use_pseudotime) {
    colData(sce)$pseudotime <- as.numeric(colData(sce)[[pseudotime_col]])
  }
  sce
}

# ---- Walk cell_population_*.rds files and build per-cell metadata -------------
#' List a run's cell-population snapshots in timepoint order
#'
#' @param cell_pop_dir Directory holding the simulator's
#'   `cell_population_<savename>_time_<t>.rds` snapshots for one run.
#' @return A data frame with `path` (full path) and `timepoint` (numeric, parsed
#'   from the file name), sorted by increasing timepoint. Errors when the
#'   directory contains no matching snapshot.
list_cell_population_files <- function(cell_pop_dir) {
  fs <- list.files(cell_pop_dir,
                   pattern = 'cell_population_.+_time_[0-9.]+\\.rds$',
                   full.names = TRUE)
  if (length(fs) == 0) {
    stop(sprintf('No cell_population_*_time_*.rds files in %s', cell_pop_dir))
  }
  tps <- as.numeric(sub('.*_time_([0-9.]+)\\.rds$', '\\1', basename(fs)))
  ord <- order(tps)
  data.frame(path = fs[ord], timepoint = tps[ord], stringsAsFactors = FALSE)
}

#' Build one metadata row per extant leaf cell at every stopping point
#'
#' Reads each snapshot and keeps the cells that are both `alive` and `terminal`,
#' i.e. the leaves of the lineage at that stopping point. Lineage depth is the
#' number of `_`-separated segments in the cell's lineage string minus one, so a
#' founder has depth 0. Because a lineage can still be a leaf at several
#' stopping points, `cell_id` repeats across rows and `sample_id` is the unique
#' per-observation identifier used for matrix column names.
#'
#' @param cellpop_files Data frame of snapshot `path`/`timepoint` pairs, as
#'   returned by `list_cell_population_files()`.
#' @param celltype_map Optional named list mapping simulator cell types
#'   (`ct1`, `ct2`, ...) to reference labels. When given, `cell_type` holds the
#'   mapped label and `cell_type_sim` preserves the simulator's own.
#' @return A data frame with one row per living leaf cell per stopping point:
#'   `cell_id` (lineage string), `cell_type`, `lineage_depth`, `birth_time`,
#'   `induced_editing`, `timepoint`, `is_terminal` (TRUE for rows from the last
#'   stopping point, not the simulator's per-cell terminal flag), `sample_id`,
#'   and `cell_type_sim` when a map was applied. Errors when no snapshot holds a
#'   living leaf cell.
build_covariate_df <- function(cellpop_files, celltype_map = NULL) {
  per_tp <- lapply(seq_len(nrow(cellpop_files)), function(i) {
    cp <- readRDS(cellpop_files$path[i])
    sampled <- vapply(
      cp,
      function(c) isTRUE(c$alive) && isTRUE(c$terminal),
      logical(1)
    )
    cp <- cp[sampled]
    if (length(cp) == 0) return(NULL)
    data.frame(
      cell_id = vapply(cp, function(c) c$linstring, character(1)),
      cell_type = vapply(cp, function(c) c$celltype, character(1)),
      lineage_depth = vapply(cp, function(c) {
        ls <- c$linstring
        if (is.null(ls) || is.na(ls) || ls == '') return(0L)
        length(strsplit(ls, '_', fixed = TRUE)[[1]]) - 1L
      }, integer(1)),
      birth_time = vapply(cp, function(c) {
        if (is.null(c$birth_time)) NA_real_ else as.numeric(c$birth_time)
      }, numeric(1)),
      induced_editing = vapply(cp, function(c) {
        ie <- c$induced_editing
        if (is.null(ie)) NA_character_ else as.character(ie)
      }, character(1)),
      timepoint = cellpop_files$timepoint[i],
      stringsAsFactors = FALSE
    )
  })
  meta <- do.call(rbind, per_tp[!vapply(per_tp, is.null, logical(1))])
  if (is.null(meta) || nrow(meta) == 0) {
    stop('No alive cells found across any stopping point.')
  }
  meta$is_terminal <- meta$timepoint == max(cellpop_files$timepoint)
  if (!is.null(celltype_map)) {
    keys <- names(celltype_map)
    new_ct <- meta$cell_type
    for (k in keys) new_ct[meta$cell_type == k] <- as.character(celltype_map[[k]])
    meta$cell_type_sim <- meta$cell_type
    meta$cell_type <- new_ct
  }
  # A lineage can be alive at several stopping points, so cell_id alone is not
  # unique in the combined output. sample_id identifies one cell-timepoint
  # observation and is safe for matrix/SingleCellExperiment column names.
  meta$sample_id <- make.unique(paste0(meta$cell_id, '_time_', meta$timepoint))
  meta
}

#' Add lineage depth rescaled to [0, 1] as a pseudotime covariate
#'
#' @param meta Cell metadata carrying a `lineage_depth` column.
#' @return `meta` with a `pseudotime` column holding depth divided by the
#'   maximum depth, or 0 for every row when the maximum depth is 0.
add_pseudotime <- function(meta) {
  max_depth <- max(meta$lineage_depth, na.rm = TRUE)
  meta$pseudotime <- if (max_depth == 0) 0 else meta$lineage_depth / max_depth
  meta
}

#' Cap the number of simulated cells per stopping point
#'
#' @param meta Cell metadata with a `timepoint` column.
#' @param max_per_tp Maximum rows to keep per timepoint, or NA for no cap.
#' @return `meta` unchanged when `max_per_tp` is NA, otherwise the rows
#'   regrouped by timepoint with each over-sized group reduced to a random
#'   subset of `max_per_tp` rows. Row order follows the split by timepoint, and
#'   the sampling depends on the seed set at script start.
maybe_downsample <- function(meta, max_per_tp) {
  if (is.na(max_per_tp)) return(meta)
  do.call(rbind, lapply(split(meta, meta$timepoint), function(sub) {
    if (nrow(sub) <= max_per_tp) return(sub)
    sub[sample(nrow(sub), max_per_tp), , drop = FALSE]
  }))
}

# ---- scDesign3 fit (cached) and simulate --------------------------------------
#' Compute the cache key identifying an scDesign3 fit of this reference
#'
#' Thin wrapper over `scdesign3_fit_cache_key()` from `scdesign3_helpers.R`,
#' which keys on the reference's resolved path, modification time and size, the
#' installed scDesign3 version, and the covariate choices below.
#'
#' @param ref_path Path to the reference dataset; must exist.
#' @param celltype_col Reference cell-type column used for the fit.
#' @param use_pseudotime Logical; whether pseudotime enters the model.
#' @param pseudotime_col Reference pseudotime column.
#' @return A character key (SHA-1 when `digest` is installed, otherwise a
#'   deterministic fallback) used to name the cached fit file.
ref_cache_key <- function(ref_path, celltype_col, use_pseudotime, pseudotime_col) {
  scdesign3_fit_cache_key(
    ref_path,
    celltype_col = celltype_col,
    pseudotime_col = pseudotime_col,
    use_pseudotime = use_pseudotime
  )
}

#' Reuse a cached scDesign3 fit, or fit the reference and cache it
#'
#' Thin wrapper over `fit_or_load_scdesign3()` from `scdesign3_helpers.R`.
#'
#' @param sce Standardized reference `SingleCellExperiment` from
#'   `load_reference()`.
#' @param cache_path Path of the cached fit `.rds` for this reference and
#'   covariate choice.
#' @param use_pseudotime Logical; adds a smooth pseudotime term to the default
#'   marginal-mean formula.
#' @param ncores Positive integer count of workers for marginal and copula
#'   fitting.
#' @return The fit list (`sce`, `data`, `marginal_list`, `copula_fit`,
#'   `mu_formula`, `family_use`, and the covariate flags).
#' @section Side effects: When no compatible cache exists, the helper fits the
#'   model and writes it to `cache_path`, creating that directory if needed.
fit_or_load <- function(sce, cache_path, use_pseudotime, ncores) {
  fit_or_load_scdesign3(
    sce,
    cache_path,
    use_pseudotime = use_pseudotime,
    ncores = ncores
  )
}

#' Simulate counts for the simulator-derived covariates
#'
#' Thin wrapper over `simulate_scdesign3_counts()` from `scdesign3_helpers.R`.
#'
#' @param fit Fit list returned by `fit_or_load()`.
#' @param new_meta Target cell metadata; must carry unique, non-missing
#'   `sample_id` values and cell types the fit knows about.
#' @param ncores Positive integer count of workers for parameter extraction and
#'   simulation.
#' @return A counts matrix with one row per reference feature and one column per
#'   row of `new_meta`, named by `sample_id`.
simulate_for_meta <- function(fit, new_meta, ncores) {
  simulate_scdesign3_counts(fit, new_meta, ncores = ncores)
}

# ---- Main --------------------------------------------------------------------
celltype_map <- if (!is.null(opt$celltype_map)) rjson::fromJSON(file = opt$celltype_map) else NULL

cellpop_files <- list_cell_population_files(cell_pop_dir)
meta <- build_covariate_df(cellpop_files, celltype_map = celltype_map)
meta <- add_pseudotime(meta)
meta <- maybe_downsample(meta, opt$max_cells_per_timepoint)
message(sprintf('Will simulate %d cells across %d stopping points (terminal: %d).',
                nrow(meta), length(unique(meta$timepoint)), sum(meta$is_terminal)))

ref_sce <- load_reference(opt$reference, opt$celltype_col, opt$pseudotime_col, opt$use_pseudotime)

cache_key <- ref_cache_key(opt$reference, opt$celltype_col,
                           opt$use_pseudotime, opt$pseudotime_col)
cache_path <- file.path(opt$cache_dir, paste0('fit_', cache_key, '.rds'))
fit <- fit_or_load(ref_sce, cache_path, opt$use_pseudotime, opt$ncores)

new_count <- simulate_for_meta(fit, meta, opt$ncores)

# ---- Save outputs ------------------------------------------------------------
suppressPackageStartupMessages(library(SingleCellExperiment))
sim_sce <- SingleCellExperiment(
  assays = list(counts = new_count),
  colData = S4Vectors::DataFrame(meta, row.names = meta$sample_id)
)
saveRDS(sim_sce, file.path(opt$output_dir, 'sim_sce_all.rds'))

for (tp in unique(meta$timepoint)) {
  sub_sce <- sim_sce[, meta$timepoint == tp]
  saveRDS(sub_sce,
          file.path(opt$output_dir, sprintf('sim_sce_time_%s.rds', as.character(tp))))
}

write.csv(meta, file.path(opt$output_dir, 'cell_metadata.csv'), row.names = FALSE)
cell_count_summary <- aggregate(cell_id ~ timepoint + cell_type + is_terminal,
                                data = meta, FUN = length)
colnames(cell_count_summary)[ncol(cell_count_summary)] <- 'n_cells'
write.csv(cell_count_summary,
          file.path(opt$output_dir, 'cell_count_summary.csv'), row.names = FALSE)

message(sprintf('Wrote outputs under %s', opt$output_dir))
