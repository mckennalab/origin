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

build_covariate_df <- function(cellpop_files, celltype_map = NULL) {
  per_tp <- lapply(seq_len(nrow(cellpop_files)), function(i) {
    cp <- readRDS(cellpop_files$path[i])
    alive <- vapply(cp, function(c) isTRUE(c$alive), logical(1))
    cp <- cp[alive]
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
  meta
}

add_pseudotime <- function(meta) {
  max_depth <- max(meta$lineage_depth, na.rm = TRUE)
  meta$pseudotime <- if (max_depth == 0) 0 else meta$lineage_depth / max_depth
  meta
}

maybe_downsample <- function(meta, max_per_tp) {
  if (is.na(max_per_tp)) return(meta)
  do.call(rbind, lapply(split(meta, meta$timepoint), function(sub) {
    if (nrow(sub) <= max_per_tp) return(sub)
    sub[sample(nrow(sub), max_per_tp), , drop = FALSE]
  }))
}

# ---- scDesign3 fit (cached) and simulate --------------------------------------
ref_cache_key <- function(ref_path, celltype_col, use_pseudotime, pseudotime_col) {
  fi <- file.info(ref_path)
  payload <- paste(normalizePath(ref_path), fi$mtime, fi$size,
                   celltype_col, use_pseudotime, pseudotime_col, sep = '|')
  if (requireNamespace('digest', quietly = TRUE)) {
    digest::digest(payload, algo = 'sha1')
  } else {
    sprintf('%x', sum(utf8ToInt(payload)))
  }
}

fit_or_load <- function(sce, cache_path, use_pseudotime, ncores) {
  if (file.exists(cache_path)) {
    message(sprintf('Loading cached scDesign3 fit: %s', cache_path))
    return(readRDS(cache_path))
  }
  suppressPackageStartupMessages(library(scDesign3))
  mu_formula <- if (use_pseudotime) {
    'cell_type + s(pseudotime, k = 4, bs = "cr")'
  } else {
    'cell_type'
  }
  pseudotime_arg <- if (use_pseudotime) 'pseudotime' else NULL
  message('Fitting scDesign3 on reference (cached afterwards) ...')
  dat <- construct_data(sce = sce, assay_use = 'counts',
                        celltype = 'cell_type',
                        pseudotime = pseudotime_arg,
                        spatial = NULL, other_covariates = NULL,
                        corr_by = '1',
                        parallelization = 'mclapply', n_cores = ncores)
  marginal_list <- fit_marginal(data = dat, predictor = 'gene',
                                mu_formula = mu_formula,
                                sigma_formula = '1',
                                family_use = 'nb',
                                n_cores = ncores, usebam = FALSE)
  copula_list <- fit_copula(sce = sce, assay_use = 'counts',
                            marginal_list = marginal_list,
                            family_use = 'nb', copula = 'gaussian',
                            n_cores = ncores, input_data = dat$dat)
  fit <- list(sce = sce, dat = dat,
              marginal_list = marginal_list,
              copula_list = copula_list,
              use_pseudotime = use_pseudotime,
              mu_formula = mu_formula)
  saveRDS(fit, cache_path)
  message(sprintf('Cached fit at %s', cache_path))
  fit
}

simulate_for_meta <- function(fit, new_meta, ncores) {
  suppressPackageStartupMessages({
    library(scDesign3)
    library(SingleCellExperiment)
  })
  ref_levels <- levels(colData(fit$sce)$cell_type)
  unknown <- setdiff(unique(new_meta$cell_type), ref_levels)
  if (length(unknown) > 0) {
    stop(sprintf('Cell types absent from reference: %s. Use --celltype_map to remap.',
                 paste(unknown, collapse = ', ')))
  }
  new_cov <- if (fit$use_pseudotime) {
    data.frame(cell_type = factor(new_meta$cell_type, levels = ref_levels),
               pseudotime = new_meta$pseudotime,
               corr_group = factor('1'),
               stringsAsFactors = FALSE)
  } else {
    data.frame(cell_type = factor(new_meta$cell_type, levels = ref_levels),
               corr_group = factor('1'),
               stringsAsFactors = FALSE)
  }

  para_new <- extract_para(sce = fit$sce, marginal_list = fit$marginal_list,
                           n_cores = ncores, family_use = 'nb',
                           new_covariate = new_cov, data = fit$dat$dat)
  new_count <- simu_new(sce = fit$sce,
                        mean_mat = para_new$mean_mat,
                        sigma_mat = para_new$sigma_mat,
                        zero_mat = para_new$zero_mat,
                        quantile_mat = NULL,
                        copula_list = fit$copula_list,
                        n_cores = ncores,
                        family_use = 'nb',
                        input_data = fit$dat$dat,
                        new_covariate = new_cov,
                        important_feature = rep(TRUE, nrow(fit$sce)))
  rownames(new_count) <- rownames(fit$sce)
  colnames(new_count) <- new_meta$cell_id
  new_count
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
  colData = S4Vectors::DataFrame(meta, row.names = meta$cell_id)
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
