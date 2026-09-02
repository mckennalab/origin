# Shared scDesign3 reference loading, model fitting, and conditional simulation.
# These helpers do not execute a CLI when sourced.

#' Emit a timestamped scDesign3 pipeline-stage message
#'
#' @param message_text One non-empty string naming the stage being entered.
#' @return `TRUE`, invisibly.
#' @section Side effects:
#' Writes `[<local timestamp>] <message_text>` to the message connection
#' (stderr).
scdesign3_log_stage <- function(message_text){
  message_text <- as.character(message_text)
  if(length(message_text) != 1 || is.na(message_text) ||
     !nzchar(message_text)){
    stop('scDesign3 stage log message must be one non-empty string.')
  }
  message(sprintf(
    '[%s] %s',
    format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z'),
    message_text
  ))
  invisible(TRUE)
}

#' Parse a comma-delimited CLI value into validated covariate column names
#'
#' Splits on commas, trims surrounding whitespace, drops empty entries, and
#' rejects names that `make.names()` would alter, because the results are pasted
#' into a model formula.
#'
#' @param value Comma-delimited string of column names. `NULL`, zero-length,
#'   `NA`, and the empty string all yield no columns.
#' @param expected_length Optional exact count of names required; the call fails
#'   when the parsed count differs. Checked before de-duplication.
#' @return A character vector of unique, syntactically valid column names,
#'   possibly zero-length.
split_scdesign3_columns <- function(value, expected_length = NULL){
  if(is.null(value) || length(value) == 0 || is.na(value) || !nzchar(value)){
    return(character())
  }
  columns <- trimws(strsplit(value, ',', fixed = TRUE)[[1]])
  columns <- columns[nzchar(columns)]
  if(!is.null(expected_length) && length(columns) != expected_length){
    stop(sprintf(
      'Expected %d comma-delimited column names, received %d.',
      expected_length,
      length(columns)
    ))
  }
  if(any(make.names(columns) != columns)){
    stop('scDesign3 covariate column names must be syntactically valid R names.')
  }
  unique(columns)
}

#' Load a reference dataset and standardize its scDesign3 predictor columns
#'
#' Reads an `.h5ad`, or an `.rds` holding a `SingleCellExperiment` or a `Seurat`
#' object, and normalizes it into the shape the rest of the pipeline assumes:
#' an assay literally named `counts` and `colData` columns named `cell_type`,
#' `pseudotime`, `spatial1`, and `spatial2`. The returned object keeps the
#' Bioconductor orientation of genes x cells — rows are features, columns are
#' cells — which is the orientation `fit_or_load_scdesign3()` and
#' `simulate_scdesign3_counts()` require.
#'
#' @details
#' When no assay is named `counts`, the first assay is aliased to `counts` with
#' a message; the raw counts the negative-binomial marginals assume are the
#' caller's responsibility. Source column names are copied into the standard
#' names rather than renamed, so the originals remain.
#'
#' @param path Reference file; `.h5ad` needs `zellkonverter` and a Seurat
#'   `.rds` needs `Seurat`. Must exist.
#' @param celltype_col `colData` column holding cell-type labels; copied to
#'   `cell_type` as a factor whose levels bound what
#'   `scdesign3_new_covariates()` will later accept. Missing labels fail.
#' @param pseudotime_col `colData` column holding pseudotime; used only when
#'   `use_pseudotime` is `TRUE`.
#' @param use_pseudotime When `TRUE`, copy `pseudotime_col` to `pseudotime`; it
#'   must be finite and hold at least two distinct values.
#' @param spatial_cols Either zero or exactly two `colData` column names, copied
#'   to `spatial1` and `spatial2`; they must be finite and span at least ten
#'   distinct locations.
#' @param other_covariates Additional `colData` columns to retain as predictors.
#'   Character columns become factors; numeric ones must be finite, and none may
#'   contain `NA`.
#' @return The standardized `SingleCellExperiment` (genes x cells).
#' @note Requires `SingleCellExperiment` and `SummarizedExperiment`, but not
#'   `scDesign3` itself.
load_scdesign3_reference <- function(path,
                                     celltype_col = 'cell_type',
                                     pseudotime_col = 'pseudotime',
                                     use_pseudotime = FALSE,
                                     spatial_cols = character(),
                                     other_covariates = character()){
  required_packages <- c('SingleCellExperiment', 'SummarizedExperiment')
  missing_packages <- required_packages[!vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )]
  if(length(missing_packages) > 0){
    stop(sprintf(
      'Reference loading requires: %s.',
      paste(missing_packages, collapse = ', ')
    ))
  }

  path <- normalizePath(path, mustWork = TRUE)
  extension <- tolower(tools::file_ext(path))
  if(extension == 'h5ad'){
    if(!requireNamespace('zellkonverter', quietly = TRUE)){
      stop('zellkonverter is required to read an h5ad reference.')
    }
    sce <- zellkonverter::readH5AD(path)
  } else if(extension == 'rds'){
    object <- readRDS(path)
    if(methods::is(object, 'SingleCellExperiment')){
      sce <- object
    } else if(inherits(object, 'Seurat')){
      if(!requireNamespace('Seurat', quietly = TRUE)){
        stop('Seurat is required to convert a Seurat RDS reference.')
      }
      sce <- Seurat::as.SingleCellExperiment(object)
    } else{
      stop(sprintf(
        'Unsupported reference R object class: %s.',
        paste(class(object), collapse = ', ')
      ))
    }
  } else{
    stop('Reference must be a SingleCellExperiment/Seurat RDS or an h5ad file.')
  }

  metadata <- as.data.frame(SummarizedExperiment::colData(sce))
  required_columns <- unique(c(
    celltype_col,
    if(isTRUE(use_pseudotime)) pseudotime_col else character(),
    spatial_cols,
    other_covariates
  ))
  missing_columns <- setdiff(required_columns, names(metadata))
  if(length(missing_columns) > 0){
    stop(sprintf(
      'Reference colData is missing: %s.',
      paste(missing_columns, collapse = ', ')
    ))
  }
  if(length(spatial_cols) != 0 && length(spatial_cols) != 2){
    stop('spatial_cols must contain exactly two reference column names.')
  }

  assay_names <- SummarizedExperiment::assayNames(sce)
  if(!('counts' %in% assay_names)){
    if(length(assay_names) == 0){
      stop('Reference has no assays.')
    }
    fallback_assay <- assay_names[1]
    message(sprintf(
      'Aliasing reference assay %s to counts.',
      fallback_assay
    ))
    SummarizedExperiment::assay(sce, 'counts') <-
      SummarizedExperiment::assay(sce, fallback_assay)
  }

  SummarizedExperiment::colData(sce)$cell_type <- factor(
    as.character(metadata[[celltype_col]])
  )
  if(anyNA(SummarizedExperiment::colData(sce)$cell_type)){
    stop('Reference cell-type labels cannot be missing.')
  }

  if(isTRUE(use_pseudotime)){
    pseudotime <- as.numeric(metadata[[pseudotime_col]])
    if(any(!is.finite(pseudotime))){
      stop('Reference pseudotime must be finite numeric data.')
    }
    if(length(unique(pseudotime)) < 2){
      stop('Reference pseudotime must contain at least two distinct values.')
    }
    SummarizedExperiment::colData(sce)$pseudotime <- pseudotime
  }
  if(length(spatial_cols) == 2){
    spatial_values <- metadata[, spatial_cols, drop = FALSE]
    if(nrow(unique(spatial_values)) < 10){
      stop('Reference spatial coordinates need at least ten distinct locations.')
    }
    for(index in seq_len(2)){
      values <- as.numeric(metadata[[spatial_cols[index]]])
      if(any(!is.finite(values))){
        stop('Reference spatial coordinates must be finite numeric data.')
      }
      SummarizedExperiment::colData(sce)[[paste0('spatial', index)]] <- values
    }
  }

  for(covariate in other_covariates){
    values <- metadata[[covariate]]
    if(is.character(values)){
      values <- factor(values)
    }
    if(is.numeric(values) && any(!is.finite(values))){
      stop(sprintf('Reference covariate %s must be finite.', covariate))
    }
    if(anyNA(values)){
      stop(sprintf('Reference covariate %s cannot contain missing values.', covariate))
    }
    SummarizedExperiment::colData(sce)[[covariate]] <- values
  }
  sce
}

#' Build the default marginal-mean formula from the standardized reference
#'
#' Assembles the right-hand side of the per-gene mean model, including only the
#' terms the reference can actually support: `cell_type` when more than one
#' level is present, a pseudotime term when requested, a spatial term when
#' requested, and any additional covariates verbatim. Pseudotime uses a cubic
#' regression spline `s(pseudotime, k, bs = "cr")` with `k` capped at 4 and at
#' one below the number of distinct values, falling back to a linear
#' `pseudotime` term when fewer than four distinct values exist. Spatial uses a
#' Gaussian-process smooth `s(spatial1, spatial2, k, bs = "gp")` with `k` set to
#' one quarter of the reference cell count, clamped to `[5, 50]`.
#'
#' @param sce Standardized reference `SingleCellExperiment` (genes x cells);
#'   `ncol(sce)` is read as the reference cell count.
#' @param use_pseudotime When `TRUE`, add the pseudotime term.
#' @param use_spatial When `TRUE`, add the two-dimensional spatial smooth.
#' @param other_covariates Character vector of extra terms appended as-is.
#' @return A one-element character formula right-hand side, or `"1"` when no
#'   term applies (an intercept-only mean model).
scdesign3_default_mu_formula <- function(sce,
                                         use_pseudotime = FALSE,
                                         use_spatial = FALSE,
                                         other_covariates = character()){
  metadata <- as.data.frame(SummarizedExperiment::colData(sce))
  terms <- character()
  if(length(unique(metadata$cell_type)) > 1){
    terms <- c(terms, 'cell_type')
  }
  if(isTRUE(use_pseudotime)){
    unique_pseudotime <- length(unique(metadata$pseudotime))
    if(unique_pseudotime >= 4){
      pseudotime_k <- min(4L, unique_pseudotime - 1L)
      terms <- c(terms, sprintf(
        's(pseudotime, k = %d, bs = "cr")',
        pseudotime_k
      ))
    } else{
      terms <- c(terms, 'pseudotime')
    }
  }
  if(isTRUE(use_spatial)){
    spatial_k <- max(5L, min(50L, floor(ncol(sce) / 4)))
    terms <- c(terms, sprintf(
      's(spatial1, spatial2, k = %d, bs = "gp")',
      spatial_k
    ))
  }
  terms <- c(terms, other_covariates)
  if(length(terms) == 0){
    '1'
  } else{
    paste(terms, collapse = ' + ')
  }
}

#' Derive a cache key covering everything a fitted scDesign3 model depends on
#'
#' Hashes the reference identity (normalized path, mtime, size), the installed
#' `scDesign3` version, and the full model specification, so that editing the
#' reference file or changing any predictor choice yields a different key and
#' forces a refit.
#'
#' @param reference_path Path to the reference dataset; must exist.
#' @param celltype_col Reference cell-type column name.
#' @param pseudotime_col Reference pseudotime column name.
#' @param use_pseudotime Whether pseudotime is a predictor.
#' @param spatial_cols Reference spatial column names, joined with commas.
#' @param other_covariates Additional predictor names, joined with commas.
#' @param mu_formula Explicit marginal-mean formula, or `NULL` to record the
#'   literal `<default>` placeholder.
#' @param family_use Marginal count family, `'nb'` (negative binomial) by
#'   default, matching `fit_or_load_scdesign3()`.
#' @return A single string: a SHA-1 digest when `digest` is installed, otherwise
#'   an 8-hex-digit fallback checksum of the payload.
#' @note Uses `scDesign3` only to read its version, recording
#'   `'not-installed'` when the package is absent, so the key is computable
#'   without it.
scdesign3_fit_cache_key <- function(reference_path,
                                    celltype_col,
                                    pseudotime_col,
                                    use_pseudotime,
                                    spatial_cols = character(),
                                    other_covariates = character(),
                                    mu_formula = NULL,
                                    family_use = 'nb'){
  reference_path <- normalizePath(reference_path, mustWork = TRUE)
  file_information <- file.info(reference_path)
  package_version <- if(requireNamespace('scDesign3', quietly = TRUE)){
    as.character(utils::packageVersion('scDesign3'))
  } else{
    'not-installed'
  }
  payload <- paste(
    reference_path,
    file_information$mtime,
    file_information$size,
    package_version,
    celltype_col,
    pseudotime_col,
    use_pseudotime,
    paste(spatial_cols, collapse = ','),
    paste(other_covariates, collapse = ','),
    if(is.null(mu_formula)) '<default>' else mu_formula,
    family_use,
    sep = '|'
  )
  if(requireNamespace('digest', quietly = TRUE)){
    digest::digest(payload, algo = 'sha1')
  } else{
    raw_values <- utf8ToInt(payload)
    sprintf('%08x', sum(raw_values * seq_along(raw_values)) %% .Machine$integer.max)
  }
}

#' Fit the scDesign3 marginal and copula models, or reuse a cached fit
#'
#' Returns a cached fit when `cache_path` holds one carrying every expected
#' field, otherwise runs the three scDesign3 stages against the reference and
#' caches the result.
#'
#' @details
#' The fit is fixed to these model choices. `construct_data()` reads the
#' `counts` assay of `sce` (genes x cells) with `celltype = 'cell_type'` and
#' `corr_by = '1'`, so all reference cells form a single correlation group.
#' `fit_marginal()` fits one `family_use` GLM/GAM per gene (`predictor =
#' 'gene'`) with mean model `mu_formula` and `sigma_formula = '1'`, so the
#' dispersion is constant across cells within a gene; `usebam = FALSE`.
#' `fit_copula()` then fits a Gaussian copula (`copula = 'gaussian'`) over the
#' top 80% of features by expression (`important_feature = 0.8`) with
#' `if_sparse = FALSE`. All three stages parallelize via `mcmapply`, which is
#' fork-based and therefore Unix-only.
#'
#' @param sce Standardized reference `SingleCellExperiment` from
#'   `load_scdesign3_reference()`, in genes x cells orientation with a `counts`
#'   assay.
#' @param cache_path RDS path read for an existing fit and written on a refit;
#'   a cache missing any required field is ignored with a message.
#' @param use_pseudotime When `TRUE`, `pseudotime` is passed to scDesign3 as the
#'   pseudotime covariate.
#' @param use_spatial When `TRUE`, `spatial1` and `spatial2` are passed as the
#'   spatial covariates.
#' @param other_covariates Additional predictor column names, or a zero-length
#'   vector for none.
#' @param mu_formula Marginal-mean formula; `NULL` or empty falls back to
#'   `scdesign3_default_mu_formula()`.
#' @param family_use Marginal count distribution passed to both `fit_marginal()`
#'   and `fit_copula()`; `'nb'` (negative binomial) by default.
#' @param ncores Worker count for all three stages; must be one positive
#'   integer.
#' @return A named list with `sce`, `data` (the `construct_data()` result),
#'   `marginal_list`, `copula_fit`, `mu_formula`, `family_use`,
#'   `use_pseudotime`, `use_spatial`, `other_covariates`, and
#'   `scdesign3_version`.
#' @section Side effects:
#' Creates `dirname(cache_path)` and writes the fit to `cache_path` after a
#' refit; forks worker processes during fitting.
#' @note Requires the `scDesign3` and `SummarizedExperiment` packages.
fit_or_load_scdesign3 <- function(sce,
                                  cache_path,
                                  use_pseudotime = FALSE,
                                  use_spatial = FALSE,
                                  other_covariates = character(),
                                  mu_formula = NULL,
                                  family_use = 'nb',
                                  ncores = 1){
  if(!requireNamespace('scDesign3', quietly = TRUE)){
    stop('The scDesign3 package is required for expression simulation.')
  }
  if(!requireNamespace('SummarizedExperiment', quietly = TRUE)){
    stop('The SummarizedExperiment package is required.')
  }
  ncores <- as.integer(ncores)
  if(length(ncores) != 1 || is.na(ncores) || ncores < 1){
    stop('ncores must be one positive integer.')
  }
  if(file.exists(cache_path)){
    fit <- readRDS(cache_path)
    required_fields <- c(
      'sce', 'data', 'marginal_list', 'copula_fit', 'mu_formula',
      'family_use', 'use_pseudotime', 'use_spatial', 'other_covariates'
    )
    if(all(required_fields %in% names(fit))){
      message(sprintf('Loading cached scDesign3 fit: %s', cache_path))
      return(fit)
    }
    message(sprintf('Ignoring incompatible scDesign3 cache: %s', cache_path))
  }

  if(is.null(mu_formula) || !nzchar(mu_formula)){
    mu_formula <- scdesign3_default_mu_formula(
      sce,
      use_pseudotime,
      use_spatial,
      other_covariates
    )
  }
  pseudotime_argument <- if(isTRUE(use_pseudotime)) 'pseudotime' else NULL
  spatial_argument <- if(isTRUE(use_spatial)){
    c('spatial1', 'spatial2')
  } else{
    NULL
  }
  other_argument <- if(length(other_covariates) > 0){
    other_covariates
  } else{
    NULL
  }

  message(sprintf('Fitting scDesign3 marginal model: %s', mu_formula))
  data <- scDesign3::construct_data(
    sce = sce,
    assay_use = 'counts',
    celltype = 'cell_type',
    pseudotime = pseudotime_argument,
    spatial = spatial_argument,
    other_covariates = other_argument,
    corr_by = '1',
    parallelization = 'mcmapply'
  )
  marginal_list <- scDesign3::fit_marginal(
    data = data,
    predictor = 'gene',
    mu_formula = mu_formula,
    sigma_formula = '1',
    family_use = family_use,
    n_cores = ncores,
    usebam = FALSE,
    parallelization = 'mcmapply'
  )
  copula_fit <- scDesign3::fit_copula(
    sce = sce,
    assay_use = 'counts',
    input_data = data$dat,
    marginal_list = marginal_list,
    family_use = family_use,
    copula = 'gaussian',
    important_feature = 0.8,
    if_sparse = FALSE,
    n_cores = ncores,
    parallelization = 'mcmapply'
  )

  fit <- list(
    sce = sce,
    data = data,
    marginal_list = marginal_list,
    copula_fit = copula_fit,
    mu_formula = mu_formula,
    family_use = family_use,
    use_pseudotime = isTRUE(use_pseudotime),
    use_spatial = isTRUE(use_spatial),
    other_covariates = other_covariates,
    scdesign3_version = as.character(utils::packageVersion('scDesign3'))
  )
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(fit, cache_path)
  message(sprintf('Cached scDesign3 fit at %s', cache_path))
  fit
}

#' Project simulation metadata onto the fitted model's covariate frame
#'
#' Selects exactly the predictor columns the fit was built on (every column of
#' `fit$data$dat` except `corr_group`) and coerces each to match the reference's
#' type: factors are rebuilt on the reference levels, and any level the
#' reference never saw is an error rather than a silent `NA`. `corr_group` is
#' set to `1` for every cell, matching the single correlation group
#' `fit_or_load_scdesign3()` fits with `corr_by = '1'`.
#'
#' @param fit Fit list from `fit_or_load_scdesign3()`; only `fit$data$dat` is
#'   read, itself a cells x covariates frame over the reference cells.
#' @param metadata Target-cell covariate data frame; must contain every
#'   predictor column, with finite values for the numeric ones.
#' @return A data frame with one row per row of `metadata`, holding the
#'   predictor columns in the reference's order and types plus `corr_group`.
scdesign3_new_covariates <- function(fit, metadata){
  reference_covariates <- fit$data$dat
  predictor_columns <- setdiff(names(reference_covariates), 'corr_group')
  missing_columns <- setdiff(predictor_columns, names(metadata))
  if(length(missing_columns) > 0){
    stop(sprintf(
      'Simulation metadata is missing scDesign3 covariates: %s.',
      paste(missing_columns, collapse = ', ')
    ))
  }

  result <- metadata[, predictor_columns, drop = FALSE]
  for(column in predictor_columns){
    reference_values <- reference_covariates[[column]]
    if(is.factor(reference_values)){
      values <- as.character(result[[column]])
      unknown <- setdiff(unique(values), levels(reference_values))
      if(length(unknown) > 0){
        stop(sprintf(
          'Covariate %s contains levels absent from the reference: %s.',
          column,
          paste(unknown, collapse = ', ')
        ))
      }
      result[[column]] <- factor(values, levels = levels(reference_values))
    } else if(is.numeric(reference_values) || is.integer(reference_values)){
      values <- as.numeric(result[[column]])
      if(any(!is.finite(values))){
        stop(sprintf('Covariate %s must be finite numeric data.', column))
      }
      result[[column]] <- values
    } else{
      result[[column]] <- result[[column]]
    }
  }
  result$corr_group <- 1
  result
}

#' Simulate a count matrix for target cells from a fitted scDesign3 model
#'
#' Evaluates the fitted per-gene marginals at the target covariates with
#' `extract_para()`, then draws correlated counts with `simu_new()` under the
#' fitted Gaussian copula.
#'
#' @details
#' The marginal family is whatever the fit was built with (`fit$family_use`,
#' `'nb'` by default), and the copula, correlation grouping, and important-
#' feature set are reused from `fit$copula_fit` — none of them can be varied
#' here. `quantile_mat` is passed as `NULL`, so counts are sampled from the
#' fitted marginals rather than quantile-matched to observed values, and
#' `fit$data$filtered_gene` carries through the genes scDesign3 excluded during
#' fitting. Work is parallelized with `mcmapply`.
#'
#' @param fit Fit list from `fit_or_load_scdesign3()`.
#' @param metadata Target-cell covariate data frame; must carry unique
#'   non-missing `sample_id` values and every predictor the fit needs.
#' @param ncores Worker count for both stages; must be one positive integer.
#' @return A count matrix in genes x cells orientation: rows are the reference
#'   genes, named by `rownames(fit$sce)`, and columns are the target cells,
#'   named by `metadata$sample_id`.
#' @note Calls `scDesign3::extract_para()` and `scDesign3::simu_new()` directly
#'   with no availability check, so the package must be installed; unlike
#'   `fit_or_load_scdesign3()` a missing package surfaces as a namespace-load
#'   error rather than an explanatory message.
simulate_scdesign3_counts <- function(fit, metadata, ncores = 1){
  ncores <- as.integer(ncores)
  if(length(ncores) != 1 || is.na(ncores) || ncores < 1){
    stop('ncores must be one positive integer.')
  }
  if(!('sample_id' %in% names(metadata)) ||
     anyNA(metadata$sample_id) ||
     anyDuplicated(metadata$sample_id)){
    stop('metadata must contain unique non-missing sample_id values.')
  }
  new_covariates <- scdesign3_new_covariates(fit, metadata)
  parameters <- scDesign3::extract_para(
    sce = fit$sce,
    assay_use = 'counts',
    marginal_list = fit$marginal_list,
    n_cores = ncores,
    family_use = fit$family_use,
    new_covariate = new_covariates,
    data = fit$data$dat,
    parallelization = 'mcmapply'
  )
  counts <- scDesign3::simu_new(
    sce = fit$sce,
    assay_use = 'counts',
    mean_mat = parameters$mean_mat,
    sigma_mat = parameters$sigma_mat,
    zero_mat = parameters$zero_mat,
    quantile_mat = NULL,
    copula_list = fit$copula_fit$copula_list,
    n_cores = ncores,
    family_use = fit$family_use,
    input_data = fit$data$dat,
    new_covariate = new_covariates,
    important_feature = fit$copula_fit$important_feature,
    parallelization = 'mcmapply',
    filtered_gene = fit$data$filtered_gene
  )
  rownames(counts) <- rownames(fit$sce)
  colnames(counts) <- metadata$sample_id
  counts
}
