# Shared scDesign3 reference loading, model fitting, and conditional simulation.
# These helpers do not execute a CLI when sourced.

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
