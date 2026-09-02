#!/usr/bin/env Rscript

# Command-line entry point that attaches simulated single-cell transcriptomes
# to a finished lineage-recording run (PhysiCell, standalone import, or
# Gillespie): it resolves the run directory, joins terminal-cell covariates,
# standardizes a reference dataset, fits or reuses a cached scDesign3 model,
# simulates counts, and writes them under lineage_recording/sc_profiles before
# printing SCDESIGN3_OUTPUT_DIR=<path>. Run with Rscript; it sources
# physicell_lineage.R, scdesign3_helpers.R, and physicell_scdesign3.R from its
# own directory and needs scDesign3 installed (jsonlite for --celltype-map).

script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- normalizePath(
  sub('^--file=', '', script_argument),
  mustWork = TRUE
)
repo_root <- dirname(script_path)
source(file.path(repo_root, 'physicell_lineage.R'))
source(file.path(repo_root, 'scdesign3_helpers.R'))
source(file.path(repo_root, 'physicell_scdesign3.R'))

#' Build the command-line usage text for this script
#'
#' A pure string builder with no arguments and no external dependencies; it is
#' emitted by `--help` and embedded in the argument-parsing error messages.
#'
#' @return A single newline-delimited character string listing the required
#'   flags, the optional flags, and their defaults.
physicell_scdesign3_usage <- function(){
  paste(
    'Usage:',
    '  Rscript generate_physicell_sc_profiles.R --run-dir DIR --reference REF.rds [options]',
    '',
    'Required:',
    '  -D, --run-dir PATH             PhysiCell pipeline run directory',
    '  -R, --reference PATH           Reference SCE/Seurat RDS or AnnData h5ad',
    '',
    'Options:',
    '  -O, --output-dir PATH          Default: RUN_DIR/lineage_recording/sc_profiles',
    '      --lineage-table PATH       Override PhysiCell lineage_table.csv',
    '  -C, --celltype-col NAME        Reference cell-type column; default: cell_type',
    '      --celltype-map PATH        Optional JSON simulated-to-reference mapping',
    '      --cell-type NAME           Override recording-manifest cell type',
    '      --use-pseudotime           Model lineage depth against pseudotime',
    '      --pseudotime-col NAME      Reference pseudotime column; default: pseudotime',
    '      --reference-spatial-cols X,Y  Optional reference coordinate columns',
    '      --other-covariates LIST    Comma-delimited matched covariates',
    '      --mu-formula FORMULA       Optional explicit marginal mean formula',
    '      --cache-dir PATH           Default: output/scdesign3_fits',
    '      --ncores N                 Default: 4',
    '      --max-cells N              Optional random testing cap',
    '      --seed N                   Default: 1',
    '      --compress-csv BOOL        Write CSV tables as .csv.gz; default: true',
    '  -h, --help',
    sep = '\n'
  )
}

#' Parse and type the scDesign3 profile-generation command line
#'
#' Walks the arguments left to right against a short/long alias table. Values
#' may be given either as a separate following argument or inline after an `=`
#' (the first `=` splits, and any later ones are kept in the value).
#' `--use-pseudotime` is the only bare switch. `--ncores`, `--max-cells`, and
#' `--seed` are coerced to integers and `--compress-csv` to a logical from a
#' case-insensitive `true` or `false`; other values stay character strings.
#'
#' @param arguments Character vector of trailing command-line arguments, as
#'   returned by `commandArgs(trailingOnly = TRUE)`.
#' @return A named list of options with the defaults already applied:
#'   `run_dir` and `reference` (both `NULL` until supplied), `output_dir`,
#'   `lineage_table`, `celltype_map`, `cell_type`, `reference_spatial_cols`,
#'   `other_covariates`, and `mu_formula` (`NULL`), `celltype_col`
#'   (`cell_type`), `use_pseudotime` (`FALSE`), `pseudotime_col`
#'   (`pseudotime`), `cache_dir` (`<repo>/output/scdesign3_fits`), `ncores`
#'   (`4`), `max_cells` (`NA_integer_`, meaning no cap), `seed` (`1`), and
#'   `compress_csv` (`TRUE`).
#' @section Side effects: On `-h` or `--help`, prints the usage text and ends
#'   the R session with status `0` rather than returning. Reads the
#'   script-level `repo_root` to build the default `cache_dir`.
#' @note Stops on an unrecognised flag, on a value-taking flag given without a
#'   value, on a non-integer value for an integer option, and on a
#'   `--compress-csv` value other than `true` or `false`. It does not check
#'   that the required options were supplied; the caller does that.
parse_physicell_scdesign3_args <- function(arguments){
  options <- list(
    run_dir = NULL,
    reference = NULL,
    output_dir = NULL,
    lineage_table = NULL,
    celltype_col = 'cell_type',
    celltype_map = NULL,
    cell_type = NULL,
    use_pseudotime = FALSE,
    pseudotime_col = 'pseudotime',
    reference_spatial_cols = NULL,
    other_covariates = NULL,
    mu_formula = NULL,
    cache_dir = file.path(repo_root, 'output', 'scdesign3_fits'),
    ncores = 4L,
    max_cells = NA_integer_,
    seed = 1L,
    compress_csv = TRUE
  )
  aliases <- c(
    '-D' = 'run_dir',
    '--run-dir' = 'run_dir',
    '-R' = 'reference',
    '--reference' = 'reference',
    '-O' = 'output_dir',
    '--output-dir' = 'output_dir',
    '--lineage-table' = 'lineage_table',
    '-C' = 'celltype_col',
    '--celltype-col' = 'celltype_col',
    '--celltype-map' = 'celltype_map',
    '--cell-type' = 'cell_type',
    '--pseudotime-col' = 'pseudotime_col',
    '--reference-spatial-cols' = 'reference_spatial_cols',
    '--other-covariates' = 'other_covariates',
    '--mu-formula' = 'mu_formula',
    '--cache-dir' = 'cache_dir',
    '--ncores' = 'ncores',
    '--max-cells' = 'max_cells',
    '--seed' = 'seed',
    '--compress-csv' = 'compress_csv'
  )
  integer_options <- c('ncores', 'max_cells', 'seed')
  argument_index <- 1L
  while(argument_index <= length(arguments)){
    argument <- arguments[argument_index]
    if(argument %in% c('-h', '--help')){
      cat(physicell_scdesign3_usage(), '\n')
      quit(save = 'no', status = 0)
    }
    if(argument == '--use-pseudotime'){
      options$use_pseudotime <- TRUE
      argument_index <- argument_index + 1L
      next
    }

    inline_parts <- strsplit(argument, '=', fixed = TRUE)[[1]]
    flag <- inline_parts[1]
    option_name <- unname(aliases[flag])
    if(length(option_name) == 0 || is.na(option_name)){
      stop(sprintf(
        'Unknown option: %s\n\n%s',
        flag,
        physicell_scdesign3_usage()
      ))
    }
    if(length(inline_parts) > 1){
      value <- paste(inline_parts[-1], collapse = '=')
    } else{
      argument_index <- argument_index + 1L
      if(argument_index > length(arguments)){
        stop(sprintf('Option %s requires a value.', flag))
      }
      value <- arguments[argument_index]
    }
    if(option_name %in% integer_options){
      value <- suppressWarnings(as.integer(value))
      if(is.na(value)){
        stop(sprintf('Option %s requires an integer value.', flag))
      }
    } else if(option_name == 'compress_csv'){
      normalized_value <- tolower(value)
      if(!(normalized_value %in% c('true', 'false'))){
        stop(sprintf('Option %s requires true or false.', flag))
      }
      value <- normalized_value == 'true'
    }
    options[[option_name]] <- value
    argument_index <- argument_index + 1L
  }
  options
}

options <- parse_physicell_scdesign3_args(commandArgs(trailingOnly = TRUE))
if(is.null(options$run_dir) || is.null(options$reference)){
  stop(sprintf(
    '--run-dir and --reference are required.\n\n%s',
    physicell_scdesign3_usage()
  ))
}
if(length(options$ncores) != 1 || is.na(options$ncores) || options$ncores < 1){
  stop('--ncores must be one positive integer.')
}
if(!is.na(options$max_cells) && options$max_cells < 1){
  stop('--max-cells must be positive when supplied.')
}
set.seed(options$seed)

run_dir <- normalizePath(options$run_dir, mustWork = TRUE)
reference_path <- normalizePath(options$reference, mustWork = TRUE)
recording_dir <- if(dir.exists(file.path(run_dir, 'lineage_recording'))){
  file.path(run_dir, 'lineage_recording')
} else if(!is.na(resolve_physicell_csv_path(
  file.path(run_dir, 'lineage_nodes.csv'),
  required = FALSE
))){
  # Standalone PhysiCell imports and Gillespie runs write directly into their
  # selected output directory.
  run_dir
} else{
  stop(sprintf('Lineage-recording tables not found under: %s.', run_dir))
}
lineage_table_path <- if(is.null(options$lineage_table)){
  gillespie_cells <- resolve_physicell_csv_path(
    file.path(recording_dir, 'gillespie_cell_states.csv'),
    required = FALSE
  )
  if(!is.na(gillespie_cells)){
    gillespie_cells
  } else{
    file.path(
      run_dir,
      'physicell_build',
      'output',
      'lineage_table.csv'
    )
  }
} else{
  options$lineage_table
}
lineage_table_path <- normalizePath(lineage_table_path, mustWork = TRUE)
output_dir <- if(is.null(options$output_dir)){
  file.path(recording_dir, 'sc_profiles')
} else{
  options$output_dir
}
if(!grepl('^/', output_dir)){
  output_dir <- file.path(repo_root, output_dir)
}
cache_dir <- options$cache_dir
if(!grepl('^/', cache_dir)){
  cache_dir <- file.path(repo_root, cache_dir)
}

spatial_columns <- split_scdesign3_columns(
  options$reference_spatial_cols,
  expected_length = if(is.null(options$reference_spatial_cols)) NULL else 2
)
other_covariates <- split_scdesign3_columns(options$other_covariates)
scdesign3_log_stage(
  'scDesign3: assembling lineage, spatial, and recording covariates.'
)
metadata <- build_physicell_sc_covariates(
  recording_dir,
  lineage_table_path,
  cell_type = options$cell_type
)
if(!is.null(options$celltype_map)){
  if(!requireNamespace('jsonlite', quietly = TRUE)){
    stop('jsonlite is required when --celltype-map is supplied.')
  }
  celltype_map_path <- normalizePath(options$celltype_map, mustWork = TRUE)
  celltype_map <- jsonlite::fromJSON(
    celltype_map_path,
    simplifyVector = TRUE
  )
  metadata <- apply_physicell_celltype_map(metadata, celltype_map)
}

scdesign3_log_stage('scDesign3: loading and standardizing reference data.')
reference_sce <- load_scdesign3_reference(
  reference_path,
  celltype_col = options$celltype_col,
  pseudotime_col = options$pseudotime_col,
  use_pseudotime = options$use_pseudotime,
  spatial_cols = spatial_columns,
  other_covariates = other_covariates
)
metadata <- align_physicell_sc_covariates(
  metadata,
  reference_sce,
  use_pseudotime = options$use_pseudotime,
  use_spatial = length(spatial_columns) == 2
)
if(!is.na(options$max_cells) && nrow(metadata) > options$max_cells){
  metadata <- metadata[
    sample(seq_len(nrow(metadata)), options$max_cells),
    ,
    drop = FALSE
  ]
}
scdesign3_log_stage(sprintf(
  'scDesign3: prepared %s target cells and %s reference cells.',
  format(nrow(metadata), big.mark = ',', scientific = FALSE, trim = TRUE),
  format(ncol(reference_sce), big.mark = ',', scientific = FALSE, trim = TRUE)
))

missing_simulated_covariates <- setdiff(other_covariates, names(metadata))
if(length(missing_simulated_covariates) > 0){
  stop(sprintf(
    'Requested covariates are absent from PhysiCell metadata: %s.',
    paste(missing_simulated_covariates, collapse = ', ')
  ))
}
mu_formula <- options$mu_formula
if(is.null(mu_formula) || !nzchar(mu_formula)){
  mu_formula <- scdesign3_default_mu_formula(
    reference_sce,
    use_pseudotime = options$use_pseudotime,
    use_spatial = length(spatial_columns) == 2,
    other_covariates = other_covariates
  )
}
cache_key <- scdesign3_fit_cache_key(
  reference_path,
  celltype_col = options$celltype_col,
  pseudotime_col = options$pseudotime_col,
  use_pseudotime = options$use_pseudotime,
  spatial_cols = spatial_columns,
  other_covariates = other_covariates,
  mu_formula = mu_formula
)
cache_path <- file.path(cache_dir, paste0('fit_', cache_key, '.rds'))
scdesign3_log_stage(
  'scDesign3: loading or fitting the reference marginal/copula model.'
)
fit <- fit_or_load_scdesign3(
  reference_sce,
  cache_path,
  use_pseudotime = options$use_pseudotime,
  use_spatial = length(spatial_columns) == 2,
  other_covariates = other_covariates,
  mu_formula = mu_formula,
  ncores = options$ncores
)
scdesign3_log_stage(
  'scDesign3 model ready; beginning transcriptome simulation.'
)
counts <- simulate_scdesign3_counts(
  fit,
  metadata,
  ncores = options$ncores
)
scdesign3_log_stage(sprintf(
  paste(
    'scDesign3 transcriptome simulation finished for %s cells and %s',
    'genes; beginning output.'
  ),
  format(ncol(counts), big.mark = ',', scientific = FALSE, trim = TRUE),
  format(nrow(counts), big.mark = ',', scientific = FALSE, trim = TRUE)
))
write_physicell_scdesign3_outputs(
  counts,
  metadata,
  output_dir,
  fit,
  reference_path,
  options$seed,
  compress_csv = options$compress_csv
)
scdesign3_log_stage('scDesign3 output complete.')

cat(sprintf(
  'Generated %d-gene scDesign3 profiles for %d PhysiCell terminal cells.\n',
  nrow(counts),
  ncol(counts)
))
cat(sprintf(
  'SCDESIGN3_OUTPUT_DIR=%s\n',
  normalizePath(output_dir, mustWork = TRUE)
))
