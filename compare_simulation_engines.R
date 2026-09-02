#!/usr/bin/env Rscript

# Command-line ENTRY POINT for the engine-comparison benchmark; the paired
# LIBRARY is `engine_comparison.R`, which defines the comparison and summary
# functions this driver calls. Sources `physicell_lineage.R`,
# `gillespie_lineage.R`, and `engine_comparison.R` from its own directory,
# runs the paired continuous-time vs. time-step benchmark on one native JSON
# parameter file, and writes four CSV tables plus one PNG to the output
# directory. Run: Rscript compare_simulation_engines.R --params PARAMS.json

script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- normalizePath(sub('^--file=', '', script_argument), mustWork = TRUE)
repo_root <- dirname(script_path)
source(file.path(repo_root, 'physicell_lineage.R'))
source(file.path(repo_root, 'gillespie_lineage.R'))
source(file.path(repo_root, 'engine_comparison.R'))

#' Build the command-line usage text for this driver
#'
#' @return One string containing the full usage message, with embedded
#'   newlines, listing the required `--params` option and every optional flag
#'   with its default.
comparison_cli_usage <- function(){
  paste(
    'Usage:',
    '  Rscript compare_simulation_engines.R --params PARAMS.json [options]',
    '',
    'Required:',
    '  -P, --params PATH          Native remote_mito JSON parameter file',
    '',
    'Options:',
    '  -O, --output-dir PATH     Default: output/engine_comparison/<JSON stem>',
    '      --replicates N        Replicates per engine; default: 50',
    '      --end-time NUMBER     Default: max(sim_length)',
    '      --time-step NUMBER    Default: numeric JSON time_inc',
    '      --max-cells N         Per-replicate safety limit; default: 1000000',
    '      --seed N              First replicate seed; default: 1',
    '      --progress BOOL       Default: true',
    '      --progress-updates N  Default: 10',
    '      --compress-csv BOOL   Default: true',
    '  -h, --help',
    sep = '\n'
  )
}

#' Parse and type-check this driver's command-line arguments
#'
#' Accepts both `--flag value` and `--flag=value` forms, plus the short
#' aliases `-P` (`--params`) and `-O` (`--output-dir`). Unknown flags,
#' missing values, and values that fail conversion to the flag's type raise
#' an error that includes the usage text where helpful.
#'
#' @param arguments Character vector of trailing command-line arguments, as
#'   returned by `commandArgs(trailingOnly = TRUE)`.
#' @return A named list of typed options: `params` and `output_dir`
#'   (character, or `NULL` when not given), `replicates`, `max_cells`,
#'   `seed`, and `progress_updates` (integers), `end_time` and `time_step`
#'   (numerics, or `NULL` to let the library derive them), and `progress` and
#'   `compress_csv` (logicals; only the case-insensitive literals `true` and
#'   `false` are accepted).
#' @section Side effects: `-h`/`--help` prints the usage text and terminates
#'   the R process via `quit()` with status 0.
parse_comparison_cli <- function(arguments){
  options <- list(
    params = NULL,
    output_dir = NULL,
    replicates = 50L,
    end_time = NULL,
    time_step = NULL,
    max_cells = 1000000L,
    seed = 1L,
    progress = TRUE,
    progress_updates = 10L,
    compress_csv = TRUE
  )
  aliases <- c(
    '-P' = 'params', '--params' = 'params',
    '-O' = 'output_dir', '--output-dir' = 'output_dir',
    '--replicates' = 'replicates', '--end-time' = 'end_time',
    '--time-step' = 'time_step', '--max-cells' = 'max_cells',
    '--seed' = 'seed', '--progress' = 'progress',
    '--progress-updates' = 'progress_updates',
    '--compress-csv' = 'compress_csv'
  )
  integer_names <- c('replicates', 'max_cells', 'seed', 'progress_updates')
  numeric_names <- c('end_time', 'time_step')
  logical_names <- c('progress', 'compress_csv')
  index <- 1L
  while(index <= length(arguments)){
    argument <- arguments[index]
    if(argument %in% c('-h', '--help')){
      cat(comparison_cli_usage(), '\n')
      quit(save = 'no', status = 0)
    }
    parts <- strsplit(argument, '=', fixed = TRUE)[[1]]
    flag <- parts[1]
    name <- unname(aliases[flag])
    if(length(name) == 0L || is.na(name)){
      stop(sprintf('Unknown option: %s\n\n%s', flag, comparison_cli_usage()))
    }
    if(length(parts) > 1L){
      raw_value <- paste(parts[-1], collapse = '=')
    } else{
      index <- index + 1L
      if(index > length(arguments)){
        stop(sprintf('Option %s requires a value.', flag))
      }
      raw_value <- arguments[index]
    }
    value <- if(name %in% integer_names){
      suppressWarnings(as.integer(raw_value))
    } else if(name %in% numeric_names){
      suppressWarnings(as.numeric(raw_value))
    } else if(name %in% logical_names){
      normalized <- tolower(raw_value)
      if(!(normalized %in% c('true', 'false'))){
        stop(sprintf('Option %s requires true or false.', flag))
      }
      normalized == 'true'
    } else{
      raw_value
    }
    if(length(value) != 1L || is.na(value)){
      stop(sprintf('Option %s has an invalid value.', flag))
    }
    options[[name]] <- value
    index <- index + 1L
  }
  options
}

options <- parse_comparison_cli(commandArgs(trailingOnly = TRUE))
if(is.null(options$params)){
  stop(sprintf('--params is required.\n\n%s', comparison_cli_usage()))
}
if(!requireNamespace('jsonlite', quietly = TRUE)){
  stop('The jsonlite package is required to read the parameter JSON.')
}
params_path <- normalizePath(options$params, mustWork = TRUE)
params <- jsonlite::fromJSON(params_path, simplifyVector = FALSE)
if(is.null(options$output_dir)){
  parameter_stem <- tools::file_path_sans_ext(basename(params_path))
  options$output_dir <- file.path(
    repo_root,
    'output',
    'engine_comparison',
    parameter_stem
  )
}
dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)

comparison <- compare_population_engines(
  params,
  replicates = options$replicates,
  end_time = options$end_time,
  time_step = options$time_step,
  seed = options$seed,
  max_cells = options$max_cells,
  show_progress = options$progress,
  progress_updates = options$progress_updates
)
comparison_summary <- summarize_engine_comparison(comparison)
cell_type_summary <- summarize_engine_cell_types(comparison)

write_physicell_csv(
  comparison$replicate_summary,
  file.path(options$output_dir, 'engine_replicate_outcomes.csv'),
  row.names = FALSE,
  compress = options$compress_csv
)
write_physicell_csv(
  comparison$cell_type_summary,
  file.path(options$output_dir, 'engine_replicate_cell_types.csv'),
  row.names = FALSE,
  compress = options$compress_csv
)
write_physicell_csv(
  comparison_summary,
  file.path(options$output_dir, 'engine_comparison_summary.csv'),
  row.names = FALSE,
  compress = options$compress_csv
)
write_physicell_csv(
  cell_type_summary,
  file.path(options$output_dir, 'engine_cell_type_summary.csv'),
  row.names = FALSE,
  compress = options$compress_csv
)
write_engine_comparison_plot(
  comparison,
  file.path(options$output_dir, 'engine_comparison.png')
)

cat(sprintf(
  'Compared %d replicates per engine through time %g with time step %g.\n',
  comparison$replicates,
  comparison$end_time,
  comparison$time_step
))
print(
  comparison_summary[
    comparison_summary$metric %in% c('final_cells', 'divisions', 'deaths'),
    c(
      'metric', 'continuous_time_mean', 'time_step_mean',
      'relative_difference_percent'
    )
  ],
  row.names = FALSE
)
cat(sprintf(
  'OUTPUT_DIR=%s\n',
  normalizePath(options$output_dir, mustWork = TRUE)
))
