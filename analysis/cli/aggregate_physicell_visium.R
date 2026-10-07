#!/usr/bin/env Rscript

# Command-line driver that pools the Visium-style spatial lineage summaries
# produced by repeated single-replicate runs into one batch-level result.
#
# Pipeline position: run_physicell_visium_replicates.sh drives
# simulate_physicell_lineage.R and analyze_physicell_visium.R once per
# replicate, leaving a replicate_*/visium_spatial directory per tumor. This
# script is the final step: it reads those per-replicate summaries, pools
# slices within a tumor, and reports across-tumor statistics.
#
# Loads physicell_lineage.R (CSV read/write helpers) and physicell_visium.R,
# which supplies `aggregate_physicell_visium_results()`, from the directory
# holding this script.
#
# Command-line arguments
#       --batch-dir PATH      Required. Searched recursively for any
#                             visium_spatial directory holding a
#                             cell_distance_summary table; the directory that
#                             contains visium_spatial names the replicate.
#   -O, --output-dir PATH     Destination directory; a relative path is
#                             resolved against the repository root. Default
#                             <batch-dir>/aggregate.
#       --compress-csv BOOL   Write CSV tables as .csv.gz. Default true.
#   -h, --help                Print the usage block and exit with status 0.
#
# Inputs read: every <batch-dir>/**/visium_spatial/cell_distance_summary.csv
# (or .csv.gz) together with the spot_distance_summary.csv beside it. At least
# one such cell-distance file must exist.
#
# Outputs, all under the resolved output directory, with .csv written as
# .csv.gz unless --compress-csv false is passed:
#   all_slice_distance_summaries.csv, replicate_distance_summary.csv,
#   aggregate_distance_summary.csv, all_slice_spot_distance_summaries.csv,
#   replicate_spot_distance_summary.csv, aggregate_spot_distance_summary.csv,
#   aggregate_lineage_correlogram.png, and aggregate_complete.txt.
#
# The last stdout line is OUTPUT_DIR=<path>, which the wrapper scripts parse.

script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- normalizePath(sub('^--file=', '', script_argument), mustWork = TRUE)
# The repository root is 2 levels above this script; every path this
# script resolves, including its outputs, is relative to that root.
repo_root <- normalizePath(file.path(dirname(script_path), '..', '..'))
source(file.path(repo_root, 'load_origin.R'))

#' Build this script's command-line usage text
#'
#' @return One string containing the newline-separated usage block: the
#'   invocation line and the batch-directory, output-directory, compression,
#'   and help options.
visium_aggregate_usage <- function(){
  paste(
    'Usage:',
    '  Rscript aggregate_physicell_visium.R --batch-dir PATH [options]',
    '',
    'Options:',
    '      --batch-dir PATH       Directory containing replicate_*/visium_spatial',
    '  -O, --output-dir PATH     Default: <batch-dir>/aggregate',
    '      --compress-csv BOOL   Default: true',
    '  -h, --help',
    sep = '\n'
  )
}

# ---- Parse the command line ----
# Accepts both `--flag value` and `--flag=value`. Unknown flags and missing
# values are fatal, and `-h`/`--help` prints the usage block and exits 0.

arguments <- commandArgs(trailingOnly = TRUE)
options <- list(batch_dir = NULL, output_dir = NULL, compress_csv = TRUE)
aliases <- c(
  '--batch-dir' = 'batch_dir',
  '-O' = 'output_dir', '--output-dir' = 'output_dir',
  '--compress-csv' = 'compress_csv'
)
index <- 1L
while(index <= length(arguments)){
  argument <- arguments[index]
  if(argument %in% c('-h', '--help')){
    cat(visium_aggregate_usage(), '\n')
    quit(save = 'no', status = 0)
  }
  parts <- strsplit(argument, '=', fixed = TRUE)[[1]]
  flag <- parts[1]
  option_name <- unname(aliases[flag])
  if(length(option_name) == 0L || is.na(option_name)){
    stop(sprintf('Unknown option: %s\n\n%s', flag, visium_aggregate_usage()))
  }
  if(length(parts) > 1L){
    value <- paste(parts[-1], collapse = '=')
  } else{
    index <- index + 1L
    if(index > length(arguments)) stop(sprintf('Option %s requires a value.', flag))
    value <- arguments[index]
  }
  if(option_name == 'compress_csv'){
    value <- tolower(value)
    if(!(value %in% c('true', 'false'))){
      stop('--compress-csv requires true or false.')
    }
    value <- value == 'true'
  }
  options[[option_name]] <- value
  index <- index + 1L
}
# ---- Resolve paths, aggregate the replicates, and report ----

if(is.null(options$batch_dir)){
  stop(sprintf('--batch-dir is required.\n\n%s', visium_aggregate_usage()))
}
batch_dir <- normalizePath(options$batch_dir, mustWork = TRUE)
if(is.null(options$output_dir)){
  options$output_dir <- file.path(batch_dir, 'aggregate')
} else if(!grepl('^/', options$output_dir)){
  options$output_dir <- file.path(repo_root, options$output_dir)
}
result <- aggregate_physicell_visium_results(
  batch_dir,
  output_dir = options$output_dir,
  compress_csv = options$compress_csv
)
cat(sprintf(
  'Aggregated %d replicate-distance rows into %s\n',
  nrow(result$replicate_summary),
  result$output_dir
))
cat(sprintf('OUTPUT_DIR=%s\n', result$output_dir))

