#!/usr/bin/env Rscript

# Standalone command-line entry point for the Gillespie lineage engine.
#
# Resolves the repository root from this script's own --file= path and sources
# prime_editing.R, physicell_lineage.R, physicell_mito.R, ecdna_lineage.R,
# gillespie_lineage.R, and gillespie_pipeline.R from it, then parses the options
# below and hands the parsed parameter JSON to run_gillespie_lineage_pipeline().
# sim5_code.R reaches the same pipeline when the JSON sets
# "simulation_engine": "gillespie"; this script is the direct route.

script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- normalizePath(sub('^--file=', '', script_argument), mustWork = TRUE)
# The repository root is 3 levels above this script; every path this
# script resolves, including its outputs, is relative to that root.
repo_root <- normalizePath(file.path(dirname(script_path), '..', '..', '..'))
source(file.path(repo_root, 'load_origin.R'))

#' Return the standalone Gillespie command-line help text
#'
#' @return One newline-joined string listing the required `--params` option and
#'   every optional flag together with its default.
gillespie_cli_usage <- function(){
  paste(
    'Usage:',
    '  Rscript simulate_gillespie_lineage.R --params PARAMS.json [options]',
    '',
    'Required:',
    '  -P, --params PATH              remote_mito JSON parameter file',
    '',
    'Options:',
    '  -O, --output-dir PATH         Default: output/gillespie/<JSON stem>',
    '      --end-time NUMBER         Default: max(sim_length)',
    '      --modalities LIST         barcode, mitochondrial, ecDNA, or both',
    '      --max-cells N             Active-cell safety limit; default: 1000000',
    '      --num-integrations N      Integrated barcode recorders',
    '      --founder-label-sites N   Stable founder labels; default: 0',
    '      --mt-genomes-per-cell N   Mitochondrial bottleneck; default: 8',
    '      --write-mt-fasta BOOL     Default: false',
    '      --compress-csv BOOL       Default: true',
    '      --progress BOOL           Default: true',
    '      --progress-updates N      Default: 20',
    '      --seed N                  Default: JSON random_seed',
    '  -h, --help',
    sep = '\n'
  )
}

#' Parse the standalone Gillespie command-line options
#'
#' Accepts both the `--flag value` and `--flag=value` forms (a value may itself
#' contain `=`; only the first is treated as the separator), plus the short
#' aliases `-P` for `--params` and `-O` for `--output-dir`. `-h`/`--help` prints
#' the usage text and exits immediately.
#'
#' @details
#' Every recognized flag and the type its value is coerced to:
#' `-P`/`--params` (character, the remote_mito JSON parameter file) and
#' `-O`/`--output-dir` (character) are the only options with dedicated list
#' slots; `--modalities` (character: `barcode`, `mitochondrial`, `ecDNA`,
#' `both`, or `all`) is also left as text. `--end-time` is numeric.
#' `--max-cells`, `--num-integrations`, `--founder-label-sites`,
#' `--mt-genomes-per-cell`, `--progress-updates`, and `--seed` are integer.
#' `--write-mt-fasta`, `--compress-csv`, and `--progress` accept only `true` or
#' `false`, compared case-insensitively. Defaults are not applied here: an
#' option the caller omits is simply absent, and
#' `run_gillespie_lineage_pipeline()` falls back to the JSON `gillespie` block
#' and then to its own defaults.
#'
#' @param arguments Character vector of trailing command-line arguments, as
#'   returned by `commandArgs(trailingOnly = TRUE)`.
#' @return A named list of parsed options. `params` and `output_dir` are always
#'   present and are `NULL` when not supplied. Stops on an unrecognized flag, a
#'   flag given without a value, a boolean flag given anything but true/false,
#'   or a value that does not coerce to a single non-`NA` value of its type.
#' @section Side effects: `-h` or `--help` prints the usage text and calls
#'   `quit(save = 'no', status = 0)`, ending the R session.
parse_gillespie_cli <- function(arguments){
  options <- list(params = NULL, output_dir = NULL)
  aliases <- c(
    '-P' = 'params', '--params' = 'params',
    '-O' = 'output_dir', '--output-dir' = 'output_dir',
    '--end-time' = 'end_time', '--modalities' = 'modalities',
    '--max-cells' = 'max_cells', '--num-integrations' = 'num_integrations',
    '--founder-label-sites' = 'founder_label_sites',
    '--mt-genomes-per-cell' = 'mt_genomes_per_cell',
    '--write-mt-fasta' = 'write_mt_fasta',
    '--compress-csv' = 'compress_csv', '--progress' = 'progress',
    '--progress-updates' = 'progress_updates', '--seed' = 'seed'
  )
  numeric_names <- c('end_time')
  integer_names <- c(
    'max_cells', 'num_integrations', 'founder_label_sites',
    'mt_genomes_per_cell', 'progress_updates', 'seed'
  )
  logical_names <- c('write_mt_fasta', 'compress_csv', 'progress')
  index <- 1L
  while(index <= length(arguments)){
    argument <- arguments[index]
    if(argument %in% c('-h', '--help')){
      cat(gillespie_cli_usage(), '\n')
      quit(save = 'no', status = 0)
    }
    parts <- strsplit(argument, '=', fixed = TRUE)[[1]]
    flag <- parts[1]
    name <- unname(aliases[flag])
    if(length(name) == 0L || is.na(name)){
      stop(sprintf('Unknown option: %s\n\n%s', flag, gillespie_cli_usage()))
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
    parsed_value <- if(name %in% numeric_names){
      as.numeric(raw_value)
    } else if(name %in% integer_names){
      as.integer(raw_value)
    } else if(name %in% logical_names){
      normalized <- tolower(raw_value)
      if(!(normalized %in% c('true', 'false'))){
        stop(sprintf('Option %s requires true or false.', flag))
      }
      normalized == 'true'
    } else{
      raw_value
    }
    if(length(parsed_value) != 1L || is.na(parsed_value)){
      stop(sprintf('Option %s has an invalid value.', flag))
    }
    options[[name]] <- parsed_value
    index <- index + 1L
  }
  options
}

# Top-level driver. --params is required; without --output-dir the run writes to
# <repo>/output/gillespie/<JSON stem>, created if it does not exist. Every
# parsed option other than params and output_dir is forwarded to
# run_gillespie_lineage_pipeline() as an override of the matching key in the
# JSON "gillespie" block.
#
# Files written into the output directory. Every CSV below gains a .gz suffix
# unless --compress-csv false is given.
#
#   Always: lineage_nodes.csv, terminal_cells.csv, physicell_lineage_full.nwk,
#     physicell_lineage_sampled.nwk, gillespie_population_events.csv,
#     gillespie_division_events.csv, gillespie_checkpoint_summary.csv,
#     gillespie_cell_states.csv, gillespie_run_manifest.csv, and
#     r_timing_summary.csv.
#   Barcode modality: barcode_alleles and barcode_binary_score_matrix (sparse
#     .rds or dense .csv depending on the recorder's compact-output setting),
#     barcode_profiles.rds, barcode_target_layout.csv, mutation_events.csv, and
#     run_manifest.csv; barcode_reference.fasta and barcode_sequences.fasta for
#     a plain barcode; the palincode_* or prime_editing_* matrices and
#     prime_editing_target_manifest.csv for those recorders.
#   Mitochondrial modality: mitochondrial_profiles.rds,
#     mitochondrial_mutation_events.csv, mitochondrial_variant_fractions.csv,
#     mitochondrial_variant_fraction_matrix.rds, mitochondrial_reference.fasta,
#     mitochondrial_manifest.csv, and mitochondrial_sampled_haplotypes.fasta
#     when --write-mt-fasta true.
#   ecDNA modality: ecdna_cell_summary.csv, ecdna_haplotypes.csv,
#     ecdna_mutation_events.csv, ecdna_species_manifest.csv, ecdna_manifest.csv,
#     ecdna_terminal_profiles.rds, and the ecdna_*_matrix_sparse.rds matrices.
#   Barcode or mitochondrial: mutation_event_descendant_matrix_sparse.rds and
#     mutation_event_descendant_manifest.csv.
#   Barcode and mitochondrial together: combined_lineage_feature_matrix.rds and
#     combined_lineage_feature_manifest.csv.
#
# The pipeline also prints a phase timing summary and the GILLESPIE_OUTPUT_DIR=
# and OUTPUT_DIR= lines that the wrapper shell scripts parse.
options <- parse_gillespie_cli(commandArgs(trailingOnly = TRUE))
if(is.null(options$params)){
  stop(sprintf('--params is required.\n\n%s', gillespie_cli_usage()))
}
if(!requireNamespace('jsonlite', quietly = TRUE)){
  stop('The jsonlite package is required to read the parameter JSON.')
}
params_path <- normalizePath(options$params, mustWork = TRUE)
params <- jsonlite::fromJSON(params_path, simplifyVector = FALSE)
if(is.null(options$output_dir)){
  stem <- tools::file_path_sans_ext(basename(params_path))
  options$output_dir <- file.path(repo_root, 'output', 'gillespie', stem)
}
overrides <- options[setdiff(names(options), c('params', 'output_dir'))]
run_gillespie_lineage_pipeline(
  params,
  params_path = params_path,
  output_dir = options$output_dir,
  overrides = overrides
)
