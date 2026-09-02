#!/usr/bin/env Rscript

# Command-line driver for the lineage-recorder benchmark sweep. It parses the
# grid definition (tree shapes x population seeds x recorder systems x sample
# sizes x integration counts or mitochondrial observation depths) and hands it
# to run_lineage_benchmark(). The sweep is simulation-only: it writes the
# ground-truth Newick tree plus sparse state/character matrices for every
# derived condition so that tree reconstruction can be run separately. The
# script sources prime_editing.R, physicell_lineage.R, gillespie_lineage.R,
# physicell_mito.R and lineage_benchmark.R from its own directory, so it must
# stay alongside them in the repository root; the Matrix and jsonlite packages
# are required at run time.
#
# Command-line arguments (every long flag also accepts --flag=value):
#   -O, --output-dir PATH         Benchmark root. A relative path is resolved
#                                 against the repository root, not the current
#                                 working directory. Default:
#                                 output/lineage_benchmark_<YYYYmmdd_HHMMSS>.
#   --shapes LIST                 Comma list of population shapes; the sweep
#                                 accepts balanced, comb, neutral, turnover and
#                                 hierarchical. Default:
#                                 balanced,comb,neutral,hierarchical.
#   --tree-sizes LIST             Comma list of nested sample sizes.
#                                 Default: 250,1000,2000,5000.
#   --integration-counts LIST     Comma list of barcode integration counts used
#                                 by the non-mitochondrial recorders.
#                                 Default: 1,2,5,10,20.
#   --mt-observation-depths LIST  Comma list of sampled mitochondrial genomes
#                                 per cell. Default: the integration counts.
#   --mt-genomes-per-cell N       Biological mtDNA copies simulated per cell;
#                                 must be at least the largest observation
#                                 depth. Default: 32.
#   --systems LIST                Comma list of recorder systems: baseline,
#                                 prime, palincode, mitochondrial.
#                                 Default: all four.
#   --seeds LIST                  Population seeds as a comma list or an
#                                 inclusive A:B range. Default: 1:10.
#   --synthetic-duration N        Simulated duration for the balanced and comb
#                                 shapes; must be positive. Default: 10.
#   --write-dense-csv BOOL        Also write dense CSV copies of the per-
#                                 condition matrices. Default: false.
#   --resume BOOL                 Reuse a matching non-empty output directory
#                                 instead of failing. Default: true.
#   --progress BOOL               Log population and recording progress.
#                                 Default: true.
#   -h, --help                    Print usage and exit with status 0.
#
# Expected inputs: none on disk for a fresh run - the populations, recorders
# and conditions are all simulated. When --resume is true an existing output
# directory is read back: benchmark_settings.rds must match the requested grid
# exactly, and completed population.rds / sample_sets.rds / recording_complete
# / condition_complete artifacts are loaded instead of resimulated.
#
# Outputs, written under <output-dir>. CSV files are gzip-compressed, so each
# .csv path below is written as .csv.gz:
#   <output-dir>/benchmark_settings.rds and benchmark_settings.json
#   <output-dir>/shape_<shape>/seed_<NNNNNN>/population/
#       population.rds, population_params.json, lineage_nodes.csv,
#       terminal_cells.csv, physicell_lineage_full.nwk,
#       physicell_lineage_sampled.nwk, division_events.csv
#   <output-dir>/shape_<shape>/seed_<NNNNNN>/samples/
#       sample_sets.rds and n_<size>/{sample_cells.csv, ground_truth_tree.nwk}
#   <output-dir>/shape_<shape>/seed_<NNNNNN>/recorder_<system>/
#       full_<NN>_integrations/ (or full_<NN>_genomes_per_cell/ for the
#       mitochondrial system): the recorder simulation at maximum depth, with
#       recording_model.rds, recording_params.json, recording_type.txt,
#       recording_complete.txt and the system-specific matrices/manifests
#       written by write_physicell_recording_outputs() or
#       write_physicell_mito_outputs()
#   <output-dir>/shape_<shape>/seed_<NNNNNN>/recorder_<system>/conditions/
#       n_<size>/k_<NN>/ (or depth_<NN>/ for the mitochondrial system):
#       recording_state_matrix_sparse.rds,
#       recording_character_matrix_sparse.rds,
#       recording_logical_target_matrix_sparse.rds, target_manifest.csv,
#       ground_truth_tree.nwk, a copy of sample_cells.csv,
#       condition_manifest.csv, condition_complete.txt, plus
#       mitochondrial_variant_fraction_matrix_sparse.rds,
#       mitochondrial_binary_variant_matrix_sparse.rds and
#       mitochondrial_variant_manifest.csv for the mitochondrial system, and
#       dense recording_*_matrix.csv copies when --write-dense-csv is true
#   <output-dir>/benchmark_manifest.csv   one row per derived condition
#   <output-dir>/benchmark_summary.csv    simulation counts and grid maxima
#   <output-dir>/benchmark_complete.txt   completion stamp
#
# Standard output ends with the machine-readable line
# BENCHMARK_OUTPUT_DIR=<path>.

script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- normalizePath(
  sub('^--file=', '', script_argument),
  mustWork = TRUE
)
repo_root <- dirname(script_path)
source(file.path(repo_root, 'prime_editing.R'))
source(file.path(repo_root, 'physicell_lineage.R'))
source(file.path(repo_root, 'gillespie_lineage.R'))
source(file.path(repo_root, 'physicell_mito.R'))
source(file.path(repo_root, 'lineage_benchmark.R'))

#' Build the command-line help text for the lineage benchmark sweep
#'
#' Assembles the usage block shown by `-h`/`--help` and embedded in parse
#' errors, including the note that the script simulates only and leaves tree
#' reconstruction to a separate step.
#'
#' @return A single character string with newline-separated usage lines.
lineage_benchmark_cli_usage <- function(){
  paste(
    'Usage:',
    '  Rscript run_lineage_benchmark.R [options]',
    '',
    'Options:',
    '  -O, --output-dir PATH          Default: output/lineage_benchmark_<timestamp>',
    '      --shapes LIST              Default: balanced,comb,neutral,hierarchical',
    '      --tree-sizes LIST          Default: 250,1000,2000,5000',
    '      --integration-counts LIST  Default: 1,2,5,10,20',
    '      --mt-observation-depths LIST  Default: same as integration counts',
    '      --mt-genomes-per-cell N    Biological mtDNA copies/cell; default: 32',
    '      --systems LIST             Default: baseline,prime,palincode,mitochondrial',
    '      --seeds LIST               Comma list or A:B range; default: 1:10',
    '      --synthetic-duration N     Balanced/comb duration; default: 10',
    '      --write-dense-csv BOOL     Also write dense CSV matrices; default: false',
    '      --resume BOOL              Resume a matching output directory; default: true',
    '      --progress BOOL            Show population/recording progress; default: true',
    '  -h, --help',
    '',
    'The script performs simulation only. Every derived condition contains a',
    'ground-truth Newick tree and sparse state/character matrices for separate',
    'tree reconstruction.',
    sep = '\n'
  )
}

#' Validate a logical command-line value
#'
#' Accepts only the literal words `true` and `false`, case-insensitively.
#'
#' @param value The raw command-line value; coerced with `as.character()`.
#' @param flag The flag name used in the error message, e.g. `--resume`.
#' @return `TRUE` or `FALSE`; anything else raises an error naming `flag`.
lineage_benchmark_cli_boolean <- function(value, flag){
  normalized <- tolower(as.character(value))
  if(length(normalized) != 1L || !(normalized %in% c('true', 'false'))){
    stop(sprintf('%s requires true or false.', flag))
  }
  normalized == 'true'
}

#' Expand a seed specification into an integer vector
#'
#' Accepts either an inclusive `A:B` range of non-negative integers, which is
#' expanded with `seq.int()`, or a comma-delimited list, which is validated by
#' `lineage_benchmark_integer_vector()` with a minimum of 0 (so the values must
#' be unique, finite, whole and non-negative, and come back sorted).
#'
#' @param value A length-one character value such as `"1:10"` or `"1,4,9"`.
#' @param name The option name used in error messages, e.g. `--seeds`.
#' @return An integer vector of seeds.
lineage_benchmark_cli_sequence <- function(value, name){
  value <- trimws(as.character(value))
  if(length(value) != 1L || !nzchar(value)){
    stop(sprintf('%s requires a non-empty integer list.', name))
  }
  if(grepl('^[0-9]+:[0-9]+$', value)){
    endpoints <- as.integer(strsplit(value, ':', fixed = TRUE)[[1]])
    return(seq.int(endpoints[1], endpoints[2]))
  }
  lineage_benchmark_integer_vector(value, name, minimum = 0L)
}

#' Parse the lineage benchmark command line
#'
#' Walks the argument vector, accepting either `--flag value` or `--flag=value`
#' (the first `=` splits the token, so values may themselves contain `=`), and
#' maps recognised flags to option names through an alias table. Only three
#' option classes are converted here: `--synthetic-duration` becomes a positive
#' finite number, `--mt-genomes-per-cell` becomes a single validated integer,
#' and the logical flags are resolved by `lineage_benchmark_cli_boolean()`.
#' Every other option is retained as its raw string and expanded later by the
#' driver code, which splits the comma lists and validates the integer vectors.
#' An unknown flag or a flag with no following value raises an error carrying
#' the usage text.
#'
#' @param arguments Character vector of trailing command-line arguments, as
#'   returned by `commandArgs(trailingOnly = TRUE)`.
#' @return A named list of 12 options: `output_dir` (`NULL` until supplied),
#'   `shapes`, `tree_sizes`, `integration_counts`, `mt_observation_depths`
#'   (`NULL` means "use the integration counts"), `mt_genomes_per_cell`,
#'   `systems`, `seeds`, `synthetic_duration`, `write_dense_csv`, `resume`, and
#'   `progress`.
#' @section Side effects: On `-h` or `--help` the usage text is printed to
#'   standard output and the R session terminates with status 0, so the
#'   function does not return in that case.
parse_lineage_benchmark_cli <- function(arguments){
  options <- list(
    output_dir = NULL,
    shapes = 'balanced,comb,neutral,hierarchical',
    tree_sizes = '250,1000,2000,5000',
    integration_counts = '1,2,5,10,20',
    mt_observation_depths = NULL,
    mt_genomes_per_cell = 32L,
    systems = 'baseline,prime,palincode,mitochondrial',
    seeds = '1:10',
    synthetic_duration = 10,
    write_dense_csv = FALSE,
    resume = TRUE,
    progress = TRUE
  )
  aliases <- c(
    '-O' = 'output_dir', '--output-dir' = 'output_dir',
    '--shapes' = 'shapes', '--tree-sizes' = 'tree_sizes',
    '--integration-counts' = 'integration_counts',
    '--mt-observation-depths' = 'mt_observation_depths',
    '--mt-genomes-per-cell' = 'mt_genomes_per_cell',
    '--systems' = 'systems', '--seeds' = 'seeds',
    '--synthetic-duration' = 'synthetic_duration',
    '--write-dense-csv' = 'write_dense_csv', '--resume' = 'resume',
    '--progress' = 'progress'
  )
  index <- 1L
  while(index <= length(arguments)){
    argument <- arguments[index]
    if(argument %in% c('-h', '--help')){
      cat(lineage_benchmark_cli_usage(), '\n')
      quit(save = 'no', status = 0)
    }
    parts <- strsplit(argument, '=', fixed = TRUE)[[1]]
    flag <- parts[1]
    option_name <- unname(aliases[flag])
    if(length(option_name) == 0L || is.na(option_name)){
      stop(sprintf(
        'Unknown option: %s\n\n%s',
        flag,
        lineage_benchmark_cli_usage()
      ))
    }
    if(length(parts) > 1L){
      value <- paste(parts[-1L], collapse = '=')
    } else{
      index <- index + 1L
      if(index > length(arguments)){
        stop(sprintf('%s requires a value.', flag))
      }
      value <- arguments[index]
    }
    if(option_name == 'synthetic_duration'){
      value <- suppressWarnings(as.numeric(value))
      if(length(value) != 1L || !is.finite(value) || value <= 0){
        stop('--synthetic-duration requires one positive number.')
      }
    } else if(option_name == 'mt_genomes_per_cell'){
      value <- lineage_benchmark_integer_vector(
        value,
        '--mt-genomes-per-cell'
      )
      if(length(value) != 1L){
        stop('--mt-genomes-per-cell requires one positive integer.')
      }
    } else if(option_name %in% c(
      'write_dense_csv', 'resume', 'progress'
    )){
      value <- lineage_benchmark_cli_boolean(value, flag)
    }
    options[[option_name]] <- value
    index <- index + 1L
  }
  options
}

options <- parse_lineage_benchmark_cli(commandArgs(trailingOnly = TRUE))
if(is.null(options$output_dir)){
  options$output_dir <- file.path(
    repo_root,
    'output',
    paste0('lineage_benchmark_', format(Sys.time(), '%Y%m%d_%H%M%S'))
  )
} else if(!grepl('^(/|[A-Za-z]:[/\\])', options$output_dir)){
  options$output_dir <- file.path(repo_root, options$output_dir)
}

result <- run_lineage_benchmark(
  output_dir = options$output_dir,
  shapes = trimws(strsplit(options$shapes, ',', fixed = TRUE)[[1]]),
  tree_sizes = lineage_benchmark_integer_vector(
    options$tree_sizes,
    '--tree-sizes'
  ),
  integration_counts = lineage_benchmark_integer_vector(
    options$integration_counts,
    '--integration-counts'
  ),
  mt_observation_depths = if(is.null(options$mt_observation_depths)){
    NULL
  } else{
    lineage_benchmark_integer_vector(
      options$mt_observation_depths,
      '--mt-observation-depths'
    )
  },
  mt_genomes_per_cell = options$mt_genomes_per_cell,
  systems = trimws(strsplit(options$systems, ',', fixed = TRUE)[[1]]),
  seeds = lineage_benchmark_cli_sequence(options$seeds, '--seeds'),
  synthetic_duration = options$synthetic_duration,
  write_dense_csv = options$write_dense_csv,
  resume = options$resume,
  show_progress = options$progress
)
cat(sprintf('BENCHMARK_OUTPUT_DIR=%s\n', result$output_dir))
