#!/usr/bin/env Rscript

# Command-line driver that renders the parameterized Visium spatial-lineage
# notebook to a standalone HTML report. It runs after analyze_physicell_visium.R
# (or run_physicell_visium_replicates.sh) has produced the analysis CSVs. The
# script copies notebooks/physicell_visium_lineage_distance.Rmd into a fresh
# temporary staging directory and renders from there, which keeps rmarkdown's
# resource paths relative and avoids absolute-output path problems, then copies
# the result to the requested destination. It loads no repository R modules but
# requires the rmarkdown package (and ggplot2, which the notebook itself
# demands).
#
# Command-line arguments (every long flag also accepts --flag=value):
#   --input-path PATH          Required. Either one visium_spatial analysis
#                              directory or a replicate-batch root containing
#                              several of them; the notebook discovers the
#                              analysis directories underneath it.
#   -O, --output-file PATH     Destination HTML file; by default the file
#                              physicell_visium_lineage_distance.html inside
#                              --input-path. A relative path is resolved
#                              against the repository root, not the current
#                              working directory.
#   --max-pairs-per-slice N    Cap on saved cell pairs plotted per slice; must
#                              be positive. Default: 50000.
#   --recent-mrca-hours N      Recent-MRCA enrichment threshold in hours; must
#                              be non-negative. Default: 6.
#   --seed N                   Seed for the notebook's display-only pair
#                              subsampling. Default: 1.
#   -h, --help                 Print usage and exit with status 0.
#
# Expected inputs:
#   <repo-root>/notebooks/physicell_visium_lineage_distance.Rmd  the rendered
#       notebook source, which must sit next to this script's directory
#   <input-path>/.../visium_spatial/cell_distance_summary.csv[.gz] and the
#       companion slice/spot summaries and saved cell-pair tables written by
#       analyze_physicell_visium.R
#
# Parameters passed to the Rmd (via rmarkdown::render(params = ...), matching
# the params block in its YAML header):
#   input_path           the normalized absolute --input-path
#   max_pairs_per_slice  --max-pairs-per-slice
#   recent_mrca_hours    --recent-mrca-hours
#   seed                 --seed
#
# Outputs:
#   <output-file>  the rendered html_notebook document, copied out of the
#                  temporary staging directory; its parent directory is created
#                  if needed
#
# Standard output ends with the machine-readable line NOTEBOOK_OUTPUT=<path>.

script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- normalizePath(sub('^--file=', '', script_argument), mustWork = TRUE)
# The repository root is 2 levels above this script; every path this
# script resolves, including its outputs, is relative to that root.
repo_root <- normalizePath(file.path(dirname(script_path), '..', '..'))

#' Build the command-line help text for the Visium notebook renderer
#'
#' Assembles the usage block shown by `-h`/`--help` and embedded in parse
#' errors. The listed defaults mirror those in
#' `parse_physicell_visium_notebook_cli()`.
#'
#' @return A single character string with newline-separated usage lines.
physicell_visium_notebook_usage <- function(){
  paste(
    'Usage:',
    '  Rscript render_physicell_visium_notebook.R --input-path PATH [options]',
    '',
    'Required:',
    '      --input-path PATH          One visium_spatial directory or replicate batch',
    '',
    'Options:',
    '  -O, --output-file PATH        Default: <input>/physicell_visium_lineage_distance.html',
    '      --max-pairs-per-slice N   Plotting cap; default: 50000',
    '      --recent-mrca-hours N     Enrichment threshold; default: 6',
    '      --seed N                  Plot-pair sampling seed; default: 1',
    '  -h, --help',
    sep = '\n'
  )
}

#' Parse and validate the Visium notebook renderer command line
#'
#' Walks the argument vector, accepting either `--flag value` or `--flag=value`
#' (the first `=` splits the token, so values may themselves contain `=`), maps
#' recognised flags to option names through an alias table, and coerces
#' `--max-pairs-per-slice` and `--seed` to integer and `--recent-mrca-hours` to
#' numeric; the two path options stay as raw strings. An unknown flag, a flag
#' with no following value, or a value that coerces to `NA` raises an error
#' carrying the usage text. After the whole line is consumed the plotting cap is
#' required to be positive and the recent-MRCA threshold non-negative.
#'
#' @param arguments Character vector of trailing command-line arguments, as
#'   returned by `commandArgs(trailingOnly = TRUE)`.
#' @return A named list of five options: `input_path` and `output_file` (both
#'   `NULL` until supplied), `max_pairs_per_slice`, `recent_mrca_hours`, and
#'   `seed`.
#' @section Side effects: On `-h` or `--help` the usage text is printed to
#'   standard output and the R session terminates with status 0, so the
#'   function does not return in that case.
parse_physicell_visium_notebook_cli <- function(arguments){
  options <- list(
    input_path = NULL,
    output_file = NULL,
    max_pairs_per_slice = 50000L,
    recent_mrca_hours = 6,
    seed = 1L
  )
  aliases <- c(
    '--input-path' = 'input_path',
    '-O' = 'output_file', '--output-file' = 'output_file',
    '--max-pairs-per-slice' = 'max_pairs_per_slice',
    '--recent-mrca-hours' = 'recent_mrca_hours',
    '--seed' = 'seed'
  )
  integer_options <- c('max_pairs_per_slice', 'seed')
  numeric_options <- 'recent_mrca_hours'
  index <- 1L
  while(index <= length(arguments)){
    argument <- arguments[index]
    if(argument %in% c('-h', '--help')){
      cat(physicell_visium_notebook_usage(), '\n')
      quit(save = 'no', status = 0)
    }
    parts <- strsplit(argument, '=', fixed = TRUE)[[1]]
    flag <- parts[1]
    option_name <- unname(aliases[flag])
    if(length(option_name) == 0L || is.na(option_name)){
      stop(sprintf(
        'Unknown option: %s\n\n%s',
        flag,
        physicell_visium_notebook_usage()
      ))
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
    value <- if(option_name %in% integer_options){
      suppressWarnings(as.integer(raw_value))
    } else if(option_name %in% numeric_options){
      suppressWarnings(as.numeric(raw_value))
    } else{
      raw_value
    }
    if(length(value) != 1L || is.na(value)){
      stop(sprintf('Option %s has an invalid value.', flag))
    }
    options[[option_name]] <- value
    index <- index + 1L
  }
  if(options$max_pairs_per_slice < 1L){
    stop('--max-pairs-per-slice must be positive.')
  }
  if(options$recent_mrca_hours < 0){
    stop('--recent-mrca-hours must be non-negative.')
  }
  options
}

options <- parse_physicell_visium_notebook_cli(commandArgs(trailingOnly = TRUE))
if(is.null(options$input_path)){
  stop(sprintf(
    '--input-path is required.\n\n%s',
    physicell_visium_notebook_usage()
  ))
}
if(!requireNamespace('rmarkdown', quietly = TRUE)){
  stop('The rmarkdown package is required to render the notebook.')
}
input_path <- normalizePath(options$input_path, mustWork = TRUE)
if(is.null(options$output_file)){
  output_file <- file.path(input_path, 'physicell_visium_lineage_distance.html')
} else if(grepl('^/', options$output_file)){
  output_file <- options$output_file
} else{
  output_file <- file.path(repo_root, options$output_file)
}
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

notebook_source <- file.path(
  repo_root,
  'notebooks',
  'physicell_visium_lineage_distance.Rmd'
)
staging_dir <- tempfile('physicell_visium_notebook_')
dir.create(staging_dir)
on.exit(unlink(staging_dir, recursive = TRUE), add = TRUE)
staged_notebook <- file.path(staging_dir, basename(notebook_source))
if(!file.copy(notebook_source, staged_notebook)){
  stop('Could not stage the Visium notebook for rendering.')
}
old_working_directory <- getwd()
on.exit(setwd(old_working_directory), add = TRUE)
setwd(staging_dir)

rendered <- rmarkdown::render(
  input = basename(staged_notebook),
  output_format = 'html_notebook',
  output_file = basename(output_file),
  params = list(
    input_path = input_path,
    max_pairs_per_slice = options$max_pairs_per_slice,
    recent_mrca_hours = options$recent_mrca_hours,
    seed = options$seed
  ),
  envir = new.env(parent = globalenv()),
  quiet = FALSE
)
if(!file.copy(rendered, output_file, overwrite = TRUE)){
  stop(sprintf('Could not copy the rendered notebook to %s.', output_file))
}

cat(sprintf(
  'NOTEBOOK_OUTPUT=%s\n',
  normalizePath(output_file, mustWork = TRUE)
))
