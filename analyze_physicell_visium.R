#!/usr/bin/env Rscript

# Command-line driver for the simulated Visium spatial-lineage analysis of a
# single completed PhysiCell run. It validates the sectioning, spot-array,
# pair-sampling and null-model options, delegates the work to
# run_physicell_visium_analysis(), and prints the per-slice occupancy and
# spatial-versus-lineage correlation summary. The script sources
# physicell_lineage.R and physicell_visium.R from its own directory, so it must
# stay alongside them in the repository root.
#
# Command-line arguments (every long flag also accepts --flag=value):
#   --run-dir PATH             Required. Completed run directory produced by
#                              run_physicell_10000_pipeline.sh or a sibling
#                              pipeline script.
#   -O, --output-dir PATH      Analysis destination. A relative path is
#                              resolved against the repository root, not the
#                              current working directory.
#                              Default: <run-dir>/visium_spatial.
#   --slice-offsets LIST       Comma-separated section-plane offsets in microns
#                              from the live-cell centroid; one slice is
#                              analyzed per value. Default: 0.
#   --plane-normal X,Y,Z       Section-plane normal vector. Default: 0,0,1.
#   --section-thickness UM     Finite section thickness in microns. Default: 5.
#   --spot-diameter UM         Capture-spot diameter in microns. Default: 55.
#   --spot-pitch UM            Center-to-center spot spacing in microns.
#                              Default: 100.
#   --randomize-alignment BOOL Randomize the array rotation and translation per
#                              slice (true/false). Default: true.
#   --distance-breaks LIST     Increasing spatial-distance bin edges in
#                              microns; Inf is accepted as the final edge.
#                              Default: 0,20,40,60,80,100,150,200,Inf.
#   --recent-mrca-hours LIST   Recent-MRCA age thresholds in hours; must be
#                              non-negative. Default: 6,12.
#   --permutations N           Spatial tip-label null replicates. Default: 20.
#   --max-pairs-per-bin N      Cell pairs saved and evaluated per distance bin.
#                              Default: 100000.
#   --max-exact-pairs N        Enumerate every unordered cell pair below this
#                              count. Default: 2000000.
#   --max-candidate-pairs N    Uniform pair draws used above the exact limit.
#                              Default: 2000000.
#   --time-units-per-hour N    Simulation time units per hour; PhysiCell emits
#                              minutes. Default: 60.
#   --seed N                   Analysis seed; each slice derives its own seed
#                              from it. Default: 1.
#   --compress-csv BOOL        Write gzip-compressed CSV. Default: true.
#   --progress BOOL            Emit per-slice progress messages. Default: true.
#   -h, --help                 Print usage and exit with status 0.
#
# Expected inputs, all read from --run-dir:
#   <run-dir>/physicell_build/output/lineage_table.csv  live cell ID, x, y, z
#   <run-dir>/lineage_recording/terminal_cells.csv      sample_id, physicell_id,
#                                                       node_id
#   <run-dir>/lineage_recording/lineage_nodes.csv       event-resolved lineage
#                                                       nodes
#
# Outputs, written under <output-dir>. Every .csv path below gains a .gz suffix
# while --compress-csv is true, and <NNN> is the 1-based slice index:
#   <output-dir>/slice_<NNN>/visium_spots.csv
#   <output-dir>/slice_<NNN>/slice_cells.csv
#   <output-dir>/slice_<NNN>/spot_cell_membership.csv
#   <output-dir>/slice_<NNN>/cell_pair_lineage_sample.csv
#   <output-dir>/slice_<NNN>/cell_distance_summary.csv
#   <output-dir>/slice_<NNN>/spot_pair_lineage.csv
#   <output-dir>/slice_<NNN>/spot_distance_summary.csv
#   <output-dir>/slice_<NNN>/slice_summary.csv
#   <output-dir>/slice_<NNN>/slice_manifest.csv
#   <output-dir>/slice_<NNN>/visium_slice.png
#   <output-dir>/slice_<NNN>/cell_lineage_correlogram.png
#   <output-dir>/cell_distance_summary.csv   all slices concatenated
#   <output-dir>/spot_distance_summary.csv   all slices concatenated
#   <output-dir>/slice_summary.csv           one row per slice
#   <output-dir>/analysis_manifest.csv       resolved analysis settings
#   <output-dir>/analysis_complete.txt       completion timestamp line
#
# Standard output ends with the machine-readable line OUTPUT_DIR=<path>.

script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- normalizePath(sub('^--file=', '', script_argument), mustWork = TRUE)
repo_root <- dirname(script_path)
source(file.path(repo_root, 'physicell_lineage.R'))
source(file.path(repo_root, 'physicell_visium.R'))

#' Build the command-line help text for the single-run Visium analysis
#'
#' Assembles the usage block shown by `-h`/`--help` and embedded in parse
#' errors. The listed defaults mirror those in `parse_physicell_visium_cli()`.
#'
#' @return A single character string with newline-separated usage lines.
physicell_visium_cli_usage <- function(){
  paste(
    'Usage:',
    '  Rscript analyze_physicell_visium.R --run-dir PHYSICELL_RUN [options]',
    '',
    'Required:',
    '      --run-dir PATH             Completed run from run_physicell_10000_pipeline.sh',
    '',
    'Options:',
    '  -O, --output-dir PATH         Default: <run-dir>/visium_spatial',
    '      --slice-offsets LIST       Plane offsets from tumor center; default: 0',
    '      --plane-normal X,Y,Z      Default: 0,0,1',
    '      --section-thickness UM    Default: 5',
    '      --spot-diameter UM        Conventional Visium default: 55',
    '      --spot-pitch UM           Conventional Visium default: 100',
    '      --randomize-alignment BOOL Random grid rotation/translation; default: true',
    '      --distance-breaks LIST    Default: 0,20,40,60,80,100,150,200,Inf',
    '      --recent-mrca-hours LIST  Default: 6,12',
    '      --permutations N          Spatial tip-label null replicates; default: 20',
    '      --max-pairs-per-bin N     Saved/evaluated pairs per distance bin; default: 100000',
    '      --max-exact-pairs N       Enumerate all pairs below this limit; default: 2000000',
    '      --max-candidate-pairs N   Uniform pair draws above exact limit; default: 2000000',
    '      --time-units-per-hour N   PhysiCell uses 60 minutes/hour; default: 60',
    '      --seed N                  Analysis seed; default: 1',
    '      --compress-csv BOOL       Default: true',
    '      --progress BOOL           Default: true',
    '  -h, --help',
    sep = '\n'
  )
}

#' Parse and type-check the Visium analysis command line
#'
#' Walks the argument vector, accepting either `--flag value` or `--flag=value`
#' (the first `=` splits the token, so values may themselves contain `=`).
#' Recognised flags are mapped to option names through an alias table and each
#' value is coerced by option class: integer for the permutation/pair/seed
#' options, numeric for the geometry and time-scale options, and `true`/`false`
#' for the logical options. Comma-delimited list options (`--slice-offsets`,
#' `--plane-normal`, `--distance-breaks`, `--recent-mrca-hours`) are left as raw
#' strings here and parsed downstream by `run_physicell_visium_analysis()`.
#' Any unknown flag, missing value, non-`true`/`false` logical, or value that
#' coerces to `NA` raises an error carrying the usage text.
#'
#' @param arguments Character vector of trailing command-line arguments, as
#'   returned by `commandArgs(trailingOnly = TRUE)`.
#' @return A named list of 18 options: `run_dir`, `output_dir` (both `NULL`
#'   until supplied), `slice_offsets`, `plane_normal`, `section_thickness`,
#'   `spot_diameter`, `spot_pitch`, `randomize_alignment`, `distance_breaks`,
#'   `recent_mrca_hours`, `permutations`, `max_pairs_per_bin`,
#'   `max_exact_pairs`, `max_candidate_pairs`, `time_units_per_hour`, `seed`,
#'   `compress_csv`, and `progress`.
#' @section Side effects: On `-h` or `--help` the usage text is printed to
#'   standard output and the R session terminates with status 0, so the
#'   function does not return in that case.
parse_physicell_visium_cli <- function(arguments){
  options <- list(
    run_dir = NULL,
    output_dir = NULL,
    slice_offsets = '0',
    plane_normal = '0,0,1',
    section_thickness = 5,
    spot_diameter = 55,
    spot_pitch = 100,
    randomize_alignment = TRUE,
    distance_breaks = '0,20,40,60,80,100,150,200,Inf',
    recent_mrca_hours = '6,12',
    permutations = 20L,
    max_pairs_per_bin = 100000L,
    max_exact_pairs = 2000000L,
    max_candidate_pairs = 2000000L,
    time_units_per_hour = 60,
    seed = 1L,
    compress_csv = TRUE,
    progress = TRUE
  )
  aliases <- c(
    '--run-dir' = 'run_dir',
    '-O' = 'output_dir', '--output-dir' = 'output_dir',
    '--slice-offsets' = 'slice_offsets',
    '--plane-normal' = 'plane_normal',
    '--section-thickness' = 'section_thickness',
    '--spot-diameter' = 'spot_diameter',
    '--spot-pitch' = 'spot_pitch',
    '--randomize-alignment' = 'randomize_alignment',
    '--distance-breaks' = 'distance_breaks',
    '--recent-mrca-hours' = 'recent_mrca_hours',
    '--permutations' = 'permutations',
    '--max-pairs-per-bin' = 'max_pairs_per_bin',
    '--max-exact-pairs' = 'max_exact_pairs',
    '--max-candidate-pairs' = 'max_candidate_pairs',
    '--time-units-per-hour' = 'time_units_per_hour',
    '--seed' = 'seed',
    '--compress-csv' = 'compress_csv',
    '--progress' = 'progress'
  )
  integer_options <- c(
    'permutations', 'max_pairs_per_bin', 'max_exact_pairs',
    'max_candidate_pairs', 'seed'
  )
  numeric_options <- c(
    'section_thickness', 'spot_diameter', 'spot_pitch', 'time_units_per_hour'
  )
  logical_options <- c('randomize_alignment', 'compress_csv', 'progress')
  index <- 1L
  while(index <= length(arguments)){
    argument <- arguments[index]
    if(argument %in% c('-h', '--help')){
      cat(physicell_visium_cli_usage(), '\n')
      quit(save = 'no', status = 0)
    }
    parts <- strsplit(argument, '=', fixed = TRUE)[[1]]
    flag <- parts[1]
    option_name <- unname(aliases[flag])
    if(length(option_name) == 0L || is.na(option_name)){
      stop(sprintf('Unknown option: %s\n\n%s', flag, physicell_visium_cli_usage()))
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
    if(option_name %in% integer_options){
      value <- suppressWarnings(as.integer(raw_value))
    } else if(option_name %in% numeric_options){
      value <- suppressWarnings(as.numeric(raw_value))
    } else if(option_name %in% logical_options){
      normalized <- tolower(raw_value)
      if(!(normalized %in% c('true', 'false'))){
        stop(sprintf('Option %s requires true or false.', flag))
      }
      value <- normalized == 'true'
    } else{
      value <- raw_value
    }
    if(length(value) != 1L || is.na(value)){
      stop(sprintf('Option %s has an invalid value.', flag))
    }
    options[[option_name]] <- value
    index <- index + 1L
  }
  options
}

options <- parse_physicell_visium_cli(commandArgs(trailingOnly = TRUE))
if(is.null(options$run_dir)){
  stop(sprintf('--run-dir is required.\n\n%s', physicell_visium_cli_usage()))
}
run_dir <- normalizePath(options$run_dir, mustWork = TRUE)
if(is.null(options$output_dir)){
  options$output_dir <- file.path(run_dir, 'visium_spatial')
} else if(!grepl('^/', options$output_dir)){
  options$output_dir <- file.path(repo_root, options$output_dir)
}

result <- run_physicell_visium_analysis(
  run_dir = run_dir,
  output_dir = options$output_dir,
  slice_offsets = options$slice_offsets,
  plane_normal = options$plane_normal,
  section_thickness = options$section_thickness,
  spot_diameter = options$spot_diameter,
  spot_pitch = options$spot_pitch,
  randomize_array_alignment = options$randomize_alignment,
  distance_breaks = options$distance_breaks,
  recent_threshold_hours = options$recent_mrca_hours,
  permutations = options$permutations,
  max_pairs_per_bin = options$max_pairs_per_bin,
  max_exact_pairs = options$max_exact_pairs,
  max_candidate_pairs = options$max_candidate_pairs,
  time_units_per_hour = options$time_units_per_hour,
  seed = options$seed,
  compress_csv = options$compress_csv,
  show_progress = options$progress
)

cat(sprintf(
  'Analyzed %d slice(s); output written to %s\n',
  nrow(result$slice_summary),
  result$output_dir
))
print(
  result$slice_summary[, c(
    'slice_id', 'section_cells', 'captured_cells', 'occupied_spots',
    'mean_cells_per_occupied_spot', 'spearman_spatial_vs_mrca'
  )],
  row.names = FALSE
)
cat(sprintf('OUTPUT_DIR=%s\n', result$output_dir))

