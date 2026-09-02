#!/usr/bin/env Rscript

# Command-line driver that turns a finished PhysiCell run into simulated
# lineage-recorder data.
#
# Pipeline position: PhysiCell writes a division-event CSV (and optionally a
# live-cell CSV); this script replays those events into an event-resolved
# lineage table, writes the ground-truth lineage, and then layers the requested
# recording modalities (barcode, mitochondrial, ecDNA) on top of that lineage.
# Downstream analysis (physicell_visium.R, the scdesign3 helpers, the benchmark
# scripts) consumes the tables written here.
#
# Loads prime_editing.R, physicell_lineage.R, physicell_mito.R, and
# ecdna_lineage.R from the directory holding this script. The jsonlite package
# is required to read the parameter JSON; the Matrix package is additionally
# required when both the barcode and mitochondrial modalities are requested.
#
# Command-line arguments
#   -L, --lineage PATH        Required. PhysiCell division CSV with time,
#                             parent_ID, and daughter_ID columns.
#   -P, --params PATH         Required. remote_mito JSON parameter file.
#   -O, --output-dir PATH     Destination directory; created if absent.
#                             Default <repo>/output/physicell/<lineage stem>.
#       --end-time NUMBER     Sampling time of the terminal cells. Default is
#                             the maximum sim_length found in the JSON.
#       --founder-time NUMBER Birth time of the founder segments. Default 0.
#       --founders PATH       Optional CSV listing every day-zero founder ID.
#       --live-cells PATH     Optional PhysiCell live-cell CSV with an ID
#                             column; restricts which terminal cells are
#                             sampled and simulated.
#       --cell-type NAME      Cell type whose rates are read from the JSON.
#                             Default founder_cell_type.
#       --editing-state STATE auto, induced, or uninduced. Default auto, which
#                             splits branches at the model's induction time.
#       --modalities LIST     Comma-separated selection of lineage, barcode,
#                             mitochondrial, ecDNA, both (barcode plus
#                             mitochondrial), or all. Default barcode; the
#                             alias mito maps to mitochondrial.
#       --num-integrations N  Barcode integrations per cell. Default is the
#                             maximum max_bc_ints_per_cell in the JSON.
#       --founder-label-sites N  Stable allele-coded founder barcode sites.
#                             Default 0.
#       --mt-genomes-per-cell N  Fixed mitochondrial bottleneck size.
#                             Default 8.
#       --write-mt-fasta BOOL Write one sampled mt haplotype per cell.
#                             Default false.
#       --compress-csv BOOL   Write CSV tables as .csv.gz. Default true.
#       --progress BOOL       Emit stage and progress logging. Default true.
#       --progress-updates N  Approximate progress updates per phase.
#                             Default 20; must be a positive integer.
#       --seed N              Base RNG seed. Default 1. The ecDNA model uses
#                             seed + 5 and the ecDNA simulation seed + 6.
#   -h, --help                Print the usage block and exit with status 0.
#
# Inputs read: the --lineage CSV, the --params JSON, the optional --founders
# and --live-cells CSVs, and any barcode-reference file the JSON names relative
# to the parameter file's own directory.
#
# Outputs, all under the resolved output directory. Every .csv below is written
# as .csv.gz unless --compress-csv false is passed.
#   always            lineage_nodes.csv, terminal_cells.csv,
#                     physicell_lineage_full.nwk,
#                     physicell_lineage_sampled.nwk, r_timing_summary.csv
#   barcode           mutation_events.csv, barcode_target_layout.csv,
#                     barcode_profiles.rds, run_manifest.csv, plus the dense or
#                     sparse allele/score/character matrices belonging to the
#                     configured recorder and any optional FASTA files
#   mitochondrial     mitochondrial_mutation_events.csv,
#                     mitochondrial_variant_fractions.csv,
#                     mitochondrial_profiles.rds,
#                     mitochondrial_reference.fasta,
#                     mitochondrial_manifest.csv, the optional sparse
#                     variant-fraction matrix, and
#                     mitochondrial_sampled_haplotypes.fasta when
#                     --write-mt-fasta true
#   ecDNA             ecdna_manifest.csv, ecdna_mutation_events.csv,
#                     ecdna_cell_summary.csv, ecdna_haplotypes.csv,
#                     ecdna_species_manifest.csv,
#                     ecdna_terminal_profiles.rds, and the ecDNA sparse
#                     static-ID / recorder matrices
#   barcode or mt     mutation_event_descendant_matrix_sparse.rds and
#                     mutation_event_descendant_manifest.csv
#   barcode and mt    combined_lineage_feature_matrix.rds and
#                     combined_lineage_feature_manifest.csv, with the ecDNA
#                     feature blocks appended when ecDNA also ran
#
# The last stdout line is OUTPUT_DIR=<path>, which the bash wrappers parse.

r_simulator_start_time <- unname(proc.time()[['elapsed']])

# ---- Bootstrap: locate the repository root and load simulator modules ----

script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- normalizePath(
  sub('^--file=', '', script_argument),
  mustWork = TRUE
)
# The repository root is 3 levels above this script; every path this
# script resolves, including its outputs, is relative to that root.
repo_root <- normalizePath(file.path(dirname(script_path), '..', '..', '..'))
source(file.path(repo_root, 'load_origin.R'))

# ---- Command-line interface ----

#' Build this script's command-line usage text
#'
#' @return One string containing the newline-separated usage block: the
#'   invocation line, the two required arguments, and every supported option
#'   with its default.
physicell_cli_usage <- function(){
  paste(
    'Usage:',
    '  Rscript simulate_physicell_lineage.R --lineage DIVISIONS.csv --params PARAMS.json [options]',
    '',
    'Required:',
    '  -L, --lineage PATH          PhysiCell time,parent_ID,daughter_ID CSV',
    '  -P, --params PATH           remote_mito JSON parameter file',
    '',
    'Options:',
    '  -O, --output-dir PATH       Default: output/physicell/<lineage stem>',
    '      --end-time NUMBER       Default: maximum sim_length in the JSON',
    '      --founder-time NUMBER   Default: 0',
    '      --founders PATH         Optional CSV containing every day-zero founder ID',
    '      --live-cells PATH       Optional PhysiCell live-cell CSV with an ID column',
    '      --cell-type NAME        Default: founder_cell_type',
    '      --editing-state STATE   auto, induced, or uninduced; default: auto',
    paste(
      '      --modalities LIST       lineage, barcode, mitochondrial, ecDNA,',
      'both, or all; default: barcode'
    ),
    '      --num-integrations N    Default: maximum max_bc_ints_per_cell',
    '      --founder-label-sites N Stable allele-coded founder barcode sites; default: 0',
    '      --mt-genomes-per-cell N Fixed mitochondrial bottleneck size; default: 8',
    '      --write-mt-fasta BOOL   Write one sampled mt haplotype/cell; default: false',
    '      --compress-csv BOOL     Write CSV tables as .csv.gz; default: true',
    '      --progress BOOL         Show replay progress; default: true',
    '      --progress-updates N    Approximate updates per phase; default: 20',
    '      --seed N                Default: 1',
    '  -h, --help',
    sep = '\n'
  )
}

#' Parse and validate the trailing command-line arguments
#'
#' Walks the arguments left to right, accepting both `--flag value` and
#' `--flag=value` (everything after the first `=` is the value). Flags are
#' resolved through a short/long alias table; an unrecognised flag, a flag with
#' no value, or a value that fails its numeric, integer, or true/false
#' coercion raises an error, and the unknown-option error embeds the usage
#' text. Later occurrences of the same flag overwrite earlier ones.
#'
#' @param arguments Character vector of trailing command-line arguments, as
#'   returned by `commandArgs(trailingOnly = TRUE)`.
#' @return A named list of options with defaults filled in: `lineage`,
#'   `params`, and `output_dir` (`NULL` until supplied), `end_time`
#'   (`NA_real_`), `founder_time` (`0`), `founders` and `live_cells` (`NULL`),
#'   `cell_type` (`NULL`), `editing_state` (`'auto'`), `modalities`
#'   (`'barcode'`), `num_integrations` (`NA_integer_`), `founder_label_sites`
#'   (`0`), `mt_genomes_per_cell` (`8`), `write_mt_fasta` (`FALSE`),
#'   `compress_csv` (`TRUE`), `progress` (`TRUE`), `progress_updates` (`20`),
#'   and `seed` (`1`). Path options are returned as given, not normalized.
#' @note Only `progress_updates >= 1` is checked here. That `--lineage` and
#'   `--params` were supplied, that the paths exist, and that `--modalities`
#'   names known modalities are checked later by the top-level script;
#'   `--editing-state` is validated further downstream, by the per-branch
#'   mutators in `physicell_lineage.R`.
#' @section Side effects: `-h`/`--help` prints the usage text and calls
#'   `quit(status = 0)`, ending the session rather than returning.
parse_physicell_cli_args <- function(arguments){
  options <- list(
    lineage = NULL,
    params = NULL,
    output_dir = NULL,
    end_time = NA_real_,
    founder_time = 0,
    founders = NULL,
    live_cells = NULL,
    cell_type = NULL,
    editing_state = 'auto',
    modalities = 'barcode',
    num_integrations = NA_integer_,
    founder_label_sites = 0L,
    mt_genomes_per_cell = 8L,
    write_mt_fasta = FALSE,
    compress_csv = TRUE,
    progress = TRUE,
    progress_updates = 20L,
    seed = 1L
  )
  aliases <- c(
    '-L' = 'lineage',
    '--lineage' = 'lineage',
    '-P' = 'params',
    '--params' = 'params',
    '-O' = 'output_dir',
    '--output-dir' = 'output_dir',
    '--end-time' = 'end_time',
    '--founder-time' = 'founder_time',
    '--founders' = 'founders',
    '--live-cells' = 'live_cells',
    '--cell-type' = 'cell_type',
    '--editing-state' = 'editing_state',
    '--modalities' = 'modalities',
    '--num-integrations' = 'num_integrations',
    '--founder-label-sites' = 'founder_label_sites',
    '--mt-genomes-per-cell' = 'mt_genomes_per_cell',
    '--write-mt-fasta' = 'write_mt_fasta',
    '--compress-csv' = 'compress_csv',
    '--progress' = 'progress',
    '--progress-updates' = 'progress_updates',
    '--seed' = 'seed'
  )
  numeric_options <- c('end_time', 'founder_time')
  integer_options <- c(
    'num_integrations', 'founder_label_sites', 'mt_genomes_per_cell',
    'progress_updates', 'seed'
  )

  argument_index <- 1L
  while(argument_index <= length(arguments)){
    argument <- arguments[argument_index]
    if(argument %in% c('-h', '--help')){
      cat(physicell_cli_usage(), '\n')
      quit(save = 'no', status = 0)
    }

    inline_parts <- strsplit(argument, '=', fixed = TRUE)[[1]]
    flag <- inline_parts[1]
    option_name <- unname(aliases[flag])
    if(length(option_name) == 0 || is.na(option_name)){
      stop(sprintf('Unknown option: %s\n\n%s', flag, physicell_cli_usage()))
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

    if(option_name %in% numeric_options){
      value <- suppressWarnings(as.numeric(value))
      if(is.na(value)){
        stop(sprintf('Option %s requires a numeric value.', flag))
      }
    } else if(option_name %in% integer_options){
      value <- suppressWarnings(as.integer(value))
      if(is.na(value)){
        stop(sprintf('Option %s requires an integer value.', flag))
      }
    } else if(option_name %in% c(
      'write_mt_fasta',
      'compress_csv',
      'progress'
    )){
      normalized_value <- tolower(value)
      if(!(normalized_value %in% c('true', 'false'))){
        stop(sprintf('Option %s requires true or false.', flag))
      }
      value <- normalized_value == 'true'
    }
    options[[option_name]] <- value
    argument_index <- argument_index + 1L
  }
  if(options$progress_updates < 1){
    stop('--progress-updates must be one positive integer.')
  }
  options
}

# ---- Phase timing ----

#' Create a wall-clock recorder for the simulator's sequential phases
#'
#' Phases are contiguous rather than independent: each `finish_phase()` call
#' closes the interval that began when the recorder was created or when the
#' previous phase finished, so the recorded phases partition the whole run and
#' their durations sum to the total. All durations are elapsed (wall) seconds,
#' not CPU seconds.
#'
#' @param total_start Elapsed-seconds timestamp marking the start of the run,
#'   on the same scale as `clock`. Must be one finite number; defaults to
#'   `clock()` evaluated at construction time.
#' @param clock Zero-argument function returning the current elapsed seconds;
#'   defaults to `proc.time()[['elapsed']]`. It must return one finite value
#'   that never decreases between calls.
#' @return A named list of two closures, `finish_phase` and `summary`, sharing
#'   the recorder's accumulated phase table.
new_physicell_timing_recorder <- function(
    total_start = NULL,
    clock = function(){
      unname(proc.time()[['elapsed']])
    }){
  if(is.null(total_start)){
    total_start <- clock()
  }
  total_start <- as.numeric(total_start)
  if(length(total_start) != 1 || !is.finite(total_start)){
    stop('Timing total_start must be one finite number.')
  }
  phase_start <- total_start
  phase_timings <- setNames(numeric(), character())

  #' Internal: close the running phase and record its elapsed time
  #'
  #' @param label One non-empty string naming the phase; must not repeat a
  #'   label already recorded.
  #' @return Invisibly, the phase duration in seconds.
  #' @section Side effects: Appends the labelled duration to the enclosing
  #'   recorder's `phase_timings` and advances its `phase_start`.
  finish_phase <- function(label){
    label <- as.character(label)
    if(length(label) != 1 || is.na(label) || !nzchar(label)){
      stop('Timing phase label must be one non-empty string.')
    }
    if(label %in% names(phase_timings)){
      stop(sprintf('Timing phase label is duplicated: %s.', label))
    }
    now <- as.numeric(clock())
    if(length(now) != 1 || !is.finite(now) || now < phase_start){
      stop('Timing clock must return a finite, non-decreasing number.')
    }
    duration <- now - phase_start
    phase_timings <<- c(
      phase_timings,
      setNames(duration, label)
    )
    phase_start <<- now
    invisible(duration)
  }

  #' Internal: summarise every recorded phase plus the total runtime
  #'
  #' Reads the clock without closing a phase, so any work done since the last
  #' `finish_phase()` call is counted only in the total row.
  #'
  #' @return A data frame with one row per recorded phase followed by a
  #'   `Total R simulator` row, and columns `phase`, `elapsed_seconds`
  #'   (rounded to six decimals), `elapsed` (a compact label from
  #'   `format_physicell_progress_duration()`), and `percent_of_total`, which
  #'   is `0` for every row when the total duration is not positive.
  summary <- function(){
    now <- as.numeric(clock())
    if(length(now) != 1 || !is.finite(now) || now < phase_start){
      stop('Timing clock must return a finite, non-decreasing number.')
    }
    durations <- c(
      phase_timings,
      'Total R simulator' = now - total_start
    )
    total_duration <- unname(durations[['Total R simulator']])
    data.frame(
      phase = names(durations),
      elapsed_seconds = round(unname(durations), 6),
      elapsed = vapply(
        unname(durations),
        format_physicell_progress_duration,
        character(1)
      ),
      percent_of_total = if(total_duration > 0){
        100 * unname(durations) / total_duration
      } else{
        rep(0, length(durations))
      },
      stringsAsFactors = FALSE
    )
  }

  list(
    finish_phase = finish_phase,
    summary = summary
  )
}

#' Format a timing summary as aligned terminal and log lines
#'
#' @param timing_summary Data frame carrying at least the `phase` and
#'   `elapsed` columns produced by a timing recorder's `summary()`, with one
#'   or more rows.
#' @return A character vector: an `R simulator timing summary:` heading
#'   followed by one indented line per row, whose `phase:` labels are padded
#'   to a common width so the elapsed values line up.
format_physicell_timing_summary <- function(timing_summary){
  required <- c('phase', 'elapsed')
  if(!is.data.frame(timing_summary) ||
     !all(required %in% names(timing_summary)) ||
     nrow(timing_summary) == 0){
    stop('timing_summary must contain phase and elapsed columns.')
  }
  labels <- paste0(as.character(timing_summary$phase), ':')
  label_width <- max(nchar(labels, type = 'width'))
  c(
    'R simulator timing summary:',
    paste0(
      '  ',
      format(labels, width = label_width, justify = 'left'),
      '  ',
      timing_summary$elapsed
    )
  )
}

# ---- Parse arguments and read the PhysiCell and parameter inputs ----

options <- parse_physicell_cli_args(commandArgs(trailingOnly = TRUE))
timing_recorder <- new_physicell_timing_recorder(
  total_start = r_simulator_start_time
)
if(is.null(options$lineage) || is.null(options$params)){
  stop(sprintf('--lineage and --params are required.\n\n%s', physicell_cli_usage()))
}
if(!requireNamespace('jsonlite', quietly = TRUE)){
  stop('The jsonlite package is required to read the parameter JSON.')
}

lineage_path <- normalizePath(options$lineage, mustWork = TRUE)
params_path <- normalizePath(options$params, mustWork = TRUE)
live_cells_path <- if(is.null(options$live_cells)){
  NULL
} else{
  normalizePath(options$live_cells, mustWork = TRUE)
}
founders_path <- if(is.null(options$founders)){
  NULL
} else{
  normalizePath(options$founders, mustWork = TRUE)
}

params <- jsonlite::fromJSON(params_path, simplifyVector = FALSE)
physicell_log_stage(
  'Reading and validating PhysiCell lineage inputs.',
  enabled = options$progress
)
division_events <- read_physicell_divisions(lineage_path)
terminal_ids <- if(is.null(live_cells_path)){
  NULL
} else{
  read_physicell_terminal_ids(live_cells_path)
}
supplied_founder_ids <- if(is.null(founders_path)){
  NULL
} else{
  read_physicell_founder_ids(founders_path)
}

# ---- Reconstruct the event-resolved lineage ----
# The founders are taken from --founders when supplied. When no divisions were
# recorded the live cells themselves become the founders; otherwise the
# founders are inferred from the division events.

end_time <- options$end_time
if(is.na(end_time)){
  end_time <- max(as.numeric(unlist(params$sim_length, use.names = FALSE)))
}
founder_ids <- if(!is.null(supplied_founder_ids)){
  supplied_founder_ids
} else if(nrow(division_events) == 0){
  terminal_ids
} else{
  NULL
}

nodes <- build_physicell_lineage(
  division_events,
  end_time = end_time,
  founder_time = options$founder_time,
  founder_ids = founder_ids,
  show_progress = options$progress,
  progress_updates = options$progress_updates
)
physicell_log_stage(
  sprintf(
    paste(
      'Lineage reconstruction finished: %s divisions became %s',
      'event-resolved nodes.'
    ),
    format(
      nrow(division_events),
      big.mark = ',',
      scientific = FALSE,
      trim = TRUE
    ),
    format(nrow(nodes), big.mark = ',', scientific = FALSE, trim = TRUE)
  ),
  enabled = options$progress
)

# ---- Resolve the output directory and the requested modalities ----
# `both` expands to barcode plus mitochondrial, `all` adds ecDNA, and `mito`
# is accepted as an alias for `mitochondrial`.

if(is.null(options$output_dir)){
  lineage_stem <- tools::file_path_sans_ext(basename(lineage_path))
  options$output_dir <- file.path(
    repo_root,
    'output',
    'physicell',
    lineage_stem
  )
}
dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)

requested_modalities <- tolower(trimws(unlist(strsplit(
  options$modalities,
  ',',
  fixed = TRUE
))))
if('both' %in% requested_modalities){
  requested_modalities <- c(
    requested_modalities[requested_modalities != 'both'],
    'barcode',
    'mitochondrial'
  )
}
if('all' %in% requested_modalities){
  requested_modalities <- c(
    requested_modalities[requested_modalities != 'all'],
    'barcode',
    'mitochondrial',
    'ecdna'
  )
}
requested_modalities[requested_modalities == 'mito'] <- 'mitochondrial'
requested_modalities[requested_modalities == 'ecdna'] <- 'ecdna'
requested_modalities <- unique(requested_modalities)
unknown_modalities <- setdiff(
  requested_modalities,
  c('lineage', 'barcode', 'mitochondrial', 'ecdna')
)
if(length(requested_modalities) == 0 || length(unknown_modalities) > 0){
  stop(paste(
    '--modalities must contain lineage, barcode, mitochondrial, ecDNA,',
    'both, or all.'
  ))
}

# ---- Select terminal cells and write the ground-truth lineage output ----
# Writes lineage_nodes.csv, terminal_cells.csv, physicell_lineage_full.nwk,
# and physicell_lineage_sampled.nwk. `--modalities lineage` alone stops the
# run after this section apart from the timing summary.

terminal_node_table <- nodes[nodes$is_terminal, , drop = FALSE]
if(!is.null(terminal_ids)){
  terminal_node_table <- terminal_node_table[
    terminal_node_table$physicell_id %in% terminal_ids,
    ,
    drop = FALSE
  ]
}
timing_recorder$finish_phase(
  'Startup, input, and lineage reconstruction'
)
write_physicell_lineage_outputs(
  nodes,
  terminal_node_table,
  options$output_dir,
  show_progress = options$progress,
  progress_updates = options$progress_updates,
  compress_csv = options$compress_csv
)
timing_recorder$finish_phase('Lineage tables and Newick output')

# ---- Barcode recording modality ----
# Prepares the recorder model from the JSON, replays it over every branch, and
# writes the barcode events, target layout, profiles, and score matrices.

barcode_simulation <- NULL
barcode_output <- NULL
if('barcode' %in% requested_modalities){
  physicell_log_stage(
    'Preparing barcode recording model.',
    enabled = options$progress
  )
  barcode_model <- prepare_physicell_recording_model(
    params,
    cell_type = options$cell_type,
    num_integrations = if(is.na(options$num_integrations)){
      NULL
    } else{
      options$num_integrations
    },
    founder_label_sites = options$founder_label_sites,
    params_dir = dirname(params_path),
    seed = options$seed
  )
  timing_recorder$finish_phase('Barcode model preparation')
  physicell_log_stage(
    sprintf(
      'Starting %s lineage simulation across %s nodes.',
      barcode_model$recorder_system,
      format(nrow(nodes), big.mark = ',', scientific = FALSE, trim = TRUE)
    ),
    enabled = options$progress
  )
  barcode_simulation <- simulate_recording_on_physicell_lineage(
    nodes,
    barcode_model,
    editing_state = tolower(options$editing_state),
    terminal_physicell_ids = terminal_ids,
    seed = options$seed,
    show_progress = options$progress,
    progress_updates = options$progress_updates
  )
  timing_recorder$finish_phase('Barcode lineage simulation')
  physicell_log_stage(
    sprintf(
      paste(
        '%s lineage simulation finished: %s events across %s terminal',
        'cells. Beginning barcode output.'
      ),
      barcode_model$recorder_system,
      format(
        nrow(barcode_simulation$mutation_events),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      ),
      format(
        nrow(barcode_simulation$terminal_nodes),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      )
    ),
    enabled = options$progress
  )
  barcode_output <- write_physicell_recording_outputs(
    barcode_simulation,
    barcode_model,
    options$output_dir,
    show_progress = options$progress,
    write_lineage = FALSE,
    progress_updates = options$progress_updates,
    compress_csv = options$compress_csv
  )
  timing_recorder$finish_phase('Barcode output')
}

# ---- Mitochondrial recording modality ----
# Propagates a fixed --mt-genomes-per-cell bottleneck along the same lineage
# and writes the mt events, variant fractions, profiles, and reference.

mitochondrial_simulation <- NULL
mitochondrial_output <- NULL
if('mitochondrial' %in% requested_modalities){
  physicell_log_stage(
    'Preparing mitochondrial lineage-tracing model.',
    enabled = options$progress
  )
  mitochondrial_model <- prepare_physicell_mito_model(
    params,
    cell_type = options$cell_type,
    genomes_per_cell = options$mt_genomes_per_cell,
    seed = options$seed
  )
  timing_recorder$finish_phase('Mitochondrial model preparation')
  physicell_log_stage(
    sprintf(
      'Starting mitochondrial lineage simulation across %s nodes.',
      format(nrow(nodes), big.mark = ',', scientific = FALSE, trim = TRUE)
    ),
    enabled = options$progress
  )
  mitochondrial_simulation <- simulate_mito_on_physicell_lineage(
    nodes,
    mitochondrial_model,
    params,
    editing_state = tolower(options$editing_state),
    terminal_physicell_ids = terminal_ids,
    seed = options$seed,
    show_progress = options$progress,
    progress_updates = options$progress_updates
  )
  timing_recorder$finish_phase('Mitochondrial lineage simulation')
  physicell_log_stage(
    sprintf(
      paste(
        'Mitochondrial lineage simulation finished: %s events across %s',
        'terminal cells. Beginning mitochondrial output.'
      ),
      format(
        nrow(mitochondrial_simulation$mutation_events),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      ),
      format(
        nrow(mitochondrial_simulation$terminal_nodes),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      )
    ),
    enabled = options$progress
  )
  mitochondrial_output <- write_physicell_mito_outputs(
    mitochondrial_simulation,
    mitochondrial_model,
    options$output_dir,
    write_fasta = options$write_mt_fasta,
    show_progress = options$progress,
    compress_csv = options$compress_csv
  )
  timing_recorder$finish_phase('Mitochondrial output')
}

# ---- ecDNA recording modality ----
# Non-Mendelian ecDNA propagation plus its CRISPR recorder. Uses offset seeds
# (--seed + 5 for the model, --seed + 6 for the simulation) and samples the
# terminal cells already chosen above.

ecdna_simulation <- NULL
ecdna_output <- NULL
if('ecdna' %in% requested_modalities){
  physicell_log_stage(
    'Preparing non-Mendelian ecDNA barcode model.',
    enabled = options$progress
  )
  ecdna_model <- prepare_ecdna_model(params, seed = options$seed + 5L)
  timing_recorder$finish_phase('ecDNA model preparation')
  physicell_log_stage(
    sprintf(
      'Starting ecDNA propagation and recorder simulation across %s nodes.',
      format(nrow(nodes), big.mark = ',', scientific = FALSE, trim = TRUE)
    ),
    enabled = options$progress
  )
  ecdna_simulation <- simulate_ecdna_on_lineage(
    nodes,
    ecdna_model,
    terminal_physicell_ids = terminal_node_table$physicell_id,
    seed = options$seed + 6L,
    show_progress = options$progress,
    progress_updates = options$progress_updates
  )
  timing_recorder$finish_phase('ecDNA lineage simulation')
  physicell_log_stage(
    sprintf(
      paste(
        'ecDNA lineage simulation finished with %s aggregated edit rows;',
        'beginning output.'
      ),
      format(
        nrow(ecdna_simulation$mutation_events),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      )
    ),
    enabled = options$progress
  )
  ecdna_output <- write_ecdna_outputs(
    ecdna_simulation,
    ecdna_model,
    options$output_dir,
    show_progress = options$progress,
    compress_csv = options$compress_csv
  )
  timing_recorder$finish_phase('ecDNA output')
}

# ---- Mutation-event descendant matrix ----
# Written whenever at least one of the barcode or mitochondrial simulations
# ran; ecDNA events are not included here.

event_tables <- list()
if(!is.null(barcode_simulation)){
  event_tables$barcode <- barcode_simulation$mutation_events
}
if(!is.null(mitochondrial_simulation)){
  event_tables$mitochondrial <- mitochondrial_simulation$mutation_events
}
if(length(event_tables) > 0L){
  write_physicell_event_descendant_outputs(
    nodes,
    terminal_node_table,
    event_tables,
    options$output_dir,
    show_progress = options$progress,
    progress_updates = options$progress_updates,
    compress_csv = options$compress_csv
  )
  timing_recorder$finish_phase('Mutation-event descendant matrix output')
}

# ---- Combined feature matrix ----
# Only when both barcode and mitochondrial data exist. The mitochondrial (and,
# if present, ecDNA) rows are reordered to the barcode row order before the
# column blocks are bound together; a terminal sample ID missing from any
# block is a fatal mismatch.

if(all(c('barcode', 'mitochondrial') %in% requested_modalities) &&
   !is.null(barcode_output) && !is.null(mitochondrial_output)){
  physicell_log_stage(
    'Combined output: aligning barcode and mitochondrial feature matrices.',
    enabled = options$progress
  )
  if(!requireNamespace('Matrix', quietly = TRUE)){
    stop('The Matrix package is required for combined lineage features.')
  }
  barcode_features <- barcode_output$binary_scores
  if(!inherits(barcode_features, 'Matrix')){
    barcode_features <- Matrix::Matrix(
      barcode_features,
      sparse = TRUE
    )
  }
  mitochondrial_features <- mitochondrial_output$sparse_scores
  if(is.null(mitochondrial_features)){
    stop('Mitochondrial sparse features were not generated.')
  }
  if(!identical(rownames(barcode_features), rownames(mitochondrial_features))){
    mitochondrial_indices <- match(
      rownames(barcode_features),
      rownames(mitochondrial_features)
    )
    if(anyNA(mitochondrial_indices)){
      stop('Barcode and mitochondrial terminal sample IDs do not match.')
    }
    mitochondrial_features <- mitochondrial_features[
      mitochondrial_indices,
      ,
      drop = FALSE
    ]
  }
  physicell_log_stage(
    'Combined output: writing aligned sparse feature matrix.',
    enabled = options$progress
  )
  combined_features <- cbind(barcode_features, mitochondrial_features)
  combined_modalities <- c(
    rep('barcode', ncol(barcode_features)),
    rep('mitochondrial', ncol(mitochondrial_features))
  )
  combined_systems <- c(
    rep(barcode_model$recorder_system, ncol(barcode_features)),
    rep('mitochondrial lineage tracing', ncol(mitochondrial_features))
  )
  if(!is.null(ecdna_output)){
    ecdna_recorder_features <- ecdna_output$recorder_edit_presence
    ecdna_static_features <- ecdna_output$static_id_copy_numbers
    ecdna_indices <- match(
      rownames(combined_features),
      rownames(ecdna_recorder_features)
    )
    if(anyNA(ecdna_indices)){
      stop('ecDNA and barcode/mitochondrial terminal sample IDs do not match.')
    }
    ecdna_recorder_features <- ecdna_recorder_features[
      ecdna_indices, , drop = FALSE
    ]
    ecdna_static_features <- ecdna_static_features[
      ecdna_indices, , drop = FALSE
    ]
    combined_features <- cbind(
      combined_features,
      ecdna_recorder_features,
      ecdna_static_features
    )
    combined_modalities <- c(
      combined_modalities,
      rep('ecdna_recorder', ncol(ecdna_recorder_features)),
      rep('ecdna_static_id_copy_number', ncol(ecdna_static_features))
    )
    combined_systems <- c(
      combined_systems,
      rep('ecDNA CRISPR recorder', ncol(ecdna_recorder_features)),
      rep('ecDNA static ID', ncol(ecdna_static_features))
    )
  }
  saveRDS(
    combined_features,
    file.path(options$output_dir, 'combined_lineage_feature_matrix.rds')
  )
  physicell_log_stage(
    'Combined output: writing feature manifest.',
    enabled = options$progress
  )
  feature_manifest <- data.frame(
    feature = colnames(combined_features),
    modality = combined_modalities,
    recorder_system = combined_systems,
    stringsAsFactors = FALSE
  )
  write_physicell_csv(
    feature_manifest,
    file.path(options$output_dir, 'combined_lineage_feature_manifest.csv'),
    row.names = FALSE,
    compress = options$compress_csv
  )
  physicell_log_stage(
    'Combined lineage output complete.',
    enabled = options$progress
  )
  timing_recorder$finish_phase('Combined lineage output')
}

# ---- Timing summary and stdout report ----
# Writes r_timing_summary.csv, prints the per-phase timing table, and ends with
# the OUTPUT_DIR=<path> line the wrapper scripts parse.

written_output_dir <- normalizePath(options$output_dir, mustWork = TRUE)
physicell_log_stage(
  'All requested lineage simulation and output phases completed.',
  enabled = options$progress
)

r_timing_summary <- timing_recorder$summary()
timing_summary_path <- write_physicell_csv(
  r_timing_summary,
  file.path(
    written_output_dir,
    'r_timing_summary.csv'
  ),
  row.names = FALSE,
  compress = options$compress_csv
)

cat(sprintf('Imported %d PhysiCell divisions into %d event-resolved nodes.\n',
            nrow(division_events), nrow(nodes)))
if(identical(requested_modalities, 'lineage')){
  cat(sprintf(
    'Generated ground-truth lineage data for %d terminal cells.\n',
    nrow(terminal_node_table)
  ))
} else{
  cat(sprintf(
    'Generated %s data for %d terminal cells.\n',
    paste(requested_modalities, collapse = ' + '),
    nrow(terminal_node_table)
  ))
}
cat(
  paste(format_physicell_timing_summary(r_timing_summary), collapse = '\n'),
  '\n'
)
cat(sprintf('R timing summary written to %s.\n', timing_summary_path))
cat(sprintf('OUTPUT_DIR=%s\n', written_output_dir))
