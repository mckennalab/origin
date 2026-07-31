#!/usr/bin/env Rscript

r_simulator_start_time <- unname(proc.time()[['elapsed']])

script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- normalizePath(
  sub('^--file=', '', script_argument),
  mustWork = TRUE
)
repo_root <- dirname(script_path)
source(file.path(repo_root, 'physicell_lineage.R'))
source(file.path(repo_root, 'physicell_mito.R'))

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
    '      --modalities LIST       barcode, mitochondrial, or both; default: barcode',
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
  requested_modalities <- c('barcode', 'mitochondrial')
}
requested_modalities[requested_modalities == 'mito'] <- 'mitochondrial'
requested_modalities <- unique(requested_modalities)
unknown_modalities <- setdiff(
  requested_modalities,
  c('barcode', 'mitochondrial')
)
if(length(requested_modalities) == 0 || length(unknown_modalities) > 0){
  stop('--modalities must contain barcode, mitochondrial, or both.')
}

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

event_tables <- list()
if(!is.null(barcode_simulation)){
  event_tables$barcode <- barcode_simulation$mutation_events
}
if(!is.null(mitochondrial_simulation)){
  event_tables$mitochondrial <- mitochondrial_simulation$mutation_events
}
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
    modality = c(
      rep('barcode', ncol(barcode_features)),
      rep('mitochondrial', ncol(mitochondrial_features))
    ),
    recorder_system = c(
      rep(barcode_model$recorder_system, ncol(barcode_features)),
      rep('mitochondrial lineage tracing', ncol(mitochondrial_features))
    ),
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
cat(sprintf(
  'Generated %s recording data for %d terminal cells.\n',
  paste(requested_modalities, collapse = ' + '),
  nrow(terminal_node_table)
))
cat(
  paste(format_physicell_timing_summary(r_timing_summary), collapse = '\n'),
  '\n'
)
cat(sprintf('R timing summary written to %s.\n', timing_summary_path))
cat(sprintf('OUTPUT_DIR=%s\n', written_output_dir))
