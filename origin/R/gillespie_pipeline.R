# Gillespie engine driver: turns one parsed remote_mito parameter set into a
# complete run directory of lineage, recording, and feature tables.
#
# Source prime_editing.R, physicell_lineage.R, physicell_mito.R,
# ecdna_lineage.R, and gillespie_lineage.R before this file; every stage below
# is an entry point defined in one of them. simulate_gillespie_lineage.R is the
# standalone CLI wrapper, and sim5_code.R sources this file when the JSON sets
# "simulation_engine": "gillespie".
#
# End-to-end stage order inside run_gillespie_lineage_pipeline():
#   1. Population. Consumes `params`; produces an exact continuous-time
#      birth/death forest: event-resolved node table, live terminal rows, and
#      population event/division/checkpoint tables.
#   2. Lineage output. Consumes that forest; produces node/tip CSVs, full and
#      sampled Newick trees, the three population event tables, and a
#      PhysiCell-shaped terminal-cell covariate table (gillespie_cell_states).
#   3. Barcode overlay (modality `barcode`). Consumes the node table plus
#      per-cell-type recording models; produces exact-time barcode mutation
#      events, the target layout, and a cell-by-position score matrix.
#   4. Mitochondrial overlay (modality `mitochondrial`). Consumes the node
#      table plus per-cell-type mt models; produces mt variant events,
#      heteroplasmy fractions, and a sparse cell-by-position score matrix.
#   5. ecDNA overlay (modality `ecdna`). Consumes the node table plus one ecDNA
#      model; produces copy-number and recorder-edit tables plus sparse
#      static-ID copy-number and recorder-edit matrices.
#   6. Event-descendant matrix. Consumes the barcode and mitochondrial event
#      tables; produces one sparse cell-by-unique-event matrix and a manifest.
#   7. Combined features. Consumes the barcode and mitochondrial score matrices
#      (ecDNA columns appended when that overlay also ran); produces
#      combined_lineage_feature_matrix.rds and its modality manifest.
#   8. Manifest and timing. Consume the run counters and phase timers; produce
#      gillespie_run_manifest.csv and r_timing_summary.csv.

#' Read one Gillespie configuration key with a fallback
#'
#' @export
#' @param configuration The `gillespie` block of a parsed parameter JSON, as a
#'   list; an empty list is acceptable.
#' @param name Character key to look up in that block.
#' @param default Value returned when the key is absent or `NULL`.
#' @return The configured value, or `default` when the key is missing.
gillespie_pipeline_option <- function(configuration, name, default){
  value <- configuration[[name]]
  if(is.null(value)) default else value
}

#' Normalise and validate the requested recording modalities
#'
#' Accepts a scalar, vector, or list; every element is lower-cased, trimmed,
#' and split on commas. The aliases `both` (barcode plus mitochondrial), `all`
#' (barcode, mitochondrial, and ecDNA), and `mito` (mitochondrial) are expanded
#' before validation, so mixed-case spellings such as `ecDNA` are accepted.
#'
#' @export
#' @param value Modality selection taken from the `gillespie` JSON block or a
#'   command-line override.
#' @return A character vector of unique modality names drawn from `barcode`,
#'   `mitochondrial`, and `ecdna`, in first-seen order.
#' @note Stops when the normalised result is empty or names anything outside
#'   those three modalities.
gillespie_modalities <- function(value){
  modalities <- tolower(trimws(as.character(unlist(value, use.names = FALSE))))
  modalities <- unlist(strsplit(modalities, ',', fixed = TRUE), use.names = FALSE)
  modalities <- trimws(modalities)
  if('both' %in% modalities){
    modalities <- c(modalities[modalities != 'both'], 'barcode', 'mitochondrial')
  }
  if('all' %in% modalities){
    modalities <- c(
      modalities[modalities != 'all'],
      'barcode',
      'mitochondrial',
      'ecdna'
    )
  }
  modalities[modalities == 'mito'] <- 'mitochondrial'
  modalities <- unique(modalities[nzchar(modalities)])
  unknown <- setdiff(modalities, c('barcode', 'mitochondrial', 'ecdna'))
  if(length(modalities) == 0L || length(unknown) > 0L){
    stop('Gillespie modalities must contain barcode, mitochondrial, ecDNA, both, or all.')
  }
  modalities
}

#' Create a wall-clock phase timer for one pipeline run
#'
#' Records elapsed `proc.time()` seconds at construction. Each phase is
#' measured as the gap since the previous `finish()` call, so the recorded
#' phases tile the run without gaps or overlap.
#'
#' @export
#' @return A named list of two closures, `finish` and `summary`, sharing one
#'   accumulating named numeric vector of phase durations.
gillespie_timing_recorder <- function(){
  total_start <- unname(proc.time()[['elapsed']])
  phase_start <- total_start
  phases <- setNames(numeric(), character())
  # Internal: Close the current phase and open the next one
  #
  # @param name Character label stored for the interval that just ended.
  # @return `NULL`, invisibly.
  # @section Side effects: Appends to the enclosing `phases` vector and moves
  #   the enclosing `phase_start` marker to the current time.
  finish <- function(name){
    now <- unname(proc.time()[['elapsed']])
    phases <<- c(phases, setNames(now - phase_start, name))
    phase_start <<- now
    invisible(NULL)
  }
  # Internal: Summarise every phase recorded so far
  #
  # @return A data frame with one row per finished phase plus a trailing
  #   `Total R simulator` row, holding `phase`, `elapsed_seconds` (rounded to
  #   six digits), and `percent_of_total` (all zero when the total is not
  #   positive).
  summary <- function(){
    now <- unname(proc.time()[['elapsed']])
    values <- c(phases, 'Total R simulator' = now - total_start)
    data.frame(
      phase = names(values),
      elapsed_seconds = round(unname(values), 6),
      percent_of_total = if(values[['Total R simulator']] > 0){
        100 * unname(values) / values[['Total R simulator']]
      } else{
        rep(0, length(values))
      },
      stringsAsFactors = FALSE
    )
  }
  list(finish = finish, summary = summary)
}

#' Run the complete Gillespie lineage-recording pipeline
#'
#' Simulates one exact continuous-time population, replays the requested
#' recording modalities on its branches, and writes every resulting table,
#' tree, and matrix into `output_dir`. Each option is resolved in the order
#' `overrides`, then the `gillespie` block of `params`, then a built-in
#' default.
#'
#' @export
#' @details Recognised options and their defaults: `progress` (`TRUE`),
#'   `progress_updates` (`20`), `compress_csv` (`TRUE`), `seed`
#'   (`params$random_seed`), `end_time` (`max(params$sim_length)`), `max_cells`
#'   (`1000000`), `modalities` (`both`), `num_integrations` (`NULL`, meaning
#'   the model's own default), `founder_label_sites` (`0`),
#'   `mt_genomes_per_cell` (`8`), and `write_mt_fasta` (`FALSE`).
#'
#'   Derived seeds are fixed offsets of the resolved base seed, so each overlay
#'   is reproducible independently: `seed + 1` and `seed + 2` for the barcode
#'   model and its branch simulation, `seed + 3` and `seed + 4` for the
#'   mitochondrial model and simulation, `seed + 5` and `seed + 6` for the
#'   ecDNA model and simulation.
#'
#'   The eight stages are those listed in the file header. Stages 3 to 5 run
#'   only for the modalities present in the normalised `modalities` vector;
#'   stage 6 runs whenever at least one of the barcode and mitochondrial
#'   overlays produced events; stage 7 runs only when both of those overlays
#'   ran. Every stage boundary is closed on the phase timer, and the phase
#'   names appear verbatim in `r_timing_summary.csv`.
#'
#' @param params Parsed remote_mito parameter list.
#'   `params$cell_type_dict$founder_cell_type` selects which per-cell-type
#'   model is treated as the reference model for output,
#'   `params$cell_type_dict$cell_type_params` supplies the cell-type names
#'   recorded in the manifest, and `params$sim_length` and `params$random_seed`
#'   supply the default end time and seed.
#' @param params_path Path the parameters were read from. Only its directory is
#'   used, as the root for resolving relative barcode-reference paths; `NULL`
#'   means the working directory.
#' @param output_dir Destination directory for every output file; created
#'   recursively when it does not exist.
#' @param overrides Named list of Gillespie options, typically from the
#'   command line, taking precedence over the `gillespie` JSON block.
#' @return Invisibly, a named list with `population` (the full population
#'   result), `barcode_simulation`, `mitochondrial_simulation`, and
#'   `ecdna_simulation` (each `NULL` when that modality was not requested),
#'   `output_dir` (the normalised destination path), and `timing_summary` (the
#'   phase timing data frame).
#' @section Side effects: Creates `output_dir`. Writes the lineage node/tip
#'   tables and Newick trees, `gillespie_population_events.csv`,
#'   `gillespie_division_events.csv`, `gillespie_checkpoint_summary.csv`,
#'   `gillespie_cell_states.csv`, the per-modality recording outputs,
#'   the event-descendant matrix and its manifest,
#'   `combined_lineage_feature_matrix.rds`,
#'   `combined_lineage_feature_manifest.csv`, `gillespie_run_manifest.csv`, and
#'   `r_timing_summary.csv` (CSV tables gzip-compressed unless `compress_csv`
#'   is false). Seeds the RNG through the simulator entry points it calls.
#'   Prints the timing table and the trailing `GILLESPIE_OUTPUT_DIR=` and
#'   `OUTPUT_DIR=` lines that the wrapper shell scripts parse.
#' @note Stops when `output_dir` cannot be created, when no cells survive the
#'   population simulation, when the `Matrix` package is missing while
#'   combining features, or when the ecDNA terminal sample IDs do not match the
#'   barcode and mitochondrial ones. Warns, without changing behaviour, when
#'   `params$fusion_events_per_mito_per_division` or
#'   `params$split_events_per_mito_per_division` is positive: this engine
#'   applies a fixed-genome bottleneck at division instead of replaying
#'   organelle fusion and split counts. The ecDNA feature columns reach the
#'   combined matrix only when the barcode and mitochondrial overlays both ran,
#'   and ecDNA events are never included in the event-descendant matrix.
run_gillespie_lineage_pipeline <- function(params,
                                            params_path = NULL,
                                            output_dir,
                                            overrides = list()){
  configuration <- params$gillespie
  if(is.null(configuration)){
    configuration <- list()
  }
  # Internal: Resolve one option from the overrides then the configuration
  #
  # @param name Character option name.
  # @param default Value used when neither `overrides` nor `configuration`
  #   supplies the option.
  # @return The override when one is present and not `NULL`, otherwise the
  #   configured value, otherwise `default`.
  value <- function(name, default){
    if(!is.null(overrides[[name]])){
      overrides[[name]]
    } else{
      gillespie_pipeline_option(configuration, name, default)
    }
  }
  progress <- isTRUE(value('progress', TRUE))
  progress_updates <- as.integer(value('progress_updates', 20L))
  compress_csv <- isTRUE(value('compress_csv', TRUE))
  seed <- as.integer(value('seed', params$random_seed))
  end_time <- as.numeric(value(
    'end_time',
    max(as.numeric(unlist(params$sim_length, use.names = FALSE)))
  ))
  max_cells <- as.integer(value('max_cells', 1000000L))
  modalities <- gillespie_modalities(value('modalities', 'both'))
  num_integrations <- value('num_integrations', NULL)
  if(!is.null(num_integrations)){
    num_integrations <- as.integer(num_integrations)
  }
  founder_label_sites <- as.integer(value('founder_label_sites', 0L))
  mt_genomes_per_cell <- as.integer(value('mt_genomes_per_cell', 8L))
  write_mt_fasta <- isTRUE(value('write_mt_fasta', FALSE))
  # Off by default: on a large run these four files cost more than the
  # simulation that produced them, and nothing downstream of tree building
  # opens any of them. See write_physicell_recording_outputs().
  write_allele_matrix <- isTRUE(value('write_allele_matrix', FALSE))
  write_mutation_events <- isTRUE(value('write_mutation_events', FALSE))
  write_profiles <- isTRUE(value('write_profiles', FALSE))
  write_barcode_fasta <- isTRUE(value('write_barcode_fasta', FALSE))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(output_dir)){
    stop(sprintf('Could not create Gillespie output directory: %s.', output_dir))
  }
  params_dir <- if(is.null(params_path)) '.' else dirname(params_path)
  timings <- gillespie_timing_recorder()

  if(progress){
    cat(sprintf(
      '[Gillespie] Starting exact population simulation through time %.6g.\n',
      end_time
    ))
  }
  population <- simulate_gillespie_population(
    params,
    end_time = end_time,
    seed = seed,
    max_cells = max_cells,
    show_progress = progress,
    progress_updates = progress_updates
  )
  timings$finish('Gillespie population simulation')
  if(nrow(population$terminal_nodes) == 0L){
    stop('No cells survived the Gillespie population simulation.')
  }
  if(progress){
    cat(sprintf(
      paste(
        '[Gillespie] Population simulation finished: %s divisions, %s',
        'deaths, %s surviving cells, and %s event-resolved nodes.\n'
      ),
      format(population$counts[['divisions']], big.mark = ','),
      format(population$counts[['deaths']], big.mark = ','),
      format(population$counts[['surviving_cells']], big.mark = ','),
      format(population$counts[['nodes']], big.mark = ',')
    ))
  }

  write_physicell_lineage_outputs(
    population$nodes,
    population$terminal_nodes,
    output_dir,
    show_progress = progress,
    progress_updates = progress_updates,
    compress_csv = compress_csv
  )
  write_physicell_csv(
    population$event_log,
    file.path(output_dir, 'gillespie_population_events.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    population$division_events,
    file.path(output_dir, 'gillespie_division_events.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    population$checkpoint_summary,
    file.path(output_dir, 'gillespie_checkpoint_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    gillespie_sc_cell_states(population),
    file.path(output_dir, 'gillespie_cell_states.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  timings$finish('Population tables and Newick output')

  barcode_simulation <- barcode_output <- barcode_model <- NULL
  if('barcode' %in% modalities){
    if(progress){
      cat('[Gillespie] Preparing cell-type-specific barcode hazards.\n')
    }
    barcode_models <- prepare_gillespie_barcode_models(
      params,
      num_integrations = num_integrations,
      founder_label_sites = founder_label_sites,
      params_dir = params_dir,
      seed = seed + 1L
    )
    founder_type <- as.character(params$cell_type_dict$founder_cell_type)
    barcode_model <- barcode_models[[founder_type]]
    timings$finish('Barcode model preparation')
    if(progress){
      cat(sprintf(
        '[Gillespie] Simulating continuous barcode recording across %s nodes.\n',
        format(nrow(population$nodes), big.mark = ',')
      ))
    }
    barcode_simulation <- simulate_recording_on_physicell_lineage(
      population$nodes,
      barcode_models,
      editing_state = 'auto',
      terminal_physicell_ids = population$terminal_nodes$physicell_id,
      seed = seed + 2L,
      show_progress = progress,
      progress_updates = progress_updates
    )
    timings$finish('Barcode lineage simulation')
    if(progress){
      cat(sprintf(
        paste(
          '[Gillespie] Barcode lineage simulation finished with %s events;',
          'beginning output.\n'
        ),
        format(nrow(barcode_simulation$mutation_events), big.mark = ',')
      ))
    }
    barcode_output <- write_physicell_recording_outputs(
      barcode_simulation,
      barcode_model,
      output_dir,
      show_progress = progress,
      write_lineage = FALSE,
      progress_updates = progress_updates,
      compress_csv = compress_csv,
      write_allele_matrix = write_allele_matrix,
      write_mutation_events = write_mutation_events,
      write_profiles = write_profiles,
      write_barcode_fasta = write_barcode_fasta
    )
    timings$finish('Barcode output')
  }

  mitochondrial_simulation <- mitochondrial_output <- mitochondrial_model <- NULL
  if('mitochondrial' %in% modalities){
    fusion_rate <- suppressWarnings(as.numeric(
      params$fusion_events_per_mito_per_division
    ))
    split_rate <- suppressWarnings(as.numeric(
      params$split_events_per_mito_per_division
    ))
    if(any(c(fusion_rate, split_rate) > 0, na.rm = TRUE)){
      warning(paste(
        'The Gillespie mitochondrial overlay uses a fixed-genome bottleneck',
        'at division; legacy organelle fusion/split counts are not replayed.'
      ), call. = FALSE)
    }
    if(progress){
      cat('[Gillespie] Preparing cell-type-specific mitochondrial hazards.\n')
    }
    mitochondrial_models <- prepare_gillespie_mito_models(
      params,
      genomes_per_cell = mt_genomes_per_cell,
      seed = seed + 3L
    )
    founder_type <- as.character(params$cell_type_dict$founder_cell_type)
    mitochondrial_model <- mitochondrial_models[[founder_type]]
    timings$finish('Mitochondrial model preparation')
    if(progress){
      cat(sprintf(
        '[Gillespie] Simulating continuous mitochondrial tracing across %s nodes.\n',
        format(nrow(population$nodes), big.mark = ',')
      ))
    }
    mitochondrial_simulation <- simulate_mito_on_physicell_lineage(
      population$nodes,
      mitochondrial_models,
      params,
      editing_state = 'auto',
      terminal_physicell_ids = population$terminal_nodes$physicell_id,
      seed = seed + 4L,
      show_progress = progress,
      progress_updates = progress_updates
    )
    timings$finish('Mitochondrial lineage simulation')
    if(progress){
      cat(sprintf(
        paste(
          '[Gillespie] Mitochondrial lineage simulation finished with %s',
          'events; beginning output.\n'
        ),
        format(nrow(mitochondrial_simulation$mutation_events), big.mark = ',')
      ))
    }
    mitochondrial_output <- write_physicell_mito_outputs(
      mitochondrial_simulation,
      mitochondrial_model,
      output_dir,
      write_fasta = write_mt_fasta,
      show_progress = progress,
      compress_csv = compress_csv
    )
    timings$finish('Mitochondrial output')
  }

  ecdna_simulation <- ecdna_output <- ecdna_model <- NULL
  if('ecdna' %in% modalities){
    if(progress){
      cat('[Gillespie] Preparing non-Mendelian ecDNA barcode model.\n')
    }
    ecdna_model <- prepare_ecdna_model(params, seed = seed + 5L)
    timings$finish('ecDNA model preparation')
    if(progress){
      cat(sprintf(
        '[Gillespie] Simulating ecDNA propagation across %s nodes.\n',
        format(nrow(population$nodes), big.mark = ',')
      ))
    }
    ecdna_simulation <- simulate_ecdna_on_lineage(
      population$nodes,
      ecdna_model,
      terminal_physicell_ids = population$terminal_nodes$physicell_id,
      seed = seed + 6L,
      show_progress = progress,
      progress_updates = progress_updates
    )
    timings$finish('ecDNA lineage simulation')
    if(progress){
      cat(sprintf(
        paste(
          '[Gillespie] ecDNA simulation finished with %s aggregated edit',
          'rows; beginning output.\n'
        ),
        format(nrow(ecdna_simulation$mutation_events), big.mark = ',')
      ))
    }
    ecdna_output <- write_ecdna_outputs(
      ecdna_simulation,
      ecdna_model,
      output_dir,
      show_progress = progress,
      compress_csv = compress_csv
    )
    timings$finish('ecDNA output')
  }

  event_tables <- list()
  if(!is.null(barcode_simulation)){
    event_tables$barcode <- barcode_simulation$mutation_events
  }
  if(!is.null(mitochondrial_simulation)){
    event_tables$mitochondrial <- mitochondrial_simulation$mutation_events
  }
  if(length(event_tables) > 0L){
    write_physicell_event_descendant_outputs(
      population$nodes,
      population$terminal_nodes,
      event_tables,
      output_dir,
      show_progress = progress,
      progress_updates = progress_updates,
      compress_csv = compress_csv
    )
    timings$finish('Mutation-event descendant matrix output')
  }

  if(!is.null(barcode_output) && !is.null(mitochondrial_output)){
    if(!requireNamespace('Matrix', quietly = TRUE)){
      stop('The Matrix package is required for combined lineage features.')
    }
    barcode_features <- barcode_output$binary_scores
    if(!inherits(barcode_features, 'Matrix')){
      barcode_features <- Matrix::Matrix(barcode_features, sparse = TRUE)
    }
    mitochondrial_features <- mitochondrial_output$sparse_scores
    mitochondrial_features <- mitochondrial_features[
      match(rownames(barcode_features), rownames(mitochondrial_features)),
      ,
      drop = FALSE
    ]
    combined_features <- cbind(barcode_features, mitochondrial_features)
    combined_modalities <- c(
      rep('barcode', ncol(barcode_features)),
      rep('mitochondrial', ncol(mitochondrial_features))
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
    }
    saveRDS(
      combined_features,
      file.path(output_dir, 'combined_lineage_feature_matrix.rds')
    )
    feature_manifest <- data.frame(
      feature = colnames(combined_features),
      modality = combined_modalities,
      stringsAsFactors = FALSE
    )
    write_physicell_csv(
      feature_manifest,
      file.path(output_dir, 'combined_lineage_feature_manifest.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
    timings$finish('Combined lineage output')
  }

  manifest <- data.frame(
    property = c(
      'simulation_engine', 'end_time', 'requested_end_time', 'stop_reason',
      'random_seed', 'founders', 'divisions', 'deaths', 'lineage_nodes',
      'surviving_cells', 'cell_types', 'modalities', 'max_cells'
    ),
    value = c(
      'gillespie', population$end_time, population$requested_end_time,
      population$stop_reason, seed, population$counts[['founders']],
      population$counts[['divisions']], population$counts[['deaths']],
      population$counts[['nodes']], population$counts[['surviving_cells']],
      paste(names(params$cell_type_dict$cell_type_params), collapse = ';'),
      paste(modalities, collapse = ';'), max_cells
    ),
    stringsAsFactors = FALSE
  )
  write_physicell_csv(
    manifest,
    file.path(output_dir, 'gillespie_run_manifest.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  timings$finish('Manifest output')
  timing_summary <- timings$summary()
  timing_path <- write_physicell_csv(
    timing_summary,
    file.path(output_dir, 'r_timing_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )

  written_output_dir <- normalizePath(output_dir, mustWork = TRUE)
  cat('R simulator timing summary:\n')
  for(row_index in seq_len(nrow(timing_summary))){
    cat(sprintf(
      '  %-45s %10.3f s\n',
      paste0(timing_summary$phase[row_index], ':'),
      timing_summary$elapsed_seconds[row_index]
    ))
  }
  cat(sprintf('R timing summary written to %s.\n', timing_path))
  cat(sprintf('GILLESPIE_OUTPUT_DIR=%s\n', written_output_dir))
  cat(sprintf('OUTPUT_DIR=%s\n', written_output_dir))

  invisible(list(
    population = population,
    barcode_simulation = barcode_simulation,
    mitochondrial_simulation = mitochondrial_simulation,
    ecdna_simulation = ecdna_simulation,
    output_dir = written_output_dir,
    timing_summary = timing_summary
  ))
}
