# Scalable mitochondrial lineage-recording replay for event-resolved PhysiCell
# trees. Source physicell_lineage.R before using these functions.

physicell_mito_rate_set <- function(params,
                                    cell_type,
                                    editing_state,
                                    mitochondrial_reference){
  cell_params <- params$cell_type_dict$cell_type_params[[cell_type]]
  if(is.null(cell_params)){
    stop(sprintf('Cell type %s is absent from the parameter JSON.', cell_type))
  }
  state_params <- cell_params[[editing_state]]
  if(is.null(state_params)){
    stop(sprintf(
      'Editing state %s is absent for cell type %s.',
      editing_state,
      cell_type
    ))
  }
  cell_cycle_length <- as.numeric(cell_params$cell_cycle_length)
  if(length(cell_cycle_length) != 1 || !is.finite(cell_cycle_length) ||
     cell_cycle_length <= 0){
    stop('Selected cell type must have a positive finite cell_cycle_length.')
  }

  substitution_matrix <- physicell_substitution_matrix(
    state_params$mt_substitution_model,
    state_params$mt_sub_model_params,
    mitochondrial_reference
  )
  insertion_probability <- as.numeric(
    state_params$mt_bg_insertion_prob_per_division
  )
  deletion_probability <- as.numeric(
    state_params$mt_bg_deletion_prob_per_division
  )
  if(length(insertion_probability) != 1 || !is.finite(insertion_probability) ||
     insertion_probability < 0 ||
     length(deletion_probability) != 1 || !is.finite(deletion_probability) ||
     deletion_probability < 0){
    stop('Mitochondrial background indel probabilities must be finite and non-negative.')
  }

  bases <- c('A', 'G', 'C', 'T')
  from_indices <- match(mitochondrial_reference, bases)
  substitution_hazards <- matrix(
    0,
    nrow = length(mitochondrial_reference),
    ncol = 4,
    dimnames = list(NULL, bases)
  )
  for(base_index in seq_along(bases)){
    positions <- which(from_indices == base_index)
    if(length(positions) > 0){
      substitution_hazards[positions, ] <- matrix(
        rep(
          physicell_probability_hazard(
            substitution_matrix[base_index, ],
            cell_cycle_length
          ),
          each = length(positions)
        ),
        nrow = length(positions),
        byrow = FALSE
      )
    }
  }
  insertion_hazard <- physicell_probability_hazard(
    insertion_probability,
    cell_cycle_length
  )
  deletion_hazard <- physicell_probability_hazard(
    deletion_probability,
    cell_cycle_length
  )
  position_hazards <- rowSums(substitution_hazards) +
    insertion_hazard + deletion_hazard

  invariant_fraction <- as.numeric(cell_params$mt_invariant_sites)
  if(length(invariant_fraction) != 1 || !is.finite(invariant_fraction) ||
     invariant_fraction < 0 || invariant_fraction > 1){
    stop('mt_invariant_sites must be a fraction between zero and one.')
  }
  num_invariants <- round(invariant_fraction * length(mitochondrial_reference))
  invariant_positions <- if(num_invariants > 0){
    sample(
      seq_along(mitochondrial_reference),
      num_invariants,
      replace = FALSE
    )
  } else{
    integer()
  }
  if(length(invariant_positions) > 0){
    substitution_hazards[invariant_positions, ] <- 0
    position_hazards[invariant_positions] <- 0
  }

  list(
    substitution_matrix = substitution_matrix,
    substitution_hazards = substitution_hazards,
    insertion_hazard = unname(insertion_hazard),
    deletion_hazard = unname(deletion_hazard),
    position_hazards = position_hazards,
    cumulative_position_hazards = cumsum(position_hazards),
    total_hazard = sum(position_hazards),
    invariant_positions = invariant_positions,
    cell_cycle_length = cell_cycle_length
  )
}

prepare_physicell_mito_model <- function(params,
                                         cell_type = NULL,
                                         genomes_per_cell = 8,
                                         seed = NULL){
  if(!is.list(params)){
    stop('params must be the parsed remote_mito JSON parameter list.')
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.null(cell_type)){
    cell_type <- as.character(params$cell_type_dict$founder_cell_type)
  }
  mitochondrial_length <- as.integer(params$mito_genome_length)
  genomes_per_cell <- as.integer(genomes_per_cell)
  if(length(mitochondrial_length) != 1 || is.na(mitochondrial_length) ||
     mitochondrial_length < 1){
    stop('mito_genome_length must be one positive integer.')
  }
  if(length(genomes_per_cell) != 1 || is.na(genomes_per_cell) ||
     genomes_per_cell < 1){
    stop('genomes_per_cell must be one positive integer.')
  }

  bases <- c('A', 'G', 'C', 'T')
  mitochondrial_reference <- sample(
    bases,
    mitochondrial_length,
    replace = TRUE
  )
  rate_sets <- list(
    uninduced_editing_params = physicell_mito_rate_set(
      params,
      cell_type,
      'uninduced_editing_params',
      mitochondrial_reference
    ),
    induced_editing_params = physicell_mito_rate_set(
      params,
      cell_type,
      'induced_editing_params',
      mitochondrial_reference
    )
  )

  editing_induction_time <- as.numeric(params$editing_induction$timepoint)
  if(length(editing_induction_time) != 1 || !is.finite(editing_induction_time)){
    editing_induction_time <- Inf
  }
  adapter <- params$physicell_adapter
  retain_internal_profiles <- if(is.null(adapter$retain_internal_profiles)){
    TRUE
  } else{
    isTRUE(adapter$retain_internal_profiles)
  }

  list(
    mitochondrial_reference = mitochondrial_reference,
    mitochondrial_length = mitochondrial_length,
    genomes_per_cell = genomes_per_cell,
    cell_type = cell_type,
    editing_induction_time = editing_induction_time,
    retain_internal_profiles = retain_internal_profiles,
    rate_sets = rate_sets
  )
}

initialize_physicell_mito_profile <- function(model, params){
  genomes <- replicate(
    model$genomes_per_cell,
    setNames(numeric(), character()),
    simplify = FALSE
  )
  site_fraction <- as.numeric(params$baseline_heteroplasmy_sites_frac)
  penetrance_parameters <- as.numeric(unlist(
    params$baseline_heteroplasmy_variant_frac_dist,
    use.names = FALSE
  ))
  transition_probability <- as.numeric(
    params$heteroplasmy_variant_transition_prob
  )
  if(length(site_fraction) != 1 || !is.finite(site_fraction) ||
     site_fraction < 0 || site_fraction > 1 ||
     length(penetrance_parameters) != 2 ||
     any(!is.finite(penetrance_parameters)) ||
     any(penetrance_parameters <= 0) ||
     length(transition_probability) != 1 ||
     !is.finite(transition_probability) ||
     transition_probability < 0 || transition_probability > 1){
    stop('Invalid baseline mitochondrial heteroplasmy parameters.')
  }

  num_heteroplasmic_sites <- round(site_fraction * model$mitochondrial_length)
  heteroplasmic_sites <- if(num_heteroplasmic_sites > 0){
    sample(
      seq_len(model$mitochondrial_length),
      num_heteroplasmic_sites,
      replace = FALSE
    )
  } else{
    integer()
  }
  bases <- c('A', 'G', 'C', 'T')
  transition_destination <- c(2L, 1L, 4L, 3L)
  transversion_destinations <- list(
    c(3L, 4L),
    c(3L, 4L),
    c(1L, 2L),
    c(1L, 2L)
  )
  initial_events <- list()
  event_index <- 0L

  for(position in heteroplasmic_sites){
    penetrance <- stats::rbeta(
      1,
      penetrance_parameters[1],
      penetrance_parameters[2]
    )
    variant_genomes <- which(
      stats::runif(model$genomes_per_cell) < penetrance
    )
    if(length(variant_genomes) == 0){
      next
    }
    from_base <- match(model$mitochondrial_reference[position], bases)
    for(genome_index in variant_genomes){
      if(stats::runif(1) < transition_probability){
        allele <- transition_destination[from_base]
        event_type <- 'initial_transition'
      } else{
        allele <- sample(transversion_destinations[[from_base]], 1)
        event_type <- 'initial_transversion'
      }
      genomes[[genome_index]][as.character(position)] <- allele
      event_index <- event_index + 1L
      initial_events[[event_index]] <- data.frame(
        genome = genome_index,
        position = position,
        event = event_type,
        reference = bases[from_base],
        alternate = bases[allele],
        allele = allele,
        stringsAsFactors = FALSE
      )
    }
  }

  events <- if(length(initial_events) == 0){
    data.frame(
      genome = integer(),
      position = integer(),
      event = character(),
      reference = character(),
      alternate = character(),
      allele = numeric(),
      stringsAsFactors = FALSE
    )
  } else{
    do.call(rbind, initial_events)
  }
  list(profile = genomes, events = events)
}

inherit_physicell_mito_profile <- function(parent_profile, genomes_per_cell){
  inherited_indices <- sample(
    seq_along(parent_profile),
    genomes_per_cell,
    replace = TRUE
  )
  lapply(inherited_indices, function(index) parent_profile[[index]])
}

sample_physicell_mito_position <- function(rate_set, excluded_positions){
  if(rate_set$total_hazard <= 0){
    return(NA_integer_)
  }
  repeat{
    position <- findInterval(
      stats::runif(1, 0, rate_set$total_hazard),
      rate_set$cumulative_position_hazards
    ) + 1L
    if(position <= length(rate_set$position_hazards) &&
       rate_set$position_hazards[position] > 0 &&
       !(position %in% excluded_positions)){
      return(position)
    }
  }
}

mutate_physicell_mito_segment <- function(profile,
                                          duration,
                                          rate_set,
                                          model,
                                          segment_start = 0){
  empty_events <- data.frame(
    genome = integer(),
    position = integer(),
    event = character(),
    reference = character(),
    alternate = character(),
    allele = numeric(),
    event_time = numeric(),
    stringsAsFactors = FALSE
  )
  if(duration <= 0 || rate_set$total_hazard <= 0){
    return(list(profile = profile, events = empty_events))
  }

  bases <- c('A', 'G', 'C', 'T')
  event_rows <- list()
  event_index <- 0L
  for(genome_index in seq_along(profile)){
    existing_positions <- as.integer(names(profile[[genome_index]]))
    existing_hazard <- if(length(existing_positions) == 0){
      0
    } else{
      sum(rate_set$position_hazards[existing_positions])
    }
    available_hazard <- max(rate_set$total_hazard - existing_hazard, 0)
    num_events <- stats::rpois(1, available_hazard * duration)
    max_events <- model$mitochondrial_length - length(existing_positions)
    num_events <- min(num_events, max_events)
    if(num_events == 0){
      next
    }

    selected_positions <- integer()
    for(unused in seq_len(num_events)){
      position <- sample_physicell_mito_position(
        rate_set,
        c(existing_positions, selected_positions)
      )
      if(is.na(position)){
        break
      }
      selected_positions <- c(selected_positions, position)
    }

    for(position in selected_positions){
      event_hazards <- c(
        setNames(
          rate_set$substitution_hazards[position, ],
          paste0('sub_', seq_along(bases))
        ),
        insertion = rate_set$insertion_hazard,
        deletion = rate_set$deletion_hazard
      )
      event_hazards <- event_hazards[event_hazards > 0]
      selected_event <- sample(
        names(event_hazards),
        1,
        prob = event_hazards,
        replace = FALSE
      )
      if(grepl('^sub_[1-4]$', selected_event)){
        allele <- as.integer(sub('^sub_', '', selected_event))
        event_type <- 'substitution'
        alternate <- bases[allele]
      } else if(selected_event == 'deletion'){
        allele <- -1
        event_type <- 'deletion'
        alternate <- '-'
      } else{
        inserted_base <- sample(seq_along(bases), 1)
        allele <- as.numeric(paste0('0.', inserted_base))
        event_type <- 'insertion'
        alternate <- paste0('+', bases[inserted_base])
      }
      profile[[genome_index]][as.character(position)] <- allele
      event_index <- event_index + 1L
      event_rows[[event_index]] <- data.frame(
        genome = genome_index,
        position = position,
        event = event_type,
        reference = model$mitochondrial_reference[position],
        alternate = alternate,
        allele = allele,
        event_time = stats::runif(
          1,
          min = segment_start,
          max = segment_start + duration
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  events <- if(length(event_rows) == 0){
    empty_events
  } else{
    do.call(rbind, event_rows)
  }
  list(profile = profile, events = events)
}

mutate_physicell_mito_branch <- function(profile,
                                         start_time,
                                         end_time,
                                         model,
                                         editing_state = 'auto'){
  allowed_states <- c('auto', 'induced', 'uninduced')
  if(!(editing_state %in% allowed_states)){
    stop(sprintf(
      'editing_state must be one of: %s.',
      paste(allowed_states, collapse = ', ')
    ))
  }
  segments <- if(editing_state == 'induced'){
    data.frame(
      start = start_time,
      end = end_time,
      state = 'induced_editing_params',
      stringsAsFactors = FALSE
    )
  } else if(editing_state == 'uninduced'){
    data.frame(
      start = start_time,
      end = end_time,
      state = 'uninduced_editing_params',
      stringsAsFactors = FALSE
    )
  } else if(end_time <= model$editing_induction_time){
    data.frame(
      start = start_time,
      end = end_time,
      state = 'uninduced_editing_params',
      stringsAsFactors = FALSE
    )
  } else if(start_time >= model$editing_induction_time){
    data.frame(
      start = start_time,
      end = end_time,
      state = 'induced_editing_params',
      stringsAsFactors = FALSE
    )
  } else{
    data.frame(
      start = c(start_time, model$editing_induction_time),
      end = c(model$editing_induction_time, end_time),
      state = c('uninduced_editing_params', 'induced_editing_params'),
      stringsAsFactors = FALSE
    )
  }

  all_events <- list()
  for(segment_index in seq_len(nrow(segments))){
    mutation <- mutate_physicell_mito_segment(
      profile,
      duration = segments$end[segment_index] - segments$start[segment_index],
      rate_set = model$rate_sets[[segments$state[segment_index]]],
      model = model,
      segment_start = segments$start[segment_index]
    )
    profile <- mutation$profile
    if(nrow(mutation$events) > 0){
      mutation$events$segment_start <- segments$start[segment_index]
      mutation$events$segment_end <- segments$end[segment_index]
      mutation$events$editing_state <- segments$state[segment_index]
      all_events[[length(all_events) + 1L]] <- mutation$events
    }
  }

  events <- if(length(all_events) == 0){
    data.frame(
      genome = integer(),
      position = integer(),
      event = character(),
      reference = character(),
      alternate = character(),
      allele = numeric(),
      event_time = numeric(),
      segment_start = numeric(),
      segment_end = numeric(),
      editing_state = character(),
      stringsAsFactors = FALSE
    )
  } else{
    do.call(rbind, all_events)
  }
  list(profile = profile, events = events)
}

simulate_mito_on_physicell_lineage <- function(nodes,
                                               model,
                                               params,
                                               editing_state = 'auto',
                                               terminal_physicell_ids = NULL,
                                               seed = 1,
                                               show_progress = FALSE,
                                               progress_updates = 20L){
  if(!is.data.frame(nodes) || !all(c(
    'node_id', 'physicell_id', 'parent_node_id', 'birth_time',
    'end_time', 'is_terminal'
  ) %in% names(nodes))){
    stop('nodes is not a valid PhysiCell lineage node table.')
  }
  set.seed(seed)
  profiles <- list()
  all_events <- list()
  remaining_children <- table(
    nodes$parent_node_id[!is.na(nodes$parent_node_id)]
  )
  report_progress <- new_physicell_progress_reporter(
    total = nrow(nodes),
    label = 'mitochondrial lineage tracing',
    enabled = show_progress,
    updates = progress_updates
  )
  report_progress(0L)

  for(node_index in seq_len(nrow(nodes))){
    node <- nodes[node_index, , drop = FALSE]
    if(is.na(node$parent_node_id)){
      initialized <- initialize_physicell_mito_profile(model, params)
      start_profile <- initialized$profile
      if(nrow(initialized$events) > 0){
        initialized$events$segment_start <- node$birth_time
        initialized$events$segment_end <- node$birth_time
        initialized$events$event_time <- node$birth_time
        initialized$events$editing_state <- 'founder_heteroplasmy'
        initialized$events$node_id <- node$node_id
        initialized$events$physicell_id <- node$physicell_id
        initialized$events$parent_node_id <- NA_character_
        initialized$events$branch_start <- node$birth_time
        initialized$events$branch_end <- node$birth_time
        all_events[[length(all_events) + 1L]] <- initialized$events
      }
    } else{
      parent_node_id <- as.character(node$parent_node_id)
      parent_profile <- profiles[[parent_node_id]]
      if(is.null(parent_profile)){
        stop(sprintf(
          'Parent mitochondrial profile %s was not available before child %s.',
          node$parent_node_id,
          node$node_id
        ))
      }
      remaining_children[parent_node_id] <-
        remaining_children[parent_node_id] - 1L
      if(!isTRUE(model$retain_internal_profiles) &&
         remaining_children[parent_node_id] <= 0){
        profiles[[parent_node_id]] <- NULL
      }
      start_profile <- inherit_physicell_mito_profile(
        parent_profile,
        model$genomes_per_cell
      )
    }

    mutation <- mutate_physicell_mito_branch(
      start_profile,
      node$birth_time,
      node$end_time,
      model,
      editing_state
    )
    profiles[[node$node_id]] <- mutation$profile
    if(nrow(mutation$events) > 0){
      mutation$events$node_id <- node$node_id
      mutation$events$physicell_id <- node$physicell_id
      mutation$events$parent_node_id <- node$parent_node_id
      mutation$events$branch_start <- node$birth_time
      mutation$events$branch_end <- node$end_time
      all_events[[length(all_events) + 1L]] <- mutation$events
    }
    report_progress(node_index)
  }

  terminal_rows <- nodes[nodes$is_terminal, , drop = FALSE]
  if(!is.null(terminal_physicell_ids)){
    terminal_physicell_ids <- unique(
      physicell_id(terminal_physicell_ids, 'terminal_physicell_ids')
    )
    unknown_ids <- setdiff(terminal_physicell_ids, terminal_rows$physicell_id)
    if(length(unknown_ids) > 0){
      stop(sprintf(
        'Live-cell ID(s) absent from the reconstructed lineage: %s.',
        paste(unknown_ids, collapse = ', ')
      ))
    }
    terminal_rows <- terminal_rows[
      terminal_rows$physicell_id %in% terminal_physicell_ids,
      ,
      drop = FALSE
    ]
  }
  if(nrow(terminal_rows) == 0){
    stop('No terminal cells remain for mitochondrial output.')
  }

  events <- if(length(all_events) == 0){
    data.frame(
      genome = integer(),
      position = integer(),
      event = character(),
      reference = character(),
      alternate = character(),
      allele = numeric(),
      event_time = numeric(),
      segment_start = numeric(),
      segment_end = numeric(),
      editing_state = character(),
      node_id = character(),
      physicell_id = character(),
      parent_node_id = character(),
      branch_start = numeric(),
      branch_end = numeric(),
      stringsAsFactors = FALSE
    )
  } else{
    do.call(rbind, all_events)
  }
  if(nrow(events) > 0){
    events <- events[
      order(
        events$event_time,
        events$node_id,
        events$genome,
        events$position
      ),
      ,
      drop = FALSE
    ]
    rownames(events) <- NULL
  }

  list(
    nodes = nodes,
    profiles = profiles,
    terminal_nodes = terminal_rows,
    mutation_events = events,
    editing_state = editing_state,
    seed = seed
  )
}

physicell_mito_variant_fractions <- function(simulation){
  terminal_nodes <- simulation$terminal_nodes
  rows <- list()
  row_index <- 0L
  for(cell_index in seq_len(nrow(terminal_nodes))){
    node <- terminal_nodes[cell_index, , drop = FALSE]
    profile <- simulation$profiles[[node$node_id]]
    genome_rows <- lapply(seq_along(profile), function(genome_index){
      variants <- profile[[genome_index]]
      if(length(variants) == 0){
        return(NULL)
      }
      data.frame(
        genome = genome_index,
        position = as.integer(names(variants)),
        allele = as.numeric(variants),
        stringsAsFactors = FALSE
      )
    })
    genome_rows <- Filter(Negate(is.null), genome_rows)
    if(length(genome_rows) == 0){
      next
    }
    variants <- do.call(rbind, genome_rows)
    counts <- stats::aggregate(
      variants$genome,
      by = list(position = variants$position, allele = variants$allele),
      FUN = length
    )
    names(counts)[3] <- 'variant_genomes'
    counts$variant_fraction <- counts$variant_genomes / length(profile)
    counts$sample_id <- paste0('cell_', node$physicell_id)
    counts$physicell_id <- node$physicell_id
    counts$node_id <- node$node_id
    row_index <- row_index + 1L
    rows[[row_index]] <- counts[, c(
      'sample_id', 'physicell_id', 'node_id', 'position', 'allele',
      'variant_genomes', 'variant_fraction'
    )]
  }

  if(length(rows) == 0){
    return(data.frame(
      sample_id = character(),
      physicell_id = character(),
      node_id = character(),
      position = integer(),
      allele = numeric(),
      variant_genomes = integer(),
      variant_fraction = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

physicell_mito_haplotype_sequence <- function(genome, reference){
  bases <- c('A', 'G', 'C', 'T')
  sequence <- reference
  if(length(genome) == 0){
    return(paste0(sequence, collapse = ''))
  }
  positions <- as.integer(names(genome))
  alleles <- as.numeric(genome)
  substitutions <- which(alleles %% 1 == 0 & alleles > 0)
  deletions <- which(alleles == -1)
  insertions <- which(alleles %% 1 != 0)
  if(length(substitutions) > 0){
    sequence[positions[substitutions]] <- bases[as.integer(alleles[substitutions])]
  }
  if(length(deletions) > 0){
    sequence[positions[deletions]] <- ''
  }
  if(length(insertions) > 0){
    inserted_codes <- as.integer(round(abs(alleles[insertions] %% 1) * 10)) %% 10
    sequence[positions[insertions]] <- paste0(
      sequence[positions[insertions]],
      bases[inserted_codes]
    )
  }
  paste0(sequence, collapse = '')
}

write_physicell_mito_outputs <- function(simulation,
                                         model,
                                         output_dir,
                                         write_fasta = FALSE,
                                         show_progress = FALSE,
                                         compress_csv = TRUE){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  terminal_nodes <- simulation$terminal_nodes
  terminal_profiles <- simulation$profiles[terminal_nodes$node_id]
  names(terminal_profiles) <- paste0('cell_', terminal_nodes$physicell_id)
  physicell_log_stage(
    sprintf(
      paste(
        'Mitochondrial output started for %s terminal cells and',
        '%s mutation events.'
      ),
      format(
        nrow(terminal_nodes),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      ),
      format(
        nrow(simulation$mutation_events),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      )
    ),
    enabled = show_progress
  )
  physicell_log_stage(
    'Mitochondrial output: summarizing terminal-cell heteroplasmy.',
    enabled = show_progress
  )
  variant_fractions <- physicell_mito_variant_fractions(simulation)

  physicell_log_stage(
    'Mitochondrial output: saving profiles and mutation-event table.',
    enabled = show_progress
  )
  saveRDS(
    terminal_profiles,
    file.path(output_dir, 'mitochondrial_profiles.rds')
  )
  write_physicell_csv(
    simulation$mutation_events,
    file.path(output_dir, 'mitochondrial_mutation_events.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  physicell_log_stage(
    'Mitochondrial output: writing variant fractions and reference FASTA.',
    enabled = show_progress
  )
  write_physicell_csv(
    variant_fractions,
    file.path(output_dir, 'mitochondrial_variant_fractions.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  writeLines(
    c(
      '>remote_mito_mitochondrial_reference',
      paste0(model$mitochondrial_reference, collapse = '')
    ),
    file.path(output_dir, 'mitochondrial_reference.fasta')
  )

  sparse_scores <- NULL
  if(requireNamespace('Matrix', quietly = TRUE)){
    physicell_log_stage(
      'Mitochondrial output: assembling sparse variant-fraction matrix.',
      enabled = show_progress
    )
    if(nrow(variant_fractions) == 0){
      sparse_scores <- Matrix::sparseMatrix(
        i = integer(),
        j = integer(),
        dims = c(nrow(terminal_nodes), 0),
        dimnames = list(
          paste0('cell_', terminal_nodes$physicell_id),
          character()
        )
      )
    } else{
      variant_fractions$variant_name <- paste0(
        'mt_', variant_fractions$position, '_', variant_fractions$allele
      )
      variant_names <- unique(variant_fractions$variant_name)
      sample_names <- paste0('cell_', terminal_nodes$physicell_id)
      sparse_scores <- Matrix::sparseMatrix(
        i = match(variant_fractions$sample_id, sample_names),
        j = match(variant_fractions$variant_name, variant_names),
        x = variant_fractions$variant_fraction,
        dims = c(length(sample_names), length(variant_names)),
        dimnames = list(sample_names, variant_names)
      )
    }
    saveRDS(
      sparse_scores,
      file.path(output_dir, 'mitochondrial_variant_fraction_matrix.rds')
    )
    physicell_log_stage(
      'Mitochondrial output: sparse variant-fraction matrix written.',
      enabled = show_progress
    )
  }

  if(isTRUE(write_fasta)){
    physicell_log_stage(
      'Mitochondrial output: writing sampled haplotype FASTA.',
      enabled = show_progress
    )
    fasta_connection <- file(
      file.path(output_dir, 'mitochondrial_sampled_haplotypes.fasta'),
      open = 'w'
    )
    on.exit(close(fasta_connection), add = TRUE)
    for(cell_name in names(terminal_profiles)){
      profile <- terminal_profiles[[cell_name]]
      sampled_genome <- sample(seq_along(profile), 1)
      writeLines(
        c(
          paste0('>', cell_name, '_mt_genome_', sampled_genome),
          physicell_mito_haplotype_sequence(
            profile[[sampled_genome]],
            model$mitochondrial_reference
          )
        ),
        fasta_connection
      )
    }
    close(fasta_connection)
    on.exit(NULL, add = FALSE)
  }

  physicell_log_stage(
    'Mitochondrial output: writing manifest.',
    enabled = show_progress
  )
  manifest <- data.frame(
    property = c(
      'num_sampled_terminal_cells',
      'mitochondrial_length',
      'genomes_per_cell',
      'num_mitochondrial_events',
      'num_variant_observations',
      'cell_type',
      'editing_state',
      'random_seed',
      'sampled_haplotype_fasta'
    ),
    value = c(
      nrow(terminal_nodes),
      model$mitochondrial_length,
      model$genomes_per_cell,
      nrow(simulation$mutation_events),
      nrow(variant_fractions),
      model$cell_type,
      simulation$editing_state,
      simulation$seed,
      isTRUE(write_fasta)
    ),
    stringsAsFactors = FALSE
  )
  write_physicell_csv(
    manifest,
    file.path(output_dir, 'mitochondrial_manifest.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  physicell_log_stage(
    'Mitochondrial output complete.',
    enabled = show_progress
  )

  invisible(list(
    terminal_profiles = terminal_profiles,
    variant_fractions = variant_fractions,
    sparse_scores = sparse_scores
  ))
}
