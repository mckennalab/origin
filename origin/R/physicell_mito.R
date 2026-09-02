# Scalable mitochondrial lineage-recording replay for event-resolved PhysiCell
# trees. Source physicell_lineage.R before using these functions.

#' Build the mitochondrial hazard set for one cell type and editing state
#'
#' Reads the per-cell-cycle (per-division) probabilities configured for
#' `cell_type` in `editing_state` and converts every one of them into a
#' continuous-time hazard with `physicell_probability_hazard()`, that is
#' `-log(1 - p) / cell_cycle_length`. Every `*_hazard` field of the result is
#' therefore a rate per unit of lineage time (the units `cell_cycle_length` is
#' expressed in); only `substitution_matrix` is left in per-division
#' probability space. Positions selected as invariant have their substitution
#' hazards and their total position hazard zeroed, so they can never be drawn.
#'
#' @export
#' @param params Parsed remote_mito JSON parameter list.
#' @param cell_type Name of a cell type present in
#'   `params$cell_type_dict$cell_type_params`; it must carry a positive finite
#'   `cell_cycle_length` and an `mt_invariant_sites` fraction in `[0, 1]`.
#' @param editing_state Name of the rate block to read, either
#'   `'uninduced_editing_params'` or `'induced_editing_params'`; it supplies
#'   `mt_substitution_model`, `mt_sub_model_params`,
#'   `mt_bg_insertion_prob_per_division`, and
#'   `mt_bg_deletion_prob_per_division`, the last two of which must be finite
#'   and non-negative.
#' @param mitochondrial_reference Character vector of reference bases, one per
#'   mt genome position, over `A`, `G`, `C`, `T`.
#' @return A named list: `substitution_matrix` (4x4 per-division destination
#'   probabilities, rows and columns ordered `A, G, C, T`),
#'   `substitution_hazards` (position-by-4 destination hazards),
#'   `insertion_hazard` and `deletion_hazard` (scalars, uniform across
#'   positions), `position_hazards` (per-position total event hazard),
#'   `cumulative_position_hazards` (its `cumsum`, cached for weighted
#'   sampling), `total_hazard`, `invariant_positions`, and `cell_cycle_length`.
#' @section Side effects: Draws the invariant-site set with `sample()`, so it
#'   consumes the random-number stream.
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

#' Prepare the mitochondrial replay model for one cell type
#'
#' Draws a uniformly random mt reference sequence of `params$mito_genome_length`
#' bases and builds both the uninduced and the induced hazard set against that
#' one reference, so the two states share a reference but each converts its own
#' per-division probabilities into hazards.
#'
#' @export
#' @param params Parsed remote_mito JSON parameter list.
#' @param cell_type Cell type whose rate blocks are used; defaults to
#'   `params$cell_type_dict$founder_cell_type`.
#' @param genomes_per_cell Fixed number of mt genomes every cell carries and
#'   that each daughter resamples at division; one positive integer.
#' @param seed Optional seed set before any random draw in this function.
#' @return A named list with `mitochondrial_reference` (character vector of
#'   bases), `mitochondrial_length`, `genomes_per_cell`, `cell_type`,
#'   `editing_induction_time` (`Inf` when
#'   `params$editing_induction$timepoint` is missing or non-finite),
#'   `retain_internal_profiles` (from
#'   `params$physicell_adapter$retain_internal_profiles`, default `TRUE`),
#'   `variant_selection_coefficient` (from
#'   `non_mendelian_selection$mitochondrial_variant_coefficient`), and
#'   `rate_sets`, a list keyed `uninduced_editing_params` and
#'   `induced_editing_params`.
#' @section Side effects: Calls `set.seed()` when `seed` is supplied, and
#'   consumes the random-number stream for the reference sequence and again for
#'   each rate set, which samples its own invariant sites independently.
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
  variant_selection_coefficient <- non_mendelian_selection_coefficient(
    params,
    'mitochondrial_variant_coefficient'
  )

  list(
    mitochondrial_reference = mitochondrial_reference,
    mitochondrial_length = mitochondrial_length,
    genomes_per_cell = genomes_per_cell,
    cell_type = cell_type,
    editing_induction_time = editing_induction_time,
    retain_internal_profiles = retain_internal_profiles,
    variant_selection_coefficient = variant_selection_coefficient,
    rate_sets = rate_sets
  )
}

#' Draw a founder cell's mt genomes and its baseline heteroplasmy
#'
#' Every genome starts identical to the reference. A fraction
#' `params$baseline_heteroplasmy_sites_frac` of positions is chosen without
#' replacement; each chosen site draws a penetrance from the Beta distribution
#' whose two shape parameters are given by
#' `params$baseline_heteroplasmy_variant_frac_dist`, and each genome
#' independently carries the variant with that probability. A
#' carried variant is a transition (`A<->G`, `C<->T`) with probability
#' `params$heteroplasmy_variant_transition_prob`, and otherwise a transversion
#' drawn uniformly from the two purine/pyrimidine alternatives.
#'
#' @export
#' @param model Prepared model from `prepare_physicell_mito_model()`; supplies
#'   `genomes_per_cell`, `mitochondrial_length`, and the reference bases.
#' @param params Parsed remote_mito JSON parameter list; the site fraction, the
#'   two Beta parameters (both positive), and the transition probability are
#'   validated here.
#' @return A list with `profile`, a list of `genomes_per_cell` sparse genomes
#'   (named numeric vectors whose names are positions as characters and whose
#'   values are the substituted base codes `1`-`4`), and `events`, a data frame
#'   with columns `genome`, `position`, `event` (`initial_transition` or
#'   `initial_transversion`), `reference`, `alternate`, `allele`; zero rows
#'   when nothing was drawn.
#' @section Side effects: Consumes the random-number stream.
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

#' Resample a daughter cell's mt genomes through the division bottleneck
#'
#' Draws `genomes_per_cell` genomes from the parent with replacement, so the
#' two division products drift independently. With a selection coefficient
#' `s > 0`, a parental genome carrying `k` recorded variants is weighted
#' `(1 - s)^k`, penalizing mutation burden; `s = 0` is neutral uniform
#' sampling. If every weight collapses to zero the draw falls back to a uniform
#' choice among the lowest-burden genomes.
#'
#' @export
#' @param parent_profile Parent cell's list of sparse mt genomes.
#' @param genomes_per_cell Number of genomes the daughter receives.
#' @param variant_selection_coefficient One finite value in `[0, 1]`.
#' @return A list of `genomes_per_cell` genomes copied from `parent_profile`;
#'   the same parental genome may be drawn more than once.
#' @section Side effects: Consumes the random-number stream.
inherit_physicell_mito_profile <- function(parent_profile,
                                           genomes_per_cell,
                                           variant_selection_coefficient = 0){
  variant_selection_coefficient <- suppressWarnings(
    as.numeric(variant_selection_coefficient)
  )
  if(length(variant_selection_coefficient) != 1L ||
     !is.finite(variant_selection_coefficient) ||
     variant_selection_coefficient < 0 ||
     variant_selection_coefficient > 1){
    stop('variant_selection_coefficient must be one value in [0, 1].')
  }
  if(variant_selection_coefficient == 0){
    inherited_indices <- sample(
      seq_along(parent_profile),
      genomes_per_cell,
      replace = TRUE
    )
  } else{
    variant_burden <- lengths(parent_profile)
    inheritance_weights <-
      (1 - variant_selection_coefficient)^variant_burden
    if(!any(inheritance_weights > 0)){
      inheritance_weights <- as.numeric(
        variant_burden == min(variant_burden)
      )
    }
    inherited_indices <- sample(
      seq_along(parent_profile),
      genomes_per_cell,
      replace = TRUE,
      prob = inheritance_weights
    )
  }
  lapply(inherited_indices, function(index) parent_profile[[index]])
}

#' Sample one mutable position in proportion to its total event hazard
#'
#' Inverse-CDF sampling against the rate set's cached
#' `cumulative_position_hazards`, retried until the draw lands on a position
#' that has non-zero hazard and is not excluded.
#'
#' @export
#' @param rate_set A rate set from `physicell_mito_rate_set()`; hazards, not
#'   per-division probabilities.
#' @param excluded_positions Integer positions to reject, normally the sites a
#'   genome has already mutated plus those already selected for this segment.
#' @return One integer position, or `NA_integer_` when the rate set has no
#'   hazard at all.
#' @note The retry loop does not terminate if `excluded_positions` already
#'   covers every position with non-zero hazard, so callers must bound the
#'   number of positions they request.
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

#' Mutate mt genomes over one segment of constant editing state
#'
#' Works entirely in hazard space. For each genome the hazard of the positions
#' it has already mutated is subtracted from the rate set's `total_hazard`,
#' because recording sites are irreversible; the number of new events is then
#' `rpois(1, available_hazard * duration)`, capped at the number of positions
#' the genome has left. Positions are drawn without replacement in proportion
#' to their hazard, and each event's type is drawn among the four destination
#' substitution hazards plus the insertion and deletion hazards at that
#' position. Every event receives a time drawn uniformly inside the segment.
#'
#' @export
#' @param profile List of sparse mt genomes to mutate.
#' @param duration Segment length, in the same time units as
#'   `cell_cycle_length`; a non-positive duration is a no-op, as is a rate set
#'   with zero total hazard.
#' @param rate_set Hazards for the editing state that holds over the segment.
#' @param model Prepared model, read for `mitochondrial_length` and for the
#'   reference base recorded on each event row.
#' @param segment_start Absolute start time of the segment; used only as the
#'   lower bound for the sampled `event_time` values.
#' @return A list with the updated `profile` and an `events` data frame with
#'   columns `genome`, `position`, `event` (`substitution`, `insertion`, or
#'   `deletion`), `reference`, `alternate`, `allele`, `event_time`.
#' @details Alleles use the repo encoding: `1`-`4` for a substitution to `A`,
#'   `G`, `C`, `T`, `-1` for a deletion, and the decimal `0.b` for an insertion
#'   of base code `b`. The `alternate` column carries the human-readable form
#'   instead: the base letter, `-`, or `+<base>`.
#' @section Side effects: Consumes the random-number stream.
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

#' Mutate one lineage branch, splitting it at editing induction
#'
#' `editing_state` `'induced'` or `'uninduced'` forces a single rate set over
#' the whole branch. `'auto'` compares the branch against
#' `model$editing_induction_time` and splits `[start_time, end_time]` into an
#' uninduced segment followed by an induced one when induction falls strictly
#' inside the branch; branches wholly before or after induction stay single
#' segments.
#'
#' @export
#' @param profile List of sparse mt genomes entering the branch.
#' @param start_time Branch start, normally the node's birth time.
#' @param end_time Branch end, normally the node's end time.
#' @param model Prepared model supplying `rate_sets` and
#'   `editing_induction_time`.
#' @param editing_state One of `'auto'`, `'induced'`, `'uninduced'`.
#' @return A list with the branch-final `profile` and an `events` data frame
#'   holding the columns of `mutate_physicell_mito_segment()` plus
#'   `segment_start`, `segment_end`, and `editing_state`, the name of the rate
#'   set that produced each event.
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

#' Replay mitochondrial evolution over a reconstructed lineage
#'
#' Walks `nodes` in row order, which must already be topological so that a
#' parent is processed before its children. A founder row (`parent_node_id` is
#' `NA`) receives a fresh baseline-heteroplasmy profile; every other row
#' inherits its parent's genomes through the division bottleneck and is then
#' mutated over its own `[birth_time, end_time]` interval. Rows whose `origin`
#' is `induction_continuation` mark a rate-state change rather than a cell
#' division, so they inherit the parent profile without a bottleneck.
#'
#' @export
#' @param nodes Lineage node table; must be a data frame with at least
#'   `node_id`, `physicell_id`, `parent_node_id`, `birth_time`, `end_time`,
#'   `is_terminal`. Optional `cell_type`, `editing_state`, and `origin` columns
#'   are honored when present.
#' @param model One prepared model, or a named collection of prepared models
#'   keyed by cell type; a collection requires a `cell_type` column on `nodes`
#'   whose every value has a model.
#' @param params Parsed remote_mito JSON parameter list, forwarded to founder
#'   initialization.
#' @param editing_state Fallback state policy (`'auto'`, `'induced'`,
#'   `'uninduced'`) applied when `nodes` has no `editing_state` column.
#' @param terminal_physicell_ids Optional PhysiCell IDs restricting the sampled
#'   terminal set; every ID must appear among the terminal nodes.
#' @param seed Seed set once before the replay begins.
#' @param show_progress Whether to emit progress and ETA messages.
#' @param progress_updates Number of progress checkpoints requested.
#' @return A named list with `nodes` (as supplied), `profiles` (genome lists
#'   keyed by `node_id`), `terminal_nodes` (the sampled terminal rows),
#'   `mutation_events` (every event row plus `node_id`, `physicell_id`,
#'   `parent_node_id`, `branch_start`, `branch_end`, sorted by event time),
#'   `editing_state`, and `seed`. Founder heteroplasmy rows appear with
#'   `editing_state` `'founder_heteroplasmy'` and all their times set to the
#'   founder's birth time.
#' @details When `retain_internal_profiles` is off, a parent's profile is
#'   dropped as soon as its last child has been processed, so `profiles` then
#'   contains only the nodes still needed downstream. For a model collection
#'   that flag is read from the first model in the collection.
#' @section Side effects: Calls `set.seed(seed)` and consumes the
#'   random-number stream.
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
  model_is_collection <- is.null(model$mitochondrial_length)
  if(model_is_collection){
    if(length(model) == 0 || is.null(names(model)) ||
       any(!nzchar(names(model))) ||
       any(vapply(
         model,
         function(value) is.null(value$mitochondrial_length),
         logical(1)
       ))){
      stop('A mitochondrial model collection must contain named prepared models.')
    }
    if(!('cell_type' %in% names(nodes))){
      stop('Lineage nodes need a cell_type column when using model collections.')
    }
    missing_models <- setdiff(unique(as.character(nodes$cell_type)), names(model))
    if(length(missing_models) > 0){
      stop(sprintf(
        'No mitochondrial model is available for cell type(s): %s.',
        paste(missing_models, collapse = ', ')
      ))
    }
    output_model <- model[[1]]
  } else{
    output_model <- model
  }
  # Internal: Select the prepared model that applies to one node
  #
  # @param node One-row slice of the node table.
  # @return The single prepared model, or the collection entry named by the
  #   node's `cell_type`.
  node_model <- function(node){
    if(model_is_collection){
      model[[as.character(node$cell_type)]]
    } else{
      model
    }
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
    branch_model <- node_model(node)
    if(is.na(node$parent_node_id)){
      initialized <- initialize_physicell_mito_profile(branch_model, params)
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
      if(!isTRUE(output_model$retain_internal_profiles) &&
         remaining_children[parent_node_id] <= 0){
        profiles[[parent_node_id]] <- NULL
      }
      start_profile <- if('origin' %in% names(nodes) &&
                          identical(
                            as.character(node$origin),
                            'induction_continuation'
                          )){
        # Induction changes a rate state but is not a cell division, so it must
        # not introduce an artificial mitochondrial bottleneck.
        parent_profile
      } else{
        inherit_physicell_mito_profile(
          parent_profile,
          branch_model$genomes_per_cell,
          branch_model$variant_selection_coefficient
        )
      }
    }

    branch_editing_state <- if('editing_state' %in% names(nodes)){
      as.character(node$editing_state)
    } else{
      editing_state
    }
    mutation <- mutate_physicell_mito_branch(
      start_profile,
      node$birth_time,
      node$end_time,
      branch_model,
      branch_editing_state
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

#' Summarize terminal-cell heteroplasmy as per-variant fractions
#'
#' For each sampled terminal cell, counts how many of its mt genomes carry each
#' distinct `(position, allele)` pair and divides by the cell's genome count.
#' Cells whose genomes are all reference contribute no rows.
#'
#' @export
#' @param simulation Result of `simulate_mito_on_physicell_lineage()`; its
#'   `terminal_nodes` and `profiles` entries are used.
#' @return A long-form data frame with columns `sample_id`
#'   (`cell_<physicell_id>`), `physicell_id`, `node_id`, `position`, `allele`,
#'   `variant_genomes` (carrier count), and `variant_fraction` (carriers over
#'   genomes per cell); a zero-row frame with the same columns when nothing
#'   mutated.
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

#' Reconstruct the nucleotide sequence of one mt genome
#'
#' Applies the sparse allele encoding on top of the reference: an integer
#' allele `1`-`4` replaces the base with `A`, `G`, `C`, `T`; `-1` removes the
#' base; and a fractional allele `0.b` appends base code `b` immediately after
#' the reference base at that position.
#'
#' @export
#' @param genome One sparse genome: a named numeric vector whose names are
#'   positions and whose values are alleles. An empty genome returns the
#'   reference unchanged.
#' @param reference Character vector of reference bases, one per position.
#' @return One collapsed sequence string. Deletions drop bases, so the result
#'   is not necessarily the same length as `reference`.
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

#' Write the mitochondrial replay outputs to disk
#'
#' Summarizes terminal-cell heteroplasmy, then writes the profiles, the event
#' table, the variant fractions, the reference sequence, and a run manifest.
#'
#' @export
#' @param simulation Result of `simulate_mito_on_physicell_lineage()`.
#' @param model Prepared model matching that replay; supplies the reference
#'   sequence and the manifest's length, genome-count, coefficient, and
#'   cell-type fields.
#' @param output_dir Destination directory, created recursively if absent.
#' @param write_fasta Whether to also emit one randomly sampled haplotype per
#'   terminal cell.
#' @param show_progress Whether to emit stage log messages.
#' @param compress_csv Whether the CSV tables are gzip-compressed, which adds a
#'   `.gz` suffix to their names.
#' @return Invisibly, a list with `terminal_profiles` (keyed
#'   `cell_<physicell_id>`), `variant_fractions`, and `sparse_scores`, which is
#'   `NULL` when the Matrix package is not installed.
#' @details The sparse score matrix has rows named `cell_<physicell_id>` and
#'   columns named `mt_<position>_<allele>`, with variant fractions as entries.
#' @section Side effects: Creates `output_dir` and writes
#'   `mitochondrial_profiles.rds`, `mitochondrial_mutation_events.csv`,
#'   `mitochondrial_variant_fractions.csv`, `mitochondrial_reference.fasta`
#'   (never compressed), and `mitochondrial_manifest.csv`; adds
#'   `mitochondrial_variant_fraction_matrix.rds` when Matrix is available and
#'   `mitochondrial_sampled_haplotypes.fasta` when `write_fasta` is `TRUE`.
#'   Choosing which haplotype to write consumes the random-number stream.
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
      'variant_selection_coefficient',
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
      model$variant_selection_coefficient,
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
