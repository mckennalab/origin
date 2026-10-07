# Non-Mendelian ecDNA inheritance with an integrated static identifier and a
# compact irreversible CRISPR recorder. Profiles are grouped by
# species/label/haplotype so copy-number-rich simulations do not allocate one R
# object per ecDNA molecule.

#' Read one ecDNA configuration key with a fallback
#'
#' @export
#' @param configuration Parsed `ecdna_adapter` list from the remote_mito
#'   parameter JSON; may be an empty list.
#' @param name Configuration key to look up.
#' @param default Value returned when the key is absent (`NULL`).
#' @return The configured value when present, otherwise `default`.
ecdna_config_value <- function(configuration, name, default){
  value <- configuration[[name]]
  if(is.null(value)) default else value
}

#' Validate one probability-valued ecDNA configuration entry
#'
#' Coerces with `as.numeric()` and requires a single finite scalar in the
#' closed interval `[0, 1]`; anything else raises an error naming the key.
#'
#' @export
#' @param value Candidate probability, in any form coercible to numeric.
#' @param name Diagnostic name used in the error message.
#' @return The validated numeric scalar.
ecdna_probability <- function(value, name){
  value <- suppressWarnings(as.numeric(value))
  if(length(value) != 1L || !is.finite(value) || value < 0 || value > 1){
    stop(sprintf('%s must be one probability in [0, 1].', name))
  }
  value
}

#' Validate one integer-valued ecDNA configuration entry
#'
#' @export
#' @param value Candidate count, coerced with `as.integer()`; must be a single
#'   non-`NA` value at or above the lower bound.
#' @param name Diagnostic name used in the error message.
#' @param allow_zero When `TRUE` the lower bound is `0`, otherwise `1`.
#' @return The validated integer scalar.
ecdna_positive_integer <- function(value, name, allow_zero = FALSE){
  value <- suppressWarnings(as.integer(value))
  lower_bound <- if(allow_zero) 0L else 1L
  if(length(value) != 1L || is.na(value) || value < lower_bound){
    stop(sprintf(
      '%s must be one %s integer.',
      name,
      if(allow_zero) 'non-negative' else 'positive'
    ))
  }
  value
}

#' Resolve the ecDNA species, recorder, and segregation model
#'
#' Reads the `ecdna_adapter` block of a parsed remote_mito parameter list and
#' turns it into the immutable model object every other function in this file
#' consumes. Founder species are named `ecDNA_<zero-padded index>`. The nearest
#' whole number to `labeled_species_fraction * num_species` species (or exactly
#' `num_labeled_species` when that key is present) is drawn uniformly at random
#' to be labeled; each labeled species gets an immutable static ID, a random
#' `A`/`C`/`G`/`T` string of `static_id_length` bases redrawn until all IDs are
#' distinct, plus an irreversible CRISPR recorder array of
#' `num_recorder_targets` binary sites. Unlabeled species carry copy number
#' only and never acquire recorder state.
#'
#' @export
#' @details
#' Units of every configured rate:
#' \itemize{
#'   \item `edit_probability_per_target_per_cell_cycle` is dimensionless, in
#'     `[0, 1)`, per recorder target per cell cycle. It is converted per cell
#'     type into the continuous-time hazard `-log(1 - p) / cell_cycle_length`,
#'     that is expected edits per unit simulation time per unedited copy per
#'     target, so that exactly `p` of the copies edit over one cycle in
#'     expectation. A cell type with infinite `cell_cycle_length` gets hazard
#'     `0` and never records.
#'   \item `recorder_start_time` is an absolute time in the same units as the
#'     lineage table's `birth_time` and `end_time`; edits before it are
#'     suppressed.
#'   \item `replication_probability` and `daughter_1_segregation_probability`
#'     are dimensionless per-copy probabilities applied once per cell division,
#'     not per unit time.
#'   \item `label_selection_coefficient` and
#'     `recorder_edit_selection_coefficient` are dimensionless coefficients in
#'     `[0, 1]` read from the `non_mendelian_selection` block. They multiply the
#'     per-copy replication probability at division only and never feed back
#'     into cellular division or death rates.
#'   \item `max_copies_per_cell` is a copy count per cell, enforced after each
#'     division.
#' }
#' Validation rejects more than 64 recorder targets, an
#' `initial_copies_per_species` vector that is neither length one nor one
#' positive value per species, a labeled-species count above `num_species`, a
#' `static_id_length` too short to give the labeled species distinct IDs, and
#' any cell type without a positive `cell_cycle_length`.
#'
#' @param params Parsed remote_mito parameter list. `params$ecdna_adapter`
#'   supplies the configuration (documented defaults are used when it is
#'   absent), `params$cell_type_dict$cell_type_params` supplies one
#'   `cell_cycle_length` per cell type, and `params$random_seed` is the fallback
#'   seed.
#' @param seed Optional non-negative integer seed; `params$random_seed` is used
#'   when this is `NULL`.
#' @return A named list holding `num_species`, `species_ids`,
#'   `initial_copies_per_species`, `labeled_species_indices`, the requested and
#'   `realized_labeled_species_fraction`, `static_ids` (`NA` for unlabeled
#'   species), `static_id_length`, `num_recorder_targets`, `target_names`,
#'   `edit_probability_per_target_per_cell_cycle`, `edit_hazards` (a per-target
#'   named numeric vector for each cell type, keyed by cell-type name),
#'   `recorder_start_time`, `replication_probability`, both selection
#'   coefficients, `daughter_1_segregation_probability`, `max_copies_per_cell`,
#'   `retain_internal_profiles`, and the resolved `seed`.
#' @section Side effects:
#' Calls `set.seed(seed)`, resetting the global RNG stream, before sampling the
#' labeled species and their static IDs.
#' @note Depends on `non_mendelian_selection_coefficient()` from
#'   `physicell_lineage.R`, which must be sourced first.
prepare_ecdna_model <- function(params, seed = NULL){
  if(!is.list(params)){
    stop('params must be a parsed remote_mito parameter list.')
  }
  configuration <- params$ecdna_adapter
  if(is.null(configuration)){
    configuration <- list()
  }
  if(is.null(seed)){
    seed <- params$random_seed
  }
  seed <- ecdna_positive_integer(seed, 'ecDNA seed', allow_zero = TRUE)
  set.seed(seed)

  num_species <- ecdna_positive_integer(
    ecdna_config_value(configuration, 'num_species', 10L),
    'ecdna_adapter.num_species'
  )
  num_targets <- ecdna_positive_integer(
    ecdna_config_value(configuration, 'num_recorder_targets', 6L),
    'ecdna_adapter.num_recorder_targets'
  )
  if(num_targets > 64L){
    stop('ecdna_adapter.num_recorder_targets cannot exceed 64.')
  }
  copies <- as.integer(unlist(
    ecdna_config_value(configuration, 'initial_copies_per_species', 5L),
    use.names = FALSE
  ))
  if(length(copies) == 1L){
    copies <- rep(copies, num_species)
  }
  if(length(copies) != num_species || anyNA(copies) || any(copies < 1L)){
    stop(paste(
      'ecdna_adapter.initial_copies_per_species must be one positive',
      'integer or one value per species.'
    ))
  }
  labeled_fraction <- ecdna_probability(
    ecdna_config_value(configuration, 'labeled_species_fraction', 0.25),
    'ecdna_adapter.labeled_species_fraction'
  )
  configured_labeled_count <- configuration$num_labeled_species
  num_labeled <- if(is.null(configured_labeled_count)){
    min(num_species, max(
      0L,
      as.integer(floor(labeled_fraction * num_species + 0.5))
    ))
  } else{
    count <- ecdna_positive_integer(
      configured_labeled_count,
      'ecdna_adapter.num_labeled_species',
      allow_zero = TRUE
    )
    if(count > num_species){
      stop('ecdna_adapter.num_labeled_species cannot exceed num_species.')
    }
    count
  }
  labeled_species <- if(num_labeled == 0L){
    integer()
  } else{
    sort(sample.int(num_species, num_labeled, replace = FALSE))
  }

  edit_probabilities <- as.numeric(unlist(
    ecdna_config_value(
      configuration,
      'edit_probability_per_target_per_cell_cycle',
      0.01
    ),
    use.names = FALSE
  ))
  if(length(edit_probabilities) == 1L){
    edit_probabilities <- rep(edit_probabilities, num_targets)
  }
  if(length(edit_probabilities) != num_targets ||
     any(!is.finite(edit_probabilities)) ||
     any(edit_probabilities < 0) || any(edit_probabilities >= 1)){
    stop(paste(
      'ecdna_adapter.edit_probability_per_target_per_cell_cycle must',
      'contain probabilities in [0, 1), one value or one per target.'
    ))
  }
  recorder_start_time <- suppressWarnings(as.numeric(ecdna_config_value(
    configuration,
    'recorder_start_time',
    0
  )))
  if(length(recorder_start_time) != 1L ||
     !is.finite(recorder_start_time) || recorder_start_time < 0){
    stop('ecdna_adapter.recorder_start_time must be finite and non-negative.')
  }
  replication_probability <- ecdna_probability(
    ecdna_config_value(configuration, 'replication_probability', 1),
    'ecdna_adapter.replication_probability'
  )
  segregation_probability <- ecdna_probability(
    ecdna_config_value(configuration, 'daughter_1_segregation_probability', 0.5),
    'ecdna_adapter.daughter_1_segregation_probability'
  )
  static_id_length <- ecdna_positive_integer(
    ecdna_config_value(configuration, 'static_id_length', 12L),
    'ecdna_adapter.static_id_length'
  )
  max_copies <- ecdna_positive_integer(
    ecdna_config_value(configuration, 'max_copies_per_cell', 10000L),
    'ecdna_adapter.max_copies_per_cell'
  )
  label_selection_coefficient <- non_mendelian_selection_coefficient(
    params,
    'ecdna_label_coefficient'
  )
  recorder_edit_selection_coefficient <- non_mendelian_selection_coefficient(
    params,
    'ecdna_recorder_edit_coefficient'
  )
  species_width <- max(3L, nchar(num_species))
  species_ids <- paste0(
    'ecDNA_',
    formatC(seq_len(num_species), width = species_width, flag = '0')
  )
  static_ids <- rep(NA_character_, num_species)
  if(num_labeled > 0L){
    if(log(num_labeled) > static_id_length * log(4) + 1e-12){
      stop('static_id_length is too short to uniquely label the selected species.')
    }
    repeat{
      candidate_ids <- vapply(
        seq_len(num_labeled),
        function(index){
          paste0(
            sample(c('A', 'C', 'G', 'T'), static_id_length, TRUE),
            collapse = ''
          )
        },
        character(1)
      )
      if(!anyDuplicated(candidate_ids)){
        break
      }
    }
    static_ids[labeled_species] <- candidate_ids
  }
  cell_cycle_lengths <- vapply(
    params$cell_type_dict$cell_type_params,
    function(type_params){
      value <- as.numeric(type_params$cell_cycle_length)
      if(length(value) != 1L || is.na(value) || value <= 0){
        stop('Every ecDNA cell type needs a positive cell_cycle_length.')
      }
      value
    },
    numeric(1)
  )
  target_names <- paste0('target_', seq_len(num_targets))
  edit_hazards <- lapply(cell_cycle_lengths, function(cycle_length){
    if(is.infinite(cycle_length)){
      setNames(rep(0, num_targets), target_names)
    } else{
      setNames(-log1p(-edit_probabilities) / cycle_length, target_names)
    }
  })

  list(
    num_species = num_species,
    species_ids = species_ids,
    initial_copies_per_species = copies,
    labeled_species_indices = labeled_species,
    labeled_species_fraction = labeled_fraction,
    realized_labeled_species_fraction = num_labeled / num_species,
    static_ids = static_ids,
    static_id_length = static_id_length,
    num_recorder_targets = num_targets,
    target_names = target_names,
    edit_probability_per_target_per_cell_cycle = edit_probabilities,
    edit_hazards = edit_hazards,
    recorder_start_time = recorder_start_time,
    replication_probability = replication_probability,
    label_selection_coefficient = label_selection_coefficient,
    recorder_edit_selection_coefficient =
      recorder_edit_selection_coefficient,
    daughter_1_segregation_probability = segregation_probability,
    max_copies_per_cell = max_copies,
    retain_internal_profiles = isTRUE(ecdna_config_value(
      configuration,
      'retain_internal_profiles',
      FALSE
    )),
    seed = seed
  )
}

#' Build the founder ecDNA profile of a root cell
#'
#' Emits one grouped row per species carrying that species'
#' `initial_copies_per_species` copies. Labeled species start with an all-`0`
#' haplotype string of length `num_recorder_targets`; unlabeled species carry
#' the empty string.
#'
#' @export
#' @param model Prepared model from `prepare_ecdna_model()`.
#' @return A data frame with `species_index` (integer), `labeled` (logical),
#'   `haplotype` (character, `0` unedited and `1` edited per target), and
#'   `count` (integer copies in that group).
initialize_ecdna_profile <- function(model){
  labeled <- seq_len(model$num_species) %in% model$labeled_species_indices
  data.frame(
    species_index = seq_len(model$num_species),
    labeled = labeled,
    haplotype = ifelse(
      labeled,
      strrep('0', model$num_recorder_targets),
      ''
    ),
    count = as.integer(model$initial_copies_per_species),
    stringsAsFactors = FALSE
  )
}

#' Merge duplicate ecDNA groups and drop empty ones
#'
#' Removes rows with no remaining copies, sums the counts of rows sharing the
#' same species/labeled/haplotype key, orders the result by species index,
#' label flag, and haplotype, and resets row names. This keeps the profile one
#' row per distinct molecule state rather than one row per molecule.
#'
#' @export
#' @param profile Grouped ecDNA profile.
#' @return The collapsed profile with the same columns; an empty or all-zero
#'   input comes back as a zero-row data frame.
collapse_ecdna_profile <- function(profile){
  if(nrow(profile) == 0L){
    return(profile)
  }
  profile <- profile[profile$count > 0L, , drop = FALSE]
  if(nrow(profile) == 0L){
    return(profile)
  }
  key <- paste(
    profile$species_index,
    as.integer(profile$labeled),
    profile$haplotype,
    sep = '|'
  )
  groups <- split(seq_len(nrow(profile)), key)
  collapsed <- lapply(groups, function(indices){
    row <- profile[indices[1], , drop = FALSE]
    row$count <- sum(profile$count[indices])
    row
  })
  result <- do.call(rbind, collapsed)
  rownames(result) <- NULL
  result[order(result$species_index, result$labeled, result$haplotype), , drop = FALSE]
}

#' Draw exact edit times conditional on editing inside a segment
#'
#' Inverse-CDF sample of an exponential waiting time with rate `hazard`,
#' truncated to `(start_time, start_time + duration]`, i.e. conditioned on the
#' edit having occurred within the segment. Draws are independent across copies.
#'
#' @export
#' @param number Number of edited copies to draw times for; `0` yields
#'   `numeric()`.
#' @param hazard Edit hazard in events per unit simulation time; a non-positive
#'   hazard yields `numeric()`.
#' @param duration Segment length in simulation time units; a non-positive
#'   duration yields `numeric()`.
#' @param start_time Absolute simulation time at which the segment begins.
#' @return A numeric vector of `number` absolute event times inside the segment.
ecdna_truncated_event_times <- function(number, hazard, duration, start_time){
  if(number == 0L || hazard <= 0 || duration <= 0){
    return(numeric())
  }
  event_probability <- -expm1(-hazard * duration)
  start_time - log1p(-stats::runif(number) * event_probability) / hazard
}

#' Apply irreversible recorder edits over one constant-hazard interval
#'
#' Walks the recorder targets in index order. For a given target, every copy in
#' every labeled group that is still `0` at that target edits independently with
#' probability `1 - exp(-hazard * duration)`, so the number of edited copies per
#' group is a binomial draw; those copies are moved out of their group into a
#' new group whose haplotype has that target flipped to `1`, and the profile is
#' collapsed before the next target. Edits are irreversible and unlabeled
#' species (empty haplotype) are never touched. Because the per-target
#' probability does not depend on the rest of the haplotype, the result is
#' equivalent to each copy editing each target independently within the
#' interval, and a copy can gain several edits in one interval.
#'
#' @export
#' @param profile Grouped ecDNA profile inherited at the start of the interval.
#' @param duration Interval length in simulation time units; a non-positive
#'   duration, an empty profile, or a model with no labeled species returns the
#'   profile unchanged.
#' @param segment_start Absolute simulation time at which the interval begins;
#'   used to place the sampled event times.
#' @param hazards Numeric vector of per-target edit hazards, in events per unit
#'   simulation time per unedited copy, for the branch's cell type; targets with
#'   hazard `0` are skipped.
#' @param model Prepared model, used for `num_recorder_targets` and
#'   `labeled_species_indices`.
#' @return A named list with `profile` (the collapsed grouped profile) and
#'   `events`: `NULL` when nothing edited, otherwise one aggregate row per
#'   source group and target holding `species_index`, `target_index`,
#'   `edited_copies`, the `first_event_time`, `mean_event_time`, and
#'   `last_event_time` of the copies edited in that draw, and the segment
#'   bounds. Individual copies are not tracked, so event rows are per group,
#'   not per molecule.
mutate_ecdna_profile_segment <- function(profile,
                                         duration,
                                         segment_start,
                                         hazards,
                                         model){
  if(duration <= 0 || nrow(profile) == 0L ||
     length(model$labeled_species_indices) == 0L){
    return(list(profile = profile, events = NULL))
  }
  event_rows <- list()
  event_count <- 0L
  for(target_index in seq_len(model$num_recorder_targets)){
    hazard <- unname(hazards[target_index])
    if(hazard <= 0){
      next
    }
    target_probability <- -expm1(-hazard * duration)
    initial_rows <- nrow(profile)
    for(row_index in seq_len(initial_rows)){
      if(!profile$labeled[row_index] || profile$count[row_index] == 0L){
        next
      }
      states <- strsplit(profile$haplotype[row_index], '', fixed = TRUE)[[1]]
      if(states[target_index] != '0'){
        next
      }
      edited_copies <- stats::rbinom(
        1L,
        size = profile$count[row_index],
        prob = target_probability
      )
      if(edited_copies == 0L){
        next
      }
      profile$count[row_index] <- profile$count[row_index] - edited_copies
      states[target_index] <- '1'
      new_haplotype <- paste0(states, collapse = '')
      profile <- rbind(
        profile,
        data.frame(
          species_index = profile$species_index[row_index],
          labeled = TRUE,
          haplotype = new_haplotype,
          count = edited_copies,
          stringsAsFactors = FALSE
        )
      )
      event_times <- ecdna_truncated_event_times(
        edited_copies,
        hazard,
        duration,
        segment_start
      )
      event_count <- event_count + 1L
      event_rows[[event_count]] <- data.frame(
        species_index = profile$species_index[row_index],
        target_index = target_index,
        edited_copies = edited_copies,
        first_event_time = min(event_times),
        mean_event_time = mean(event_times),
        last_event_time = max(event_times),
        segment_start = segment_start,
        segment_end = segment_start + duration,
        stringsAsFactors = FALSE
      )
    }
    profile <- collapse_ecdna_profile(profile)
  }
  list(
    profile = profile,
    events = if(length(event_rows) == 0L) NULL else do.call(rbind, event_rows)
  )
}

#' Edit a profile over the recorder-active part of one lineage branch
#'
#' Clips the branch to start no earlier than `model$recorder_start_time` and
#' applies the branch cell type's hazards to whatever span remains. A branch
#' that ends before the recorder switches on is returned unchanged.
#'
#' @export
#' @param profile Grouped ecDNA profile inherited at `start_time`.
#' @param start_time Branch birth time in simulation time units.
#' @param end_time Branch end time in the same units.
#' @param cell_type Cell-type key; it must name an entry of
#'   `model$edit_hazards` or an error is raised.
#' @param model Prepared model from `prepare_ecdna_model()`.
#' @return The same `profile`/`events` list shape as
#'   `mutate_ecdna_profile_segment()`.
mutate_ecdna_profile_branch <- function(profile,
                                        start_time,
                                        end_time,
                                        cell_type,
                                        model){
  mutation_start <- max(start_time, model$recorder_start_time)
  if(end_time <= mutation_start){
    return(list(profile = profile, events = NULL))
  }
  hazards <- model$edit_hazards[[as.character(cell_type)]]
  if(is.null(hazards)){
    stop(sprintf('No ecDNA edit hazards are configured for cell type %s.', cell_type))
  }
  mutate_ecdna_profile_segment(
    profile,
    duration = end_time - mutation_start,
    segment_start = mutation_start,
    hazards = hazards,
    model = model
  )
}

#' Enforce the per-cell ecDNA copy ceiling
#'
#' When a cell's total copy number exceeds `maximum`, the whole pool is
#' resampled as a single multinomial draw of exactly `maximum` copies with group
#' probabilities proportional to the current counts. Composition is therefore
#' preserved in expectation, but a rare group can be lost outright.
#'
#' @export
#' @param profile Grouped ecDNA profile.
#' @param maximum Copy ceiling per cell, in copies.
#' @return The profile unchanged when the total is at or below the cap,
#'   otherwise the downsampled and collapsed profile.
cap_ecdna_profile <- function(profile, maximum){
  total <- sum(profile$count)
  if(total <= maximum){
    return(profile)
  }
  profile$count <- as.integer(stats::rmultinom(
    1L,
    size = maximum,
    prob = profile$count
  )[, 1])
  collapse_ecdna_profile(profile)
}

#' Replicate and randomly segregate ecDNA between two daughters
#'
#' The non-Mendelian division step, applied once per binary division. Each copy
#' present in the parent first retains itself and independently produces one
#' additional copy with probability `p_rep`, so a group of `n` copies becomes
#' `n + Binomial(n, p_rep)` copies and can at most double. The replicated pool
#' of each group is then split by an independent
#' `Binomial(replicated, daughter_1_segregation_probability)` draw assigned to
#' daughter 1, with daughter 2 receiving the exact complement. Daughters
#' therefore get complementary, generally unequal copy numbers, total copy
#' number is conserved across the pair, and either daughter can lose a species
#' entirely. Each daughter is finally collapsed and capped at
#' `model$max_copies_per_cell`.
#'
#' @export
#' @details
#' Selection enters only through `p_rep`, the per-copy replication probability.
#' When both coefficients are zero, `p_rep` is the scalar
#' `model$replication_probability` for every group. Otherwise it is computed per
#' group as
#' `min(1, replication_probability * (1 - s_label)^labeled * (1 - s_edit)^b)`,
#' where `s_label` is `label_selection_coefficient`, `s_edit` is
#' `recorder_edit_selection_coefficient`, `labeled` is `1` for a labeled group
#' and `0` otherwise, and `b` is that haplotype's edit burden, its count of
#' non-`0` characters. The label penalty is applied once to every labeled copy
#' regardless of burden; the edit penalty compounds once per edited target.
#' Segregation itself is never selection-weighted, and neither coefficient
#' alters the cell's division or death rate: the fitness coupling is
#' intracellular copy replication only.
#'
#' @param profile Parent grouped ecDNA profile; a zero-row profile is handed to
#'   both daughters unchanged.
#' @param model Prepared model supplying `replication_probability`, both
#'   selection coefficients, `daughter_1_segregation_probability`, and
#'   `max_copies_per_cell`.
#' @return An unnamed length-two list of grouped daughter profiles.
partition_ecdna_profile <- function(profile, model){
  if(nrow(profile) == 0L){
    return(list(profile, profile))
  }
  effective_replication_probability <- if(
    model$label_selection_coefficient == 0 &&
    model$recorder_edit_selection_coefficient == 0
  ){
    model$replication_probability
  } else{
    edit_burden <- nchar(gsub('0', '', profile$haplotype, fixed = TRUE))
    relative_replication <- ifelse(
      profile$labeled,
      1 - model$label_selection_coefficient,
      1
    ) * (1 - model$recorder_edit_selection_coefficient)^edit_burden
    pmin(1, model$replication_probability * relative_replication)
  }
  replicated <- profile$count + stats::rbinom(
    nrow(profile),
    size = profile$count,
    prob = effective_replication_probability
  )
  daughter_1_counts <- stats::rbinom(
    nrow(profile),
    size = replicated,
    prob = model$daughter_1_segregation_probability
  )
  daughter_1 <- profile
  daughter_2 <- profile
  daughter_1$count <- daughter_1_counts
  daughter_2$count <- replicated - daughter_1_counts
  list(
    cap_ecdna_profile(
      collapse_ecdna_profile(daughter_1),
      model$max_copies_per_cell
    ),
    cap_ecdna_profile(
      collapse_ecdna_profile(daughter_2),
      model$max_copies_per_cell
    )
  )
}

#' Propagate ecDNA profiles down an event-resolved lineage
#'
#' Sweeps the node table once in parent-before-child order. A root node (`NA`
#' parent) starts from the founder profile; every other node inherits the
#' profile its parent handed it. Each node then edits its profile over the
#' recorder-active part of its own branch and passes it on according to its
#' child count: a node with one child, such as an induction-continuation node,
#' hands the profile through unchanged with no replication or segregation; a
#' node with two children splits it via `partition_ecdna_profile()`; a node with
#' no children keeps it as terminal output. More than two children is an error,
#' as is a table that is not sorted parent-before-child.
#'
#' @export
#' @param nodes Event-resolved lineage data frame. It must contain `node_id`,
#'   `physicell_id`, `parent_node_id`, `birth_time`, `end_time`, and
#'   `is_terminal`, list every parent before its children, and use `NA` parents
#'   for roots. An optional `cell_type` column selects the per-branch hazard
#'   set; without it the first configured cell type is used on every branch.
#' @param model Prepared model from `prepare_ecdna_model()`.
#' @param terminal_physicell_ids Optional sampled-cell filter; terminal nodes
#'   whose `physicell_id` is outside this set are dropped from the output and at
#'   least one must remain.
#' @param seed Integer seed for the recording RNG stream.
#' @param show_progress Whether to report sweep progress to the console.
#' @param progress_updates Number of progress checkpoints to report.
#' @return A named list with `nodes` (the input table), `terminal_nodes` (the
#'   sampled terminal subset), `profiles` (grouped profiles named by `node_id`,
#'   leaves only unless `model$retain_internal_profiles` is `TRUE`),
#'   `mutation_events` (the concatenated aggregate edit rows, extended with
#'   `species_id`, `static_id`, `node_id`, `physicell_id`, `parent_node_id`, and
#'   `cell_type`, or an empty typed frame when nothing edited), and `seed`.
#' @section Side effects:
#' Calls `set.seed(seed)`, resetting the global RNG stream, and prints progress
#' lines when `show_progress` is `TRUE`.
#' @note An inherited profile is released as soon as its node consumes it, so
#'   internal nodes are absent from `profiles` unless
#'   `model$retain_internal_profiles` is set. Depends on
#'   `new_physicell_progress_reporter()` from `physicell_lineage.R`.
simulate_ecdna_on_lineage <- function(nodes,
                                      model,
                                      terminal_physicell_ids = NULL,
                                      seed = 1L,
                                      show_progress = FALSE,
                                      progress_updates = 20L){
  required <- c(
    'node_id', 'physicell_id', 'parent_node_id', 'birth_time', 'end_time',
    'is_terminal'
  )
  if(!is.data.frame(nodes) || !all(required %in% names(nodes))){
    stop('nodes is not a valid event-resolved lineage table.')
  }
  set.seed(seed)
  num_nodes <- nrow(nodes)
  parent_indices <- match(nodes$parent_node_id, nodes$node_id)
  non_roots <- which(!is.na(nodes$parent_node_id))
  if(anyNA(parent_indices[non_roots]) ||
     any(parent_indices[non_roots] >= non_roots)){
    stop('ecDNA lineage nodes must be in parent-before-child order.')
  }
  children <- split(non_roots, parent_indices[non_roots])
  start_profiles <- vector('list', num_nodes)
  profiles <- vector('list', num_nodes)
  names(profiles) <- nodes$node_id
  all_events <- vector('list', num_nodes)
  event_count <- 0L
  report_progress <- new_physicell_progress_reporter(
    total = num_nodes,
    label = 'ecDNA lineage recording',
    enabled = show_progress,
    updates = progress_updates
  )
  report_progress(0L)

  for(node_index in seq_len(num_nodes)){
    node <- nodes[node_index, , drop = FALSE]
    profile <- if(is.na(node$parent_node_id)){
      initialize_ecdna_profile(model)
    } else{
      start_profiles[[node_index]]
    }
    if(is.null(profile)){
      stop(sprintf('No inherited ecDNA profile was available for %s.', node$node_id))
    }
    start_profiles[node_index] <- list(NULL)
    cell_type <- if('cell_type' %in% names(nodes)){
      as.character(node$cell_type)
    } else{
      names(model$edit_hazards)[1]
    }
    mutation <- mutate_ecdna_profile_branch(
      profile,
      node$birth_time,
      node$end_time,
      cell_type,
      model
    )
    profile <- mutation$profile
    if(!is.null(mutation$events)){
      mutation$events$species_id <- model$species_ids[
        mutation$events$species_index
      ]
      mutation$events$static_id <- model$static_ids[
        mutation$events$species_index
      ]
      mutation$events$node_id <- node$node_id
      mutation$events$physicell_id <- node$physicell_id
      mutation$events$parent_node_id <- node$parent_node_id
      mutation$events$cell_type <- cell_type
      event_count <- event_count + 1L
      all_events[[event_count]] <- mutation$events
    }

    child_indices <- children[[as.character(node_index)]]
    if(is.null(child_indices)){
      profiles[[node_index]] <- profile
    } else if(length(child_indices) == 1L){
      start_profiles[[child_indices]] <- profile
      if(isTRUE(model$retain_internal_profiles)){
        profiles[[node_index]] <- profile
      }
    } else if(length(child_indices) == 2L){
      daughters <- partition_ecdna_profile(profile, model)
      start_profiles[[child_indices[1]]] <- daughters[[1]]
      start_profiles[[child_indices[2]]] <- daughters[[2]]
      if(isTRUE(model$retain_internal_profiles)){
        profiles[[node_index]] <- profile
      }
    } else{
      stop(sprintf(
        'Node %s has %d children; ecDNA segregation requires at most two.',
        node$node_id,
        length(child_indices)
      ))
    }
    report_progress(node_index)
  }

  terminal_nodes <- nodes[nodes$is_terminal, , drop = FALSE]
  if(!is.null(terminal_physicell_ids)){
    terminal_nodes <- terminal_nodes[
      terminal_nodes$physicell_id %in% as.character(terminal_physicell_ids),
      ,
      drop = FALSE
    ]
  }
  if(nrow(terminal_nodes) == 0L){
    stop('No sampled terminal cells remain for ecDNA output.')
  }
  events <- if(event_count == 0L){
    data.frame(
      species_index = integer(), target_index = integer(),
      edited_copies = integer(), first_event_time = numeric(),
      mean_event_time = numeric(), last_event_time = numeric(),
      segment_start = numeric(), segment_end = numeric(),
      species_id = character(), static_id = character(), node_id = character(),
      physicell_id = character(), parent_node_id = character(),
      cell_type = character(), stringsAsFactors = FALSE
    )
  } else{
    do.call(rbind, all_events[seq_len(event_count)])
  }
  profiles <- profiles[!vapply(profiles, is.null, logical(1))]
  list(
    nodes = nodes,
    terminal_nodes = terminal_nodes,
    profiles = profiles,
    mutation_events = events,
    seed = seed
  )
}

#' Summarize terminal ecDNA profiles into per-cell and long haplotype tables
#'
#' Walks the sampled terminal nodes and reduces each retained profile to one
#' burden summary row plus its grouped haplotype rows. Sample identifiers are
#' `cell_<physicell_id>`.
#'
#' @export
#' @param simulation Completed simulation from `simulate_ecdna_on_lineage()`.
#' @param model Prepared model, used for `species_ids`, `static_ids`, and
#'   `num_recorder_targets`.
#' @return A named list with `cell_summary` (one row per sampled terminal cell:
#'   identifiers, `total_ecdna_copies`, `labeled_ecdna_copies`, the number of
#'   distinct retained species overall and among labeled species, and
#'   `ecdna_recorder_edit_fraction`, the copy-weighted count of edited target
#'   slots divided by `labeled_ecdna_copies * num_recorder_targets`, reported as
#'   `0` when no labeled copies survive) and `haplotypes` (one row per
#'   cell/species/haplotype group with `species_id`, `static_id`, `labeled`,
#'   `haplotype`, `edit_count`, and `count`; an empty typed frame when no cell
#'   retains any copies).
ecdna_terminal_tables <- function(simulation, model){
  summary_rows <- vector('list', nrow(simulation$terminal_nodes))
  haplotype_rows <- vector('list', nrow(simulation$terminal_nodes))
  haplotype_count <- 0L
  for(cell_index in seq_len(nrow(simulation$terminal_nodes))){
    node <- simulation$terminal_nodes[cell_index, , drop = FALSE]
    profile <- simulation$profiles[[as.character(node$node_id)]]
    labeled_rows <- profile$labeled
    total_copies <- sum(profile$count)
    labeled_copies <- sum(profile$count[labeled_rows])
    edit_count <- 0
    if(any(labeled_rows)){
      edit_count <- sum(vapply(
        which(labeled_rows),
        function(row_index){
          states <- strsplit(profile$haplotype[row_index], '', fixed = TRUE)[[1]]
          sum(states == '1') * profile$count[row_index]
        },
        numeric(1)
      ))
    }
    summary_rows[[cell_index]] <- data.frame(
      sample_id = paste0('cell_', node$physicell_id),
      physicell_id = node$physicell_id,
      node_id = node$node_id,
      total_ecdna_copies = total_copies,
      labeled_ecdna_copies = labeled_copies,
      observed_ecdna_species = length(unique(profile$species_index)),
      observed_labeled_species = length(unique(profile$species_index[labeled_rows])),
      ecdna_recorder_edit_fraction = if(labeled_copies > 0){
        edit_count / (labeled_copies * model$num_recorder_targets)
      } else{
        0
      },
      stringsAsFactors = FALSE
    )
    if(nrow(profile) > 0L){
      profile$sample_id <- paste0('cell_', node$physicell_id)
      profile$physicell_id <- node$physicell_id
      profile$node_id <- node$node_id
      profile$species_id <- model$species_ids[profile$species_index]
      profile$static_id <- model$static_ids[profile$species_index]
      profile$edit_count <- vapply(
        profile$haplotype,
        function(haplotype){
          if(!nzchar(haplotype)) 0L else sum(strsplit(haplotype, '', TRUE)[[1]] == '1')
        },
        integer(1)
      )
      haplotype_count <- haplotype_count + 1L
      haplotype_rows[[haplotype_count]] <- profile[, c(
        'sample_id', 'physicell_id', 'node_id', 'species_id', 'static_id',
        'labeled', 'haplotype', 'edit_count', 'count'
      )]
    }
  }
  list(
    cell_summary = do.call(rbind, summary_rows),
    haplotypes = if(haplotype_count == 0L){
      data.frame(
        sample_id = character(), physicell_id = character(),
        node_id = character(), species_id = character(), static_id = character(),
        labeled = logical(), haplotype = character(), edit_count = integer(),
        count = integer(), stringsAsFactors = FALSE
      )
    } else{
      do.call(rbind, haplotype_rows[seq_len(haplotype_count)])
    }
  )
}

#' Build sparse cell-by-feature ecDNA matrices from terminal inheritance
#'
#' Feature columns cover the labeled species only. `static_id_copy_numbers` has
#' one column per labeled species, named `<species_id>_<static_id>`, holding
#' that cell's copy number of the species. The two recorder matrices have one
#' column per species and target, named `<species_id>_<static_id>_<target>` and
#' ordered species-major so that species `k`'s target `t` is column
#' `(k - 1) * num_recorder_targets + t`; the fraction matrix stores the share of
#' that species' copies in the cell whose target is edited, and the presence
#' matrix is the same sparsity pattern with every stored value set to `1`.
#' Entries reflect the copies a cell actually inherited, not every topological
#' descendant of the branch on which an edit arose, because random segregation
#' breaks that Mendelian assumption.
#'
#' @export
#' @param simulation Completed simulation from `simulate_ecdna_on_lineage()`.
#' @param model Prepared model supplying `labeled_species_indices`,
#'   `species_ids`, `static_ids`, `target_names`, and `num_recorder_targets`.
#' @return A named list of three sparse matrices with rows named
#'   `cell_<physicell_id>`: `static_id_copy_numbers`,
#'   `recorder_edit_fractions`, and `recorder_edit_presence`.
#' @note Requires the `Matrix` package; its absence raises an error.
ecdna_feature_matrices <- function(simulation, model){
  if(!requireNamespace('Matrix', quietly = TRUE)){
    stop('The Matrix package is required for ecDNA feature matrices.')
  }
  sample_ids <- paste0('cell_', simulation$terminal_nodes$physicell_id)
  labeled_species <- model$labeled_species_indices
  static_feature_names <- paste0(
    model$species_ids[labeled_species],
    '_',
    model$static_ids[labeled_species]
  )
  target_feature_names <- as.vector(t(outer(
    static_feature_names,
    model$target_names,
    paste,
    sep = '_'
  )))
  num_cells <- nrow(simulation$terminal_nodes)
  static_i_chunks <- static_j_chunks <- static_x_chunks <- vector('list', num_cells)
  target_i_chunks <- target_j_chunks <- target_x_chunks <- vector('list', num_cells)
  for(cell_index in seq_len(nrow(simulation$terminal_nodes))){
    node <- simulation$terminal_nodes[cell_index, , drop = FALSE]
    profile <- simulation$profiles[[as.character(node$node_id)]]
    for(species_position in seq_along(labeled_species)){
      species_index <- labeled_species[species_position]
      species_rows <- which(profile$species_index == species_index & profile$labeled)
      if(length(species_rows) == 0L){
        next
      }
      total_copies <- sum(profile$count[species_rows])
      if(total_copies == 0L){
        next
      }
      static_i_chunks[[cell_index]] <- c(
        static_i_chunks[[cell_index]],
        cell_index
      )
      static_j_chunks[[cell_index]] <- c(
        static_j_chunks[[cell_index]],
        species_position
      )
      static_x_chunks[[cell_index]] <- c(
        static_x_chunks[[cell_index]],
        total_copies
      )
      edited_counts <- numeric(model$num_recorder_targets)
      for(row_index in species_rows){
        states <- strsplit(profile$haplotype[row_index], '', fixed = TRUE)[[1]]
        edited_counts <- edited_counts +
          (states == '1') * profile$count[row_index]
      }
      edited_targets <- which(edited_counts > 0)
      if(length(edited_targets) > 0L){
        target_i_chunks[[cell_index]] <- c(
          target_i_chunks[[cell_index]],
          rep(cell_index, length(edited_targets))
        )
        target_j_chunks[[cell_index]] <- c(
          target_j_chunks[[cell_index]],
          (species_position - 1L) * model$num_recorder_targets + edited_targets
        )
        target_x_chunks[[cell_index]] <- c(
          target_x_chunks[[cell_index]],
          edited_counts[edited_targets] / total_copies
        )
      }
    }
  }
  static_i <- as.integer(unlist(static_i_chunks, use.names = FALSE))
  static_j <- as.integer(unlist(static_j_chunks, use.names = FALSE))
  static_x <- as.numeric(unlist(static_x_chunks, use.names = FALSE))
  target_i <- as.integer(unlist(target_i_chunks, use.names = FALSE))
  target_j <- as.integer(unlist(target_j_chunks, use.names = FALSE))
  target_x <- as.numeric(unlist(target_x_chunks, use.names = FALSE))
  static_matrix <- Matrix::sparseMatrix(
    i = static_i,
    j = static_j,
    x = static_x,
    dims = c(length(sample_ids), length(static_feature_names)),
    dimnames = list(sample_ids, static_feature_names)
  )
  edit_fraction_matrix <- Matrix::sparseMatrix(
    i = target_i,
    j = target_j,
    x = target_x,
    dims = c(length(sample_ids), length(target_feature_names)),
    dimnames = list(sample_ids, target_feature_names)
  )
  edit_presence_matrix <- edit_fraction_matrix
  if(length(edit_presence_matrix@x) > 0L){
    edit_presence_matrix@x[] <- 1
  }
  list(
    static_id_copy_numbers = static_matrix,
    recorder_edit_fractions = edit_fraction_matrix,
    recorder_edit_presence = edit_presence_matrix
  )
}

#' Write every ecDNA artifact for one completed simulation
#'
#' Summarizes the terminal profiles, builds the feature matrices, and persists
#' both alongside the species manifest and a run manifest of the model settings.
#'
#' @export
#' @param simulation Completed simulation from `simulate_ecdna_on_lineage()`.
#' @param model Prepared model from `prepare_ecdna_model()`.
#' @param output_dir Destination directory, created recursively when missing.
#' @param show_progress Whether to emit stage log lines.
#' @param compress_csv Whether the tabular outputs are gzip-compressed.
#' @return Invisibly, the terminal tables and feature matrices concatenated into
#'   one list: `cell_summary`, `haplotypes`, `static_id_copy_numbers`,
#'   `recorder_edit_fractions`, and `recorder_edit_presence`.
#' @section Side effects:
#' Creates `output_dir` and writes `ecdna_cell_summary.csv`,
#' `ecdna_haplotypes.csv`, `ecdna_mutation_events.csv`,
#' `ecdna_species_manifest.csv`, and `ecdna_manifest.csv` through
#' `write_physicell_csv()` (gzipped when `compress_csv` is `TRUE`), plus the
#' serialized `ecdna_terminal_profiles.rds`,
#' `ecdna_static_id_copy_number_matrix_sparse.rds`,
#' `ecdna_recorder_edit_fraction_matrix_sparse.rds`, and
#' `ecdna_recorder_character_matrix_sparse.rds`.
#' @note Depends on `write_physicell_csv()` and `physicell_log_stage()` from
#'   `physicell_lineage.R`.
write_ecdna_outputs <- function(simulation,
                                model,
                                output_dir,
                                show_progress = FALSE,
                                compress_csv = TRUE){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  physicell_log_stage('ecDNA output: summarizing terminal profiles.', show_progress)
  tables <- ecdna_terminal_tables(simulation, model)
  matrices <- ecdna_feature_matrices(simulation, model)
  write_physicell_csv(
    tables$cell_summary,
    file.path(output_dir, 'ecdna_cell_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    tables$haplotypes,
    file.path(output_dir, 'ecdna_haplotypes.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    simulation$mutation_events,
    file.path(output_dir, 'ecdna_mutation_events.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  species_table <- data.frame(
    species_index = seq_len(model$num_species),
    species_id = model$species_ids,
    initially_labeled = seq_len(model$num_species) %in%
      model$labeled_species_indices,
    static_id = model$static_ids,
    initial_copies = model$initial_copies_per_species,
    stringsAsFactors = FALSE
  )
  write_physicell_csv(
    species_table,
    file.path(output_dir, 'ecdna_species_manifest.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  saveRDS(
    simulation$profiles[simulation$terminal_nodes$node_id],
    file.path(output_dir, 'ecdna_terminal_profiles.rds')
  )
  saveRDS(
    matrices$static_id_copy_numbers,
    file.path(output_dir, 'ecdna_static_id_copy_number_matrix_sparse.rds')
  )
  saveRDS(
    matrices$recorder_edit_fractions,
    file.path(output_dir, 'ecdna_recorder_edit_fraction_matrix_sparse.rds')
  )
  saveRDS(
    matrices$recorder_edit_presence,
    file.path(output_dir, 'ecdna_recorder_character_matrix_sparse.rds')
  )
  manifest <- data.frame(
    property = c(
      'num_sampled_terminal_cells', 'num_species', 'num_initially_labeled_species',
      'requested_labeled_species_fraction',
      'realized_labeled_species_fraction', 'num_recorder_targets',
      'recorder_start_time',
      'replication_probability', 'label_selection_coefficient',
      'recorder_edit_selection_coefficient',
      'daughter_1_segregation_probability',
      'max_copies_per_cell', 'num_aggregated_mutation_rows', 'random_seed'
    ),
    value = c(
      nrow(simulation$terminal_nodes), model$num_species,
      length(model$labeled_species_indices), model$labeled_species_fraction,
      model$realized_labeled_species_fraction,
      model$num_recorder_targets, model$recorder_start_time,
      model$replication_probability,
      model$label_selection_coefficient,
      model$recorder_edit_selection_coefficient,
      model$daughter_1_segregation_probability, model$max_copies_per_cell,
      nrow(simulation$mutation_events), simulation$seed
    ),
    stringsAsFactors = FALSE
  )
  write_physicell_csv(
    manifest,
    file.path(output_dir, 'ecdna_manifest.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  physicell_log_stage('ecDNA output complete.', show_progress)
  invisible(c(tables, matrices))
}
