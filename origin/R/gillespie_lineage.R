# Exact continuous-time population simulation for remote_mito.
#
# The population process is sampled with Gillespie's direct method. Barcode
# and mitochondrial mutations are conditionally independent Poisson processes
# on the resulting branches and are replayed by physicell_lineage.R and
# physicell_mito.R. Keeping those processes separate is exact while mutations
# do not alter division, death, or differentiation propensities.

#' Coerce one value to a validated numeric or integer scalar
#'
#' Used to normalize parameters that arrive from JSON or the command line, where
#' a value may be a string, a length-zero list, or `NA`.
#'
#' @export
#' @param value Candidate value of any type accepted by `as.numeric()` or
#'   `as.integer()`; coercion warnings are suppressed.
#' @param name Diagnostic name used in the error message.
#' @param mode Either `'numeric'` (the default) or `'integer'`, matched with
#'   `match.arg()`.
#' @return The coerced length-one value. Stops when the coercion does not yield
#'   exactly one non-`NA` value.
gillespie_scalar <- function(value, name, mode = c('numeric', 'integer')){
  mode <- match.arg(mode)
  converted <- if(mode == 'integer'){
    suppressWarnings(as.integer(value))
  } else{
    suppressWarnings(as.numeric(value))
  }
  if(length(converted) != 1L || is.na(converted)){
    stop(sprintf('%s must be one %s value.', name, mode))
  }
  converted
}

#' Convert a per-interval probability into a constant continuous-time hazard
#'
#' Solves `p = 1 - exp(-h * interval)` exactly for `h`, so integrating the
#' returned hazard over one interval reproduces the requested probability. This
#' is how per-cell-cycle probabilities in the parameter JSON become
#' per-unit-time propensities for the Gillespie draw.
#'
#' @export
#' @param probability Probability of at least one event during one interval;
#'   must be finite and in `[0, 1)`.
#' @param interval Interval length in simulation time units; must be positive,
#'   or `Inf` to mean "never".
#' @return The hazard `-log1p(-probability) / interval`, in events per unit
#'   time. An infinite interval returns `0` without inspecting `probability`.
gillespie_probability_hazard <- function(probability, interval){
  probability <- gillespie_scalar(probability, 'probability')
  interval <- gillespie_scalar(interval, 'interval')
  if(!is.finite(interval) || interval <= 0){
    if(is.infinite(interval) && interval > 0){
      return(0)
    }
    stop('interval must be positive (or Inf for a zero hazard).')
  }
  if(!is.finite(probability) || probability < 0 || probability >= 1){
    stop('probability must be in [0, 1).')
  }
  -log1p(-probability) / interval
}

#' Draw the cells that one induction event acts on
#'
#' Mirrors the selection semantics of the native timestep engine: a fixed count
#' is capped at the number of available cells, while a fraction is realized as a
#' binomial draw rather than a rounded count.
#'
#' @export
#' @param cell_indices Node indices of the currently active cells, coerced with
#'   `as.integer()`.
#' @param num_cells Optional non-negative whole number of cells to induce; takes
#'   precedence over `frac_cells` whenever it is supplied.
#' @param frac_cells Optional fraction in `[0, 1]`, consulted only when
#'   `num_cells` is `NULL`; the realized count is
#'   `rbinom(1, length(cell_indices), frac_cells)`.
#' @return An integer vector of selected indices sampled without replacement,
#'   possibly empty. Returns `integer()` when no cells are available.
gillespie_induction_selection <- function(cell_indices,
                                           num_cells = NULL,
                                           frac_cells = NULL){
  cell_indices <- as.integer(cell_indices)
  num_available <- length(cell_indices)
  if(num_available == 0L){
    return(integer())
  }
  if(!is.null(num_cells)){
    num_cells <- suppressWarnings(as.numeric(num_cells))
    if(length(num_cells) != 1L || !is.finite(num_cells) ||
       num_cells < 0 || num_cells %% 1 != 0){
      stop('Induction num_cells must be one non-negative integer.')
    }
    num_selected <- min(as.integer(num_cells), num_available)
  } else{
    frac_cells <- suppressWarnings(as.numeric(frac_cells))
    if(length(frac_cells) != 1L || !is.finite(frac_cells) ||
       frac_cells < 0 || frac_cells > 1){
      stop('Induction frac_cells must be one fraction between zero and one.')
    }
    num_selected <- stats::rbinom(1L, num_available, frac_cells)
  }
  if(num_selected == 0L){
    integer()
  } else{
    sample(cell_indices, num_selected, replace = FALSE)
  }
}

#' Normalize one induction block from the parameter JSON
#'
#' @export
#' @param params Parsed remote_mito parameter list.
#' @param name Name of the induction block, such as `'editing_induction'` or
#'   `'differentiation_induction'`.
#' @return A list with `timepoint` (a simulation time; `Inf` when the block is
#'   absent or its time does not parse, which disables the induction),
#'   `num_cells`, and `frac_cells`. A block carrying only a `timepoint` is
#'   normalized to `frac_cells = 1`, so every active cell is induced.
gillespie_induction_spec <- function(params, name){
  specification <- params[[name]]
  if(is.null(specification)){
    return(list(timepoint = Inf, num_cells = 0L, frac_cells = NULL))
  }
  timepoint <- suppressWarnings(as.numeric(specification$timepoint))
  if(length(timepoint) != 1L || is.na(timepoint)){
    timepoint <- Inf
  }
  if(is.null(specification$num_cells) && is.null(specification$frac_cells)){
    # An induction block containing only a time means all cells, matching the
    # PhysiCell adapter configuration used by the organoid examples.
    specification$frac_cells <- 1
  }
  list(
    timepoint = timepoint,
    num_cells = specification$num_cells,
    frac_cells = specification$frac_cells
  )
}

#' Build the induced and uninduced daughter cell-type transition matrices
#'
#' Row `i` of a transition matrix gives the probabilities that a dividing parent
#' of type `i` produces a daughter of each type. The two daughters of a division
#' are drawn independently from that row.
#'
#' @export
#' @param params Parsed parameter list; the matrices are read from
#'   `cell_type_dict$uninduced_transition_matrix` and
#'   `cell_type_dict$induced_transition_matrix`.
#' @param cell_types Character vector of cell-type names that fixes the row and
#'   column order of both matrices.
#' @return A list with `uninduced` and `induced`, each a square numeric matrix
#'   whose dimnames are `cell_types`. A block omitted from the JSON becomes the
#'   identity matrix, i.e. daughters keep the parent's type.
gillespie_transition_matrices <- function(params, cell_types){
  cell_type_dict <- params$cell_type_dict
  num_types <- length(cell_types)
  identity_matrix <- diag(num_types)
  dimnames(identity_matrix) <- list(cell_types, cell_types)

  # Internal: validate one square, row-stochastic transition matrix
  #
  # @param value Raw JSON value, unlisted and filled row by row
  #   (`byrow = TRUE`); `NULL` yields the identity matrix.
  # @param label Diagnostic name used in the error messages.
  # @return A `num_types` by `num_types` numeric matrix named by cell type.
  #   Stops when an entry is not finite and non-negative, or when a row does
  #   not sum to one within `1e-8`.
  parse_matrix <- function(value, label){
    if(is.null(value)){
      return(identity_matrix)
    }
    matrix_value <- matrix(
      as.numeric(unlist(value, use.names = FALSE)),
      nrow = num_types,
      byrow = TRUE
    )
    if(!identical(dim(matrix_value), c(num_types, num_types)) ||
       any(!is.finite(matrix_value)) || any(matrix_value < 0)){
      stop(sprintf(
        '%s must be a finite, non-negative %d by %d matrix.',
        label,
        num_types,
        num_types
      ))
    }
    row_totals <- rowSums(matrix_value)
    if(any(abs(row_totals - 1) > 1e-8)){
      stop(sprintf('Every row of %s must sum to one.', label))
    }
    dimnames(matrix_value) <- list(cell_types, cell_types)
    matrix_value
  }

  list(
    uninduced = parse_matrix(
      cell_type_dict$uninduced_transition_matrix,
      'cell_type_dict.uninduced_transition_matrix'
    ),
    induced = parse_matrix(
      cell_type_dict$induced_transition_matrix,
      'cell_type_dict.induced_transition_matrix'
    )
  )
}

#' Compute the per-cell division and death propensities of every cell type
#'
#' Both propensities are constant hazards in events per unit simulation time.
#' Division uses the reciprocal of the mean cell-cycle length, so the waiting
#' time to a division is exponential with that mean rather than deterministic as
#' in the timestep engine. Death converts a per-cell-cycle probability with
#' `gillespie_probability_hazard()`, so absent any competing division a cell
#' would die with exactly that probability over one cycle length.
#'
#' @export
#' @param params Parsed parameter list. Rates come from
#'   `cell_type_dict$cell_type_params[[type]]$cell_cycle_length` (must be
#'   positive, or `Inf` for a non-dividing type) and `death_per_cell_cycle_prob`
#'   (treated as `0` when absent).
#' @param cell_types Character vector of cell-type names.
#' @return A list of three numeric vectors named by cell type: `division`,
#'   `death`, and `total`, their sum. `total` is the single-cell propensity the
#'   Gillespie loop sums over the population.
gillespie_cell_type_rates <- function(params, cell_types){
  cell_type_params <- params$cell_type_dict$cell_type_params
  division <- death <- setNames(numeric(length(cell_types)), cell_types)
  for(cell_type in cell_types){
    type_params <- cell_type_params[[cell_type]]
    cycle_length <- suppressWarnings(as.numeric(type_params$cell_cycle_length))
    if(length(cycle_length) != 1L || is.na(cycle_length) || cycle_length <= 0){
      stop(sprintf(
        'cell_type_params.%s.cell_cycle_length must be positive or Inf.',
        cell_type
      ))
    }
    division[cell_type] <- if(is.infinite(cycle_length)) 0 else 1 / cycle_length
    death_probability <- type_params$death_per_cell_cycle_prob
    if(is.null(death_probability)){
      death_probability <- 0
    }
    death[cell_type] <- gillespie_probability_hazard(
      death_probability,
      cycle_length
    )
  }
  list(division = division, death = death, total = division + death)
}

#' Prepare one barcode recording model per cell type on a shared layout
#'
#' Barcode edits are replayed on the branches after the population process has
#' been simulated, so every cell type must agree on one recorder reference,
#' target layout, and integration count and may differ only in its editing
#' hazards. This builds one `prepare_physicell_recording_model()` result per
#' configured cell type and rejects any set that is not layout-compatible.
#'
#' @export
#' @param params Parsed parameter list; cell types are taken from
#'   `cell_type_dict$cell_type_params`.
#' @param num_integrations Optional number of integrated barcode recorders;
#'   `NULL` leaves the value the recording model derives from `params`.
#' @param founder_label_sites Number of stable founder-label sites; default `0`.
#' @param params_dir Directory used to resolve paths referenced by the
#'   parameters, such as a fixed barcode sequence file.
#' @param seed Seed handed to every per-type model so they all draw the same
#'   reference and layout.
#' @return A named list of prepared barcode models, one entry per cell type.
#' @note Requires `physicell_lineage.R` to have been sourced. Stops when the
#'   per-type models disagree on barcode sequence, base-editor targets, nuclease
#'   target classes, integration count, or profile storage.
prepare_gillespie_barcode_models <- function(params,
                                              num_integrations = NULL,
                                              founder_label_sites = 0L,
                                              params_dir = '.',
                                              seed = 1L){
  if(!exists('prepare_physicell_recording_model', mode = 'function')){
    stop('Source physicell_lineage.R before preparing barcode models.')
  }
  cell_types <- names(params$cell_type_dict$cell_type_params)
  models <- setNames(lapply(cell_types, function(cell_type){
    prepare_physicell_recording_model(
      params,
      cell_type = cell_type,
      num_integrations = num_integrations,
      founder_label_sites = founder_label_sites,
      params_dir = params_dir,
      seed = seed
    )
  }), cell_types)
  reference_model <- models[[1]]
  compatible <- vapply(models, function(model){
    identical(model$barcode_sequence, reference_model$barcode_sequence) &&
      identical(model$be_targets, reference_model$be_targets) &&
      identical(model$nuc_target_classes, reference_model$nuc_target_classes) &&
      identical(model$num_integrations, reference_model$num_integrations) &&
      identical(model$profile_storage, reference_model$profile_storage)
  }, logical(1))
  if(any(!compatible)){
    stop('Cell-type barcode models do not share one compatible recorder layout.')
  }
  models
}

#' Prepare one mitochondrial mutation model per cell type on a shared reference
#'
#' As with the barcode models, mitochondrial mutations are replayed on branches
#' afterwards, so all cell types must share one reference genome and bottleneck
#' size and differ only in their mutation hazards.
#'
#' @export
#' @param params Parsed parameter list; cell types are taken from
#'   `cell_type_dict$cell_type_params`.
#' @param genomes_per_cell Number of mitochondrial genomes carried through each
#'   division bottleneck; default `8`.
#' @param seed Seed handed to every per-type model so they all draw the same
#'   mitochondrial reference.
#' @return A named list of prepared mitochondrial models, one entry per cell
#'   type.
#' @note Requires `physicell_mito.R` to have been sourced. Stops when the
#'   per-type models disagree on the mitochondrial reference or on
#'   `genomes_per_cell`.
prepare_gillespie_mito_models <- function(params,
                                           genomes_per_cell = 8L,
                                           seed = 1L){
  if(!exists('prepare_physicell_mito_model', mode = 'function')){
    stop('Source physicell_mito.R before preparing mitochondrial models.')
  }
  cell_types <- names(params$cell_type_dict$cell_type_params)
  models <- setNames(lapply(cell_types, function(cell_type){
    prepare_physicell_mito_model(
      params,
      cell_type = cell_type,
      genomes_per_cell = genomes_per_cell,
      seed = seed
    )
  }), cell_types)
  reference_model <- models[[1]]
  compatible <- vapply(models, function(model){
    identical(
      model$mitochondrial_reference,
      reference_model$mitochondrial_reference
    ) && identical(model$genomes_per_cell, reference_model$genomes_per_cell)
  }, logical(1))
  if(any(!compatible)){
    stop('Cell-type mitochondrial models do not share one compatible reference.')
  }
  models
}

#' Flatten the active-cell slots into a plain node-index vector
#'
#' @export
#' @param active_nodes Preallocated integer vector whose first `active_count`
#'   entries hold the active node indices. A list of index groups is also
#'   accepted, for compatibility with the earlier grouped representation.
#' @param active_count Number of slots currently in use; ignored for the list
#'   form.
#' @return An integer vector of active node indices, empty when no slots are in
#'   use.
gillespie_active_indices <- function(active_nodes, active_count = NULL){
  if(is.list(active_nodes)){
    return(as.integer(unlist(active_nodes, use.names = FALSE)))
  }
  if(is.null(active_count) || active_count == 0L){
    return(integer())
  }
  as.integer(active_nodes[seq_len(active_count)])
}

#' Simulate an exact continuous-time birth, death, and differentiation lineage
#'
#' Samples the population process with Gillespie's direct method and returns the
#' complete lineage as a node table whose branch lengths are in simulation time
#' units. No mutations are simulated here; the barcode and mitochondrial Poisson
#' processes are replayed on these branches afterwards.
#'
#' @export
#' @details
#' Each active cell carries exactly two reactions, both constant-hazard Poisson
#' processes measured in events per unit simulation time and both supplied by
#' `gillespie_cell_type_rates()`: division at `1 / cell_cycle_length` and death
#' at `-log(1 - death_per_cell_cycle_prob) / cell_cycle_length`. A cell's total
#' propensity is their sum and depends only on its cell type, so the population
#' total changes only when cells are added or removed, never as time passes.
#'
#' One step of the loop draws the waiting time to the next event as
#' `rexp(1, total_propensity)`, where `total_propensity` is the sum over every
#' active cell, maintained incrementally next to a Fenwick (binary indexed) tree
#' over the active slots. When the drawn time falls at or beyond the next
#' boundary the clock instead advances to that boundary with no event, so
#' induction times and `sim_length` checkpoints are hit exactly rather than
#' being straddled. Otherwise a uniform draw on `[0, total_propensity]` picks
#' the reacting cell from the Fenwick tree in logarithmic time, and a Bernoulli
#' draw with probability `division / (division + death)` for that cell's type
#' decides whether the event is a division or a death.
#'
#' A division retires the parent node at the event time and appends two daughter
#' nodes whose types are drawn independently from the parent's row of the
#' induced or uninduced transition matrix. A death retires the node and flags it
#' `died`. An induction likewise retires the node and appends a single
#' continuation node that keeps the same cell id, type, and generation but
#' carries the new induction flags, so the induced state is inherited by every
#' later descendant.
#'
#' @param params Parsed remote_mito parameter list. Reads `random_seed`,
#'   `sim_length`, `num_init_cells`, `cell_type_dict` (founder type, per-type
#'   rates, transition matrices), `editing_induction`, and
#'   `differentiation_induction`.
#' @param end_time Time at which the simulation stops; defaults to
#'   `max(sim_length)`. Must be finite and non-negative.
#' @param seed Integer random seed; defaults to `params$random_seed`.
#' @param max_cells Hard cap on simultaneously active cells. Reaching it ends
#'   the run with a warning and `stop_reason = 'max_cells'`. It also sizes the
#'   preallocated active-slot and Fenwick vectors. Default `1000000`.
#' @param show_progress Whether to print throttled progress lines.
#' @param progress_updates Number of evenly spaced times at which progress is
#'   reported; must be positive.
#' @return A list with `nodes` (one row per lineage node, carrying `birth_time`,
#'   `end_time`, `branch_length`, `division_event`, `origin`, `is_terminal`,
#'   `alive_at_end`, `died`, `cell_type`, `generation`, `editing_state`, and
#'   `differentiation_induced`), `terminal_nodes` (the rows alive at the final
#'   time), `event_log`, `division_events`, `checkpoint_summary`, `end_time`
#'   (the time actually reached), `requested_end_time`, `stop_reason`, `seed`, a
#'   `counts` vector, and the `rates` list.
#' @section Side effects: Calls `set.seed(seed)`, replacing the caller's RNG
#'   stream. Prints progress lines when `show_progress` is set, and warns when
#'   `max_cells` is reached or no cell survives to the final time.
#' @note Heteroplasmy-dependent survival is rejected up front: it would couple
#'   mitochondrial mutations to the population propensities and invalidate the
#'   separation that lets mutations be replayed after the fact.
simulate_gillespie_population <- function(params,
                                          end_time = NULL,
                                          seed = NULL,
                                          max_cells = 1000000L,
                                          show_progress = TRUE,
                                          progress_updates = 20L){
  if(!is.list(params)){
    stop('params must be a parsed remote_mito parameter list.')
  }
  if(isTRUE(params$consider_cell_heteroplasmy_scores)){
    stop(paste(
      'The Gillespie engine does not yet support heteroplasmy-dependent',
      'survival because that couples mitochondrial mutations to population',
      'propensities. Set consider_cell_heteroplasmy_scores to false.'
    ))
  }
  if(is.null(seed)){
    seed <- params$random_seed
  }
  seed <- gillespie_scalar(seed, 'seed', 'integer')
  set.seed(seed)
  if(is.null(end_time)){
    end_time <- max(as.numeric(unlist(params$sim_length, use.names = FALSE)))
  }
  end_time <- gillespie_scalar(end_time, 'end_time')
  if(!is.finite(end_time) || end_time < 0){
    stop('end_time must be one finite non-negative number.')
  }
  max_cells <- gillespie_scalar(max_cells, 'max_cells', 'integer')
  if(max_cells < 1L){
    stop('max_cells must be positive.')
  }
  progress_updates <- gillespie_scalar(
    progress_updates,
    'progress_updates',
    'integer'
  )
  if(progress_updates < 1L){
    stop('progress_updates must be positive.')
  }

  cell_type_params <- params$cell_type_dict$cell_type_params
  cell_types <- names(cell_type_params)
  if(length(cell_types) == 0L || any(!nzchar(cell_types))){
    stop('cell_type_dict.cell_type_params must name at least one cell type.')
  }
  founder_type <- as.character(params$cell_type_dict$founder_cell_type)
  if(length(founder_type) != 1L || !(founder_type %in% cell_types)){
    stop('cell_type_dict.founder_cell_type must name a configured cell type.')
  }
  rates <- gillespie_cell_type_rates(params, cell_types)
  transitions <- gillespie_transition_matrices(params, cell_types)
  num_founders <- gillespie_scalar(
    params$num_init_cells,
    'num_init_cells',
    'integer'
  )
  if(num_founders < 1L){
    stop('num_init_cells must be positive.')
  }
  if(num_founders > max_cells){
    stop('num_init_cells exceeds the Gillespie max_cells safety limit.')
  }

  edit_spec <- gillespie_induction_spec(params, 'editing_induction')
  diff_spec <- gillespie_induction_spec(params, 'differentiation_induction')

  # An induction specified as a fixed num_cells induces that many cells and no
  # more. With one founder that covers the whole population, so the same
  # parameter file behaves completely differently once num_init_cells is
  # raised: the uninduced founders and all their descendants simply never
  # record. Nothing fails, the run just carries a diluted signal, so warn
  # rather than error -- inducing a deliberate subset is legitimate.
  warn_partial_founder_induction <- function(spec, name){
    if(is.null(spec$num_cells) || !is.finite(spec$timepoint)){
      return(invisible(NULL))
    }
    requested <- suppressWarnings(as.integer(spec$num_cells))
    if(is.na(requested) || requested >= num_founders){
      return(invisible(NULL))
    }
    # Only a concern while the population is still just the founders; after
    # they divide, inducing a subset is an ordinary experimental design.
    if(spec$timepoint > 0){
      return(invisible(NULL))
    }
    warning(sprintf(
      paste0(
        '%s.num_cells is %d but num_init_cells is %d, so %d of the %d ',
        'founders and all of their descendants will never be induced. Set ',
        '%s.frac_cells = 1 (and num_cells = NULL) to induce every founder.'
      ),
      name, requested, num_founders, num_founders - requested, num_founders,
      name
    ), call. = FALSE)
  }
  warn_partial_founder_induction(edit_spec, 'editing_induction')
  warn_partial_founder_induction(diff_spec, 'differentiation_induction')
  stopping_points <- sort(unique(as.numeric(unlist(
    params$sim_length,
    use.names = FALSE
  ))))
  stopping_points <- stopping_points[
    is.finite(stopping_points) & stopping_points >= 0 & stopping_points <= end_time
  ]
  boundaries <- sort(unique(c(
    stopping_points,
    edit_spec$timepoint[
      is.finite(edit_spec$timepoint) & edit_spec$timepoint >= 0 &
        edit_spec$timepoint <= end_time
    ],
    diff_spec$timepoint[
      is.finite(diff_spec$timepoint) & diff_spec$timepoint >= 0 &
        diff_spec$timepoint <= end_time
    ],
    end_time
  )))

  # Node arrays grow geometrically. Node indices are also stable handles in the
  # active-cell vectors, avoiding hash lookup and linear deletion costs.
  node_capacity <- max(1024L, min(max_cells * 2L, num_founders * 4L))
  node_id <- rep(NA_character_, node_capacity)
  physicell_id_value <- rep(NA_character_, node_capacity)
  parent_node_id <- rep(NA_character_, node_capacity)
  birth_time <- numeric(node_capacity)
  end_time_value <- rep(NA_real_, node_capacity)
  division_event <- rep(NA_integer_, node_capacity)
  origin <- rep(NA_character_, node_capacity)
  cell_type_value <- rep(NA_character_, node_capacity)
  generation <- integer(node_capacity)
  editing_induced <- logical(node_capacity)
  differentiation_induced <- logical(node_capacity)
  died <- logical(node_capacity)
  child_count <- integer(node_capacity)
  active_position <- integer(node_capacity)
  node_count <- 0L

  # Internal: grow the preallocated node columns to hold `required` nodes
  #
  # Grows geometrically (1.5x, never past `.Machine$integer.max`) and
  # explicitly initializes every newly added slot of all fourteen parallel node
  # vectors so no slot is left as the recycled value of an extended vector.
  #
  # @param required Number of node slots that must exist.
  # @return `invisible(NULL)`; `node_capacity` and the node vectors are updated
  #   in the enclosing frame. Stops when the capacity cannot be grown further.
  grow_nodes <- function(required){
    if(required <= node_capacity){
      return(invisible(NULL))
    }
    new_capacity <- node_capacity
    while(new_capacity < required){
      new_capacity <- min(
        .Machine$integer.max,
        max(required, ceiling(new_capacity * 1.5))
      )
      if(new_capacity <= node_capacity){
        stop('Gillespie node capacity exceeded.')
      }
    }
    old_capacity <- node_capacity
    length(node_id) <<- new_capacity
    length(physicell_id_value) <<- new_capacity
    length(parent_node_id) <<- new_capacity
    length(birth_time) <<- new_capacity
    length(end_time_value) <<- new_capacity
    length(division_event) <<- new_capacity
    length(origin) <<- new_capacity
    length(cell_type_value) <<- new_capacity
    length(generation) <<- new_capacity
    length(editing_induced) <<- new_capacity
    length(differentiation_induced) <<- new_capacity
    length(died) <<- new_capacity
    length(child_count) <<- new_capacity
    length(active_position) <<- new_capacity
    new_indices <- seq.int(old_capacity + 1L, new_capacity)
    node_id[new_indices] <<- NA_character_
    physicell_id_value[new_indices] <<- NA_character_
    parent_node_id[new_indices] <<- NA_character_
    birth_time[new_indices] <<- 0
    end_time_value[new_indices] <<- NA_real_
    division_event[new_indices] <<- NA_integer_
    origin[new_indices] <<- NA_character_
    cell_type_value[new_indices] <<- NA_character_
    generation[new_indices] <<- 0L
    editing_induced[new_indices] <<- FALSE
    differentiation_induced[new_indices] <<- FALSE
    died[new_indices] <<- FALSE
    child_count[new_indices] <<- 0L
    active_position[new_indices] <<- 0L
    node_capacity <<- new_capacity
    invisible(NULL)
  }

  active_nodes <- integer(max_cells)
  propensity_tree <- numeric(max_cells)
  active_count <- 0L
  active_total_propensity <- 0
  next_cell_id <- 0L
  # Internal: fold a propensity delta into one Fenwick-tree slot
  #
  # @param position Active slot index, in `1:max_cells`.
  # @param difference Signed propensity change added to every tree node that
  #   covers `position`.
  # @return `invisible(NULL)`; `propensity_tree` is updated in the enclosing
  #   frame in `O(log max_cells)` time.
  update_propensity_tree <- function(position, difference){
    while(position <= max_cells){
      propensity_tree[position] <<- propensity_tree[position] + difference
      position <- position + bitwAnd(position, -position)
    }
    invisible(NULL)
  }
  # Internal: find the active slot that owns one point of the propensity mass
  #
  # Descends the Fenwick tree by decreasing powers of two, so the weighted
  # choice costs `O(log max_cells)` instead of a linear scan of the population.
  #
  # @param target Uniform draw on the total active propensity.
  # @return The one-based active slot index whose cumulative propensity spans
  #   `target`. Stops when that slot falls outside `1:active_count`, which
  #   would mean the tree and the running total have diverged.
  select_propensity_slot <- function(target){
    tree_index <- 0L
    bit <- as.integer(2 ^ floor(log(max_cells, base = 2)))
    while(bit >= 1L){
      candidate <- tree_index + bit
      if(candidate <= max_cells && propensity_tree[candidate] <= target){
        tree_index <- candidate
        target <- target - propensity_tree[candidate]
      }
      bit <- bit %/% 2L
    }
    selected_slot <- tree_index + 1L
    if(selected_slot < 1L || selected_slot > active_count){
      stop('Internal Gillespie propensity selection is inconsistent.')
    }
    selected_slot
  }
  # Internal: make one node an active cell
  #
  # Appends the node to the active slots, records its reverse position, and
  # adds its cell type's total propensity to both the Fenwick tree and the
  # running population total.
  #
  # @param node_index Node index to activate.
  # @return The node index, invisibly. Stops when all `max_cells` slots are in
  #   use.
  add_active <- function(node_index){
    if(active_count >= max_cells){
      stop('Internal active-cell capacity exceeded max_cells.')
    }
    active_count <<- active_count + 1L
    active_nodes[active_count] <<- node_index
    active_position[node_index] <<- active_count
    propensity <- rates$total[cell_type_value[node_index]]
    update_propensity_tree(active_count, propensity)
    active_total_propensity <<- active_total_propensity + propensity
    invisible(node_index)
  }
  # Internal: retire one active cell by swapping the last slot into its place
  #
  # Keeps `active_nodes`, `active_position`, `propensity_tree`, and
  # `active_total_propensity` mutually consistent, and clamps a residual total
  # whose magnitude has fallen below `1e-12` to exactly zero.
  #
  # @param node_index Node index to deactivate.
  # @return The node index, invisibly. Stops when the node is not recorded as
  #   active at its stored slot.
  remove_active <- function(node_index){
    position <- active_position[node_index]
    if(position < 1L || position > active_count ||
       active_nodes[position] != node_index){
      stop('Internal active-cell index is inconsistent.')
    }
    removed_propensity <- rates$total[cell_type_value[node_index]]
    last_position <- active_count
    last_node <- active_nodes[last_position]
    update_propensity_tree(position, -removed_propensity)
    if(position < last_position){
      last_propensity <- rates$total[cell_type_value[last_node]]
      update_propensity_tree(last_position, -last_propensity)
      active_nodes[position] <<- last_node
      active_position[last_node] <<- position
      update_propensity_tree(position, last_propensity)
    }
    active_nodes[last_position] <<- 0L
    active_position[node_index] <<- 0L
    active_count <<- active_count - 1L
    active_total_propensity <<- active_total_propensity - removed_propensity
    if(abs(active_total_propensity) < 1e-12){
      active_total_propensity <<- 0
    }
    invisible(node_index)
  }
  # Internal: append one lineage node and immediately activate it
  #
  # Nodes are appended parent before child, so the finished table is already in
  # topological order.
  #
  # @param cell_id Cell identifier stored as `physicell_id`; a continuation
  #   node deliberately reuses its parent's identifier.
  # @param parent_index Parent node index, or `NA_integer_` for a founder.
  #   Supplying it also increments the parent's child count.
  # @param time Birth time of the new node.
  # @param type Cell-type name.
  # @param node_origin Provenance label: `'founder'`, `'daughter_1'`,
  #   `'daughter_2'`, or `'induction_continuation'`.
  # @param node_generation Number of divisions separating the node from its
  #   founder.
  # @param is_editing_induced Whether the cell carries the editing induction.
  # @param is_differentiation_induced Whether the cell carries the
  #   differentiation induction.
  # @param division_index Index of the division that created the node, or
  #   `NA_integer_` for founders and continuations.
  # @return The index of the new node.
  add_node <- function(cell_id,
                       parent_index = NA_integer_,
                       time,
                       type,
                       node_origin,
                       node_generation,
                       is_editing_induced,
                       is_differentiation_induced,
                       division_index = NA_integer_){
    node_count <<- node_count + 1L
    grow_nodes(node_count)
    node_id[node_count] <<- paste0('node_', node_count)
    physicell_id_value[node_count] <<- as.character(cell_id)
    parent_node_id[node_count] <<- if(is.na(parent_index)){
      NA_character_
    } else{
      child_count[parent_index] <<- child_count[parent_index] + 1L
      node_id[parent_index]
    }
    birth_time[node_count] <<- time
    division_event[node_count] <<- division_index
    origin[node_count] <<- node_origin
    cell_type_value[node_count] <<- type
    generation[node_count] <<- node_generation
    editing_induced[node_count] <<- is_editing_induced
    differentiation_induced[node_count] <<- is_differentiation_induced
    add_active(node_count)
    node_count
  }

  event_rows <- list()
  event_count <- 0L
  # Internal: append one typed row to the population event log
  #
  # @param time Time of the event.
  # @param event Event label: `'division'`, `'death'`, `'editing_induction'`,
  #   or `'differentiation_induction'`.
  # @param node_index Node the event is attributed to; its id, cell id, parent
  #   id, and cell type are copied into the row.
  # @param detail Free-text detail column, empty by default.
  # @return The appended one-row data frame, as the value of the assignment;
  #   `event_rows` and `event_count` are updated in the enclosing frame.
  append_event <- function(time, event, node_index, detail = ''){
    event_count <<- event_count + 1L
    event_rows[[event_count]] <<- data.frame(
      time = time,
      event = event,
      node_id = node_id[node_index],
      physicell_id = physicell_id_value[node_index],
      parent_node_id = parent_node_id[node_index],
      cell_type = cell_type_value[node_index],
      detail = detail,
      stringsAsFactors = FALSE
    )
  }
  division_rows <- list()
  division_row_count <- 0L
  checkpoint_rows <- list()
  checkpoint_row_count <- 0L

  for(founder_index in seq_len(num_founders)){
    next_cell_id <- next_cell_id + 1L
    add_node(
      cell_id = next_cell_id,
      time = 0,
      type = founder_type,
      node_origin = 'founder',
      node_generation = 0L,
      is_editing_induced = FALSE,
      is_differentiation_induced = FALSE
    )
  }

  # Internal: choose which currently active cells one induction applies to
  #
  # @param specification Normalized induction spec from
  #   `gillespie_induction_spec()`.
  # @return Integer node indices chosen by `gillespie_induction_selection()`.
  select_active <- function(specification){
    gillespie_induction_selection(
      gillespie_active_indices(active_nodes, active_count),
      specification$num_cells,
      specification$frac_cells
    )
  }
  # Internal: apply every induction whose timepoint coincides with `time`
  #
  # Editing and differentiation inductions that land on the same boundary are
  # resolved together, so a cell selected by both gains one continuation node
  # carrying both flags rather than a zero-length chain of two nodes. Each
  # selected cell is retired at `time` and replaced by a continuation node
  # holding the union of its old and newly applied induction flags.
  #
  # @param time Boundary time, matched to an induction timepoint within
  #   `1e-10`.
  # @return `invisible(NULL)`; the node table, active slots, and event log are
  #   updated in the enclosing frame.
  apply_inductions <- function(time){
    do_edit <- is.finite(edit_spec$timepoint) &&
      abs(time - edit_spec$timepoint) <= 1e-10
    do_diff <- is.finite(diff_spec$timepoint) &&
      abs(time - diff_spec$timepoint) <= 1e-10
    if(!do_edit && !do_diff){
      return(invisible(NULL))
    }
    edit_selected <- if(do_edit) select_active(edit_spec) else integer()
    diff_selected <- if(do_diff) select_active(diff_spec) else integer()
    selected <- union(edit_selected, diff_selected)
    if(length(selected) == 0L){
      return(invisible(NULL))
    }
    # A single continuation node handles coincident induction boundaries and
    # prevents artificial zero-length chains.
    for(old_index in selected){
      old_cell_id <- physicell_id_value[old_index]
      old_type <- cell_type_value[old_index]
      old_generation <- generation[old_index]
      new_edit <- editing_induced[old_index] || old_index %in% edit_selected
      new_diff <- differentiation_induced[old_index] || old_index %in% diff_selected
      remove_active(old_index)
      end_time_value[old_index] <<- time
      continuation_index <- add_node(
        cell_id = old_cell_id,
        parent_index = old_index,
        time = time,
        type = old_type,
        node_origin = 'induction_continuation',
        node_generation = old_generation,
        is_editing_induced = new_edit,
        is_differentiation_induced = new_diff
      )
      if(old_index %in% edit_selected){
        append_event(time, 'editing_induction', continuation_index)
      }
      if(old_index %in% diff_selected){
        append_event(time, 'differentiation_induction', continuation_index)
      }
    }
    invisible(NULL)
  }
  # Internal: record the active-cell counts per cell type at one time
  #
  # @param time Checkpoint time, one of the `sim_length` stopping points.
  # @return `invisible(NULL)`; appends one row per cell type to
  #   `checkpoint_rows`, each carrying that type's active count alongside the
  #   total over all types, the node count, and the population-event count.
  record_checkpoint <- function(time){
    active_indices <- gillespie_active_indices(active_nodes, active_count)
    counts <- tabulate(
      match(cell_type_value[active_indices], cell_types),
      nbins = length(cell_types)
    )
    for(type_index in seq_along(cell_types)){
      checkpoint_row_count <<- checkpoint_row_count + 1L
      checkpoint_rows[[checkpoint_row_count]] <<- data.frame(
        time = time,
        cell_type = cell_types[type_index],
        active_cells = counts[type_index],
        total_active_cells = sum(counts),
        total_nodes = node_count,
        population_events = event_count,
        stringsAsFactors = FALSE
      )
    }
    invisible(NULL)
  }

  # Time-zero induction is applied before any stochastic event.
  apply_inductions(0)
  if(any(abs(stopping_points) <= 1e-10)){
    record_checkpoint(0)
  }
  boundary_index <- which(boundaries > 0)[1]
  if(is.na(boundary_index)){
    boundary_index <- length(boundaries) + 1L
  }
  current_time <- 0
  stochastic_events <- 0L
  divisions <- 0L
  deaths <- 0L
  stop_reason <- 'end_time'
  if(active_count >= max_cells && end_time > 0){
    stop_reason <- 'max_cells'
    warning(sprintf(
      'Gillespie simulation started at the max_cells=%d safety limit.',
      max_cells
    ), call. = FALSE)
  }
  progress_thresholds <- if(end_time > 0){
    seq(0, end_time, length.out = progress_updates + 1L)[-1]
  } else{
    numeric()
  }
  progress_index <- 1L
  # Internal: emit any progress lines whose time threshold has been passed
  #
  # @return `invisible(NULL)`; prints one throttled line per crossed threshold,
  #   reporting the current time, active-cell count, and stochastic event
  #   count. Does nothing unless `show_progress` is set and `end_time` is
  #   positive.
  report_time_progress <- function(){
    if(!isTRUE(show_progress) || end_time <= 0){
      return(invisible(NULL))
    }
    while(progress_index <= length(progress_thresholds) &&
          current_time >= progress_thresholds[progress_index] - 1e-10){
      cat(sprintf(
        '[Gillespie] time %.6g / %.6g (%.0f%%); active=%s; events=%s\n',
        current_time,
        end_time,
        100 * current_time / end_time,
        format(active_count, big.mark = ',', scientific = FALSE),
        format(stochastic_events, big.mark = ',', scientific = FALSE)
      ))
      flush.console()
      progress_index <<- progress_index + 1L
    }
    invisible(NULL)
  }

  while(current_time < end_time - 1e-12 && stop_reason == 'end_time'){
    next_boundary <- if(boundary_index <= length(boundaries)){
      boundaries[boundary_index]
    } else{
      end_time
    }
    total_propensity <- active_total_propensity
    candidate_time <- if(total_propensity > 0){
      current_time + stats::rexp(1L, total_propensity)
    } else{
      Inf
    }

    if(candidate_time >= next_boundary){
      current_time <- next_boundary
      apply_inductions(current_time)
      if(any(abs(stopping_points - current_time) <= 1e-10)){
        record_checkpoint(current_time)
      }
      boundary_index <- boundary_index + 1L
      report_time_progress()
      next
    }

    current_time <- candidate_time
    selected_slot <- select_propensity_slot(
      stats::runif(1L, min = 0, max = total_propensity)
    )
    selected_node <- active_nodes[selected_slot]
    selected_type <- cell_type_value[selected_node]
    division_probability <- rates$division[selected_type] /
      rates$total[selected_type]
    is_division <- stats::runif(1L) < division_probability
    remove_active(selected_node)
    end_time_value[selected_node] <- current_time
    stochastic_events <- stochastic_events + 1L

    if(is_division){
      divisions <- divisions + 1L
      append_event(current_time, 'division', selected_node)
      transition_matrix <- if(differentiation_induced[selected_node]){
        transitions$induced
      } else{
        transitions$uninduced
      }
      daughter_types <- sample(
        cell_types,
        2L,
        replace = TRUE,
        prob = transition_matrix[selected_type, ]
      )
      for(daughter_number in 1:2){
        next_cell_id <- next_cell_id + 1L
        daughter_index <- add_node(
          cell_id = next_cell_id,
          parent_index = selected_node,
          time = current_time,
          type = daughter_types[daughter_number],
          node_origin = paste0('daughter_', daughter_number),
          node_generation = generation[selected_node] + 1L,
          is_editing_induced = editing_induced[selected_node],
          is_differentiation_induced = differentiation_induced[selected_node],
          division_index = divisions
        )
        division_row_count <- division_row_count + 1L
        division_rows[[division_row_count]] <- data.frame(
          time = current_time,
          parent_ID = physicell_id_value[selected_node],
          daughter_ID = physicell_id_value[daughter_index],
          parent_node_id = node_id[selected_node],
          daughter_node_id = node_id[daughter_index],
          daughter_cell_type = cell_type_value[daughter_index],
          stringsAsFactors = FALSE
        )
      }
    } else{
      deaths <- deaths + 1L
      died[selected_node] <- TRUE
      append_event(current_time, 'death', selected_node)
    }

    if(active_count >= max_cells){
      stop_reason <- 'max_cells'
      warning(sprintf(
        'Gillespie simulation stopped at %.6g after reaching max_cells=%d.',
        current_time,
        max_cells
      ), call. = FALSE)
      break
    }
    report_time_progress()
  }

  final_time <- if(stop_reason == 'max_cells') current_time else end_time
  final_active <- gillespie_active_indices(active_nodes, active_count)
  end_time_value[final_active] <- final_time
  indices <- seq_len(node_count)
  branch_length <- end_time_value[indices] - birth_time[indices]
  if(any(!is.finite(branch_length)) || any(branch_length < -1e-10)){
    invalid_index <- which(!is.finite(branch_length) | branch_length < -1e-10)[1]
    stop(sprintf(
      paste(
        'The Gillespie engine produced an invalid branch length for %s',
        '(birth=%s, end=%s, active=%s, origin=%s).'
      ),
      node_id[invalid_index],
      birth_time[invalid_index],
      end_time_value[invalid_index],
      invalid_index %in% final_active,
      origin[invalid_index]
    ))
  }
  branch_length[abs(branch_length) < 1e-12] <- 0
  is_leaf <- child_count[indices] == 0L
  alive_at_end <- indices %in% final_active
  nodes <- data.frame(
    node_id = node_id[indices],
    physicell_id = physicell_id_value[indices],
    parent_node_id = parent_node_id[indices],
    birth_time = birth_time[indices],
    end_time = end_time_value[indices],
    branch_length = branch_length,
    division_event = division_event[indices],
    origin = origin[indices],
    is_terminal = is_leaf,
    alive_at_end = alive_at_end,
    died = died[indices],
    cell_type = cell_type_value[indices],
    generation = generation[indices],
    editing_state = ifelse(
      editing_induced[indices],
      'induced',
      'uninduced'
    ),
    differentiation_induced = differentiation_induced[indices],
    stringsAsFactors = FALSE
  )
  terminal_nodes <- nodes[nodes$alive_at_end, , drop = FALSE]
  if(nrow(terminal_nodes) == 0L){
    warning('No cells survived to the final Gillespie time.', call. = FALSE)
  }
  event_log <- if(length(event_rows) == 0L){
    data.frame(
      time = numeric(), event = character(), node_id = character(),
      physicell_id = character(), parent_node_id = character(),
      cell_type = character(), detail = character(),
      stringsAsFactors = FALSE
    )
  } else{
    do.call(rbind, event_rows)
  }
  division_events <- if(length(division_rows) == 0L){
    data.frame(
      time = numeric(), parent_ID = character(), daughter_ID = character(),
      parent_node_id = character(), daughter_node_id = character(),
      daughter_cell_type = character(), stringsAsFactors = FALSE
    )
  } else{
    do.call(rbind, division_rows)
  }
  checkpoint_summary <- if(length(checkpoint_rows) == 0L){
    data.frame(
      time = numeric(), cell_type = character(), active_cells = integer(),
      total_active_cells = integer(), total_nodes = integer(),
      population_events = integer(), stringsAsFactors = FALSE
    )
  } else{
    do.call(rbind, checkpoint_rows)
  }

  list(
    nodes = nodes,
    terminal_nodes = terminal_nodes,
    event_log = event_log,
    division_events = division_events,
    checkpoint_summary = checkpoint_summary,
    end_time = final_time,
    requested_end_time = end_time,
    stop_reason = stop_reason,
    seed = seed,
    counts = c(
      founders = num_founders,
      divisions = divisions,
      deaths = deaths,
      stochastic_events = stochastic_events,
      nodes = node_count,
      surviving_cells = nrow(terminal_nodes)
    ),
    rates = rates
  )
}

#' Build a PhysiCell-compatible covariate table for the surviving cells
#'
#' Supplies the per-cell covariates the downstream single-cell helpers expect
#' from the PhysiCell adapter. The Gillespie engine has no spatial model, so the
#' coordinate and neighbour columns exist but are constant.
#'
#' @export
#' @param simulation Completed result of `simulate_gillespie_population()`.
#' @return A data frame with one row per node alive at the final time and the
#'   columns `ID`, `parent_ID` (the parent node's cell id, `NA` for founders),
#'   `x`, `y`, `z` (all zero), `neighbor_IDs` (empty strings), `cell_type`,
#'   `developmental_pseudotime` (the node's generation divided by the largest
#'   generation anywhere in the lineage, or zero when nothing divided),
#'   `generation`, `birth_time`, and `end_time`. A simulation with no surviving
#'   cells returns those columns with zero rows.
gillespie_sc_cell_states <- function(simulation){
  terminal_nodes <- simulation$terminal_nodes
  if(nrow(terminal_nodes) == 0L){
    return(data.frame(
      ID = character(), parent_ID = character(), x = numeric(), y = numeric(),
      z = numeric(), neighbor_IDs = character(), cell_type = character(),
      developmental_pseudotime = numeric(), generation = integer(),
      birth_time = numeric(), end_time = numeric(), stringsAsFactors = FALSE
    ))
  }
  node_parent_index <- match(
    terminal_nodes$parent_node_id,
    simulation$nodes$node_id
  )
  maximum_generation <- max(simulation$nodes$generation)
  pseudotime <- if(maximum_generation > 0){
    terminal_nodes$generation / maximum_generation
  } else{
    rep(0, nrow(terminal_nodes))
  }
  data.frame(
    ID = terminal_nodes$physicell_id,
    parent_ID = ifelse(
      is.na(node_parent_index),
      NA_character_,
      simulation$nodes$physicell_id[node_parent_index]
    ),
    x = 0,
    y = 0,
    z = 0,
    neighbor_IDs = '',
    cell_type = terminal_nodes$cell_type,
    developmental_pseudotime = pseudotime,
    generation = terminal_nodes$generation,
    birth_time = terminal_nodes$birth_time,
    end_time = terminal_nodes$end_time,
    stringsAsFactors = FALSE
  )
}
