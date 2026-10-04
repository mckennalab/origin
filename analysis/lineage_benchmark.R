# Simulation-only benchmark for comparing integrated and mitochondrial lineage
# recorders.
#
# Population trees are generated once at the largest requested size. Recorder
# profiles are generated once at the largest requested integration count or
# mitochondrial observation depth. Nested cell and recording subsets then
# produce paired reconstruction inputs without rerunning biology or mutation.

# ---- Validation and logging helpers ----

#' Coerce a benchmark grid argument to a sorted unique integer vector
#'
#' Accepts a numeric vector, a list of numbers, or a single comma-delimited
#' string such as `'1,2,5'`. Every element must be finite, whole, at least
#' `minimum`, and distinct; duplicates are rejected rather than collapsed.
#'
#' @param value Numeric vector, list of numbers, or one comma-delimited
#'   character scalar.
#' @param name Argument label used in the error message.
#' @param minimum Smallest permitted value; defaults to `1L`.
#' @return An integer vector sorted in increasing order.
lineage_benchmark_integer_vector <- function(value, name, minimum = 1L){
  if(length(value) == 1L && is.character(value)){
    value <- trimws(unlist(strsplit(value, ',', fixed = TRUE)))
  }
  value <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
  if(length(value) == 0L || any(!is.finite(value)) ||
     any(value < minimum) || any(value %% 1 != 0) || anyDuplicated(value)){
    stop(sprintf('%s must contain unique integers >= %d.', name, minimum))
  }
  sort(as.integer(value))
}

#' Print one timestamped benchmark progress line
#'
#' @param ... Message components pasted together without a separator.
#' @param enabled When not `TRUE` the call is a no-op; defaults to `TRUE`.
#' @return `NULL`, invisibly.
#' @section Side effects: Writes one line to standard output and flushes the
#'   console.
lineage_benchmark_log <- function(..., enabled = TRUE){
  if(!isTRUE(enabled)){
    return(invisible(NULL))
  }
  cat(sprintf(
    '[%s] %s\n',
    format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z'),
    paste0(..., collapse = '')
  ))
  flush.console()
  invisible(NULL)
}

# ---- Population parameter construction ----

#' Build one native editing-state block for the benchmark
#'
#' Every background process is switched off, so the only mutation sources are
#' the base-editing rate and the mitochondrial substitution rate supplied here.
#' Both substitution models are Jukes-Cantor; barcode substitution, and
#' insertion and deletion of both barcode and mitochondrial sequence, are held
#' at zero.
#'
#' @param be_probability Base-edit probability per target per division.
#' @param mt_probability Mitochondrial substitution probability per base per
#'   division under the `'JC'` model; defaults to `0`.
#' @return A named list of native editing parameters covering the
#'   mitochondrial, barcode, base-editing, and nuclease channels.
lineage_benchmark_editing_state <- function(be_probability,
                                            mt_probability = 0){
  list(
    mt_substitution_model = 'JC',
    mt_sub_model_params = list(mt_probability),
    mt_bg_insertion_prob_per_division = 0,
    mt_bg_deletion_prob_per_division = 0,
    bc_substitution_model = 'JC',
    bc_sub_model_params = list(0),
    bc_bg_insertion_prob_per_division = 0,
    bc_bg_deletion_prob_per_division = 0,
    be_mutations_per_target_per_division = be_probability,
    nuc_insertions_per_target_per_division = 0,
    nuc_deletions_per_target_per_division = 0
  )
}

#' Build one benchmark cell-type parameter block
#'
#' The uninduced state records nothing and the induced state records at
#' `be_probability`; no invariant barcode or mitochondrial sites are reserved.
#'
#' @param cell_cycle_length Mean cell-cycle length in simulation time units.
#' @param death_probability Death probability per cell cycle; defaults to `0`.
#' @param be_probability Induced base-edit probability per target per division;
#'   defaults to `0.00952`.
#' @return A named list with `cell_cycle_length`,
#'   `death_per_cell_cycle_prob`, zero invariant-site counts, and the
#'   `uninduced_editing_params` and `induced_editing_params` blocks.
lineage_benchmark_cell_type <- function(cell_cycle_length,
                                        death_probability = 0,
                                        be_probability = 0.00952){
  list(
    cell_cycle_length = cell_cycle_length,
    death_per_cell_cycle_prob = death_probability,
    bc_invariant_sites = 0,
    mt_invariant_sites = 0,
    uninduced_editing_params = lineage_benchmark_editing_state(0),
    induced_editing_params = lineage_benchmark_editing_state(be_probability)
  )
}

#' Assemble native population parameters for one benchmark topology
#'
#' `balanced` and `comb` share a single cycling type and exist so the synthetic
#' division generator can produce the two topology extremes; `neutral` is that
#' same cycling type grown by the Gillespie engine; `turnover` adds a death
#' probability of 0.3 per cell cycle; `hierarchical` uses stem, progenitor, and
#' terminal types with a fixed transition matrix and an effectively
#' non-dividing terminal type (cell-cycle length 1e9). Editing and
#' differentiation are both induced at time 0.
#'
#' @param shape One of `'balanced'`, `'comb'`, `'neutral'`, `'turnover'`, or
#'   `'hierarchical'`; matched case-insensitively, anything else is an error.
#' @param maximum_cells Cell cap for the population; the first validated
#'   integer is used.
#' @param seed Random seed stored as `random_seed`; may be `0` or greater.
#' @return A named list of native population parameters for the Gillespie
#'   engine, including a `benchmark` element carrying the resolved `shape` and
#'   `maximum_cells`.
lineage_benchmark_population_params <- function(shape,
                                                maximum_cells,
                                                seed){
  shape <- tolower(as.character(shape))
  maximum_cells <- lineage_benchmark_integer_vector(
    maximum_cells,
    'maximum_cells'
  )[1]
  seed <- lineage_benchmark_integer_vector(seed, 'seed', minimum = 0L)[1]
  if(shape %in% c('balanced', 'comb')){
    cell_types <- list(
      cycling = lineage_benchmark_cell_type(1, 0)
    )
    transition_matrix <- list(list(1))
    founder_type <- 'cycling'
  } else if(shape == 'neutral'){
    cell_types <- list(
      cycling = lineage_benchmark_cell_type(1, 0)
    )
    transition_matrix <- list(list(1))
    founder_type <- 'cycling'
  } else if(shape == 'turnover'){
    cell_types <- list(
      cycling = lineage_benchmark_cell_type(1, 0.3)
    )
    transition_matrix <- list(list(1))
    founder_type <- 'cycling'
  } else if(shape == 'hierarchical'){
    cell_types <- list(
      stem = lineage_benchmark_cell_type(1.5, 0),
      progenitor = lineage_benchmark_cell_type(0.75, 0),
      terminal = lineage_benchmark_cell_type(1e9, 0)
    )
    transition_matrix <- list(
      list(0.75, 0.25, 0),
      list(0, 0.70, 0.30),
      list(0, 0, 1)
    )
    founder_type <- 'stem'
  } else{
    stop(sprintf(
      'Unknown benchmark shape %s; use balanced, comb, neutral, turnover, or hierarchical.',
      shape
    ))
  }
  list(
    simulation_engine = 'gillespie',
    num_init_cells = 1L,
    sim_length = list(1000),
    random_seed = seed,
    consider_cell_heteroplasmy_scores = FALSE,
    editing_induction = list(timepoint = 0),
    differentiation_induction = list(timepoint = 0),
    cell_type_dict = list(
      founder_cell_type = founder_type,
      cell_type_params = cell_types,
      uninduced_transition_matrix = transition_matrix,
      induced_transition_matrix = transition_matrix
    ),
    benchmark = list(shape = shape, maximum_cells = maximum_cells)
  )
}

# ---- Synthetic topology generation ----

#' Generate a synthetic division table with an extreme tree topology
#'
#' Emits exactly `terminal_cells - 1` divisions in PhysiCell persistent-ID
#' form: the parent keeps its ID and each division mints one new integer
#' daughter ID. `'comb'` divides cell `1` repeatedly, giving a fully asymmetric
#' caterpillar with evenly spaced division times; `'balanced'` divides in
#' breadth-first queue order, so terminal depths differ by at most one and
#' division times are set by generation level scaled across `duration`.
#'
#' @param shape Either `'balanced'` or `'comb'`; any other shape is an error.
#' @param terminal_cells Exact number of live cells to end with; `1` yields a
#'   zero-row table.
#' @param duration Positive total simulated time spanned by the division
#'   times; defaults to `10`.
#' @return A data frame with `time`, `parent_ID`, and `daughter_ID` columns,
#'   one row per division.
lineage_benchmark_synthetic_divisions <- function(shape,
                                                  terminal_cells,
                                                  duration = 10){
  shape <- tolower(as.character(shape))
  if(!(shape %in% c('balanced', 'comb'))){
    stop('Synthetic benchmark divisions require balanced or comb shape.')
  }
  terminal_cells <- lineage_benchmark_integer_vector(
    terminal_cells,
    'terminal_cells'
  )[1]
  duration <- suppressWarnings(as.numeric(duration))
  if(length(duration) != 1L || !is.finite(duration) || duration <= 0){
    stop('Synthetic benchmark duration must be one positive number.')
  }
  if(terminal_cells == 1L){
    return(data.frame(
      time = numeric(), parent_ID = character(), daughter_ID = character(),
      stringsAsFactors = FALSE
    ))
  }
  number_divisions <- terminal_cells - 1L
  parent_ids <- daughter_ids <- character(number_divisions)
  event_levels <- integer(number_divisions)
  next_cell_id <- 1L
  if(shape == 'comb'){
    for(event_index in seq_len(number_divisions)){
      next_cell_id <- next_cell_id + 1L
      parent_ids[event_index] <- '1'
      daughter_ids[event_index] <- as.character(next_cell_id)
      event_levels[event_index] <- event_index
    }
    event_times <- seq_len(number_divisions) * duration / terminal_cells
  } else{
    queue_ids <- integer(2L * terminal_cells)
    queue_depths <- integer(2L * terminal_cells)
    queue_ids[1] <- 1L
    queue_start <- 1L
    queue_end <- 1L
    for(event_index in seq_len(number_divisions)){
      parent_id <- queue_ids[queue_start]
      parent_depth <- queue_depths[queue_start]
      queue_start <- queue_start + 1L
      next_cell_id <- next_cell_id + 1L
      parent_ids[event_index] <- as.character(parent_id)
      daughter_ids[event_index] <- as.character(next_cell_id)
      event_levels[event_index] <- parent_depth + 1L
      queue_end <- queue_end + 1L
      queue_ids[queue_end] <- parent_id
      queue_depths[queue_end] <- parent_depth + 1L
      queue_end <- queue_end + 1L
      queue_ids[queue_end] <- next_cell_id
      queue_depths[queue_end] <- parent_depth + 1L
    }
    event_times <- event_levels * duration / (max(event_levels) + 1L)
  }
  data.frame(
    time = event_times,
    parent_ID = parent_ids,
    daughter_ID = daughter_ids,
    stringsAsFactors = FALSE
  )
}

#' Add the node fields the recorders expect to a synthetic lineage
#'
#' Generation is filled in by reading each node's parent generation, which
#' assumes the node table lists every parent before its daughters. All nodes
#' are marked induced and differentiated, terminal nodes are the survivors, and
#' nothing dies, which matches the deterministic synthetic topologies.
#'
#' @param nodes Event-resolved node table from `build_physicell_lineage()`,
#'   carrying `node_id`, `parent_node_id`, and `is_terminal` columns.
#' @param cell_type Type label written to every node; defaults to `'cycling'`.
#' @return `nodes` with added `cell_type`, `generation` (0 for founders),
#'   `editing_state`, `differentiation_induced`, `alive_at_end`, and `died`
#'   columns.
lineage_benchmark_annotate_nodes <- function(nodes, cell_type = 'cycling'){
  parent_indices <- match(nodes$parent_node_id, nodes$node_id)
  generation <- integer(nrow(nodes))
  if(nrow(nodes) > 1L){
    for(node_index in seq_len(nrow(nodes))){
      if(!is.na(parent_indices[node_index])){
        generation[node_index] <- generation[parent_indices[node_index]] + 1L
      }
    }
  }
  nodes$cell_type <- cell_type
  nodes$generation <- generation
  nodes$editing_state <- 'induced'
  nodes$differentiation_induced <- TRUE
  nodes$alive_at_end <- nodes$is_terminal
  nodes$died <- FALSE
  nodes
}

# ---- Population simulation and nested sampling ----

#' Simulate one benchmark population with an exact live-cell count
#'
#' `'balanced'` and `'comb'` are built directly from a synthetic division table
#' so their topology is exact; every other shape is grown with
#' `simulate_gillespie_population()` capped at `terminal_cells`. Either way the
#' function errors unless the surviving terminal count equals `terminal_cells`
#' exactly, so the nested sample panel always has the full population to draw
#' from.
#'
#' @param shape Benchmark topology name accepted by
#'   `lineage_benchmark_population_params()`.
#' @param terminal_cells Exact number of live terminal cells required.
#' @param seed Random seed for the population.
#' @param synthetic_duration Total time spanned by synthetic divisions, also
#'   used as the lineage end time; defaults to `10`.
#' @param show_progress Whether to emit progress lines; defaults to `TRUE`.
#' @return A named list with `population` (nodes, terminal nodes, division
#'   events, counts, and stop reason) and `params` (the native population
#'   parameters).
simulate_lineage_benchmark_population <- function(shape,
                                                  terminal_cells,
                                                  seed,
                                                  synthetic_duration = 10,
                                                  show_progress = TRUE){
  shape <- tolower(as.character(shape))
  params <- lineage_benchmark_population_params(shape, terminal_cells, seed)
  lineage_benchmark_log(
    sprintf('Generating %s population with %s live cells.',
            shape, format(terminal_cells, big.mark = ',')),
    enabled = show_progress
  )
  if(shape %in% c('balanced', 'comb')){
    divisions <- lineage_benchmark_synthetic_divisions(
      shape,
      terminal_cells,
      duration = synthetic_duration
    )
    nodes <- build_physicell_lineage(
      divisions,
      end_time = synthetic_duration,
      founder_ids = '1',
      show_progress = FALSE
    )
    nodes <- lineage_benchmark_annotate_nodes(nodes, 'cycling')
    terminal_nodes <- nodes[nodes$is_terminal, , drop = FALSE]
    population <- list(
      nodes = nodes,
      terminal_nodes = terminal_nodes,
      division_events = divisions,
      event_log = data.frame(),
      checkpoint_summary = data.frame(),
      end_time = synthetic_duration,
      requested_end_time = synthetic_duration,
      stop_reason = 'terminal_cell_target',
      seed = seed,
      counts = c(
        founders = 1L,
        divisions = nrow(divisions),
        deaths = 0L,
        nodes = nrow(nodes),
        surviving_cells = nrow(terminal_nodes)
      ),
      rates = NULL
    )
  } else{
    population <- suppressWarnings(simulate_gillespie_population(
      params,
      end_time = 1000,
      seed = seed,
      max_cells = terminal_cells,
      show_progress = FALSE
    ))
  }
  if(nrow(population$terminal_nodes) != terminal_cells){
    stop(sprintf(
      '%s population produced %d live cells instead of %d (stop reason: %s).',
      shape,
      nrow(population$terminal_nodes),
      terminal_cells,
      population$stop_reason
    ))
  }
  list(population = population, params = params)
}

#' Draw nested terminal-cell samples for every requested tree size
#'
#' The terminal IDs are shuffled once under `seed` and each requested size is
#' the leading prefix of that single order, so a smaller sample is always a
#' strict subset of every larger one and the whole panel is reproducible.
#'
#' @param population Population list carrying a `terminal_nodes` table with a
#'   `physicell_id` column.
#' @param tree_sizes Requested sample sizes; none may exceed the live terminal
#'   population.
#' @param seed Integer seed for the shuffle.
#' @return A list of character vectors of PhysiCell IDs, named by tree size.
#' @section Side effects: Calls `set.seed()`, which resets the global RNG
#'   stream.
lineage_benchmark_sample_sets <- function(population,
                                          tree_sizes,
                                          seed){
  tree_sizes <- lineage_benchmark_integer_vector(tree_sizes, 'tree_sizes')
  terminal_ids <- as.character(population$terminal_nodes$physicell_id)
  if(max(tree_sizes) > length(terminal_ids)){
    stop('A requested tree size exceeds the live terminal population.')
  }
  set.seed(as.integer(seed))
  terminal_order <- sample(terminal_ids, length(terminal_ids), replace = FALSE)
  sets <- lapply(tree_sizes, function(number){
    terminal_order[seq_len(number)]
  })
  names(sets) <- as.character(tree_sizes)
  sets
}

# ---- Recorder chemistry presets ----

#' Build a valid zero-target native target specification
#'
#' Switches off the base-editing or nuclease channel of a recorder while
#' keeping the parameter block structurally complete.
#'
#' @return A named list target specification whose `num_targets` is `NULL` and
#'   whose editing window has size `0`.
lineage_benchmark_empty_target_spec <- function(){
  list(
    num_targets = NULL,
    edit_rate_class_fractions = list(high = 1, medium = 0, low = 0),
    config = 'U',
    editing_window = list(size = 0L, decaying = FALSE, close_after_edit = TRUE),
    prime_editing_system = FALSE
  )
}

#' Build the PEtracer mark pool used by the benchmark panel
#'
#' PEtracer (doi:10.1126/science.adx3800) places three edit sites per cassette
#' and gives each site its own alphabet of eight predefined 5nt marks, for 24
#' distinct marks in total. The alphabet is what keeps independent edits at one
#' site distinguishable: with a single outcome per site, two cells that edit it
#' separately look identical to two cells sharing an ancestor, and that
#' homoplasy misleads reconstruction.
#'
#' Marks are emitted in per-site blocks of eight, which is the layout
#' `expand_prime_editing_marks()` expects. Within a site every mark carries the
#' same efficiency, so the installed mark is uniform over the eight; the
#' efficiency differs BETWEEN sites, taking the in vitro modification rates the
#' paper reports for its three edit sites (73.4%, 47.7%, 86.5%).
#'
#' @return A list of 24 pegRNA definitions, ordered site-major.
lineage_benchmark_petracer_pool <- function(){
  site_efficiency <- c(0.734, 0.477, 0.865)
  bases <- c('A', 'C', 'G', 'T')
  pool <- list()
  for(site in seq_along(site_efficiency)){
    for(mark in 1:8){
      index <- (site - 1L) * 8L + mark
      # Base-4 enumeration of the mark index, so all 24 sequences are distinct
      # by construction. An ad-hoc arithmetic mix collides: marks that share a
      # sequence are the same character, which silently shrinks the alphabet.
      digits <- (index - 1L) %/% 4L^(4:0) %% 4L
      sequence <- paste(bases[digits + 1L], collapse = '')
      pool[[index]] <- list(
        pegRNA_id = sprintf('petracer_site%d_mark%02d', site, mark),
        edit_sequence = sequence,
        editing_efficiency = site_efficiency[site],
        description = sprintf(
          'PEtracer site %d mark %d (site efficiency %.3f)',
          site, mark, site_efficiency[site]
        )
      )
    }
  }
  pool
}

#' Return the fixed benchmark pegRNA pool
#'
#' @return A list of six pegRNA definitions, each with `pegRNA_id`, a
#'   five-base `edit_sequence`, an `editing_efficiency` stepping 0.90, 0.75,
#'   0.60, 0.45, 0.30, 0.15, and a `description`.
lineage_benchmark_prime_pool <- function(){
  list(
    list(
      pegRNA_id = 'benchmark_peg_01', edit_sequence = 'AGGCT',
      editing_efficiency = 0.90,
      description = 'benchmark high-efficiency pegRNA'
    ),
    list(
      pegRNA_id = 'benchmark_peg_02', edit_sequence = 'CTTGA',
      editing_efficiency = 0.75,
      description = 'benchmark high-medium-efficiency pegRNA'
    ),
    list(
      pegRNA_id = 'benchmark_peg_03', edit_sequence = 'GACCA',
      editing_efficiency = 0.60,
      description = 'benchmark medium-efficiency pegRNA'
    ),
    list(
      pegRNA_id = 'benchmark_peg_04', edit_sequence = 'TGGAT',
      editing_efficiency = 0.45,
      description = 'benchmark medium-low-efficiency pegRNA'
    ),
    list(
      pegRNA_id = 'benchmark_peg_05', edit_sequence = 'ACGGT',
      editing_efficiency = 0.30,
      description = 'benchmark low-efficiency pegRNA'
    ),
    list(
      pegRNA_id = 'benchmark_peg_06', edit_sequence = 'GTCCA',
      editing_efficiency = 0.15,
      description = 'benchmark rare-efficiency pegRNA'
    )
  )
}

#' Report how many logical targets one integration carries
#'
#' @param system Recorder name: `'baseline'`, `'prime'`, `'palincode'`, or
#'   `'mitochondrial'`; matched case-insensitively.
#' @return `50L`, `6L`, or `2L` for the integrated recorders, and
#'   `NA_integer_` for `'mitochondrial'`, whose observation unit is sampled
#'   genome depth rather than an integration count. Any other name is an
#'   error.
lineage_benchmark_targets_per_integration <- function(system){
  switch(
    tolower(as.character(system)),
    baseline = 50L,
    prime = 6L,
    palincode = 2L,
    mitochondrial = NA_integer_,
    stop(sprintf('Unknown benchmark recorder system: %s.', system))
  )
}

#' Derive mitochondrial recording parameters from population parameters
#'
#' Adds a 16,569-base mitochondrial reference with no baseline heteroplasmy and
#' rewrites every cell type so the uninduced state accumulates nothing and the
#' induced state substitutes at `substitution_probability` per base per
#' division under a Jukes-Cantor model. Mitochondrial insertions and deletions
#' stay at zero, profiles are stored sparsely, and internal-node profiles are
#' discarded.
#'
#' @param population_params Native population parameters to extend.
#' @param substitution_probability Induced mitochondrial substitution
#'   probability per base per division, which must lie in [0, 1]; defaults to
#'   `2e-6`.
#' @return The population parameters with mitochondrial genome, heteroplasmy,
#'   adapter, and per-cell-type editing fields set.
lineage_benchmark_mito_params <- function(population_params,
                                          substitution_probability = 2e-6){
  substitution_probability <- suppressWarnings(
    as.numeric(substitution_probability)
  )
  if(length(substitution_probability) != 1L ||
     !is.finite(substitution_probability) ||
     substitution_probability < 0 || substitution_probability > 1){
    stop('Mitochondrial substitution probability must be in [0, 1].')
  }
  params <- population_params
  params$mito_genome_length <- 16569L
  params$baseline_heteroplasmy_sites_frac <- 0
  params$baseline_heteroplasmy_variant_frac_dist <- list(1, 9)
  params$heteroplasmy_variant_transition_prob <- 2 / 3
  params$physicell_adapter <- list(
    recorder_system = 'mitochondrial',
    profile_storage = 'sparse',
    compact_output = TRUE,
    retain_internal_profiles = FALSE
  )
  for(cell_type in names(params$cell_type_dict$cell_type_params)){
    cell_params <- params$cell_type_dict$cell_type_params[[cell_type]]
    cell_params$uninduced_editing_params$mt_substitution_model <- 'JC'
    cell_params$uninduced_editing_params$mt_sub_model_params <- list(0)
    cell_params$uninduced_editing_params$mt_bg_insertion_prob_per_division <- 0
    cell_params$uninduced_editing_params$mt_bg_deletion_prob_per_division <- 0
    cell_params$induced_editing_params$mt_substitution_model <- 'JC'
    cell_params$induced_editing_params$mt_sub_model_params <- list(
      substitution_probability
    )
    cell_params$induced_editing_params$mt_bg_insertion_prob_per_division <- 0
    cell_params$induced_editing_params$mt_bg_deletion_prob_per_division <- 0
    params$cell_type_dict$cell_type_params[[cell_type]] <- cell_params
  }
  params
}

#' Assemble native recording parameters for one integrated recorder preset
#'
#' `'baseline'` is a 1,435-base barcode with 50 base-editing targets in a
#' decaying seven-base window and an `A --> G` conversion; `'prime'` is an
#' 80-base barcode with six prime-editing targets driven by the fixed pegRNA
#' pool at 0.08 induced edits per cell cycle; `'palincode'` is a two-base
#' barcode with two cbits (`PalT7`, `PalRNF2`) and their left, right, and both
#' outcome fractions. All three use uniform barcode nucleotide composition,
#' sparse profile storage, and compact output.
#'
#' @param population_params Native population parameters to extend.
#' @param system One of `'baseline'`, `'prime'`, or `'palincode'`;
#'   mitochondrial recording uses `lineage_benchmark_mito_params()` instead.
#' @param maximum_integrations Barcode integrations per cell; the first
#'   validated integer is used.
#' @return The population parameters extended with barcode composition and
#'   length, integration count, base-editing and nuclease target
#'   specifications, and the adapter block for the chosen system.
lineage_benchmark_recording_params <- function(population_params,
                                               system,
                                               maximum_integrations){
  system <- tolower(as.character(system))
  if(!(system %in% c('baseline', 'prime', 'palincode'))){
    stop('Recorder system must be baseline, prime, or palincode.')
  }
  maximum_integrations <- lineage_benchmark_integer_vector(
    maximum_integrations,
    'maximum_integrations'
  )[1]
  params <- population_params
  params$bc_nuc_composition <- list(
    frac_a = 0.25, frac_g = 0.25, frac_c = 0.25, frac_t = 0.25
  )
  params$barcode_sequence <- NULL
  params$max_bc_ints_per_cell <- list(maximum_integrations)
  params$be_conversion_pattern <- 'A --> G'
  params$physicell_adapter <- list(
    recorder_system = system,
    profile_storage = 'sparse',
    compact_output = TRUE,
    retain_internal_profiles = FALSE
  )
  params$prime_editing_backend <- NULL
  params$palincode_adapter <- NULL

  if(system == 'baseline'){
    params$bc_length <- 1435L
    params$physicell_adapter$recorder_system <-
      'BASELINE-like hyperdCas12a-ABE8e'
    params$be_targets <- list(
      num_targets = 50L,
      edit_rate_class_fractions = list(high = 0.4, medium = 0, low = 0.6),
      config = 'S:10:24',
      editing_window = list(
        size = 7L, decaying = TRUE, close_after_edit = FALSE
      )
    )
    params$nuclease_targets <- lineage_benchmark_empty_target_spec()
  } else if(system == 'prime'){
    params$bc_length <- 80L
    params$physicell_adapter$recorder_system <- 'prime editing'
    params$be_targets <- lineage_benchmark_empty_target_spec()
    # PEtracer geometry: three edit sites per cassette, each with its own
    # alphabet of eight marks. Targets are NOT bound to a single pegRNA here --
    # target_pegRNA_ids would pin one outcome per site and reintroduce the
    # homoplasy the alphabet exists to avoid -- so the pool is blocked instead
    # and marks_per_target carves it into per-site alphabets.
    params$nuclease_targets <- list(
      num_targets = 3L,
      edit_rate_class_fractions = list(high = 1, medium = 0, low = 0),
      config = 'S:4:10',
      editing_window = list(
        size = 0L, decaying = FALSE, close_after_edit = TRUE
      ),
      prime_editing_system = TRUE
    )
    params$prime_editing_backend <- list(
      enabled = TRUE,
      pegRNAs = lineage_benchmark_petracer_pool(),
      assignment = 'cycle',
      marks_per_target = 8L,
      induced_edit_probability_per_cell_cycle = 0.18,
      uninduced_edit_probability_per_cell_cycle = 0,
      static_id_length = 12L
    )
  } else{
    params$bc_length <- 2L
    params$physicell_adapter$recorder_system <- 'PALINCODE'
    params$be_targets <- lineage_benchmark_empty_target_spec()
    params$nuclease_targets <- lineage_benchmark_empty_target_spec()
    params$palincode_adapter <- list(
      enabled = TRUE,
      num_cbits_per_integration = 2L,
      cbit_names = list('PalT7', 'PalRNF2'),
      static_id_length = 12L,
      uninduced_edit_probability_per_cbit_per_cell_cycle = 0,
      induced_edit_probability_per_cbit_per_cell_cycle = list(0.08, 0.05),
      left_edit_fraction = list(0.495, 0.35),
      right_edit_fraction = list(0.495, 0.64),
      both_edit_fraction = list(0.01, 0.01)
    )
  }
  params
}

# ---- Recording matrix subsetting and logical targets ----

#' Parse the integration number out of recording-matrix column names
#'
#' @param column_names Character vector of native matrix column names, every
#'   one of which must begin with `int_<number>_`.
#' @return An integer vector of integration numbers, one per column.
lineage_benchmark_column_integrations <- function(column_names){
  matches <- regexec('^int_([0-9]+)_', column_names)
  parts <- regmatches(column_names, matches)
  integrations <- suppressWarnings(as.integer(vapply(parts, function(value){
    if(length(value) < 2L) NA_character_ else value[2L]
  }, character(1))))
  if(anyNA(integrations)){
    stop('A recording matrix column does not begin with int_<number>_.')
  }
  integrations
}

#' Take a paired cell and integration submatrix of a recording matrix
#'
#' Rows come back in `sample_ids` order, and only columns whose parsed
#' integration number is at or below `integrations` are kept. This is what
#' makes every nested integration panel a subset of the single maximum
#' recording rather than an independent simulation.
#'
#' @param matrix_value Sparse matrix whose rows are named `cell_<physicell_id>`
#'   and whose columns are named `int_<number>_...`.
#' @param sample_ids Ordered sample row names; every one must be present.
#' @param integrations Highest integration number to retain.
#' @return The row- and column-subset matrix, with dimensions preserved even
#'   for a single row or column.
lineage_benchmark_subset_matrix <- function(matrix_value,
                                            sample_ids,
                                            integrations){
  row_indices <- match(sample_ids, rownames(matrix_value))
  if(anyNA(row_indices)){
    stop('A benchmark sample is absent from the full recording matrix.')
  }
  column_integrations <- lineage_benchmark_column_integrations(
    colnames(matrix_value)
  )
  retained_columns <- which(column_integrations <= integrations)
  matrix_value[row_indices, retained_columns, drop = FALSE]
}

#' Take a cell subset of a matrix, keeping every feature column
#'
#' @param matrix_value Matrix whose rows are named by sample ID.
#' @param sample_ids Ordered sample row names; every one must be present.
#' @return The row-subset matrix in `sample_ids` order with all columns
#'   retained.
lineage_benchmark_subset_rows <- function(matrix_value, sample_ids){
  row_indices <- match(sample_ids, rownames(matrix_value))
  if(anyNA(row_indices)){
    stop('A benchmark sample is absent from the full recording matrix.')
  }
  matrix_value[row_indices, , drop = FALSE]
}

#' Collapse a native state matrix to one column per logical target
#'
#' PALINCODE and prime-editing recorders already carry one state per physical
#' target, so the matching `int_<n>_pos_<p>` columns are selected and renamed.
#' A BASELINE target instead spans an editing window of several barcode
#' positions, so each window is collapsed into one categorical value: an edited
#' position contributes bit `2^(j - 1)` for its rank `j` within the window,
#' giving `0` for an unedited target and a distinct nonzero integer for each
#' combination of edited positions. Window membership is taken from
#' `model$windows`, falling back to the target position itself when a window is
#' empty.
#'
#' @param state_matrix Sparse native state matrix with `int_<n>_pos_<p>`
#'   columns and cells as rows.
#' @param model Prepared recording model; `is_palincode` and
#'   `is_prime_editing` select the pass-through branches, otherwise
#'   `be_targets` and `windows` drive the bitmask.
#' @param integrations Number of integrations to emit columns for; targets are
#'   emitted for integrations `1` through this value.
#' @return A sparse matrix with one column per integration and target, named
#'   `int_<n>_target_<i>`, rows in `state_matrix` order.
#' @note Requires the Matrix package, and rejects any window wider than 52
#'   positions because the bitmask would exceed exact double precision.
#' Build the character matrix tree building will read
#'
#' The dispatch point for recorder-specific character encodings. Tree building
#' consumes `character_matrix`, so whatever this returns is what a reconstruction
#' method actually sees -- which makes it the place a recorder either keeps its
#' allele alphabet or loses it. A WT-Cas9 array is the case in point: its
#' deletions carry distinct breakpoints, and returning binary scores collapses
#' every one of them to "edited", leaving independent events indistinguishable
#' from shared ancestry.
#'
#' The default is the prime-editing rule, which returns the mark identities when
#' a recorder has a mark alphabet and binary scores otherwise. Recorders injected
#' by a harness override this.
#'
#' @param raw_alleles Cell-by-position allele matrix.
#' @param binary_scores Cell-by-position binary score matrix.
#' @param model Prepared recording model.
#' @param integrations Integrations to encode.
#' @return A cell-by-character matrix.
lineage_benchmark_character_matrix <- function(raw_alleles,
                                               binary_scores,
                                               model,
                                               integrations){
  hook <- origin_recorder_hook(model, 'character_matrix')
  if(!is.null(hook)){
    return(hook(raw_alleles, binary_scores, model, integrations))
  }
  physicell_prime_character_matrix(raw_alleles, binary_scores, model)
}

lineage_benchmark_logical_target_matrix <- function(state_matrix,
                                                    model,
                                                    integrations){
  if(!requireNamespace('Matrix', quietly = TRUE)){
    stop('The Matrix package is required for logical target matrices.')
  }
  integrations <- as.integer(integrations)
  if(isTRUE(model$is_palincode)){
    target_positions <- model$palincode$cbit_positions
    target_names <- model$palincode$cbit_names
    retained <- unlist(lapply(seq_len(integrations), function(integration){
      match(
        paste0('int_', integration, '_pos_', target_positions),
        colnames(state_matrix)
      )
    }))
    if(anyNA(retained)){
      stop('PALINCODE logical targets are absent from the state matrix.')
    }
    result <- state_matrix[, retained, drop = FALSE]
    colnames(result) <- unlist(lapply(seq_len(integrations), function(integration){
      paste0('int_', integration, '_target_', seq_along(target_names))
    }))
    return(result)
  }
  if(isTRUE(model$is_prime_editing)){
    target_positions <- model$prime_editing$targets$target_position
    retained <- unlist(lapply(seq_len(integrations), function(integration){
      match(
        paste0('int_', integration, '_pos_', target_positions),
        colnames(state_matrix)
      )
    }))
    if(anyNA(retained)){
      stop('Prime-editing logical targets are absent from the state matrix.')
    }
    result <- state_matrix[, retained, drop = FALSE]
    colnames(result) <- unlist(lapply(seq_len(integrations), function(integration){
      paste0('int_', integration, '_target_', seq_along(target_positions))
    }))
    return(result)
  }

  target_positions <- sort(as.integer(names(model$be_targets)))
  original_target_positions <- as.integer(names(model$be_targets))
  num_targets <- length(target_positions)
  num_rows <- nrow(state_matrix)
  num_columns <- integrations * num_targets
  row_chunks <- vector('list', num_columns)
  value_chunks <- vector('list', num_columns)
  column_names <- character(num_columns)
  output_column <- 0L
  for(integration in seq_len(integrations)){
    for(target_index in seq_len(num_targets)){
      output_column <- output_column + 1L
      window_index <- match(
        target_positions[target_index],
        original_target_positions
      )
      window_name <- paste0('be_window_', window_index)
      positions <- model$windows[[window_name]]
      if(is.null(positions) || length(positions) == 0L){
        positions <- target_positions[target_index]
      }
      position_columns <- match(
        paste0('int_', integration, '_pos_', positions),
        colnames(state_matrix)
      )
      position_columns <- position_columns[!is.na(position_columns)]
      if(length(position_columns) == 0L){
        stop(sprintf(
          'BASELINE target %d integration %d has no state-matrix positions.',
          target_index,
          integration
        ))
      }
      if(length(position_columns) > 52L){
        stop('A BASELINE logical target exceeds exact double bitmask capacity.')
      }
      edited <- state_matrix[, position_columns, drop = FALSE]
      if(length(edited@x) > 0L){
        edited@x[] <- 1
      }
      weights <- 2 ^ (seq_along(position_columns) - 1L)
      states <- as.numeric(edited %*% weights)
      nonzero <- which(states != 0)
      row_chunks[[output_column]] <- nonzero
      value_chunks[[output_column]] <- states[nonzero]
      column_names[output_column] <- paste0(
        'int_', integration, '_target_', target_index
      )
    }
  }
  counts <- lengths(row_chunks)
  Matrix::sparseMatrix(
    i = as.integer(unlist(row_chunks, use.names = FALSE)),
    j = as.integer(rep.int(seq_len(num_columns), counts)),
    x = as.numeric(unlist(value_chunks, use.names = FALSE)),
    dims = c(num_rows, num_columns),
    dimnames = list(rownames(state_matrix), column_names)
  )
}

#' Annotate a target layout with its logical positions and state encoding
#'
#' PALINCODE and prime-editing layouts map one logical target to one physical
#' position and are labelled `0=wild_type;1=left;2=right;3=both` and
#' `0=unedited;1=assigned_pegRNA_edit` respectively. A BASELINE target lists
#' every position in its editing window, semicolon-separated in the order the
#' bitmask bits are assigned by `lineage_benchmark_logical_target_matrix()`.
#'
#' @param target_layout Native target layout with `position` and `target_index`
#'   columns; a zero-row layout is returned unchanged.
#' @param model Prepared recording model supplying the window definitions.
#' @return The layout with added `logical_positions` and
#'   `logical_state_encoding` character columns.
lineage_benchmark_logical_target_manifest <- function(target_layout, model){
  if(nrow(target_layout) == 0L){
    return(target_layout)
  }
  if(isTRUE(model$is_palincode)){
    target_layout$logical_positions <- as.character(target_layout$position)
    target_layout$logical_state_encoding <-
      '0=wild_type;1=left;2=right;3=both'
    return(target_layout)
  }
  if(isTRUE(model$is_prime_editing)){
    target_layout$logical_positions <- as.character(target_layout$position)
    target_layout$logical_state_encoding <-
      '0=unedited;1=assigned_pegRNA_edit'
    return(target_layout)
  }
  target_indices <- sort(unique(target_layout$target_index))
  sorted_target_positions <- target_layout$position[
    match(target_indices, target_layout$target_index)
  ]
  original_target_positions <- as.integer(names(model$be_targets))
  positions_by_target <- setNames(vapply(
    target_indices,
    function(target_index){
      primary_position <- sorted_target_positions[
        match(target_index, target_indices)
      ]
      window_index <- match(primary_position, original_target_positions)
      positions <- model$windows[[paste0('be_window_', window_index)]]
      if(is.null(positions) || length(positions) == 0L){
        positions <- target_layout$position[
          target_layout$target_index == target_index
        ][1]
      }
      paste(positions, collapse = ';')
    },
    character(1)
  ), as.character(target_indices))
  target_layout$logical_positions <- unname(
    positions_by_target[as.character(target_layout$target_index)]
  )
  target_layout$logical_state_encoding <- paste(
    'nonzero integer bitmask over logical_positions in listed order;',
    '0=unedited'
  )
  target_layout
}

# ---- Population and sample output ----

#' Write a value as pretty reproducibility JSON
#'
#' @param value Any jsonlite-serializable value; length-one elements are
#'   unboxed and full numeric precision is kept.
#' @param path Destination file path.
#' @return The normalized path of the written file, invisibly.
#' @section Side effects: Writes `path`.
lineage_benchmark_write_json <- function(value, path){
  if(!requireNamespace('jsonlite', quietly = TRUE)){
    stop('The jsonlite package is required for benchmark parameter output.')
  }
  jsonlite::write_json(
    value,
    path,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = 'null',
    digits = NA
  )
  invisible(normalizePath(path, mustWork = TRUE))
}

#' Write the reusable population bundle for one shape and seed
#'
#' The saved `population.rds` is what a resumed run reloads instead of
#' regenerating the population.
#'
#' @param population_bundle List with `population` and `params`, as returned by
#'   `simulate_lineage_benchmark_population()`.
#' @param output_dir Destination directory, created if it does not exist.
#' @param show_progress Whether to emit progress lines; defaults to `TRUE`.
#' @return `output_dir`, invisibly.
#' @section Side effects: Creates `output_dir` and writes `population.rds`,
#'   `population_params.json`, the compressed lineage tables and Newick trees
#'   produced by `write_physicell_lineage_outputs()`, and
#'   `division_events.csv.gz`.
write_lineage_benchmark_population <- function(population_bundle,
                                               output_dir,
                                               show_progress = TRUE){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(population_bundle, file.path(output_dir, 'population.rds'))
  lineage_benchmark_write_json(
    population_bundle$params,
    file.path(output_dir, 'population_params.json')
  )
  population <- population_bundle$population
  write_physicell_lineage_outputs(
    population$nodes,
    population$terminal_nodes,
    output_dir,
    show_progress = show_progress,
    compress_csv = TRUE
  )
  write_physicell_csv(
    population$division_events,
    file.path(output_dir, 'division_events.csv'),
    row.names = FALSE,
    compress = TRUE
  )
  invisible(output_dir)
}

#' Write the sampled-cell table and ground-truth tree for every tree size
#'
#' Each size gets its own `n_<size>` subdirectory holding the sampled terminal
#' cells in draw order (`sample_rank`) and the ground-truth Newick tree pruned
#' to exactly those cells. Sample IDs are written as `cell_<physicell_id>`,
#' which is also how the recording matrices name their rows.
#'
#' @param population Population list with `nodes` and `terminal_nodes`.
#' @param sample_sets Named list of PhysiCell ID vectors from
#'   `lineage_benchmark_sample_sets()`.
#' @param output_dir Destination directory, created if it does not exist.
#' @param show_progress Whether to emit progress lines; defaults to `TRUE`.
#' @return A character vector of normalized ground-truth tree paths, named by
#'   tree size.
#' @section Side effects: Creates `output_dir` and one `n_<size>` subdirectory
#'   per size; writes `sample_sets.rds`, `sample_cells.csv.gz`, and
#'   `ground_truth_tree.nwk`.
write_lineage_benchmark_samples <- function(population,
                                            sample_sets,
                                            output_dir,
                                            show_progress = TRUE){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sample_sets, file.path(output_dir, 'sample_sets.rds'))
  paths <- setNames(character(length(sample_sets)), names(sample_sets))
  for(size_name in names(sample_sets)){
    size_dir <- file.path(output_dir, paste0('n_', size_name))
    dir.create(size_dir, recursive = TRUE, showWarnings = FALSE)
    physicell_ids <- sample_sets[[size_name]]
    terminal_indices <- match(
      physicell_ids,
      population$terminal_nodes$physicell_id
    )
    if(anyNA(terminal_indices)){
      stop('A nested sample ID is absent from the population terminal table.')
    }
    sample_table <- population$terminal_nodes[
      terminal_indices,
      intersect(
        c('node_id', 'physicell_id', 'cell_type', 'generation',
          'birth_time', 'end_time', 'branch_length'),
        names(population$terminal_nodes)
      ),
      drop = FALSE
    ]
    sample_table$sample_id <- paste0('cell_', sample_table$physicell_id)
    sample_table$sample_rank <- seq_len(nrow(sample_table))
    sample_table <- sample_table[, c(
      'sample_id', 'sample_rank',
      setdiff(names(sample_table), c('sample_id', 'sample_rank'))
    ), drop = FALSE]
    write_physicell_csv(
      sample_table,
      file.path(size_dir, 'sample_cells.csv'),
      row.names = FALSE,
      compress = TRUE
    )
    tree_path <- file.path(size_dir, 'ground_truth_tree.nwk')
    write_physicell_lineage_newick(
      population$nodes,
      tree_path,
      terminal_physicell_ids = physicell_ids,
      show_progress = FALSE,
      progress_label = paste0('n=', size_name, ' truth tree')
    )
    paths[[size_name]] <- normalizePath(tree_path, mustWork = TRUE)
    lineage_benchmark_log(
      sprintf('Wrote nested %s-cell ground-truth tree.', size_name),
      enabled = show_progress
    )
  }
  paths
}

# ---- Mitochondrial observation panels ----

#' Draw a reproducible genome sampling order for every cell
#'
#' Each cell gets one without-replacement draw of `maximum_depth` genome
#' indices out of its biological copy number; a depth-`d` panel then reads the
#' first `d` entries of that order, so a shallower panel is a strict subset of
#' every deeper one.
#'
#' @param terminal_profiles Named list of per-cell mitochondrial profiles, each
#'   a list of genomes; the names must be present and non-empty.
#' @param maximum_depth Deepest observation depth to support; it may not exceed
#'   the smallest per-cell copy number.
#' @param seed Integer seed for the draws.
#' @return A list of integer genome-index vectors of length `maximum_depth`,
#'   one per cell, in `terminal_profiles` order.
#' @section Side effects: Calls `set.seed()`, which resets the global RNG
#'   stream.
lineage_benchmark_mito_sampling_orders <- function(terminal_profiles,
                                                   maximum_depth,
                                                   seed){
  maximum_depth <- lineage_benchmark_integer_vector(
    maximum_depth,
    'maximum_depth'
  )[1]
  if(length(terminal_profiles) == 0L || is.null(names(terminal_profiles)) ||
     any(!nzchar(names(terminal_profiles)))){
    stop('Named terminal mitochondrial profiles are required for sampling.')
  }
  profile_sizes <- lengths(terminal_profiles)
  if(any(profile_sizes < maximum_depth)){
    stop(sprintf(
      paste(
        'Mitochondrial observation depth %d exceeds the smallest biological',
        'copy number (%d).'
      ),
      maximum_depth,
      min(profile_sizes)
    ))
  }
  set.seed(as.integer(seed))
  lapply(profile_sizes, function(profile_size){
    sample.int(profile_size, maximum_depth, replace = FALSE)
  })
}

#' Build the shared variant manifest for a mitochondrial observation panel
#'
#' Scans the first `maximum_depth` sampled genomes of every cell and returns
#' the sorted union of the variants they carry, so every depth in the panel
#' shares one feature universe. Alleles use the repository encoding: `-1` is a
#' deletion, `1`-`4` are the substitutions `A`, `G`, `C`, `T`, and a decimal
#' value is an insertion whose fractional part carries the inserted base code.
#'
#' @param terminal_profiles Named list of per-cell profiles, each a list of
#'   genomes whose variants are numeric alleles named by position.
#' @param sampling_orders Genome orders from
#'   `lineage_benchmark_mito_sampling_orders()`, keyed by the same cell names.
#' @param maximum_depth Number of leading genomes per cell to scan.
#' @param reference Mitochondrial reference bases indexed by position.
#' @return A data frame with `feature` (`mt_<position>_<allele>`), `position`,
#'   `allele`, `reference`, `alternate` (a base, `-` for a deletion, or
#'   `+<base>` for an insertion), and `logical_state_encoding`
#'   (`0=not_observed;1=observed`), ordered by position then allele; zero rows
#'   when nothing was observed.
lineage_benchmark_mito_variant_manifest <- function(terminal_profiles,
                                                     sampling_orders,
                                                     maximum_depth,
                                                     reference){
  observed_positions <- list()
  observed_alleles <- list()
  observation_index <- 0L
  for(sample_id in names(terminal_profiles)){
    profile <- terminal_profiles[[sample_id]]
    selected <- sampling_orders[[sample_id]][seq_len(maximum_depth)]
    for(genome_index in selected){
      variants <- profile[[genome_index]]
      if(length(variants) == 0L){
        next
      }
      observation_index <- observation_index + 1L
      observed_positions[[observation_index]] <- as.integer(names(variants))
      observed_alleles[[observation_index]] <- as.numeric(variants)
    }
  }
  if(observation_index == 0L){
    return(data.frame(
      feature = character(), position = integer(), allele = numeric(),
      reference = character(), alternate = character(),
      logical_state_encoding = character(), stringsAsFactors = FALSE
    ))
  }
  variants <- unique(data.frame(
    position = as.integer(unlist(observed_positions, use.names = FALSE)),
    allele = as.numeric(unlist(observed_alleles, use.names = FALSE)),
    stringsAsFactors = FALSE
  ))
  variants <- variants[order(variants$position, variants$allele), , drop = FALSE]
  rownames(variants) <- NULL
  bases <- c('A', 'G', 'C', 'T')
  variants$feature <- paste0('mt_', variants$position, '_', variants$allele)
  variants$reference <- reference[variants$position]
  variants$alternate <- vapply(variants$allele, function(allele){
    if(allele == -1){
      '-'
    } else if(allele %% 1 == 0 && allele >= 1 && allele <= 4){
      bases[as.integer(allele)]
    } else{
      inserted_code <- as.integer(round(abs(allele %% 1) * 10)) %% 10
      paste0('+', bases[inserted_code])
    }
  }, character(1))
  variants$logical_state_encoding <- '0=not_observed;1=observed'
  variants[, c(
    'feature', 'position', 'allele', 'reference', 'alternate',
    'logical_state_encoding'
  ), drop = FALSE]
}

#' Build per-depth heteroplasmy and presence matrices on a shared feature set
#'
#' For each requested depth the first `depth` genomes of every cell are read
#' and a variant's value is the number of those sampled genomes carrying it
#' divided by `depth`, so `state` entries are heteroplasmy fractions in (0, 1]
#' and reach `1` only when every sampled genome carries the variant. The
#' `character` matrix is the same matrix binarised to `1` for any observation.
#' All depths use the `variant_manifest` features as columns, so matrices from
#' different depths are directly comparable.
#'
#' @param terminal_profiles Named list of per-cell mitochondrial profiles.
#' @param sampling_orders Per-cell genome orders covering the largest depth.
#' @param depths Observation depths to build, validated as unique integers.
#' @param variant_manifest Manifest whose `feature` column defines the columns.
#' @return A list keyed by depth as a character string, each element a list
#'   with `state` (heteroplasmy fractions) and `character` (binary presence)
#'   sparse matrices whose rows are the cells in `terminal_profiles` order.
#' @note Requires the Matrix package.
lineage_benchmark_mito_observation_matrices <- function(terminal_profiles,
                                                        sampling_orders,
                                                        depths,
                                                        variant_manifest){
  if(!requireNamespace('Matrix', quietly = TRUE)){
    stop('The Matrix package is required for mitochondrial observation matrices.')
  }
  depths <- lineage_benchmark_integer_vector(depths, 'depths')
  sample_ids <- names(terminal_profiles)
  features <- variant_manifest$feature
  results <- setNames(vector('list', length(depths)), as.character(depths))
  for(depth in depths){
    row_chunks <- vector('list', length(sample_ids))
    column_chunks <- vector('list', length(sample_ids))
    value_chunks <- vector('list', length(sample_ids))
    for(cell_index in seq_along(sample_ids)){
      sample_id <- sample_ids[cell_index]
      profile <- terminal_profiles[[sample_id]]
      selected <- sampling_orders[[sample_id]][seq_len(depth)]
      sampled_variants <- unlist(profile[selected], use.names = TRUE)
      if(length(sampled_variants) == 0L){
        next
      }
      sampled_features <- paste0(
        'mt_', names(sampled_variants), '_', as.numeric(sampled_variants)
      )
      feature_indices <- match(sampled_features, features)
      feature_indices <- feature_indices[!is.na(feature_indices)]
      if(length(feature_indices) == 0L){
        next
      }
      counts <- tabulate(feature_indices, nbins = length(features))
      nonzero <- which(counts > 0L)
      row_chunks[[cell_index]] <- rep.int(cell_index, length(nonzero))
      column_chunks[[cell_index]] <- nonzero
      value_chunks[[cell_index]] <- counts[nonzero] / depth
    }
    state_matrix <- Matrix::sparseMatrix(
      i = as.integer(unlist(row_chunks, use.names = FALSE)),
      j = as.integer(unlist(column_chunks, use.names = FALSE)),
      x = as.numeric(unlist(value_chunks, use.names = FALSE)),
      dims = c(length(sample_ids), length(features)),
      dimnames = list(sample_ids, features)
    )
    character_matrix <- state_matrix
    if(length(character_matrix@x) > 0L){
      character_matrix@x[] <- 1
    }
    results[[as.character(depth)]] <- list(
      state = state_matrix,
      character = character_matrix
    )
  }
  results
}

#' Simulate the mitochondrial overlay once and derive every depth panel
#'
#' Mutation is simulated at the full biological copy number
#' (`genomes_per_cell`) and observation depth is a sampling decision made
#' afterwards, so every depth in the panel shares one simulated history. The
#' sampling seed is derived deterministically from `seed`, and all cell-type
#' models must agree on a single mitochondrial reference.
#'
#' @param population_bundle Bundle with `population` and `params`.
#' @param observation_depths Depths to build panels for; the maximum may not
#'   exceed `genomes_per_cell`.
#' @param genomes_per_cell Biological mitochondrial copy number per cell.
#' @param seed Integer seed for model preparation and mutation simulation.
#' @param output_dir Destination directory, created if it does not exist.
#' @param show_progress Whether to emit progress lines; defaults to `TRUE`.
#' @return A named list with `type = 'mitochondrial'`, `model`,
#'   `observation_matrices`, `variant_manifest`, `observation_depths`,
#'   `genomes_per_cell`, and `params`.
#' @section Side effects: Creates `output_dir` and writes the exact
#'   mutation-event output of `write_physicell_mito_outputs()` plus
#'   `recording_model.rds`, `mitochondrial_sampling_orders.rds`,
#'   `mitochondrial_observation_matrices.rds`,
#'   `mitochondrial_variant_manifest.csv.gz`, `recording_params.json`,
#'   `recording_type.txt`, and the `recording_complete.txt` resume marker.
simulate_lineage_benchmark_mitochondrial <- function(population_bundle,
                                                     observation_depths,
                                                     genomes_per_cell,
                                                     seed,
                                                     output_dir,
                                                     show_progress = TRUE){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  observation_depths <- lineage_benchmark_integer_vector(
    observation_depths,
    'observation_depths'
  )
  genomes_per_cell <- lineage_benchmark_integer_vector(
    genomes_per_cell,
    'genomes_per_cell'
  )[1]
  if(max(observation_depths) > genomes_per_cell){
    stop('Mitochondrial observation depth cannot exceed genomes_per_cell.')
  }
  params <- lineage_benchmark_mito_params(population_bundle$params)
  population <- population_bundle$population
  cell_types <- unique(as.character(population$nodes$cell_type))
  models <- setNames(lapply(cell_types, function(cell_type){
    prepare_physicell_mito_model(
      params,
      cell_type = cell_type,
      genomes_per_cell = genomes_per_cell,
      seed = seed
    )
  }), cell_types)
  references <- vapply(models, function(model){
    paste0(model$mitochondrial_reference, collapse = '')
  }, character(1))
  if(length(unique(references)) != 1L){
    stop('Prepared mitochondrial cell-type models have different references.')
  }
  output_model <- models[[1L]]
  lineage_benchmark_log(
    sprintf(
      paste(
        'Simulating mitochondrial recording with %d biological genomes per',
        'cell across %s nodes.'
      ),
      genomes_per_cell,
      format(nrow(population$nodes), big.mark = ',')
    ),
    enabled = show_progress
  )
  simulation <- simulate_mito_on_physicell_lineage(
    population$nodes,
    models,
    params,
    terminal_physicell_ids = population$terminal_nodes$physicell_id,
    seed = seed,
    show_progress = show_progress
  )
  output <- write_physicell_mito_outputs(
    simulation,
    output_model,
    output_dir,
    show_progress = show_progress,
    compress_csv = TRUE
  )
  sampling_seed <- as.integer(
    (as.double(seed) * 104729 + 7919) %% .Machine$integer.max
  )
  sampling_orders <- lineage_benchmark_mito_sampling_orders(
    output$terminal_profiles,
    max(observation_depths),
    sampling_seed
  )
  names(sampling_orders) <- names(output$terminal_profiles)
  variant_manifest <- lineage_benchmark_mito_variant_manifest(
    output$terminal_profiles,
    sampling_orders,
    max(observation_depths),
    output_model$mitochondrial_reference
  )
  observation_matrices <- lineage_benchmark_mito_observation_matrices(
    output$terminal_profiles,
    sampling_orders,
    observation_depths,
    variant_manifest
  )
  saveRDS(output_model, file.path(output_dir, 'recording_model.rds'))
  saveRDS(sampling_orders, file.path(output_dir, 'mitochondrial_sampling_orders.rds'))
  saveRDS(
    observation_matrices,
    file.path(output_dir, 'mitochondrial_observation_matrices.rds')
  )
  write_physicell_csv(
    variant_manifest,
    file.path(output_dir, 'mitochondrial_variant_manifest.csv'),
    row.names = FALSE,
    compress = TRUE
  )
  lineage_benchmark_write_json(
    params,
    file.path(output_dir, 'recording_params.json')
  )
  writeLines('mitochondrial', file.path(output_dir, 'recording_type.txt'))
  writeLines('complete', file.path(output_dir, 'recording_complete.txt'))
  list(
    type = 'mitochondrial',
    model = output_model,
    observation_matrices = observation_matrices,
    variant_manifest = variant_manifest,
    observation_depths = observation_depths,
    genomes_per_cell = genomes_per_cell,
    params = params
  )
}

# ---- Recorder simulation and resumption ----

#' Simulate the maximum recording overlay for one recorder system
#'
#' Recording is simulated once at the largest level the sweep will ask for, so
#' every smaller condition can be produced by subsetting instead of by
#' resimulating mutation. `'mitochondrial'` is delegated to
#' `simulate_lineage_benchmark_mitochondrial()`; the integrated systems prepare
#' barcode models and record along the existing lineage.
#'
#' @param population_bundle Bundle with `population` and `params`.
#' @param system `'baseline'`, `'prime'`, `'palincode'`, or `'mitochondrial'`.
#' @param maximum_integrations Highest integration count to simulate; for
#'   `'mitochondrial'` it is used only as the fallback observation depth when
#'   `observation_depths` is `NULL`.
#' @param seed Integer seed for model preparation and recording.
#' @param output_dir Destination directory, created if it does not exist.
#' @param show_progress Whether to emit progress lines; defaults to `TRUE`.
#' @param observation_depths Mitochondrial observation depths; ignored by the
#'   integrated systems.
#' @param mitochondrial_genomes_per_cell Biological copy number used by the
#'   mitochondrial recorder; defaults to `32L`.
#' @return For the integrated systems, a named list with `type = 'integrated'`,
#'   `model`, `raw_alleles` (native allele matrix), `binary_scores`,
#'   `target_layout`, and `params`; for `'mitochondrial'`, the bundle returned
#'   by `simulate_lineage_benchmark_mitochondrial()`.
#' @section Side effects: Creates `output_dir` and writes the recording output
#'   of `write_physicell_recording_outputs()` plus `recording_model.rds`,
#'   `recording_params.json`, `recording_type.txt`, and the
#'   `recording_complete.txt` resume marker.
simulate_lineage_benchmark_recorder <- function(population_bundle,
                                                system,
                                                maximum_integrations,
                                                seed,
                                                output_dir,
                                                show_progress = TRUE,
                                                observation_depths = NULL,
                                                mitochondrial_genomes_per_cell = 32L){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  system <- tolower(as.character(system))
  if(system == 'mitochondrial'){
    if(is.null(observation_depths)){
      observation_depths <- maximum_integrations
    }
    return(simulate_lineage_benchmark_mitochondrial(
      population_bundle = population_bundle,
      observation_depths = observation_depths,
      genomes_per_cell = mitochondrial_genomes_per_cell,
      seed = seed,
      output_dir = output_dir,
      show_progress = show_progress
    ))
  }
  params <- lineage_benchmark_recording_params(
    population_bundle$params,
    system,
    maximum_integrations
  )
  models <- prepare_gillespie_barcode_models(
    params,
    num_integrations = maximum_integrations,
    founder_label_sites = 0L,
    seed = seed
  )
  output_model <- models[[1L]]
  population <- population_bundle$population
  lineage_benchmark_log(
    sprintf(
      'Simulating %s recording with %d integrations across %s nodes.',
      system,
      maximum_integrations,
      format(nrow(population$nodes), big.mark = ',')
    ),
    enabled = show_progress
  )
  simulation <- simulate_recording_on_physicell_lineage(
    population$nodes,
    models,
    terminal_physicell_ids = population$terminal_nodes$physicell_id,
    seed = seed,
    show_progress = show_progress
  )
  output <- write_physicell_recording_outputs(
    simulation,
    output_model,
    output_dir,
    show_progress = show_progress,
    write_lineage = FALSE,
    compress_csv = TRUE,
    # This benchmark reads the allele matrix back and records the mutation
    # event file in its condition manifest, so it asks for both rather than
    # taking the writer's defaults.
    write_allele_matrix = TRUE,
    write_mutation_events = TRUE
  )
  saveRDS(output_model, file.path(output_dir, 'recording_model.rds'))
  lineage_benchmark_write_json(
    params,
    file.path(output_dir, 'recording_params.json')
  )
  writeLines('complete', file.path(output_dir, 'recording_complete.txt'))
  writeLines('integrated', file.path(output_dir, 'recording_type.txt'))
  list(
    type = 'integrated',
    model = output_model,
    raw_alleles = output$raw_alleles,
    binary_scores = output$binary_scores,
    target_layout = physicell_baseline_target_layout(output_model),
    params = params
  )
}

#' Reload a completed maximum-recording directory
#'
#' Reads `recording_type.txt` to decide which artifacts are required, treating
#' a missing marker as an integrated recording, and errors when any required
#' file is absent so a partially written directory is never resumed.
#'
#' @param output_dir Directory previously written by
#'   `simulate_lineage_benchmark_recorder()`.
#' @return The same list shape the simulating function returns, except that
#'   `params` is `NULL` because the parameter block is not reloaded.
load_lineage_benchmark_recorder <- function(output_dir){
  recording_type_path <- file.path(output_dir, 'recording_type.txt')
  recording_type <- if(file.exists(recording_type_path)){
    trimws(readLines(recording_type_path, warn = FALSE)[1L])
  } else{
    'integrated'
  }
  required <- if(recording_type == 'mitochondrial'){
    c(
      'recording_complete.txt', 'recording_model.rds',
      'mitochondrial_observation_matrices.rds',
      'mitochondrial_variant_manifest.csv.gz'
    )
  } else{
    c(
      'recording_complete.txt', 'recording_model.rds',
      'barcode_alleles_sparse.rds', 'barcode_binary_score_matrix_sparse.rds'
    )
  }
  missing <- required[!file.exists(file.path(output_dir, required))]
  if(length(missing) > 0L){
    stop(sprintf(
      'Incomplete resumed recording directory %s; missing: %s.',
      output_dir,
      paste(missing, collapse = ', ')
    ))
  }
  model <- readRDS(file.path(output_dir, 'recording_model.rds'))
  if(recording_type == 'mitochondrial'){
    observation_matrices <- readRDS(file.path(
      output_dir,
      'mitochondrial_observation_matrices.rds'
    ))
    return(list(
      type = 'mitochondrial',
      model = model,
      observation_matrices = observation_matrices,
      variant_manifest = read_physicell_csv(file.path(
        output_dir,
        'mitochondrial_variant_manifest.csv'
      )),
      observation_depths = as.integer(names(observation_matrices)),
      genomes_per_cell = model$genomes_per_cell,
      params = NULL
    ))
  }
  list(
    type = 'integrated',
    model = model,
    raw_alleles = readRDS(file.path(output_dir, 'barcode_alleles_sparse.rds')),
    binary_scores = readRDS(file.path(
      output_dir,
      'barcode_binary_score_matrix_sparse.rds'
    )),
    target_layout = physicell_baseline_target_layout(model),
    params = NULL
  )
}

# ---- Condition output ----

#' Express a path relative to the benchmark root when it lies inside it
#'
#' Keeps manifest paths portable when the benchmark directory is moved.
#'
#' @param path Candidate path; it need not exist.
#' @param root Benchmark root directory, which must exist.
#' @return The root-relative path when `path` is inside `root`, otherwise the
#'   normalized `path` unchanged.
lineage_benchmark_relative_path <- function(path, root){
  path <- normalizePath(path, mustWork = FALSE)
  root <- normalizePath(root, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if(startsWith(path, prefix)){
    substring(path, nchar(prefix) + 1L)
  } else{
    path
  }
}

#' Write one reconstruction-input condition
#'
#' Subsets the single maximum recording down to this condition's cells and
#' recording level and writes everything a downstream reconstruction job needs.
#' For the integrated systems `integrations` is an integration count and
#' columns are filtered by integration number, giving native alleles in the
#' state matrix, binary scores in the character matrix, and the collapsed
#' per-target encoding in the logical target matrix. For the mitochondrial
#' recorder `integrations` is an observation depth: the matching depth panel is
#' used, state values are heteroplasmy fractions in (0, 1], the character and
#' logical target matrices are the same binary presence matrix, and features
#' are restricted to those observed in these cells at the deepest simulated
#' depth so the whole depth series shares one feature universe.
#'
#' @details The one-row manifest separates the two observation units:
#'   `integrations` is `NA` for the mitochondrial recorder and
#'   `observation_depth` is `NA` for the integrated ones, while
#'   `observations_per_cell` always carries the level and
#'   `observation_dimension` names the unit as `'integrations'` or
#'   `'sampled_mitochondrial_genomes'`. Matrix paths are stored relative to
#'   `benchmark_root`.
#'
#' @param recording Recording bundle from
#'   `simulate_lineage_benchmark_recorder()` or
#'   `load_lineage_benchmark_recorder()`.
#' @param system Recorder system name, also recorded in the manifest.
#' @param shape Population shape name recorded in the manifest.
#' @param population_seed Seed used for the population.
#' @param recorder_seed Seed used for the recording.
#' @param physicell_ids Sampled terminal PhysiCell IDs; matrix rows are taken
#'   as `cell_<id>` in this order.
#' @param integrations Integration count, or observation depth when the
#'   recording is mitochondrial.
#' @param truth_tree_path Ground-truth Newick file copied into the condition.
#' @param sample_cells_path Sampled-cell table copied into the condition.
#' @param full_recording_dir Directory holding the maximum recording, recorded
#'   in the manifest as the source of the exact mutation events.
#' @param output_dir Destination directory, created if it does not exist.
#' @param benchmark_root Root used to make the manifest paths portable.
#' @param write_dense_csv Whether to also write dense CSV copies of the three
#'   matrices; defaults to `FALSE`.
#' @return A one-row condition manifest data frame.
#' @section Side effects: Creates `output_dir` and writes the state, character,
#'   and logical-target sparse matrices (plus `mitochondrial_`-named copies for
#'   that recorder), the optional dense CSVs, `target_manifest.csv.gz`,
#'   `ground_truth_tree.nwk`, a copy of the sampled-cell table,
#'   `condition_manifest.csv.gz`, and the `condition_complete.txt` resume
#'   marker.
write_lineage_benchmark_condition <- function(recording,
                                              system,
                                              shape,
                                              population_seed,
                                              recorder_seed,
                                              physicell_ids,
                                              integrations,
                                              truth_tree_path,
                                              sample_cells_path,
                                              full_recording_dir,
                                              output_dir,
                                              benchmark_root,
                                              write_dense_csv = FALSE){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  sample_ids <- paste0('cell_', physicell_ids)
  is_mitochondrial <- identical(recording$type, 'mitochondrial') ||
    identical(system, 'mitochondrial')
  if(is_mitochondrial){
    depth_key <- as.character(integrations)
    depth_matrices <- recording$observation_matrices[[depth_key]]
    if(is.null(depth_matrices)){
      stop(sprintf('Mitochondrial observation depth %s was not simulated.', depth_key))
    }
    maximum_depth_key <- as.character(max(recording$observation_depths))
    maximum_depth_character <- lineage_benchmark_subset_rows(
      recording$observation_matrices[[maximum_depth_key]]$character,
      sample_ids
    )
    retained_features <- which(Matrix::colSums(maximum_depth_character) > 0)
    state_matrix <- lineage_benchmark_subset_rows(
      depth_matrices$state,
      sample_ids
    )[, retained_features, drop = FALSE]
    character_matrix <- lineage_benchmark_subset_rows(
      depth_matrices$character,
      sample_ids
    )[, retained_features, drop = FALSE]
    logical_target_matrix <- character_matrix
  } else{
    state_matrix <- lineage_benchmark_subset_matrix(
      recording$raw_alleles,
      sample_ids,
      integrations
    )
    # A PEtracer-style run carries which mark each site received, and that
    # identity IS the character state. Binarising it here would discard the
    # alphabet before tree building ever sees it, leaving independent edits at
    # one site indistinguishable from shared ancestry. Every other recorder,
    # and prime editing with a single mark, still gets the binary scores.
    character_matrix <- lineage_benchmark_subset_matrix(
      lineage_benchmark_character_matrix(
        recording$raw_alleles, recording$binary_scores, recording$model,
        integrations
      ),
      sample_ids,
      integrations
    )
    logical_target_matrix <- lineage_benchmark_logical_target_matrix(
      state_matrix,
      recording$model,
      integrations
    )
  }
  saveRDS(
    state_matrix,
    file.path(output_dir, 'recording_state_matrix_sparse.rds')
  )
  saveRDS(
    character_matrix,
    file.path(output_dir, 'recording_character_matrix_sparse.rds')
  )
  saveRDS(
    logical_target_matrix,
    file.path(output_dir, 'recording_logical_target_matrix_sparse.rds')
  )
  if(is_mitochondrial){
    saveRDS(
      state_matrix,
      file.path(
        output_dir,
        'mitochondrial_variant_fraction_matrix_sparse.rds'
      )
    )
    saveRDS(
      character_matrix,
      file.path(output_dir, 'mitochondrial_binary_variant_matrix_sparse.rds')
    )
  }
  if(isTRUE(write_dense_csv)){
    write_physicell_csv(
      as.matrix(state_matrix),
      file.path(output_dir, 'recording_state_matrix.csv'),
      row.names = TRUE,
      compress = TRUE
    )
    write_physicell_csv(
      as.matrix(character_matrix),
      file.path(output_dir, 'recording_character_matrix.csv'),
      row.names = TRUE,
      compress = TRUE
    )
    write_physicell_csv(
      as.matrix(logical_target_matrix),
      file.path(output_dir, 'recording_logical_target_matrix.csv'),
      row.names = TRUE,
      compress = TRUE
    )
  }
  target_layout <- if(is_mitochondrial){
    manifest_indices <- match(
      colnames(state_matrix),
      recording$variant_manifest$feature
    )
    if(anyNA(manifest_indices)){
      stop('A mitochondrial matrix feature is absent from its variant manifest.')
    }
    recording$variant_manifest[manifest_indices, , drop = FALSE]
  } else{
    retained_layout <- recording$target_layout[
      recording$target_layout$integration <= integrations,
      ,
      drop = FALSE
    ]
    lineage_benchmark_logical_target_manifest(
      retained_layout,
      recording$model
    )
  }
  write_physicell_csv(
    target_layout,
    file.path(output_dir, 'target_manifest.csv'),
    row.names = FALSE,
    compress = TRUE
  )
  if(is_mitochondrial){
    write_physicell_csv(
      target_layout,
      file.path(output_dir, 'mitochondrial_variant_manifest.csv'),
      row.names = FALSE,
      compress = TRUE
    )
  }
  file.copy(
    truth_tree_path,
    file.path(output_dir, 'ground_truth_tree.nwk'),
    overwrite = TRUE
  )
  resolved_sample_path <- resolve_physicell_csv_path(sample_cells_path)
  file.copy(
    resolved_sample_path,
    file.path(output_dir, basename(resolved_sample_path)),
    overwrite = TRUE
  )

  targets_per_integration <- lineage_benchmark_targets_per_integration(system)
  condition_id <- if(is_mitochondrial){
    sprintf(
      '%s_seed_%06d_%s_n_%05d_depth_%02d',
      shape,
      population_seed,
      system,
      length(physicell_ids),
      integrations
    )
  } else{
    sprintf(
      '%s_seed_%06d_%s_n_%05d_k_%02d',
      shape,
      population_seed,
      system,
      length(physicell_ids),
      integrations
    )
  }
  logical_targets <- ncol(logical_target_matrix)
  full_events_name <- if(is_mitochondrial){
    'mitochondrial_mutation_events.csv'
  } else{
    'mutation_events.csv'
  }
  condition <- data.frame(
    condition_id = condition_id,
    shape = shape,
    population_seed = as.integer(population_seed),
    recorder_seed = as.integer(recorder_seed),
    system = system,
    tree_size = length(physicell_ids),
    integrations = if(is_mitochondrial) NA_integer_ else integrations,
    observation_depth = if(is_mitochondrial) integrations else NA_integer_,
    observations_per_cell = integrations,
    observation_dimension = if(is_mitochondrial){
      'sampled_mitochondrial_genomes'
    } else{
      'integrations'
    },
    biological_mitochondrial_genomes = if(is_mitochondrial){
      recording$genomes_per_cell
    } else{
      NA_integer_
    },
    targets_per_integration = targets_per_integration,
    logical_targets = logical_targets,
    logical_target_matrix_columns = ncol(logical_target_matrix),
    logical_target_matrix_nonzero = length(logical_target_matrix@x),
    state_matrix_columns = ncol(state_matrix),
    character_matrix_columns = ncol(character_matrix),
    character_matrix_nonzero = length(character_matrix@x),
    ground_truth_tree = lineage_benchmark_relative_path(
      file.path(output_dir, 'ground_truth_tree.nwk'),
      benchmark_root
    ),
    state_matrix = lineage_benchmark_relative_path(
      file.path(output_dir, 'recording_state_matrix_sparse.rds'),
      benchmark_root
    ),
    character_matrix = lineage_benchmark_relative_path(
      file.path(output_dir, 'recording_character_matrix_sparse.rds'),
      benchmark_root
    ),
    logical_target_matrix = lineage_benchmark_relative_path(
      file.path(output_dir, 'recording_logical_target_matrix_sparse.rds'),
      benchmark_root
    ),
    target_manifest = lineage_benchmark_relative_path(
      physicell_csv_output_path(
        file.path(output_dir, 'target_manifest.csv'),
        TRUE
      ),
      benchmark_root
    ),
    full_recording_events = lineage_benchmark_relative_path(
      physicell_csv_output_path(
        file.path(full_recording_dir, full_events_name),
        TRUE
      ),
      benchmark_root
    ),
    stringsAsFactors = FALSE
  )
  write_physicell_csv(
    condition,
    file.path(output_dir, 'condition_manifest.csv'),
    row.names = FALSE,
    compress = TRUE
  )
  writeLines('complete', file.path(output_dir, 'condition_complete.txt'))
  condition
}

#' Reload a finished condition manifest
#'
#' @param output_dir Condition directory.
#' @return The one-row condition manifest, or `NULL` when
#'   `condition_complete.txt` is absent, which is what lets the sweep skip
#'   conditions that already finished.
load_lineage_benchmark_condition <- function(output_dir){
  if(!file.exists(file.path(output_dir, 'condition_complete.txt'))){
    return(NULL)
  }
  read_physicell_csv(file.path(output_dir, 'condition_manifest.csv'))
}

# ---- Benchmark sweep driver ----

#' Run the resumable paired population and recorder benchmark sweep
#'
#' Sweeps shapes, seeds, recorder systems, tree sizes, and recording levels
#' while simulating biology only once per cell: one population per shape and
#' seed at the largest tree size, and one recording per shape, seed, and system
#' at the largest integration count or mitochondrial observation depth. Every
#' condition is then a nested cell and recording subset of those, so conditions
#' sharing a population or a recording are paired rather than independently
#' simulated. Sample and recorder seeds are derived deterministically from the
#' population seed together with the shape and system indices.
#'
#' @details Resuming reuses any `population.rds`, `sample_sets.rds`, completed
#'   recording directory, and completed condition already on disk; a stored
#'   `benchmark_settings.rds` must match the requested grid exactly or the run
#'   stops. This function writes reconstruction inputs only. It does not
#'   reconstruct trees and computes no accuracy or tree-distance metric, so no
#'   value it returns is a reconstruction score.
#'
#' @param output_dir Benchmark root directory; if it exists and is not empty,
#'   `resume` must be `TRUE`.
#' @param shapes Population shapes to sweep, from `balanced`, `comb`,
#'   `neutral`, `turnover`, and `hierarchical`.
#' @param tree_sizes Sampled cell counts; the largest sets the simulated
#'   population size.
#' @param integration_counts Integration counts for the integrated recorders;
#'   the largest is simulated and the rest are column subsets.
#' @param systems Recorder systems to sweep, from `baseline`, `prime`,
#'   `palincode`, and `mitochondrial`.
#' @param seeds Population seeds; must be unique and non-negative.
#' @param synthetic_duration Duration passed to the synthetic `balanced` and
#'   `comb` generators; defaults to `10`.
#' @param write_dense_csv Whether each condition also writes dense CSV copies
#'   of its matrices; defaults to `FALSE`.
#' @param resume Whether an existing non-empty `output_dir` may be reused;
#'   defaults to `TRUE`.
#' @param show_progress Whether to emit progress lines; defaults to `TRUE`.
#' @param mt_observation_depths Mitochondrial observation depths; defaults to
#'   `integration_counts` and may not exceed `mt_genomes_per_cell`.
#' @param mt_genomes_per_cell Biological mitochondrial copy number per cell,
#'   one positive integer; defaults to `32L`.
#' @return Invisibly, a named list with `output_dir`, the row-bound condition
#'   `manifest`, and the `summary` data frame of simulation counts.
#' @section Side effects: Creates the whole benchmark directory tree and writes
#'   `benchmark_settings.rds`, `benchmark_settings.json`,
#'   `benchmark_manifest.csv.gz`, `benchmark_summary.csv.gz`, and
#'   `benchmark_complete.txt`, plus every population, sample, recording, and
#'   condition artifact beneath it. The nested sampling and mitochondrial
#'   sampling steps call `set.seed()`.
#' @note Requires the Matrix package.
run_lineage_benchmark <- function(output_dir,
                                  shapes = c(
                                    'balanced', 'comb', 'neutral',
                                    'hierarchical'
                                  ),
                                  tree_sizes = c(250L, 1000L, 2000L, 5000L),
                                  integration_counts = c(1L, 2L, 5L, 10L, 20L),
                                  systems = c(
                                    'baseline', 'prime', 'palincode',
                                    'mitochondrial'
                                  ),
                                  seeds = 1:10,
                                  synthetic_duration = 10,
                                  write_dense_csv = FALSE,
                                  resume = TRUE,
                                  show_progress = TRUE,
                                  mt_observation_depths = NULL,
                                  mt_genomes_per_cell = 32L){
  if(!requireNamespace('Matrix', quietly = TRUE)){
    stop('The Matrix package is required for the lineage benchmark.')
  }
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  shapes <- unique(tolower(trimws(as.character(shapes))))
  allowed_shapes <- c(
    'balanced', 'comb', 'neutral', 'turnover', 'hierarchical'
  )
  unknown_shapes <- setdiff(shapes, allowed_shapes)
  if(length(shapes) == 0L || length(unknown_shapes) > 0L){
    stop(sprintf(
      'Unknown shapes: %s.',
      paste(unknown_shapes, collapse = ', ')
    ))
  }
  systems <- unique(tolower(trimws(as.character(systems))))
  unknown_systems <- setdiff(
    systems,
    c('baseline', 'prime', 'palincode', 'mitochondrial')
  )
  if(length(systems) == 0L || length(unknown_systems) > 0L){
    stop(sprintf(
      'Unknown recorder systems: %s.',
      paste(unknown_systems, collapse = ', ')
    ))
  }
  tree_sizes <- lineage_benchmark_integer_vector(tree_sizes, 'tree_sizes')
  integration_counts <- lineage_benchmark_integer_vector(
    integration_counts,
    'integration_counts'
  )
  if(is.null(mt_observation_depths)){
    mt_observation_depths <- integration_counts
  }
  mt_observation_depths <- lineage_benchmark_integer_vector(
    mt_observation_depths,
    'mt_observation_depths'
  )
  mt_genomes_per_cell <- lineage_benchmark_integer_vector(
    mt_genomes_per_cell,
    'mt_genomes_per_cell'
  )
  if(length(mt_genomes_per_cell) != 1L){
    stop('mt_genomes_per_cell must be one positive integer.')
  }
  if('mitochondrial' %in% systems &&
     max(mt_observation_depths) > mt_genomes_per_cell){
    stop(paste(
      'The maximum mitochondrial observation depth cannot exceed',
      'mt_genomes_per_cell.'
    ))
  }
  seeds <- lineage_benchmark_integer_vector(seeds, 'seeds', minimum = 0L)
  maximum_cells <- max(tree_sizes)
  maximum_integrations <- max(integration_counts)
  settings <- list(
    benchmark_version = 2L,
    shapes = shapes,
    tree_sizes = tree_sizes,
    integration_counts = integration_counts,
    mt_observation_depths = mt_observation_depths,
    mt_genomes_per_cell = mt_genomes_per_cell,
    systems = systems,
    seeds = seeds,
    synthetic_duration = synthetic_duration,
    write_dense_csv = isTRUE(write_dense_csv)
  )

  if(dir.exists(output_dir)){
    contents <- list.files(output_dir, all.files = TRUE, no.. = TRUE)
    if(length(contents) > 0L && !isTRUE(resume)){
      stop('Benchmark output directory exists and is not empty; enable resume.')
    }
  } else{
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  settings_path <- file.path(output_dir, 'benchmark_settings.rds')
  if(file.exists(settings_path)){
    previous_settings <- readRDS(settings_path)
    if(!isTRUE(all.equal(previous_settings, settings, check.attributes = TRUE))){
      stop('Existing benchmark settings do not match this requested run.')
    }
  } else{
    saveRDS(settings, settings_path)
    lineage_benchmark_write_json(
      settings,
      file.path(output_dir, 'benchmark_settings.json')
    )
  }

  condition_rows <- list()
  condition_count <- 0L
  for(shape_index in seq_along(shapes)){
    shape <- shapes[shape_index]
    for(population_seed in seeds){
      run_dir <- file.path(
        output_dir,
        paste0('shape_', shape),
        sprintf('seed_%06d', population_seed)
      )
      population_dir <- file.path(run_dir, 'population')
      population_path <- file.path(population_dir, 'population.rds')
      if(file.exists(population_path)){
        lineage_benchmark_log(
          sprintf('Resuming %s seed %d population.', shape, population_seed),
          enabled = show_progress
        )
        population_bundle <- readRDS(population_path)
      } else{
        population_bundle <- simulate_lineage_benchmark_population(
          shape,
          maximum_cells,
          population_seed,
          synthetic_duration = synthetic_duration,
          show_progress = show_progress
        )
        write_lineage_benchmark_population(
          population_bundle,
          population_dir,
          show_progress = show_progress
        )
      }

      samples_dir <- file.path(run_dir, 'samples')
      sample_sets_path <- file.path(samples_dir, 'sample_sets.rds')
      if(file.exists(sample_sets_path)){
        sample_sets <- readRDS(sample_sets_path)
        expected_names <- as.character(tree_sizes)
        if(!identical(names(sample_sets), expected_names)){
          stop('Existing nested sample sizes do not match benchmark settings.')
        }
      } else{
        sample_seed <- as.integer(
          (as.double(population_seed) * 1009 + shape_index * 9173) %%
            .Machine$integer.max
        )
        sample_sets <- lineage_benchmark_sample_sets(
          population_bundle$population,
          tree_sizes,
          sample_seed
        )
      }
      truth_paths <- write_lineage_benchmark_samples(
        population_bundle$population,
        sample_sets,
        samples_dir,
        show_progress = show_progress
      )

      for(system_index in seq_along(systems)){
        system <- systems[system_index]
        system_dir <- file.path(run_dir, paste0('recorder_', system))
        system_levels <- if(system == 'mitochondrial'){
          mt_observation_depths
        } else{
          integration_counts
        }
        maximum_system_level <- max(system_levels)
        full_recording_dir <- if(system == 'mitochondrial'){
          file.path(
            system_dir,
            sprintf('full_%02d_genomes_per_cell', mt_genomes_per_cell)
          )
        } else{
          file.path(
            system_dir,
            sprintf('full_%02d_integrations', maximum_integrations)
          )
        }
        recorder_seed <- as.integer(
          (as.double(population_seed) * 10007 +
             shape_index * 1009 + system_index * 97) %%
            .Machine$integer.max
        )
        if(file.exists(file.path(
          full_recording_dir,
          'recording_complete.txt'
        ))){
          lineage_benchmark_log(
            sprintf(
              'Resuming %s recorder for %s seed %d.',
              system, shape, population_seed
            ),
            enabled = show_progress
          )
          recording <- load_lineage_benchmark_recorder(full_recording_dir)
        } else{
          recording <- simulate_lineage_benchmark_recorder(
            population_bundle,
            system,
            maximum_system_level,
            recorder_seed,
            full_recording_dir,
            show_progress = show_progress,
            observation_depths = mt_observation_depths,
            mitochondrial_genomes_per_cell = mt_genomes_per_cell
          )
        }

        for(size_name in names(sample_sets)){
          physicell_ids <- sample_sets[[size_name]]
          truth_path <- truth_paths[[size_name]]
          sample_cells_path <- file.path(
            samples_dir,
            paste0('n_', size_name),
            'sample_cells.csv'
          )
          for(integrations in system_levels){
            condition_dir <- file.path(
              system_dir,
              'conditions',
              paste0('n_', size_name),
              if(system == 'mitochondrial'){
                sprintf('depth_%02d', integrations)
              } else{
                sprintf('k_%02d', integrations)
              }
            )
            condition <- load_lineage_benchmark_condition(condition_dir)
            if(is.null(condition)){
              condition <- write_lineage_benchmark_condition(
                recording = recording,
                system = system,
                shape = shape,
                population_seed = population_seed,
                recorder_seed = recorder_seed,
                physicell_ids = physicell_ids,
                integrations = integrations,
                truth_tree_path = truth_path,
                sample_cells_path = sample_cells_path,
                full_recording_dir = full_recording_dir,
                output_dir = condition_dir,
                benchmark_root = output_dir,
                write_dense_csv = write_dense_csv
              )
            }
            condition_count <- condition_count + 1L
            condition_rows[[condition_count]] <- condition
          }
        }
      }
    }
  }

  manifest <- do.call(rbind, condition_rows)
  rownames(manifest) <- NULL
  write_physicell_csv(
    manifest,
    file.path(output_dir, 'benchmark_manifest.csv'),
    row.names = FALSE,
    compress = TRUE
  )
  summary <- data.frame(
    property = c(
      'population_simulations', 'full_recording_simulations',
      'derived_reconstruction_conditions', 'maximum_tree_size',
      'maximum_integrations', 'maximum_mitochondrial_observation_depth',
      'mitochondrial_genomes_per_cell'
    ),
    value = c(
      length(shapes) * length(seeds),
      length(shapes) * length(seeds) * length(systems),
      nrow(manifest),
      maximum_cells,
      maximum_integrations,
      max(mt_observation_depths),
      mt_genomes_per_cell
    ),
    stringsAsFactors = FALSE
  )
  write_physicell_csv(
    summary,
    file.path(output_dir, 'benchmark_summary.csv'),
    row.names = FALSE,
    compress = TRUE
  )
  writeLines('complete', file.path(output_dir, 'benchmark_complete.txt'))
  lineage_benchmark_log(
    sprintf(
      'Benchmark complete with %s reconstruction-input conditions.',
      format(nrow(manifest), big.mark = ',')
    ),
    enabled = show_progress
  )
  invisible(list(
    output_dir = normalizePath(output_dir, mustWork = TRUE),
    manifest = manifest,
    summary = summary
  ))
}
