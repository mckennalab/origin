#!/usr/bin/env Rscript

# Lightweight regression tests for pure helpers and static call wiring. These
# tests intentionally avoid sourcing the CLI entry points, which execute a full
# simulation and require the complete conda environment.

# ---- Repository root resolution ----

#' Resolve the repository root from this script's own path
#'
#' The suite is launched as `Rscript tests/regression_tests.R` from an arbitrary
#' working directory, so every source file and fixture below is addressed
#' relative to `repo_root` rather than to `getwd()`.
script_argument <- commandArgs(trailingOnly = FALSE)[
  grepl('^--file=', commandArgs(trailingOnly = FALSE))
][1]
script_path <- sub('^--file=', '', script_argument)
repo_root <- normalizePath(file.path(dirname(script_path), '..'), mustWork = TRUE)

# ---- Test harness helpers ----

#' Load one named top-level function out of a repository source file
#'
#' Parses `relative_path` and evaluates only the top-level assignment
#' `function_name <- function(...)`, skipping every other statement in the file.
#' This is what keeps the suite dependency-light: the CLI entry points are never
#' executed, so no simulation runs and no package outside the loaded body is
#' needed. Any free variable the loaded function refers to must be pre-seeded in
#' `envir` by the caller; several blocks below stub `sparseMatrix`,
#' `str_extract`, `parLapply` and the `parallel` cluster verbs this way.
#'
#' @param relative_path Source file path relative to `repo_root`.
#' @param function_name Name of the top-level function assignment to evaluate.
#' @param envir Environment the assignment is evaluated into; defaults to a
#'   fresh environment whose parent is the global environment.
#' @return Invisibly, the loaded function. Stops when the file contains no
#'   matching top-level function assignment.
#' @section Side effects: Binds `function_name` inside `envir`.
load_named_function <- function(relative_path, function_name, envir = new.env(parent = globalenv())){
  expressions <- parse(file.path(repo_root, relative_path))
  for(expression in expressions){
    if(is.call(expression) &&
       identical(expression[[1]], as.name('<-')) &&
       identical(expression[[2]], as.name(function_name)) &&
      is.call(expression[[3]]) &&
       identical(expression[[3]][[1]], as.name('function'))){
      eval(expression, envir = envir)
      return(invisible(envir[[function_name]]))
    }
  }
  stop(sprintf('Function %s not found in %s.', function_name, relative_path))
}

#' Assert that a value is exactly TRUE
#'
#' @param value Value under test; only a length-one, non-NA `TRUE` passes.
#' @param message Text reported when the assertion fails.
#' @return Invisible `NULL` on success; stops with `message` otherwise.
expect_true <- function(value, message){
  if(!isTRUE(value)){
    stop(message, call. = FALSE)
  }
}

#' Assert two values match, reporting both structures on failure
#'
#' @param actual Value produced by the code under test.
#' @param expected Reference value.
#' @param message Text prefixed to the failure report.
#' @param tolerance Numeric tolerance handed to `all.equal()`; the 1e-12 default
#'   is effectively exact for the integer and small-double fixtures used here.
#' @return Invisible `NULL` on success; otherwise stops with `message` followed
#'   by `str()` renderings of `expected` and `actual`. Attributes such as names
#'   and dimnames are compared, not just values.
expect_equal <- function(actual, expected, message, tolerance = 1e-12){
  if(!isTRUE(all.equal(actual, expected, tolerance = tolerance, check.attributes = TRUE))){
    stop(
      sprintf(
        '%s\nExpected: %s\nActual: %s',
        message,
        paste(capture.output(str(expected)), collapse = ' '),
        paste(capture.output(str(actual)), collapse = ' ')
      ),
      call. = FALSE
    )
  }
}

#' Assert that an expression signals an error
#'
#' @param expression Expression to evaluate; it is forced inside `tryCatch()`,
#'   so lazy evaluation defers it until the handler is installed.
#' @param message Text reported when the expression completes without error.
#' @return Invisible `NULL` when an error was raised; stops otherwise.
expect_error <- function(expression, message){
  failed <- FALSE
  tryCatch(
    force(expression),
    error = function(error) failed <<- TRUE
  )
  expect_true(failed, message)
}

#' Collect every call to a named function inside a parsed expression
#'
#' Backs the static wiring checks further down: instead of running the
#' simulator, the tests parse `sim5_code.R` and inspect the syntax tree for the
#' argument names actually supplied at each call site.
#'
#' @param expression A parsed R expression, e.g. one element of a `parse()`
#'   result.
#' @param function_name Name of the called function to look for; only calls
#'   whose head is that bare symbol match.
#' @return A list of matching call objects in traversal order.
collect_calls <- function(expression, function_name){
  matches <- list()

  #' Internal: Recursively visit one syntax node, recording matching calls
  #'
  #' Descends into a `function` node's body only, so default expressions
  #' attached to formal arguments are not searched.
  #'
  #' @param node A language object, symbol, or constant.
  #' @return Invisibly `NULL`; matches are appended to `matches` in the
  #'   enclosing frame with `<<-`.
  walk <- function(node){
    if(!is.call(node)) return()
    if(is.symbol(node[[1]]) && identical(as.character(node[[1]]), function_name)){
      matches[[length(matches) + 1L]] <<- node
    }
    if(identical(node[[1]], as.name('function'))){
      walk(node[[3]])
      return()
    }
    if(length(node) > 1){
      for(index in 2:length(node)){
        walk(node[[index]])
      }
    }
  }
  walk(expression)
  matches
}

# ---- Barcode composition and target layout ----

# Barcode composition and target placement.
barcode_env <- new.env(parent = globalenv())
generate_non_be_target_sequence <- load_named_function(
  file.path('legacy', 'add_intervening_be_targets_to_seq.r'),
  'generate_non_be_target_sequence',
  barcode_env
)
generate_target_indices <- load_named_function(
  file.path('legacy', 'add_intervening_be_targets_to_seq.r'),
  'generate_target_indices',
  barcode_env
)

#' Guard nucleotide composition for a barcode with no base-editing targets
#'
#' Requests a ten-base backbone with A/G/C/T weights 0.1/0.2/0.3/0.4 and zero BE
#' targets, then tabulates the realized bases. Asserts the exact allocation
#' 1/2/3/4.
#'
#' Regression: a zero-BE-target barcode ignored the requested nucleotide
#' fractions, and the rounding adjustment moved counts in the wrong direction
#' (BUG_AUDIT.md, "Barcode composition"). The fix allocates by largest
#' remainder, which this exact-count assertion pins down.
set.seed(1)
sequence_without_targets <- generate_non_be_target_sequence(
  barcode_length = 10,
  nuc_fracs = c(0.1, 0.2, 0.3, 0.4),
  target_from = 'A',
  be_target_count = 0
)
expect_equal(
  as.integer(table(factor(sequence_without_targets, levels = c('A', 'G', 'C', 'T')))),
  c(1L, 2L, 3L, 4L),
  'Barcode generation must honor nucleotide composition when there are no BE targets.'
)

#' Guard that BE target bases are reserved out of the returned backbone
#'
#' Asks for a ten-base barcode carrying four BE targets edited from `A`, with a
#' backbone composition that excludes `A` entirely. Asserts the returned
#' backbone holds exactly `barcode_length - be_target_count` bases and contains
#' no target base, i.e. the four target slots are reserved rather than emitted.
#'
#' Regression: same "Barcode composition" entry in BUG_AUDIT.md, whose fix
#' reserves target bases before returning the shuffled backbone.
set.seed(2)
target_reserved_sequence <- generate_non_be_target_sequence(
  barcode_length = 10,
  nuc_fracs = c(0, 0.5, 0.5, 0),
  target_from = 'A',
  be_target_count = 4
)
expect_true(
  length(target_reserved_sequence) == 6 && !('A' %in% target_reserved_sequence),
  'Reserved targets must leave exactly barcode_length - be_target_count backbone bases.'
)

#' Reject a spaced target layout that runs past the end of the barcode
#'
#' Three spaced targets starting at position 8 with one intervening base would
#' need positions 8, 10 and 12 inside a ten-base barcode, which is impossible.
#'
#' Regression: invalid configurations could leave `all_inds` undefined or place
#' spaced targets beyond the barcode (BUG_AUDIT.md, "Target layout"); the layout
#' is now validated and fails explicitly.
expect_error(
  generate_target_indices('S', 3, 8, 10, 1),
  'An out-of-bounds spaced target configuration must fail explicitly.'
)

# ---- Deletion run grouping in score matrices ----

# Deletion runs.
score_env <- new.env(parent = globalenv())
group_deletions <- load_named_function(file.path('legacy', 'mut_to_scoremat.r'), 'group_deletions', score_env)

#' Condense only genuinely adjacent deletions into one run
#'
#' Feeds `group_deletions()` four deleted positions: 1 and 2 are adjacent within
#' cell1/integration1, 5 is isolated in the same cell/integration, and 9 belongs
#' to cell2. Asserts exactly three runs come back - one two-base deletion in
#' cell1 plus single-base deletions in cell1 and cell2 - encoded as `d2`/`d1`
#' in the mutation-value column.
#'
#' Regression: a misplaced parenthesis in the adjacency test collapsed every
#' deletion in one cell/integration into a single run, and boundary rows were
#' attributed to the wrong run (BUG_AUDIT.md, "Score matrices").
deletions <- data.frame(
  linstring = c('cell1', 'cell1', 'cell1', 'cell2'),
  ints_mutated = c(1, 1, 1, 1),
  positions_mutated = c(1, 2, 5, 9),
  mut_vals = -1
)
grouped <- group_deletions(deletions)
expected_grouped <- rbind(
  c('cell1', '1', '2', 'd2'),
  c('cell1', '1', '5', 'd1'),
  c('cell2', '1', '9', 'd1')
)
colnames(expected_grouped) <- colnames(deletions)
expect_equal(grouped, expected_grouped, 'Only adjacent deletions may be condensed into one run.')

# ---- Prime-editing guide translation ----

# Prime-editing guide translation.
prime_env <- new.env(parent = globalenv())
prime_env$int_to_nuc_list <- list('1' = 'A', '2' = 'G', '3' = 'C', '4' = 'T')
create_prime_editing_basepos_seqs <- load_named_function(
  file.path('legacy', 'sim5_code.R'),
  'create_prime_editing_basepos_seqs',
  prime_env
)

#' Translate each integer prime-editing guide through its own lookup entry
#'
#' Generates three five-base guides for three target indices, then
#' independently re-translates the returned integer guides through
#' `int_to_nuc_list` and compares against the nucleotide guides the function
#' produced. Nucleotides are encoded `A, G, C, T` as `1, 2, 3, 4`.
#'
#' Regression: integer guide translation returned the entire lookup table for
#' every base, so every guide collapsed to the same malformed string
#' (BUG_AUDIT.md, "Prime editing").
set.seed(3)
prime_guides <- create_prime_editing_basepos_seqs(1:3, guide_length = 5, num_unique_guides = 3)
translated <- vapply(prime_guides$ind_to_prime_seq_int_map, function(ints){
  paste(vapply(ints, function(int) prime_env$int_to_nuc_list[[as.character(int)]], character(1)),
        collapse = '')
}, character(1))
expect_equal(
  unname(unlist(prime_guides$ind_to_prime_seq_nuc_map)),
  unname(translated),
  'Prime-editing nucleotide guides must translate their corresponding integer guides.'
)

# ---- Known-pegRNA prime-editing backend ----

# The shared backend assigns known pegRNAs to recorder targets and treats each
# pegRNA efficiency as a multiplier of the base editing hazard.
prime_backend_env <- new.env(parent = globalenv())
sys.source(file.path(repo_root, 'origin', 'R', 'prime_editing.R'), envir = prime_backend_env)

#' Explicit pegRNA assignments keep target order and permit pool reuse
#'
#' Wires three recorder targets to the pegRNA IDs `peg_high`, `peg_zero`,
#' `peg_high` drawn from a two-entry pool. Asserts the prepared backend reports
#' the pegRNA IDs in target order, so a pegRNA used at more than one target is
#' neither deduplicated nor reordered.
known_prime_params <- list(
  nuclease_targets = list(prime_editing_system = FALSE),
  prime_editing_backend = list(
    enabled = TRUE,
    target_pegRNA_ids = list('peg_high', 'peg_zero', 'peg_high'),
    pegRNAs = list(
      list(
        pegRNA_id = 'peg_high',
        edit_sequence = 'AGGCT',
        editing_efficiency = 1
      ),
      list(
        pegRNA_id = 'peg_zero',
        edit_sequence = 'CTTGA',
        editing_efficiency = 0
      )
    )
  )
)
known_prime_backend <- prime_backend_env$prepare_prime_editing_backend(
  known_prime_params,
  target_positions = c(2L, 4L, 6L),
  seed = 5L
)
expect_equal(
  known_prime_backend$targets$pegRNA_id,
  c('peg_high', 'peg_zero', 'peg_high'),
  'Explicit pegRNA assignments must retain target order and permit pool reuse.'
)

#' A pegRNA efficiency scales editing on the complementary-survival hazard scale
#'
#' Asserts `prime_editing_scale_probability(0.36, 0.5)` is `0.2`: the efficiency
#' acts as an exponent on the survival probability, 1 - (1 - 0.36)^0.5, rather
#' than as a linear multiplier of the per-cell-cycle edit probability.
expect_equal(
  prime_backend_env$prime_editing_scale_probability(0.36, 0.5),
  0.2,
  'A pegRNA efficiency must scale editing on the complementary-survival hazard scale.'
)

#' Known pegRNA pools resolve relative to the parameter JSON directory
#'
#' Loads `../data/example_pegRNA_pool.csv` with `params_dir` pointed at
#' `example_json_params`, mirroring how a parameter file names a pool by a path
#' relative to itself. Asserts four pegRNAs load with efficiencies inside
#' [0, 1].
example_prime_pool <- prime_backend_env$read_prime_editing_pool(
  list(prime_editing_backend = list(
    pegRNA_pool_path = '../data/example_pegRNA_pool.csv'
  )),
  params_dir = file.path(repo_root, 'example_json_params')
)
expect_true(
  nrow(example_prime_pool) == 4L &&
    all(example_prime_pool$editing_efficiency >= 0) &&
    all(example_prime_pool$editing_efficiency <= 1),
  'Known pegRNA pools must load relative to the parameter JSON directory.'
)

#' Reject a pegRNA whose editing efficiency lies outside [0, 1]
#'
#' An inline pegRNA definition with `editing_efficiency = 1.1` must be refused
#' when the pool is read, rather than yielding an out-of-range hazard later.
expect_error(
  prime_backend_env$read_prime_editing_pool(list(
    prime_editing_backend = list(pegRNAs = list(list(
      pegRNA_id = 'invalid',
      edit_sequence = 'ACGT',
      editing_efficiency = 1.1
    )))
  )),
  'PegRNA editing efficiencies outside [0, 1] must be rejected.'
)

# ---- Editing windows, target classes, induction, and founder setup ----

window_env <- new.env(parent = globalenv())
parse_target_count_arguments <- load_named_function(
  file.path('legacy', 'sim5_code.R'),
  'parse_target_count_arguments',
  window_env
)
gcd <- load_named_function(file.path('legacy', 'sim5_code.R'), 'gcd', window_env)
get_new_be_targets <- load_named_function(file.path('legacy', 'sim5_code.R'), 'get_new_be_targets', window_env)
get_new_nuc_targets <- load_named_function(file.path('legacy', 'sim5_code.R'), 'get_new_nuc_targets', window_env)
sample_induced_cells <- load_named_function(
  file.path('legacy', 'sim5_code.R'),
  'sample_induced_cells',
  window_env
)
initialize_founder_population <- load_named_function(
  file.path('legacy', 'sim5_code.R'),
  'initialize_founder_population',
  window_env
)
create_ground_truth_tree <- load_named_function(
  file.path('legacy', 'sim5_code.R'),
  'create_ground_truth_tree',
  window_env
)
setup_sim <- load_named_function(
  file.path('legacy', 'sim5_code.R'),
  'setup_sim',
  window_env
)

#' GCD of rates written in scientific notation
#'
#' Asserts `gcd(1e-5, 3e-5)` is `1e-5`. The simulator derives its time increment
#' from the GCD of the configured rates, which the parameter JSON commonly
#' writes in scientific notation.
expect_equal(
  gcd(1e-5, 3e-5),
  1e-5,
  'GCD calculation must support values printed in scientific notation.'
)

#' Target-class rounding must never produce a negative residual class
#'
#' Splits three targets across high/medium/low fractions of 0.5/0.5/0 and
#' asserts the exact allocation 2/1/0.
#'
#' Regression: the high and medium counts were rounded independently and the low
#' count taken as whatever remained, which could go negative (BUG_AUDIT.md,
#' "Target classes"). All three classes are now allocated by largest remainder.
expect_equal(
  parse_target_count_arguments(
    3,
    list(high = 0.5, medium = 0.5, low = 0)
  ),
  list(num_h = 2, num_m = 1, num_l = 0),
  'Target-class rounding must never produce a negative residual class.'
)

#' A zero-width editing window returns initialized empty collections
#'
#' Calls `get_new_nuc_targets()` with an editing window of zero and asserts both
#' returned elements are empty lists rather than unbound names.
#'
#' Regression: a direct zero-width call returned an uninitialized
#' `growing_window_editrates` (BUG_AUDIT.md, "Editing windows").
expect_equal(
  get_new_nuc_targets(0, list('2' = 'High'), FALSE, 1:4),
  list(growing_window_editrates = list(), nuc_window_to_target_ind_list = list()),
  'A zero-width editing window must return initialized empty collections.'
)

#' One window name must hold both the target and its expanded positions
#'
#' Expands a single BE target at position 2 with an editing window of one over a
#' three-base barcode. Asserts the reverse mapping lists the target together
#' with the neighboring eligible position under the one key `be_window_1`.
#'
#' Regression: expanded positions were stored under target-position keys while
#' the target itself was stored under `be_window_*`/`nuc_window_*`, so the
#' reverse mappings described separate, incomplete windows (BUG_AUDIT.md,
#' "Editing windows").
be_window <- get_new_be_targets(
  be_editing_window = 1,
  basepos_erc_be_list = list('2' = 'High'),
  decaying_editing = FALSE,
  baseline_seq_ints_bc = c(1, 1, 2)
)
expect_equal(
  be_window$be_window_to_target_ind_list,
  list(be_window_1 = c(2L, 1L)),
  'Editing-window mappings must keep the target and expanded positions under one window name.'
)

#' Fractional induction works for a population larger than one cell
#'
#' Draws a 50% editing induction from five living cells and asserts the result
#' is a character subset of those cells.
#'
#' Regression: fractional induction drew one Bernoulli value per living cell and
#' passed the resulting vector as `sample(size = ...)`, which crashed for any
#' population above one cell (BUG_AUDIT.md, "Editing induction"). A single
#' binomial count is now drawn and that many cells sampled.
set.seed(10)
fractionally_induced <- sample_induced_cells(
  letters[1:5],
  frac_cells = 0.5
)
expect_true(
  is.character(fractionally_induced) &&
    length(fractionally_induced) <= 5 &&
    all(fractionally_induced %in% letters[1:5]),
  'Fractional induction must return a valid subset when multiple cells are alive.'
)

#' Fixed-count induction caps selection at the number of living cells
#'
#' Requests ten induced cells from a population of five and asserts all five are
#' returned, with no error and no padding.
expect_equal(
  sort(sample_induced_cells(letters[1:5], num_cells = 10)),
  letters[1:5],
  'Fixed-count induction must cap selection at the number of living cells.'
)

#' Every configured founder gets its own record, schedule, and induction state
#'
#' Initializes three founders with distinct division points, marking founders 1
#' and 3 for time-zero editing induction and founder 2 for differentiation.
#' Asserts the population is named `1`, `2`, `3`, that each founder keeps its
#' own `elig_div_points`, and that `induced_editing` is assigned per founder.
#'
#' Regression: `num_init_cells` was accepted by `setup_sim()` but population
#' initialization always created only cell `1` (BUG_AUDIT.md, "Founder
#' population").
founder_population <- initialize_founder_population(
  init_pop_size = 3,
  founder_cell_type = 'ct1',
  division_points_by_founder = list(1, 2, 3),
  init_incoming_mt_profile = matrix(0, nrow = 2, ncol = 2),
  init_incoming_bc_profile = matrix(0, nrow = 2, ncol = 2),
  mito_to_genome_map = list('1' = 1, '2' = 2),
  initial_heteroplasmy_score = 0,
  init_heteroplasmy_survive_prob = 1,
  editing_induced_founders = c('1', '3'),
  differentiation_induced_founders = '2'
)
expect_equal(
  names(founder_population),
  c('1', '2', '3'),
  'num_init_cells must create one independently named record per founder.'
)
expect_equal(
  unname(vapply(founder_population, `[[`, numeric(1), 'elig_div_points')),
  c(1, 2, 3),
  'Each founder must retain its independently drawn division schedule.'
)
expect_equal(
  unname(vapply(founder_population, `[[`, character(1), 'induced_editing')),
  c('induced_editing_params', 'uninduced_editing_params', 'induced_editing_params'),
  'Time-zero editing induction must be assigned per founder.'
)

#' Multiple founders join under one synthetic time-zero root in Newick output
#'
#' Writes ground-truth trees for the three-founder population and for a single
#' undivided founder into a temporary output root. Asserts the multi-founder
#' tree carries three tips under one internal node, that the exported Newick
#' contains every initialized cell, and that a lone undivided founder still
#' yields valid Newick (`1;`).
#'
#' Regression: same "Founder population" entry in BUG_AUDIT.md - multiple roots
#' are now joined beneath a synthetic time-zero root.
#'
#' @section Side effects: Creates and then removes a temporary directory holding
#'   `processed_newicks/test/*.newick`.
founder_tree_output <- tempfile('remote_mito_founder_tree_')
dir.create(founder_tree_output)
founder_tree <- create_ground_truth_tree(
  founder_population,
  urid = 'test',
  save_path_stem = 'three_founders',
  output_root = founder_tree_output
)
expect_equal(
  c(length(founder_tree$tip.label), founder_tree$Nnode),
  c(3L, 1L),
  'Multiple founders must be joined by one synthetic time-zero tree root.'
)
expect_equal(
  sort(ape::read.tree(file.path(
    founder_tree_output,
    'processed_newicks',
    'test',
    'three_founders.newick'
  ))$tip.label),
  c('1', '2', '3'),
  'The exported multi-founder Newick must contain every initialized cell.'
)
single_founder_tree <- create_ground_truth_tree(
  founder_population[1],
  urid = 'test',
  save_path_stem = 'one_founder',
  output_root = founder_tree_output
)
expect_equal(
  readLines(file.path(
    founder_tree_output,
    'processed_newicks',
    'test',
    'one_founder.newick'
  )),
  '1;',
  'An undivided single-founder population must still produce valid Newick.'
)
unlink(founder_tree_output, recursive = TRUE)

# ---- setup_sim founder count and worker RNG wiring ----

# Exercise setup_sim with a non-networked cluster double so the configured
# founder count and worker seed wiring are covered without opening sockets.

#' setup_sim honors the configured founder count and seeds worker RNG streams
#'
#' Stubs `makeCluster`, `clusterSetRNGStream`, `clusterEvalQ` and
#' `clusterExport` inside the loaded function's environment so no socket is
#' opened, then calls `setup_sim()` with every formal supplied by name. The
#' `clusterSetRNGStream` stub records both the seed it was handed and the
#' main-process `.Random.seed`, then deliberately calls `set.seed(999)` so any
#' failure to restore the main RNG state becomes observable.
#'
#' Asserts three founders are created, that the fixed-count time-zero editing
#' induction selects exactly two of them, that the configured `random_seed`
#' reaches `clusterSetRNGStream()` as `iseed`, and that seeding the workers
#' leaves the main simulation RNG state unchanged.
#'
#' Regression: `num_init_cells` was ignored by population initialization
#' (BUG_AUDIT.md, "Founder population"), and the configured seed initialized
#' only the main R process, leaving PSOCK worker streams uncontrolled
#' (BUG_AUDIT.md, "Parallel reproducibility").
window_env$sparseMatrix <- Matrix::sparseMatrix
window_env$makeCluster <- function(num_clusters){
  list(num_clusters = num_clusters)
}
window_env$clusterSetRNGStream <- function(cl, iseed){
  window_env$observed_worker_seed <- iseed
  window_env$master_state_before_worker_seed <-
    get('.Random.seed', envir = .GlobalEnv)
  set.seed(999)
  invisible(cl)
}
window_env$clusterEvalQ <- function(...) invisible(NULL)
window_env$clusterExport <- function(...) invisible(NULL)
window_env$input_args <- list(
  random_seed = 91,
  editing_induction = list(timepoint = 0, num_cells = 2, frac_cells = NULL),
  differentiation_induction = list(timepoint = 1, num_cells = 1, frac_cells = NULL)
)
setup_arguments <- setNames(
  vector('list', length(formals(setup_sim))),
  names(formals(setup_sim))
)
setup_arguments[c(
  'num_clusters', 'init_pop_size', 'sim_length', 'cell_type_cell_cycle_length',
  'num_rows_mt', 'num_cols_mt', 'num_rows_bc', 'num_cols_bc', 'time_inc',
  'init_incoming_mt_profile', 'founder_cell_type', 'mito_to_genome_map',
  'initial_heteroplasmy_score', 'init_heteroplasmy_survive_prob'
)] <- list(
  2, 3, 2, list(ct1 = 1),
  2, 2, 2, 2, 1,
  matrix(0, nrow = 2, ncol = 2), 'ct1', list('1' = 1, '2' = 2),
  0, 1
)
set.seed(91)
setup_population <- do.call(setup_sim, setup_arguments)
expect_equal(
  names(setup_population),
  c('1', '2', '3'),
  'setup_sim must honor num_init_cells rather than hard-coding one founder.'
)
expect_equal(
  sum(vapply(
    setup_population,
    function(cell) cell$induced_editing == 'induced_editing_params',
    logical(1)
  )),
  2L,
  'Time-zero fixed-count editing induction must select the requested number of founders.'
)
expect_equal(
  window_env$observed_worker_seed,
  91,
  'setup_sim must pass random_seed to the parallel RNG stream initializer.'
)
expect_equal(
  get('.Random.seed', envir = .GlobalEnv),
  window_env$master_state_before_worker_seed,
  'Initializing worker RNG streams must not reset the main simulation RNG state.'
)

# ---- Substitution models and rate heterogeneity ----

# Substitution and heterogeneity helpers.
sub_env <- new.env(parent = globalenv())
sub_env$harmonic.mean <- function(values) length(values) / sum(1 / values)
hky_sub_rate_mat <- load_named_function(file.path('legacy', 'substitution_models.r'), 'hky_sub_rate_mat', sub_env)
simplify_target_rates <- load_named_function(
  file.path('legacy', 'substitution_models.r'),
  'SIMPLIFY_target_site_gamma_based_sub_rates',
  sub_env
)
scale_nontarget_rates <- load_named_function(
  file.path('legacy', 'substitution_models.r'),
  'nontarget_scale_gamma_heterogeneity',
  sub_env
)

#' HKY applies the transition multiplier to A<->G and C<->T only
#'
#' Builds an HKY rate matrix with uniform base fractions and a
#' transition/transversion ratio of three, then checks two ratios in the
#' `A, G, C, T` row/column order: A->G over A->C, and C->T over C->A. Both must
#' equal the configured ratio.
#'
#' Regression: the transition multiplier was applied to the wrong cells - mostly
#' destination G - and multiplied an already transition-scaled baseline a second
#' time (BUG_AUDIT.md, "HKY").
hky <- hky_sub_rate_mat(
  frac_a = 0.25,
  frac_g = 0.25,
  frac_c = 0.25,
  frac_t = 0.25,
  transition_to_transversion_ratio = 3,
  baseline_transition_rate = 3e-6,
  baseline_transversion_rate = 1e-6
)
expect_true(
  abs(hky[1, 2] / hky[1, 3] - 3) < 1e-12 &&
    abs(hky[3, 4] / hky[3, 1] - 3) < 1e-12,
  'HKY must apply the transition multiplier to A<->G and C<->T only.'
)

#' A target-rate request with no targets returns an empty list
#'
#' Calls the target-site gamma rate generator with zero targets and empty index
#' vectors, asserting an empty list rather than an error.
#'
#' Regression: the no-target test checked `length(total_num_targets)`, which is
#' always one for a scalar, and then attempted an invalid gamma cut
#' (BUG_AUDIT.md, "Target-rate generation").
expect_equal(
  simplify_target_rates(10, integer(), integer(), integer()),
  list(),
  'A target-rate request with no targets must return an empty list.'
)

#' A gamma shape of zero means no heterogeneity, not zero mutation rates
#'
#' Passes a nested background-rate list - one scalar position and one
#' per-destination-base position - through the heterogeneity scaler with
#' `shape_param = 0` and asserts the rates return unchanged.
#'
#' Regression: `shape_param = 0` replaced every background mutation rate with
#' zero, and a degenerate gamma draw replaced rates instead of scaling them
#' (BUG_AUDIT.md, "Heterogeneity").
base_rates <- list('1' = 0.1, '2' = list(A = 0.2, G = 0.3))
expect_equal(
  scale_nontarget_rates(base_rates, shape_param = 0),
  base_rates,
  'shape_param = 0 must mean no heterogeneity, not zero mutation rates.'
)

# ---- Nonuniform mutation sampling, editing windows, and transversions ----

# Target transversions must receive a destination-base vector.
mutation_env <- new.env(parent = globalenv())
mutation_env$sparseMatrix <- Matrix::sparseMatrix
filter_elig_ints_by_edit_window <- load_named_function(
  file.path('legacy', 'nonuniform_muts_heterogeneous.R'),
  'filter_elig_ints_by_edit_window',
  mutation_env
)
load_named_function(file.path('legacy', 'nonuniform_muts_heterogeneous.R'), 'non_uniform_editing', mutation_env)
load_named_function(file.path('legacy', 'nonuniform_muts_heterogeneous.R'), 'get_background_edit_inds', mutation_env)
transversion_func <- load_named_function(
  file.path('legacy', 'nonuniform_muts_heterogeneous.R'),
  'transversion_func',
  mutation_env
)

#' A rate-one target edits every eligible integration once, keeping its name
#'
#' Runs target editing for a single position `'1'` at probability one across two
#' eligible integrations. Asserts both integrations appear exactly once in
#' `i_coords` and that the position name survives into `j_coords`.
#'
#' Regression: `sapply()` simplification discarded target-position names for
#' multi-integration results, and sampling with replacement produced fewer
#' unique edits than the binomial draw called for (BUG_AUDIT.md, "Target
#' mutation sampling").
set.seed(4)
target_edit_coordinates <- mutation_env$non_uniform_editing(
  pos_er_list = list('1' = 1),
  num_integrations = 2,
  eligible_ints = list('1' = 1:2)
)
expect_true(
  !identical(target_edit_coordinates$i_coords, FALSE) &&
    identical(sort(as.integer(target_edit_coordinates$i_coords)), 1:2) &&
    identical(as.integer(target_edit_coordinates$j_coords), c(1L, 1L)),
  paste(
    'A rate-one target must edit every eligible integration exactly once while',
    'retaining its position name.'
  )
)

#' Background edits are drawn per position and never duplicate a coordinate
#'
#' Uses a four-integration, two-position background rate list with position `1`
#' at probability one and position `2` at probability zero. Asserts four edits
#' are drawn, that every integration is edited once at position 1, and that the
#' probability-zero position is never touched.
#'
#' Regression: position-specific probabilities were averaged before drawing one
#' global edit count, then coordinates were sampled with replacement - which
#' distorted heterogeneous rates and lost edits to duplicate coordinates
#' (BUG_AUDIT.md, "Background mutation sampling").
background_edit_coordinates <- mutation_env$get_background_edit_inds(
  num_rows = 4,
  num_cols = 2,
  bg_pos_er_list = list('1' = 1, '2' = 0),
  mut_type = 'transition'
)
expect_equal(
  background_edit_coordinates$num_edits,
  4,
  'Independent background sampling must honor a per-position probability of one.'
)
expect_equal(
  sort(as.integer(background_edit_coordinates$i_coords)),
  1:4,
  'Every integration must be edited at a probability-one background position.'
)
expect_equal(
  as.integer(background_edit_coordinates$j_coords),
  rep(1L, 4),
  'A probability-zero background position must never be edited.'
)

#' A saturated position closes its entire editing window
#'
#' Two positions share `window1`, and position `2` has no unedited integrations
#' left. Asserts both positions come back with an empty eligible set.
#'
#' Regression: an empty eligible set at one position was discarded before the
#' intersection, which reopened a window that should have been closed
#' (BUG_AUDIT.md, "Editing windows").
closed_window <- filter_elig_ints_by_edit_window(
  pos_to_window_inds_list = list('1' = 'window1', '2' = 'window1'),
  window_to_pos_inds_list = list('window1' = c(1, 2)),
  pos_to_unedited_int_list = list('1' = 1:2, '2' = integer(0))
)
expect_equal(
  closed_window,
  list('1' = integer(0), '2' = integer(0)),
  'An empty eligible set at one edited position must close its entire editing window.'
)

#' A forced target transversion writes the configured BE destination base
#'
#' Runs the transversion path over a one-integration, one-position matrix with
#' the target transversion rate at one, `force_target_transversions = TRUE`, and
#' destination base `4` (T, in the `A, G, C, T` = `1, 2, 3, 4` encoding).
#' Asserts the mutation matrix holds `4`.
#'
#' Regression: the nonuniform transversion path omitted the required destination
#' bases and crashed, and forced conversions did not use `be_conversion_pattern`
#' (BUG_AUDIT.md, "Target transversions").
set.seed(4)
target_transversion <- transversion_func(
  mut_mat = Matrix::Matrix(0, nrow = 1, ncol = 1, sparse = TRUE),
  num_rows = 1,
  num_cols = 1,
  bg_transversion_pos_er_list = list('1' = list(C = 0, T = 0)),
  baseline_ints = 1L,
  bg_sub_prob_mat = matrix(c(0, 0, 0.5, 0.5,
                             0, 0, 0.5, 0.5,
                             0.5, 0.5, 0, 0,
                             0.5, 0.5, 0, 0), nrow = 4, byrow = TRUE),
  target_transversion_pos_er_list = list('1' = 1),
  force_target_transversions = TRUE,
  target_transversion_to_base = 4L
)
expect_equal(
  as.numeric(target_transversion[1, 1]),
  4,
  'A forced target transversion must produce the base from be_conversion_pattern.'
)

# ---- Robinson-Foulds result filename parsing ----

# RF filename parsing must preserve distinct per-cell-type sampling fractions.
results_env <- new.env(parent = globalenv())
results_env$str_extract <- stringr::str_extract
results_env$str_match <- stringr::str_match
results_env$str_match_all <- stringr::str_match_all
extract_subrun_details <- load_named_function(
  file.path('legacy', 'process_results_from_bash.r'),
  'extract_subrun_details',
  results_env
)

#' Each cell type keeps its own sampling fraction when RF filenames are parsed
#'
#' Parses a result filename whose sampling block encodes `ct1-0.25_ct2-0.75` and
#' asserts both fractions survive against their own cell types.
#'
#' Regression: every parsed cell type received the first sampling fraction found
#' in the filename (BUG_AUDIT.md, "Result aggregation"); the contiguous
#' `cell-type-rate` block is now parsed in one pass.
sampling_result <- extract_subrun_details(
  'fasta_proc_bc_list_3_ints_RP_0.5_samp_ct1-0.25_ct2-0.75_res_time_10_TERM_rf.txt',
  'bc',
  'bc_cell_rec_fracs',
  'fasta'
)
expect_equal(
  unname(sampling_result),
  'ct1: 0.25, ct2: 0.75',
  'Each cell type must retain its own sampling fraction when RF filenames are parsed.'
)

# ---- Single-cell snapshot covariates ----

# Repeated lineages across single-cell snapshots need unique observation IDs.
sc_env <- new.env(parent = globalenv())
build_covariate_df <- load_named_function(
  file.path('legacy', 'generate_sc_profiles_from_bash.r'),
  'build_covariate_df',
  sc_env
)

#' Snapshot observations get unique IDs and exclude nonterminal ancestors
#'
#' Writes the same two-record population - one cell that is both `alive` and
#' `terminal`, plus a still-`alive` but nonterminal parent - to two snapshot RDS
#' files taken at timepoints 1 and 2, then builds the covariate frame. Asserts
#' the leaf appears once per timepoint as `1_time_1` and `1_time_2`, and that
#' the historical parent contributes no row.
#'
#' Regression: one lineage sampled at several stopping points produced duplicate
#' SCE column/row names, and historical parents that remain `alive` but are no
#' longer terminal leaves were resampled at later snapshots (BUG_AUDIT.md,
#' "Single-cell profiles").
#'
#' @section Side effects: Writes and then removes two RDS files under
#'   `tempdir()`.
snapshot_paths <- file.path(tempdir(), c('remote_mito_snapshot_1.rds', 'remote_mito_snapshot_2.rds'))
cell_record <- list(
  alive = TRUE,
  terminal = TRUE,
  linstring = '1',
  celltype = 'ct1',
  birth_time = 0,
  induced_editing = 'uninduced_editing_params'
)
internal_record <- cell_record
internal_record$terminal <- FALSE
internal_record$linstring <- 'parent'
saveRDS(list('1' = cell_record, parent = internal_record), snapshot_paths[1])
saveRDS(list('1' = cell_record, parent = internal_record), snapshot_paths[2])
snapshot_meta <- build_covariate_df(
  data.frame(path = snapshot_paths, timepoint = c(1, 2), stringsAsFactors = FALSE)
)
unlink(snapshot_paths)
expect_equal(
  snapshot_meta$sample_id,
  c('1_time_1', '1_time_2'),
  paste(
    'Repeated lineage IDs at different stopping points must receive unique sample IDs,',
    'and historical nonterminal ancestors must not be resampled.'
  )
)

# ---- Static call-wiring checks against sim5_code.R ----

# Keep the main simulator's long argument lists wired to their definitions.
sim_expressions <- parse(file.path(repo_root, 'legacy', 'sim5_code.R'))

#' Retrieve one parsed top-level assignment from the simulator source
#'
#' Filters the already-parsed `sim_expressions` for the first top-level
#' `name <- ...` assignment, so formals and call arguments can be compared
#' without evaluating any simulator code.
#'
#' @param name Name of the assigned symbol to look up.
#' @return The matching call object; subscripting fails when the simulator has
#'   no such top-level assignment.
top_level_assignment <- function(name){
  Filter(function(expression){
    is.call(expression) &&
      identical(expression[[1]], as.name('<-')) &&
      identical(expression[[2]], as.name(name))
  }, as.list(sim_expressions))[[1]]
}

#' sim_arglist must exactly match the setup_sim signature
#'
#' Compares the names in the `sim_arglist` list literal against the formals of
#' `setup_sim`, both read from the parsed syntax tree rather than from a run.
#'
#' Regression: `setup_sim()` required `cold_startup` and `cell_population`,
#' which `sim_arglist` never supplied, so the run failed before simulation began
#' (BUG_AUDIT.md, "Main call wiring").
setup_formal_names <- names(as.list(top_level_assignment('setup_sim')[[3]][[2]]))
sim_arg_names <- names(as.list(top_level_assignment('sim_arglist')[[3]])[-1])
expect_equal(
  sort(sim_arg_names),
  sort(setup_formal_names),
  'sim_arglist names must exactly match setup_sim formals.'
)

#' The single multi_core_func call must supply every formal exactly once
#'
#' Walks the parsed simulator for `multi_core_func()` call sites, asserts there
#' is exactly one, and compares its supplied argument names against the
#' function's formals.
#'
#' Regression: the sole call omitted the score-matrix controls and the initial
#' mt-genome count (BUG_AUDIT.md, "Main call wiring").
multi_formal_names <- names(as.list(top_level_assignment('multi_core_func')[[3]][[2]]))
multi_calls <- unlist(
  lapply(sim_expressions, collect_calls, function_name = 'multi_core_func'),
  recursive = FALSE
)
expect_true(length(multi_calls) == 1, 'Expected exactly one multi_core_func invocation.')
expect_equal(
  sort(names(as.list(multi_calls[[1]])[-1])),
  sort(multi_formal_names),
  'The main multi_core_func invocation must supply every formal exactly once.'
)

#' The configured random seed must reach the parallel worker RNG streams
#'
#' Asserts the parsed simulator contains exactly one `clusterSetRNGStream()`
#' call site.
#'
#' Regression: the configured seed initialized only the main R process, leaving
#' PSOCK worker streams uncontrolled (BUG_AUDIT.md, "Parallel reproducibility").
cluster_rng_calls <- unlist(
  lapply(sim_expressions, collect_calls, function_name = 'clusterSetRNGStream'),
  recursive = FALSE
)
expect_true(
  length(cluster_rng_calls) == 1,
  'The configured random seed must initialize the parallel worker RNG streams.'
)

# ---- Molecular recovery and FASTA export ----

# Zero molecular recovery must produce a well-formed empty profile.
recovery_env <- new.env(parent = globalenv())
get_profiles_ints_and_umis <- load_named_function(
  file.path('legacy', 'mut_to_fasta_difflen_ints.r'),
  'get_profiles_ints_and_umis',
  recovery_env
)
load_named_function(
  file.path('legacy', 'mut_to_fasta_difflen_ints.r'),
  'get_one_cell_sequence',
  recovery_env
)
write_all_cell_sequences <- load_named_function(
  file.path('legacy', 'mut_to_fasta_difflen_ints.r'),
  'write_all_cell_sequences',
  recovery_env
)

#' Zero recovery probability yields a well-formed empty profile
#'
#' Recovers barcode integrations from a one-cell population whose profiles hold
#' three integrations, with `int_rec_prob = 0`. Asserts a zero-row mutation
#' matrix of the right column count, no recovered integration indices, and no
#' reported UMIs.
#'
#' Regression: recovery used `max(binomial_draw, 1)`, which made zero recovery
#' impossible and inflated low recovery rates (BUG_AUDIT.md, "Molecular
#' recovery").
recovery_population <- list(
  cell1 = list(
    incoming_bc_profiles = matrix(1:6, nrow = 3),
    incoming_mt_profiles = matrix(1:6, nrow = 3)
  )
)
zero_recovery <- get_profiles_ints_and_umis(
  cell_pop = recovery_population,
  bc_or_mt = 'bc',
  int_rec_prob = 0,
  num_ints = 3,
  umis = c('u1', 'u2', 'u3')
)[[1]]
expect_equal(
  dim(zero_recovery$mut_mat),
  c(0L, 2L),
  'A zero recovery probability must return a zero-row mutation matrix.'
)
expect_equal(
  zero_recovery$which_ints_recovered,
  integer(0),
  'A zero recovery probability must not force an observed integration.'
)
expect_equal(
  zero_recovery$recovered_umis,
  character(0),
  'No UMI may be reported when no integration was recovered.'
)

#' All-zero recovery writes an empty FASTA instead of failing
#'
#' Stubs the worker cluster with a serial `parLapply` and replaces `write.fasta`
#' with a function that errors, so any attempt to emit a record is caught, then
#' writes the all-empty population. Asserts the `_TERM.fasta` file exists and is
#' zero bytes.
#'
#' @section Side effects: Writes and then removes one FASTA file under
#'   `tempdir()`.
recovery_env$one_cluster <- NULL
recovery_env$parLapply <- function(cl, X, fun) lapply(X, fun)
recovery_env$write.fasta <- function(...){
  stop('write.fasta must not be called when every profile is empty.')
}
empty_fasta_stem <- tempfile(fileext = '.fasta')
empty_fasta_path <- sub(
  '\\.fasta$',
  '_TERM.fasta',
  empty_fasta_stem
)
suppressWarnings(write_all_cell_sequences(
  cell_mutmats = list(cell1 = zero_recovery$mut_mat),
  reference = c('A', 'C'),
  output_fasta_name = empty_fasta_stem,
  fasta_type = 'TERM',
  bc_integration_umis = list(cell1 = character())
))
expect_true(
  file.exists(empty_fasta_path) && file.info(empty_fasta_path)$size == 0,
  'All-zero recovery must produce an empty FASTA instead of failing.'
)
unlink(empty_fasta_path)

# ---- PhysiCell lineage import and dual-modality replay ----

# PhysiCell persistent-ID divisions must become an event-resolved binary tree,
# and recording mutations must be inherited along that tree.
physicell_env <- new.env(parent = globalenv())
sys.source(file.path(repo_root, 'origin', 'R', 'prime_editing.R'), envir = physicell_env)
sys.source(file.path(repo_root, 'origin', 'R', 'physicell_lineage.R'), envir = physicell_env)
sys.source(file.path(repo_root, 'origin', 'R', 'gillespie_lineage.R'), envir = physicell_env)
sys.source(file.path(repo_root, 'analysis', 'lineage_benchmark.R'), envir = physicell_env)
sys.source(file.path(repo_root, 'analysis', 'physicell_visium.R'), envir = physicell_env)

#' Marker-specific non-Mendelian selection coefficients override the global one
#'
#' Configures a global coefficient of 0.2 with an ecDNA-label override of 0.4.
#' Asserts the labeled marker resolves to its override, that an unlisted marker
#' falls back to the global value, and that a coefficient above one is rejected.
selection_params <- list(
  non_mendelian_selection = list(
    coefficient = 0.2,
    ecdna_label_coefficient = 0.4
  )
)
expect_equal(
  physicell_env$non_mendelian_selection_coefficient(
    selection_params,
    'ecdna_label_coefficient'
  ),
  0.4,
  'A marker-specific non-Mendelian selection coefficient must override the global coefficient.'
)
expect_equal(
  physicell_env$non_mendelian_selection_coefficient(
    selection_params,
    'mitochondrial_variant_coefficient'
  ),
  0.2,
  'Unspecified non-Mendelian markers must use the global selection coefficient.'
)
expect_error(
  physicell_env$non_mendelian_selection_coefficient(
    list(non_mendelian_selection = list(coefficient = 1.01)),
    'mitochondrial_variant_coefficient'
  ),
  'Non-Mendelian selection coefficients above one must be rejected.'
)

#' PhysiCell CSV output is gzipped by default and read back transparently
#'
#' Writes a two-row frame through the shared CSV writer and asserts the returned
#' path is a real `.csv.gz` with no plain `.csv` left behind, that the reader
#' round-trips the frame from the compressed file, that `compress = FALSE`
#' retains the legacy `.csv` path, and that the path resolver prefers the
#' compressed form when both exist.
#'
#' @section Side effects: Writes and then removes `<tempfile>.csv` and
#'   `<tempfile>.csv.gz`.
csv_fixture <- data.frame(
  sample_id = c('cell_1', 'cell_2'),
  event_time = c(0.25, 1.5),
  stringsAsFactors = FALSE
)
csv_base_path <- tempfile('physicell_csv_', fileext = '.csv')
compressed_csv_path <- physicell_env$write_physicell_csv(
  csv_fixture,
  csv_base_path
)
expect_true(
  identical(compressed_csv_path, normalizePath(
    paste0(csv_base_path, '.gz'),
    mustWork = TRUE
  )) &&
    !file.exists(csv_base_path),
  'PhysiCell CSV output must use a real .csv.gz file by default.'
)
expect_equal(
  physicell_env$read_physicell_csv(csv_base_path),
  csv_fixture,
  'PhysiCell CSV readers must resolve and read compressed output transparently.'
)
plain_csv_path <- physicell_env$write_physicell_csv(
  csv_fixture,
  csv_base_path,
  compress = FALSE
)
expect_true(
  identical(plain_csv_path, normalizePath(csv_base_path, mustWork = TRUE)),
  'The CSV compression opt-out must retain the legacy .csv path.'
)
expect_equal(
  physicell_env$resolve_physicell_csv_path(csv_base_path),
  paste0(csv_base_path, '.gz'),
  'Compressed CSV output must be preferred when both file forms exist.'
)
unlink(c(csv_base_path, paste0(csv_base_path, '.gz')))

#' The R timing recorder preserves phase and total wall-clock durations
#'
#' Drives the recorder with a stubbed clock so the durations are exact rather
#' than machine-dependent, finishing two phases and then taking the summary.
#' Asserts the elapsed seconds are 2.5, 2.5 and 6, and that the rendered summary
#' lines align the phase labels and format each duration to one decimal.
load_named_function(
  file.path('origin', 'inst', 'scripts', 'simulate_physicell_lineage.R'),
  'new_physicell_timing_recorder',
  physicell_env
)
load_named_function(
  file.path('origin', 'inst', 'scripts', 'simulate_physicell_lineage.R'),
  'format_physicell_timing_summary',
  physicell_env
)

#' Fixture: a deterministic wall-clock stub for the timing recorder
#'
#' Yields 12.5, 15 and 16 on successive calls, which against a `total_start` of
#' 10 gives two 2.5-second phases inside a 6-second total.
#'
#' @return A zero-argument function returning the next stubbed clock reading.
fake_clock <- local({
  values <- c(12.5, 15, 16)
  value_index <- 0L
  function(){
    value_index <<- value_index + 1L
    values[value_index]
  }
})
timing_recorder <- physicell_env$new_physicell_timing_recorder(
  total_start = 10,
  clock = fake_clock
)
timing_recorder$finish_phase('Lineage simulation')
timing_recorder$finish_phase('Output')
timing_fixture <- timing_recorder$summary()
expect_equal(
  timing_fixture$elapsed_seconds,
  c(2.5, 2.5, 6),
  'The R timing recorder must preserve phase and total wall-clock durations.'
)
expect_equal(
  physicell_env$format_physicell_timing_summary(timing_fixture),
  c(
    'R simulator timing summary:',
    '  Lineage simulation:  2.5s',
    '  Output:              2.5s',
    '  Total R simulator:   6.0s'
  ),
  'The final timing summary must render aligned phase and total durations.'
)

#' PhysiCell persistent-ID divisions become an event-resolved binary tree
#'
#' Reads the fixture division log - parent 0 divides at t=1 and again at t=2,
#' and the resulting cell 1 divides at t=2.5 - with an end time of 3, plus the
#' fixture live-cell table listing IDs 0, 2 and 3. Asserts the reconstruction
#' has seven nodes of which four are terminal, that a parent ID retained across
#' repeated divisions produces a continuing-parent tip, and that the validated
#' division frame is marked so it is not revalidated on reuse.
physicell_fixture_dir <- file.path(repo_root, 'tests', 'fixtures')
physicell_divisions <- physicell_env$read_physicell_divisions(
  file.path(physicell_fixture_dir, 'physicell_divisions.csv')
)
physicell_live_ids <- physicell_env$read_physicell_terminal_ids(
  file.path(physicell_fixture_dir, 'physicell_live_cells.csv')
)
physicell_nodes <- physicell_env$build_physicell_lineage(
  physicell_divisions,
  end_time = 3
)
expect_equal(
  c(nrow(physicell_nodes), sum(physicell_nodes$is_terminal)),
  c(7L, 4L),
  paste(
    'Three persistent-parent PhysiCell divisions must become three internal',
    'segments and four terminal event-resolved branches.'
  )
)
expect_equal(
  physicell_nodes$physicell_id[physicell_nodes$is_terminal],
  c('0', '2', '1', '3'),
  'Repeated division by a retained PhysiCell parent ID must create a continuing-parent tip.'
)
expect_true(
  isTRUE(attr(physicell_divisions, 'physicell_validated')),
  'Validated PhysiCell divisions must be marked to avoid duplicate validation.'
)

#' Each mutation-event column marks exactly its sampled terminal descendants
#'
#' Places three literal events on a founder branch, an internal branch and a
#' terminal branch of the fixture lineage, then builds the event-by-cell
#' descendant matrix. Asserts the incidence pattern for the four terminal cells
#' and that the reported descendant counts equal the matrix column sums.
event_descendant_terminal_nodes <- physicell_nodes[
  physicell_nodes$is_terminal,
  ,
  drop = FALSE
]
event_descendant_fixture <- data.frame(
  node_id = physicell_nodes$node_id[c(1, 2, 7)],
  event = c('founder_edit', 'internal_edit', 'terminal_edit'),
  stringsAsFactors = FALSE
)
event_descendant_result <-
  physicell_env$physicell_event_descendant_matrix(
    physicell_nodes,
    event_descendant_terminal_nodes,
    event_descendant_fixture,
    event_ids = c('event_a', 'event_b', 'event_c')
  )
expected_event_descendants <- matrix(
  c(
    1, 1, 0,
    1, 1, 0,
    1, 0, 0,
    1, 0, 1
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(
    c('cell_0', 'cell_2', 'cell_1', 'cell_3'),
    c('event_a', 'event_b', 'event_c')
  )
)
expect_equal(
  as.matrix(event_descendant_result$matrix),
  expected_event_descendants,
  paste(
    'Each literal mutation-event column must mark exactly the sampled',
    'terminal descendants of the branch where that event arose.'
  )
)
expect_equal(
  event_descendant_result$descendant_counts,
  c(4L, 2L, 1L),
  'Mutation-event descendant counts must match the matrix column sums.'
)

#' The event-descendant writer keeps a one-to-one column-to-event mapping
#'
#' Splits the same three events across a `barcode` and a `mitochondrial`
#' modality and writes the outputs. Asserts the sparse matrix and the compressed
#' manifest both exist, that the matrix column names match the manifest
#' `event_id` order, and that the manifest records each event's modality.
#'
#' @section Side effects: Writes a temporary output directory containing
#'   `mutation_event_descendant_matrix_sparse.rds` and
#'   `mutation_event_descendant_manifest.csv.gz`.
event_descendant_output_path <- tempfile('event_descendant_output_')
event_descendant_output <- physicell_env$write_physicell_event_descendant_outputs(
  physicell_nodes,
  event_descendant_terminal_nodes,
  list(
    barcode = event_descendant_fixture[1:2, , drop = FALSE],
    mitochondrial = event_descendant_fixture[3, , drop = FALSE]
  ),
  event_descendant_output_path
)
written_event_descendant_manifest <- physicell_env$read_physicell_csv(
  file.path(
    event_descendant_output_path,
    'mutation_event_descendant_manifest.csv'
  )
)
expect_true(
  file.exists(file.path(
    event_descendant_output_path,
    'mutation_event_descendant_matrix_sparse.rds'
  )) &&
    file.exists(file.path(
      event_descendant_output_path,
      'mutation_event_descendant_manifest.csv.gz'
    )) &&
    identical(
      colnames(event_descendant_output$matrix),
      written_event_descendant_manifest$event_id
    ) &&
    identical(
      written_event_descendant_manifest$modality,
      c('barcode', 'barcode', 'mitochondrial')
    ),
  paste(
    'The event-descendant writer must preserve a one-to-one mapping from',
    'matrix columns to modality-specific mutation-event rows.'
  )
)

#' A daughter ID that already existed as a parent must be rejected
#'
#' Feeds the validator a division log in which ID 0 acts as a parent at t=1 and
#' is then re-emitted as a daughter at t=2, which would silently corrupt the
#' persistent-ID lineage.
expect_error(
  physicell_env$validate_physicell_divisions(data.frame(
    time = c(1, 2),
    parent_ID = c(0, 1),
    daughter_ID = c(1, 0)
  )),
  'A daughter ID that already existed as a parent must still be rejected.'
)

#' A 10,000-division retained-parent lineage stays exact and linear-time
#'
#' Builds a lineage in which one persistent parent divides 10,000 times. Asserts
#' the reconstruction has `2n + 1` nodes with `n + 1` terminal branches and
#' completes well under ten seconds, which fails if validation or reconstruction
#' is superlinear.
large_division_count <- 10000L
large_physicell_divisions <- data.frame(
  time = seq_len(large_division_count),
  parent_ID = rep(0L, large_division_count),
  daughter_ID = seq_len(large_division_count)
)
large_lineage_elapsed <- system.time({
  large_physicell_nodes <- physicell_env$build_physicell_lineage(
    large_physicell_divisions,
    end_time = large_division_count + 1
  )
})[['elapsed']]
expect_true(
  nrow(large_physicell_nodes) == 2L * large_division_count + 1L &&
    sum(large_physicell_nodes$is_terminal) == large_division_count + 1L,
  'A large retained-parent lineage must preserve its exact binary-tree dimensions.'
)
expect_true(
  large_lineage_elapsed < 10,
  paste(
    'A 10,000-division lineage should complete comfortably with linear-time',
    'validation and reconstruction.'
  )
)

#' The Newick writer streams a 10,000-generation lineage without recursion
#'
#' Renders the same deep retained-parent lineage to a file. Asserts the file is
#' written non-empty and inside the ten-second budget, which a recursive writer
#' or one that re-renders subtrees would not meet.
#'
#' @section Side effects: Writes and then removes one `.nwk` file under
#'   `tempdir()`.
large_newick_path <- tempfile(fileext = '.nwk')
large_newick_elapsed <- system.time(
  physicell_env$write_physicell_lineage_newick(
    large_physicell_nodes,
    large_newick_path
  )
)[['elapsed']]
expect_true(
  file.exists(large_newick_path) &&
    file.info(large_newick_path)$size > 0 &&
    large_newick_elapsed < 10,
  paste(
    'The iterative Newick writer must handle a retained-parent lineage with',
    '10,000 generations without recursion or superlinear subtree rendering.'
  )
)
unlink(large_newick_path)

#' Lineage progress is throttled and reports count, percent, rate, and ETA
#'
#' Requests two updates over ten nodes and calls the reporter at 0, 1, 5 and 10.
#' Asserts exactly three messages are emitted - the throttled call at 1 is
#' suppressed - and that they carry the starting notice, the halfway mark, and a
#' completion line with a nodes/second rate and an ETA.
progress_messages <- capture.output(
  {
    progress_reporter <- physicell_env$new_physicell_progress_reporter(
      total = 10,
      label = 'test recorder',
      enabled = TRUE,
      updates = 2
    )
    progress_reporter(0)
    progress_reporter(1)
    progress_reporter(5)
    progress_reporter(10)
  },
  type = 'message'
)
expect_true(
  length(progress_messages) == 3 &&
    grepl('0/10 nodes \\(0.0%\\): starting', progress_messages[1]) &&
    grepl('5/10 nodes \\(50.0%\\)', progress_messages[2]) &&
    grepl('10/10 nodes \\(100.0%\\)', progress_messages[3]) &&
    grepl('nodes/s, ETA', progress_messages[3]),
  'Lineage progress must be throttled and report counts, percent, rate, and ETA.'
)

#' Disabled lineage progress emits nothing
#'
#' Repeats the reporter with `enabled = FALSE` and asserts the message stream is
#' empty.
quiet_progress <- capture.output(
  {
    progress_reporter <- physicell_env$new_physicell_progress_reporter(
      total = 10,
      label = 'quiet recorder',
      enabled = FALSE,
      updates = 2
    )
    progress_reporter(0)
    progress_reporter(10)
  },
  type = 'message'
)
expect_equal(
  quiet_progress,
  character(),
  'Disabled lineage progress must not emit messages.'
)

#' Stage logging is timestamped, and silent when disabled
#'
#' Asserts an enabled stage log emits one `[YYYY-MM-DD ...] <message>` line and
#' that a disabled one emits nothing.
stage_messages <- capture.output(
  physicell_env$physicell_log_stage(
    'Test output phase complete.',
    enabled = TRUE
  ),
  type = 'message'
)
expect_true(
  length(stage_messages) == 1 &&
    grepl(
      '^\\[[0-9]{4}-[0-9]{2}-[0-9]{2} .+\\] Test output phase complete\\.$',
      stage_messages
    ),
  'Stage logging must include a timestamp and the supplied phase transition.'
)
quiet_stage_messages <- capture.output(
  physicell_env$physicell_log_stage(
    'Hidden output phase.',
    enabled = FALSE
  ),
  type = 'message'
)
expect_equal(
  quiet_stage_messages,
  character(),
  'Disabled stage logging must remain quiet.'
)

#' scDesign3 stage logging emits the supplied timestamped transition
#'
#' Sources the scDesign3 helpers into their own environment and asserts one
#' message is emitted, ending in the supplied phase text.
scdesign3_helper_env <- new.env(parent = globalenv())
sys.source(
  file.path(repo_root, 'analysis', 'scdesign3_helpers.R'),
  envir = scdesign3_helper_env
)
scdesign3_stage_messages <- capture.output(
  scdesign3_helper_env$scdesign3_log_stage(
    'scDesign3 test phase complete.'
  ),
  type = 'message'
)
expect_true(
  length(scdesign3_stage_messages) == 1 &&
    grepl(
      'scDesign3 test phase complete\\.$',
      scdesign3_stage_messages
    ),
  'scDesign3 stage logging must emit the supplied timestamped transition.'
)

#' Recording replay filters to live cells and inherits founder-branch edits
#'
#' Prepares a recording model from the fixture parameter JSON, then replays it
#' over the fixture lineage while passing the live-cell IDs. Asserts the
#' terminal set is narrowed to the imported live IDs `0`, `2`, `3`, and that an
#' A-to-G edit arising on the founder branch appears in every sampled
#' descendant's profile.
physicell_params <- jsonlite::fromJSON(
  file.path(physicell_fixture_dir, 'physicell_params.json'),
  simplifyVector = FALSE
)
physicell_model <- physicell_env$prepare_physicell_recording_model(
  physicell_params,
  params_dir = physicell_fixture_dir,
  seed = 7
)
physicell_simulation <- physicell_env$simulate_recording_on_physicell_lineage(
  physicell_nodes,
  physicell_model,
  terminal_physicell_ids = physicell_live_ids,
  seed = 7
)
expect_equal(
  physicell_simulation$terminal_nodes$physicell_id,
  c('0', '2', '3'),
  'The optional PhysiCell live-cell table must filter extinct terminal branches.'
)
for(node_id in physicell_simulation$terminal_nodes$node_id){
  expect_true(
    all(physicell_simulation$profiles[[node_id]][, 2] == 2),
    'A founder-branch A-to-G recording edit must be inherited by every sampled descendant.'
  )
}

#' Newick rendering matches for the full tree, the sampled tree, and the stream
#'
#' Compares the in-memory renderer against literal expected strings for the full
#' lineage and for the live-cell-pruned lineage, then writes the pruned tree
#' through the streaming writer with progress on. Asserts pruning removes only
#' unrequested tips while retaining unary ancestors and elapsed branch lengths,
#' that the streamed file matches the in-memory string exactly, and that a
#' completion message is reported.
#'
#' @section Side effects: Writes and then removes one `.nwk` file under
#'   `tempdir()`.
expected_full_physicell_newick <- paste0(
  '(((cell_0:1,cell_2:1)node_000002:1,',
  '(cell_1:0.5,cell_3:0.5)node_000003:1.5)',
  'node_000001:1)remote_mito_root;'
)
expected_sampled_physicell_newick <- paste0(
  '(((cell_0:1,cell_2:1)node_000002:1,',
  '(cell_3:0.5)node_000003:1.5)',
  'node_000001:1)remote_mito_root;'
)
expect_equal(
  physicell_env$physicell_lineage_to_newick(physicell_nodes),
  expected_full_physicell_newick,
  'Iterative Newick rendering must preserve the established full-tree syntax.'
)
expect_equal(
  physicell_env$physicell_lineage_to_newick(
    physicell_nodes,
    terminal_physicell_ids = physicell_live_ids
  ),
  expected_sampled_physicell_newick,
  paste(
    'Sampled Newick rendering must prune unrequested tips while retaining',
    'unary ancestors and elapsed branch lengths.'
  )
)
streamed_newick_path <- tempfile(fileext = '.nwk')
newick_progress <- capture.output(
  physicell_env$write_physicell_lineage_newick(
    physicell_nodes,
    streamed_newick_path,
    terminal_physicell_ids = physicell_live_ids,
    show_progress = TRUE,
    progress_updates = 2
  ),
  type = 'message'
)
expect_equal(
  readLines(streamed_newick_path, warn = FALSE),
  expected_sampled_physicell_newick,
  'Buffered Newick output must match the compatibility in-memory renderer.'
)
expect_true(
  any(grepl('Newick rendering.*100.0%', newick_progress)),
  'Streaming Newick output must report completion when progress is enabled.'
)
unlink(streamed_newick_path)

#' Barcode mutation events carry exact times inside their own branch
#'
#' Asserts every recorded event has an `event_time` between its branch start and
#' branch end, i.e. events are placed in continuous time on the branch that
#' produced them rather than snapped to a node.
expect_true(
  nrow(physicell_simulation$mutation_events) > 0 &&
    all(
      physicell_simulation$mutation_events$event_time >=
        physicell_simulation$mutation_events$branch_start
    ) &&
    all(
      physicell_simulation$mutation_events$event_time <=
        physicell_simulation$mutation_events$branch_end
    ),
  'Dynamic barcode events must have exact times inside their lineage branches.'
)

# ---- Mitochondrial replay on PhysiCell lineages ----

sys.source(file.path(repo_root, 'origin', 'R', 'physicell_mito.R'), envir = physicell_env)
physicell_mito_model <- physicell_env$prepare_physicell_mito_model(
  physicell_params,
  genomes_per_cell = 4,
  seed = 7
)

#' Mitochondrial daughter sampling: neutral path and full purifying selection
#'
#' Builds a four-genome parent profile with two wild-type genomes and two
#' carrying a variant. Under `variant_selection_coefficient = 0` the daughter
#' draw must reproduce, seed for seed, the legacy `sample(..., replace = TRUE)`
#' path computed independently in the test. Under a coefficient of one, a
#' 100-genome bottleneck must contain no variant-bearing genome while wild-type
#' genomes remain available.
mitochondrial_selection_fixture <- list(
  setNames(numeric(), character()),
  setNames(numeric(), character()),
  setNames(2, '10'),
  setNames(3, '11')
)
set.seed(23L)
mitochondrial_neutral_expected_indices <- sample(
  seq_along(mitochondrial_selection_fixture),
  25L,
  replace = TRUE
)
set.seed(23L)
mitochondrial_neutral_profile <- physicell_env$inherit_physicell_mito_profile(
  mitochondrial_selection_fixture,
  genomes_per_cell = 25L,
  variant_selection_coefficient = 0
)
expect_equal(
  mitochondrial_neutral_profile,
  mitochondrial_selection_fixture[mitochondrial_neutral_expected_indices],
  'Neutral mitochondrial selection must preserve the legacy seeded sampling path.'
)
mitochondrial_selected_profile <- physicell_env$inherit_physicell_mito_profile(
  mitochondrial_selection_fixture,
  genomes_per_cell = 100L,
  variant_selection_coefficient = 1
)
expect_true(
  all(lengths(mitochondrial_selected_profile) == 0L),
  paste(
    'A mitochondrial selection coefficient of one must exclude variant-bearing',
    'genomes when wild-type genomes are available at the daughter bottleneck.'
  )
)

#' Mitochondrial replay shares the live-cell filter, bottleneck, and event times
#'
#' Replays the mitochondrial model over the same fixture lineage. Asserts the
#' terminal set matches the barcode replay's live-cell filter, that every
#' imported terminal cell carries exactly the configured number of genomes, that
#' mutation events fall inside their branch intervals, and that per-cell
#' heteroplasmy fractions land in (0, 1].
physicell_mito_simulation <- physicell_env$simulate_mito_on_physicell_lineage(
  physicell_nodes,
  physicell_mito_model,
  physicell_params,
  terminal_physicell_ids = physicell_live_ids,
  seed = 7
)
expect_equal(
  physicell_mito_simulation$terminal_nodes$physicell_id,
  c('0', '2', '3'),
  'Mitochondrial replay must use the same PhysiCell live-cell filter as barcode replay.'
)
expect_true(
  all(vapply(
    physicell_mito_simulation$profiles[
      physicell_mito_simulation$terminal_nodes$node_id
    ],
    length,
    integer(1)
  ) == 4),
  'Every imported terminal cell must retain the configured mitochondrial bottleneck size.'
)
expect_true(
  nrow(physicell_mito_simulation$mutation_events) > 0 &&
    all(
      physicell_mito_simulation$mutation_events$event_time >=
        physicell_mito_simulation$mutation_events$branch_start
    ) &&
    all(
      physicell_mito_simulation$mutation_events$event_time <=
        physicell_mito_simulation$mutation_events$branch_end
    ),
  'Dynamic mitochondrial events must have exact times inside their lineage branches.'
)
physicell_mito_fractions <- physicell_env$physicell_mito_variant_fractions(
  physicell_mito_simulation
)
expect_true(
  nrow(physicell_mito_fractions) > 0 &&
    all(physicell_mito_fractions$variant_fraction > 0) &&
    all(physicell_mito_fractions$variant_fraction <= 1),
  'Mitochondrial replay must emit valid per-cell heteroplasmy fractions.'
)

#' PhysiCell scDesign3 covariates join lineage, spatial, and recorder features
#'
#' Writes the barcode and mitochondrial outputs plus a sparse binary score
#' matrix into one directory, then builds the covariate frame against the
#' fixture live-cell table. Asserts sample IDs match the shared PhysiCell
#' terminal IDs, that lineage, spatial, barcode and mitochondrial covariates are
#' all present, that lineage pseudotime stays in the unit interval, and that
#' barcode edit counts are read from the compact sparse BASELINE output.
#'
#' @section Side effects: Creates a temporary directory of PhysiCell outputs
#'   that later blocks reuse; it is removed after the organoid covariate test.
sys.source(file.path(repo_root, 'analysis', 'physicell_scdesign3.R'), envir = physicell_env)
physicell_sc_output <- tempfile('physicell_sc_covariates_')
dir.create(physicell_sc_output)
physicell_barcode_output <- physicell_env$write_physicell_recording_outputs(
  physicell_simulation,
  physicell_model,
  physicell_sc_output
)
physicell_env$write_physicell_mito_outputs(
  physicell_mito_simulation,
  physicell_mito_model,
  physicell_sc_output
)
saveRDS(
  Matrix::Matrix(physicell_barcode_output$binary_scores, sparse = TRUE),
  file.path(
    physicell_sc_output,
    'barcode_binary_score_matrix_sparse.rds'
  )
)
physicell_sc_metadata <- physicell_env$build_physicell_sc_covariates(
  physicell_sc_output,
  file.path(physicell_fixture_dir, 'physicell_live_cells.csv')
)
expect_equal(
  physicell_sc_metadata$sample_id,
  paste0('cell_', physicell_live_ids),
  'scDesign3 covariates must preserve the shared PhysiCell terminal sample IDs.'
)
expect_true(
  all(c(
    'lineage_depth', 'lineage_pseudotime', 'x', 'y', 'z',
    'neighbor_count', 'barcode_edit_count', 'mt_variant_count',
    'mt_heteroplasmy_burden'
  ) %in% names(physicell_sc_metadata)),
  'PhysiCell scDesign3 metadata must join lineage, spatial, barcode, and mitochondrial covariates.'
)
expect_true(
  all(is.finite(physicell_sc_metadata$lineage_pseudotime)) &&
    all(physicell_sc_metadata$lineage_pseudotime >= 0) &&
    all(physicell_sc_metadata$lineage_pseudotime <= 1),
  'PhysiCell lineage pseudotime must remain in the unit interval.'
)
expect_true(
  all(physicell_sc_metadata$barcode_edit_count > 0),
  'scDesign3 covariates must read compact sparse BASELINE barcode output.'
)

# Multi-founder organoid imports must retain undivided founders and may assign
# stable allele-coded founder labels before branch recording begins.

#' Undivided founders survive import and receive stable founder labels
#'
#' Imports one division under two declared founders, so founder 1 never divides.
#' Asserts multi-founder Newick keeps the divided clade first and the undivided
#' founder as its own root-level tip, that all three IDs are terminal, and that
#' a founder-label site is inherited within a founder clone while differing
#' between founders.
organoid_nodes <- physicell_env$build_physicell_lineage(
  data.frame(time = 1, parent_ID = 0, daughter_ID = 2),
  end_time = 2,
  founder_ids = c(0, 1)
)
expect_equal(
  physicell_env$physicell_lineage_to_newick(organoid_nodes),
  paste0(
    '((cell_0:1,cell_2:1)node_000001:1,',
    'cell_1:2)remote_mito_root;'
  ),
  'Iterative Newick rendering must preserve multi-founder root ordering.'
)
expect_equal(
  sort(organoid_nodes$physicell_id[organoid_nodes$is_terminal]),
  c('0', '1', '2'),
  'An explicitly supplied undivided founder must remain a terminal lineage root.'
)
organoid_model <- physicell_env$prepare_physicell_recording_model(
  physicell_params,
  founder_label_sites = 1,
  params_dir = physicell_fixture_dir,
  seed = 9
)
organoid_recording <- physicell_env$simulate_recording_on_physicell_lineage(
  organoid_nodes,
  organoid_model,
  seed = 9
)
organoid_terminal_profiles <- setNames(
  organoid_recording$profiles[organoid_recording$terminal_nodes$node_id],
  organoid_recording$terminal_nodes$physicell_id
)
founder_label_position <- organoid_model$founder_label_positions[1]
expect_true(
  organoid_terminal_profiles[['0']][1, founder_label_position] ==
    organoid_terminal_profiles[['2']][1, founder_label_position] &&
    organoid_terminal_profiles[['0']][1, founder_label_position] !=
      organoid_terminal_profiles[['1']][1, founder_label_position],
  'Founder barcodes must be inherited within a founder clone and differ across founders.'
)

# ---- BASELINE recorder preset ----

#' The BASELINE preset defines five sparse recorders with 50 targets each
#'
#' Prepares a recording model from the neural-organoid BASELINE parameter file
#' with five integrations and no founder-label sites. Asserts the target count,
#' integration count and sparse profile storage, and that rate preparation
#' caches one event hazard and one total hazard per barcode position, with
#' positive totals at every active position.
baseline_params <- jsonlite::fromJSON(
  file.path(
    repo_root,
    'example_json_params',
    'physicell_neural_organoid_baseline.json'
  ),
  simplifyVector = FALSE
)
baseline_model <- physicell_env$prepare_physicell_recording_model(
  baseline_params,
  num_integrations = 5,
  founder_label_sites = 0,
  seed = 2
)
expect_true(
  length(baseline_model$be_targets) == 50 &&
    baseline_model$num_integrations == 5 &&
    identical(baseline_model$profile_storage, 'sparse'),
  'The BASELINE preset must define five sparse recorders with 50 targets each.'
)
baseline_induced_rates <-
  baseline_model$rate_sets[['induced_editing_params']]
expect_true(
  length(baseline_induced_rates$event_hazards) ==
    baseline_model$barcode_length &&
    length(baseline_induced_rates$total_hazards) ==
      baseline_model$barcode_length &&
    all(
      baseline_induced_rates$total_hazards[
        baseline_induced_rates$active_positions
      ] > 0
    ),
  'BASELINE rate preparation must cache per-position event and total hazards.'
)

#' Sparse BASELINE replay records events, frees internal profiles, and repeats
#'
#' Replays the BASELINE model over the fixture lineage with a long end time,
#' then replays it again under the same seed. Asserts mutations are recorded,
#' that the retained profiles are exactly the terminal nodes - internal profiles
#' having been released - and that the second run reproduces the first run's
#' event table exactly.
baseline_nodes <- physicell_env$build_physicell_lineage(
  physicell_divisions,
  end_time = 30240
)
baseline_recording <- physicell_env$simulate_recording_on_physicell_lineage(
  baseline_nodes,
  baseline_model,
  seed = 2
)
expect_true(
  nrow(baseline_recording$mutation_events) > 0 &&
    setequal(
      names(baseline_recording$profiles),
      baseline_recording$terminal_nodes$node_id
    ),
  paste(
    'Sparse BASELINE replay must record mutations while releasing',
    'non-terminal profiles.'
  )
)
baseline_recording_repeat <-
  physicell_env$simulate_recording_on_physicell_lineage(
    baseline_nodes,
    baseline_model,
    seed = 2
  )
expect_equal(
  baseline_recording_repeat$mutation_events,
  baseline_recording$mutation_events,
  'Vectorized BASELINE replay must remain reproducible for a fixed seed.'
)

#' The growing mutation-event buffer keeps every event past its initial size
#'
#' Replays eight founders with no divisions over a very long branch, so every
#' active position of every integration must fire and the event count exceeds
#' the columnar buffer's initial allocation. Asserts the recorded event count
#' equals the analytically expected total.
buffer_growth_nodes <- physicell_env$build_physicell_lineage(
  data.frame(
    time = numeric(),
    parent_ID = integer(),
    daughter_ID = integer()
  ),
  end_time = 1e14,
  founder_ids = 0:7
)
buffer_growth_recording <-
  physicell_env$simulate_recording_on_physicell_lineage(
    buffer_growth_nodes,
    baseline_model,
    seed = 4
  )
expected_buffer_events <- 8L * baseline_model$num_integrations *
  length(baseline_induced_rates$active_positions)
expect_true(
  expected_buffer_events > 4096L &&
    nrow(buffer_growth_recording$mutation_events) ==
      expected_buffer_events,
  paste(
    'The growing columnar mutation-event buffer must preserve every event',
    'after exceeding its initial allocation.'
  )
)

#' A 1,000-founder BASELINE replay uses vectorized coordinate draws
#'
#' Replays 1,000 founders with no divisions and asserts all 1,000 profiles come
#' back inside a ten-second budget, which a scalar loop over integrations and
#' positions would not meet.
vectorized_recording_nodes <- physicell_env$build_physicell_lineage(
  data.frame(
    time = numeric(),
    parent_ID = integer(),
    daughter_ID = integer()
  ),
  end_time = 1,
  founder_ids = 0:999
)
vectorized_recording_elapsed <- system.time(
  vectorized_recording <- physicell_env$simulate_recording_on_physicell_lineage(
    vectorized_recording_nodes,
    baseline_model,
    seed = 5
  )
)[['elapsed']]
expect_true(
  length(vectorized_recording$profiles) == 1000L &&
    vectorized_recording_elapsed < 10,
  paste(
    'A 1,000-founder BASELINE replay must use vectorized coordinate draws',
    'rather than scalar integration-position loops.'
  )
)

#' The sequential fallback still enforces close-after-edit windows
#'
#' Turns on `close_be_window` and mutates one profile over a very long duration,
#' which forces the sequential path rather than the vectorized draw. Asserts no
#' integration ever holds more than one edited position within any single BE
#' window.
closed_window_model <- baseline_model
closed_window_model$close_be_window <- TRUE
closed_window_mutation <- physicell_env$mutate_physicell_barcode_segment(
  physicell_env$initialize_physicell_barcode_profile(closed_window_model),
  duration = 1e14,
  rate_set = baseline_induced_rates,
  model = closed_window_model
)
be_window_names <- names(closed_window_model$windows)[
  grepl('^be_window_', names(closed_window_model$windows))
]
expect_true(
  nrow(closed_window_mutation$events) > 0 &&
    all(vapply(
      seq_len(closed_window_model$num_integrations),
      function(integration){
        edited_positions <- as.integer(names(
          closed_window_mutation$profile[[integration]]
        ))
        all(vapply(
          be_window_names,
          function(window_name){
            sum(
              edited_positions %in%
                closed_window_model$windows[[window_name]]
            ) <= 1L
          },
          logical(1)
        ))
      },
      logical(1)
    )),
  'The sequential fallback must still enforce close-after-edit windows.'
)

#' Compact BASELINE output writes exact alleles, sparse matrices, and progress
#'
#' Writes the BASELINE replay to a temporary directory. Asserts the binary score
#' matrix is a `Matrix` object, that the target layout has 250 rows (five
#' integrations of 50 targets), and that the compressed layout, event, manifest
#' and sparse-matrix files all exist. The raw allele matrix is then rebuilt
#' position by position from the retained terminal profiles and compared
#' exactly, which pins the preallocated compact triplet assembly. Finally the
#' sparse assembly is re-run with progress enabled and must report completion.
#'
#' @section Side effects: Creates and then removes a temporary BASELINE output
#'   directory.
baseline_output_path <- tempfile('physicell_baseline_sparse_')
dir.create(baseline_output_path)
baseline_output <- physicell_env$write_physicell_recording_outputs(
  baseline_recording,
  baseline_model,
  baseline_output_path
)
baseline_layout <- physicell_env$read_physicell_csv(
  file.path(baseline_output_path, 'barcode_target_layout.csv')
)
expect_true(
  inherits(baseline_output$binary_scores, 'Matrix') &&
    nrow(baseline_layout) == 250 &&
    file.exists(file.path(
      baseline_output_path,
      'barcode_target_layout.csv.gz'
    )) &&
    file.exists(file.path(
      baseline_output_path,
      'mutation_events.csv.gz'
    )) &&
    file.exists(file.path(
      baseline_output_path,
      'run_manifest.csv.gz'
    )) &&
    file.exists(file.path(
      baseline_output_path,
      'barcode_binary_score_matrix_sparse.rds'
    )),
  'Compact BASELINE output must contain 250 integrated target definitions and a sparse matrix.'
)
baseline_terminal_profiles <- baseline_recording$profiles[
  baseline_recording$terminal_nodes$node_id
]
names(baseline_terminal_profiles) <- paste0(
  'cell_',
  baseline_recording$terminal_nodes$physicell_id
)
manual_baseline_alleles <- matrix(
  0,
  nrow = length(baseline_terminal_profiles),
  ncol = baseline_model$num_integrations *
    length(baseline_model$output_positions),
  dimnames = dimnames(as.matrix(baseline_output$raw_alleles))
)
for(cell_index in seq_along(baseline_terminal_profiles)){
  profile <- baseline_terminal_profiles[[cell_index]]
  for(integration in seq_len(baseline_model$num_integrations)){
    integration_values <- profile[[integration]]
    retained_positions <- match(
      names(integration_values),
      as.character(baseline_model$output_positions),
      nomatch = 0L
    )
    retained <- retained_positions > 0L
    if(any(retained)){
      manual_baseline_alleles[
        cell_index,
        (integration - 1L) * length(baseline_model$output_positions) +
          retained_positions[retained]
      ] <- unname(integration_values[retained])
    }
  }
}
expect_equal(
  as.matrix(baseline_output$raw_alleles),
  manual_baseline_alleles,
  'Preallocated compact barcode triplets must preserve exact allele values.'
)
sparse_matrix_progress <- capture.output(
  invisible(
    physicell_env$physicell_sparse_recording_matrix(
      baseline_terminal_profiles,
      baseline_model,
      show_progress = TRUE,
      progress_updates = 2
    )
  ),
  type = 'message'
)
expect_true(
  any(grepl('Barcode sparse matrix assembly.*100.0%', sparse_matrix_progress)),
  'Compact barcode matrix assembly must report completion when enabled.'
)
unlink(baseline_output_path, recursive = TRUE)

#' An alive column in the PhysiCell state table excludes dead observations
#'
#' Writes a two-row state table with one dead cell and asserts only the living
#' ID is returned as a terminal ID.
#'
#' @section Side effects: Writes and then removes one CSV under `tempdir()`.
alive_filter_path <- tempfile('physicell_alive_filter_', fileext = '.csv')
utils::write.csv(
  data.frame(ID = c(0, 1), alive = c(TRUE, FALSE)),
  alive_filter_path,
  row.names = FALSE
)
expect_equal(
  physicell_env$read_physicell_terminal_ids(alive_filter_path),
  '0',
  'An alive column in the PhysiCell state table must exclude dead observations.'
)
unlink(alive_filter_path)

#' Organoid covariates use per-cell developmental state and ecDNA summaries
#'
#' Extends the fixture live-cell table with organoid columns - cell type,
#' founder, culture day, time in state, developmental pseudotime, oxygen and
#' nutrient - and writes an ecDNA per-cell summary next to the earlier PhysiCell
#' outputs. Asserts the covariate frame takes its cell types from the PhysiCell
#' state, that developmental pseudotime replaces lineage depth as the expression
#' pseudotime, that founder, developmental, niche and ecDNA covariates are all
#' retained, and that ecDNA copy numbers align to the terminal sample IDs.
#'
#' @section Side effects: Writes a temporary organoid state CSV and an
#'   `ecdna_cell_summary.csv.gz`, then removes both temporary directories.
organoid_state_path <- tempfile('physicell_organoid_state_', fileext = '.csv')
organoid_state <- utils::read.csv(
  file.path(physicell_fixture_dir, 'physicell_live_cells.csv'),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
organoid_state$cell_type <- c(
  'radial_glia', 'neural_progenitor', 'neuron'
)
organoid_state$founder_ID <- c(0, 0, 0)
organoid_state$culture_day <- 30
organoid_state$state_start_time <- c(10, 20, 25)
organoid_state$time_in_state <- c(20, 10, 5)
organoid_state$transition_count <- c(2, 3, 4)
organoid_state$developmental_pseudotime <- c(0.4, 0.65, 0.9)
organoid_state$oxygen <- c(20, 15, 10)
organoid_state$nutrient <- c(1, 0.9, 0.8)
organoid_state$alive <- TRUE
utils::write.csv(organoid_state, organoid_state_path, row.names = FALSE)
physicell_env$write_physicell_csv(
  data.frame(
    sample_id = paste0('cell_', organoid_state$ID),
    total_ecdna_copies = c(12, 18, 7),
    labeled_ecdna_copies = c(4, 5, 2),
    observed_ecdna_species = c(3, 4, 2),
    observed_labeled_species = c(1, 2, 1),
    ecdna_recorder_edit_fraction = c(0.1, 0.2, 0),
    stringsAsFactors = FALSE
  ),
  file.path(physicell_sc_output, 'ecdna_cell_summary.csv'),
  row.names = FALSE,
  compress = TRUE
)
organoid_sc_metadata <- physicell_env$build_physicell_sc_covariates(
  physicell_sc_output,
  organoid_state_path
)
expect_equal(
  organoid_sc_metadata$cell_type,
  organoid_state$cell_type,
  'scDesign3 covariates must use per-cell PhysiCell developmental states.'
)
expect_equal(
  organoid_sc_metadata$pseudotime,
  organoid_state$developmental_pseudotime,
  'Organoid developmental pseudotime must replace lineage depth as expression pseudotime.'
)
expect_true(
  all(c(
    'founder_id', 'culture_day', 'time_in_state', 'oxygen', 'nutrient',
    'radial_position', 'ecdna_copy_number', 'ecdna_labeled_copy_number',
    'ecdna_species_count', 'ecdna_recorder_edit_fraction'
  ) %in% names(organoid_sc_metadata)),
  paste(
    'Organoid scDesign3 metadata must retain founder, developmental, niche,',
    'and ecDNA covariates.'
  )
)
expect_equal(
  organoid_sc_metadata$ecdna_copy_number,
  c(12L, 18L, 7L),
  'scDesign3 ecDNA covariates must align to terminal sample IDs.'
)
unlink(organoid_state_path)
unlink(physicell_sc_output, recursive = TRUE)

# ---- Exact Gillespie population engine ----

# Exact Gillespie population engine.
gillespie_env <- new.env(parent = globalenv())
sys.source(file.path(repo_root, 'origin', 'R', 'gillespie_lineage.R'), envir = gillespie_env)

#' Gillespie population replay is reproducible, ordered, and continuous-time
#'
#' Runs a 200-cell birth-death population twice under the same seed. Asserts the
#' two node tables are identical, that every branch length is non-negative with
#' each parent appearing before its children, and that division and death times
#' are genuinely off-grid rather than snapped to a fixed timestep.
gillespie_test_params <- list(
  num_init_cells = 200L,
  sim_length = list(0.75),
  random_seed = 101L,
  consider_cell_heteroplasmy_scores = FALSE,
  cell_type_dict = list(
    founder_cell_type = 'cycling',
    cell_type_params = list(
      cycling = list(
        cell_cycle_length = 1,
        death_per_cell_cycle_prob = 0.1
      )
    ),
    uninduced_transition_matrix = list(list(1)),
    induced_transition_matrix = list(list(1))
  )
)
gillespie_first <- gillespie_env$simulate_gillespie_population(
  gillespie_test_params,
  seed = 44L,
  show_progress = FALSE
)
gillespie_second <- gillespie_env$simulate_gillespie_population(
  gillespie_test_params,
  seed = 44L,
  show_progress = FALSE
)
expect_equal(
  gillespie_first$nodes,
  gillespie_second$nodes,
  'Gillespie population simulation must be reproducible for a fixed seed.'
)
expect_true(
  all(gillespie_first$nodes$branch_length >= 0) &&
    all(
      is.na(gillespie_first$nodes$parent_node_id) |
        match(
          gillespie_first$nodes$parent_node_id,
          gillespie_first$nodes$node_id
        ) < seq_len(nrow(gillespie_first$nodes))
    ),
  'Gillespie nodes must have non-negative branches in parent-before-child order.'
)
stochastic_times <- gillespie_first$event_log$time[
  gillespie_first$event_log$event %in% c('division', 'death')
]
expect_true(
  length(stochastic_times) > 0 && any(abs(stochastic_times - round(stochastic_times)) > 1e-8),
  'Gillespie population events must occur at continuous, non-grid times.'
)

#' Mean Gillespie growth agrees with the birth-death expectation
#'
#' Grows 10,000 cells for a quarter of a cell cycle with a per-cycle death
#' probability of 0.1, converting that probability to the continuous death rate
#' `-log1p(-0.1)`. Asserts the realized terminal count is within 3% of
#' `n0 * exp((birth - death) * t)`.
gillespie_statistical_params <- gillespie_test_params
gillespie_statistical_params$num_init_cells <- 10000L
gillespie_statistical_params$sim_length <- list(0.25)
gillespie_statistical <- gillespie_env$simulate_gillespie_population(
  gillespie_statistical_params,
  seed = 73L,
  max_cells = 50000L,
  show_progress = FALSE
)
birth_rate <- 1
death_rate <- -log1p(-0.1)
expected_population <- 10000 * exp((birth_rate - death_rate) * 0.25)
expect_true(
  abs(nrow(gillespie_statistical$terminal_nodes) - expected_population) /
    expected_population < 0.03,
  'Gillespie mean population growth must agree with the birth-death expectation.'
)

#' Geometric node-buffer growth initializes every appended row
#'
#' Runs a pure-birth population long enough to force several buffer doublings.
#' Asserts the node table grows past 1,600 rows with no `NA` left in the
#' `is_terminal` or `died` columns, i.e. newly allocated rows are initialized
#' rather than left missing.
gillespie_growth_params <- gillespie_test_params
gillespie_growth_params$num_init_cells <- 400L
gillespie_growth_params$sim_length <- list(1.5)
gillespie_growth_params$cell_type_dict$cell_type_params$cycling$
  death_per_cell_cycle_prob <- 0
gillespie_growth <- gillespie_env$simulate_gillespie_population(
  gillespie_growth_params,
  seed = 11L,
  max_cells = 10000L,
  show_progress = FALSE
)
expect_true(
  nrow(gillespie_growth$nodes) > 1600L &&
    !anyNA(gillespie_growth$nodes$is_terminal) &&
    !anyNA(gillespie_growth$nodes$died),
  'Geometrically grown Gillespie node columns must initialize every appended row.'
)

#' Partial editing induction selects an exact count at an exact time
#'
#' Induces exactly ten of the living cells at t=0.2. Asserts ten
#' `editing_induction` events are logged and that every continuation branch
#' created by the induction is born precisely at 0.2.
gillespie_induction_params <- gillespie_test_params
gillespie_induction_params$num_init_cells <- 50L
gillespie_induction_params$sim_length <- list(0.5)
gillespie_induction_params$editing_induction <- list(
  timepoint = 0.2,
  num_cells = 10L
)
gillespie_induction <- gillespie_env$simulate_gillespie_population(
  gillespie_induction_params,
  seed = 91L,
  show_progress = FALSE
)
expect_equal(
  sum(gillespie_induction$event_log$event == 'editing_induction'),
  10L,
  'A fixed-size editing induction must select exactly the requested active cells.'
)
expect_true(
  all(
    gillespie_induction$nodes$birth_time[
      gillespie_induction$nodes$origin == 'induction_continuation'
    ] == 0.2
  ),
  'Partial induction must create exact-time continuation branches.'
)

# ---- Continuous-time versus fixed-timestep engine comparison ----

# Lightweight continuous-time versus timestep comparison.
comparison_env <- new.env(parent = gillespie_env)
sys.source(file.path(repo_root, 'analysis', 'engine_comparison.R'), envir = comparison_env)

#' The discrete comparison engine matches its pure-birth expectation
#'
#' Runs 10,000 cells for one cell cycle in four steps of 0.25 with no death, so
#' each step divides a cell with probability `1 - exp(-0.25)`. Asserts the final
#' population is within 2% of `10000 * (1 + (1 - exp(-0.25)))^4` and that the
#' per-cell-type counts sum to the reported final population.
comparison_params <- gillespie_test_params
comparison_params$num_init_cells <- 10000L
comparison_params$sim_length <- list(1)
comparison_params$cell_type_dict$cell_type_params$cycling$
  death_per_cell_cycle_prob <- 0
timestep_outcome <- comparison_env$simulate_timestep_population_outcome(
  comparison_params,
  end_time = 1,
  time_step = 0.25,
  seed = 212L,
  max_cells = 100000L
)
expected_timestep_population <- 10000 *
  (1 + (1 - exp(-0.25))) ^ 4
expect_true(
  abs(timestep_outcome$final_cells - expected_timestep_population) /
    expected_timestep_population < 0.02,
  'The discrete comparison engine must match its pure-birth expectation.'
)
expect_equal(
  sum(timestep_outcome$cell_type_counts),
  timestep_outcome$final_cells,
  'Time-step cell-type counts must sum to the final population.'
)

# ---- Non-Mendelian ecDNA propagation and recorder characters ----

# Non-Mendelian ecDNA propagation and recorder characters.
ecdna_env <- new.env(parent = physicell_env)
sys.source(file.path(repo_root, 'origin', 'R', 'ecdna_lineage.R'), envir = ecdna_env)

#' Fixture: a four-species ecDNA model and a unary-then-dividing lineage
#'
#' The lineage holds a founder segment, an induction continuation (a unary node
#' that is not a division), and two daughters of one true division. The model
#' starts with four species at ten copies each, half of them labeled, perfect
#' replication, and fair 50/50 segregation.
ecdna_params <- list(
  random_seed = 5L,
  cell_type_dict = list(
    founder_cell_type = 'ct1',
    cell_type_params = list(ct1 = list(cell_cycle_length = 1))
  ),
  ecdna_adapter = list(
    num_species = 4L,
    initial_copies_per_species = 10L,
    labeled_species_fraction = 0.5,
    static_id_length = 8L,
    num_recorder_targets = 3L,
    edit_probability_per_target_per_cell_cycle = 0.4,
    recorder_start_time = 0,
    replication_probability = 1,
    daughter_1_segregation_probability = 0.5,
    max_copies_per_cell = 1000L
  )
)
ecdna_nodes <- data.frame(
  node_id = c('ec_root', 'ec_continuation', 'ec_daughter_1', 'ec_daughter_2'),
  physicell_id = c('1', '1', '2', '3'),
  parent_node_id = c(NA, 'ec_root', 'ec_continuation', 'ec_continuation'),
  birth_time = c(0, 0.5, 1, 1),
  end_time = c(0.5, 1, 2, 2),
  branch_length = c(0.5, 0.5, 1, 1),
  division_event = c(NA, NA, 1, 1),
  origin = c(
    'founder', 'induction_continuation', 'daughter_1', 'daughter_2'
  ),
  is_terminal = c(FALSE, FALSE, TRUE, TRUE),
  cell_type = 'ct1',
  stringsAsFactors = FALSE
)
ecdna_model <- ecdna_env$prepare_ecdna_model(ecdna_params, seed = 17L)

#' Neutral ecDNA selection preserves the legacy seeded partition draws
#'
#' Recomputes the legacy sequence - one binomial replication draw per species
#' row followed by one binomial daughter-1 draw - under a fixed seed, then runs
#' `partition_ecdna_profile()` from the same seed. Asserts both daughters match
#' the hand-computed counts, so introducing selection did not perturb the RNG
#' stream of the neutral path.
ecdna_selection_profile <- data.frame(
  species_index = c(1L, 2L),
  labeled = c(FALSE, TRUE),
  haplotype = c('', '000'),
  count = c(10000L, 10000L),
  stringsAsFactors = FALSE
)
ecdna_neutral_model <- ecdna_model
ecdna_neutral_model$max_copies_per_cell <- 100000L
set.seed(31L)
ecdna_neutral_replicated <- ecdna_selection_profile$count + stats::rbinom(
  nrow(ecdna_selection_profile),
  size = ecdna_selection_profile$count,
  prob = ecdna_neutral_model$replication_probability
)
ecdna_neutral_daughter_1 <- stats::rbinom(
  nrow(ecdna_selection_profile),
  size = ecdna_neutral_replicated,
  prob = ecdna_neutral_model$daughter_1_segregation_probability
)
set.seed(31L)
ecdna_neutral_daughters <- ecdna_env$partition_ecdna_profile(
  ecdna_selection_profile,
  ecdna_neutral_model
)
expect_equal(
  ecdna_neutral_daughters[[1]]$count,
  ecdna_neutral_daughter_1,
  'Neutral ecDNA selection must preserve the legacy seeded daughter-1 draw.'
)
expect_equal(
  ecdna_neutral_daughters[[2]]$count,
  ecdna_neutral_replicated - ecdna_neutral_daughter_1,
  'Neutral ecDNA selection must preserve the legacy seeded daughter-2 draw.'
)

#' An ecDNA label selection coefficient of one blocks labeled replication
#'
#' Partitions a profile holding 10,000 unlabeled and 10,000 labeled copies with
#' the label coefficient at one and the edit coefficient at zero. Asserts the
#' unlabeled pool doubles to 20,000 while the labeled pool stays at 10,000
#' across both daughters.
ecdna_label_selection_model <- ecdna_model
ecdna_label_selection_model$label_selection_coefficient <- 1
ecdna_label_selection_model$recorder_edit_selection_coefficient <- 0
ecdna_label_selection_model$max_copies_per_cell <- 100000L
ecdna_selected_daughters <- ecdna_env$partition_ecdna_profile(
  ecdna_selection_profile,
  ecdna_label_selection_model
)
ecdna_joint_selected_profile <- do.call(rbind, ecdna_selected_daughters)
expect_equal(
  sum(ecdna_joint_selected_profile$count[
    !ecdna_joint_selected_profile$labeled
  ]),
  20000,
  'Unlabeled ecDNA copies must retain their neutral replication probability.'
)
expect_equal(
  sum(ecdna_joint_selected_profile$count[
    ecdna_joint_selected_profile$labeled
  ]),
  10000,
  'An ecDNA label selection coefficient of one must block labeled-copy replication.'
)

#' An ecDNA edit selection coefficient of one blocks edited-copy replication
#'
#' Partitions two labeled haplotypes of the same species - `000` unedited and
#' `100` edited - with the edit coefficient at one. Asserts the unedited
#' haplotype doubles to 20,000 while the edited haplotype stays at 10,000.
ecdna_edit_selection_profile <- data.frame(
  species_index = c(1L, 1L),
  labeled = c(TRUE, TRUE),
  haplotype = c('000', '100'),
  count = c(10000L, 10000L),
  stringsAsFactors = FALSE
)
ecdna_edit_selection_model <- ecdna_model
ecdna_edit_selection_model$label_selection_coefficient <- 0
ecdna_edit_selection_model$recorder_edit_selection_coefficient <- 1
ecdna_edit_selection_model$max_copies_per_cell <- 100000L
ecdna_edit_selected_daughters <- ecdna_env$partition_ecdna_profile(
  ecdna_edit_selection_profile,
  ecdna_edit_selection_model
)
ecdna_joint_edit_profile <- do.call(rbind, ecdna_edit_selected_daughters)
expect_equal(
  sum(ecdna_joint_edit_profile$count[
    ecdna_joint_edit_profile$haplotype == '000'
  ]),
  20000,
  'Unedited ecDNA recorder copies must retain neutral replication.'
)
expect_equal(
  sum(ecdna_joint_edit_profile$count[
    ecdna_joint_edit_profile$haplotype == '100'
  ]),
  10000,
  'An edit selection coefficient of one must block edited ecDNA-copy replication.'
)

#' Unary continuations do not replicate ecDNA; true divisions conserve each pool
#'
#' Replays ecDNA over the fixture lineage. Asserts the joint terminal copy
#' number is 80 - the 40 founder copies doubled by the single true division,
#' with the induction continuation contributing no replication - that each of
#' the four species conserves exactly 20 joint copies, that the labeled-species
#' fraction selected exactly two species, and that recorder edits carry
#' branch-event times inside the simulated interval.
ecdna_simulation <- ecdna_env$simulate_ecdna_on_lineage(
  ecdna_nodes,
  ecdna_model,
  terminal_physicell_ids = c('2', '3'),
  seed = 19L,
  show_progress = FALSE
)
ecdna_terminal_copy_numbers <- vapply(
  ecdna_simulation$profiles,
  function(profile) sum(profile$count),
  numeric(1)
)
expect_equal(
  sum(ecdna_terminal_copy_numbers),
  80,
  paste(
    'A unary induction continuation must not replicate ecDNA, while one true',
    'division with perfect replication must double the joint daughter pool.'
  )
)
for(species_index in seq_len(ecdna_model$num_species)){
  species_copy_total <- sum(vapply(
    ecdna_simulation$profiles,
    function(profile){
      sum(profile$count[profile$species_index == species_index])
    },
    numeric(1)
  ))
  expect_equal(
    species_copy_total,
    20,
    'Joint daughter segregation must conserve each perfectly replicated species pool.'
  )
}
expect_equal(
  length(ecdna_model$labeled_species_indices),
  2L,
  'The configured labeled-species fraction must select an exact species count.'
)
expect_true(
  nrow(ecdna_simulation$mutation_events) > 0L &&
    all(ecdna_simulation$mutation_events$first_event_time >= 0) &&
    all(ecdna_simulation$mutation_events$last_event_time <= 2),
  'ecDNA recorder edits must retain valid continuous branch-event times.'
)

#' The ecDNA character matrix has one column per labeled-species recorder target
#'
#' Asserts the recorder-edit presence matrix is two terminal cells by six
#' columns: two labeled species times three recorder targets.
ecdna_matrices <- ecdna_env$ecdna_feature_matrices(
  ecdna_simulation,
  ecdna_model
)
expect_equal(
  dim(ecdna_matrices$recorder_edit_presence),
  c(2L, 6L),
  'The ecDNA character matrix must contain one column per labeled static-ID target.'
)

# ---- Known-template prime editing on PhysiCell lineages ----

# Known prime-editing templates resolve once, retain exact event sequences,
# and respect pegRNA-specific zero/nonzero efficiencies.

#' Prime-editing preparation keeps known pegRNAs and unique integration IDs
#'
#' Configures three spaced nuclease targets across two integrations, wired to
#' pegRNAs with efficiencies 1, 0.5 and 0. Asserts the model is typed as prime
#' editing, that the pegRNA IDs stay in target order, and that the two
#' integrations receive distinct static IDs.
prime_recording_params <- physicell_params
prime_recording_params$max_bc_ints_per_cell <- list(2L)
prime_recording_params$nuclease_targets$num_targets <- 3L
prime_recording_params$nuclease_targets$config <- 'S:1:1'
prime_recording_params$nuclease_targets$prime_editing_system <- TRUE
prime_recording_params$physicell_adapter <- list(
  recorder_system = 'prime editing',
  profile_storage = 'sparse',
  compact_output = TRUE,
  retain_internal_profiles = FALSE
)
prime_recording_params$prime_editing_backend <- list(
  enabled = TRUE,
  target_pegRNA_ids = list('peg_one', 'peg_half', 'peg_zero'),
  induced_edit_probability_per_cell_cycle = 1,
  uninduced_edit_probability_per_cell_cycle = 0,
  static_id_length = 10L,
  pegRNAs = list(
    list(
      pegRNA_id = 'peg_one',
      edit_sequence = 'AGGCT',
      editing_efficiency = 1
    ),
    list(
      pegRNA_id = 'peg_half',
      edit_sequence = 'CTTGA',
      editing_efficiency = 0.5
    ),
    list(
      pegRNA_id = 'peg_zero',
      edit_sequence = 'GACCA',
      editing_efficiency = 0
    )
  )
)
prime_recording_model <- physicell_env$prepare_physicell_recording_model(
  prime_recording_params,
  num_integrations = 2L,
  seed = 79L
)
expect_true(
  isTRUE(prime_recording_model$is_prime_editing) &&
    identical(
      prime_recording_model$prime_editing$targets$pegRNA_id,
      c('peg_one', 'peg_half', 'peg_zero')
    ) &&
    length(unique(
      prime_recording_model$prime_editing$integration_static_ids
    )) == 2L,
  'Prime-editing preparation must retain known pegRNAs and unique integration IDs.'
)

#' Prime-editing events report their template and skip zero-efficiency pegRNAs
#'
#' Mutates the initial profile over a ten-unit branch starting at time 4.
#' Asserts the reported alternate sequences are exactly two copies each of the
#' `peg_one` and `peg_half` templates - one per integration - that every event
#' is a `prime_edit` timed inside [4, 14], and that the target driven by the
#' zero-efficiency pegRNA stays unedited in both integrations.
prime_initial_profile <- physicell_env$initialize_physicell_barcode_profile(
  prime_recording_model
)
set.seed(83L)
prime_mutation <- physicell_env$mutate_physicell_barcode_segment(
  prime_initial_profile,
  duration = 10,
  rate_set = prime_recording_model$rate_sets$induced_editing_params,
  model = prime_recording_model,
  segment_start = 4
)
expect_equal(
  sort(prime_mutation$events$alternate),
  sort(rep(c('AGGCT', 'CTTGA'), each = 2L)),
  'Prime-editing events must report each assigned edit template exactly.'
)
expect_true(
  all(prime_mutation$events$event == 'prime_edit') &&
    all(prime_mutation$events$event_time >= 4) &&
    all(prime_mutation$events$event_time <= 14) &&
    all(vapply(prime_mutation$profile, function(integration){
      identical(as.numeric(integration[c('1', '3', '5')]), c(1, 1, NA_real_))
    }, logical(1))),
  paste(
    'Prime-editing events must use continuous branch times and a zero-efficiency',
    'pegRNA must remain unedited.'
  )
)

#' A resolved prime-editing target stays locked against re-editing
#'
#' Runs a second branch segment over the already-edited profile and asserts no
#' further events are produced.
prime_locked <- physicell_env$mutate_physicell_barcode_segment(
  prime_mutation$profile,
  duration = 10,
  rate_set = prime_recording_model$rate_sets$induced_editing_params,
  model = prime_recording_model,
  segment_start = 14
)
expect_equal(
  nrow(prime_locked$events),
  0L,
  'A resolved prime-editing target must remain locked against re-editing.'
)

#' Prime-editing output writes state, character, and target-manifest files
#'
#' Replays prime editing over the fixture lineage and writes the outputs.
#' Asserts both sparse matrices exist, that the manifest lists all six targets
#' (three per integration) with static ID, pegRNA ID, edit sequence, efficiency
#' and the derived effective per-cell-cycle probability, and that the shared
#' barcode character output equals the written prime-editing character matrix.
#'
#' @section Side effects: Creates and then removes a temporary prime-editing
#'   output directory.
prime_recording_nodes <- physicell_env$build_physicell_lineage(
  physicell_divisions,
  end_time = 3
)
prime_recording_simulation <-
  physicell_env$simulate_recording_on_physicell_lineage(
    prime_recording_nodes,
    prime_recording_model,
    terminal_physicell_ids = physicell_live_ids,
    seed = 89L
  )
prime_output_dir <- tempfile('prime_editing_output_')
dir.create(prime_output_dir)
prime_recording_output <- physicell_env$write_physicell_recording_outputs(
  prime_recording_simulation,
  prime_recording_model,
  prime_output_dir,
  show_progress = FALSE
)
prime_layout <- physicell_env$read_physicell_csv(file.path(
  prime_output_dir,
  'prime_editing_target_manifest.csv'
))
expect_true(
  file.exists(file.path(
    prime_output_dir,
    'prime_editing_state_matrix_sparse.rds'
  )) &&
    file.exists(file.path(
      prime_output_dir,
      'prime_editing_character_matrix_sparse.rds'
    )) &&
    nrow(prime_layout) == 6L &&
    all(c(
      'static_id', 'pegRNA_id', 'edit_sequence', 'editing_efficiency',
      'induced_effective_probability_per_cell_cycle'
    ) %in% names(prime_layout)),
  paste(
    'Prime-editing output must include state/character matrices and the',
    'complete known-pegRNA target manifest.'
  )
)
expect_equal(
  prime_recording_output$binary_scores,
  readRDS(file.path(
    prime_output_dir,
    'prime_editing_character_matrix_sparse.rds'
  )),
  'The common barcode character output must preserve prime-editing states.'
)
unlink(prime_output_dir, recursive = TRUE)

# ---- PALINCODE cBit recorder ----

# PALINCODE cBits resolve once to left, right, or rare simultaneous-both states.

#' PALINCODE preparation yields a typed model with unique static integration IDs
#'
#' Configures three cBits per integration across two integrations, each cBit
#' forced to a single outcome (left, right, or both). Asserts the model is typed
#' as PALINCODE and that the two integrations carry distinct static IDs of the
#' configured length.
palincode_params <- physicell_params
palincode_params$max_bc_ints_per_cell <- list(2L)
palincode_params$physicell_adapter <- list(
  recorder_system = 'PALINCODE',
  profile_storage = 'sparse',
  compact_output = TRUE,
  retain_internal_profiles = FALSE
)
palincode_params$palincode_adapter <- list(
  num_cbits_per_integration = 3L,
  cbit_names = list('left_forced', 'right_forced', 'both_forced'),
  static_id_length = 12L,
  induced_edit_probability_per_cbit_per_cell_cycle = list(1, 1, 1),
  uninduced_edit_probability_per_cbit_per_cell_cycle = 0,
  left_edit_fraction = list(1, 0, 0),
  right_edit_fraction = list(0, 1, 0),
  both_edit_fraction = list(0, 0, 1)
)
palincode_model <- physicell_env$prepare_physicell_recording_model(
  palincode_params,
  num_integrations = 2L,
  seed = 71L
)
expect_true(
  isTRUE(palincode_model$is_palincode) &&
    identical(palincode_model$recorder_system, 'PALINCODE') &&
    length(unique(palincode_model$palincode$integration_static_ids)) == 2L &&
    all(nchar(palincode_model$palincode$integration_static_ids) == 12L),
  'PALINCODE preparation must create a typed model and unique static integration IDs.'
)

#' PALINCODE supports left, right, and simultaneous-both cBit outcomes
#'
#' Mutates the initial profile over a ten-unit branch starting at time 4 with
#' each cBit's conditional outcome fractions forced. Asserts exactly two of each
#' outcome are produced - one per integration - and that every event time falls
#' inside [4, 14].
palincode_initial_profile <- physicell_env$initialize_physicell_barcode_profile(
  palincode_model
)
palincode_mutation <- physicell_env$mutate_physicell_barcode_segment(
  palincode_initial_profile,
  duration = 10,
  rate_set = palincode_model$rate_sets$induced_editing_params,
  model = palincode_model,
  segment_start = 4
)
expect_equal(
  as.integer(table(factor(
    palincode_mutation$events$event,
    levels = c('palincode_left', 'palincode_right', 'palincode_both')
  ))),
  c(2L, 2L, 2L),
  'PALINCODE must support left, right, and simultaneous-both cBit outcomes.'
)
expect_true(
  all(palincode_mutation$events$event_time >= 4) &&
    all(palincode_mutation$events$event_time <= 14),
  'PALINCODE events must retain exact continuous times inside the branch segment.'
)

#' A resolved PALINCODE cBit stays locked and keeps its recorded state
#'
#' Runs a second branch segment over the already-edited profile. Asserts no
#' further events are produced and that the profile is unchanged.
palincode_locked <- physicell_env$mutate_physicell_barcode_segment(
  palincode_mutation$profile,
  duration = 10,
  rate_set = palincode_model$rate_sets$induced_editing_params,
  model = palincode_model,
  segment_start = 14
)
expect_equal(
  nrow(palincode_locked$events),
  0L,
  'A resolved PALINCODE cBit must remain locked against sequential re-editing.'
)
expect_equal(
  palincode_locked$profile,
  palincode_mutation$profile,
  'PALINCODE state locking must preserve every resolved cBit state.'
)

#' PALINCODE character output allocates left/right/both columns per cBit
#'
#' Lifts one edited profile into the sparse recording matrix and then into the
#' categorical character matrix. Asserts the matrix is one cell by 18 columns
#' (two integrations times three cBits times three outcomes) and that the six
#' edited cBits occupy exactly one column each.
palincode_raw_matrix <- physicell_env$physicell_sparse_recording_matrix(
  list(cell_test = palincode_mutation$profile),
  palincode_model
)
palincode_character_matrix <-
  physicell_env$physicell_palincode_character_matrix(
    palincode_raw_matrix,
    palincode_model
  )
expect_equal(
  dim(palincode_character_matrix),
  c(1L, 18L),
  'PALINCODE character output must allocate left/right/both columns per cBit.'
)
expect_equal(
  sum(palincode_character_matrix),
  6,
  'Each edited PALINCODE cBit must occupy exactly one categorical character column.'
)

#' PALINCODE conditional outcome fractions must sum to one
#'
#' Overrides `both_edit_fraction` so each cBit's left/right/both fractions no
#' longer sum to one, and asserts model preparation refuses the configuration.
expect_error(
  physicell_env$prepare_physicell_recording_model(
    within(palincode_params, {
      palincode_adapter$both_edit_fraction <- list(0.1, 0.1, 0.1)
    }),
    num_integrations = 2L,
    seed = 71L
  ),
  'PALINCODE conditional outcome fractions that do not sum to one must be rejected.'
)

#' PALINCODE output writes state, character, and calibrated cBit layout files
#'
#' Replays PALINCODE over the fixture lineage and writes the outputs. Asserts
#' both sparse matrices exist, that the target layout lists all six cBits with
#' their static ID, name and calibrated left/right/both fractions, and that the
#' shared barcode character output equals the written PALINCODE character
#' matrix.
#'
#' @section Side effects: Creates and then removes a temporary PALINCODE output
#'   directory.
palincode_nodes <- physicell_env$build_physicell_lineage(
  physicell_divisions,
  end_time = 3
)
palincode_simulation <- physicell_env$simulate_recording_on_physicell_lineage(
  palincode_nodes,
  palincode_model,
  terminal_physicell_ids = physicell_live_ids,
  seed = 73L
)
palincode_output_dir <- tempfile('palincode_output_')
dir.create(palincode_output_dir)
palincode_output <- physicell_env$write_physicell_recording_outputs(
  palincode_simulation,
  palincode_model,
  palincode_output_dir,
  show_progress = FALSE
)
palincode_layout <- physicell_env$read_physicell_csv(file.path(
  palincode_output_dir,
  'barcode_target_layout.csv'
))
expect_true(
  file.exists(file.path(
    palincode_output_dir,
    'palincode_state_matrix_sparse.rds'
  )) &&
    file.exists(file.path(
      palincode_output_dir,
      'palincode_character_matrix_sparse.rds'
    )) &&
    nrow(palincode_layout) == 6L &&
    all(c(
      'static_id', 'target_name', 'left_edit_fraction',
      'right_edit_fraction', 'both_edit_fraction'
    ) %in% names(palincode_layout)),
  'PALINCODE output must include state/character matrices and a calibrated cBit layout.'
)
expect_equal(
  palincode_output$binary_scores,
  readRDS(file.path(
    palincode_output_dir,
    'palincode_character_matrix_sparse.rds'
  )),
  'The shared barcode character output must preserve PALINCODE outcome orientation.'
)
unlink(palincode_output_dir, recursive = TRUE)

# ---- Visium array geometry, spatial capture, and lineage queries ----

# Conventional Visium geometry, spatial capture, and binary-lifting lineage
# queries must preserve the exact simulated genealogy without a quadratic
# cophenetic matrix.

#' The conventional 6.5-mm Visium array is 4,992 staggered spots in 78 rows
#'
#' Builds the array at a 100-unit pitch with no rotation or translation. Asserts
#' 78 rows of 64 spots each and the even/odd staggering, where even array rows
#' hold only even array columns and odd rows only odd columns.
visium_array <- physicell_env$make_visium_6_5mm_array(
  spot_pitch = 100,
  rotation_degrees = 0,
  translation = c(0, 0),
  slice_id = 'test'
)
expect_true(
  nrow(visium_array) == 4992L &&
    length(unique(visium_array$array_row)) == 78L &&
    all(table(visium_array$array_row) == 64L) &&
    all(visium_array$array_col[visium_array$array_row == 0L] %% 2L == 0L) &&
    all(visium_array$array_col[visium_array$array_row == 1L] %% 2L == 1L),
  'The conventional 6.5-mm Visium array must use 4,992 staggered spots in 78 rows.'
)

#' Every horizontal or diagonal Visium neighbor is exactly one pitch apart
#'
#' Measures the in-plane distance from the first spot to its same-row neighbor
#' and to its diagonal neighbor in the next row, asserting both equal the
#' configured 100-unit pitch.
same_row_distance <- sqrt(sum(
  (as.numeric(visium_array[1, c('plane_u', 'plane_v')]) -
     as.numeric(visium_array[2, c('plane_u', 'plane_v')]))^2
))
diagonal_index <- which(
  visium_array$array_row == 1L & visium_array$array_col == 1L
)[1]
diagonal_distance <- sqrt(sum(
  (as.numeric(visium_array[1, c('plane_u', 'plane_v')]) -
     as.numeric(visium_array[diagonal_index, c('plane_u', 'plane_v')]))^2
))
expect_equal(
  c(same_row_distance, diagonal_distance),
  c(100, 100),
  'Every horizontal or diagonal neighboring Visium spot must be one pitch apart.'
)

#' Visium lineage queries return exact MRCAs, distances, and founder separation
#'
#' Builds a binary-lifting index over a six-node fixture containing two separate
#' founders, then queries three pairs. Asserts the MRCAs are `n3`, `n1` and `NA`
#' for the cross-founder pair; that patristic distances are the elapsed-time
#' sums 2, 4 and 6; and that division distances count lineage edges within a
#' founder while staying `NA` across founders rather than inventing a count. The
#' index answers pairs directly, avoiding a quadratic cophenetic matrix.
visium_test_nodes <- data.frame(
  node_id = paste0('n', 1:6),
  parent_node_id = c(NA, 'n1', 'n1', 'n3', 'n3', NA),
  birth_time = c(0, 1, 1, 2, 2, 0),
  end_time = c(1, 3, 2, 3, 3, 3),
  branch_length = c(1, 2, 1, 1, 1, 3),
  stringsAsFactors = FALSE
)
visium_lca <- physicell_env$build_visium_lca_index(visium_test_nodes)
visium_queries <- physicell_env$query_visium_lineage_pairs(
  visium_lca,
  c(4L, 2L, 2L),
  c(5L, 4L, 6L)
)
expect_equal(
  visium_queries$mrca_node_id,
  c('n3', 'n1', NA_character_),
  'Visium lineage queries must return exact MRCAs and identify separate founders.'
)
expect_equal(
  visium_queries$patristic_distance,
  c(2, 4, 6),
  'Visium lineage queries must match elapsed-time patristic distance.'
)
expect_equal(
  visium_queries$division_distance,
  c(2L, 3L, NA_integer_),
  paste(
    'Visium lineage queries must report within-founder lineage-edge separation',
    'without inventing a cross-founder division count.'
  )
)

#' Section capture respects finite thickness and circular spot footprints
#'
#' Slices four cells against a z-normal plane of thickness 5 centered at the
#' origin, so only the three cells within 2.5 units of the plane are sampled.
#' The retained cells are then assigned to two spots 200 units apart with a
#' 55-unit spot diameter. Asserts all three section cells fall inside a
#' footprint and that spot cellularity is updated to 2 and 1.
visium_basis <- physicell_env$visium_plane_basis(c(0, 0, 1))
visium_test_cells <- data.frame(
  sample_id = paste0('cell_', 1:4),
  physicell_id = as.character(1:4),
  node_id = c('n2', 'n4', 'n5', 'n6'),
  x = c(0, 10, 100, 200),
  y = c(0, 0, 0, 0),
  z = c(0, 1, 3, 0),
  stringsAsFactors = FALSE
)
visium_slice <- physicell_env$slice_physicell_cells(
  visium_test_cells,
  visium_basis,
  center = c(0, 0, 0),
  offset = 0,
  thickness = 5
)
expect_equal(
  visium_slice$sample_id,
  c('cell_1', 'cell_2', 'cell_4'),
  'Only cell centers inside the configured finite section thickness may be sampled.'
)
small_spots <- data.frame(
  barcode = c('a', 'b'), spot_id = c('s1', 's2'), in_tissue = 0L,
  array_row = 0:1, array_col = 0:1, plane_u = c(0, 200),
  plane_v = c(0, 0), n_cells = 0L, stringsAsFactors = FALSE
)
assigned_slice <- physicell_env$assign_cells_to_visium_spots(
  visium_slice,
  small_spots,
  spot_diameter = 55
)
expect_true(
  identical(assigned_slice$cells$captured, c(TRUE, TRUE, TRUE)) &&
    identical(assigned_slice$spots$n_cells, c(2L, 1L)),
  'Section cells must map to circular Visium footprints and update spot cellularity.'
)

#' Empty or sparsely sampled sections keep typed zero-count outputs
#'
#' Assigns an empty section to the same spots and summarizes an empty spot-pair
#' frame. Asserts no cells are produced, that every spot reports zero
#' cellularity, and that the pair-lineage summary still returns its full
#' seven-row typed skeleton.
empty_assigned_slice <- physicell_env$assign_cells_to_visium_spots(
  visium_slice[FALSE, , drop = FALSE],
  small_spots,
  spot_diameter = 55
)
empty_spot_summary <- physicell_env$summarize_visium_spot_pair_lineage(
  data.frame(),
  spot_pitch = 100
)
expect_true(
  nrow(empty_assigned_slice$cells) == 0L &&
    sum(empty_assigned_slice$spots$n_cells) == 0L &&
    nrow(empty_spot_summary$summary) == 7L,
  'Empty or sparsely sampled sections must retain typed zero-count Visium outputs.'
)

# ---- Lineage benchmark grid ----

# The simulation-only benchmark must derive paired cell/integration subsets
# with one logical target column per physical target or cBit.

#' The benchmark grid covers every condition with one column per logical target
#'
#' Runs the simulation-only benchmark over four recorder systems, two tree
#' sizes, two integration counts and two mitochondrial observation depths.
#' Asserts the manifest holds all 16 conditions, that each non-mitochondrial
#' system reports its expected targets per integration (50 BASELINE, six
#' prime-editing, two PALINCODE), and that every logical target matrix has one
#' column per logical target.
#'
#' @section Side effects: Creates a temporary benchmark output directory that
#'   the resume test below reuses; it is removed afterwards.
benchmark_output_dir <- tempfile('lineage_benchmark_output_')
benchmark_result <- physicell_env$run_lineage_benchmark(
  output_dir = benchmark_output_dir,
  shapes = 'balanced',
  tree_sizes = c(4L, 8L),
  integration_counts = c(1L, 2L),
  mt_observation_depths = c(1L, 2L),
  mt_genomes_per_cell = 4L,
  systems = c('baseline', 'prime', 'palincode', 'mitochondrial'),
  seeds = 1L,
  show_progress = FALSE
)
benchmark_manifest <- benchmark_result$manifest
expect_equal(
  nrow(benchmark_manifest),
  16L,
  'The benchmark grid must contain every system/tree-size/integration condition.'
)
expected_benchmark_targets <- c(
  baseline = 50L,
  prime = 6L,
  palincode = 2L
)
expect_true(
  all(
    benchmark_manifest$targets_per_integration[
      benchmark_manifest$system != 'mitochondrial'
    ] == expected_benchmark_targets[
      benchmark_manifest$system[
        benchmark_manifest$system != 'mitochondrial'
      ]
    ]
  ) &&
    all(
      benchmark_manifest$logical_target_matrix_columns ==
        benchmark_manifest$logical_targets
    ),
  paste(
    'Benchmark logical matrices must preserve 50 BASELINE, six prime-editing,',
    'or two PALINCODE targets per integration.'
  )
)

#' Mitochondrial conditions separate sampled-genome depth from integration count
#'
#' Asserts the mitochondrial rows carry no integration count, vary only along
#' observation depth, declare `sampled_mitochondrial_genomes` as their
#' observation dimension, and record the configured biological genome count per
#' cell.
mitochondrial_conditions <- benchmark_manifest[
  benchmark_manifest$system == 'mitochondrial',
  ,
  drop = FALSE
]
expect_true(
  all(is.na(mitochondrial_conditions$integrations)) &&
    identical(
      sort(unique(mitochondrial_conditions$observation_depth)),
      c(1L, 2L)
    ) &&
    all(
      mitochondrial_conditions$observation_dimension ==
        'sampled_mitochondrial_genomes'
    ) &&
    all(mitochondrial_conditions$biological_mitochondrial_genomes == 4L),
  paste(
    'Mitochondrial conditions must distinguish sampled-genome depth from',
    'integrated-recorder counts.'
  )
)

#' Every benchmark logical matrix matches its requested cells and targets
#'
#' Loads each condition's stored logical target matrix and asserts its
#' dimensions equal the manifest's tree size by logical target count.
for(condition_index in seq_len(nrow(benchmark_manifest))){
  logical_matrix <- readRDS(file.path(
    benchmark_output_dir,
    benchmark_manifest$logical_target_matrix[condition_index]
  ))
  expect_equal(
    dim(logical_matrix),
    c(
      benchmark_manifest$tree_size[condition_index],
      benchmark_manifest$logical_targets[condition_index]
    ),
    'Every benchmark logical matrix must match its requested cells and targets.'
  )
}

#' Nested mitochondrial sampling keeps a common variant universe
#'
#' Compares the eight-cell mitochondrial character matrices observed at depth 1
#' and depth 2. Asserts they share dimensions, that no variant seen at the
#' shallower depth is lost at the deeper one, and that the exact recording-event
#' files exist for every mitochondrial condition.
mitochondrial_n8 <- mitochondrial_conditions[
  mitochondrial_conditions$tree_size == 8L,
  ,
  drop = FALSE
]
mitochondrial_depth1 <- readRDS(file.path(
  benchmark_output_dir,
  mitochondrial_n8$character_matrix[
    mitochondrial_n8$observation_depth == 1L
  ]
))
mitochondrial_depth2 <- readRDS(file.path(
  benchmark_output_dir,
  mitochondrial_n8$character_matrix[
    mitochondrial_n8$observation_depth == 2L
  ]
))
expect_true(
  identical(dim(mitochondrial_depth1), dim(mitochondrial_depth2)) &&
    all(as.matrix(mitochondrial_depth1) <= as.matrix(mitochondrial_depth2)) &&
    all(file.exists(file.path(
      benchmark_output_dir,
      mitochondrial_conditions$full_recording_events
    ))),
  paste(
    'Nested mitochondrial sampling must retain a common variant universe,',
    'preserve variants seen at shallower depth, and expose exact events.'
  )
)

#' A resumed benchmark reuses completed simulations, manifest unchanged
#'
#' Re-runs the identical benchmark against the same output directory and asserts
#' the manifest is unchanged, i.e. resume detects finished conditions instead of
#' re-simulating them under a fresh RNG stream.
resumed_benchmark <- physicell_env$run_lineage_benchmark(
  output_dir = benchmark_output_dir,
  shapes = 'balanced',
  tree_sizes = c(4L, 8L),
  integration_counts = c(1L, 2L),
  mt_observation_depths = c(1L, 2L),
  mt_genomes_per_cell = 4L,
  systems = c('baseline', 'prime', 'palincode', 'mitochondrial'),
  seeds = 1L,
  show_progress = FALSE
)
expect_equal(
  resumed_benchmark$manifest,
  benchmark_manifest,
  'A resumed benchmark must reuse completed simulations without changing its manifest.'
)
unlink(benchmark_output_dir, recursive = TRUE)

# ---- Completion ----

#' Reaching this line means every assertion above held
cat('All regression tests passed.\n')
