# Building a parameter list.
#
# The simulator takes a deeply nested list, and callers have been assembling it
# by hand: roughly thirty lines of nesting for a minimal run, with every
# editing rate repeated across the induced and uninduced blocks of every cell
# type. Mistakes surface as validation errors raised deep in the call stack,
# where the message names a field but not which of several copies is wrong.
#
# origin_params() fills the structure from flat arguments and validates on the
# way in. It is a convenience over the list, not a replacement for it: the
# result is an ordinary parameter list and can be modified afterwards like any
# other.

#' Build a parameter list for a Gillespie run
#'
#' @export
#' @param sim_length Simulated time units.
#' @param founders Starting cells.
#' @param cell_cycle_length Mean time between divisions.
#' @param death Death probability per cell cycle. Growth is steeply non-linear
#'   in this near the critical point; see `origin_growth_response()`.
#' @param barcode_length Positions per integration.
#' @param integrations Barcode integrations per cell.
#' @param be_rate Base-editing probability per target per division.
#' @param nuc_insertion_rate,nuc_deletion_rate Nuclease rates per target per
#'   division.
#' @param targets Targets per integration, or `NULL` for none.
#' @param target_config Target layout string, e.g. `"S:15:20"`.
#' @param conversion Base-editing conversion pattern.
#' @param composition Barcode base composition, as a named list of fractions.
#'   The default is all-A, which makes every position an eligible A to G target;
#'   a uniform composition leaves only about a quarter of them editable.
#' @param seed Random seed.
#' @param mito_genomes Mitochondrial genome copies per cell.
#' @param mito_rate Mitochondrial substitution probability per genome per
#'   division.
#' @return A parameter list suitable for `run_gillespie_lineage_pipeline()`.
origin_params <- function(sim_length = 10,
                          founders = 10L,
                          cell_cycle_length = 1,
                          death = 0.4,
                          barcode_length = 100L,
                          integrations = 5L,
                          be_rate = 0,
                          nuc_insertion_rate = 0,
                          nuc_deletion_rate = 0,
                          targets = NULL,
                          target_config = 'S:1:0',
                          conversion = 'A --> G',
                          composition = list(frac_a = 1, frac_g = 0,
                                             frac_c = 0, frac_t = 0),
                          seed = 1L,
                          mito_genomes = 0L,
                          mito_rate = 0){
  positive <- function(value, name, allow_zero = FALSE){
    value <- suppressWarnings(as.numeric(value))
    if(length(value) != 1L || !is.finite(value) ||
       (allow_zero && value < 0) || (!allow_zero && value <= 0)){
      stop(sprintf('%s must be one finite %s number.', name,
                   if(allow_zero) 'non-negative' else 'positive'))
    }
    value
  }
  probability <- function(value, name){
    value <- suppressWarnings(as.numeric(value))
    if(length(value) != 1L || !is.finite(value) || value < 0 || value > 1){
      stop(sprintf('%s must be one probability in [0, 1].', name))
    }
    value
  }
  sim_length <- positive(sim_length, 'sim_length')
  cell_cycle_length <- positive(cell_cycle_length, 'cell_cycle_length')
  death <- probability(death, 'death')
  if(death >= 1){
    stop('death must be below 1, or no lineage survives.')
  }
  be_rate <- probability(be_rate, 'be_rate')
  nuc_insertion_rate <- probability(nuc_insertion_rate, 'nuc_insertion_rate')
  nuc_deletion_rate <- probability(nuc_deletion_rate, 'nuc_deletion_rate')
  mito_rate <- probability(mito_rate, 'mito_rate')
  founders <- as.integer(positive(founders, 'founders'))
  barcode_length <- as.integer(positive(barcode_length, 'barcode_length'))
  integrations <- as.integer(positive(integrations, 'integrations'))
  mito_genomes <- as.integer(positive(mito_genomes, 'mito_genomes',
                                      allow_zero = TRUE))
  if(!is.null(targets)){
    targets <- as.integer(positive(targets, 'targets'))
    if(targets > barcode_length){
      stop('targets cannot exceed barcode_length.')
    }
  }
  if(be_rate == 0 && nuc_insertion_rate == 0 && nuc_deletion_rate == 0 &&
     mito_rate == 0){
    stop('Every editing rate is zero, so the run would record nothing.')
  }

  # One editing block, shared by the induced and uninduced states of every cell
  # type. Hand-built lists repeat these fourteen fields per state per cell type,
  # which is where they drift apart.
  editing <- list(
    be_mutations_per_target_per_division = be_rate,
    nuc_insertions_per_target_per_division = nuc_insertion_rate,
    nuc_deletions_per_target_per_division = nuc_deletion_rate,
    mt_mutations_per_genome_per_division = mito_rate,
    num_mt_genomes = mito_genomes,
    bc_bg_insertion_prob_per_division = 0,
    bc_bg_deletion_prob_per_division = 0,
    bc_substitution_model = 'JC',
    bc_sub_model_params = list(0),
    mt_substitution_model = 'JC',
    mt_sub_model_params = list(0),
    mt_bg_insertion_prob_per_division = 0,
    mt_bg_deletion_prob_per_division = 0
  )
  params <- list(
    num_init_cells = founders,
    sim_length = list(sim_length),
    random_seed = seed,
    bc_length = barcode_length,
    mito_genome_length = 16600,
    max_bc_ints_per_cell = list(integrations),
    bc_nuc_composition = composition,
    be_conversion_pattern = conversion,
    editing_induction = list(timepoint = 0, num_cells = NULL, frac_cells = 1),
    differentiation_induction = list(timepoint = 0, num_cells = NULL,
                                     frac_cells = 1),
    cell_type_dict = list(
      founder_cell_type = 'progenitor',
      cell_type_params = list(progenitor = list(
        cell_cycle_length = cell_cycle_length,
        death_per_cell_cycle_prob = death,
        bc_invariant_sites = 0,
        mt_invariant_sites = 0,
        induced_editing_params = editing,
        uninduced_editing_params = editing
      )),
      uninduced_transition_matrix = list(list(1)),
      induced_transition_matrix = list(list(1))
    )
  )
  if(!is.null(targets)){
    params$be_targets <- list(
      num_targets = targets,
      config = target_config,
      edit_rate_class_fractions = list(high = 1, medium = 0, low = 0),
      editing_window = list(size = 0, decaying = FALSE,
                            close_after_edit = FALSE)
    )
  }
  params
}
