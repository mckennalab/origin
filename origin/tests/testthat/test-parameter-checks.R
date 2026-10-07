# Guards for three configuration mistakes that used to fail unhelpfully or not
# at all. Each test asserts the message actually names the offending parameter,
# because a check that fires with a generic message is barely better than the
# original failure.
#
# Note when editing these: build the invalid cases by direct assignment, not
# modifyList(). modifyList merges nested lists recursively, so replacing
# bc_nuc_composition with a partial list leaves the original fields in place and
# the "invalid" case silently becomes valid.

editing_state <- function() {
  list(
    be_mutations_per_target_per_division = 0.02,
    nuc_insertions_per_target_per_division = 0,
    nuc_deletions_per_target_per_division = 0,
    bc_bg_insertion_prob_per_division = 0,
    bc_bg_deletion_prob_per_division = 0,
    bc_substitution_model = "JC", bc_sub_model_params = list(0),
    mt_substitution_model = "JC", mt_sub_model_params = list(0),
    mt_bg_insertion_prob_per_division = 0,
    mt_bg_deletion_prob_per_division = 0
  )
}

recording_params <- function() {
  state <- editing_state()
  list(
    num_init_cells = 1, sim_length = list(3), random_seed = 1,
    bc_length = 300, max_bc_ints_per_cell = list(1),
    bc_nuc_composition = list(frac_a = 0.25, frac_g = 0.25,
                              frac_c = 0.25, frac_t = 0.25),
    be_conversion_pattern = "A --> G",
    be_targets = list(
      num_targets = 8, config = "S:10:24",
      edit_rate_class_fractions = list(high = 0.5, medium = 0, low = 0.5),
      editing_window = list(size = 0, decaying = FALSE, close_after_edit = FALSE)
    ),
    editing_induction = list(timepoint = 0, num_cells = NULL, frac_cells = 1),
    differentiation_induction = list(timepoint = 0, num_cells = NULL,
                                     frac_cells = 1),
    cell_type_dict = list(
      founder_cell_type = "p",
      cell_type_params = list(p = c(list(
        cell_cycle_length = 1, death_per_cell_cycle_prob = 0,
        bc_invariant_sites = 0, mt_invariant_sites = 0,
        induced_editing_params = state, uninduced_editing_params = state
      ))),
      uninduced_transition_matrix = list(list(1)),
      induced_transition_matrix = list(list(1))
    )
  )
}

run_pipeline <- function(params) {
  run_gillespie_lineage_pipeline(
    params, params_path = NULL,
    output_dir = file.path(tempdir(), paste0("checks-", sample.int(1e6, 1))),
    overrides = list(progress = FALSE, seed = 1, modalities = list("barcode"))
  )
}

test_that("a missing bc_nuc_composition is reported by name", {
  params <- recording_params()
  params$bc_nuc_composition <- NULL
  expect_error(run_pipeline(params), "bc_nuc_composition is required")
})

test_that("a partly specified bc_nuc_composition names the absent fractions", {
  params <- recording_params()
  params$bc_nuc_composition <- list(frac_a = 0.25, frac_g = 0.25)
  expect_error(run_pipeline(params), "frac_c, frac_t")
})

test_that("bc_nuc_composition fractions that cannot form a sequence are rejected", {
  params <- recording_params()
  params$bc_nuc_composition <- list(frac_a = 0, frac_g = 0, frac_c = 0,
                                    frac_t = 0)
  expect_error(run_pipeline(params), "sum above zero")
})

test_that("a target layout longer than the barcode reports the arithmetic", {
  params <- recording_params()
  params$be_targets$num_targets <- 40
  # The message should carry the last target position, the configured length,
  # and a target count that would fit, so the caller does not have to rederive
  # which of the four numbers to change.
  expect_error(run_pipeline(params), "last target at 985")
  expect_error(run_pipeline(params), "bc_length is 300")
  expect_error(run_pipeline(params), "at most 12")
})

test_that("inducing fewer cells than founders warns instead of diluting silently", {
  params <- recording_params()
  params$num_init_cells <- 5
  params$editing_induction <- list(timepoint = 0, num_cells = 1,
                                   frac_cells = NULL)
  expect_warning(run_pipeline(params), "4 of the 5 founders")
})

test_that("valid configurations raise nothing", {
  # Not expect_silent(): the pipeline always prints a timing summary and its
  # OUTPUT_DIR line to stdout. What matters is that it raises no condition.
  expect_warning(
    {
      result <- run_pipeline(recording_params())
      unlink(result$output_dir, recursive = TRUE)
    },
    regexp = NA
  )

  # Several founders with frac_cells = 1 induces everything, so no warning.
  params <- recording_params()
  params$num_init_cells <- 5
  expect_warning(
    {
      result <- run_pipeline(params)
      unlink(result$output_dir, recursive = TRUE)
    },
    regexp = NA
  )
})

test_that("inducing a subset after the founders have divided is not warned about", {
  # Inducing a fraction of a grown population is an ordinary design, so the
  # warning must be limited to inductions at time zero.
  params <- recording_params()
  params$num_init_cells <- 5
  params$editing_induction <- list(timepoint = 1, num_cells = 1,
                                   frac_cells = NULL)
  expect_warning(
    {
      result <- run_pipeline(params)
      unlink(result$output_dir, recursive = TRUE)
    },
    regexp = NA
  )
})
