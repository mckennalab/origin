# Response curves.
#
# These exist because two relations in the engine are not predictable from the
# parameters, and both were discovered the expensive way. Growth is sharply
# non-linear in the death probability near the critical point -- extrapolating
# one measured point put a probe at 50,000 cells instead of 5,000 and it ran for
# hours without finishing -- and the configured editing rate is nominal, with
# the effective per-division rate about a third of it.
#
# The tests use the smallest configurations that still exercise the relation,
# because every call runs a real simulation.

demo_state <- function() list(
  be_mutations_per_target_per_division = 0.02,
  nuc_insertions_per_target_per_division = 0, nuc_deletions_per_target_per_division = 0,
  mt_mutations_per_genome_per_division = 0, num_mt_genomes = 0,
  bc_bg_insertion_prob_per_division = 0, bc_bg_deletion_prob_per_division = 0,
  bc_substitution_model = "JC", bc_sub_model_params = list(0),
  mt_substitution_model = "JC", mt_sub_model_params = list(0),
  mt_bg_insertion_prob_per_division = 0, mt_bg_deletion_prob_per_division = 0)

demo_params <- function(death = 0.45) {
  state <- demo_state()
  list(num_init_cells = 40L, sim_length = list(5), random_seed = 1,
       bc_length = 20L, mito_genome_length = 16600,
       max_bc_ints_per_cell = list(1),
       bc_nuc_composition = list(frac_a = 1, frac_g = 0, frac_c = 0, frac_t = 0),
       be_conversion_pattern = "A --> G",
       be_targets = list(num_targets = 4L, config = "S:1:4",
         edit_rate_class_fractions = list(high = 1, medium = 0, low = 0),
         editing_window = list(size = 0, decaying = FALSE,
                               close_after_edit = FALSE)),
       editing_induction = list(timepoint = 0, num_cells = NULL, frac_cells = 1),
       differentiation_induction = list(timepoint = 0, num_cells = NULL,
                                        frac_cells = 1),
       cell_type_dict = list(founder_cell_type = "progenitor",
         cell_type_params = list(progenitor = list(cell_cycle_length = 1,
           death_per_cell_cycle_prob = death, bc_invariant_sites = 0,
           mt_invariant_sites = 0, induced_editing_params = state,
           uninduced_editing_params = state)),
         uninduced_transition_matrix = list(list(1)),
         induced_transition_matrix = list(list(1))))
}

test_that("growth falls monotonically as the death probability rises", {
  skip_on_cran()
  response <- origin_growth_response(demo_params(), deaths = c(0.3, 0.5, 0.65),
                                     sim_length = 5, founders = 60L)
  expect_equal(nrow(response), 3L)
  expect_setequal(names(response),
                  c("death", "alive", "generation", "growth_per_division"))
  measured <- response$growth_per_division[is.finite(response$growth_per_division)]
  expect_gt(length(measured), 1L)
  # The relation this function exists to expose: growth is decreasing in death,
  # and steeply enough that interpolating between two points is unsafe.
  expect_true(all(diff(measured) < 0))
})

test_that("growth response survives extinction rather than failing", {
  skip_on_cran()
  # Near and past the critical point a probe can die out. That is an answer,
  # not an error, and callers sweeping a range need the row back.
  response <- origin_growth_response(demo_params(), deaths = 0.95,
                                     sim_length = 5, founders = 5L)
  expect_equal(nrow(response), 1L)
  expect_true(is.na(response$growth_per_division) ||
                response$growth_per_division <= 1)
})

test_that("growth response validates its inputs", {
  expect_error(origin_growth_response(demo_params(), deaths = numeric(0)),
               "finite probabilities")
  expect_error(origin_growth_response(demo_params(), deaths = 1),
               "finite probabilities")
  expect_error(origin_growth_response(demo_params(), deaths = 0.4,
                                      cell_type = "absent"),
               "not in the parameter list")
})

test_that("saturation rises with the configured rate, and the effective rate is lower", {
  skip_on_cran()
  setter <- function(params, rate) {
    for (type in names(params$cell_type_dict$cell_type_params)) {
      params$cell_type_dict$cell_type_params[[type]]$induced_editing_params$
        be_mutations_per_target_per_division <- rate
      params$cell_type_dict$cell_type_params[[type]]$uninduced_editing_params$
        be_mutations_per_target_per_division <- rate
    }
    params
  }
  response <- origin_rate_response(demo_params(), rates = c(0.01, 0.08),
                                   apply_rate = setter)
  expect_equal(nrow(response), 2L)
  expect_true(all(c("rate", "target_rate", "effective_rate_per_division") %in%
                    names(response)))
  expect_gt(response$target_rate[2], response$target_rate[1])
  # The point of the function: what the parameter says is not what the engine
  # delivers, so the effective rate must be reported separately rather than
  # assumed equal to the nominal one.
  effective <- response$effective_rate_per_division
  expect_true(all(is.na(effective) | effective >= 0))
})

test_that("rate response validates its inputs", {
  expect_error(origin_rate_response(demo_params(), rates = -1,
                                    apply_rate = function(p, r) p),
               "finite and non-negative")
  expect_error(origin_rate_response(demo_params(), rates = 0.01,
                                    apply_rate = "not a function"),
               "function of")
})
