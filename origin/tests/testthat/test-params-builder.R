# Building a parameter list.
#
# The structure is deep enough that callers were writing thirty lines of nesting
# per run, repeating the fourteen-field editing block across the induced and
# uninduced states of every cell type. The two things worth pinning are that the
# result actually runs, and that the mistakes which used to surface deep in the
# simulator are caught at the boundary instead.

test_that("the built list runs through the pipeline", {
  skip_on_cran()
  params <- origin_params(sim_length = 4, founders = 5L, barcode_length = 20L,
                          integrations = 2L, be_rate = 0.05, targets = 4L)
  output <- file.path(tempdir(), "origin-params-builder")
  dir.create(output, showWarnings = FALSE)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  expect_error(
    run_gillespie_lineage_pipeline(params, params_path = NULL,
                                   output_dir = output,
                                   overrides = list(progress = FALSE, seed = 1,
                                                    modalities = "barcode")),
    NA)
  expect_true(file.exists(file.path(output,
                                    "barcode_binary_score_matrix.csv.gz")))
})

test_that("the editing block is shared by both states, so they cannot drift", {
  params <- origin_params(be_rate = 0.02, targets = 4L)
  cell <- params$cell_type_dict$cell_type_params$progenitor
  expect_identical(cell$induced_editing_params, cell$uninduced_editing_params)
  expect_equal(cell$induced_editing_params$be_mutations_per_target_per_division,
               0.02)
})

test_that("scalars are boxed where the reader expects a length-one list", {
  # sim_length and max_bc_ints_per_cell are read as boxed scalars; a bare number
  # there fails validation inside the simulator rather than here.
  params <- origin_params(sim_length = 7, integrations = 3L, be_rate = 0.01)
  expect_true(is.list(params$sim_length))
  expect_equal(params$sim_length[[1L]], 7)
  expect_true(is.list(params$max_bc_ints_per_cell))
  expect_equal(params$max_bc_ints_per_cell[[1L]], 3L)
})

test_that("a run that would record nothing is rejected", {
  # Every rate zero is a configuration error, not a simulation worth running.
  expect_error(origin_params(be_rate = 0, nuc_insertion_rate = 0,
                             nuc_deletion_rate = 0, mito_rate = 0),
               "record nothing")
})

test_that("invalid values are caught at the boundary", {
  expect_error(origin_params(death = 1, be_rate = 0.01), "below 1")
  expect_error(origin_params(death = -0.1, be_rate = 0.01), "probability")
  expect_error(origin_params(be_rate = 1.5), "probability")
  expect_error(origin_params(sim_length = 0, be_rate = 0.01), "positive")
  expect_error(origin_params(founders = 0, be_rate = 0.01), "positive")
  expect_error(origin_params(be_rate = 0.01, barcode_length = 10L,
                             targets = 20L),
               "cannot exceed")
})

test_that("targets are optional and omitted cleanly", {
  expect_null(origin_params(be_rate = 0.01)$be_targets)
  expect_equal(origin_params(be_rate = 0.01, targets = 6L,
                             target_config = "S:15:20")$be_targets$config,
               "S:15:20")
})

test_that("the default composition makes every position editable", {
  # All-A with an A to G conversion; a uniform composition would leave only
  # about a quarter of positions eligible and silently quarter the recorder.
  params <- origin_params(be_rate = 0.01)
  expect_equal(params$bc_nuc_composition$frac_a, 1)
  expect_equal(params$be_conversion_pattern, "A --> G")
})
