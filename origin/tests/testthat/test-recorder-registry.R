# Recorder registration.
#
# These pin the property the registry exists for: a recorder is either
# dispatched or absent, and the check for that reflects the wiring rather than
# the parameters. The previous arrangement installed recorders by assigning over
# the simulator's functions, so configuring the parameters set the model flag a
# guard would test while leaving the dispatch untouched -- a run then produced
# plausible output from the wrong model, and the guard passed.

teardown_registry <- function() {
  for (name in origin_recorders()) unregister_origin_recorder(name)
}

test_that("a registered recorder is reported, an unregistered one is not", {
  on.exit(teardown_registry(), add = TRUE)
  teardown_registry()
  expect_null(origin_active_recorder(list(is_demo = TRUE)))

  register_origin_recorder(
    "demo", applies = function(model) isTRUE(model$is_demo),
    mutate = function(profile, duration, rate_set, model, segment_start = 0) {
      list(profile = "demo", events = NULL)
    })
  expect_equal(origin_recorders(), "demo")
  expect_equal(origin_active_recorder(list(is_demo = TRUE)), "demo")
  expect_null(origin_active_recorder(list(is_demo = FALSE)))

  expect_true(unregister_origin_recorder("demo"))
  expect_null(origin_active_recorder(list(is_demo = TRUE)))
  expect_false(unregister_origin_recorder("demo"))
})

test_that("the flag alone cannot make a recorder look dispatched", {
  # The failure this replaces: a model carrying the recorder's own flag, set
  # from the parameters, with nothing registered. The old guard asked the model
  # and passed; this one asks the registry and does not.
  on.exit(teardown_registry(), add = TRUE)
  teardown_registry()
  configured_but_unwired <- list(is_demo = TRUE)
  expect_null(origin_active_recorder(configured_but_unwired))
  expect_null(origin_recorder_hook(configured_but_unwired, "mutate"))
})

test_that("mutation dispatches to the registered recorder", {
  on.exit(teardown_registry(), add = TRUE)
  teardown_registry()
  register_origin_recorder(
    "demo", applies = function(model) isTRUE(model$is_demo),
    mutate = function(profile, duration, rate_set, model, segment_start = 0) {
      list(profile = "dispatched", events = NULL)
    })
  model <- list(is_demo = TRUE, num_integrations = 1L)
  result <- mutate_physicell_barcode_segment(
    profile = matrix(0, 1, 4), duration = 1, rate_set = list(), model = model)
  expect_equal(result$profile, "dispatched")

  # A model the recorder does not claim takes the built-in path untouched.
  unregister_origin_recorder("demo")
  expect_null(origin_recorder_hook(model, "mutate"))
})

test_that("layout and character-matrix hooks are fetched independently", {
  on.exit(teardown_registry(), add = TRUE)
  teardown_registry()
  register_origin_recorder(
    "demo", applies = function(model) isTRUE(model$is_demo),
    layout = function(model) data.frame(position = 1L),
    character_matrix = function(raw_alleles, binary_scores, model, integrations) {
      "characters"
    })
  model <- list(is_demo = TRUE)
  expect_equal(physicell_baseline_target_layout(model)$position, 1L)
  hook <- origin_recorder_hook(model, "character_matrix")
  expect_equal(hook(NULL, NULL, model, 1L), "characters")
  # A hook the recorder leaves unset reads as absent, not as an error.
  expect_null(origin_recorder_hook(model, "mutate"))
})

test_that("prepare hooks all run and must return a model", {
  on.exit(teardown_registry(), add = TRUE)
  teardown_registry()
  register_origin_recorder("one", applies = function(model) FALSE,
                           prepare = function(model, params) {
                             model$seen <- c(model$seen, "one"); model
                           })
  register_origin_recorder("two", applies = function(model) FALSE,
                           prepare = function(model, params) {
                             model$seen <- c(model$seen, "two"); model
                           })
  out <- origin_apply_recorder_prepare(list(), list())
  expect_setequal(out$seen, c("one", "two"))

  register_origin_recorder("bad", applies = function(model) FALSE,
                           prepare = function(model, params) "not a model")
  expect_error(origin_apply_recorder_prepare(list(), list()),
               "did not return a model")
})

test_that("registration validates its arguments", {
  on.exit(teardown_registry(), add = TRUE)
  teardown_registry()
  expect_error(register_origin_recorder("", applies = function(m) TRUE,
                                        mutate = identity),
               "non-empty string")
  expect_error(register_origin_recorder("x", applies = "not a function",
                                        mutate = identity),
               "applies")
  expect_error(register_origin_recorder("x", applies = function(m) TRUE,
                                        mutate = "not a function"),
               "functions or NULL")
  expect_error(register_origin_recorder("x", applies = function(m) TRUE),
               "at least one hook")
})

test_that("a recorder whose predicate errors is skipped, not fatal", {
  # A half-built model reaching applies() should not take down the run.
  on.exit(teardown_registry(), add = TRUE)
  teardown_registry()
  register_origin_recorder("explodes",
                           applies = function(model) stop("boom"),
                           mutate = identity)
  expect_null(origin_active_recorder(list()))
})
