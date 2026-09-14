# Per-target editing-rate dispersion.
#
# The gamma shape controls how unequal target rates are (CV = 1/sqrt(shape)).
# It was previously fixed at 0.5, which could not reproduce a measured BASELINE
# construct: at 0.5 the simulator produced no saturating targets against 17%
# observed. These tests pin the default for backwards compatibility, check that
# the parameter is honoured end to end, and check that it is validated rather
# than silently coerced.

test_that("the dispersion default stays 0.5 for existing parameter files", {
  expect_equal(physicell_target_dispersion(NULL), 0.5)
  expect_equal(physicell_target_dispersion(list()), 0.5)
  expect_equal(
    physicell_target_dispersion(list(num_targets = 10)), 0.5
  )
  # A block that sets one gets it back, including through JSON's list boxing.
  expect_equal(
    physicell_target_dispersion(list(edit_rate_dispersion_shape = 0.09)),
    0.09
  )
  expect_equal(
    physicell_target_dispersion(list(edit_rate_dispersion_shape = list(0.09))),
    0.09
  )
})

test_that("an invalid dispersion shape is rejected, not coerced", {
  for (bad in list(0, -1, NA_real_, Inf, c(0.1, 0.2), "0.1")) {
    expect_error(
      physicell_target_dispersion(list(edit_rate_dispersion_shape = bad)),
      "edit_rate_dispersion_shape"
    )
  }
  classes <- setNames(rep("High", 20), as.character(seq_len(20)))
  expect_error(
    draw_physicell_target_rates(classes, 0.01, shape = 0),
    "dispersion shape"
  )
  expect_error(
    draw_physicell_target_rates(classes, 0.01, shape = NA_real_),
    "dispersion shape"
  )
})

test_that("a smaller shape spreads target rates further apart", {
  classes <- setNames(rep(c("High", "Low"), each = 250),
                      as.character(seq_len(500)))
  mean_rate <- 0.02
  set.seed(1)
  wide <- draw_physicell_target_rates(classes, mean_rate, shape = 0.05)
  set.seed(1)
  narrow <- draw_physicell_target_rates(classes, mean_rate, shape = 5)

  # Dispersion is the point, so compare spread rather than level.
  expect_gt(stats::sd(wide) / mean(wide), stats::sd(narrow) / mean(narrow))
  # A heavy-tailed draw puts more targets near zero AND more far above the mean;
  # a single-tailed check would pass on a distribution that merely shifted.
  expect_gt(mean(wide < mean_rate / 10), mean(narrow < mean_rate / 10))
  expect_gt(max(wide), max(narrow))
})

test_that("the default call path is unchanged by the new argument", {
  classes <- setNames(rep(c("High", "Low"), each = 30),
                      as.character(seq_len(60)))
  set.seed(42)
  implicit <- draw_physicell_target_rates(classes, 0.01)
  set.seed(42)
  explicit <- draw_physicell_target_rates(classes, 0.01, shape = 0.5)
  expect_identical(implicit, explicit)
})

test_that("zero and empty inputs still short-circuit before the shape is used", {
  classes <- setNames(rep("High", 5), as.character(seq_len(5)))
  # Zero mean returns zeros whatever the shape, and must not divide by it.
  expect_equal(
    unname(draw_physicell_target_rates(classes, 0, shape = 0.09)),
    rep(0, 5)
  )
  expect_length(
    draw_physicell_target_rates(setNames(character(), character()), 0.01), 0
  )
})

test_that("the shape reaches the rate set from the target block", {
  # physicell_rate_set() reads be_targets and nuclease_targets separately, so a
  # shape set on one block must not leak into the other.
  spec <- list(
    be_targets = list(edit_rate_dispersion_shape = 0.09),
    nuclease_targets = list()
  )
  expect_equal(physicell_target_dispersion(spec$be_targets), 0.09)
  expect_equal(physicell_target_dispersion(spec$nuclease_targets), 0.5)
})
