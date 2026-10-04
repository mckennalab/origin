# Empirical per-target editing rates.
#
# The Gamma dispersion model cannot reproduce a measured BASELINE construct:
# fitted at its best it still read 0.236 mean site rate against 0.305 observed,
# over-producing dead targets and under-producing saturated ones, because a
# one-parameter family has no way to be both. Supplying the measured rates
# directly removes the family from the question. These tests check that the
# rates are honoured, that the Gamma path is untouched when they are absent,
# and that the class stratification really is dropped.

test_that("an absent empirical block leaves the gamma path alone", {
  expect_null(physicell_target_empirical_rates(NULL))
  expect_null(physicell_target_empirical_rates(list()))
  expect_null(physicell_target_empirical_rates(list(num_targets = 10)))

  classes <- setNames(rep("High", 50), as.character(seq_len(50)))
  set.seed(1)
  gamma_rates <- draw_physicell_target_rates(classes, 0.02, shape = 0.5)
  set.seed(1)
  unchanged <- draw_physicell_target_rates(classes, 0.02, shape = 0.5,
                                           empirical = NULL)
  expect_equal(gamma_rates, unchanged)
})

test_that("supplied rates are the ones that come back", {
  # Three distinct values: every drawn rate must be one of them, and over many
  # targets all three must appear, which a mean-and-shape model cannot promise.
  pool <- c(0.001, 0.05, 0.4)
  spec <- list(edit_rate_empirical = pool)
  expect_equal(physicell_target_empirical_rates(spec), pool)
  # JSON list boxing survives, as it does for the dispersion shape.
  expect_equal(
    physicell_target_empirical_rates(list(edit_rate_empirical = as.list(pool))),
    pool
  )

  classes <- setNames(rep("High", 500), as.character(seq_len(500)))
  set.seed(2)
  rates <- draw_physicell_target_rates(classes, 0.02, shape = 0.5,
                                       empirical = pool)
  expect_length(rates, 500)
  expect_true(all(rates %in% pool))
  expect_setequal(unique(rates), pool)
  expect_equal(names(rates), names(classes))
})

test_that("empirical draws ignore the class stratification", {
  # The gamma path sends Low and High targets to different quantile bins, so
  # their rates differ systematically. Empirical draws must not: the whole point
  # is independent sampling, because stratifying flattens the spread between
  # integrations that the measured construct shows.
  pool <- stats::runif(2000, 0, 0.3)
  classes <- setNames(rep(c("Low", "High"), each = 1000),
                      as.character(seq_len(2000)))
  set.seed(3)
  rates <- draw_physicell_target_rates(classes, 0.02, shape = 0.5,
                                       empirical = pool)
  low <- mean(rates[classes == "Low"])
  high <- mean(rates[classes == "High"])
  expect_lt(abs(low - high), 0.02)

  # Two integrations drawn independently differ; under stratification they
  # would be forced to the same composition.
  draw_one <- function() {
    mean(draw_physicell_target_rates(
      setNames(rep("High", 272), as.character(seq_len(272))),
      0.02, shape = 0.5, empirical = pool))
  }
  set.seed(4)
  means <- replicate(30, draw_one())
  expect_gt(stats::sd(means), 0)
})

test_that("invalid empirical rates are rejected, not coerced", {
  for (bad in list(numeric(0), c(0.1, NA), c(0.1, Inf), c(-0.1, 0.2),
                   c(0.1, 1.5), "0.1", c("0.1", "0.2"))) {
    expect_error(
      physicell_target_empirical_rates(list(edit_rate_empirical = bad)),
      "edit_rate_empirical"
    )
  }
  classes <- setNames("High", "1")
  expect_error(
    draw_physicell_target_rates(classes, 0.02, empirical = numeric(0)),
    "finite and non-negative"
  )
  expect_error(
    draw_physicell_target_rates(classes, 0.02, empirical = c(0.1, -1)),
    "finite and non-negative"
  )
})

# ---- optional outputs --------------------------------------------------------
#
# The output phase cost 83% of a large run and every byte of it went to files
# the analyses never opened. These checks pin the defaults, because a default
# that silently flips back would not fail anything else here -- the runs would
# just get slow again.

test_that("the expensive outputs are off unless asked for", {
  arguments <- formals(write_physicell_recording_outputs)
  for (flag in c("write_allele_matrix", "write_mutation_events",
                 "write_profiles", "write_barcode_fasta")) {
    expect_false(isTRUE(eval(arguments[[flag]])), label = flag)
  }
  # The outputs that are read downstream keep their old behaviour.
  expect_true(isTRUE(eval(arguments$write_lineage)))
  expect_true(isTRUE(eval(arguments$compress_csv)))
})
