# Deletions as intervals.
#
# A Cas9 deletion removes a contiguous run of sequence. The recorder used to
# flip the cut position to -1 and leave every base between targets untouched,
# so a deletion had no length and no endpoints, and two events that removed
# different sequence were indistinguishable. The GSM8791703 HL60 recording
# carries 4,579 distinct deletion alleles over 223 start coordinates and 190
# lengths, none of which that model can express. These tests pin the ported
# behaviour: bases are consumed outwards from the cut, on both sides.

profile_of <- function(values) matrix(values, nrow = 1L)

test_that("a position reports the bases it can still lose", {
  expect_equal(physicell_deletable_bases(0), 1L)      # reference only
  expect_equal(physicell_deletable_bases(3), 1L)      # substitution
  expect_equal(physicell_deletable_bases(0.2), 2L)    # reference + insert
  expect_equal(physicell_deletable_bases(-1), 0L)     # already gone
  for (bad in list(NA_real_, Inf, numeric(0), c(1, 2))) {
    expect_error(physicell_deletable_bases(bad), "one finite encoded allele")
  }
})

test_that("a deletion consumes a contiguous run on both sides of the cut", {
  profile <- profile_of(rep(0, 21))
  result <- apply_physicell_deletion(profile, 1L, 11L, left_bases = 3L,
                                     right_bases = 2L, barcode_length = 21L)
  expect_equal(result$positions, 8:13)
  expect_true(all(result$profile[1L, 8:13] == -1))
  # Everything outside the interval is untouched, which is the whole point:
  # the old behaviour left the spacer bases alone and only marked targets.
  expect_true(all(result$profile[1L, c(1:7, 14:21)] == 0))
})

test_that("the cut position goes even when both extents are zero", {
  profile <- profile_of(rep(0, 5))
  result <- apply_physicell_deletion(profile, 1L, 3L, 0L, 0L, 5L)
  expect_equal(result$positions, 3L)
  expect_equal(result$profile[1L, ], c(0, 0, -1, 0, 0))
})

test_that("deletions stop at the ends of the barcode", {
  profile <- profile_of(rep(0, 5))
  result <- apply_physicell_deletion(profile, 1L, 1L, left_bases = 10L,
                                     right_bases = 10L, barcode_length = 5L)
  expect_equal(result$positions, 1:5)
  expect_true(all(result$profile[1L, ] == -1))
})

test_that("an insertion costs two bases, and a partial hit keeps the reference", {
  # Budget of 1 against a position holding reference + insert: the inserted
  # base goes, the reference survives, and the run stops there.
  profile <- profile_of(c(0, 0.3, 0, 0, 0))
  result <- apply_physicell_deletion(profile, 1L, 3L, left_bases = 1L,
                                     right_bases = 0L, barcode_length = 5L)
  expect_equal(result$profile[1L, 2], 0)
  expect_false(2L %in% result$positions)
  # Budget of 2 empties it.
  profile <- profile_of(c(0, 0.3, 0, 0, 0))
  result <- apply_physicell_deletion(profile, 1L, 3L, left_bases = 2L,
                                     right_bases = 0L, barcode_length = 5L)
  expect_equal(result$profile[1L, 2], -1)
  expect_true(2L %in% result$positions)
})

test_that("already-deleted positions are free to cross", {
  # A later deletion reaching across an earlier one must not be stopped by it,
  # otherwise spans could never accumulate over divisions.
  profile <- profile_of(c(0, 0, -1, 0, 0, 0, 0))
  result <- apply_physicell_deletion(profile, 1L, 6L, left_bases = 3L,
                                     right_bases = 0L, barcode_length = 7L)
  expect_true(all(result$profile[1L, 2:6] == -1))
  expect_false(3L %in% result$positions)   # it was already gone
})

test_that("both extents are drawn, and drawn independently", {
  set.seed(1)
  draws <- t(replicate(400, draw_physicell_deletion_extent()))
  expect_setequal(colnames(draws), c("left", "right"))
  expect_true(all(draws >= 1))
  # Legacy defaults are mean-one exponential, so the median lands at 1 and the
  # two sides must not be locked together.
  expect_lt(abs(stats::median(draws[, "left"]) - 1), 1.5)
  expect_false(all(draws[, "left"] == draws[, "right"]))

  # A construct that needs long deletions asks for them.
  set.seed(2)
  long <- t(replicate(400, draw_physicell_deletion_extent(
    list(left = list(shape = 4, rate = 0.05),
         right = list(shape = 4, rate = 0.05)))))
  expect_gt(mean(long[, "left"]), 40)
  expect_gt(mean(long[, "right"]), 40)
})

test_that("an invalid extent block is rejected, not silently defaulted", {
  for (bad in list(list(left = list(shape = 0, rate = 1)),
                   list(left = list(shape = 1, rate = -1)),
                   list(right = list(shape = NA_real_, rate = 1)),
                   list(right = list(shape = "1", rate = 1)))) {
    expect_error(draw_physicell_deletion_extent(bad), "shape and rate")
  }
})

test_that("a local component makes the draw a two-class mixture", {
  # Measured lengths are bimodal: a third of cuts repair locally and never
  # reach the next target, the rest resect. One gamma fitted across both sits
  # between the modes and describes neither.
  spec <- list(left = list(shape = 4.1475, rate = 0.06407),
               right = list(shape = 4.1475, rate = 0.06407),
               local = list(probability = 0.335, shape = 0.5931,
                            rate = 0.15473))
  set.seed(11)
  totals <- replicate(4000, sum(draw_physicell_deletion_extent(spec)) + 1L)
  # Both modes must be present, and in roughly the measured proportions.
  expect_gt(mean(totals <= 26), 0.20)
  expect_lt(mean(totals <= 26), 0.50)
  expect_gt(mean(totals > 26), 0.50)
  # A probability of zero must reproduce the pure long class.
  spec$local$probability <- 0
  set.seed(11)
  long_only <- replicate(2000, sum(draw_physicell_deletion_extent(spec)) + 1L)
  expect_lt(mean(long_only <= 26), 0.05)
})

test_that("an invalid local probability is rejected", {
  for (bad in list(-0.1, 1.5, NA_real_, "0.3", c(0.2, 0.3))) {
    expect_error(
      draw_physicell_deletion_extent(list(
        left = list(shape = 1, rate = 1), right = list(shape = 1, rate = 1),
        local = list(probability = bad, shape = 1, rate = 1))),
      "local probability")
  }
})
