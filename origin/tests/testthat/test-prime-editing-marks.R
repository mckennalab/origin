# PEtracer-style mark alphabets for prime editing.
#
# PEtracer installs one of eight predefined 5nt marks at each edit site
# (doi:10.1126/science.adx3800). Modelling that site as binary edited/unedited
# throws away the information that distinguishes two independent edits from two
# cells sharing an ancestor, so these tests pin that an edit records WHICH mark
# it installed, and that the single-mark path is unchanged.

pool <- data.frame(
  pool_index = 1:8,
  pegRNA_id = sprintf("mark_%d", 1:8),
  edit_sequence = c("ACGTA", "AGTCG", "ATCGT", "CAGTA",
                    "CTAGC", "GATCA", "GCTAG", "TGCAT"),
  editing_efficiency = rep(1, 8),
  spacer_sequence = rep("GACTGACTGACTGACTGACT", 8),
  pbs_sequence = rep("CGTACGTAC", 8),
  rtt_sequence = c("ACGTA", "AGTCG", "ATCGT", "CAGTA",
                   "CTAGC", "GATCA", "GCTAG", "TGCAT"),
  description = rep("test mark", 8),
  stringsAsFactors = FALSE
)

test_that("one mark per target reproduces the single-outcome recorder", {
  marks <- expand_prime_editing_marks(pool, c(5L, 11L, 17L), marks_per_target = 1L)
  expect_equal(nrow(marks), 3L)
  expect_true(all(marks$allele == 1))
  expect_true(all(marks$mark_probability == 1))
  # One pegRNA per site, taken in pool order, as the flat assignment did.
  expect_equal(marks$pegRNA_id, c("mark_1", "mark_2", "mark_3"))
})

test_that("each target gets its own block of distinct marks", {
  marks <- expand_prime_editing_marks(pool, c(5L, 11L), marks_per_target = 4L)
  expect_equal(nrow(marks), 8L)
  by_target <- split(marks, marks$target_position)
  for (block in by_target) {
    expect_equal(nrow(block), 4L)
    # A site must not be able to install the same mark twice: duplicates would
    # silently halve its real alphabet.
    expect_equal(anyDuplicated(block$edit_sequence), 0L)
    expect_equal(sort(block$allele), 1:4)
    expect_equal(sum(block$mark_probability), 1)
  }
  # Different sites draw different blocks, so a mark is not shared across sites
  # by construction.
  expect_false(identical(by_target[[1]]$edit_sequence,
                         by_target[[2]]$edit_sequence))
})

test_that("an alphabet larger than the pool is refused, not silently cycled", {
  expect_error(
    expand_prime_editing_marks(pool, 5L, marks_per_target = 9L),
    "marks_per_target"
  )
  for (bad in list(0L, -1L, NA_integer_, c(2L, 3L))) {
    expect_error(
      expand_prime_editing_marks(pool, 5L, marks_per_target = bad),
      "marks_per_target"
    )
  }
})

test_that("mark probability follows relative efficiency, not absolute", {
  skewed <- pool[1:2, ]
  skewed$editing_efficiency <- c(0.75, 0.25)
  marks <- expand_prime_editing_marks(skewed, 5L, marks_per_target = 2L)
  expect_equal(marks$mark_probability, c(0.75, 0.25))

  # Scaling every efficiency changes how often the SITE edits, which is the
  # site's own rate, and must not change which mark wins.
  scaled <- skewed
  scaled$editing_efficiency <- scaled$editing_efficiency / 10
  expect_equal(
    expand_prime_editing_marks(scaled, 5L, marks_per_target = 2L)$mark_probability,
    marks$mark_probability
  )
})

test_that("an all-zero efficiency block falls back to a uniform alphabet", {
  # Otherwise the weights would be 0/0 and sampling would fail at run time,
  # long after the pool was configured.
  dead <- pool[1:4, ]
  dead$editing_efficiency <- rep(0, 4)
  marks <- expand_prime_editing_marks(dead, 5L, marks_per_target = 4L)
  expect_equal(marks$mark_probability, rep(0.25, 4))
})

test_that("the backend reports an encoding that names every mark", {
  params <- list(
    nuclease_targets = list(prime_editing_system = TRUE),
    prime_editing_backend = list(
      enabled = TRUE, marks_per_target = 8,
      pegRNAs = pool[, setdiff(names(pool), "pool_index")]
    )
  )
  backend <- prepare_prime_editing_backend(params, c(5L, 11L, 17L), seed = 1)
  skip_if(is.null(backend), "prime editing backend unavailable")
  expect_equal(backend$marks_per_target, 8L)
  expect_equal(nrow(backend$marks), 24L)
  expect_equal(unname(backend$state_encoding), 0:8)
  expect_equal(names(backend$state_encoding),
               c("unedited", sprintf("mark_%d", 1:8)))
  # The per-site table stays one row per site so existing consumers still work.
  expect_equal(nrow(backend$targets), 3L)
})

test_that("each site's rate comes from its own mark block, not the first", {
  # assign_prime_editing_targets() walks the pool one entry per site, so with a
  # blocked mark pool every site would otherwise inherit block one's
  # efficiency. That silently gives all sites the same rate, which is invisible
  # in the output until per-site saturation is inspected.
  blocked <- do.call(rbind, lapply(1:3, function(site) {
    data.frame(
      pegRNA_id = sprintf("s%d_m%d", site, 1:4),
      edit_sequence = sprintf("%s%s", c("AA", "AC", "AG", "AT"),
                              c("A", "C", "G", "T")[site]),
      editing_efficiency = c(0.9, 0.5, 0.2)[site],
      spacer_sequence = "GACTGACTGACTGACTGACT",
      pbs_sequence = "CGTACGTAC",
      rtt_sequence = "AAAAA",
      description = "blocked mark",
      stringsAsFactors = FALSE
    )
  }))
  params <- list(
    nuclease_targets = list(prime_editing_system = TRUE),
    prime_editing_backend = list(
      enabled = TRUE, marks_per_target = 4, pegRNAs = blocked
    )
  )
  backend <- prepare_prime_editing_backend(params, c(4L, 15L, 26L), seed = 1)
  skip_if(is.null(backend), "prime editing backend unavailable")
  expect_equal(backend$targets$editing_efficiency, c(0.9, 0.5, 0.2))
  # And the marks each site can install come from its own block.
  by_site <- split(backend$marks$editing_efficiency,
                   backend$marks$target_position)
  expect_equal(unname(vapply(by_site, unique, numeric(1))), c(0.9, 0.5, 0.2))
})
