source("load_origin.R")

# The previous implementation, kept here only to prove the ranged one agrees.
old_apply <- function(profile, integration, position, left_bases, right_bases,
                      barcode_length) {
  left_bases <- max(0L, as.integer(left_bases)); right_bases <- max(0L, as.integer(right_bases))
  emptied <- integer()
  consume <- function(profile, coordinate, budget) {
    available <- physicell_deletable_bases(physicell_barcode_value(profile, integration, coordinate))
    if (available == 0) return(list(profile = profile, budget = budget, emptied = FALSE))
    if (budget >= available) {
      profile <- set_physicell_barcode_value(profile, integration, coordinate, -1)
      return(list(profile = profile, budget = budget - available, emptied = TRUE))
    }
    profile <- set_physicell_barcode_value(profile, integration, coordinate, 0)
    list(profile = profile, budget = 0L, emptied = FALSE)
  }
  step <- consume(profile, position, .Machine$integer.max); profile <- step$profile
  if (step$emptied) emptied <- c(emptied, position)
  budget <- left_bases; coordinate <- position - 1L
  while (budget > 0 && coordinate >= 1L) {
    step <- consume(profile, coordinate, budget); profile <- step$profile
    budget <- step$budget; if (step$emptied) emptied <- c(emptied, coordinate)
    coordinate <- coordinate - 1L
  }
  budget <- right_bases; coordinate <- position + 1L
  while (budget > 0 && coordinate <= barcode_length) {
    step <- consume(profile, coordinate, budget); profile <- step$profile
    budget <- step$budget; if (step$emptied) emptied <- c(emptied, coordinate)
    coordinate <- coordinate + 1L
  }
  list(profile = profile, positions = sort(emptied))
}

make_sparse <- function(n_int, len, seed) {
  set.seed(seed)
  structure(lapply(seq_len(n_int), function(i) {
    k <- sample(seq_len(len), 12)                 # a few prior edits
    v <- sample(c(-1, 0.1, 0.2, 0.3, 0.4), 12, TRUE)
    stats::setNames(v, as.character(k))
  }), class = c("physicell_sparse_barcode", "list"))
}

LEN <- 215L
cat("=== equivalence over randomised cuts ===\n")
set.seed(7); mismatch <- 0
for (trial in 1:400) {
  p <- make_sparse(2, LEN, trial)
  pos <- sample(LEN, 1); lb <- sample(0:90, 1); rb <- sample(0:90, 1)
  a <- old_apply(p, 1L, pos, lb, rb, LEN)
  b <- apply_physicell_deletion(p, 1L, pos, lb, rb, LEN)
  ra <- physicell_barcode_row(a$profile, 1L, LEN)
  rb2 <- physicell_barcode_row(b$profile, 1L, LEN)
  if (!identical(a$positions, b$positions) || !isTRUE(all.equal(ra, rb2))) mismatch <- mismatch + 1
}
cat(sprintf("  %d of 400 trials disagree\n", mismatch))

cat("\n=== timing, 2000 resected cuts (mean 65bp/side) ===\n")
bench <- function(fn, label) {
  set.seed(3); p <- make_sparse(2, LEN, 1)
  t0 <- Sys.time()
  for (i in 1:2000) {
    pos <- sample(LEN, 1)
    r <- fn(p, 1L, pos, 65L, 65L, LEN)
    if (i %% 50 == 0) p <- make_sparse(2, LEN, i)   # reset so it stays realistic
  }
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  %-10s %6.2fs\n", label, el)); el
}
o <- bench(old_apply, "per-base"); n <- bench(apply_physicell_deletion, "ranged")
cat(sprintf("\nspeedup: %.1fx\n", o / n))
