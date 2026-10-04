# Score each local-class share on the span distribution, which is what the
# parameter actually controls.
#
# Span 0 is dropped and both sides renormalised: 11.4% of real molecules are
# pristine because those cells suppressed Cas9, which is a population the model
# has no way to express and which we are deliberately not fitting. Leaving it in
# would charge every candidate the same fixed penalty and bias the comparison
# toward whichever one edits least.
coll <- read.delim("analysis/results/flare_real/deletion_span_collapsed.tsv")
real <- coll$molecules[coll$largest_shared_span >= 1]
real <- real / sum(real)
targets <- 15L + 26L * (0:7)

rows <- list()
prefix <- Sys.getenv("PREFIX", "flare_local")
for (p in c("0", "0.15", "0.335", "0.5")) {
  f <- sprintf("analysis/results/%s_%s/calibrated_molecules.rds", prefix, p)
  if (!file.exists(f)) { cat(sprintf("local %-6s : not finished\n", p)); next }
  m <- readRDS(f)$molecules
  span <- apply(m, 1, function(r) { v <- r[r != 0]; if (!length(v)) 0L else max(table(v)) })
  sim <- table(factor(span[span >= 1], levels = 1:8))
  sim <- as.numeric(sim) / sum(sim)
  rows[[p]] <- data.frame(
    local_share = as.numeric(p),
    rate = read.csv(sprintf("analysis/results/%s_%s/calibrated_settings.csv", prefix, p))$rate_multiplier,
    target_rate = mean(m != 0),
    span_ge7 = sum(sim[7:8]),
    span_2to6 = sum(sim[2:6]),
    distinct_calls = mean(apply(m, 1, function(r) length(unique(r)))),
    L1_distance = sum(abs(sim - real)),
    frac_8of8 = mean(rowSums(m != 0) == 8),
    frac_7or8 = mean(rowSums(m != 0) >= 7),
    frac_1to6 = mean(rowSums(m != 0) >= 1 & rowSums(m != 0) <= 6))
}
out <- do.call(rbind, rows)
if (!is.null(out)) {
  cat(sprintf("\nreal (span>=1, renormalised): span>=7 %.3f, span 2-6 %.3f\n\n",
              sum(real[7:8]), sum(real[2:6])))
  print(out, row.names = FALSE, digits = 4)
  cat(sprintf("\nbest by L1 on the span distribution: local share %s\n",
              out$local_share[which.min(out$L1_distance)]))
}
