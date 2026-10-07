# Methods — Lineage recorder editing rates

## Overview

Editing rates for the two recorders with published or in-house measurements —
BASELINE base editing and prime editing — were converted from cumulative
observed editing into per-cell-division rates, so that they are directly
comparable with the per-division parameters the simulator takes. BASELINE rates
were estimated from a clone 84 recording; prime editing rates were taken from
Choi et al., Supplementary Table 1.

## Converting cumulative editing to a per-division rate

Both recorders write irreversible edits, so a target is unedited at sampling
only if it escaped editing at every division. For a constant per-division
editing probability *r* over *n* divisions, the cumulative edited fraction is

    p = 1 - (1 - r)^n,   hence   r = 1 - (1 - p)^(1/n)

All per-division rates below use this inversion, and each was checked by
round-trip. The linear approximation *r ≈ p/n* is adequate only for small *p*:
at p = 0.059 it is 2% low, at p = 0.1456 it is 7% low, and at the cumulative
rates reached by saturated BASELINE targets it is meaningless.

Two properties of this transform matter for how the results are summarised.
Because it is monotone, the **median** per-division rate equals the transform of
the median cumulative rate, so no correction is needed. The **mean** is not
preserved, and is dominated by targets near saturation, which is why medians are
reported throughout.

Above roughly 95% cumulative editing the per-division rate is close to
unidentifiable: many different rates all saturate before generation *n* and land
at the same *p*. At p = 1 the inversion returns r = 1 regardless of the true
rate. Estimates are therefore reported both with and without the saturated
fraction, and the saturated targets are flagged rather than silently included.

## BASELINE: empirical rates from clone 84

### Data and encoding

The input is a cell-by-position matrix of 9,349 cells against 9,248 positions,
comprising 34 barcodes each carrying 272 target sites. Entries are `0`
(observed, unedited), `1` (observed, edited), `?` (site not called in that cell)
and `+` (barcode not captured in that cell). The four categories were verified to
account for every entry, with counts summing exactly to 9,349 × 9,248.

Capture is the dominant source of missing data: 25.4% of entries are `+` and
0.8% are `?`, leaving a median of 7,363 of 9,349 cells (78.8%) observed per
position.

### Per-position editing rate

Editing rate per position was computed over observed entries only, as
`edited / (edited + unedited)`. Both `?` and `+` were excluded from the
denominator rather than being counted as unedited: scoring an uncalled site as
"no edit" would understate every rate by the local dropout, and dropout is not
uniform across positions. Positions observed in no cell carry no rate and were
held as missing rather than as zero.

The mean editing rate across positions is 0.307, and the rate pooled over all
observed entries is 0.3055. Their near-identity indicates that no position is
distorting the average through low coverage.

### The distribution is bimodal

Mean editing rate (0.307) and median (0.031) differ by an order of magnitude,
with an interquartile range of 0.003–0.765. The distribution is not a single
mode with a low tail but two opposing peaks: 33% of positions sit between 0.05%
and 1% editing, while 1,799 positions (19.5%) are edited in at least 99% of
cells, 18 of them in every cell. The saturated peak is what carries the mean
above the median.

Per-barcode means are tight — 0.261 to 0.389 across the 34 barcodes — so this
bimodality is a property of sites within a barcode rather than of some
integrations editing faster than others. Applying a 1% threshold retains between
156 and 179 of 272 sites in every barcode (median 165, 60.9% overall), a spread
of only 23 sites.

### Filtering

Two thresholds are reported. A 0.05% floor removes 528 of 9,248 positions (5.7%)
and moves the mean editing rate from 0.307 to 0.326. A 1% floor removes 3,614
positions and moves it to 0.503.

Both are one-sided, and this is a limitation worth stating: positions edited in
essentially every cell are as uninformative for reconstruction as positions never
edited, since a character shared by all cells resolves no split. A two-sided
0.05%–95% band retains 6,744 positions (72.9%); a 1%–95% band retains 3,658
(about 108 per barcode). The analyses below report the one-sided and two-sided
bands separately rather than choosing between them.

### Per-division rate

Assuming 30 cell generations:

| Position set | n | Median per-division rate |
|---|---|---|
| All positions | 9,248 | 0.103% |
| ≥ 0.05% cumulative | 8,720 | 0.149% |
| ≥ 1% cumulative | 5,634 | 1.776% |
| ≥ 1% cumulative, p < 95% | 3,658 | 0.323% |

The choice of threshold moves the estimate by more than tenfold, entirely
because of how much of the saturated peak it admits. Mean per-division rates are
correspondingly unstable — 4.70% over the 0.05%-filtered set, falling to 0.72%
when saturated positions are excluded — and are not used.

Saturated targets are informative despite being saturated, because they edit
early. Converting each per-division rate into the division at which an unedited
target first edits (geometric, median = ln0.5 / ln(1−r)) shows the bands tiling
lineage depth:

| Cumulative p | Positions | Median per-division rate | Median division of first edit |
|---|---|---|---|
| 1–10% | 1,886 | 0.096% | > 30 |
| 10–50% | 1,070 | 0.872% | > 30 |
| 50–90% | 612 | 4.22% | 16.1 |
| 90–99% | 267 | 10.6% | 6.2 |
| ≥ 99% | 1,799 | 17.5% | 3.6 |

The saturated band edits at a median of division 3.6 of 30, resolving the first
two or three splits of the tree and nothing thereafter. Its rate estimate is a
lower bound, for the identifiability reason given above.

Sensitivity to the generation count is mild and approximately proportional to
1/n. Over the 1%–95% band the median per-division rate runs 0.484% at 20
generations, 0.323% at 30, and 0.242% at 40.

## Prime editing: rates from Choi et al.

Supplementary Table 1 of Choi et al. reports a first-site target editing rate of
0.059 over four days. Taking four days as four cell cycles gives

    r = 1 - (1 - 0.059)^(1/4) = 0.01509

that is, **1.51% per cell cycle**. Because the simulator's prime editing
parameter (`induced_edit_probability_per_cell_cycle`) is itself a per-division
probability, this is directly comparable without further conversion.

A second estimate is available from a mean editing rate of 0.1456 across targets
over 20 divisions, giving r = 0.784%. The two do not reconcile: they differ by
1.93×, and no constant per-division rate reproduces both. Projecting the
four-cycle rate forward to 20 divisions predicts 26.2% cumulative editing against
14.6% observed; projecting the 20-division rate back to four cycles predicts 3.1%
against 5.9% observed.

Three explanations are compatible with the data and are distinguishable with
further measurement:

1. **The rate decays.** Transient or silenced editor expression makes editing
   fastest shortly after delivery, so a four-day window samples the early phase
   while a 20-division window averages over the decay. This would appear as
   cumulative editing plateauing before the end of the experiment.
2. **The two numbers describe different targets.** 0.059 is specifically the
   first site, whereas 0.1456 is a mean across targets. If site 1 edits faster
   than average, these are two rates rather than one rate measured twice, and no
   reconciliation is required. Comparing the first-site rate against the
   across-target mean over the same window would settle this.
3. **Generation time is not one day.** The four-cycle figure assumes 4 days = 4
   divisions. At six cycles r falls to 1.01% and at eight to 0.76%, at which
   point the two estimates agree almost exactly.

The estimate is more sensitive to the assumed cycle count than to anything else:
r is 3.00% at two cycles, 1.51% at four, 1.01% at six and 0.76% at eight.

## Relation to the simulation parameter grid

The benchmark sweeps five editing tiers per technology on an exact doubling
ladder (`generate-recorder-parameter-grid.R`). For the two recorders considered
here:

| Tier | BASELINE (per target per division) | Prime (per cell cycle) |
|---|---|---|
| very_low | 0.00238 | 0.01 |
| low | 0.00476 | 0.02 |
| mid | 0.00952 | 0.04 |
| high | 0.01904 | 0.08 |
| very_high | 0.03808 | 0.16 |

Both recorders used `mid` in the large-scale runs, selected for reconstruction
accuracy rather than to match a measured construct.

The prime editing estimate from Choi et al. falls **inside** the swept range,
between `very_low` and `low`, 59.3% of the way along the doubling ladder and only
1.07× above their geometric midpoint. Expected performance at the measured rate
can therefore be obtained by interpolating between two simulated tiers rather
than extrapolating beyond the grid. The alternative 20-division estimate
(0.784%) falls 1.28× below `very_low` and outside the grid; so does the
four-cycle estimate if generation time is shorter than about 18 hours.

The BASELINE estimates straddle the grid depending on which band is used. The
1%-filtered median (1.776%) sits between `mid` and `high`, 1.87× the `mid` tier
actually simulated. The 1%–95% band median (0.323%) sits between `very_low` and
`low`, with `mid` 2.95× higher. The 0.05%-filtered median (0.149%) falls 1.60×
below `very_low`, outside the grid.

In both cases the simulated `mid` tier runs hotter than the empirical median —
2.65× for prime and 1.87–2.95× for BASELINE — which is expected given that the
tier was chosen to optimise reconstruction accuracy. A benchmark intended to
speak to these specific constructs should report from the bracketing tiers rather
than from `mid`.

## Software availability

Per-position rates, filtering and the per-division conversion are implemented in
`analysis/cli/clone84_editing_rates.R`, which writes a per-position table
(including a `per_division_rate` column) and the distribution figure. The
generation count is a command-line argument (`--generations`, default 30), as are
the rate floor (`--min-rate`) and a coverage flag (`--min-observed`). The
simulation tiers are defined in the benchmark harness at
`rust_cmd/test_harness/tree_building/generate-recorder-parameter-grid.R`.

---

## Notes for the authors (remove before submission)

- **The Choi et al. citation is incomplete.** The value is recorded here as
  "Supplementary Table 1" only; the full reference was not available when this
  was written and must be filled in and checked against the table.
- **Confirm the 4 days = 4 cell cycles assumption.** It is the single most
  influential assumption in the prime editing estimate, and it determines
  whether the measured rate falls inside or outside the simulated grid. If the
  doubling time is known for that experiment, use it.
- **Confirm the 30-generation assumption for clone 84.** It is less critical —
  the estimate scales roughly as 1/n — but it is currently an assumption rather
  than a measurement.
- **Decide which BASELINE band to lead with.** The three filtered medians span
  0.149% to 1.776% and place the construct in different parts of the grid. The
  1%–95% band is the most defensible for reconstruction purposes, since it
  excludes both uninformative extremes, but it is also the most aggressive
  filter and this should be stated plainly rather than buried.
- **The prime editing discrepancy is unresolved.** If it cannot be resolved
  before submission, report both estimates and the three candidate explanations
  rather than selecting one. Reporting only the estimate that lands inside the
  simulated grid would be the wrong choice.
- **Single-rate summaries understate the construct.** BASELINE per-target rates
  span two orders of magnitude, and no single tier reproduces the bimodal
  distribution. If the simulator supports per-target rates, fitting the observed
  distribution is preferable to matching a median.
- Numbers quoted are from the runs recorded in `analysis/results/clone84/` and
  `analysis/results/clone84_min1pct/`.
