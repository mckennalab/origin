# Methods — Single-cell transcriptional simulation

## Overview

Single-cell expression profiles are simulated on top of an already-simulated
cell lineage, so that every profile has a known position in a known tree.
Terminal (sampled) cells receive profiles; internal nodes remain latent
ancestors that are never observed but that carry the state descendants inherit.
All work used R 4.6.0 with SymSim (YosefLab/SymSim, commit `76a674b`, version
string 0.0.0.9000), ape 5.8.1, and Seurat 5.5.0.

## Lineage input and tree preparation

Lineages were produced by the continuous-time Gillespie engine described above.
Only cells alive at the sampling horizon became tips; every division point
separating them was retained as an internal node.

Two transformations were required before SymSim would accept a lineage:

**Collapsing unary nodes.** Restricting a lineage to a subset of terminal cells
leaves internal nodes with a single child. SymSim's recursive edge sampler
indexes children positionally and fails on these with an uninformative error, so
single-child nodes were collapsed with `ape::collapse.singles()`. This preserves
topology and summed branch lengths; only non-branching waypoints are removed.

**Replacing non-positive branch lengths.** SymSim allocates cells to edges in
proportion to branch length and divides by the total, so zero-length branches are
degenerate. These were set to 10^-6.

## Expression simulation

Expression was generated with `SymSim::SimulateTrueCounts()` in continuous mode
(`evf_type = "continuous"`), which propagates extrinsic variation factors (EVFs)
down the tree as a Brownian motion, so cells sharing ancestry share expression
state. Unless stated otherwise: 10 EVFs of which 6 were differential
(`nevf = 10`, `n_de_evf = 6`), sigma = 0.4, `vary = "s"`, and `min_popsize` set
to one tenth of the tip count with a floor of 5. Figures used 400 genes. Gene
counts below roughly 100 cause SymSim to fail with an opaque array-dimension
error, so a minimum of 100 is enforced. Counts were log(1 + x)-transformed before
correlation or clustering.

Note that SymSim must be **attached** rather than referenced through `::`,
because its simulation functions read bundled parameter tables by bare name.

## Assigning profiles to terminal cells

SymSim distributes cells along tree branches in proportion to branch length,
recursing on a shared budget; it does **not** place one cell per tip. Measured on
a 12-tip tree, 10 of 12 terminal edges received a cell at 200 requested cells and
11 of 12 at 600, with the missed edges of ordinary length — this is the sampler's
budget, not a short-branch artifact.

Profiles were therefore assigned explicitly. Each terminal cell took the deepest
cell sampled on its own terminal edge; where that edge was empty, it took the
deepest cell from the nearest ancestral edge. Every profile is labelled with
which rule applied and the proportion assigned from a cell's own terminal edge is
reported per simulation. Cells were over-requested at 100 per tree edge; at 64
tips this gave 62 of 64 assignments from the tip's own edge.

Runs in which two terminal cells resolved to the same sampled cell were rejected.
Such cells receive byte-identical profiles and appear as a block of correlation
exactly 1 that is visually indistinguishable from genuine clade structure; at
lower over-sampling (8 cells per edge) this affected 4 of 64 cells.

## Terminal-cell sampling

Which terminal cells become tips materially affects the lineage signal, and this
is reported rather than left implicit.

A pure-birth lineage sampled **uniformly** at a single horizon is close to
star-shaped: 55% of terminal-cell pairs sat at the maximum cophenetic distance,
mean correlation was 0.693 for the closest decile of pairs against 0.680 for
maximally distant pairs, and the association between cophenetic distance and
expression correlation was rho = -0.15. Because SymSim accumulates EVF variance
along branches, tips coalescing near the root are nearly independent.

Sampling equal numbers of cells from a small number of **disjoint clades** gives
tips shared internal branches and recovers strong structure (rho = -0.64). Clades
were selected as disjoint subtrees containing between one and three times the
per-clade cell count; an upper bound is necessary, since without one a "clade"
can encompass most of the tree and the groups cease to be distinct
(rho = -0.22). Unless stated otherwise, terminal cells were drawn from four
disjoint clades.

## Validation

Across 102 terminal cells, cophenetic distance and expression-profile distance
were positively associated (Spearman rho = 0.32), exceeding all 100
label-permuted replicates (permutation p = 0.01, null mean rho = -0.004). This
test is sensitive to tip count: at 45 tips it has so little power that it fails
roughly one run in three on genuine signal (p = 0.030, 0.069, 0.020 across three
seeds), whereas at 120 tips rho = 0.316-0.320 with p at the permutation floor
across every seed tested.

## Clustering

Simulated profiles were clustered with the standard Seurat workflow:
`NormalizeData` (LogNormalize, scale factor 10,000), `FindVariableFeatures`
(vst), `ScaleData`, `RunPCA`, `FindNeighbors`, `FindClusters` (Louvain), and
`RunUMAP` (uwot 0.2.4).

Three defaults are inappropriate for simulated matrices of this size and were
adjusted: PCA cannot return more components than cells, so `npcs` was capped at
min(genes, cells) - 1; requesting 2,000 variable features from 400 genes is not a
selection, so the request was capped at the gene count; and UMAP's `n.neighbors`
was capped below the cell count. Fifteen principal components were used.

Agreement between clusters and lineage was quantified by adjusted Rand index
(mclust 6.1.2) against lineage groups obtained by cutting the tree, via
average-linkage hierarchical clustering on cophenetic distance, into the same
number of groups the clustering found. Matching granularity in this way scores
assignment rather than rewarding a clustering for choosing the right number of
groups — but it also means a clustering at the wrong scale is not penalised for
that.

On 64 cells drawn from four clades, resolution 0.5 gave 2 clusters matching the
tree's deepest split exactly (ARI = 1.000). Finer granularity degraded: 4
clusters at resolution 1.2 (ARI = 0.560) and 9 at resolution 2.0 (ARI = 0.289).
The major bifurcation is recoverable from expression alone; shallower splits are
not, at this signal level.

## Software availability

Implemented in `analysis/symsim_profiles.R` with drivers
`analysis/cli/generate_symsim_profiles.R`, `cluster_symsim_seurat.R`, and
`plot_seurat_vs_lineage.R`. SymSim is deliberately not a dependency of the
`origin` package; tests in `tests/symsim_tests.R` skip cleanly when it is absent.

---

## Notes for the authors (remove before submission)

- **SymSim citation needs checking.** The installed `DESCRIPTION` gives only
  `0.0.0.9000`, so the commit SHA is pinned here instead. Confirm this matches
  the reference intended (Zhang et al., *Nat Commun* 2019).
- **Keep the terminal-cell sampling subsection.** A reviewer could reasonably
  read clade sampling as choosing the analysis that produces the desired result.
  Reporting the uniform-sampling numbers alongside shows the weak case is
  understood and attributable to tree shape rather than hidden. If space forces a
  cut, drop the parenthetical statistics before dropping the fact that both modes
  exist.
- **scDesign3 is not covered here.** It is a separate generator with a different
  mechanism — covariate-conditioned, fit to a real reference dataset, with no
  tree topology reaching it — and warrants its own subsection if those runs
  appear in the same paper.
- Numbers quoted are from the runs recorded in `analysis/figures/` and
  `analysis/results/`; regenerating with a different seed will shift them.
