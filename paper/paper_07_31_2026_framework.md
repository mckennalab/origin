# ORIGIN paper: proposed simulation framework and results plan

Companion planning document for `paper/paper_07_31_2026.md`. Drafted 2026-07-31.

This document proposes (1) the four results sections, (2) the exact simulations
to run for each, (3) draft prose to surround those simulations, and (4) the
tooling that must be built first. Every parameter named below was verified
against the current source; file:line references point at the code that reads
the parameter. Where a proposed experiment is **not** runnable today, that is
stated explicitly rather than assumed.

---

## 0. Framing: what the paper argues

The existing draft sells ORIGIN as an architecture — separating the
lineage-generating process from the molecular record written on it from the
observation model that recovers it. That separation is the right spine, but
architecture alone is a methods note. The four aims below turn it into a
results paper by using the separation to answer questions that a single-layer
simulator cannot:

1. **Recorder chemistry is not the main variable.** Comparing CRISPR recorder
   technologies on one fixed lineage isolates *architecture* (integrations,
   target count, window, irreversibility, inter-target deletion) from *chemistry*
   (nuclease vs. ABE vs. CBE vs. prime), and we expect architecture to dominate.
2. **Prospective and retrospective records resolve different parts of the same
   tree.** Mitochondrial variants need many divisions of drift to become
   detectable, so they mark early splits; engineered edits accrue at a constant
   hazard and mark recent ones. Their combination should be complementary rather
   than merely additive.
3. **Mitochondrial inheritance model choice is a first-order confound.** Most
   published mtDNA lineage inference implicitly assumes one segregation model.
   We can ask which models produce heteroplasmy statistics compatible with real
   human single-cell data, and whether the choice changes lineage conclusions.
4. **Lineage tools have a detection limit for altered division patterns.** In a
   3D organoid with a disease that skews lineage, we can ask which specific
   division-pattern changes survive reconstruction and sampling — and which are
   invisible.

A unifying claim for the abstract: *the observation model, not the recorder
chemistry, sets the ceiling on what developmental biology a lineage experiment
can recover.*

---

## 1. Ground rules: the two engines and their limits

Everything below has to respect a hard architectural fact — ORIGIN has two
independent recorder implementations that read the same JSON schema but are not
the same model.

| | Native forward simulator | Fixed-tree replay |
|---|---|---|
| Entry | `sim5_code.R -P <json>` via `bash_wrapper_all_combos.sh -d <dir>` | `simulate_physicell_lineage.R`, wrapped by `run_physicell_organoid_pipeline.sh` |
| Tree | simulated (birth/death/differentiation) | imported PhysiCell division log |
| Practical scale | ~2<sup>12</sup> cells (the scale used in the 07/27 draft) | 104,154 cells measured on the day-21 organoid |
| Mitochondria | two-level organelle model (mitochondrion → genomes), fusion/fission, mitophagy, cell-level selection | flat list, fixed copy number, with-replacement resample at division |
| Inter-target deletion | implemented (`interdeletion_dropout_*`) | **not implemented** |
| Prime editing | implemented but crashes (see §2) | **not implemented** (warns, falls back to 1 bp insertions) |
| Cell-type-specific rates | yes, per `cell_type_params` | **no** — one rate block for every branch |
| Observation model | cell sampling + integration/genome recovery sweeps | **none** — ground-truth profiles only |

**Consequence for the paper.** Aims 1 and 3 must run on the native simulator,
because that is the only engine with inter-target dropout, prime editing, and
the organelle model. Aim 4 must run on the replay engine, because it is the only
one that reaches 10<sup>5</sup> cells in 3D. Aim 2 should run on **both** and the
agreement between them becomes a validation figure. Any sentence comparing a
native result directly to a replay result must say which engine produced it;
they are different generative models, not two settings of one.

### Measured runtime and scale anchors

From `output/physicell_organoid_20260730_214332/` (day 21, 2,500 founders,
BASELINE preset, seed 1):

- PhysiCell C++ growth stage: **2 min 36 s** → 104,154 live cells, 101,761 divisions.
- R replay: **40 min 48 s** total (barcode 13 m 37 s, mitochondrial 18 m 21 s,
  mitochondrial output 6 m 37 s).
- **Budget ~45 min and ~700 MB per day-21 replicate.** A 12-run design is ~9 h
  and ~9 GB; with 87 GB free on this machine, keep sweeps under ~40 full runs
  or prune intermediate outputs.

**Reproducibility hazard.** Runs `20260730_171714` and `20260730_200008` both
record `random_seed,1` yet produced 98,487 vs 104,789 cells, because `--jobs`
sets `omp_num_threads` and the OpenMP thread count perturbs the RNG stream
(`run_physicell_organoid_pipeline.sh:260,266`). **Pin `--jobs` explicitly for
every paper run and record it in the manifest.**

**Environment gap.** IQ-TREE is not installed on this machine. It is declared in
`heavy_sim_env.yml` (`iqtree=3.0.1`, `beast2=2.6.3`, `cassiopeia-lineage==2.1.0`,
env name `babette_2_backup`) and invoked with `-nt 40`
(`bash_wrapper_all_combos.sh:132`), so the reconstruction stage is a cluster
job. Build that environment before committing to a sweep size.

---

## 2. Capability audit: what can actually be run

This is the honest core of the plan. The 07/27 draft describes analyses
(MMPC/GTTC clades, F1 purity, ARI, the induced-vs-uninduced EM, the α-noise
sensitivity sweep) that **are not in this repository** — `README.md:651-652`
says so directly. Planning around them as if they exist is the main risk to the
schedule.

### Runnable today, no new code

- Nuclease indel recorders; ABE **and** CBE (`be_conversion_pattern: "A --> G"` /
  `"C --> T"` are the same code path, classified at `sim5_code.R:257`).
- Target array layout: `config` = `"U"` (uniform), `"R"` (random),
  `"S:<first>:<gap>"` (spaced), parsed at `sim5_code.R:161`.
- Editing windows with positional decay (`size`, `decaying`, `close_after_edit`).
- H/M/L rate classes with within-class gamma heterogeneity.
- Irreversible target saturation; one-way inducible recording
  (`editing_induction.timepoint`).
- Inter-target deletion collapse (`interdeletion_dropout_radius`,
  `interdeletion_dropout_prob`) — **native only, same-timestep cut pairs only**.
- Swept axes as JSON arrays: `max_bc_ints_per_cell`,
  `bc_integration_recovery_prob`, `mt_genome_recovery_prob`,
  `mt_allelic_fraction_thresholds`, `sampling_fractions`, `sim_length`.
- Mitochondrial: random partitioning, fusion/fission rates, uniform mitophagy
  (`post_mitotic_mt_deletion_frac`), copy number
  (`starting_mito_per_cell` × `average_genomes_per_mito`), cell-level
  heteroplasmy selection.
- Organoid disease perturbations via XML `user_parameters` (see §6).
- Normalized Robinson–Foulds accuracy (`compare_trees_call_from_bash.r:204`).

### Small fixes (hours)

| Fix | Location | Needed by |
|---|---|---|
| Prime editing crashes: `basepos_erc_nuc_list` used at `sim5_code.R:722`, defined at `:1228` | move the block below `:1228` | Aim 1 |
| `interdeletion_dropout_*` passed as `NULL` when keys absent → `argument is of length zero` on the first deletion timestep | `sim5_code.R:3780-3781`, add defaults | Aim 1 |
| Author a FLARE preset JSON — **no FLARE config exists in the repo or its git history**; the draft's "parameterized by observed edit counts in FLARE" is currently unsupported | new `example_json_params/flare_native.json` | Aim 1 |
| Pin `--jobs` and record it | `run_physicell_organoid_pipeline.sh` | all |

### New code required (the critical path)

Ordered by how many aims they unblock:

1. **Metrics module.** Only normalized RF exists, and it is untested. Add
   triplet/quartet accuracy, clade precision/recall, ARI, F1 purity, and
   **depth-stratified accuracy** (accuracy as a function of node age — this is
   the load-bearing metric for Aim 2). *Unblocks all four aims.*
2. **Run reconstruction on PhysiCell outputs at all.** Grep of the three
   `run_physicell_*.sh` scripts for `iqtree|muscle|compare_trees` returns
   nothing — the replay engine currently produces character matrices that are
   never turned into a tree. *Unblocks Aims 2 and 4.*
3. **Replicate harness and cross-run aggregation.** There is no replicate
   mechanism; `random_seed` is a scalar and `process_results_from_bash.r`
   aggregates per run-id only. Replicates currently mean duplicating JSON files.
   *Unblocks all four aims.*
4. **Missing-data encoding.** Unrecovered integrations read as `0` (unedited),
   not `NA` — the NA-masking logic lives in the dead legacy function
   (`mut_to_scoremat.r:508-520`); the live path is `new_create_one_score_mat`.
   Recovery is a headline axis of Aim 1, and it is currently modeled as
   *false-unedited* rather than *missing*, which biases every recovery sweep.
5. **Mitochondrial observation model in the replay engine** (recovery,
   coverage, VAF thresholding). Replay currently emits ground-truth fractions,
   so mitochondria are evaluated under an optimistic observation model relative
   to barcodes — an unfair comparison in Aim 2 unless fixed or disclosed.
6. **Mitochondrial inheritance models** in `mito_dynamics()` (§5).
7. **Organoid disease knobs**: CLI flags, symmetric-neurogenic division, and a
   mosaic genotype (§6).
8. **State-specific recording hazards.** `cell_state_transitions.csv` is
   written but never read by any R script, so branches are not split at fate
   transitions and one rate block applies to every branch. Any disease that
   changes cycle length will therefore not change the per-branch edit hazard
   correctly. *Needed for Aim 4 to be interpretable.*

---

## 3. Aim 1 — Comparing CRISPR lineage recorder technologies

### Question

On an identical lineage, how much of the tree does each recorder architecture
recover, and what mechanism limits each one?

### Design

Native simulator, single founder to 2,048 terminal cells, 3 seeds per
condition. One shared population model across all recorders so that only the
record differs.

**Table 1 — recorder panel.** Chemistry and architecture are deliberately
crossed so the paper can separate them:

| ID | System | Key parameters | Purpose |
|---|---|---|---|
| R1 | Cas9 nuclease array (FLARE-like) | `nuclease_targets.num_targets: 10`, `config "S:15:20"`, window `size 0`, `close_after_edit: true`, indels on, `interdeletion_dropout_prob: 0.3`, `radius: 200` | canonical multi-cut recorder |
| R2 | R1 with dropout disabled | identical but `interdeletion_dropout_prob: 0` | **isolates inter-target deletion as the mechanism** |
| R3 | Single-cut, many integrations | `num_targets: 1`, `max_bc_ints_per_cell` swept high | dropout-immune nuclease control |
| R4 | ABE (BASELINE-like) | `be_targets.num_targets: 50`, `config "S:10:24"`, window `size 7, decaying true, close_after_edit false`, `be_conversion_pattern "A --> G"` | base-editing array |
| R5 | CBE | R4 with `"C --> T"` and backbone composition rebalanced | does conversion chemistry matter given backbone base content? |
| R6 | Prime editing | `prime_editing_system: true`, `num_unique_prime_editing_guides` swept 4/16/64 | insertion-library recorder (needs the `:722` fix) |

**Swept axes** (all are already JSON arrays, so one file yields the grid):
`max_bc_ints_per_cell` [1, 5, 10, 20, 50]; `bc_integration_recovery_prob`
[0.25, 0.5, 0.75, 1.0]; `sampling_fractions` [0.05, 0.25, 1.0]; `sim_length` at
several stopping points for accuracy-vs-time; per-target rate scaled ×0.25/×1/×4
around each system's default to probe saturation.

Full crossing is 6 × 5 × 4 × 3 × 3 seeds = 1,080 runs. **Prune to a staged
design**: sweep integrations × recovery at one rate and one time (360 runs),
then a rate × time sweep at the best-performing integration count.

### Expected figure (Figure 2)

- **2A** Accuracy (1 − normalized RF) vs. integrations, faceted by recovery
  probability, one line per recorder.
- **2B** Accuracy vs. time with the per-target rate scaling, showing saturation.
- **2C** The R1-vs-R2 contrast alone — the mechanistic panel.
- **2D** Character-level information: fraction of targets edited, number of
  distinct alleles per target, per-character entropy. *(New code.)*

### Draft prose

> **Recorder architecture, not chemistry, sets reconstruction accuracy.**
> We simulated six engineered recorder designs on a common population model,
> holding the lineage-generating process and the observation model fixed so that
> only the molecular record differed. Across the sweep, the number of
> independently inherited integrations accounted for the largest share of
> variance in accuracy, followed by integration recovery probability; the
> identity of the editing chemistry was comparatively minor. Base-editing arrays
> (R4, R5) and single-cut multi-integration nuclease designs (R3) converged to
> similar accuracy once [N] integrations were available, whereas the multi-target
> nuclease array (R1) plateaued below them at every integration count tested.
> Adenine and cytosine base editors were statistically indistinguishable once the
> barcode backbone was rebalanced for target base content, indicating that the
> conversion pattern matters only through the number of editable positions it
> creates.

> **Inter-target deletion, not information content per edit, limits multi-cut
> nuclease recorders.** Disabling inter-target deletion collapse in an otherwise
> identical nuclease recorder (R2) recovered [Δ] accuracy units, closing [X]% of
> the gap to the base-editing array. Because R2 differs from R1 only in
> `interdeletion_dropout_prob`, this attributes the nuclease deficit to the
> physical loss of intervening targets when two cuts resolve together, rather
> than to any difference in per-edit information. The residual gap is consistent
> with irreversible saturation: nuclease targets close permanently after the
> first edit, whereas a decaying base-editing window continues to accept edits at
> nearby positions.

> **Saturation defines an experiment-duration window.** Because target edits are
> irreversible, the expected fraction of edited targets grows as
> 1 − (1 − p)<sup>T/τ</sup> for per-division edit probability p and cycle length
> τ. Recorders tuned for a short experiment exhaust their targets and lose
> resolution on late divisions, while recorders tuned for a long experiment carry
> too little signal early. [Panel 2B] locates this optimum as a function of the
> intended experiment length, giving a design rule: choose p ≈ [value] / expected
> number of divisions.

**Caveats to state in Methods.** Cas9 and Cas12a are not distinguished by any
code — `physicell_adapter.recorder_system` is a free-text label only, so any
Cas12a-specific claim must be encoded through window size, spacing, and rate,
and described that way. hgRNA self-targeting and integrase/recombinase recorders
are absent from the codebase and should appear only in the Discussion as
extensions. Inter-target deletion pairs only cuts occurring in the *same*
timestep, which under-models sequential collapse — report `time_inc`.

---

## 4. Aim 2 — Mitochondrial vs. CRISPR lineage tracing

### Question

Do prospective and retrospective records resolve the same tree, or different
parts of it — and is their combination complementary?

### Design

Both engines, same trees, three character sets per tree: barcode-only,
mitochondria-only, and joint (`combine_mt_bc: true` natively;
`combined_lineage_feature_matrix.rds` on the replay path).

Axes: VAF threshold [0, 0.01, 0.05, 0.10, 0.20]; `mt_genome_recovery_prob`
[0.1, 0.25, 0.5, 1.0]; copy number via `starting_mito_per_cell` ×
`average_genomes_per_mito` [50×5 = 250, and a reduced 20×5 = 100];
`bc_integration_recovery_prob` as in Aim 1; several `sim_length` stopping points.

**The distinguishing analysis is depth-stratified accuracy**: score each
internal node of the reconstruction by the age of the corresponding true split,
and plot accuracy against node age separately for each modality. This is the
figure that makes the complementarity argument concrete rather than assertive.

### Preliminary evidence (measured, this session)

From the completed day-21 organoid replay (104,154 cells; 32 mt genomes/cell;
5 × 50 BASELINE targets), clustering a random 2,000-cell subsample against the
2,499 ground-truth founder clones:

| Character set | Callable characters | ARI at coarse (k=50) | ARI at clone resolution (k=991) |
|---|---|---|---|
| Mitochondrial, VAF ≥ 0 | 2,769 | 0.098 | **0.924** |
| Mitochondrial, VAF ≥ 0.10 | 1,022 | 0.001 | 0.011 |
| BASELINE barcode | 573 | 0.000 | 0.215 |

Mean callable variants per cell: 12.85 at VAF ≥ 0, 7.30 at ≥ 0.05, 4.85 at
≥ 0.10, 3.47 at ≥ 0.20. Barcode: 107.8 of 640 characters edited per cell (17%),
92% of characters edited somewhere in the population.

Two things follow. First, **low-frequency mitochondrial variants carry nearly
all of the clonal signal** — a conventional 10% VAF threshold destroys it
(0.924 → 0.011). That is a strong, quantitative, publishable result and it
directly contradicts standard mtscATAC practice. Second, these numbers are from
a crude average-linkage binary clustering, not the paper's IQ-TREE pipeline;
they are a pilot that justifies the experiment, not a result.

### Expected figure (Figure 3)

- **3A** Accuracy vs. VAF threshold, one line per modality. The headline panel.
- **3B** Depth-stratified accuracy: accuracy vs. true node age, barcode vs.
  mitochondria vs. joint. *(New code.)*
- **3C** Joint-vs-best-single across the recovery grid — is combination
  additive or complementary?
- **3D** Callable characters per cell vs. threshold, with the published
  mtscATAC range overlaid.

### Draft prose

> **The two modalities resolve different epochs of the same tree.** Replaying
> both recorders along identical branches isolates the difference between them to
> the record itself. Engineered edits accumulate under a constant per-unit-time
> hazard and therefore mark divisions uniformly until targets saturate;
> mitochondrial variants must first drift to a detectable heteroplasmy, which
> requires many divisions of segregation, so they preferentially mark early
> splits. Stratifying reconstruction accuracy by the age of the true split
> [Figure 3B] shows this directly: mitochondrial characters recovered [X]% of
> splits older than [T] but only [Y]% of splits in the final [n] divisions, while
> the engineered recorder showed the opposite profile. Combining the two
> modalities improved accuracy over the better single modality by [Δ], and the
> improvement was concentrated at intermediate node ages where neither record is
> individually informative — the signature of complementary rather than additive
> information.

> **Variant-calling thresholds, not variant abundance, limit mitochondrial
> lineage resolution.** [Figure 3A] Lowering the allele-fraction threshold from
> 0.10 to 0 increased clone-level agreement from [X] to [Y] while increasing the
> callable variant count from [a] to [b] per cell. Because low-frequency variants
> are precisely those that have not yet been homogenized by segregation, they
> retain the most recent lineage information; discarding them is discarding the
> signal. This suggests that mitochondrial lineage studies are limited less by
> the mutation rate than by the confidence with which a low-heteroplasmy variant
> can be called against sequencing noise, and it motivates explicit error
> modeling rather than a fixed threshold.

**Caveats to state.** The replay engine applies no recovery or thresholding to
mitochondria (ground-truth fractions only) while barcodes go through an
integration-recovery model, so the comparison is currently unfair to barcodes;
build tooling item 5 before running this aim, or report both under a
ground-truth observation model. Barcode integration dropout currently reads as
false-unedited rather than missing (tooling item 4), which flatters barcodes at
low recovery. Both must be fixed for Figure 3 to mean anything.

---

## 5. Aim 3 — Mitochondrial inheritance models against real data

### Question

Which mitochondrial segregation models produce heteroplasmy statistics
compatible with observed human single-cell mtDNA data, and does the choice
change lineage conclusions?

### Status: this aim is currently blocked

The maintained simulator implements exactly **one** inheritance model. The
native engine hard-stops on anything else:

```r
# sim5_code.R:2006-2008
if(inheritance_pattern != 'random'){
  stop("Only inheritance_pattern = 'random' is currently implemented.")
}
```

The replay engine offers a single fixed-size, with-replacement resample
(`physicell_mito.R:276-283`). The Models I/J/K/L and the "seven unique
mitochondrial inheritance regimes" in the 07/27 draft were never in this
repository — a search of the full git history for a biased or bottleneck branch
in `mito_dynamics` returns nothing. **They must be written, not restored.**

### Implementation plan

All changes are local to `mito_dynamics()` (`sim5_code.R:1992`), branching at
`:2079` and relaxing the guard at `:2006`:

| ID | Model | Change |
|---|---|---|
| M0 | Random binomial partition of mitochondria | exists |
| M1 | Biased partition, p ≠ 0.5 | partition probability is hard-coded 0.5 at `sim5_code.R:2082` — expose it |
| M2 | Fixed bottleneck to *k*, then re-expand | new branch |
| M3 | Progressive/developmental bottleneck, *k* varying by cell type | new branch, per-`cell_type_params` *k* |
| M4 | Clustered/network segregation of fused units | new branch over `mito_to_genome_map` |
| M5 | Relaxed replication | replication is exact doubling at `sim5_code.R:2062`; replace with random replication to a target copy number |
| M6 | Selective mitophagy | `post_mitotic_mt_deletion_frac` removal is currently uniform (`:2152-2195`); weight by variant severity |

Estimated 1–2 days for M1–M4, plus a day for M5–M6 and tests.

### Validation against real data — the actual contribution

This is what makes the aim more than a parameter sweep. Compute, for each
model, summary statistics that have published human counterparts, and ask which
models fall inside the observed range:

1. **Informative variants per cell** at standard thresholds, against mtscATAC
   observations (Ludwig 2019; Lareau 2021). *Current simulation: 12.85 at
   VAF ≥ 0, 4.85 at ≥ 0.10 — plausible in magnitude.*
2. **Heteroplasmy spectrum shape** — the distribution of VAF across cells and
   its bimodality (Weng 2024).
3. **Normalized heteroplasmy variance growth**,
   V = Var(f) / (f̄(1 − f̄)), which grows as 1 − (1 − 1/N<sub>e</sub>)<sup>g</sup>
   over *g* generations. **This is the discriminating statistic**: fitting an
   effective bottleneck size N<sub>e</sub> from simulated cells and comparing it
   with published human estimates gives a single, quantitative,
   model-discriminating number. M0 and M2 will give very different N<sub>e</sub>
   for the same nominal copy number.
4. **Variant-sharing structure** — the clone-size distribution of
   variant-defined groups, compared with reported clonal substructure.

**Required input the repo does not have.** There is no mtDNA data anywhere in
this project (`data/` contains only `LARRY_adata_preprocessed.h5ad`, which is
referenced by nothing). A public mtscATAC dataset must be obtained and reduced
to these four summary statistics before this aim can be written.

**Known model gaps to disclose.** The simulated mitochondrial genome is
uniform random ACGT, not rCRS; there is no strand asymmetry and no
mtDNA-specific mutational spectrum; the only mtDNA-flavored parameter is a high
HKY κ = 25 and `heteroplasmy_variant_transition_prob = 0.95`. Mutation rates
(`mt_sub_model_params [25, 5e-06, 2e-07]`) are assumed, not fitted —
`README.md:243` says so. If Aim 3 claims compatibility with real data, the
mutation spectrum is the weakest link and should be addressed or scoped.

### Expected figure (Figure 4)

- **4A** Schematic of M0–M6 segregation.
- **4B** Simulated vs. observed: the four summary statistics, with the published
  range shaded. The compatibility verdict panel.
- **4C** Fitted effective bottleneck size N<sub>e</sub> by model, against
  published estimates.
- **4D** Lineage accuracy by inheritance model — does the confound matter?

### Draft prose

> **Segregation model choice changes heteroplasmy statistics by more than the
> mutation rate does.** We implemented seven mitochondrial inheritance regimes
> spanning random partitioning, biased partitioning, fixed and progressive
> bottlenecks, clustered network segregation, relaxed replication, and selective
> mitophagy, and simulated each under an otherwise identical population model and
> mutation process. Across models, the number of callable variants per cell
> varied by [X]-fold and the normalized between-cell heteroplasmy variance by
> [Y]-fold, exceeding the variation produced by an order-of-magnitude change in
> the mutation rate. Because most mitochondrial lineage inference implicitly
> assumes one segregation model, this variation is a confound rather than a
> nuisance parameter.

> **Only [n] of seven models are compatible with observed human single-cell
> mitochondrial data.** Comparing simulated summary statistics with published
> mtscATAC observations [Figure 4B], models [list] reproduced both the observed
> per-cell variant count and the shape of the heteroplasmy spectrum, whereas
> [list] produced [too many low-frequency variants / insufficient between-cell
> variance]. Fitting an effective bottleneck size to each simulation recovered
> N<sub>e</sub> = [range], bracketing published human estimates only for models
> [list]. Notably, the models that best matched the data were not the ones that
> gave the best lineage reconstruction accuracy [Figure 4D], indicating that
> realistic mitochondrial segregation is *less* informative than the idealized
> random-partitioning assumption commonly used in simulation studies — and that
> published benchmarks based on random partitioning are optimistic.

That last sentence is the paper's sharpest claim if the result holds. It is also
falsifiable, which is the point.

---

## 6. Aim 4 — 3D neural organoid with a lineage-skewing disease

### Question

For a disease that alters progenitor division patterns, which specific changes
are detectable from lineage data, with which recorder, at what sampling depth?

### Measured wild-type baseline (from the completed day-21 run)

| Quantity | Value |
|---|---|
| Live cells at day 21 | 104,154 (101,761 divisions) |
| Composition | 57.4% neural progenitor, 24.8% neuron, 17.7% radial glia, 0.1% epithelial |
| Founder clones represented | 2,499 of 2,500 |
| Clone size | median 22, mean 41.7, max 599 |
| Clone-size skew | top 10% of clones hold 44.7% of cells |
| Fate transitions | RG→NPC 36,188; NPC→neuron 25,808; NE→RG 14,210; iPSC→NE 8,280 |

This is a strong baseline: the clone-size distribution is already highly skewed
under the null, so a disease effect has to be measured against real variance
rather than against a uniform expectation.

### Disease models

The division logic is three lines that matter (`custom.cpp:394-423`): only
radial glia divide asymmetrically, the parent **always** retains `radial_glia`,
and the daughter becomes an NPC with probability
`radial_glia_differentiating_daughter_probability` (default **1.0**) after
`neurogenesis_start` (default 7200 min = day 5).

**Runnable today with a one-line XML edit** (these are `user_parameters`, read at
runtime, and the pipeline's perl patch does *not* overwrite them):

| ID | Phenotype | Parameter change |
|---|---|---|
| D1 | Premature neurogenic switch (microcephaly-like) | `neurogenesis_start` 7200 → 2880 (day 2) |
| D2 | Progenitor expansion (macrocephaly-like) | `radial_glia_differentiating_daughter_probability` 1.0 → 0.5, and/or `neurogenesis_start` → 14400 |
| D3 | Delayed differentiation | `neuroepithelial_to_radial_glia_rate`, `ipsc_to_neuroepithelial_rate` reduced |
| D4 | Increased progenitor death | `base_apoptosis_rate` 1e-7 → 1e-5 |

I verified the pipeline's perl patch (`run_physicell_organoid_pipeline.sh:264-279`)
rewrites only `max_time`, `omp_num_threads`, `random_seed`, `initial_cells`,
`initial_organoid_radius`, `ipsc_epithelial_probability`,
`max_epithelial_cells`, the three `phase_transition_rates` blocks, and the four
cycle/differentiation rates. **None of D1–D4 is overwritten.**

One operational wrinkle: the pipeline copies the project XML into the run's
build directory and patches it there, and there is no per-run override flag, so
the values must be changed in
`physicell_projects/neural_organoid_lineage/config/PhysiCell_settings.xml`
before each run — which would silently affect every later run. For a
reproducible sweep, either keep one project copy per disease model, or add
pass-through flags to the pipeline following the pattern already used for
`ORGANOID_*` environment variables. **Adding the flags is the better option**
and takes about an hour; it also puts the disease parameters into
`pipeline_manifest.csv`, which is where the paper's provenance should come from.

**Two important structural limits.**

1. **True microcephaly is not reachable without C++.** Because the parent always
   stays radial glia, the RG pool can never *shrink* through division — p = 0
   doubles it, p = 1 holds it constant. There is no symmetric-neurogenic mode, so
   the canonical microcephaly mechanism (premature symmetric neurogenic division
   depleting the progenitor pool) cannot be simulated. **This is a ~3-line
   change** at `custom.cpp:422`: under a new probability, also call
   `transition_cell(pParent, daughter_state)`. Given that the canonical
   lineage-skewing disease is exactly this mechanism, I recommend making this
   change; D1 is otherwise a proxy for it, not the thing itself.
2. **A mosaic mutant subclone is not supported.** Every rate is a global read
   through `parameter_or()`; `custom_data` holds only `founder_ID`,
   `state_start_time`, `transition_count`, none of which is consulted by any
   rate function. Adding a `genotype` field (`conserved="false"` so it is not
   halved at division), marking a fraction of founders in `setup_tissue()`, and
   gating the cycle rate and self-renewal probability on it is ~20–40 lines
   confined to `custom.cpp` plus XML. This unlocks **D5**, a somatic
   self-renewal-advantage subclone (tuberous-sclerosis/focal-cortical-dysplasia-like),
   which is the most interesting model because it is detectable specifically as
   a clone-size outlier and is invisible to bulk composition assays.

### Detection tasks — the actual science

For each disease model, four ground-truth division-pattern statistics, each
recomputed from progressively more degraded observations:

- **T1 Clone-size distribution.** Gini coefficient, top-decile share, maximum
  clone size, recovered from founder-label barcodes.
- **T2 Within-clone fate composition.** Neuron fraction per clone — does the
  disease change the fate ratio a clone produces?
- **T3 Sister-cell fate concordance.** For each division, are the two products
  the same fate? This is the *direct* readout of symmetric vs. asymmetric
  division and the statistic most tightly coupled to the disease mechanism.
- **T4 Divisions-from-founder to neuron.** Lineage depth distribution; a
  premature switch shortens it.

Each task is evaluated on: (a) the ground-truth tree, (b) a barcode
reconstruction, (c) mitochondrial clades, (d) the joint matrix — crossed with
sampled cells [2,000 / 5,000 / 20,000 of ~104,000] and integration recovery
[0.5, 1.0]. Note that 10<sup>5</sup> leaves is far beyond IQ-TREE; subsampling is
required and is itself realistic, since organoid scRNA-seq recovers ~10<sup>4</sup>
cells.

**Design: 5 disease models (incl. WT) × 3 seeds = 15 runs ≈ 11 h, ~11 GB.**
Pin `--jobs`. Recording replay can be reused across sampling depths without
re-running PhysiCell.

### Expected figure (Figure 5)

- **5A** 3D renderings, WT vs. disease, colored by fate and by clone.
- **5B** Ground-truth division-pattern statistics T1–T4, WT vs. each disease —
  the effect sizes to be detected.
- **5C** **The detectability matrix**: disease model × detection task ×
  observation modality, colored by statistical power at a fixed sample size.
  This is the paper's practical deliverable.
- **5D** Power vs. number of cells sampled, for the hardest detectable task.

### Draft prose

> **A three-dimensional organoid model with ground-truth fate transitions.**
> We simulated cortical organoid development from 2,500 individually tracked
> iPSCs in PhysiCell, through neuroepithelial, radial-glial, and neural-progenitor
> states to postmitotic neurons, with oxygen and nutrient fields and
> mechanically realistic 3D growth. At day 21 the wild-type organoid contained
> 104,154 cells (57.4% neural progenitors, 24.8% neurons, 17.7% radial glia)
> derived from 2,499 surviving founder clones with a median clone size of 22 and a
> maximum of 599; the top decile of clones accounted for 44.7% of all cells.
> Because the simulator records every division, every fate transition with its
> time and position, and the founder of every cell, the true division pattern is
> known exactly, and lineage-based estimates of it can be scored rather than
> merely compared.

> **Disease-driven changes in division mode are detectable, but changes in
> division timing are not.** We perturbed the balance between self-renewing and
> differentiating radial-glial divisions and the onset of neurogenesis, producing
> organoids whose day-21 cell-type proportions differed by less than [X]
> percentage points from wild type — differences that would be undetectable by
> composition assays alone. Lineage reconstruction from the engineered recorder
> recovered the shift in [T1/T3] at [n] sampled cells with power [p], while
> [the timing perturbation] remained undetectable at every sampling depth tested.
> The distinguishing feature of the detectable perturbations is that they change
> the *distribution* of clone sizes and sister-cell fate concordance, whereas
> timing perturbations shift all clones together and leave the distribution's
> shape intact. Lineage recording therefore reports on division *mode* far more
> sensitively than on division *rate*.

> **Mosaic clonal advantage is the easiest phenotype to detect and the hardest
> to see any other way.** A somatic subclone with a self-renewal advantage
> introduced into a single founder produced [X] outlier clones detectable at [n]
> sampled cells, with no change in bulk cell-type composition. Because the
> phenotype is defined by the size distribution of individual clones rather than
> by any population average, it is invisible to composition-based assays and
> visible to lineage recording at modest sampling depth — the clearest case in
> our simulations where lineage tracing is not merely confirmatory but necessary.

**Caveats to state.** Recording hazards are not state-specific on the replay
path (`cell_state_transitions.csv` is written but never read), so a disease that
changes cycle length will not change the per-branch edit hazard correctly —
tooling item 8 should land before D1/D2 are interpreted quantitatively. Death
times are not exported (only a boolean `alive`), so extinct lineages have
unrecoverable terminal branch lengths. `state_start_time` and `transition_count`
are inherited verbatim at division, so `time_in_state` can exceed a cell's own
age — do not use it as a covariate without correcting. The organoid project is
untracked in git and was created after the 2026-07-29 bug audit, so it is
entirely un-audited; it needs a regression test before it carries paper claims.

Also worth correcting in the current draft: `paper_07_31_2026.md:58` states the
PhysiCell export lacks per-cell phenotype and substrate exposure. That is true
of the `tumor_3D_lineage` project but **false for the organoid**, which exports
`cell_type`, `oxygen`, `nutrient`, `developmental_pseudotime`, `founder_ID`,
`time_in_state`, and `transition_count`. Scope that sentence to the tumor
project and to "cell-cycle phase and death time", which really are missing.

---

## 7. Suggested paper structure

| Section | Content | Status |
|---|---|---|
| Fig 1 | Architecture: lineage generation / recording / observation | prose exists |
| Fig 2 | Aim 1 — recorder technology comparison | runnable after small fixes |
| Fig 3 | Aim 2 — mitochondrial vs. engineered, depth-stratified | needs metrics + mito observation model |
| Fig 4 | Aim 3 — inheritance models vs. real data | needs new models + external data |
| Fig 5 | Aim 4 — organoid disease detectability matrix | needs XML sweep + optional C++ |
| Fig 6 | scDesign3 covariate-linked expression | validated at smoke scale |

**Recommended order of execution**, because it front-loads the shared
dependencies and puts the most novel result where it can still be cut:

1. Tooling items 1–4 (metrics, PhysiCell reconstruction, replicates, NA masking).
2. Aim 1 — cheapest, and the first real result.
3. Aim 4 XML sweep — long wall-clock, start it early and let it run.
4. Aim 2 — needs tooling item 5.
5. Aim 3 — largest new-code and external-data burden; the most likely cut.

Aim 3 is the highest-risk section. If the external mtDNA comparison cannot be
assembled, it degrades gracefully into "inheritance model choice is a confound
for lineage inference" (Figure 4D alone), which still supports the paper's
thesis without claiming real-world compatibility.

---

## 8. Threats to validity to address before submission

- **No calibration anywhere.** Every biological rate in the project is an
  assumed default; `README.md:243` and `ORGANOID_SIMULATION.md:252` say so. The
  BASELINE edit rate 0.00952 is a constant carried over from an unrelated smoke
  test, and `ORGANOID_SIMULATION.md:107` states it is "not a fitted estimate".
  Either fit the headline parameters or state prominently that all results are
  relative comparisons under matched assumptions.
- **FLARE has no configuration file.** The claim that recorders are
  "parameterized by observed edit counts in FLARE and BASELINE" is currently
  unsupported for FLARE — no such config exists in the repo or its history.
- **One untested metric.** Normalized RF is the only accuracy measure and
  `compare_trees_call_from_bash.r` has no test coverage. It also treats both
  trees as unrooted while the ground truth is rooted and may contain unary
  nodes; state the normalization convention explicitly in Methods.
- **Insertions are encoded as floating-point decimal digits**
  (`BUG_AUDIT.md:68-71`), so long insertions lose precision. This bounds how
  long a prime-editing guide can be simulated safely — relevant to R6.
- **Sequencing error, amplicon capture, and read depth are not modeled at
  all.** Since Aim 2's headline claim concerns low-frequency variant calling,
  the absence of an error model is a direct threat to that claim and should be
  either implemented or scoped as an explicit assumption.
