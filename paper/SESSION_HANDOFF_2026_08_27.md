# Session handoff — ORIGIN lineage simulation paper framework

**Documented:** 2026-08-27
**Work performed:** 2026-07-31, against repo state at commit `5dd66c4` (branch `aaron_3d`)
**Measurements taken from:** `output/physicell_organoid_20260730_214332/`

> ⚠️ **Read §9 first if you are resuming after 2026-07-31.** The repository has
> changed substantially since this analysis. Several gaps identified below may
> already be closed, and **all `file:line` references are stale** for
> `sim5_code.R`, `physicell_lineage.R`, and `physicell_mito.R`.

---

## 1. What was asked

> "Propose a framework for the lineage simulation paper in
> `papers/paper_07_31_2026.md` detailing simulations results to run and
> paragraphs around those simulations, focusing on (1) comparing various CRISPR
> lineage recorder technologies (2) comparing mitochondrial and CRISPR lineage
> tracing technologies, (3) comparing various models of mitochondrial
> inheritance against mitochondrial lineage tracing to understand which are
> compatible with simulation and real world results, and (4) simulating 3d
> development of a neuronal organoid with the lineage simulation tools for a
> disease that skews lineage and determine what cell division patterns are
> detectable with lineage tools."

The named path was `papers/` (plural); the actual file is `paper/paper_07_31_2026.md`.

## 2. What was produced

| File | Status |
|---|---|
| `paper/paper_07_31_2026_framework.md` | **Created** — 697-line framework: 4 aims, simulation matrices, draft prose, capability audit, critical path |
| `paper/SESSION_HANDOFF_2026_08_27.md` | This document |
| `paper/paper_07_31_2026.md` | **Untouched** — the active manuscript draft |
| `paper/paper_2026_07_27_markdown.md` | **Untouched** — earlier draft with full Methods |

**Decision made:** the framework was written as a *companion* file rather than
into `paper_07_31_2026.md`, because that file is an active manuscript (abstract,
intro, results, references) written the same day and overwriting it would have
destroyed work. The offer to merge it in — or use it to replace the Results
section — was left open and **has not been acted on**.

## 3. How the analysis was done

Four parallel `Explore` subagents mapped the codebase (reconstruction/benchmarking,
PhysiCell organoid, mitochondrial models, CRISPR recorder models), plus direct
verification of every load-bearing claim. The subagent reports are gone; their
substance is preserved in the framework doc and below.

---

## 4. Key findings — capability audit (as of 2026-07-31)

The single most important conclusion: **three of the four aims were blocked on
missing code, not on compute.**

### 4.1 The evaluation half barely existed

- **Normalized Robinson–Foulds was the only accuracy metric**, via
  `phangorn::RF.dist(..., normalize = TRUE)` in `compare_trees_call_from_bash.r`.
  It had **no test coverage**. `accuracy = 1 - rf_dist`.
- **Absent:** quartet distance, triplets, kNN purity, clade precision/recall,
  ARI, cophenetic correlation, AUROC, F1 purity.
- **Absent reconstruction methods:** neighbor-joining, BIONJ, UPGMA, maximum
  parsimony, Camin-Sokal, Cassiopeia (greedy or ILP), hierarchical clustering.
  IQ-TREE was the only wired method, in two modes: `-m HKY` on MUSCLE-aligned
  nucleotide FASTA, and `-st BIN -m MF` on binary character matrices.
- BEAST2 via `make_babette_tree.r` exists but is **disconnected and untested**
  (the `source()` is commented out; nothing calls it).
- **The MMPC/GTTC clade analysis, F1 purity, ARI, and the induced-vs-uninduced
  EM classifier described at length in `paper_2026_07_27_markdown.md`
  (lines 288–373) are NOT in this repository.** `README.md:651-652` stated this
  outright. `compare_celltype_probs.r` is not the EM — it is a hard-coded
  plotting script whose inputs are never produced.
- **No replicate mechanism.** `random_seed` is a scalar; replicates meant
  duplicating JSON files. `process_results_from_bash.r` aggregates per run-id
  only, writing `output/param_results_files/<run_id>/stacked_results.csv`.

### 4.2 Nothing reconstructed trees from PhysiCell output

Grep of `run_physicell_*.sh` for `iqtree|muscle|compare_trees|rf` returned
nothing. The replay engine produced character matrices that never became trees.

### 4.3 Only one mitochondrial inheritance model existed

```r
# sim5_code.R (line ~2006 at the time)
if(inheritance_pattern != 'random'){
  stop("Only inheritance_pattern = 'random' is currently implemented.")
}
```

The replay engine (`physicell_mito.R`) offered a single fixed-size,
**with-replacement** resample at division. The Models I/J/K/L and "seven unique
mitochondrial inheritance regimes" from the 07/27 draft **were never in the
repository** — a search of the full git history for a biased/bottleneck/clustered
branch in `mito_dynamics` returned nothing. They must be written, not restored.

Configurable at the time: fusion/fission rates, uniform mitophagy
(`post_mitotic_mt_deletion_frac`), copy number (`starting_mito_per_cell` ×
`average_genomes_per_mito`), cell-level heteroplasmy selection, and (replay)
`--mt-genomes-per-cell`.

**No calibration to real data anywhere.** `data/` held only
`LARRY_adata_preprocessed.h5ad`, referenced by nothing. The mt genome is uniform
random ACGT (not rCRS); no strand asymmetry; no mtDNA-specific spectrum.

### 4.4 CRISPR recorders

Two independent implementations reading the same JSON schema: the native
forward simulator and the fixed-tree replay adapter. They are **not the same
model**.

Available: nuclease indels; ABE **and** CBE (same code path, selected by
`be_conversion_pattern`); target layouts `"U"`/`"R"`/`"S:<first>:<gap>"`;
editing windows with positional decay and `close_after_edit`; H/M/L rate
classes; irreversible saturation; one-way induction; inter-target deletion
collapse (`interdeletion_dropout_radius`/`_prob`, **native only, same-timestep
cut pairs only**).

Absent: hgRNA self-targeting; integrase/recombinase recorders; per-integration
heterogeneous target layouts; dox on/off windows; sequencing error / amplicon
capture; character-entropy metrics.

**Cas9 and Cas12a are not distinguished by any code** —
`physicell_adapter.recorder_system` is a free-text label only.

**No FLARE preset existed** in the repo or its git history, despite the draft
claiming recorders "parameterized by observed edit counts in FLARE and BASELINE."

### 4.5 New bugs found, not in `BUG_AUDIT.md` at the time

1. Prime editing crashed natively: `basepos_erc_nuc_list` was read ~506 lines
   before its assignment. (**May now be fixed** — see §9.)
2. `interdeletion_dropout_{radius,prob}` passed as `NULL` when the JSON keys are
   absent → `argument is of length zero` on the first timestep with a deletion.
   All shipped `physicell_*.json` presets omitted these keys.
3. Missing-data encoding: unrecovered integrations read as `0` (unedited), not
   `NA`. The NA-masking logic lives in a **dead** legacy function
   (`create_one_score_mat`); the live path is `new_create_one_score_mat`. This
   biases every integration-recovery sweep toward *false-unedited*.
4. `neural_induction_start` inconsistency: XML sets `0`, the C++ fallback is
   `2880`, and the non-BASELINE recording JSON sets
   `editing_induction.timepoint = 2880` — model onset and recorder onset are
   silently decoupled.
5. Reproducibility: `--jobs` sets `omp_num_threads`, and the OpenMP thread count
   perturbs the RNG. Two runs both recorded `random_seed,1` yet produced 98,487
   vs 104,789 cells. **Pin `--jobs` on every paper run.**

### 4.6 The organoid path was in the best shape

Disease-relevant knobs are XML `user_parameters`, read at runtime, editable
**without recompiling**: `neurogenesis_start`, `neuron_production_start`,
`gliogenesis_start`, `radial_glia_differentiating_daughter_probability`,
`astrocyte_daughter_probability`, `ipsc_to_neuroepithelial_rate`,
`neuroepithelial_to_radial_glia_rate`, `neuroepithelial_cycle_rate`,
`base_apoptosis_rate`, `hypoxic_apoptosis_rate`, `hypoxia_threshold`.

I verified the pipeline's perl patch rewrites only `max_time`,
`omp_num_threads`, `random_seed`, `initial_cells`, `initial_organoid_radius`,
`ipsc_epithelial_probability`, `max_epithelial_cells`, the three
`phase_transition_rates` blocks, and the four cycle/differentiation rates —
**so none of the disease knobs are clobbered.**

**Two structural limits (both requiring C++):**

- `organoid_division()` always retains the parent as `radial_glia`. So
  `radial_glia_differentiating_daughter_probability` = 1 holds the RG pool
  constant and 0 doubles it; **the pool can never shrink.** There is no
  symmetric-neurogenic mode, so canonical microcephaly (premature symmetric
  neurogenic division depleting the progenitor pool) is unreachable. Fix is
  ~3 lines: under a new probability, also `transition_cell(pParent, daughter_state)`.
- **Mosaic mutant subclone unsupported.** Every rate is a global read through
  `parameter_or()`; `custom_data` holds only `founder_ID`, `state_start_time`,
  `transition_count`, none consulted by any rate function. Adding a `genotype`
  field (`conserved="false"`), marking founders in `setup_tissue()`, and gating
  cycle rate / self-renewal on it is ~20–40 lines.

Also: `cell_state_transitions.csv` is **written but never read** by any R
script, so recording hazards are not state-specific and branches are not split
at fate transitions. Death times are not exported (only boolean `alive`).
`state_start_time`/`transition_count` are inherited verbatim at division, so
`time_in_state` can exceed a cell's own age.

---

## 5. Empirical measurements (reproducible)

All from `output/physicell_organoid_20260730_214332/` — day 21, 2,500 founders,
BASELINE preset (5 integrations × 50 ABE targets), 32 mt genomes/cell, seed 1.

### Organoid baseline

| Quantity | Value |
|---|---|
| Live cells at day 21 | 104,154 (101,761 divisions) |
| Composition | 57.4% neural progenitor, 24.8% neuron, 17.7% radial glia, 0.1% epithelial |
| Founder clones represented | 2,499 of 2,500 |
| Clone size | median 22, mean 41.7, max 599 |
| Clone-size skew | top 10% of clones hold 44.7% of cells |
| Fate transitions | RG→NPC 36,188; NPC→neuron 25,808; NE→RG 14,210; iPSC→NE 8,280; iPSC→epithelial 100 |

### Recorder information content

| Metric | Value |
|---|---|
| Mito VAF matrix | 104,154 cells × 71,196 variants, 1,338,788 nonzero |
| Mito variants/cell | 12.85 at VAF ≥ 0; 7.30 at ≥ 0.05; 4.85 at ≥ 0.10; 3.47 at ≥ 0.20 |
| Mito variants shared by ≥ 5 cells | 24,925 (VAF ≥ 0) |
| Founder-private variants | 46.9% in exactly 1 founder; 53.1% in ≥ 2; max 50 founders share one |
| Barcode matrix | 104,154 cells × 640 characters |
| Barcode burden | 107.8 edited characters/cell (17%); 92% of characters edited somewhere |
| Target layout | 250 primary targets (5 × 50) |

### Runtime

| Stage | Time |
|---|---|
| PhysiCell C++ growth | 2 min 36 s |
| Barcode replay | 13 min 37 s |
| Mitochondrial replay | 18 min 21 s |
| Mitochondrial output | 6 min 37 s |
| **Total R replay** | **40 min 48 s** |

**Budget ≈ 45 min and ≈ 700 MB per day-21 replicate.** 87 GB free at the time.

### Pilot analysis — clone recovery, mito vs barcode

2,000-cell random subsample, restricted to clones with ≥ 5 cells (991 clones),
average-linkage clustering on binary Jaccard distance, scored by ARI against
ground-truth `founder_ID`:

| Character set | Callable characters | ARI at k=50 | ARI at k=991 |
|---|---|---|---|
| Mitochondrial, VAF ≥ 0 | 2,769 | 0.098 | **0.924** |
| Mitochondrial, VAF ≥ 0.10 | 1,022 | 0.001 | 0.011 |
| BASELINE barcode | 573 | 0.000 | 0.215 |

**Interpretation:** low-frequency mitochondrial variants carry nearly all the
clonal signal; a conventional 10% VAF threshold destroys it (0.924 → 0.011).
This contradicts standard mtscATAC practice and is the strongest candidate
result for Aim 2. **Caveat: this is a crude clustering, not the paper's IQ-TREE
pipeline — a pilot that justifies the experiment, not a result.**

### Pilot script (verbatim — scratchpad has since been cleared)

```r
suppressMessages({library(Matrix); library(data.table)})
d <- "output/physicell_organoid_20260730_214332/lineage_recording"
st <- fread("output/physicell_organoid_20260730_214332/physicell_build/output/organoid_cell_states.csv",
            select=c("ID","founder_ID","cell_type"))
tc <- fread(file.path(d,"terminal_cells.csv.gz"))
tc <- merge(tc, st, by.x="physicell_id", by.y="ID", all.x=TRUE)
set.seed(1)
# restrict to clones with >=5 cells so clonal structure is learnable
big <- tc[, .N, by=founder_ID][N>=5]$founder_ID
sub <- tc[founder_ID %in% big][sample(.N, 2000)]
cat("subsample:", nrow(sub), "cells from", length(unique(sub$founder_ID)), "clones\n")
cat("cell types:\n"); print(table(sub$cell_type))

mt <- readRDS(file.path(d,"mitochondrial_variant_fraction_matrix.rds"))
bc <- readRDS(file.path(d,"barcode_binary_score_matrix_sparse.rds"))
idx <- match(sub$sample_id, rownames(mt)); stopifnot(!any(is.na(idx)))
MT <- mt[idx,,drop=FALSE]; BC <- bc[match(sub$sample_id, rownames(bc)),,drop=FALSE]
truth <- factor(sub$founder_ID)

ari <- function(a,b){ t<-table(a,b); s<-function(x) sum(choose(x,2))
  i<-s(t); ea<-s(rowSums(t)); eb<-s(colSums(t)); n<-s(sum(t))
  (i-ea*eb/n)/((ea+eb)/2-ea*eb/n) }

run <- function(M, thr, lab){
  B <- M; if(thr>0){ B@x[B@x<thr]<-0; B<-drop0(B) }
  keep <- Matrix::colSums(B>0) >= 3; B <- B[,keep,drop=FALSE]
  B@x[] <- 1
  D <- dist(as.matrix(B), method="binary")
  hc <- hclust(D, method="average")
  for (k in c(50, length(unique(truth)))) {
    cl <- cutree(hc, k=k)
    cat(sprintf("%-28s thr=%.2f chars=%5d k=%4d  ARI=%.4f\n", lab, thr, ncol(B), k, ari(truth, cl)))
  }
}
run(MT, 0.00, "mito VAF>=0")
run(MT, 0.10, "mito VAF>=0.10")
run(BC, 0.00, "BASELINE barcode")
```

---

## 6. Environment and feasibility

- R 4.6.0 at `/usr/local/bin/R`. Present: `ape`, `phangorn`, `Matrix`,
  `scDesign3`. **Missing: `TreeDist`, `aricode`.**
- **IQ-TREE is not installed on this machine.** Declared in `heavy_sim_env.yml`
  (env name `babette_2_backup`): `iqtree=3.0.1`, `beast2=2.6.3`,
  `cassiopeia-lineage==2.1.0`. Invoked with `-nt 40`, so reconstruction is a
  **cluster job** — this laptop has 12 cores.
- Conda envs present: `metal`, `miso`, `py39`, `simulation` — none contains IQ-TREE.
- Machine: 12 cores, 32 GB RAM, 87 GB free disk.
- `../PhysiCell` checkout exists and builds.
- 104k leaves is far beyond IQ-TREE; **subsampling to ~2,000–5,000 cells is
  required** and is realistic, since organoid scRNA-seq recovers ~10⁴ cells.

---

## 7. Corrections owed to the manuscript

1. `paper_07_31_2026.md:58` states the PhysiCell export lacks per-cell
   phenotype, substrate exposure, cell-cycle state, and death time. **True for
   `tumor_3D_lineage`, false for the organoid**, which exports `cell_type`,
   `oxygen`, `nutrient`, `developmental_pseudotime`, `founder_ID`,
   `time_in_state`, `transition_count`, plus a full time-resolved fate-transition
   log. Scope the sentence to the tumor project and to "cell-cycle phase and
   death time", which really are missing.
2. The FLARE parameterization claim is unsupported — no FLARE config exists.
3. Methods should state that `RF.dist` treats both trees as **unrooted** while
   the ground truth is rooted and may contain unary nodes.
4. The abstract's closing line is still a placeholder:
   "a NEW GENERATION OF LINEAGE POWER (fix this last line)".

---

## 8. Recommended critical path (from the framework doc)

1. **Tooling first:** metrics module (triplets/quartets, clade precision/recall,
   ARI, F1 purity, **depth-stratified accuracy**); wire reconstruction onto
   PhysiCell outputs; replicate harness + cross-run aggregation; NA masking for
   unrecovered integrations.
2. **Aim 1** (recorder comparison) — cheapest, first real result.
3. **Aim 4** (organoid disease sweep) — long wall-clock, start early and let it run.
4. **Aim 2** (mito vs CRISPR) — needs a mitochondrial observation model in replay.
5. **Aim 3** (inheritance models) — largest new-code + external-data burden;
   **most likely to be cut.** Degrades gracefully to "inheritance model choice
   is a confound for lineage inference" without the real-world comparison.

---

## 9. ⚠️ What has changed since this analysis

Checked 2026-08-27. `git status` shows substantial new work at HEAD `5dd66c4`
that did **not** exist during the audit.

**Modified since the audit** (so all line references above are stale):
`sim5_code.R`, `physicell_lineage.R`, `physicell_mito.R`,
`simulate_physicell_lineage.R`, `physicell_scdesign3.R`,
`generate_physicell_sc_profiles.R`, `README.md`, `PHYSICELL_INTEGRATION.md`,
`FUNCTION_REFERENCE.md`, `BUG_AUDIT.md`, `bash_wrapper_all_combos.sh`.

**New, entirely un-audited modules:**

| File | Size | Note |
|---|---|---|
| `lineage_benchmark.R` | 59 KB | Benchmark **data-generation** harness — synthetic populations, recorder + mito simulation, sampling, observation matrices, conditions. Includes `lineage_benchmark_mito_observation_matrices()` and `lineage_benchmark_mito_sampling_orders()`, which look like they close the **missing mitochondrial observation model** (framework tooling item 5). Does **not** appear to contain tree-accuracy metrics. |
| `run_lineage_benchmark.R` | — | driver for the above |
| `gillespie_lineage.R` | 30 KB | new continuous-time engine |
| `gillespie_pipeline.R`, `simulate_gillespie_lineage.R` | — | drivers |
| `ecdna_lineage.R` | 28 KB | extrachromosomal DNA module |
| `prime_editing.R` | — | standalone — **the prime-editing crash may be fixed** |
| `engine_comparison.R`, `compare_simulation_engines.R` | 14 KB, 6 KB | cross-engine validation |
| `physicell_visium.R`, `analyze_physicell_visium.R`, `aggregate_physicell_visium.R`, `render_physicell_visium_notebook.R`, `run_physicell_visium_replicates.sh` | — | spatial transcriptomics path |
| `run_physicell_organoid_ecdna_pipeline.sh` | — | new organoid variant |
| `example_json_params/palincode_gillespie.json`, `prime_editing_gillespie.json` | — | new presets — "palincode" is a recorder not covered by the audit |
| `notebooks/`, `phylip/` | — | new directories |

**New in `paper/` since the session** — figure and table production has started,
which the framework doc does not account for:

- `paper/figures/figure_1/` — `origin_figure1a.pdf`, `new_figure1a.pdf/.jpg`
  (Figure 1A exists; the framework's figure numbering assumed none did)
- `paper/tables/all_tables.xlsx` (plus an open Excel lock file `~$all_tables.xlsx`)
- `paper/2026.04.16.718941v1.full.pdf` — a bioRxiv preprint added 2026-08-11,
  presumably a reference or a competing/related simulator worth citing

**Before acting on §4, re-verify:**

- [ ] Does `lineage_benchmark.R` supply the mitochondrial observation model
      (recovery / coverage / VAF thresholding)? If so, tooling item 5 is done.
- [ ] Are tree-accuracy metrics beyond normalized RF now available anywhere?
- [ ] Is prime editing fixed (`prime_editing.R`)? The `basepos_erc_nuc_list`
      ordering bug reference no longer matches current line numbers.
- [ ] Does the Gillespie engine change the "two engines" framing in §1 of the
      framework doc? It may now be **three**, which would need a rewrite of that
      section and of the Aim-1/Aim-2 engine assignments.
- [ ] What is "palincode"? It is a recorder technology the Aim-1 panel does not
      include and probably should.
- [ ] Has `BUG_AUDIT.md` absorbed the five bugs listed in §4.5?

The four **aims**, the **draft prose**, the **empirical measurements** in §5, and
the **organoid disease design** in the framework doc remain valid regardless —
they depend on model structure and measured outputs, not on line numbers.
