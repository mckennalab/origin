# Conversation and development handoff — 2026-08-27

## Purpose of this document

This is a dated pickup guide for the extended development conversation that
turned `remote_mito_clean` from a primarily time-step mitochondrial/barcode
simulator into a broader multimodal lineage-simulation framework. It records
the user's requests, the resulting architecture and runnable workflows, the
important modeling decisions, current limitations, verification commands, and
the repository state that a future session should preserve.

This document is a snapshot, not a substitute for the maintained interfaces in
[`README.md`](README.md), [`FUNCTION_REFERENCE.md`](FUNCTION_REFERENCE.md),
[`PHYSICELL_INTEGRATION.md`](PHYSICELL_INTEGRATION.md), and
[`ORGANOID_SIMULATION.md`](ORGANOID_SIMULATION.md). When this file and code
disagree, inspect the current code and regression tests before changing
behavior.

## Executive summary

The repository now has three related population/lineage paths:

1. The legacy native discrete time-step simulator in `sim5_code.R`.
2. An exact, well-mixed continuous-time Gillespie population simulator.
3. A fixed-tree adapter that imports a lineage-enabled PhysiCell division log
   and replays recording systems over exact branch durations.

The event-resolved Gillespie and PhysiCell paths can simulate BASELINE-like
Cas12a base editing, a unified prime-editing recorder, PALINCODE, mitochondrial
lineage tracing, and non-Mendelian ecDNA recording. They emit ground-truth
Newick trees with branch lengths, sparse reconstruction matrices, exact event
logs, literal event-by-descendant ancestry matrices, compressed CSV output,
progress messages, and timing summaries. PhysiCell terminal states can also be
used as covariates for scDesign3 expression simulation.

Two larger demonstration workflows were built:

- a staged 10,000-cell PhysiCell tumor and virtual Visium spatial-lineage
  experiment; and
- a day-21 neural organoid starting from 2,500 iPSCs, with explicit radial glia,
  neural progenitor cells (NPCs), neurons, and a capped epithelial side fate.

The newest addition is a simulation-only lineage-recorder benchmark with
balanced, comb, neutral, turnover, and hierarchical population shapes; cell
samples of 250/1,000/2,000/5,000; recorder integration panels of
1/2/5/10/20; and a separate mitochondrial observation-depth panel.

## Chronological request log

The following list preserves the sequence and intent of the user's requests.
It is useful context for why some interfaces and output aliases exist.

1. Search the codebase for bugs and document undocumented functions.
2. Integrate the repository with an updated `../PhysiCell/` that tracks and
   exports cell lineage/division data, then simulate lineage recording on the
   imported tree.
3. Create an end-to-end script that grows a 10,000-cell PhysiCell cancer model
   and runs mitochondrial and base-editing lineage simulation.
4. Explore and then implement standard covariate-linked scDesign3 expression
   simulation.
5. Update the paper markdown to reflect the simulator, PhysiCell, and
   transcriptomics extensions.
6. Repeat the bug search, fix a specifically discussed bug number 2, report
   remaining defects, and work through six high-priority bugs.
7. Identify prerequisites missing from a `sim5_code.R` CLR run.
8. Design a neuronal-organoid simulation starting from 2,500 iPSCs and begin
   implementing it.
9. Explain longitudinal organoid reference-data format, running PhysiCell to
   day 90, and expected day-90 scale.
10. Reassess whether the early organoid cell count was biologically plausible,
    then focus on a realistic day-21 expansion with higher proliferation and
    turnover.
11. Tune the model toward approximately 25% neurons.
12. Add explicit NPCs that generate neurons and a limited epithelial branch
    originating from iPSCs.
13. Wrap the organoid PhysiCell tree with a BASELINE-like
    hyperdCas12a-ABE8e system: five integrations, 50 targets per integration,
    dynamic recording, plus mitochondrial lineage tracing.
14. Add progress reporting for large lineage simulations.
15. Diagnose slow lineage replay and make validation/tree reconstruction
    linear using hash-map parent lookups and preallocated node columns.
16. Add R stage logs announcing the end of lineage simulation and the start and
    completion of every output phase.
17. Diagnose slow Newick writing and implement a linear, iterative, buffered
    writer.
18. Diagnose slow BASELINE-like steps and implement the five discussed
    performance fixes, including compact sparse state handling and output
    improvements.
19. Add a total R runtime summary at the end of a run.
20. Compress mutation-event and other tabular output as `.csv.gz`.
21. Clarify how imported PhysiCell branch durations produce continuous-time
    mutations versus the legacy simulator's discrete time grid.
22. Confirm that neurons are post-mitotic in the current organoid model.
23. Determine whether a cell-by-event reconstruction matrix exists, then add a
    literal matrix with one column per mutation-event row and `1` for every
    sampled topological descendant.
24. Add a Gillespie framework alongside the existing time-step engine and
    execute the implementation plan.
25. Add a quick paired script comparing continuous-time and time-step
    population outcomes.
26. Add non-Mendelian ecDNA barcode simulation with partial labeling, an
    immutable static ID, a limited CRISPR recorder, copy propagation, and
    random daughter segregation.
27. Add selection coefficients against mitochondrial variant burden, ecDNA
    labels, and ecDNA recorder edits.
28. Add PALINCODE support from `paper/2026.04.16.718941v1.full.pdf`, including
    left, right, and rare simultaneous left-and-right outcomes.
29. Confirm Newick output and branch-length behavior across simulation paths.
30. Assess prime-editing support, then create one known-pegRNA backend shared
    by the time-step, Gillespie, and PhysiCell workflows, with per-pegRNA
    efficiency.
31. Design a realistic benchmark over population sizes and tree shapes.
32. Restrict benchmark sizes to 250, 1,000, 2,000, and 5,000 cells and all
    integrated systems to 1, 2, 5, 10, and 20 integrations. Retain six targets
    per prime-editing integration and two targets per PALINCODE integration.
33. Create a simulation-only benchmark script that writes ground truth and
    recording matrices but leaves tree reconstruction for a separate stage.
34. Explain how to run `run_lineage_benchmark.R` and describe its output
    directory.
35. Design a spatial tumor analysis that takes planar sections and asks how
    ground-truth lineage distance changes with physical distance.
36. Add replicated, realistic conventional-Visium sampling and an R notebook
    plotting physical distance against lineage distance.
37. Add mitochondrial mutations to `run_lineage_benchmark.R` as a first-class
    benchmark system.
38. Create this dated conversation handoff.

## Current architecture

### Native discrete time-step simulation

- Entry point: `Rscript sim5_code.R -P parameters.json`.
- Population changes and recording are evaluated on the configured time grid.
- Existing JSON files continue to select this path unless
  `simulation_engine` is set to `gillespie`.
- It retains the original mitochondrial organelle/genome dynamics and sparse
  mutation representation.
- The top-level multi-parameter wrapper is `bash_wrapper_all_combos.sh`.

### Exact Gillespie simulation

- Entry points: `simulate_gillespie_lineage.R` or automatic dispatch from
  `sim5_code.R` when `simulation_engine: "gillespie"`.
- Core code: `gillespie_lineage.R` and `gillespie_pipeline.R`.
- Division hazard is `1 / cell_cycle_length`.
- A death probability `p` per cell cycle becomes
  `-log(1-p) / cell_cycle_length`.
- Cell-type transitions occur at division using the configured transition
  matrix. Editing and differentiation induction occur at their exact times.
- Barcode, mitochondrial, and ecDNA overlays are conditionally simulated on
  the resulting tree. This is equivalent to a joint Gillespie system only
  while recorder state does not feed back into birth/death propensities.
- The current Gillespie implementation rejects legacy
  `consider_cell_heteroplasmy_scores: true` because that would violate this
  separation.

### PhysiCell fixed-tree replay

- Entry point: `simulate_physicell_lineage.R`.
- Core code: `physicell_lineage.R`, `physicell_mito.R`, and `ecdna_lineage.R`.
- Required division input:

  ```csv
  time,parent_ID,daughter_ID
  60,0,12
  90,0,19
  120,12,31
  ```

- Events can be unsorted. IDs must be non-negative integers and each daughter
  must be born once.
- The updated PhysiCell convention retains the parent's ID at division and
  gives a new ID only to the new daughter. The adapter converts every event
  into a continuing-parent branch and a new-daughter branch.
- The optional final current-cell table needs an `ID` column. It is important
  because division logs alone do not report cell removal.
- Each recording hazard is integrated over exact branch duration, so recording
  events have continuous times even though the cellular tree came from a
  PhysiCell simulation with its own numerical time step.

## Major implemented capabilities

### Lineage import, validation, and Newick

- Division validation and tree construction are linear in event count.
- Active persistent IDs are resolved with hashed cell-to-node lookups.
- Node columns are preallocated and grown geometrically when necessary.
- Newick writing uses integer parent/child arrays, an explicit stack, reverse
  ancestor marking for sampled trees, and buffered output.
- Full and sampled trees preserve elapsed-time branch lengths. Sampled trees
  retain unary ancestors so pruning does not silently change elapsed time.

### BASELINE-like Cas12a recording

- HyperdCas12a-ABE8e-like irreversible A-to-G editing.
- The organoid preset uses five integrations and 50 primary targets per
  integration, with 20 high-rate and 30 low-rate targets.
- Seven-base decaying editing windows are supported.
- The working organoid preset uses a mean induced probability of `0.00952` per
  target per 24-hour iPSC cycle. This is a simulation assumption, not a fitted
  experimental estimate.
- Sparse profiles, compact output, released internal profiles, buffered event
  accumulation, and sparse-triplet preallocation are used to keep large runs
  tractable.

### Unified prime editing

- Core code: `prime_editing.R`.
- One backend is used by the native, Gillespie, and PhysiCell paths.
- A known pegRNA pool requires `pegRNA_id`, `edit_sequence`, and
  `editing_efficiency`; optional fields include spacer, PBS, RTT, and
  description.
- Pools may be embedded in JSON or loaded from CSV. An example is
  `data/example_pegRNA_pool.csv`.
- Each integration has an immutable static ID; targets can be assigned pegRNAs
  exactly or by cycling/sampling.
- Effective probability is `1 - (1 - p)^e`, where `p` is base editing
  probability and `e` is pegRNA efficiency.
- Benchmark prime integrations contain six targets.

### PALINCODE

- Supported in event-resolved Gillespie and PhysiCell paths.
- Each cBit locks at its first event to `0=WT`, `1=left`, `2=right`, or
  `3=both`; resolved sites do not edit again.
- Left/right/both fractions are configurable per cBit and must sum to one.
- The benchmark retains two cBits per integration.
- The legacy native time-step nucleotide engine is not currently a PALINCODE
  implementation.

### Mitochondrial lineage tracing

- Mitochondrial genomes are sparse lists of variants within each cell.
- At true cell division, each child independently bottleneck-samples a fixed
  number of genomes from the parent, allowing heteroplasmy drift and loss.
- Exact mutation times and genome IDs are recorded.
- A selection coefficient can down-weight genomes according to variant burden,
  but the current replay fixes total genomes per cell and therefore does not
  model selection on mitochondrial copy number or cellular fitness.

### Non-Mendelian ecDNA

- Core code: `ecdna_lineage.R`.
- A configurable fraction or number of founder ecDNA species is labeled.
- Labeled species have an immutable static ID and limited irreversible CRISPR
  target array.
- Copies replicate and the joint pool partitions binomially between daughters,
  so inheritance is not Mendelian and species can be lost.
- Selection coefficients can penalize labeled copies and recorder-edited
  copies during intracellular propagation.
- The ecDNA reconstruction character matrix is based on actual terminal copy
  inheritance, not all topological descendants of the branch on which an edit
  arose.

### Literal mutation-event ancestry matrix

- `mutation_event_descendant_matrix_sparse.rds` has one column for every
  mutation-event table row.
- A sampled terminal cell gets `1` when it is a topological descendant of that
  event's branch.
- `mutation_event_descendant_manifest.csv.gz` maps every column to modality,
  original event row, branch, event fields, and sampled-descendant count.
- For mitochondrial events, this is ground-truth branch ancestry and does not
  claim that the particular mutated genome survived every bottleneck.

### Progress, output compression, and timing

- Large reconstruction/recording phases report processed work, percentage,
  elapsed time, throughput, and ETA.
- Timestamped stage logs announce the transition from simulation to output and
  each subsequent artifact group.
- R-generated tabular outputs are gzip-compressed by default; readers accept
  `.csv` and `.csv.gz` and prefer compressed files.
- `r_timing_summary.csv.gz` and the terminal log report each major phase and
  total R simulator wall time.

## PhysiCell workflows

### 10,000-cell tumor

The standard workflow expects the lineage-enabled PhysiCell checkout at
`../PhysiCell/` and does not modify it. It stages a copy of
`user_projects/tumor_3D_lineage` under the selected output directory.

```bash
bash run_physicell_10000_pipeline.sh
```

Useful explicit form:

```bash
bash run_physicell_10000_pipeline.sh \
  --target-cells 10000 \
  --num-integrations 1 \
  --mt-genomes-per-cell 8 \
  --seed 1 \
  --jobs 8 \
  --output-dir output/my_physicell_run
```

Use `--modalities lineage` when only the event-resolved truth tree is needed.
Use `--sc-reference` and related options to append scDesign3 expression.

### Neural organoid

The current day-21 model starts from 2,500 individually tracked iPSCs in a
150-micron-radius spheroid. Its seven states are iPSC, neuroepithelial, radial
glia, NPC, neuron, astrocyte, and epithelial. Neurons and epithelial cells are
terminal/non-dividing in the current model. Astrocytes are gated until day 45
and are therefore not expected in the day-21 preset.

```bash
bash run_physicell_organoid_pipeline.sh
```

The calibrated `early21` seed-2 run produced 104,804 cells at day 21:

- 25,986 neurons (24.79%);
- 60,213 NPCs (57.45%);
- 18,499 radial glia (17.65%);
- 99 surviving epithelial cells (0.094%);
- 7 residual iPSC/neuroepithelial cells.

The epithelial branch uses a default 2% commitment probability and a cap of
100 commitments. The output is an order-of-magnitude calibration result, not a
general biological prediction. The conservative `long90` preset remains
available, but discussion shifted toward day 21 because day-90 cell number and
composition require substantially more empirical calibration and compute.

BASELINE plus mitochondrial replay:

```bash
bash run_physicell_organoid_baseline_pipeline.sh \
  --output-dir output/organoid_baseline_day21
```

BASELINE, mitochondrial, and ecDNA replay:

```bash
bash run_physicell_organoid_ecdna_pipeline.sh \
  --output-dir output/organoid_baseline_ecdna_day21
```

Small installation/interface test:

```bash
bash run_physicell_organoid_pipeline.sh \
  --initial-cells 8 \
  --initial-radius 40 \
  --target-cells 100 \
  --days 3 \
  --founder-label-sites 3 \
  --mt-genomes-per-cell 4 \
  --output-dir output/organoid_smoke_test
```

The organoid PhysiCell project is under
`physicell_projects/neural_organoid_lineage/`.

## scDesign3 integration and longitudinal data contract

The standard covariate-linked implementation fits scDesign3 to a reference and
simulates new counts at one row of covariates per simulated terminal cell.

Accepted reference containers:

- `SingleCellExperiment` RDS;
- Seurat RDS; or
- AnnData `.h5ad`.

The reference should contain raw counts, a cell-type label, and—when a
developmental trajectory is desired—a numeric pseudotime column. Longitudinal
organoid references should retain a culture-day/time field as well. Any field
requested through `--sc-other-covariates` must exist with the same name in the
reference metadata. A JSON mapping can reconcile simulator fate names with
reference labels.

Typical organoid call:

```bash
bash run_physicell_organoid_pipeline.sh \
  --sc-reference data/cortical_organoid_reference.rds \
  --sc-celltype-col cell_type \
  --sc-pseudotime-col pseudotime \
  --sc-other-covariates culture_day,oxygen \
  --sc-celltype-map config/organoid_celltype_map.json
```

The simulated metadata can include cell state, developmental pseudotime,
culture day, time in state, transition count, x/y/z position, radial position,
oxygen, nutrient, neighbor count, barcode edit burden, mitochondrial
heteroplasmy burden, and ecDNA burden. The organoid pseudotime comes from fate
and time in state, not merely lineage depth.

Current limitation: scDesign3 conditions on the supplied covariates but does
not add a latent sibling/clone expression correlation after conditioning.

## Spatial tumor/Visium workflow

The spatial workflow tests whether physically close cells are more closely
related on average. Proximity is expected to enrich relatedness because of
local expansion, but close cells are not guaranteed to be close relatives;
mixing, mechanics, turnover, and independent clones can break that relation.

Run replicated tumors and conventional-Visium sections:

```bash
bash run_physicell_visium_replicates.sh \
  --replicates 5 \
  --target-cells 10000 \
  --slice-offsets=-50,0,50 \
  --section-thickness 5 \
  --spot-diameter 55 \
  --spot-pitch 100 \
  --jobs 8
```

The default conventional grid has 4,992 spots in 78 staggered rows, with
64 spots per row, 55-micron spot diameter, and 100-micron center spacing.
Metrics include MRCA age, patristic elapsed-time distance, lineage-edge
distance, and shared-founder status. Shuffled tip locations provide a spatial
null. Tumors, not individual cell pairs, are the replicate unit in aggregate
analysis.

Render the parameterized notebook with:

```bash
Rscript render_physicell_visium_notebook.R \
  --input-path output/physicell_visium_YYYYMMDD_HHMMSS
```

Notebook source:
`notebooks/physicell_visium_lineage_distance.Rmd`.

## Paired recorder benchmark

### Design

`run_lineage_benchmark.R` is simulation-only. It deliberately does not run
tree reconstruction.

Default cell samples:

- 250;
- 1,000;
- 2,000; and
- 5,000 cells.

Default integrated-recorder panel:

- 1, 2, 5, 10, and 20 integrations;
- BASELINE: 50 targets/integration;
- prime editing: 6 targets/integration; and
- PALINCODE: 2 cBits/integration.

Population shapes are balanced, comb, neutral asynchronous, and hierarchical;
turnover is optional. There are ten default seeds.

Mitochondrial recording is a separate observation-depth experiment rather
than an “integration” experiment:

- 32 modeled biological genomes per cell by default;
- nested without-replacement observations of 1, 2, 5, 10, and 20 genomes per
  cell;
- one shared observed-variant universe across depths for a given sampled tree;
- heteroplasmy state matrix and binary presence matrix;
- exact mitochondrial event table and ground-truth tree; and
- `integrations=NA` plus an explicit `observation_depth` in manifests.

The mitochondrial benchmark preset uses a 16,569-base reference, no founder
heteroplasmy or indels, and a Jukes-Cantor substitution probability of `2e-6`
per alternate base per cell cycle.

### Running it

Complete default benchmark:

```bash
Rscript run_lineage_benchmark.R \
  --output-dir output/lineage_benchmark
```

Small pilot:

```bash
Rscript run_lineage_benchmark.R \
  --output-dir output/lineage_benchmark_pilot \
  --shapes balanced,neutral \
  --tree-sizes 250,1000 \
  --integration-counts 1,5,20 \
  --mt-observation-depths 1,5,20 \
  --mt-genomes-per-cell 32 \
  --seeds 1:2
```

Mitochondrial-only example:

```bash
Rscript run_lineage_benchmark.R \
  --output-dir output/mitochondrial_benchmark \
  --systems mitochondrial \
  --tree-sizes 250,1000,2000,5000 \
  --mt-observation-depths 1,2,5,10,20 \
  --mt-genomes-per-cell 32 \
  --seeds 1:10
```

The default grid performs 40 population simulations and 160 full recorder
overlays, producing 3,200 reconstruction-input conditions. Runs resume only
when `benchmark_settings.rds` exactly matches. Benchmark version 2 added
mitochondrial support; use a new output directory rather than attempting to
resume a version-1 benchmark.

### Condition outputs

Every condition contains:

- `ground_truth_tree.nwk`;
- `recording_logical_target_matrix_sparse.rds`;
- `recording_state_matrix_sparse.rds`;
- `recording_character_matrix_sparse.rds`;
- `target_manifest.csv.gz`;
- `sample_cells.csv.gz`; and
- `condition_manifest.csv.gz`.

Mitochondrial conditions also contain
`mitochondrial_variant_fraction_matrix_sparse.rds`,
`mitochondrial_binary_variant_matrix_sparse.rds`, and
`mitochondrial_variant_manifest.csv.gz`. The full mitochondrial directory
contains exact mutation events, terminal genome profiles, reference FASTA,
the observation matrices, and reproducible genome sampling orders.

## Important output semantics

- `recording_state_matrix_sparse.rds` contains the native state. For
  mitochondrial benchmark conditions this is heteroplasmy among sampled
  genomes.
- `recording_character_matrix_sparse.rds` contains reconstruction characters.
  For mitochondrial benchmark conditions this is binary variant presence.
- `recording_logical_target_matrix_sparse.rds` preserves one physical logical
  target for cross-system comparisons. BASELINE windows become integer
  bitmasks; prime is binary; PALINCODE is categorical 0–3; mitochondrial is
  binary variant presence.
- Exact event logs are not equivalent to observed characters for
  non-Mendelian systems because mitochondrial/ecDNA copies can be lost.
- All matrices should be joined by row name (`cell_<PhysiCell_ID>`), never by
  implicit row order alone.

## Bug-fix and performance work

[`BUG_AUDIT.md`](BUG_AUDIT.md) is the authoritative detailed defect table. It
records fixes across startup portability, function-call wiring, parallel
workers and RNG reproducibility, barcode construction and target allocation,
editing-window logic, substitution models, heterogeneous mutation sampling,
prime editing, mitochondrial initialization/inheritance, molecular recovery,
score matrices, result aggregation, tree comparison/plotting, shell wrappers,
IGV conversion, scDesign3, and fractional induction.

Important performance outcomes from the conversation:

- linear PhysiCell validation and event-tree construction;
- linear iterative Newick rendering with bounded buffering;
- sparse recorder profiles and sparse matrix assembly;
- release of internal node profiles after their final child is processed;
- compact output for large BASELINE runs;
- compressed event/output tables;
- progress/ETA reporting around expensive phases; and
- phase plus total R timing summaries.

## Known limitations and modeling boundaries

These are especially important when interpreting results or planning the next
implementation session.

- Native insertion alleles still use floating-point decimal encoding; long or
  interacting insertions need a structured representation.
- The legacy `max_mito_per_cell` saturation behavior remains undefined and is
  not enforced.
- The old heatmap helper expects an obsolete schema and
  `fit_plot_parameters.R` depends on an untracked calibration CSV.
- Organoid developmental rates are working assumptions, not fitted biology.
- The organoid lacks rosettes, lumen/apical polarity, radial fibers, cortical
  layers, electrophysiology, and region-specific morphogen patterning.
- PhysiCell fate transitions are logged, but imported recording branches are
  not yet split at every state-transition time; recorder kinetics currently
  use one configured rate block per imported branch.
- PhysiCell division logs do not contain exact death times; dead cells are
  excluded from final observation with the current-cell table.
- The organoid founders do not share a simulated pre-aggregation expansion
  tree.
- BASELINE simulation does not yet include integration dropout, sequencing
  errors, amplicon capture, or target-specific empirical calibration.
- PALINCODE does not expand left/right outcomes into every possible
  within-window nucleotide combination.
- Mitochondrial copy number is fixed in event-resolved replay.
- ecDNA and mitochondrial selection act on intracellular propagation on an
  already generated cellular tree; they do not currently change cell fitness
  or the tree itself.
- The recorder benchmark does not reconstruct trees. A separate reconstruction
  workflow still needs to consume its manifest and report accuracy/runtime.
- Full-scale 5,000-cell, all-system benchmark execution was not part of the
  mitochondrial implementation smoke tests.

## Paper and narrative artifacts

- `paper/paper_2026_07_27_markdown.md` is the broad manuscript draft covering
  native simulation, fixed-tree replay, mitochondrial recording, PhysiCell,
  scDesign3, and spatial analysis.
- `paper/paper_07_31_2026.md` is a shorter ORIGIN-oriented draft.
- `paper/paper_07_31_2026_framework.md` is a detailed aims/results plan with
  capability and threat-to-validity sections.
- `paper/2026.04.16.718941v1.full.pdf` is the local PALINCODE source paper.

Before editing the paper again, reconcile claims against completed runs. In
particular, distinguish measured smoke/calibration results from planned
benchmark results and from biological expectations.

## Key file map

| Area | Primary files |
| --- | --- |
| Legacy time-step simulator | `sim5_code.R`, `nonuniform_muts_heterogeneous.R`, `substitution_models.r` |
| PhysiCell tree/replay | `simulate_physicell_lineage.R`, `physicell_lineage.R`, `physicell_mito.R` |
| Exact Gillespie | `simulate_gillespie_lineage.R`, `gillespie_lineage.R`, `gillespie_pipeline.R` |
| Prime editing | `prime_editing.R`, `data/example_pegRNA_pool.csv` |
| ecDNA | `ecdna_lineage.R` |
| Benchmark | `run_lineage_benchmark.R`, `lineage_benchmark.R` |
| scDesign3 | `scdesign3_helpers.R`, `physicell_scdesign3.R`, `generate_physicell_sc_profiles.R` |
| 10k PhysiCell tumor | `run_physicell_10000_pipeline.sh`, `example_json_params/physicell_10000.json` |
| Neural organoid | `run_physicell_organoid_pipeline.sh`, `run_physicell_organoid_baseline_pipeline.sh`, `run_physicell_organoid_ecdna_pipeline.sh`, `physicell_projects/neural_organoid_lineage/` |
| Spatial/Visium | `run_physicell_visium_replicates.sh`, `physicell_visium.R`, `analyze_physicell_visium.R`, `aggregate_physicell_visium.R`, `notebooks/physicell_visium_lineage_distance.Rmd` |
| Engine comparison | `compare_simulation_engines.R`, `engine_comparison.R` |
| Tests | `tests/regression_tests.R`, `tests/fixtures/` |
| Maintained docs | `README.md`, `FUNCTION_REFERENCE.md`, `PHYSICELL_INTEGRATION.md`, `ORGANOID_SIMULATION.md`, `BUG_AUDIT.md` |

## Repository state at handoff

The worktree is intentionally dirty and contains substantial user/session work.
Do not reset, clean, or overwrite it wholesale.

At the time this handoff was written:

- tracked files including the main documentation, `sim5_code.R`, PhysiCell
  adapters, scDesign3 code, and shell utilities had modifications;
- many major additions—including the Gillespie, benchmark, ecDNA, Visium,
  organoid, paper, and regression-test files—were still untracked by Git; and
- `output/`, `.DS_Store`, notebook checkpoints, and other generated artifacts
  were also present and untracked.

Before committing, deliberately separate source/documentation from large
generated outputs and incidental files. Do not infer that an untracked file is
disposable.

## Verification and reproducibility

Primary dependency-light regression command:

```bash
Rscript tests/regression_tests.R
```

This command passed on 2026-08-27 after the handoff document was created.

Useful additional checks:

```bash
Rscript -e "invisible(parse(file='lineage_benchmark.R')); invisible(parse(file='run_lineage_benchmark.R'))"
bash -n bash_wrapper_all_combos.sh
bash -n run_physicell_10000_pipeline.sh
bash -n run_physicell_organoid_pipeline.sh
bash -n run_physicell_organoid_baseline_pipeline.sh
bash -n run_physicell_organoid_ecdna_pipeline.sh
bash -n run_physicell_visium_replicates.sh
git diff --check
```

The bug-audit session also reports a compiled eight-cell PhysiCell smoke run
and a live scDesign3 1.10 conditional simulation. The 10,000-cell workflow,
full all-system benchmark, BEAST, MUSCLE, IQ-TREE, BWA, and samtools were not
all rerun as part of every subsequent change.

## Recommended pickup sequence

1. Read this file and the current `git status --short`; preserve unrelated
   worktree changes.
2. Run `Rscript tests/regression_tests.R` before making new edits.
3. Decide which population path is in scope: native time-step, well-mixed
   Gillespie, or spatial PhysiCell.
4. For PhysiCell work, verify that `../PhysiCell/` is the lineage-enabled
   checkout and confirm the division/current-cell CSV headers.
5. For organoid/scDesign3 work, select the actual longitudinal reference and
   calibrate state labels, pseudotime, culture day, and requested covariates.
6. For benchmark work, use a new version-2 output directory, start with a
   small pilot, inspect `benchmark_manifest.csv.gz`, and add reconstruction as
   a separate resumable stage.
7. Treat output matrices according to their modality semantics, especially
   the distinction between topological event ancestry and observed
   non-Mendelian allele retention.
8. Update the maintained documentation and regression tests with every public
   interface change.

## High-value next steps

- Build the separate tree-reconstruction driver for the benchmark manifest,
  with per-method timeouts and accuracy/runtime aggregation.
- Calibrate organoid growth and fate-transition parameters against a selected
  longitudinal dataset rather than a literature-scale anchor alone.
- Split imported PhysiCell branches at state-transition times so recorder and
  mitochondrial kinetics can change immediately with fate.
- Add explicit observation noise: integration dropout, capture efficiency,
  sequencing depth/error, allelic dropout, and mitochondrial sampling noise.
- Decide whether marker selection should feed back into cellular birth/death
  propensities; if so, simulate population and recorder state jointly.
- Replace floating-point insertion encoding with structured alleles.
- Run and archive full-scale performance profiles for 5,000-cell benchmark
  conditions and the complete day-21 organoid recording stack.
