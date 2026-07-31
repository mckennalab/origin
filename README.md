# remote_mito_clean

This document describes the repository’s simulation and analysis workflow.
See [`FUNCTION_REFERENCE.md`](FUNCTION_REFERENCE.md) for the complete contract
of every named function, including internal helpers. See
[`BUG_AUDIT.md`](BUG_AUDIT.md) for the repository-wide defect audit, fixes, and
remaining model-level limitations. See
[`PHYSICELL_INTEGRATION.md`](PHYSICELL_INTEGRATION.md) for growing a
lineage-enabled PhysiCell tumor and simulating barcode plus mitochondrial
recording data on its branches. See
[`ORGANOID_SIMULATION.md`](ORGANOID_SIMULATION.md) for the 2,500-founder
neural-organoid workflow.

## 1. Purpose

`remote_mito_clean` is an R-based forward-time simulation of a growing cell
population in which each cell carries:

1. a **mitochondrial (mt) genome population** — many copies per cell, undergoing
   inheritance, fusion/fission, drift, and heteroplasmy, and
2. one or more **integrated lineage barcodes (bc)** — synthetic loci that
   accumulate mutations under either base-editor (BE) or nuclease editing
   regimes.

The simulation runs an arbitrary number of cell types, each with its own cell
cycle, death rate, substitution model, target/non-target mutation rates, and
transition matrices that govern fate decisions before and after editing
induction. At chosen stopping points, the cell population is downsampled and
exported in two forms:

- nucleotide **FASTA** sequences for each cell (per integration / mt genome
  sample), suitable for IQ-TREE phylogenetic reconstruction, and
- character-state / binary **score matrices** (RDS, FASTA, PHYLIP) summarising
  per-target mutational profiles for tree reconstruction by score-based methods.

The same wrapper script then aligns sequences (`muscle`), reconstructs trees
(`iqtree`), compares them to the simulator's ground-truth lineage tree
(Robinson–Foulds distance via `phangorn`/`ape`), and aggregates results across
JSON parameter sweeps.

## 2. Pipeline at a glance

```
JSON params dir
        │
        ▼
bash_wrapper_all_combos.sh             ┐
        │                              │
        ▼                              │
Rscript sim5_code.R -P <params.json>   │  per-JSON loop
        │                              │
        ├─► output/processed_fastas/   │
        ├─► output/score_mats/         │
        ├─► output/processed_newicks/  │  (ground-truth trees)
        ├─► output/cell_populations/   │
        └─► output/run_logs|run_specs/ │
                                       │
muscle  -align *_TERM.fasta            │
iqtree  -s *_msa.fasta -m HKY ...      │  → output/recon_trees/<run_id>/fasta_*
iqtree  -s *.fasta     -st BIN -m MF   │  → output/recon_trees/<run_id>/score_*
                                       │
Rscript compare_trees_call_from_bash.r │  → output/tree_images/, rf_dist_files/
Rscript process_results_from_bash.r    │  → output/param_results_files/
Rscript convert_scoremats_to_csvs.r    │  → output/score_mats/<run_id>/csvs/
                                       ┘

parse_rf_results.ipynb (optional, post-hoc)  → heatmaps
```

## 3. Environment

The conda environment is defined by `heavy_sim_env.yml` (env name in the file:
`babette_2_backup`). It pulls a large stack from `conda-forge` / `bioconda`,
including:

- R (with `ape`, `phangorn`, `Matrix`, `data.table`, `rjson`, `optparse`,
  `seqinr`, `Biostrings`, `msa`, `parallel`, `ggplot2`, `ggmuller`, `babette`,
  etc. — see the `library()` block at the top of `sim5_code.R`),
- BEAST 2 + `beagle-lib` (for `babette`-driven Bayesian inference; currently the
  `make_babette_tree.r` source line is commented out in `sim5_code.R`),
- `muscle`, `iqtree`, `samtools`, `bwa` for the bash wrapper and IGV utility.

Activate the environment before running anything:

```bash
conda env create -f heavy_sim_env.yml
conda activate babette_2_backup
```

## 4. Running a simulation

The single user-facing command is the bash wrapper:

```bash
bash bash_wrapper_all_combos.sh                   # interactive: prompts for params dir
bash bash_wrapper_all_combos.sh -d <params_dir>   # non-interactive

# With scDesign3 single-cell profile generation enabled:
bash bash_wrapper_all_combos.sh -d <params_dir> \
     -r <reference.rds|h5ad> \
     [--celltype_col <name>] [--celltype_map <map.json>] \
     [--use_pseudotime] [--sc_ncores 4] [--sc_max_cells 5000]
```

`<params_dir>` must be a directory containing one or more JSON parameter files.
For each file in the directory, the wrapper:

1. Invokes `Rscript sim5_code.R -P <file>`. The simulator generates a
   numeric `unique_run_id` and writes its outputs under several `output/<...>/<run_id>/`
   subdirectories (see §6).
2. Reads the completed run id from the simulator’s `UNIQUE_RUN_ID=<id>` output,
   avoiding ambiguity when multiple runs share an output directory.
3. For every `*_TERM.fasta` (terminal-timepoint sequences, excluding any
   pre-existing `*_msa.fasta`), runs `muscle -align` and then
   `iqtree -m HKY -bb 1000 -alrt 1000 -nstep 80 -nt 40` to produce a
   nucleotide-based tree, writing into
   `output/recon_trees/<run_id>/fasta_<...>/`.
4. If `output/score_mats/<run_id>/phylips/` exists, runs `iqtree -st BIN -m MF`
   on each score-matrix file to produce a binary character-state tree, writing
   into `output/recon_trees/<run_id>/score_<...>/`.
   - Note: the wrapper currently expects the score-matrix files to be saved as
     `.fasta` inside the `phylips/` directory (a comment in the script flags
     this as a "temporary solution"); the `.phy`-handling lines are still
     present but commented out.
5. For each reconstructed treefile, calls
   `Rscript compare_trees_call_from_bash.r -R <treefile> -I <run_id> -S <savename> -P <json> -T <treedir>`.
6. Calls `Rscript process_results_from_bash.r -I <run_id>` to aggregate RF
   distance files into a stacked results CSV.
7. Calls `Rscript convert_scoremats_to_csvs.r --urid <run_id> --score_mat_path ...`
   to dump score matrices as human-readable CSVs.
8. **(Optional)** If `-r <reference>` was supplied, calls
   `Rscript generate_sc_profiles_from_bash.r -I <run_id> -R <reference> ...`
   to fit scDesign3 on the supplied reference scRNA-seq dataset and synthesise
   per-cell expression profiles for every living terminal cell at every stopping point
   (see §9). The fit is cached under `output/scdesign3_fits/` keyed by reference
   path + mtime + size, so subsequent runs against the same reference reuse it.

Each iteration prints the per-file wall-time via the bash `SECONDS` builtin.

`sim5_code.R` resolves its project directory from its own script path, so the
wrapper can be launched from another working directory.

### 4.1 PhysiCell lineage replay

With the lineage-enabled PhysiCell checkout at `../PhysiCell`, the complete
10,000-cell workflow is:

```bash
bash run_physicell_10000_pipeline.sh
```

The runner stages and compiles
`user_projects/tumor_3D_lineage`, grows the tumor to 10,000 current cells,
imports `cell_lineage.csv` and `lineage_table.csv`, and writes base-editing
barcode plus mitochondrial observations under a timestamped
`output/physicell_10000_*` directory. It does not modify the PhysiCell
checkout. The default recording parameters are in
`example_json_params/physicell_10000.json`.

Supplying a reference dataset adds covariate-linked scDesign3 expression:

```bash
bash run_physicell_10000_pipeline.sh \
  --sc-reference data/tumor_reference.rds \
  --sc-celltype-col cell_type \
  --sc-use-pseudotime
```

The optional expression stage uses the same `cell_<PhysiCell_ID>` sample names
as the barcode, mitochondrial, and tree outputs. Lineage depth is mapped to the
reference pseudotime range. `--sc-spatial-cols ref_x,ref_y` additionally maps
PhysiCell x/y positions to reference spatial coordinates.

To replay an existing lineage without running PhysiCell:

```bash
Rscript simulate_physicell_lineage.R \
  --lineage path/to/cell_lineage.csv \
  --live-cells path/to/lineage_table.csv \
  --params example_json_params/physicell_10000.json \
  --end-time 1440 \
  --modalities both \
  --output-dir output/physicell/my_run
```

The adapter converts PhysiCell's retained-parent-ID division events into a
binary event tree and evolves heritable barcode profiles and sparse
mitochondrial haplotypes over its exact branch durations. The current-cell
table is optional but recommended because division events alone do not
identify cells removed before the final snapshot. Full input contracts,
outputs, assumptions, and a fixture are in
[`PHYSICELL_INTEGRATION.md`](PHYSICELL_INTEGRATION.md).

### 4.2 Neural organoid development

The cortical-organoid MVP starts from 2,500 individually tracked iPSCs. Its
default `early21` preset uses a 24-hour iPSC cycle, 8-day radial-glia and NPC
cycles, and a reduced NPC-to-neuron transition hazard to target approximately
25% neurons at day 21. Radial glia produce explicit neural progenitor cells,
which produce neurons. A separate direct iPSC-to-epithelial branch defaults to
a 2% commitment probability and a 100-cell cap. The model simulates seven
developmental states in a spatial oxygen/nutrient environment, then adds
founder-resolving barcodes,
mitochondrial recording, and optional state-linked scDesign3 counts:

```bash
bash run_physicell_organoid_pipeline.sh
```

Use `--epithelial-probability` and `--max-epithelial-cells` to tune or disable
the limited epithelial side population.

The full seed-2 calibration produced 104,804 cells at day 21: 25,986 neurons
(24.79%), 60,213 NPCs, 18,499 radial glia, and 99 surviving epithelial cells
from 100 epithelial commitments.

To add a BASELINE-like Cas12a recorder and mitochondrial lineage tracing to
the same PhysiCell tree, use:

```bash
bash run_physicell_organoid_baseline_pipeline.sh \
  --output-dir output/organoid_baseline_day21
```

This preset creates five independently inherited recorder integrations with
50 primary targets each, adds founder-resolving positions, and emits aligned
sparse Cas12a, mitochondrial, and combined lineage-feature matrices. See
[`ORGANOID_SIMULATION.md`](ORGANOID_SIMULATION.md) for the model assumptions
and output contract.

With a longitudinal organoid reference:

```bash
bash run_physicell_organoid_pipeline.sh \
  --sc-reference data/cortical_organoid_reference.rds \
  --sc-celltype-col cell_type \
  --sc-pseudotime-col pseudotime
```

The staged model preserves all day-zero founders, records cell-state
transitions, and supplies final cell type, developmental pseudotime, culture
day, spatial position, oxygen, nutrient, and recording burdens to scDesign3.
Its default biological rates are working assumptions that require calibration;
the model contract, smoke-test command, outputs, and current boundaries are in
[`ORGANOID_SIMULATION.md`](ORGANOID_SIMULATION.md).

## 5. JSON parameter file

`example_json_params/short_test.json` is a complete reference example. The
top-level keys fall into the following groups.

### 5.1 Population dynamics
- `num_init_cells` — number of founder cells (typically `1`). Multiple founders
  receive independent division schedules and share a synthetic time-zero root
  in ground-truth Newick output.
- `sim_length` — array of stopping points; the simulation captures output at
  each entry.
- `time_inc` — simulation timestep. `"auto"` derives a GCD across all per-cell-type
  cell-cycle lengths.
- `num_cores` — parallel workers used by `mclapply` inside `multi_core_func`.

### 5.2 Mitochondrial genome dynamics
- `mito_genome_length` (e.g. 16600 for human mtDNA),
- `average_genomes_per_mito`, `starting_mito_per_cell`,
- `max_mito_per_cell` — retained as a saturation setting in the parameter
  schema, but not currently enforced by `mito_dynamics()`; see the bug audit,
- `mito_inheritance_pattern` — currently `"random"` only (governs how
  mitochondria are distributed to daughter cells in `mito_dynamics`; other
  values fail explicitly because those models are not implemented),
- `fusion_events_per_mito_per_division`,
  `split_events_per_mito_per_division` — Poisson rates for fusion / fission
  events per cell cycle,
- `post_mitotic_mt_deletion_frac` — fraction of mitos randomly removed after
  division.

### 5.3 Heteroplasmy initialisation
- `baseline_heteroplasmy_sites_frac` — fraction of mt positions seeded with a
  variant in the founder cell,
- `baseline_heteroplasmy_variant_frac_dist` — `[shape1, shape2]` of a Beta
  distribution drawn per site to give the per-genome penetrance,
- `heteroplasmy_variant_transition_prob` — fraction of variant genomes that
  carry the transition base; the remainder split equally between the two
  transversion bases,
- `heteroplasmy_standard_deviation` — σ of the bimodal Gaussian used in
  `draw_severity_scores`,
- `fraction_deleterious_heteroplasmy_variants` — weight on the negative-fitness
  Gaussian mode,
- `init_heteroplasmy_survive_prob` — anchor probability for
  `logistic_prob_survive_given_score`,
- `consider_cell_heteroplasmy_scores` — boolean; when false, heteroplasmy is
  tracked but does not modulate cell survival.

### 5.4 Barcode (lineage) layout
- `bc_length`, `barcode_sequence` (path to a fixed sequence; `null` to
  randomise), `bc_nuc_composition.frac_{a,g,c,t}`,
- `max_bc_ints_per_cell` — array of integration counts per cell to evaluate,
- `include_bc_umis` — prepends a 15 bp random UMI per integration (used as an
  alignment scaffold; reference sequence is written to
  `output/processed_fastas/<run_id>/reference_seqs/`),
- `force_transversions` — global flag passed through the simulation.

### 5.5 Editing systems

Two parallel blocks with the same shape, one for nuclease editing and one for
base editing:

```
"nuclease_targets" / "be_targets": {
    "num_targets": <int|null>,
    "edit_rate_class_fractions": { "high": 0.4, "medium": 0, "low": 0.6 },
    "editing_window": {
        "size": <int>,
        "decaying": <bool>,
        "close_after_edit": <bool>
    },
    "config": "S:<first_pos>:<bases_btwn>",   // S = spaced; R/U also accepted
    "interdeletion_dropout_radius": <int>,    // nuc only
    "interdeletion_dropout_prob": <float>,    // nuc only
    "prime_editing_system": <bool>,           // nuc only
    "num_unique_prime_editing_guides": <int>, // nuc only
    "prime_editing_guide_length": <int>       // nuc only
}
```

Targets are split into High/Medium/Low edit-rate classes; editing-window
expansion either propagates the same rate or *decays* it by 1–2 "degrees"
(High→Medium→Low→background) according to `drop_editrate`.

`be_conversion_pattern` (e.g. `"A --> G"`) declares which base→base substitution
the BE produces; `classify_be_mutation_type` then dispatches the BE rate into
the transition or transversion probability list.

### 5.6 Induction events
- `differentiation_induction.{timepoint, num_cells, frac_cells}` — when (and
  to how many cells) the *induced* cell-type transition matrix is applied,
- `editing_induction.{timepoint, num_cells, frac_cells}` — when editing rates
  switch from `uninduced_editing_params` to `induced_editing_params`.

### 5.7 Cell type dictionary

```
"cell_type_dict": {
    "founder_cell_type": "ct1",
    "cell_type_params": { "ct1": {...}, "ct2": {...}, ... },
    "uninduced_transition_matrix": [[...], ...],   // square, |cell_types| x |cell_types|
    "induced_transition_matrix":   [[...], ...]
}
```

Each cell type entry contains:

- `cell_cycle_length`, `death_per_cell_cycle_prob`, `sampling_fractions`
  (per-stopping-point downsampling), `mt_invariant_sites`, `bc_invariant_sites`,
- two parameter blocks `induced_editing_params` and `uninduced_editing_params`,
  each declaring:
  - `mt_substitution_model`, `mt_sub_model_params` (and the same for `bc_*`) —
    one of `JC`, `K80`, `K81`, `F81`, `HKY`, `GTR`, with model-specific
    semicolon-separated parameter strings parsed by
    `parse_sub_model_params` (substitutions models defined in
    `substitution_models.r`),
  - `bc_bg_{insertion,deletion}_prob_per_division` and the matching `mt_*` —
    background indel rates per cell cycle (converted to per-timestep by
    `estimate_prob_per_timept`),
  - `be_mutations_per_target_per_division`,
    `nuc_{insertions,deletions}_per_target_per_division` — target-class mean
    rates fed into `SIMPLIFY_target_site_gamma_based_sub_rates` in
    `nonuniform_muts_heterogeneous.R` (gamma shape 0.5; per-site rates drawn
    by bootstrap),
  - `mt_nontarget_heterogeneity_gamma` and
    `bc_nontarget_heterogeneity_gamma` — `{shape_param, num_bins, agg_metric}`
    for `nontarget_scale_gamma_heterogeneity`, which adds across-position
    heterogeneity to the non-target mutation probabilities.

Transition matrices are converted to nested lists by
`make_cell_type_transition_lists` in `parse_cell_type_specific_args.r`.

### 5.8 Output / reconstruction controls
- `scoremat_collapse_deletions`, `binarize_mutation_scores`,
  `mt_allelic_fraction_thresholds` — score-matrix construction options
  (`mut_to_scoremat.r`),
- `reconstruction_method` — `"score"`, `"fasta"`, etc. (drives which artefacts
  are written),
- `recon_modality` — `"bc"`, `"mt"`, or both,
- `combine_mt_bc` — whether to concatenate barcode and mt score matrices,
- `mt_genome_recovery_prob`, `bc_integration_recovery_prob` — per-cell sampling
  probabilities applied at recovery time,
- `fasta_type` — e.g. `"terminal"` to export only end-of-simulation sequences,
- `include_var_pos_fasta`, `plot_heatmaps`, `savename`, `random_seed`.

## 6. Output layout

All artefacts live under `output/`. Each simulation creates a fresh
`<run_id>` directory under most subtrees:

| Path | Producer | Contents |
| ---- | -------- | -------- |
| `output/run_logs/<run_id>/runlog_*.txt` | `sim5_code.R` | per-run log file |
| `output/run_specs/<run_id>/` | `sim5_code.R` | copy of the input JSON / parsed CLI args |
| `output/cell_populations/<run_id>/cell_population_*.rds` | `sim5_code.R` | full cell population snapshot at each stopping point |
| `output/mut_profiles/<run_id>/` | `sim5_code.R` | RDS of mt and bc mutation profile matrices |
| `output/processed_fastas/<run_id>/` | `mut_to_fasta_difflen_ints.r` | per-cell FASTA sequences, plus `reference_seqs/` |
| `output/processed_newicks/<run_id>/` | `sim5_code.R` | ground-truth lineage trees in Newick |
| `output/processed_lists/<run_id>/` | `sim5_code.R` | RDS of subsetted profile lists used to produce FASTAs |
| `output/score_mats/<run_id>/matrices/` | `mut_to_scoremat.r` | RDS sparse score matrices (incl. `af/` for allelic-fraction variants) |
| `output/score_mats/<run_id>/phylips/` | `mut_to_scoremat.r` | score-matrix files consumed by IQ-TREE in `BIN/MF` mode (currently saved as `.fasta`) |
| `output/score_mats/<run_id>/csvs/` | `convert_scoremats_to_csvs.r` | CSV exports of the RDS score matrices |
| `output/recon_trees/<run_id>/{fasta_*,score_*}/` | bash wrapper (`muscle` + `iqtree`) | reconstructed trees |
| `output/tree_images/<run_id>/` | `compare_trees_call_from_bash.r` | rendered tree comparison PNG/PDFs |
| `output/rf_dist_files/<run_id>/*.txt` | `compare_trees_call_from_bash.r` | one-line RF metrics per reconstruction |
| `output/param_results_files/<run_id>/stacked_results.csv` | `process_results_from_bash.r` | parsed RF results merged with input params |
| `output/induction_details/<run_id>/` | `sim5_code.R` | which cells were induced (differentiation / editing) at each event |
| `output/timing_obj/<run_id>/`, `output/plots/runtime_plots/<run_id>/` | `sim5_code.R` | per-timepoint runtime CSVs and plots |
| `output/lineplots/<run_id>/` | `make_lineplot()` | sim-length and cell-count line plots per modality |
| `output/sc_profiles/<run_id>/` | `generate_sc_profiles_from_bash.r` | one `sim_sce_all.rds` plus per-stopping-point `sim_sce_time_<t>.rds` (`SingleCellExperiment`s with simulated counts), `cell_metadata.csv`, `cell_count_summary.csv`. Only written when `-r` is supplied. |
| `output/scdesign3_fits/` | `generate_sc_profiles_from_bash.r` | cached scDesign3 fits keyed by SHA1 of reference path + mtime + size + formula bits. Shared across runs. |
| `output/physicell/<lineage_stem>/` | `simulate_physicell_lineage.R` | imported event tree, terminal metadata, barcode and mitochondrial profiles/events, allele/heteroplasmy matrices, and FASTA/Newick exports |
| `output/physicell_10000_<timestamp>/` | `run_physicell_10000_pipeline.sh` | staged PhysiCell build, raw lineage export, logs/manifests, and combined lineage-recording output; optional scDesign3 output is under `lineage_recording/sc_profiles/` |

## 7. Source files

### Entry points
- **`bash_wrapper_all_combos.sh`** — top-level orchestrator described in §4.
- **`run_physicell_10000_pipeline.sh`** — end-to-end staged PhysiCell build,
  10,000-cell tumor run, dual-modality lineage-recording replay, and optional
  scDesign3 expression generation.
- **`simulate_physicell_lineage.R`** — dependency-light PhysiCell adapter.
  Reads `time,parent_ID,daughter_ID`, converts persistent parent IDs into
  event-resolved binary branches, and writes simulated barcode and/or
  mitochondrial observations using the selected remote_mito JSON parameters.
- **`generate_physicell_sc_profiles.R`** — fits/reuses scDesign3 on a real
  reference and simulates terminal-cell counts using PhysiCell-linked
  covariates.
- **`sim5_code.R`** *(~3.9k lines)* — the simulator. Parses the JSON, builds
  per-cell-type substitution probability matrices and per-position target /
  non-target mutation probability lists, initialises the heteroplasmy sparse
  matrix (`Matrix::sparseMatrix`), then iterates cell divisions across the
  requested stopping points using `multi_core_func` (wraps
  `parallel::mclapply` over cells). Output is written through
  `all_processes_at_stopping_point` and friends. Top-level helpers include
  `parse_target_config`, `parse_target_count_arguments`, `parse_be_example`,
  `classify_be_mutation_type`, `create_bc_sequence`, `parse_sub_model_params`,
  `generate_transition_basepos_list`, `generate_transversion_basepos_list`,
  `get_new_be_targets` / `get_new_nuc_targets` (editing-window expansion),
  `estimate_prob_per_timept`, `draw_severity_scores`,
  `logistic_prob_survive_given_score`, `mito_dynamics`,
  `reassign_genome_inds`, `sample_induced_cells`,
  `initialize_founder_population`, `setup_sim`, `multi_core_func`,
  `create_ground_truth_tree`, `all_processes_at_stopping_point`, `join_endpoint_results`,
  `make_lineplot`. The script ends by calling `setup_sim` once and then
  iterating `multi_core_func` over `poss_times`.

### Sourced modules (loaded by `sim5_code.R`)
- **`fit_plot_parameters.R`** — standalone legacy utility (not sourced by the
  simulator) that fits smoothing splines from
  `imported_heatmap_plotval_dat.csv` so that `get_heatmap_params(num_cells)`
  returns suitable y-position, height, font size, and figure dimensions for
  cell-population heatmaps. Writes a diagnostic plot
  `diagnostic_fitvals_plots.png`.
- **`nonuniform_muts_heterogeneous.R`** *(~1.1k lines)* — the position- and
  cell-type-aware mutation engine. Key functions:
  `filter_elig_ints_by_edit_window` (drops integrations once any base in their
  editing window has been edited), `non_uniform_editing` (applies position-
  specific target rates respecting saturation), `get_background_edit_inds`
  (samples non-target mutations across the matrix), `transition_func` /
  related routines for converting integer-encoded mutation outcomes into
  `(i, j, x)` triplets for sparse matrix construction. Includes
  `SIMPLIFY_target_site_gamma_based_sub_rates` (gamma-based per-target rate
  bootstrap) and `nontarget_scale_gamma_heterogeneity` used by `sim5_code.R`.
- **`substitution_models.r`** — pure-function definitions of
  `jc_sub_rate_mat`, `k80_sub_rate_mat`, `k81_sub_rate_mat`,
  `f81_sub_rate_mat`, `hky_sub_rate_mat`, `gtr_sub_rate_mat`. Returns 4×4 rate
  matrices over (A, G, C, T). Selected by name in
  `parse_sub_model_params`.
- **`add_intervening_be_targets_to_seq.r`** — `generate_non_be_target_sequence`
  builds a barcode backbone with the requested nucleotide composition and a
  reserved count of target sites; companion logic inserts the targets at
  spaced positions.
- **`mut_to_fasta_difflen_ints.r`** — converts mutation matrices to FASTA.
  Handles the integer encoding (0 = WT, –1 = deletion, 1–4 = substitution,
  decimals = insertions) via `ins_to_charvec`, samples integrations per cell
  with `get_profiles_ints_and_umis`, and writes per-cell sequences with
  `write_all_cell_sequences` (parallelised over cells).
- **`mut_to_scoremat.r`** — turns mutation matrices into score matrices.
  `group_deletions` collapses contiguous deletions, `score_mat_to_phylip`
  emits PHYLIP format, `new_scoremat_to_fasta` writes character-state FASTA,
  and `new_create_one_score_mat` is the main entry point that produces
  binarised and/or allelic-fraction sparse score matrices and saves them under
  `output/score_mats/<run_id>/matrices/{,af/}`.
- **`parse_cell_type_specific_args.r`** —
  `make_cell_type_transition_lists(cell_type_names)` lifts the JSON's
  `induced_transition_matrix` and `uninduced_transition_matrix` into nested
  named lists keyed by source/target cell type.
- **`make_babette_tree.r`** — `babette` (BEAST 2) wrapper:
  `extract_state_file_path`, `extract_newick_data`,
  `phylo_obj_from_newick`, `fasta_to_phylo`. **Currently not loaded** — the
  `source(...)` line in `sim5_code.R` is commented out.

### Post-processing scripts
- **`compare_trees_call_from_bash.r`** — `optparse` CLI invoked once per
  reconstructed treefile. Loads the reconstructed tree and the simulator's
  ground-truth tree, computes the Robinson–Foulds distance with `phangorn`,
  renders the comparison via `ape`, and writes:
  - `output/tree_images/<run_id>/...` (PNG/PDF),
  - `output/rf_dist_files/<run_id>/<savename>.txt` (one-line RF value plus
    a path back to the JSON params).
- **`process_results_from_bash.r`** — aggregator. `extract_subrun_details`
  parses the savename pattern (timepoint, modality, integration counts,
  recovery rates), `process_one_subrun` joins each RF file with its source
  JSON parameters; the script stacks everything into
  `output/param_results_files/<run_id>/stacked_results.csv`.
- **`convert_scoremats_to_csvs.r`** — small `optparse` utility; iterates RDS
  files in `--score_mat_path` and writes CSVs under
  `output/score_mats/<urid>/csvs/`.
- **`compare_celltype_probs.r`** — diagnostic. Reads
  `celltype_prob_concordance/<urid>/*.rds` (only produced when the
  currently-commented "saveRDS" block in `sim5_code.R` is re-enabled) and
  builds barplots + pivot tables comparing transition / indel rates across
  cell types.

### Single-cell profile generation
- **`scdesign3_helpers.R`** — corrected shared scDesign3 reference loader,
  cache key, marginal/copula fitting, and `new_covariate` simulation helpers.
- **`physicell_scdesign3.R`** — joins PhysiCell terminal/tree/spatial data with
  barcode edit burden and mitochondrial heteroplasmy summaries, aligns
  pseudotime/spatial ranges to the reference, and writes linked SCE outputs.
- **`generate_sc_profiles_from_bash.r`** — `optparse` CLI that bolts scDesign3
  onto the simulator. Loads every `cell_population_*_time_*.rds` snapshot for a
  given run id, builds a per-cell covariate frame (`cell_id` = simulator
  `linstring`, `cell_type` from the simulator's `celltype` field, `lineage_depth`
  = number of underscores in the linstring, `pseudotime` = `lineage_depth /
  max_depth`, `birth_time`, `induced_editing`, `timepoint`, `is_terminal`),
  with a unique `sample_id` for each cell-timepoint observation,
  loads a reference scRNA-seq dataset (SCE `.rds`, Seurat `.rds`, or `.h5ad`),
  fits scDesign3 once (`construct_data` -> `fit_marginal` -> `fit_copula`,
  cached under `output/scdesign3_fits/`), then calls `extract_para` + `simu_new`
  with the simulator-derived `new_covariate` to produce a counts matrix. The
  result is written as one `SingleCellExperiment` covering every cell at every
  stopping point plus per-timepoint splits. Required CLI flags are `-I/--run_id`
  and `-R/--reference`. Cell-type label mismatches between simulator (`ct1`,
  `ct2`, …) and reference can be resolved with `--celltype_map <json>`. Pseudotime
  is opt-in via `--use_pseudotime`; the reference must then carry a
  `pseudotime` colData column (or whatever `--pseudotime_col` names).

### Notebook
- **`parse_rf_results.ipynb`** — Jupyter notebook for post-hoc visualisation.
  Reads `output/rf_dist_files/<run_id>/`, joins to the JSON params, and emits
  faceted heatmaps of accuracy across integration count, recovery rate, and
  timepoint. Helper functions include `make_results_df` and `make_heatmap`.

### Auxiliary shell utility
- **`make_igv_input.sh`** — independent of the main pipeline. Takes a
  reference FASTA and a simulated sequence FASTA, indexes the reference with
  `samtools faidx`, converts the simulated FASTA to FASTQ (uniform quality
  scores via `awk`), aligns with `bwa mem`, and produces a sorted, indexed
  `_IGV_input.bam` for visual inspection.

## 9. scDesign3 single-cell profile step (optional)

For a PhysiCell run, pass `--sc-reference` to
`run_physicell_10000_pipeline.sh`, or invoke the stage independently:

```bash
Rscript generate_physicell_sc_profiles.R \
  --run-dir output/physicell_10000_<timestamp> \
  --reference data/tumor_reference.rds \
  --celltype-col cell_type \
  --use-pseudotime \
  --other-covariates barcode_edit_fraction
```

Every requested model covariate must exist in both the reference `colData` and
the generated PhysiCell metadata. Available generated fields include
`lineage_depth`, `lineage_pseudotime`, `birth_time`, `branch_length`, `x`, `y`,
`z`, `tumor_radius`, `neighbor_count`, `barcode_edit_count`,
`barcode_edit_fraction`, `mt_variant_count`, and
`mt_heteroplasmy_burden`. Covariates are always retained as metadata but only
affect expression when requested through pseudotime, spatial options,
`--other-covariates`, or an explicit `--mu-formula`.

PhysiCell expression outputs are:

- `sim_sce_final.rds` — feature-by-terminal-cell `SingleCellExperiment`;
- `simulated_counts.rds` — the same count assay, sparse when `Matrix` is
  available;
- `cell_metadata.csv.gz` and `cell_count_summary.csv.gz`; and
- `scdesign3_manifest.csv.gz` — reference, model formula, package version, and
  seed.

PhysiCell-linked R stages gzip CSV tables by default. Use
`--compress-csv false` with either R CLI for legacy `.csv` output; downstream
readers accept both forms.

The PhysiCell pipeline also writes
`mutation_event_descendant_matrix_sparse.rds`, a literal sampled-cell by
unique-event ground-truth ancestry matrix, and
`mutation_event_descendant_manifest.csv.gz`, which maps every column back to
its barcode or mitochondrial mutation-log row.

Behaviour: every current leaf cell in `cell_population_*_time_*.rds` for which
`alive == TRUE` and `terminal == TRUE` is treated as a sampled cell. Each gets a
row of simulated counts. Cells observed at earlier stopping points later become
internal ancestors; `is_terminal` in the expression output is TRUE only for
observations from the last entry of the JSON `sim_length` array.

What scDesign3 needs that the simulator does not provide:

- A real reference dataset (`.rds` SCE/Seurat or `.h5ad`) to fit gene-level
  marginals from. Without one, no expression can be synthesised.
- A reference `colData` column matching the simulator's cell-type labels — or a
  `--celltype_map` JSON of the form
  `{ "ct1": "stem", "ct2": "progenitor", ... }`.

Pseudotime: when `--use_pseudotime` is set, the reference fit uses a smooth
formula `cell_type + s(pseudotime, k = 4, bs = "cr")` and the simulator-side
pseudotime is `lineage_depth / max(lineage_depth)` per run. If the reference
does not have a `pseudotime` colData column, the script will refuse to run
with this flag.

Required R packages (NOT in `heavy_sim_env.yml`; install separately):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("scDesign3", "SingleCellExperiment", "zellkonverter"))
install.packages(c("optparse", "rjson", "digest"))
# only if you'll feed Seurat .rds files:
# install.packages("Seurat")
```

Outputs (under `output/sc_profiles/<run_id>/`):

- `sim_sce_all.rds` — `SingleCellExperiment` covering every cell at every
  stopping point.
- `sim_sce_time_<t>.rds` — same, split per stopping-point timepoint.
- `cell_metadata.csv` — `cell_id, sample_id, cell_type, lineage_depth, birth_time,
  induced_editing, timepoint, is_terminal, pseudotime`.
- `cell_count_summary.csv` — `(timepoint, cell_type, is_terminal) -> n_cells`.

The fit is cached at `output/scdesign3_fits/fit_<sha1>.rds`. The hash combines
the reference's normalised path, mtime, size, the chosen `--celltype_col`, the
`--use_pseudotime` flag, and the `--pseudotime_col` name; change any of these
and a fresh fit is produced.

## 10. Notes and limitations

- Score-matrix files in `phylips/` are presently saved as `.fasta`; the
  commented `.phy` lines in `bash_wrapper_all_combos.sh` are kept around for
  future reuse.
- Several diagnostic outputs are gated behind `consider_cell_heteroplasmy_scores`
  or behind currently-commented-out `saveRDS` blocks (e.g. the inputs that
  `compare_celltype_probs.r` consumes).
- Per the original `README.md`, the EM-related work referenced elsewhere in
  the project has not been committed to this repository.
