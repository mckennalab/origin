# remote_mito_clean

This document describes the repository’s simulation and analysis workflow.
The code is split into an installable simulator package (`origin/`), a sourced
analysis tier (`analysis/`), and the legacy sim5 engine (`legacy/`); see
[`REPO_LAYOUT.md`](REPO_LAYOUT.md) for what lives where and which entry point
to use.
See [`FUNCTION_REFERENCE.md`](FUNCTION_REFERENCE.md) for the complete contract
of every named function, including internal helpers. Every R source file also
carries its documentation inline: a header block stating the file’s purpose,
prerequisites, command-line arguments and outputs, and a roxygen (`#'`) block
above each function giving its parameters, return shape, and side effects. See
[`BUG_AUDIT.md`](BUG_AUDIT.md) for the repository-wide defect audit, fixes, and
remaining model-level limitations. See
[`PHYSICELL_INTEGRATION.md`](PHYSICELL_INTEGRATION.md) for growing a
lineage-enabled PhysiCell tumor and simulating barcode plus mitochondrial
recording data on its branches. See
[`ORGANOID_SIMULATION.md`](ORGANOID_SIMULATION.md) for the 2,500-founder
neural-organoid workflow.

## 0. Install and run

### The simulator package

`origin/` is a standard R package and installs with no dependencies beyond
`Matrix`.

The canonical repository is <https://github.com/mckennalab/origin>, which is
private; you need access to it, and `remotes` needs a GitHub token (see
`?gitcreds::gitcreds_set`). The package lives in the `origin/` subdirectory of
that repository, not at its root.

Straight from GitHub, without cloning:

```r
install.packages("remotes")
remotes::install_github("mckennalab/origin", subdir = "origin")
```

Or from a clone, where the commands below assume you are in the repository
root:

```bash
git clone https://github.com/mckennalab/origin.git
cd origin
R CMD INSTALL origin
```

Running `R CMD INSTALL origin` from anywhere other than the repository root, or
against a checkout of a branch that predates the package, gives

```
In normalizePath(path) : path[1]="origin": No such file or directory
```

Then:

```r
library(origin)
params <- origin_params(sim_length = 10, founders = 20, death = 0.45,
                        barcode_length = 100, integrations = 5,
                        be_rate = 0.01, targets = 20)
run_gillespie_lineage_pipeline(params, params_path = NULL,
                               output_dir = "run",
                               overrides = list(modalities = "barcode"))
```

`vignette("gillespie-quickstart", package = "origin")` is the place to start.
It covers building parameters, the two relations the parameter values do not
predict, which of the three output matrices tree building reads, and how to
register a recorder.

### Without installing

The analysis tier is sourced rather than installed, and the clique harness
replaces simulator functions at run time, which only works in the global
environment. For either, use the loader instead of `library(origin)`:

```r
source("load_origin.R")            # core + analysis
origin_include_analysis <- FALSE
source("load_origin.R")            # core only
```

This needs a clone, and R's working directory must be the repository root.
`load_origin.R` resolves its own path, so `source("/full/path/load_origin.R")`
works from anywhere.

### What else you need, and only when

The package itself needs `Matrix`. Everything beyond that is for a specific
tier, so install it when you reach that tier rather than up front:

| for | packages |
|-----|----------|
| the analysis tier (`analysis/`) | `ape`, `phangorn`, `mclust`, `data.table`, `jsonlite` |
| figures | `ggplot2`, `ragg`, `patchwork` |
| tests | `testthat` |

```r
install.packages(c("ape", "phangorn", "mclust", "data.table", "jsonlite",
                   "ggplot2", "ragg", "patchwork", "testthat"))
```

Tree reconstruction is a separate checkout, `clique_2025_12_10`, holding the
`cliqueR` package and the benchmark harness. Its methods shell out to external
tools: IQ-TREE 2 (`iqtree2`), PHYLIP (`mix`), and Cassiopeia as a Python module
reached through `reticulate`. Neighbour joining, parsimony and VINE need none of
those, so a run limited to them has nothing external to install.

The conda environment in section 3 belongs to the legacy `sim5` engine in
`legacy/`, not to this package. It pulls BEAST 2, `babette`, `muscle` and
`samtools`, none of which `origin` uses.

### Reproducing the large-scale runs

`analysis/runs/` holds one driver per family of simulation, each running
simulate, reconstruct, summarise and plot end to end. Start with the source
benchmark every other run reads from:

```bash
bash analysis/runs/run_0_source_benchmark.sh
export SOURCE_BENCHMARK=$PWD/output/lineage_benchmark_<timestamp>
bash analysis/runs/run_a_recorder_parameter_grid.sh
```

See `analysis/runs/README.md` for scale, overrides and what each run produces.
Outputs land in `output/`, which is ignored.

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

## 3. Environment (legacy sim5 engine)

This section covers the legacy `sim5` engine in `legacy/`, not the `origin`
package; see section 0 for installing and running `origin`.

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

## 4. Running a simulation (legacy sim5 engine)

The command below drives the legacy `sim5` engine. For the `origin` package see
section 0; for the large-scale benchmark runs see `analysis/runs/`.

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

### 4.1.1 Replicated Visium spatial-lineage experiment

`run_physicell_visium_replicates.sh` runs independent PhysiCell tumors, takes
virtual tissue sections, overlays a conventional 6.5-mm Visium array, and
measures ground-truth lineage distance as a function of spatial distance. The
default array follows the documented conventional Visium geometry: 4,992
spots in 78 staggered rows with 64 spots per row, 55-micron spots, and
100-micron center-to-center spacing. It models a 5-micron section and randomly
rotates/translates the array for every section.

Run five 10,000-cell tumor replicates with three sections per tumor:

```bash
bash run_physicell_visium_replicates.sh \
  --replicates 5 \
  --target-cells 10000 \
  --slice-offsets=-50,0,50 \
  --jobs 8
```

The wrapper requests `--modalities lineage` by default. This writes the exact
event tree and skips barcode/mitochondrial replay, which is unnecessary for a
ground-truth spatial-lineage analysis. Set `--modalities both` to retain those
simulated observations as well. Completed replicates and analyses are reused
when the command is restarted.

Analyze an already completed tumor without rerunning PhysiCell:

```bash
Rscript analyze_physicell_visium.R \
  --run-dir output/physicell_10000_20260726_171340 \
  --slice-offsets=-50,0,50
```

Every `visium_spatial/slice_*` directory contains:

| File | Contents |
| --- | --- |
| `visium_spots.csv.gz` | All 4,992 synthetic spots with Space Ranger-style array row/column coordinates, physical coordinates, tissue flag, and cell count. |
| `slice_cells.csv.gz` | Every cell center intersecting the section, including cells falling in the gaps between spots. |
| `spot_cell_membership.csv.gz` | Cells captured by the 55-micron spot footprints and their spot assignments. |
| `cell_pair_lineage_sample.csv.gz` | Cell pairs with 2-D/3-D distance, MRCA age, patristic distance, division distance, founder identity, and spot identity. |
| `cell_distance_summary.csv.gz` | Cell-distance correlogram and tip-label-permutation null enrichment. |
| `spot_pair_lineage.csv.gz` | Mean lineage relationship for each sampled within-spot or between-spot cell mixture. |
| `spot_distance_summary.csv.gz` | Equal-spot-pair lineage summary as a function of spot-center distance. |
| `visium_slice.png` | Section cells and occupied spot footprints. |
| `cell_lineage_correlogram.png` | Observed and spatially permuted mean MRCA age curves. |

The batch-level `aggregate/` directory first pools slices within each tumor
replicate, then reports the mean, standard deviation, and standard error across
independent tumors. This avoids treating the many cell pairs from one tumor as
independent biological replicates. Pair queries use a preprocessed parent
table rather than a complete quadratic tree-distance matrix. Above the exact
pair limit, uniformly sampled pairs are retained separately within each
distance bin; the sampling mode and evaluated counts are written to the
summary.

Conventional Visium is not single-cell resolution. The exported membership
table therefore records the simulated cells contributing to each spot; it
does not relabel a mixed spot as one cell. This stage models spatial capture
geometry but not molecule-level RNA diffusion, UMI sampling, or histology-
based tissue detection.

A 10,000-cell PhysiCell spheroid is only a small region relative to a 6.5-mm
capture area; the existing example occupies roughly 7--9 spots per central
section. It is suitable for pipeline testing, but `--target-cells 100000` or
larger is preferable when the goal is a richer multi-spot spatial curve. The
spot table still contains all 4,992 array positions so the physical scale is
never silently rescaled to make a small tumor appear slide-sized.

Open [`notebooks/physicell_visium_lineage_distance.Rmd`](notebooks/physicell_visium_lineage_distance.Rmd)
in RStudio for an interactive notebook of physical distance versus MRCA age,
patristic/lineage-edge distance, recent-relative enrichment, within-spot versus
between-spot relationships, projection effects, and per-section statistics.
Set the `input_path` parameter to either one `visium_spatial/` directory or a
complete replicate batch. Render the same notebook from the command line with:

```bash
Rscript render_physicell_visium_notebook.R \
  --input-path output/physicell_visium_YYYYMMDD_HHMMSS \
  --output-file output/physicell_visium_report.html
```

For a replicate batch, sections are averaged within each tumor before the
notebook calculates across-tumor means and standard errors. The dense pair
plots are reproducibly downsampled only for display; their binned curves use
the full summaries saved by the analysis.

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

## Exact Gillespie lineage simulation

The native, non-spatial population simulator can now run either on the legacy
time grid or as an exact continuous-time birth/death process. The timestep
engine remains the default. To select Gillespie in an existing parameter file,
add:

```json
"simulation_engine": "gillespie",
"gillespie": {
  "modalities": ["barcode", "mitochondrial"],
  "max_cells": 1000000,
  "num_integrations": 5,
  "founder_label_sites": 0,
  "mt_genomes_per_cell": 8,
  "progress": true,
  "progress_updates": 20,
  "compress_csv": true
}
```

The usual entry point dispatches from the JSON automatically:

```bash
Rscript sim5_code.R -P parameters.json
```

It can also be run directly, with command-line values overriding the optional
`gillespie` block:

```bash
Rscript simulate_gillespie_lineage.R \
  --params parameters.json \
  --modalities both \
  --end-time 21 \
  --output-dir output/gillespie/day21
```

Division has hazard `1 / cell_cycle_length`. A configured death probability
`p` per cycle becomes the exact hazard `-log(1-p) / cell_cycle_length`.
Daughter types are sampled from the existing induced or uninduced transition
matrix at division. Editing and differentiation induction occur at their exact
configured times; selected cells get explicit continuation nodes so partial
induction states are inherited without moving the boundary to a timestep.

Conditional on this continuous-time tree, barcode and mitochondrial mutation
processes are simulated at exact event times with cell-type-specific hazards.
This separation is mathematically equivalent to one joint Gillespie process
because recording mutations do not change population propensities. For that
reason, Gillespie currently rejects `consider_cell_heteroplasmy_scores: true`.
The mitochondrial overlay uses a fixed-genome division bottleneck and does not
replay the legacy organelle fusion/fission counts.

In addition to the standard lineage, recorder, sparse feature, literal
event-descendant, Newick, compressed CSV, and timing artifacts, the engine
writes `gillespie_population_events.csv.gz`,
`gillespie_division_events.csv.gz`, `gillespie_checkpoint_summary.csv.gz`, and
`gillespie_cell_states.csv.gz`. The latter is accepted directly by the
covariate-linked scDesign3 entry point:

```bash
Rscript generate_physicell_sc_profiles.R \
  --run-dir output/gillespie/day21 \
  --reference data/reference.rds \
  --use-pseudotime
```

### Quick continuous-time versus time-step comparison

Use the population-only comparison script to quantify discretization effects
without paying the cost of barcode, mitochondrial, tree, or transcriptome
output for every replicate:

```bash
Rscript compare_simulation_engines.R \
  --params example_json_params/short_test.json \
  --replicates 100 \
  --time-step 1 \
  --output-dir output/engine_comparison/short_test
```

Both engines use the same division/death hazards, differentiation induction,
and daughter cell-type transition matrices. The time-step version permits at
most one competing division/death event per cell per interval, with event
probability `1-exp(-(division_hazard+death_hazard)*dt)`; this is the discrete
population approximation used to expose grid-size bias. The comparison does
not run either recording overlay.

Outputs include per-replicate outcomes and cell-type counts, an aggregate table
of continuous/time-step means and differences, and
`engine_comparison.png`. Repeat with successively smaller `--time-step` values
to check convergence toward continuous time.

## Paired lineage-recorder benchmark simulation

`run_lineage_benchmark.R` creates simulation inputs for a separate tree-
reconstruction benchmark. It does not infer trees. The default design uses:

- ground-truth samples of 250, 1,000, 2,000, and 5,000 cells;
- 1, 2, 5, 10, and 20 integrations;
- BASELINE-like recording with 50 primary targets per integration;
- prime editing with six fixed pegRNAs per integration and efficiencies from
  0.90 to 0.15;
- PALINCODE with two cBits per integration;
- mitochondrial lineage tracing with a fixed 32-genome intracellular pool and
  nested observations of 1, 2, 5, 10, or 20 sampled genomes per cell;
- balanced, comb-like, neutral asynchronous, and hierarchical
  stem/progenitor/terminal population shapes; and
- ten population seeds.

Run the complete default benchmark with:

```bash
Rscript run_lineage_benchmark.R \
  --output-dir output/lineage_benchmark
```

The run is resumable by default. To run a small pilot:

```bash
Rscript run_lineage_benchmark.R \
  --output-dir output/lineage_benchmark_pilot \
  --shapes balanced,neutral \
  --tree-sizes 250,1000 \
  --integration-counts 1,5,20 \
  --mt-observation-depths 1,5,20 \
  --seeds 1:2
```

The largest population and integration panel are simulated only once for each
shape, seed, and recorder. Mitochondrial evolution is likewise simulated once
with 32 modeled genomes per cell; observation-depth conditions use the first
1/2/5/10/20 genomes from a reproducibly shuffled, without-replacement sampling
order. Smaller conditions are nested deterministic subsets of the same cells,
integrations, and mitochondrial observations. Change the modeled mitochondrial
pool with `--mt-genomes-per-cell`; every `--mt-observation-depths` value must be
no larger than that pool. The default performs 40 population simulations and
160 full recording simulations to produce 3,200 paired reconstruction-input
conditions. An optional
`turnover` shape is available through `--shapes`; use `--help` for all controls.
The mitochondrial benchmark preset uses a 16,569-base reference, no founder
heteroplasmy or indels, and a Jukes-Cantor substitution probability of
`2e-6` per alternate base per cell cycle.

Every condition directory contains:

| File | Contents |
| --- | --- |
| `ground_truth_tree.nwk` | Exact sampled ground-truth topology with elapsed-time branch lengths. |
| `recording_logical_target_matrix_sparse.rds` | Recommended cross-system matrix: categorical integrated-recorder targets or binary mitochondrial variant presence. |
| `recording_state_matrix_sparse.rds` | Native recorder state matrix; for mitochondria, sampled variant fractions (heteroplasmy). |
| `recording_character_matrix_sparse.rds` | Native reconstruction characters; for mitochondria, binary variant presence. |
| `target_manifest.csv.gz` | Integration/target definitions and logical-state encoding. |
| `sample_cells.csv.gz` | Ordered tree tips and population metadata. |
| `condition_manifest.csv.gz` | Seeds, dimensions, logical target counts, and paths to shared event output. |

For BASELINE, the logical matrix encodes the within-target editing-window
pattern as a stable nonzero integer bitmask, so one physical target remains one
logical character. Prime-editing logical states are `0/1`; PALINCODE logical
states are `0=WT`, `1=left`, `2=right`, and `3=both`. The native matrices remain
available for reconstruction methods that explicitly model base-level or
one-hot outcomes. Mitochondrial condition directories also provide explicit
aliases: `mitochondrial_variant_fraction_matrix_sparse.rds`,
`mitochondrial_binary_variant_matrix_sparse.rds`, and
`mitochondrial_variant_manifest.csv.gz`. Their manifests report
`observation_depth` and leave `integrations` unset.

Full 20-integration mutation-event tables and recorder parameters are stored
once under each `recorder_<system>/full_20_integrations/` directory.
Mitochondrial events, terminal genome profiles, reference FASTA, observation
matrices, and sampling orders are stored under
`recorder_mitochondrial/full_32_genomes_per_cell/`. The root
`benchmark_manifest.csv.gz` is the reconstruction job table and links every
condition to its ground truth and recording matrices. Dense CSV matrices are
disabled because of their size; add `--write-dense-csv true` only when needed.

## Unified prime-editing recorder

The native time-step simulator, exact Gillespie simulator, and imported
PhysiCell lineage replay now share one prime-editing backend. Recorder targets
are defined by `nuclease_targets`; the backend assigns one known pegRNA from a
pool to each target. Every pool entry has its own editing efficiency, and every
integration inherits the same target-to-pegRNA layout plus a unique static ID.

By default each target installs a single fixed edit. Setting
`marks_per_target` gives every site an alphabet of that many marks instead, one
of which is drawn when the site edits — the PEtracer design ([Science 2025](https://www.science.org/doi/10.1126/science.adx3800)),
where each of 3 edit sites per cassette receives one of 8 predefined 5nt marks.
This matters for reconstruction rather than realism: with one outcome per site,
two cells that edit the same site independently are indistinguishable from two
cells sharing an ancestor, and that homoplasy misleads tree building. With `N`
marks they differ `(N - 1) / N` of the time. Marks are taken from the pool in
blocks, so the pool needs at least `marks_per_target` entries. `pegRNA`
efficiency then sets which mark wins when the site fires, not whether it fires,
which the site's own rate governs.

See `example_json_params/petracer_gillespie.json` for a runnable PEtracer
configuration — 12 cassettes, 3 sites each, 8 marks per site, at the published
0.05–0.1 edits per site per day. Note that editing accrues per unit of simulated
time rather than per division, so expected saturation over `sim_length` T is
`1 - (1 - p)^(T / cell_cycle_length)`; that preset reaches 0.68, inside the
60–80% band the paper reports as maximising reconstruction accuracy.

The pool can be a CSV referenced relative to the parameter JSON:

```json
"nuclease_targets": {
  "num_targets": 4,
  "config": "S:5:5",
  "edit_rate_class_fractions": {"high": 1, "medium": 0, "low": 0},
  "editing_window": {"size": 0, "decaying": false, "close_after_edit": true},
  "prime_editing_system": true
},
"prime_editing_backend": {
  "enabled": true,
  "pegRNA_pool_path": "../data/example_pegRNA_pool.csv",
  "target_pegRNA_ids": ["peg_A", "peg_B", "peg_C", "peg_D"],
  "induced_edit_probability_per_cell_cycle": 0.4,
  "uninduced_edit_probability_per_cell_cycle": 0.0,
  "static_id_length": 12
}
```

The required pool columns are `pegRNA_id`, `edit_sequence`, and
`editing_efficiency`. IDs must be unique, edit sequences must contain only
A/C/G/T, and efficiencies must be in `[0,1]`. Optional columns are
`spacer_sequence`, `pbs_sequence`, `rtt_sequence`, and `description`. The same
rows may instead be embedded as an array of objects under
`prime_editing_backend.pegRNAs`.

`target_pegRNA_ids` gives an exact assignment and must contain one pool ID per
target. If it is omitted, `assignment` may be `cycle` (the default),
`sample_with_replacement`, or `sample_without_replacement`. Scalar induced and
uninduced base probabilities apply to all targets; a vector can set one base
probability per target. For base probability `p` and pegRNA efficiency `e`, the
effective per-cycle probability is `1 - (1 - p)^e`. Thus efficiency zero
disables a pegRNA, efficiency one preserves the base rate, and intermediate
values scale the event hazard.

Run the included continuous-time example with:

```bash
Rscript simulate_gillespie_lineage.R \
  --params example_json_params/prime_editing_gillespie.json
```

The same JSON block works with `simulate_physicell_lineage.R`, or with
`sim5_code.R -P parameters.json` when `simulation_engine` is `timestep` or
`gillespie`. The old `num_unique_prime_editing_guides` and
`prime_editing_guide_length` fields remain as a compatibility fallback that
generates a random, unit-efficiency pool when no known pool is supplied.

Event-resolved PhysiCell and Gillespie output uses locked `0=unedited` and
`1=edited` states, retains the exact programmed sequence in
`mutation_events.csv.gz`, and writes:

| File | Contents |
| --- | --- |
| `prime_editing_target_manifest.csv.gz` | Static integration ID, pegRNA assignment, exact edit/template fields, raw efficiency, and base/effective probabilities. |
| `prime_editing_state_matrix_sparse.rds` | Cell-by-integrated-target categorical state matrix. |
| `prime_editing_character_matrix_sparse.rds` | Binary cell-by-integrated-target character matrix for reconstruction. |
| `barcode_binary_score_matrix_sparse.rds` | Compatibility alias of the prime-editing character matrix. |

Dense runs write the corresponding `.csv.gz` matrices. The legacy time-step
path writes the exact target assignment to its run-spec
`prime_editing_target_manifest.csv` and uses the same efficiency-adjusted
target probabilities while retaining its existing insertion-profile encoding.

## PALINCODE palindromic cBits

The event-resolved PhysiCell and Gillespie pipelines support the PALINCODE
recorder described in
[`paper/2026.04.16.718941v1.full.pdf`](paper/2026.04.16.718941v1.full.pdf).
Each palindromic cBit begins wild type and irreversibly resolves at its first
event to a left edit, right edit, or rare simultaneous edit of both sides. A
resolved site cannot later acquire the opposite-side edit. This follows the
paper's evidence that the dual class is primarily a simultaneous event and
that pre-edited targets have much lower subsequent activity.

PALINCODE is selected by adding a `palincode_adapter` block:

```json
"physicell_adapter": {
  "recorder_system": "PALINCODE",
  "profile_storage": "sparse",
  "compact_output": true,
  "retain_internal_profiles": false
},
"palincode_adapter": {
  "enabled": true,
  "num_cbits_per_integration": 2,
  "cbit_names": ["PalT7", "PalRNF2"],
  "static_id_length": 12,
  "uninduced_edit_probability_per_cbit_per_cell_cycle": 0.0,
  "induced_edit_probability_per_cbit_per_cell_cycle": [0.15, 0.08],
  "left_edit_fraction": [0.495, 0.35],
  "right_edit_fraction": [0.495, 0.64],
  "both_edit_fraction": [0.01, 0.01]
}
```

Probability fields accept either one value shared by every cBit or one value
per cBit. For each target, the three outcome fractions must sum to one. The
per-cell-cycle edit probability is converted to a continuous-time hazard using
the branch's cell-type-specific cycle length. Conditional on the first event,
one left/right/both outcome is drawn and locked. The paper explored editing
rates from 0.1% to 75% per generation and reported its highest reconstruction
accuracy around 5–25%; the example intentionally uses target-specific values
in that range rather than treating them as experimentally fitted rates.

Run the eight-generation, 30-cBit example with:

```bash
Rscript simulate_gillespie_lineage.R \
  --params example_json_params/palincode_gillespie.json
```

or replay PALINCODE on a PhysiCell tree with the same parameter schema:

```bash
Rscript simulate_physicell_lineage.R \
  --lineage divisions.csv \
  --live-cells live_cells.csv \
  --params parameters_with_palincode.json \
  --modalities barcode \
  --num-integrations 15
```

State values are `0=wild type`, `1=left`, `2=right`, and `3=both`. The main
PALINCODE outputs are:

| File | Contents |
| --- | --- |
| `barcode_target_layout.csv.gz` | Static 12-nt integration IDs, cBit names, edit rates, and left/right/both fractions. |
| `mutation_events.csv.gz` | Exact event times and `palincode_left`, `palincode_right`, or `palincode_both` outcomes. |
| `palincode_state_matrix_sparse.rds` | Cell-by-cBit categorical state matrix using values 0–3. |
| `palincode_character_matrix_sparse.rds` | One-hot left/right/both characters; wild type is all zero. |
| `barcode_binary_score_matrix_sparse.rds` | Compatibility alias of the PALINCODE one-hot character matrix. |

The one-hot matrix preserves edit orientation in combined lineage features and
can be passed to reconstruction methods that expect binary characters. The
simulator intentionally records which side was edited, matching the simple
outcome encoding used for the paper's trees; it does not currently expand each
side into its individual within-window adenine-to-guanine combinations. The
legacy native timestep simulator still uses its original nucleotide barcode
mutation engine; PALINCODE currently runs through the imported-PhysiCell or
exact-Gillespie event-resolved pipelines.

## Non-Mendelian ecDNA barcodes

Both fixed PhysiCell trees and Gillespie trees can carry a third recording
modality consisting of randomly segregating ecDNA species. The nearest whole
number to a configured fraction of founder ecDNA species is labeled (or set
`num_labeled_species` for an exact count). Every labeled species receives
an immutable static ID and a small irreversible CRISPR target array; unlabeled
species propagate copy number but produce no recorder state.

Add an `ecdna_adapter` block to the native parameter JSON:

```json
"ecdna_adapter": {
  "num_species": 10,
  "initial_copies_per_species": 5,
  "labeled_species_fraction": 0.25,
  "static_id_length": 12,
  "num_recorder_targets": 6,
  "edit_probability_per_target_per_cell_cycle": 0.01,
  "recorder_start_time": 0,
  "replication_probability": 1.0,
  "daughter_1_segregation_probability": 0.5,
  "max_copies_per_cell": 10000
}
```

Selection against non-Mendelian markers is configured independently of the
cellular tree:

```json
"non_mendelian_selection": {
  "coefficient": 0.0,
  "ecdna_label_coefficient": 0.05,
  "ecdna_recorder_edit_coefficient": 0.01,
  "mitochondrial_variant_coefficient": 0.02
}
```

All coefficients use the standard range `0 <= s <= 1`; zero is neutral. The
global `coefficient` is a fallback for any omitted marker-specific value.
Labeled ecDNA copies receive relative extra-replication propensity `1-s`, and
each recorder edit contributes another factor of `1-s`. Mitochondrial genomes
are sampled into a daughter with relative weight `(1-s)^k`, where `k` is that
genome's variant count. For ecDNA, a coefficient of one blocks additional
replication but existing copies still segregate; for mitochondria, it excludes
variant-bearing genomes from the bottleneck when wild-type genomes are present.

The mitochondrial replay currently fixes total `genomes_per_cell`, so this
controls selection against mitochondrial *variant burden*, not selection on
total mitochondrial copy number. Copy-number-dependent cellular fitness would
require variable organelle counts coupled back into the population simulator.

At a true cell division, each parental ecDNA copy retains itself and produces
one additional copy with probability `replication_probability`. The joint pool
is then divided between the two daughters by binomial segregation. Thus the
daughters receive complementary, generally unequal copy numbers and can lose
a species entirely. Induction-continuation nodes do not replicate or partition
ecDNA. Recorder targets edit continuously according to the branch cell type's
cycle length; static IDs and edited haplotypes follow ecDNA copies rather than
all descendants of the cellular node.

Run ecDNA alone or with the other modalities:

```bash
Rscript simulate_physicell_lineage.R \
  --lineage divisions.csv \
  --live-cells live_cells.csv \
  --params parameters.json \
  --modalities ecdna

Rscript simulate_gillespie_lineage.R \
  --params parameters.json \
  --modalities all
```

Key outputs are:

| File | Contents |
| --- | --- |
| `ecdna_species_manifest.csv.gz` | Species, labeling status, static IDs, and founder copy numbers. |
| `ecdna_cell_summary.csv.gz` | Total/labeled copy number, retained species, and recorder burden per sampled cell. |
| `ecdna_haplotypes.csv.gz` | Per-cell species/static-ID/CRISPR-haplotype copy counts. |
| `ecdna_mutation_events.csv.gz` | Aggregated exact-time CRISPR edit occurrences with cellular branch provenance. |
| `ecdna_static_id_copy_number_matrix_sparse.rds` | Cells by labeled static-ID copy number. |
| `ecdna_recorder_edit_fraction_matrix_sparse.rds` | Cells by static-ID/target edited-copy fraction. |
| `ecdna_recorder_character_matrix_sparse.rds` | Binary cell-by-static-ID/target character matrix. |

The ecDNA character matrix is computed from actual terminal copy inheritance.
It deliberately does not mark every cellular descendant of an edit-bearing
branch, because random ecDNA segregation violates that Mendelian assumption.
The covariate-linked scDesign3 loader also adds ecDNA copy number, labeled copy
number, retained-species counts, and recorder edit fraction when these outputs
are present.

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
    "edit_rate_dispersion_shape": <float>,    // optional; default 0.5
    "editing_window": {
        "size": <int>,
        "decaying": <bool>,
        "close_after_edit": <bool>
    },
    "config": "S:<first_pos>:<bases_btwn>",   // S = spaced; R/U also accepted
    "interdeletion_dropout_radius": <int>,    // nuc only
    "interdeletion_dropout_prob": <float>,    // nuc only
    "prime_editing_system": <bool>,           // nuc only; enables shared backend
    "num_unique_prime_editing_guides": <int>, // legacy random-pool fallback
    "prime_editing_guide_length": <int>       // legacy random-pool fallback
}
```

Targets are split into High/Medium/Low edit-rate classes; editing-window
expansion either propagates the same rate or *decays* it by 1–2 "degrees"
(High→Medium→Low→background) according to `drop_editrate`.

Within that, each target's rate is drawn from a gamma distribution with the
configured mean and `edit_rate_dispersion_shape`, whose coefficient of variation
is `1 / sqrt(shape)`. Smaller values spread target rates further apart, pushing
mass towards both a dead tail and a fast tail that saturates early in the
lineage. The default of 0.5 (CV 1.41) is kept for backwards compatibility, but
measured constructs can be considerably more dispersed: a 34-barcode BASELINE
recording was reproduced at about 0.09 (CV 3.4), where 0.5 produced no saturated
targets at all against 17% observed. See
`analysis/cli/clone84_simulation_match.R`.

`be_conversion_pattern` (e.g. `"A --> G"`) declares which base→base substitution
the BE produces; `classify_be_mutation_type` then dispatches the BE rate into
the transition or transversion probability list.

Known prime-editing pools and target-specific efficiencies are configured in
the top-level `prime_editing_backend` block described under
“Unified prime-editing recorder.”

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
