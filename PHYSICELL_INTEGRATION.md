# PhysiCell lineage integration

`simulate_physicell_lineage.R` imports the division-event CSV produced by the
lineage-enabled PhysiCell checkout and simulates heritable barcode and
mitochondrial recording data along that fixed lineage. It does not resimulate
population growth: PhysiCell supplies the topology and division times, while
the remote_mito JSON supplies barcode layout, editing targets, cell-cycle
scale, heteroplasmy initialization, and mutation rates.

For the 2,500-iPSC neural-organoid model and its richer founder, fate-transition,
and terminal-state exports, see
[`ORGANOID_SIMULATION.md`](ORGANOID_SIMULATION.md).

## Complete 10,000-cell run

The repository includes an end-to-end driver for the lineage-enabled
`user_projects/tumor_3D_lineage` project in `../PhysiCell`:

```bash
bash run_physicell_10000_pipeline.sh
```

By default it:

1. creates a staged PhysiCell build under
   `output/physicell_10000_<timestamp>/physicell_build`;
2. changes only the staged copy to stop at 10,000 current cells, use the
   requested seed/thread count, and suppress intermediate full/SVG snapshots;
3. compiles and runs the 3-D tumor;
4. validates its division and current-cell CSVs; and
5. replays both recording modalities with
   `example_json_params/physicell_10000.json`.

The PhysiCell source checkout is not modified. Run `--help` for all controls;
common overrides include:

```bash
bash run_physicell_10000_pipeline.sh \
  --target-cells 10000 \
  --num-integrations 1 \
  --mt-genomes-per-cell 8 \
  --seed 1 \
  --jobs 8 \
  --output-dir output/my_physicell_run
```

Use a new or empty output directory. The pipeline writes its resolved inputs,
cell/event counts, seed, and final simulation time to
`pipeline_manifest.csv`.

## Optional scDesign3 expression

Add a real single-cell reference to generate an expression profile for every
terminal PhysiCell cell:

```bash
bash run_physicell_10000_pipeline.sh \
  --sc-reference data/tumor_reference.rds \
  --sc-celltype-col cell_type \
  --sc-use-pseudotime \
  --sc-ncores 8
```

The reference may be a `SingleCellExperiment`/Seurat RDS or an AnnData h5ad
file and must contain a count assay plus the selected cell-type column. A JSON
mapping supplied through `--sc-celltype-map` resolves labels such as
`cancer` to a different reference label.

The generated metadata joins:

- terminal sample/PhysiCell/node IDs;
- lineage depth, normalized lineage pseudotime, birth/sampling times, and
  terminal branch length;
- PhysiCell x/y/z coordinates, radius, and neighbor count;
- base-editing edit count/fraction; and
- mitochondrial variant count and heteroplasmy summaries.

`--sc-use-pseudotime` fits a smooth reference pseudotime effect and maps
normalized lineage depth into the reference pseudotime range.
`--sc-spatial-cols ref_x,ref_y` fits a 2-D spatial effect after mapping
PhysiCell x/y into those reference coordinate ranges. Additional fields can be
modeled with `--sc-other-covariates`; each must have the same name in the
reference `colData`. These links are conditional effects, not an additional
cell-to-cell phylogenetic expression covariance model.

Outputs are written under `lineage_recording/sc_profiles/`:

| File | Contents |
| --- | --- |
| `sim_sce_final.rds` | Synthetic count assay and all linked metadata, keyed by `cell_<PhysiCell_ID>`. |
| `simulated_counts.rds` | Feature-by-cell count matrix, sparse when `Matrix` is available. |
| `cell_metadata.csv.gz` | Joined lineage, spatial, barcode, mitochondrial, and model covariates. |
| `cell_count_summary.csv.gz` | Simulated terminal-cell counts by modeled cell type. |
| `scdesign3_manifest.csv.gz` | Reference path, formula, scDesign3 version, covariates, and seed. |

The fitted marginal/copula model is cached under `output/scdesign3_fits/`.
The cache key includes the reference identity, scDesign3 version, covariate
specification, and formula.

## PhysiCell inputs

The required division CSV is the output of
`PhysiCell::save_cell_lineage(...)`:

```csv
time,parent_ID,daughter_ID
60,0,12
90,0,19
120,12,31
```

All three columns are required. Events may be unsorted; the importer sorts by
time while retaining file order for simultaneous events. IDs must be
non-negative integers, every daughter ID must be born exactly once, and a
division cannot name the parent as its own daughter.

The updated PhysiCell implementation retains the parent's ID at division and
assigns a new ID only to the newly allocated daughter. A single persistent cell
can therefore appear as `parent_ID` more than once. The importer turns every
division into a binary split:

```text
active parent segment
├── continuing-parent segment (same PhysiCell ID)
└── new-daughter segment      (new PhysiCell ID)
```

This event-resolved representation gives every branch one birth time and one
end/division time, which is the form needed to inherit recording mutations
without turning repeated divisions into a polytomy.

Division validation and event-tree construction scale linearly with the
division log: validation uses vectorized ID matching, active parent segments
are resolved through a hashed cell-to-node map, and node columns are
preallocated. A division table returned by `read_physicell_divisions()` is
marked as validated so `build_physicell_lineage()` does not repeat the same
validation pass.

PhysiCell's division log does not record cell removal. Supplying the optional
output of `PhysiCell::save_cell_lineage_table(...)` with `--live-cells` limits
the generated observations and sampled tree to IDs still present at the final
PhysiCell snapshot:

```csv
ID,parent_ID,x,y,z,neighbor_IDs,branch_length_to_parent
0,-1,0,0,0,12;19,-1
12,0,1,0,0,0,60
```

Only the `ID` column is consumed. Without this table, every cell ID that remains
active at the end of the division log is treated as a sampled terminal cell.

If the output calls are not already in the PhysiCell project driver, call them
at the desired final snapshot:

```cpp
PhysiCell::save_cell_lineage("output/cell_divisions.csv");
PhysiCell::save_cell_lineage_table("output/live_cells.csv");
```

The complete division log must be saved; a live-cell table alone cannot recover
intermediate divisions.

## Run the adapter

From this repository:

```bash
Rscript simulate_physicell_lineage.R \
  --lineage ../PhysiCell/output/cell_divisions.csv \
  --live-cells ../PhysiCell/output/live_cells.csv \
  --params example_json_params/short_test.json \
  --end-time 1440 \
  --output-dir output/physicell/my_run \
  --modalities both \
  --mt-genomes-per-cell 8 \
  --seed 1
```

`--lineage` and `--params` are required. `--end-time` defaults to the largest
`sim_length` value in the JSON and must be at least the final division time.
Other options are:

- `--founder-time`: birth time assigned to parent IDs first seen in the log;
  default `0`.
- `--cell-type`: cell-type parameter block used for every imported branch;
  default `cell_type_dict.founder_cell_type`.
- `--editing-state auto|induced|uninduced`: `auto` switches all branches at
  `editing_induction.timepoint`; the other choices force one rate block.
- `--modalities barcode|mitochondrial|both`: output modality; default
  `barcode`. A comma-delimited `barcode,mitochondrial` value is also accepted.
- `--num-integrations`: integrations simulated per terminal cell; default is
  the largest value in `max_bc_ints_per_cell`.
- `--mt-genomes-per-cell`: fixed mitochondrial genome bottleneck/resampling
  size per cell; default `8`.
- `--write-mt-fasta true|false`: optionally write one sampled mitochondrial
  haplotype per terminal cell; default `false`.
- `--compress-csv true|false`: write R-generated tabular outputs as
  gzip-compressed `.csv.gz` files; default `true`.
- `--progress true|false`: print progress and timestamped stage-transition
  logs; default `true` for command-line runs.
- `--progress-updates`: approximate number of progress messages per phase;
  default `20`.
- `--live-cells`: optional final PhysiCell live-cell table.

The command requires base R and `jsonlite`; `Matrix` is used when available to
write a sparse mitochondrial variant-fraction matrix. It does not load the
full remote_mito/BEAST analysis environment.

Progress lines cover event-tree reconstruction, barcode recording, and
mitochondrial tracing, and full and sampled Newick rendering. They report
processed and total work units or nodes, percentage, elapsed time, throughput,
and estimated time remaining. The PhysiCell pipeline wrappers merge these
messages into `recording_run.log`.
Library calls to `build_physicell_lineage()`,
`simulate_recording_on_physicell_lineage()`, and
`simulate_mito_on_physicell_lineage()` remain quiet by default; set
`show_progress = TRUE` to enable the same reporter.

Newick export uses integer parent/child arrays and an explicit traversal stack,
so runtime is linear in the retained tree size and is not limited by R's
recursion depth. Tokens are written through a bounded buffer instead of
assembling and repeatedly copying complete subtree strings. The sampled tree
uses the same traversal after a reverse ancestor-marking pass; it retains unary
ancestors so pruning does not alter elapsed branch time.

After each simulation loop, a timestamped message explicitly announces that
simulation has finished and output is beginning. Subsequent messages identify
lineage-table/Newick writing, barcode matrix and event output, mitochondrial
heteroplasmy summarization and matrix output, combined-feature output, and
final completion. A successful CLI run then prints phase wall-clock durations
and total R runtime regardless of the `--progress` setting, and writes the same
values to `r_timing_summary.csv.gz`. For example:

```text
R simulator timing summary:
  Startup, input, and lineage reconstruction:  1.2s
  Lineage tables and Newick output:             6.4s
  Barcode model preparation:                    0.8s
  Barcode lineage simulation:                   4m 03s
  Barcode output:                               18.7s
  Total R simulator:                            4m 30s
```

## Outputs

The output directory contains:

| File | Contents |
| --- | --- |
| `lineage_nodes.csv.gz` | Event-resolved branches, original PhysiCell IDs, parent nodes, birth/end times, and branch lengths. |
| `terminal_cells.csv.gz` | PhysiCell IDs and node IDs actually emitted as observations. |
| `physicell_lineage_full.nwk` | Complete event tree implied by the division log. |
| `physicell_lineage_sampled.nwk` | Tree restricted to observed/live terminal IDs; unary nodes are retained to preserve elapsed time. |
| `mutation_events.csv.gz` | Newly acquired edits with exact event time, branch, integration, position, allele, time interval, and editing-state provenance. |
| `mutation_event_descendant_matrix_sparse.rds` | Literal sampled-cell by unique mutation-event matrix. Every mutation-log row has its own column, with `1` for each sampled terminal descendant of the event's branch. |
| `mutation_event_descendant_manifest.csv.gz` | One-to-one mapping from event-matrix columns to modality, modality-specific event-row number, descendant count, and the complete original mutation-event fields. |
| `barcode_profiles.rds` | Named per-cell numeric mutation matrices using the repository's `0/-1/1..4/decimal` encoding. |
| `barcode_alleles.csv.gz` | Flattened integration-by-position allele matrix. |
| `barcode_binary_score_matrix.csv.gz` | The same observations encoded as WT `0` versus edited `1`. |
| `barcode_alleles_sparse.rds` | Compact sparse allele matrix when `physicell_adapter.compact_output` is enabled. |
| `barcode_binary_score_matrix_sparse.rds` | Compact sparse WT/edited matrix when compact output is enabled. |
| `barcode_target_layout.csv.gz` | Primary base-editing targets, integration indices, rate classes, and reference/alternate bases. |
| `barcode_reference.fasta` | Generated or configured barcode reference. |
| `barcode_sequences.fasta` | One reconstructed sequence per terminal cell and integration. |
| `run_manifest.csv.gz` | Counts and configuration summary. |
| `mitochondrial_profiles.rds` | Sparse per-cell list of inherited mitochondrial genome variants. |
| `mitochondrial_mutation_events.csv.gz` | Founder and branch mt substitutions/indels with exact event time, genome, position, branch, and rate-state provenance. |
| `mitochondrial_variant_fractions.csv.gz` | Long-form terminal-cell heteroplasmy fractions by position and allele. |
| `mitochondrial_variant_fraction_matrix.rds` | Optional sparse cell-by-variant heteroplasmy matrix, written when `Matrix` is installed. |
| `mitochondrial_reference.fasta` | Generated mitochondrial reference with `mito_genome_length` bases. |
| `mitochondrial_sampled_haplotypes.fasta` | One sampled mt genome per terminal cell when `--write-mt-fasta true`. |
| `mitochondrial_manifest.csv.gz` | Mitochondrial dimensions, observation count, state, and seed. |
| `combined_lineage_feature_matrix.rds` | Aligned sparse barcode-plus-mitochondrial feature matrix when both modalities are requested. |
| `combined_lineage_feature_manifest.csv.gz` | Feature modality and recorder-system annotations for the combined matrix. |
| `r_timing_summary.csv.gz` | Wall-clock seconds, readable duration, and percentage of total R runtime for each executed phase and the complete simulator. |

The dense barcode CSV and per-cell barcode FASTA are intentionally omitted
when compact output is enabled. `profile_storage: "sparse"` additionally
avoids retaining an entire dense barcode for every lineage node, and
`retain_internal_profiles: false` releases completed internal profiles during
the replay. These settings are used by the day-21 BASELINE-like organoid
preset.

All R-generated CSV tables are compressed by default. Pass
`--compress-csv false` to either R CLI to retain legacy uncompressed `.csv`
names. Downstream PhysiCell/scDesign3 readers accept either form and prefer
`.csv.gz` when both exist.

The event-descendant matrix is a ground-truth ancestry matrix, rather than an
observed-allele matrix. Its rows are the sampled terminal cells in lineage
depth-first order and should be aligned to other matrices by row name. An
event on a branch with no sampled descendants retains its manifest/matrix
column with all zeros. For mitochondrial events, `1` means that the cell is a
topological descendant of the event-bearing branch; it does not assert that
the mutated mitochondrial genome survived every stochastic daughter-cell
bottleneck.

## Recording model

The barcode replay uses:

- nucleotide composition or `barcode_sequence`;
- BE and nuclease target layouts, H/M/L classes, editing windows, and
  close-after-edit behavior;
- JC/K80/K81/F81/HKY/GTR barcode substitution parameters;
- background substitution/insertion/deletion probabilities;
- BE substitution and nuclease insertion/deletion probabilities;
- induced and uninduced parameter blocks.

Per-division probabilities are converted to continuous-time hazards using the
selected cell type's `cell_cycle_length` and cached with the prepared rate set.
For ordinary non-closing editing windows, all unedited coordinates in an
integration are drawn and updated as vectors; close-after-edit configurations
retain a sequential fallback so edits acquired earlier in the segment can
close later coordinates. Competing edit types are sampled on each branch for
its actual PhysiCell duration. Recording sites are irreversible after their
first edit, and both children inherit the complete parental profile at every
division. Mutation-event columns grow geometrically rather than being
row-bound per node, and compact barcode matrices use a count-then-fill pass
with exact sparse-triplet preallocation.

Replay remains reproducible for a fixed seed under this implementation.
Because vectorized draws consume the random-number stream in a different
order, a seed does not reproduce the exact event realization emitted by older
scalar-loop versions.

The mitochondrial replay uses `mito_genome_length`, baseline heteroplasmy
fraction/Beta-distribution/transition probability, `mt_invariant_sites`,
HKY/other substitution parameters, and mitochondrial background indel
probabilities. Each cell holds a configurable fixed number of sparse mt genome
profiles. Both division products independently resample the parental genomes
with replacement, providing a genetic bottleneck and heteroplasmy drift;
low-rate branch events are drawn from continuous hazards.

The random seed controls reference/barcode construction, initial
heteroplasmy, target-class rate draws, mitochondrial resampling, and branch
mutations.

## Current boundaries

- Mitochondrial replay models mutation, a fixed-size division bottleneck, and
  heteroplasmy drift. It does not replay organelle-level fusion/fission,
  variable copy number, mitochondrial dropout, or fitness-dependent cell
  survival from the native remote_mito population simulator.
- The division CSV contains no cell type. One selected remote_mito cell-type
  rate block is therefore applied to all branches.
- An explicit `--founders` CSV preserves founders that never appear as a parent
  in the division log. `--founder-label-sites` can assign stable, allele-coded
  founder barcodes before branch recording.
- `editing_state=auto` models one global induction time. The JSON's
  per-cell induction count/fraction and differentiation transition matrices
  cannot be reconstructed without corresponding per-cell PhysiCell metadata.
- Targeted insertions are emitted as ordinary one-base insertions. When
  `prime_editing_system` is enabled, the adapter warns because guide-specific
  insertion sequences are not replayed.
- A current-cell table identifies final retained cells but does not provide death
  times. Extinct tips are omitted from sampled data; their exact terminal
  branch lengths cannot be recovered from the current PhysiCell files.

The small files under `tests/fixtures/physicell_*` provide a runnable example
of the accepted formats.
