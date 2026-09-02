# Neural organoid simulation

`run_physicell_organoid_pipeline.sh` provides a first end-to-end cortical
organoid simulation that starts from 2,500 individually tracked iPSCs. It
stages a build from the lineage-enabled PhysiCell checkout, overlays the
organoid-specific model in `physicell_projects/neural_organoid_lineage`, and
does not modify `../PhysiCell`.

## Quick start

From this repository:

```bash
bash run_physicell_organoid_pipeline.sh
```

The default run:

- seeds 2,500 iPSCs in a packed 150-micron-radius spheroid;
- uses the `early21` growth preset and simulates up to 21 culture days or a
  150,000-cell safety threshold;
- targets roughly 25% neurons at day 21 while keeping the population below the
  150,000-cell safety threshold;
- tracks seven states: `iPSC`, `neuroepithelial`, `radial_glia`,
  `neural_progenitor` (NPC), `neuron`, `astrocyte`, and `epithelial`, although
  astrocytes are not expected before the current day-45 gliogenesis gate;
- models time-gated stochastic state transitions and radial-glia asymmetric
  divisions that retain the parent while producing NPC/late astrocyte daughters;
- sends a configurable minority of induced iPSCs directly to a non-proliferating
  epithelial state, with a default hard cap of 100 epithelial commitments;
- uses effective cycle lengths of 24 hours for iPSCs and 8 days for both radial
  glia and NPCs;
- uses an NPC-to-neuron transition hazard of `4.5e-5 /min` after day 6
  (a mean waiting time of approximately 15.4 days);
- exposes cells to diffusing oxygen and nutrient fields;
- writes division, founder, transition, and terminal-state tables;
- replays base-editing and mitochondrial recording; and
- reserves 12 allele-coded barcode positions to distinguish the 2,500
  day-zero founders.

The full run can be computationally expensive. A useful installation and
interface check is:

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

The output directory must be new or empty.

### Spatial PhysiCell versus Gillespie

The organoid wrapper above should remain the primary model when oxygen,
nutrient gradients, mechanics, spatial neighborhoods, or the staged
seven-state developmental program affect cell behavior. The Gillespie engine
is a complementary well-mixed model for large exact continuous-time lineage
experiments. It uses the native JSON cell-cycle, death, induction, and daughter
transition parameters, then applies the same BASELINE-like and mitochondrial
recorders to its event tree.

For a calibrated non-spatial approximation, add `"simulation_engine":
"gillespie"` and a `gillespie` block to a native parameter JSON, then run
`Rscript sim5_code.R -P <parameters.json>`. Its
`gillespie_cell_states.csv.gz` and lineage-recording outputs can be passed to
`generate_physicell_sc_profiles.R --run-dir <gillespie-output> ...`. This does
not turn the spatial PhysiCell mechanics themselves into a Gillespie process.

The earlier conservative 90-day behavior remains available explicitly:

```bash
bash run_physicell_organoid_pipeline.sh \
  --growth-preset long90 \
  --output-dir output/organoid_day90
```

`early21` defaults to 21 days and a 150,000-cell stop. `long90` defaults to 90
days, a 100,000-cell stop, a 12-day radial-glia cycle, and a 6-day NPC cycle.
Explicit `--days` and `--target-cells` values override the preset duration and
safety stop.

The epithelial side population is controlled independently:

```bash
bash run_physicell_organoid_pipeline.sh \
  --epithelial-probability 0.02 \
  --max-epithelial-cells 100
```

The probability is evaluated when an iPSC undergoes neural induction. Once the
cap is reached, subsequent induced iPSCs continue into the neuroepithelial
lineage. Epithelial cells are terminal in the current model. Set the cap to
zero to disable this branch.

## BASELINE-like Cas12a and mitochondrial recording

Use the dedicated wrapper to run the organoid model and replay two independent
recording systems on the exact PhysiCell division tree:

```bash
bash run_physicell_organoid_baseline_pipeline.sh \
  --output-dir output/organoid_baseline_day21
```

For a large lineage, event-tree reconstruction and each recording modality
print approximately 20 progress updates, including percent complete, elapsed
time, throughput, and ETA. These messages are also saved in
`recording_run.log`. Run `simulate_physicell_lineage.R` with
`--progress-updates N` to change the frequency or `--progress false` to
disable them. Timestamped stage logs announce when each lineage simulation
finishes, which output artifact group is being generated, and when combined
output is complete.

The preset is modeled after
[BASELINE](https://mckennalab.org/papers/baseline-crispr-lineage-tracing/):
five independently inherited hyperdCas12a-ABE8e recorder integrations, each
with 50 primary targets (250 integrated target sites total), irreversible
A-to-G edits, and recording active from time zero. It uses 20 high-rate and 30
low-rate targets per integration, a seven-base decaying editing window, and a
mean induced edit probability of `0.00952` per target per 24-hour iPSC cycle.
These are simulation defaults, not a fitted estimate for this organoid system.

The wrapper also adds 12 founder-label positions per integration and tracks 32
mitochondrial genomes per cell. Founder labels identify the 2,500 starting
iPSCs; dynamic Cas12a edits and mitochondrial variants are acquired at
continuous event times along each PhysiCell branch and inherited by both
division products. Mitochondrial genomes are independently bottleneck-sampled
by each child, allowing heteroplasmy drift on the same tree.

### ecDNA plus BASELINE and mitochondrial recording

The ecDNA wrapper enables all three recording modalities using the
`ecdna_adapter` block in
`example_json_params/physicell_neural_organoid_baseline.json`:

```bash
bash run_physicell_organoid_ecdna_pipeline.sh \
  --output-dir output/organoid_baseline_ecdna_day21
```

It starts each founder with ten ecDNA species at five copies/species, labels an
exact 30% of species with a 12-base static ID and six-target CRISPR recorder,
and uses unbiased joint daughter segregation. These are demonstrative defaults
that should be calibrated to an ecDNA system. The ordinary BASELINE wrapper
still defaults to barcode plus mitochondrial recording; pass `--modalities
all` to any organoid wrapper to enable ecDNA explicitly.

The preset's top-level `non_mendelian_selection` block is neutral by default.
Set `ecdna_label_coefficient`, `ecdna_recorder_edit_coefficient`, or
`mitochondrial_variant_coefficient` between zero and one to penalize propagation
of the corresponding marker. `coefficient` supplies a shared fallback. This
selection changes intracellular copy inheritance on the already simulated
PhysiCell tree; it does not change cell birth, death, or differentiation rates.

The BASELINE preset uses sparse in-memory profiles and compact output so a
full day-21 run does not create a dense
`cells x integrations x 1,435 bases` object. The key recording files are:

| File | Contents |
| --- | --- |
| `barcode_target_layout.csv.gz` | The 250 primary target/integration combinations and their rate classes. |
| `barcode_alleles_sparse.rds` | Sparse terminal-cell matrix of edited nucleotide positions and founder labels. |
| `barcode_binary_score_matrix_sparse.rds` | Sparse WT/edited version of the barcode matrix. |
| `mutation_events.csv.gz` | Cas12a and founder-label events with exact event times and branch provenance. |
| `mutation_event_descendant_matrix_sparse.rds` | Literal sampled-cell by unique barcode/mitochondrial event ancestry matrix, with one column per event-log row. |
| `mutation_event_descendant_manifest.csv.gz` | Matrix-column mapping to modality, original event row, branch, mutation details, and sampled-descendant count. |
| `mitochondrial_variant_fraction_matrix.rds` | Sparse terminal-cell mitochondrial heteroplasmy matrix. |
| `mitochondrial_mutation_events.csv.gz` | Mitochondrial events with exact event times and genome IDs. |
| `combined_lineage_feature_matrix.rds` | Aligned Cas12a plus mitochondrial features; ecDNA recorder/static-ID features are appended when `--modalities all` is used. |
| `combined_lineage_feature_manifest.csv.gz` | Feature-to-modality and recorder-system mapping. |

One 50-target integration can contribute more than 50 matrix columns because
the seven-base editing windows can contain several editable adenines. Founder
label columns are also included; `barcode_target_layout.csv.gz` is the
authoritative table of the 250 primary integrated targets.

This is BASELINE-like rather than a full experimental recovery model. It
simulates editing and inheritance but does not yet model integration dropout,
amplicon capture, sequencing errors, or target-specific rates fitted from
empirical BASELINE observations.

The literal event-descendant matrix describes ground-truth tree ancestry.
Barcode event columns also describe irreversible inherited edits.
Mitochondrial event columns deliberately mark every topological descendant,
even though bottleneck sampling can remove that mitochondrial allele from an
observed descendant.

## scDesign3 expression

Supply a longitudinal organoid reference containing a raw count assay,
cell-type labels, and developmental pseudotime:

```bash
bash run_physicell_organoid_pipeline.sh \
  --sc-reference data/cortical_organoid_reference.rds \
  --sc-celltype-col cell_type \
  --sc-pseudotime-col pseudotime \
  --sc-celltype-map config/organoid_celltype_map.json
```

The model always supplies per-cell developmental pseudotime to the scDesign3
stage. It is computed from the current PhysiCell fate and time spent in that
state, rather than from lineage depth. Optional linked predictors include
`culture_day`, `time_in_state`, `transition_count`, `oxygen`, `nutrient`,
`radial_position`, `neighbor_count`, `barcode_edit_fraction`, and
`mt_heteroplasmy_burden`. A predictor requested with
`--sc-other-covariates` must also exist under the same name in the reference
`colData`.

For example:

```bash
bash run_physicell_organoid_pipeline.sh \
  --sc-reference data/cortical_organoid_reference.rds \
  --sc-celltype-col cell_type \
  --sc-pseudotime-col pseudotime \
  --sc-other-covariates culture_day,oxygen
```

If the reference uses different labels, provide a mapping such as:

```json
{
  "iPSC": "pluripotent",
  "neuroepithelial": "neuroepithelium",
  "radial_glia": "RG",
  "neural_progenitor": "NPC",
  "neuron": "excitatory_neuron",
  "astrocyte": "astroglia",
  "epithelial": "epithelial"
}
```

## PhysiCell model and outputs

The model source is:

- `physicell_projects/neural_organoid_lineage/config/PhysiCell_settings.xml`;
- `physicell_projects/neural_organoid_lineage/custom_modules/custom.cpp`; and
- `physicell_projects/neural_organoid_lineage/custom_modules/custom.h`.

The XML contains all initial dimensions, cycle rates, transition start times,
transition hazards, hypoxia response, and substrate settings. These are
working simulation defaults, not fitted biological estimates. The radial-glia
division probability defaults to fully asymmetric after neurogenesis begins;
change the corresponding XML user parameter to model symmetric self-renewal.

The day-21 population scale is anchored to a
[reported expansion](https://www.nature.com/articles/s41598-019-48347-2)
from 9,000 plated hiPSCs to approximately 180,000 cells over about 20 days.
Scaling the same 20-fold expansion to 2,500 founders gives 50,000 cells. The
default preset uses a faster 24-hour iPSC cycle and is additionally calibrated
to a day-21 neuron fraction near 25%. The full 2,500-founder seed-2 validation
produced 104,804 cells: 25,986 neurons (24.79%), 60,213 NPCs (57.45%), 18,499
radial glia (17.65%), 99 surviving epithelial cells (0.094%), and 7 residual
iPSC/neuroepithelial cells. The epithelial branch made exactly 100 commitments.
Of all cells, 104,801 were alive; the maximum radius was 669 microns, and the
oxygen range was 15.3--27.4 mmHg. The resulting population and composition
remain order-of-magnitude calibration results, not an assertion that growth is
linear across protocols, cell lines, or organoid sizes.

PhysiCell writes:

| File | Contents |
| --- | --- |
| `cell_lineage.csv` | Persistent-parent division events. |
| `founder_cells.csv` | Every day-zero iPSC ID, including founders that never divide. |
| `cell_state_transitions.csv` | Time, cell/founder IDs, old/new state, position, oxygen, and nutrient at each transition. |
| `organoid_cell_states.csv` | Final live/dead state, founder, fate, position, neighborhood, state age, pseudotime, and niche values. |

The lineage-recording and scDesign3 outputs use the same
`cell_<PhysiCell_ID>` sample identifiers.

## Founder labels

`--founder-label-sites N` encodes each founder with two possible non-reference
alleles at `N` positions outside the configured recording targets/windows. At least
`ceiling(log2(number_of_founders))` sites are required when founder labeling is
enabled, so 12 sites distinguish up to 4,096 founders. The sites are initialized
before branch recording and are irreversible.

Founder identity is represented in `barcode_alleles.csv` and the nucleotide
FASTA. The binary WT-versus-edited score matrix intentionally discards allele
identity and therefore cannot distinguish these founder codes by itself.

## Current model boundaries

This first implementation is a runnable integration scaffold. Before treating
its results as biological predictions, it needs calibration against a selected
organoid protocol and longitudinal reference.

In particular:

- Fate changes are time-gated stochastic transitions, not a fitted gene
  regulatory network.
- The epithelial branch is a capped terminal off-target fate, not a mechanistic
  epithelial differentiation program.
- The model does not yet represent neuroepithelial rosettes, apical polarity,
  lumen formation, radial fibers, cortical layers, electrical activity, or
  region-specific morphogen patterning.
- The final state table is linked to scDesign3, but recording replay currently
  uses one configured remote_mito rate block for every branch. State-specific
  recorder and mitochondrial kinetics are a planned extension.
- The transition log records fate changes, but the lineage adapter does not yet
  split a branch at those transition times.
- Founder mitochondrial profiles are initialized independently. A shared
  pre-aggregation iPSC expansion tree is not yet simulated.
- scDesign3 conditions expression on state, pseudotime, spatial, and optional
  niche covariates. It does not add an unobserved sibling/clone expression
  covariance beyond those predictors.
- Exact death times are not included in the imported tree. Dead cells are
  filtered from final observations through the `alive` field.

Recommended calibration targets are cell counts over time, state proportions,
cycle-length distributions, organoid radius, hypoxic-core size, transition
timing, and gene-expression trajectories.
