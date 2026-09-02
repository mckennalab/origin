#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage:
  bash run_physicell_organoid_pipeline.sh [options]

Stages a lineage-enabled PhysiCell neural-organoid model, seeds individually
tracked iPSCs, simulates cortical differentiation, replays barcode and
mitochondrial/ecDNA recording, and optionally generates scDesign3 transcriptomes.

Options:
  --physicell-dir PATH       PhysiCell checkout [default: ../PhysiCell]
  --params PATH              Recording JSON [default: example_json_params/physicell_neural_organoid.json]
  --output-dir PATH          Run directory [default: output/physicell_organoid_<timestamp>]
  --growth-preset NAME       early21 or long90 [default: early21]
  --initial-cells N          Day-zero iPSC founders [default: 2500]
  --initial-radius MICRONS   Initial packed-spheroid radius [default: 150]
  --target-cells N           Safety stop [preset default: 150000 or 100000]
  --days N                   Duration [preset default: 21 or 90 days]
  --epithelial-probability P Chance an induced iPSC becomes epithelial [default: 0.02]
  --max-epithelial-cells N   Hard cap on direct iPSC-derived epithelial cells [default: 100]
  --founder-label-sites N    Stable founder barcode sites [default: 12]
  --num-integrations N       Barcodes per cell [default: 1]
  --mt-genomes-per-cell N    Fixed mt bottleneck size [default: 32]
  --modalities LIST          barcode, mitochondrial, ecDNA, both, or all [default: both]
  --write-mt-fasta BOOL      Write one mt haplotype per cell [default: false]
  --seed N                   PhysiCell and recording seed [default: 1]
  --jobs N                   Parallel compiler/simulation workers
  --sc-reference PATH        Optional SCE/Seurat RDS or h5ad
  --sc-celltype-col NAME     Reference cell-type column [default: cell_type]
  --sc-celltype-map PATH     Optional simulated-to-reference type map
  --sc-pseudotime-col NAME   Reference pseudotime column [default: pseudotime]
  --sc-spatial-cols X,Y      Optional reference spatial-coordinate columns
  --sc-other-covariates LIST Additional matched covariates
  --sc-mu-formula FORMULA    Optional explicit scDesign3 mean formula
  --sc-cache-dir PATH        Optional shared scDesign3 fit cache
  --sc-ncores N              scDesign3 workers [default: --jobs]
  -h, --help

The PhysiCell checkout is never modified. Model assumptions in this first MVP
are intentionally configurable and require calibration to a chosen organoid
protocol/reference before biological interpretation.
EOF
}

physicell_dir="$script_dir/../PhysiCell"
params_path="$script_dir/example_json_params/physicell_neural_organoid.json"
model_overlay="$script_dir/physicell_projects/neural_organoid_lineage"
output_dir=""
growth_preset=early21
initial_cells=2500
initial_radius=150
target_cells=""
culture_days=""
epithelial_probability=0.02
max_epithelial_cells=100
founder_label_sites=12
num_integrations=1
mt_genomes_per_cell=32
recording_modalities=both
write_mt_fasta=false
seed=1
sc_reference=""
sc_celltype_col=cell_type
sc_celltype_map=""
sc_pseudotime_col=pseudotime
sc_spatial_cols=""
sc_other_covariates=""
sc_mu_formula=""
sc_cache_dir=""
sc_ncores=""
if command -v sysctl >/dev/null 2>&1; then
  jobs="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
elif command -v nproc >/dev/null 2>&1; then
  jobs="$(nproc)"
else
  jobs=4
fi

while (($# > 0)); do
  case "$1" in
    --physicell-dir) physicell_dir="$2"; shift 2 ;;
    --params) params_path="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --growth-preset) growth_preset="$2"; shift 2 ;;
    --initial-cells) initial_cells="$2"; shift 2 ;;
    --initial-radius) initial_radius="$2"; shift 2 ;;
    --target-cells) target_cells="$2"; shift 2 ;;
    --days) culture_days="$2"; shift 2 ;;
    --epithelial-probability) epithelial_probability="$2"; shift 2 ;;
    --max-epithelial-cells) max_epithelial_cells="$2"; shift 2 ;;
    --founder-label-sites) founder_label_sites="$2"; shift 2 ;;
    --num-integrations) num_integrations="$2"; shift 2 ;;
    --mt-genomes-per-cell) mt_genomes_per_cell="$2"; shift 2 ;;
    --modalities) recording_modalities="$2"; shift 2 ;;
    --write-mt-fasta)
      write_mt_fasta="$(printf '%s' "$2" | tr '[:upper:]' '[:lower:]')"
      shift 2
      ;;
    --seed) seed="$2"; shift 2 ;;
    --jobs) jobs="$2"; shift 2 ;;
    --sc-reference) sc_reference="$2"; shift 2 ;;
    --sc-celltype-col) sc_celltype_col="$2"; shift 2 ;;
    --sc-celltype-map) sc_celltype_map="$2"; shift 2 ;;
    --sc-pseudotime-col) sc_pseudotime_col="$2"; shift 2 ;;
    --sc-spatial-cols) sc_spatial_cols="$2"; shift 2 ;;
    --sc-other-covariates) sc_other_covariates="$2"; shift 2 ;;
    --sc-mu-formula) sc_mu_formula="$2"; shift 2 ;;
    --sc-cache-dir) sc_cache_dir="$2"; shift 2 ;;
    --sc-ncores) sc_ncores="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      printf 'Unknown option: %s\n\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

ipsc_cycle_rate=0.000694444
case "$growth_preset" in
  early21)
    culture_days="${culture_days:-21}"
    target_cells="${target_cells:-150000}"
    radial_glia_cycle_rate=0.0000868056
    neural_progenitor_cycle_rate=0.0000868056
    neural_progenitor_to_neuron_rate=0.000045
    ;;
  long90)
    culture_days="${culture_days:-90}"
    target_cells="${target_cells:-100000}"
    radial_glia_cycle_rate=0.0000578704
    neural_progenitor_cycle_rate=0.000115741
    neural_progenitor_to_neuron_rate=0.00035
    ;;
  *)
    printf -- '--growth-preset must be early21 or long90.\n' >&2
    exit 2
    ;;
esac

if [[ -z "$sc_ncores" ]]; then
  sc_ncores="$jobs"
fi
for integer_value in \
  "$initial_cells" "$target_cells" "$culture_days" "$founder_label_sites" \
  "$num_integrations" "$mt_genomes_per_cell" "$max_epithelial_cells" \
  "$seed" "$jobs" "$sc_ncores"; do
  if [[ ! "$integer_value" =~ ^[0-9]+$ ]]; then
    printf 'Expected a non-negative integer, received: %s\n' "$integer_value" >&2
    exit 2
  fi
done
if [[ ! "$initial_radius" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  printf -- '--initial-radius must be a positive number.\n' >&2
  exit 2
fi
if [[ ! "$epithelial_probability" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   ! awk -v value="$epithelial_probability" \
     'BEGIN { exit !(value >= 0 && value <= 1) }'; then
  printf -- '--epithelial-probability must be between 0 and 1.\n' >&2
  exit 2
fi
if ((initial_cells < 1 || target_cells < initial_cells || culture_days < 1 ||
     num_integrations < 1 || mt_genomes_per_cell < 1 || jobs < 1 ||
     sc_ncores < 1)); then
  printf 'Counts, culture days, and workers must be positive; target cells must not be below founders.\n' >&2
  exit 2
fi
if [[ "$write_mt_fasta" != "true" && "$write_mt_fasta" != "false" ]]; then
  printf -- '--write-mt-fasta must be true or false.\n' >&2
  exit 2
fi

physicell_dir="$(cd "$physicell_dir" && pwd)"
params_path="$(cd "$(dirname "$params_path")" && pwd)/$(basename "$params_path")"
project_source="$physicell_dir/user_projects/tumor_3D_lineage"
for required_path in \
  "$physicell_dir/core" "$physicell_dir/modules" "$physicell_dir/BioFVM" \
  "$project_source/main.cpp" "$project_source/Makefile" \
  "$model_overlay/custom_modules/custom.cpp" \
  "$model_overlay/custom_modules/custom.h" \
  "$model_overlay/config/PhysiCell_settings.xml" "$params_path"; do
  if [[ ! -e "$required_path" ]]; then
    printf 'Required input not found: %s\n' "$required_path" >&2
    exit 1
  fi
done
for executable in make g++ Rscript perl awk; do
  if ! command -v "$executable" >/dev/null 2>&1; then
    printf 'Required executable not found: %s\n' "$executable" >&2
    exit 1
  fi
done

if [[ -n "$sc_reference" ]]; then
  sc_reference="$(cd "$(dirname "$sc_reference")" && pwd)/$(basename "$sc_reference")"
  if [[ ! -f "$sc_reference" ]]; then
    printf 'scDesign3 reference not found: %s\n' "$sc_reference" >&2
    exit 1
  fi
fi
if [[ -n "$sc_celltype_map" ]]; then
  sc_celltype_map="$(cd "$(dirname "$sc_celltype_map")" && pwd)/$(basename "$sc_celltype_map")"
  if [[ ! -f "$sc_celltype_map" ]]; then
    printf 'scDesign3 cell-type map not found: %s\n' "$sc_celltype_map" >&2
    exit 1
  fi
fi

if [[ -z "$output_dir" ]]; then
  output_dir="$script_dir/output/physicell_organoid_$(date '+%Y%m%d_%H%M%S')"
elif [[ "$output_dir" != /* ]]; then
  output_dir="$script_dir/$output_dir"
fi
if [[ -e "$output_dir" ]] &&
   [[ -n "$(find "$output_dir" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  printf 'Output directory already exists and is not empty: %s\n' "$output_dir" >&2
  exit 1
fi

build_dir="$output_dir/physicell_build"
physicell_output="$build_dir/output"
recording_output="$output_dir/lineage_recording"
mkdir -p "$build_dir/config" "$build_dir/custom_modules" "$recording_output"

cp "$project_source/main.cpp" "$build_dir/main.cpp"
cp "$project_source/Makefile" "$build_dir/Makefile"
if [[ -f "$project_source/VERSION.txt" ]]; then
  cp "$project_source/VERSION.txt" "$build_dir/VERSION.txt"
fi
cp "$model_overlay/custom_modules/custom.cpp" "$build_dir/custom_modules/custom.cpp"
cp "$model_overlay/custom_modules/custom.h" "$build_dir/custom_modules/custom.h"
cp "$model_overlay/config/PhysiCell_settings.xml" "$build_dir/config/PhysiCell_settings.xml"
ln -s "$physicell_dir/core" "$build_dir/core"
ln -s "$physicell_dir/modules" "$build_dir/modules"
ln -s "$physicell_dir/BioFVM" "$build_dir/BioFVM"

ORGANOID_TARGET_CELLS="$target_cells" perl -0pi -e \
  's/int target_number_of_cells = [0-9]+;/int target_number_of_cells = $ENV{ORGANOID_TARGET_CELLS};/' \
  "$build_dir/main.cpp"
perl -0pi -e '
  s{(\tsave_cell_lineage_table\( filename \);)}{$1\n\n\tsave_organoid_outputs( PhysiCell_settings.folder );};
  s{(\t// write the lineage outputs:)}{\tstd::cout << "Finished organoid simulation at t = "\n\t\t<< PhysiCell_globals.current_time << " min with "\n\t\t<< (*all_cells).size() << " cells." << std::endl;\n\n$1};
' "$build_dir/main.cpp"
if ! grep -q "save_organoid_outputs" "$build_dir/main.cpp"; then
  printf 'Could not add organoid output calls to the staged PhysiCell main.\n' >&2
  exit 1
fi

max_time=$((culture_days * 1440))
ORGANOID_INITIAL_CELLS="$initial_cells" \
ORGANOID_INITIAL_RADIUS="$initial_radius" \
ORGANOID_EPITHELIAL_PROBABILITY="$epithelial_probability" \
ORGANOID_MAX_EPITHELIAL_CELLS="$max_epithelial_cells" \
ORGANOID_MAX_TIME="$max_time" \
ORGANOID_SEED="$seed" \
ORGANOID_THREADS="$jobs" \
ORGANOID_IPSC_CYCLE_RATE="$ipsc_cycle_rate" \
ORGANOID_RADIAL_GLIA_CYCLE_RATE="$radial_glia_cycle_rate" \
ORGANOID_NPC_CYCLE_RATE="$neural_progenitor_cycle_rate" \
ORGANOID_NPC_TO_NEURON_RATE="$neural_progenitor_to_neuron_rate" perl -0pi -e '
  s{<max_time units="min">[^<]*</max_time>}{<max_time units="min">$ENV{ORGANOID_MAX_TIME}</max_time>};
  s{<omp_num_threads>[^<]*</omp_num_threads>}{<omp_num_threads>$ENV{ORGANOID_THREADS}</omp_num_threads>};
  s{<random_seed([^>]*)>[^<]*</random_seed>}{<random_seed$1>$ENV{ORGANOID_SEED}</random_seed>}g;
  s{<initial_cells([^>]*)>[^<]*</initial_cells>}{<initial_cells$1>$ENV{ORGANOID_INITIAL_CELLS}</initial_cells>};
  s{<initial_organoid_radius([^>]*)>[^<]*</initial_organoid_radius>}{<initial_organoid_radius$1>$ENV{ORGANOID_INITIAL_RADIUS}</initial_organoid_radius>};
  s{<ipsc_epithelial_probability([^>]*)>[^<]*</ipsc_epithelial_probability>}{<ipsc_epithelial_probability$1>$ENV{ORGANOID_EPITHELIAL_PROBABILITY}</ipsc_epithelial_probability>};
  s{<max_epithelial_cells([^>]*)>[^<]*</max_epithelial_cells>}{<max_epithelial_cells$1>$ENV{ORGANOID_MAX_EPITHELIAL_CELLS}</max_epithelial_cells>};
  s{(<cell_definition name="iPSC".*?<phase_transition_rates[^>]*>.*?<rate[^>]*>)[^<]*}{$1$ENV{ORGANOID_IPSC_CYCLE_RATE}}s;
  s{(<cell_definition name="radial_glia".*?<phase_transition_rates[^>]*>.*?<rate[^>]*>)[^<]*}{$1$ENV{ORGANOID_RADIAL_GLIA_CYCLE_RATE}}s;
  s{(<cell_definition name="neural_progenitor".*?<phase_transition_rates[^>]*>.*?<rate[^>]*>)[^<]*}{$1$ENV{ORGANOID_NPC_CYCLE_RATE}}s;
  s{<ipsc_cycle_rate([^>]*)>[^<]*</ipsc_cycle_rate>}{<ipsc_cycle_rate$1>$ENV{ORGANOID_IPSC_CYCLE_RATE}</ipsc_cycle_rate>};
  s{<radial_glia_cycle_rate([^>]*)>[^<]*</radial_glia_cycle_rate>}{<radial_glia_cycle_rate$1>$ENV{ORGANOID_RADIAL_GLIA_CYCLE_RATE}</radial_glia_cycle_rate>};
  s{<neural_progenitor_cycle_rate([^>]*)>[^<]*</neural_progenitor_cycle_rate>}{<neural_progenitor_cycle_rate$1>$ENV{ORGANOID_NPC_CYCLE_RATE}</neural_progenitor_cycle_rate>};
  s{<neural_progenitor_to_neuron_rate([^>]*)>[^<]*</neural_progenitor_to_neuron_rate>}{<neural_progenitor_to_neuron_rate$1>$ENV{ORGANOID_NPC_TO_NEURON_RATE}</neural_progenitor_to_neuron_rate>};
' "$build_dir/config/PhysiCell_settings.xml"

cp "$params_path" "$output_dir/recording_params.json"

printf 'Building staged PhysiCell neural-organoid model with %s job(s)...\n' "$jobs"
(
  cd "$build_dir"
  make -j "$jobs"
) 2>&1 | tee "$output_dir/physicell_build.log"

printf 'Simulating %s iPSC founders for up to %s day(s)...\n' \
  "$initial_cells" "$culture_days"
(
  cd "$build_dir"
  ./project ./config/PhysiCell_settings.xml
) 2>&1 | tee "$output_dir/physicell_run.log"

lineage_path="$physicell_output/cell_lineage.csv"
state_path="$physicell_output/organoid_cell_states.csv"
founder_path="$physicell_output/founder_cells.csv"
transition_path="$physicell_output/cell_state_transitions.csv"
for expected_output in \
  "$lineage_path" "$state_path" "$founder_path" "$transition_path"; do
  if [[ ! -s "$expected_output" ]]; then
    printf 'PhysiCell did not produce expected output: %s\n' "$expected_output" >&2
    exit 1
  fi
done

current_cell_count="$(
  awk -F, 'NR > 1 && $17 == "true" { count++ } END { print count + 0 }' \
    "$state_path"
)"
division_count="$(
  awk 'NR > 1 { count++ } END { print count + 0 }' "$lineage_path"
)"
end_time="$(
  sed -nE 's/.*Finished organoid simulation at t = ([0-9.eE+-]+) min.*/\1/p' \
    "$output_dir/physicell_run.log" | tail -n 1
)"
if [[ -z "$end_time" ]]; then
  printf 'Could not recover final PhysiCell time from its run log.\n' >&2
  exit 1
fi

printf 'Replaying recording on %s divisions from %s founders...\n' \
  "$division_count" "$initial_cells"
Rscript "$script_dir/origin/inst/scripts/simulate_physicell_lineage.R" \
  --lineage "$lineage_path" \
  --founders "$founder_path" \
  --live-cells "$state_path" \
  --params "$params_path" \
  --output-dir "$recording_output" \
  --end-time "$end_time" \
  --editing-state auto \
  --modalities "$recording_modalities" \
  --num-integrations "$num_integrations" \
  --founder-label-sites "$founder_label_sites" \
  --mt-genomes-per-cell "$mt_genomes_per_cell" \
  --write-mt-fasta "$write_mt_fasta" \
  --seed "$seed" 2>&1 | tee "$output_dir/recording_run.log"

if [[ -n "$sc_reference" ]]; then
  printf 'Generating cell-state-linked scDesign3 transcriptomes...\n'
  scdesign3_args=(
    --run-dir "$output_dir"
    --reference "$sc_reference"
    --lineage-table "$state_path"
    --celltype-col "$sc_celltype_col"
    --use-pseudotime
    --pseudotime-col "$sc_pseudotime_col"
    --ncores "$sc_ncores"
    --seed "$seed"
  )
  if [[ -n "$sc_celltype_map" ]]; then
    scdesign3_args+=(--celltype-map "$sc_celltype_map")
  fi
  if [[ -n "$sc_spatial_cols" ]]; then
    scdesign3_args+=(--reference-spatial-cols "$sc_spatial_cols")
  fi
  if [[ -n "$sc_other_covariates" ]]; then
    scdesign3_args+=(--other-covariates "$sc_other_covariates")
  fi
  if [[ -n "$sc_mu_formula" ]]; then
    scdesign3_args+=(--mu-formula "$sc_mu_formula")
  fi
  if [[ -n "$sc_cache_dir" ]]; then
    scdesign3_args+=(--cache-dir "$sc_cache_dir")
  fi
  Rscript "$script_dir/analysis/cli/generate_physicell_sc_profiles.R" \
    "${scdesign3_args[@]}" 2>&1 | tee "$output_dir/scdesign3_run.log"
fi

{
  printf 'property,value\n'
  printf 'initial_ipsc_founders,%s\n' "$initial_cells"
  printf 'initial_organoid_radius_microns,%s\n' "$initial_radius"
  printf 'growth_preset,%s\n' "$growth_preset"
  printf 'target_cell_safety_stop,%s\n' "$target_cells"
  printf 'requested_culture_days,%s\n' "$culture_days"
  printf 'ipsc_cycle_rate_per_min,%s\n' "$ipsc_cycle_rate"
  printf 'radial_glia_cycle_rate_per_min,%s\n' "$radial_glia_cycle_rate"
  printf 'neural_progenitor_cycle_rate_per_min,%s\n' "$neural_progenitor_cycle_rate"
  printf 'neural_progenitor_to_neuron_rate_per_min,%s\n' \
    "$neural_progenitor_to_neuron_rate"
  printf 'ipsc_epithelial_probability,%s\n' "$epithelial_probability"
  printf 'max_epithelial_cells,%s\n' "$max_epithelial_cells"
  printf 'physicell_end_time_minutes,%s\n' "$end_time"
  printf 'final_live_cells,%s\n' "$current_cell_count"
  printf 'division_events,%s\n' "$division_count"
  printf 'founder_barcode_sites,%s\n' "$founder_label_sites"
  printf 'barcode_integrations_per_cell,%s\n' "$num_integrations"
  printf 'mitochondrial_genomes_per_cell,%s\n' "$mt_genomes_per_cell"
  printf 'recording_modalities,%s\n' "$recording_modalities"
  printf 'random_seed,%s\n' "$seed"
  printf 'physicell_source,%s\n' "$physicell_dir"
  printf 'recording_params,%s\n' "$params_path"
  printf 'scdesign3_enabled,%s\n' "$([[ -n "$sc_reference" ]] && printf true || printf false)"
  printf 'scdesign3_reference,%s\n' "$sc_reference"
} > "$output_dir/pipeline_manifest.csv"

printf '\nNeural-organoid pipeline complete.\n'
printf 'PhysiCell outputs: %s\n' "$physicell_output"
printf 'Recording outputs: %s\n' "$recording_output"
if [[ -n "$sc_reference" ]]; then
  printf 'scDesign3 outputs: %s\n' "$recording_output/sc_profiles"
fi
printf 'PIPELINE_OUTPUT_DIR=%s\n' "$output_dir"
