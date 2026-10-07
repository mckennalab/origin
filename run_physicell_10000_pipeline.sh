#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage:
  bash run_physicell_10000_pipeline.sh [options]

Runs the lineage-enabled PhysiCell 3-D tumor project to a cell-count threshold,
then writes ground truth and replays the requested recording modalities on the
exported division tree.

Options:
  --physicell-dir PATH       PhysiCell checkout [default: ../PhysiCell]
  --params PATH              Recording JSON [default: example_json_params/physicell_10000.json]
  --output-dir PATH          Run directory [default: output/physicell_10000_<timestamp>]
  --target-cells N           PhysiCell stopping threshold [default: 10000]
  --num-integrations N       Barcodes per cell [default: 1]
  --mt-genomes-per-cell N    Fixed mt bottleneck size [default: 8]
  --write-mt-fasta BOOL      Write one 16.6 kb mt haplotype per cell [default: false]
  --modalities LIST          lineage, barcode, mitochondrial, ecDNA, both, or all [default: both]
  --seed N                   PhysiCell and recording seed [default: 1]
  --jobs N                   Parallel compiler jobs [default: detected CPU count]
  --sc-reference PATH        Optional SCE/Seurat RDS or h5ad; enables scDesign3
  --sc-celltype-col NAME     Reference cell-type column [default: cell_type]
  --sc-celltype-map PATH     Optional JSON simulated-to-reference type mapping
  --sc-use-pseudotime        Model normalized lineage depth as pseudotime
  --sc-pseudotime-col NAME   Reference pseudotime column [default: pseudotime]
  --sc-spatial-cols X,Y      Optional reference spatial-coordinate columns
  --sc-other-covariates LIST Comma-delimited matched reference/PhysiCell covariates
  --sc-mu-formula FORMULA    Optional explicit scDesign3 mean formula
  --sc-cache-dir PATH        Optional shared scDesign3 fit-cache directory
  --sc-ncores N              scDesign3 worker count [default: --jobs]
  -h, --help

The PhysiCell source checkout is never modified. A self-contained staged build
is created inside the run directory.
EOF
}

physicell_dir="$script_dir/../PhysiCell"
params_path="$script_dir/example_json_params/physicell_10000.json"
output_dir=""
target_cells=10000
num_integrations=1
mt_genomes_per_cell=8
write_mt_fasta=false
modalities=both
seed=1
sc_reference=""
sc_celltype_col=cell_type
sc_celltype_map=""
sc_use_pseudotime=false
sc_pseudotime_col=pseudotime
sc_spatial_cols=""
sc_other_covariates=""
sc_mu_formula=""
sc_cache_dir=""
sc_ncores=""
scdesign3_enabled=false
if command -v sysctl >/dev/null 2>&1; then
  jobs="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
elif command -v nproc >/dev/null 2>&1; then
  jobs="$(nproc)"
else
  jobs=4
fi

while (($# > 0)); do
  case "$1" in
    --physicell-dir)
      physicell_dir="$2"
      shift 2
      ;;
    --params)
      params_path="$2"
      shift 2
      ;;
    --output-dir)
      output_dir="$2"
      shift 2
      ;;
    --target-cells)
      target_cells="$2"
      shift 2
      ;;
    --num-integrations)
      num_integrations="$2"
      shift 2
      ;;
    --mt-genomes-per-cell)
      mt_genomes_per_cell="$2"
      shift 2
      ;;
    --write-mt-fasta)
      write_mt_fasta="$(printf '%s' "$2" | tr '[:upper:]' '[:lower:]')"
      shift 2
      ;;
    --modalities)
      modalities="$2"
      shift 2
      ;;
    --seed)
      seed="$2"
      shift 2
      ;;
    --jobs)
      jobs="$2"
      shift 2
      ;;
    --sc-reference)
      sc_reference="$2"
      shift 2
      ;;
    --sc-celltype-col)
      sc_celltype_col="$2"
      shift 2
      ;;
    --sc-celltype-map)
      sc_celltype_map="$2"
      shift 2
      ;;
    --sc-use-pseudotime)
      sc_use_pseudotime=true
      shift
      ;;
    --sc-pseudotime-col)
      sc_pseudotime_col="$2"
      shift 2
      ;;
    --sc-spatial-cols)
      sc_spatial_cols="$2"
      shift 2
      ;;
    --sc-other-covariates)
      sc_other_covariates="$2"
      shift 2
      ;;
    --sc-mu-formula)
      sc_mu_formula="$2"
      shift 2
      ;;
    --sc-cache-dir)
      sc_cache_dir="$2"
      shift 2
      ;;
    --sc-ncores)
      sc_ncores="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown option: %s\n\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "$sc_ncores" ]]; then
  sc_ncores="$jobs"
fi
for integer_value in \
  "$target_cells" \
  "$num_integrations" \
  "$mt_genomes_per_cell" \
  "$seed" \
  "$jobs" \
  "$sc_ncores"; do
  if [[ ! "$integer_value" =~ ^[0-9]+$ ]]; then
    printf 'Expected a non-negative integer, received: %s\n' "$integer_value" >&2
    exit 2
  fi
done
if ((target_cells < 1 || num_integrations < 1 || mt_genomes_per_cell < 1 ||
     jobs < 1 || sc_ncores < 1)); then
  printf 'Cell, integration, mt-genome, and worker counts must be positive.\n' >&2
  exit 2
fi
if [[ "$write_mt_fasta" != "true" && "$write_mt_fasta" != "false" ]]; then
  printf -- '--write-mt-fasta must be true or false.\n' >&2
  exit 2
fi

physicell_dir="$(cd "$physicell_dir" && pwd)"
params_path="$(cd "$(dirname "$params_path")" && pwd)/$(basename "$params_path")"
if [[ -n "$sc_reference" ]]; then
  scdesign3_enabled=true
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
project_source="$physicell_dir/user_projects/tumor_3D_lineage"
for required_path in \
  "$physicell_dir/core" \
  "$physicell_dir/modules" \
  "$physicell_dir/BioFVM" \
  "$project_source/main.cpp" \
  "$project_source/Makefile" \
  "$project_source/config/PhysiCell_settings.xml" \
  "$params_path"; do
  if [[ ! -e "$required_path" ]]; then
    printf 'Required input not found: %s\n' "$required_path" >&2
    exit 1
  fi
done
if ! command -v make >/dev/null 2>&1 || ! command -v g++ >/dev/null 2>&1; then
  printf 'A working make and g++ installation is required.\n' >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  printf 'Rscript is required for lineage-recording replay.\n' >&2
  exit 1
fi

if [[ -z "$output_dir" ]]; then
  output_dir="$script_dir/output/physicell_10000_$(date '+%Y%m%d_%H%M%S')"
elif [[ "$output_dir" != /* ]]; then
  output_dir="$script_dir/$output_dir"
fi
if [[ -e "$output_dir" ]] && [[ -n "$(find "$output_dir" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  printf 'Output directory already exists and is not empty: %s\n' "$output_dir" >&2
  exit 1
fi

build_dir="$output_dir/physicell_build"
physicell_output="$build_dir/output"
recording_output="$output_dir/lineage_recording"
mkdir -p "$build_dir" "$recording_output"

cp "$project_source/main.cpp" "$build_dir/main.cpp"
cp "$project_source/Makefile" "$build_dir/Makefile"
cp "$project_source/VERSION.txt" "$build_dir/VERSION.txt"
cp -R "$project_source/config" "$build_dir/config"
cp -R "$project_source/custom_modules" "$build_dir/custom_modules"
ln -s "$physicell_dir/core" "$build_dir/core"
ln -s "$physicell_dir/modules" "$build_dir/modules"
ln -s "$physicell_dir/BioFVM" "$build_dir/BioFVM"

TARGET_CELLS="$target_cells" perl -0pi -e \
  's/int target_number_of_cells = [0-9]+;/int target_number_of_cells = $ENV{TARGET_CELLS};/' \
  "$build_dir/main.cpp"
if ! grep -q "int target_number_of_cells = ${target_cells};" "$build_dir/main.cpp"; then
  printf 'Could not configure the staged PhysiCell target cell count.\n' >&2
  exit 1
fi

PHYSICELL_SEED="$seed" PHYSICELL_THREADS="$jobs" perl -0pi -e '
  s{(<full_data>.*?<enable>)[^<]*(</enable>)}{$1false$2}s;
  s{(<SVG>.*?<enable>)[^<]*(</enable>)}{$1false$2}s;
  s{<random_seed>[^<]*</random_seed>}{<random_seed>$ENV{PHYSICELL_SEED}</random_seed>};
  s{<omp_num_threads>[^<]*</omp_num_threads>}{<omp_num_threads>$ENV{PHYSICELL_THREADS}</omp_num_threads>};
' "$build_dir/config/PhysiCell_settings.xml"

cp "$params_path" "$output_dir/recording_params.json"

printf 'Building staged PhysiCell tumor project with %s job(s)...\n' "$jobs"
(
  cd "$build_dir"
  make -j "$jobs"
) 2>&1 | tee "$output_dir/physicell_build.log"

printf 'Growing the PhysiCell tumor to at least %s current cells...\n' "$target_cells"
(
  cd "$build_dir"
  ./project ./config/PhysiCell_settings.xml
) 2>&1 | tee "$output_dir/physicell_run.log"

lineage_path="$physicell_output/cell_lineage.csv"
live_cells_path="$physicell_output/lineage_table.csv"
if [[ ! -s "$lineage_path" || ! -s "$live_cells_path" ]]; then
  printf 'PhysiCell did not produce the expected lineage CSV files.\n' >&2
  exit 1
fi

current_cell_count="$(awk 'NR > 1 { count++ } END { print count + 0 }' "$live_cells_path")"
division_count="$(awk 'NR > 1 { count++ } END { print count + 0 }' "$lineage_path")"
if ((current_cell_count < target_cells)); then
  printf 'PhysiCell stopped with %s current cells, below target %s.\n' \
    "$current_cell_count" "$target_cells" >&2
  exit 1
fi

end_time="$(sed -nE 's/.* at t = ([0-9.eE+-]+) min.*/\1/p' \
  "$output_dir/physicell_run.log" | tail -n 1)"
if [[ -z "$end_time" ]]; then
  printf 'Could not recover the final PhysiCell time from its run log.\n' >&2
  exit 1
fi

printf 'Writing lineage output and replaying requested modalities (%s) on %s division events...\n' \
  "$modalities" "$division_count"
Rscript "$script_dir/origin/inst/scripts/simulate_physicell_lineage.R" \
  --lineage "$lineage_path" \
  --live-cells "$live_cells_path" \
  --params "$params_path" \
  --output-dir "$recording_output" \
  --end-time "$end_time" \
  --editing-state auto \
  --modalities "$modalities" \
  --num-integrations "$num_integrations" \
  --mt-genomes-per-cell "$mt_genomes_per_cell" \
  --write-mt-fasta "$write_mt_fasta" \
  --seed "$seed" 2>&1 | tee "$output_dir/recording_run.log"

if [[ -n "$sc_reference" ]]; then
  printf 'Generating scDesign3 expression profiles for terminal cells...\n'
  scdesign3_args=(
    --run-dir "$output_dir"
    --reference "$sc_reference"
    --celltype-col "$sc_celltype_col"
    --pseudotime-col "$sc_pseudotime_col"
    --ncores "$sc_ncores"
    --seed "$seed"
  )
  if [[ -n "$sc_celltype_map" ]]; then
    scdesign3_args+=(--celltype-map "$sc_celltype_map")
  fi
  if [[ "$sc_use_pseudotime" == "true" ]]; then
    scdesign3_args+=(--use-pseudotime)
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
    "${scdesign3_args[@]}" | tee "$output_dir/scdesign3_run.log"
fi

{
  printf 'property,value\n'
  printf 'target_cells,%s\n' "$target_cells"
  printf 'final_current_cells,%s\n' "$current_cell_count"
  printf 'division_events,%s\n' "$division_count"
  printf 'physicell_end_time,%s\n' "$end_time"
  printf 'barcode_integrations_per_cell,%s\n' "$num_integrations"
  printf 'mitochondrial_genomes_per_cell,%s\n' "$mt_genomes_per_cell"
  printf 'recording_modalities,%s\n' "$modalities"
  printf 'random_seed,%s\n' "$seed"
  printf 'physicell_source,%s\n' "$physicell_dir"
  printf 'recording_params,%s\n' "$params_path"
  printf 'scdesign3_enabled,%s\n' "$scdesign3_enabled"
  printf 'scdesign3_reference,%s\n' "$sc_reference"
  printf 'scdesign3_use_pseudotime,%s\n' "$sc_use_pseudotime"
  printf 'scdesign3_workers,%s\n' "$sc_ncores"
} > "$output_dir/pipeline_manifest.csv"

printf '\nPipeline complete.\n'
printf 'PhysiCell outputs: %s\n' "$physicell_output"
printf 'Recording outputs: %s\n' "$recording_output"
if [[ -n "$sc_reference" ]]; then
  printf 'scDesign3 outputs: %s\n' "$recording_output/sc_profiles"
fi
printf 'PIPELINE_OUTPUT_DIR=%s\n' "$output_dir"
