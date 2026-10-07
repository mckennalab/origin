#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<'EOF'
Usage:
  bash run_physicell_visium_replicates.sh [options]

Runs independent lineage-enabled PhysiCell tumors, takes conventional
Visium-style sections, and aggregates cell-cell lineage distance as a function
of distance within each section.

Options:
  --physicell-dir PATH         PhysiCell checkout [default: ../PhysiCell]
  --params PATH                Recording JSON [default: example_json_params/physicell_10000.json]
  --output-dir PATH            Batch directory [default: output/physicell_visium_<timestamp>]
  --replicates N               Independent tumor runs [default: 5]
  --first-seed N               Seed for replicate 1 [default: 1]
  --target-cells N             Terminal tumor population [default: 10000]
  --jobs N                     PhysiCell compiler/OpenMP workers [default: detected CPU count]
  --modalities LIST            Recording replay [default: lineage]
  --slice-offsets LIST         Microns from tumor center [default: -50,0,50]
  --plane-normal X,Y,Z        Section normal [default: 0,0,1]
  --section-thickness UM      Tissue-section thickness [default: 5]
  --spot-diameter UM          Conventional Visium spot diameter [default: 55]
  --spot-pitch UM             Conventional Visium center pitch [default: 100]
  --randomize-alignment BOOL  Random array rotation/translation [default: true]
  --distance-breaks LIST      Cell-distance bin edges [default: 0,20,40,60,80,100,150,200,Inf]
  --recent-mrca-hours LIST    Recent-ancestor thresholds [default: 6,12]
  --permutations N            Spatial tip-label permutations [default: 20]
  --max-pairs-per-bin N       Evaluated/output pairs per bin [default: 100000]
  --max-exact-pairs N         Exact enumeration limit [default: 2000000]
  --max-candidate-pairs N     Pair draws above exact limit [default: 2000000]
  --resume BOOL               Reuse completed replicates/analyses [default: true]
  -h, --help

The default --modalities lineage writes only the ground-truth lineage needed
for this analysis and avoids the more expensive barcode/mitochondrial replay.
Use --modalities both when those observations are also wanted.
EOF
}

physicell_dir="$script_dir/../PhysiCell"
params_path="$script_dir/example_json_params/physicell_10000.json"
output_dir=""
replicates=5
first_seed=1
target_cells=10000
modalities=lineage
slice_offsets="-50,0,50"
plane_normal="0,0,1"
section_thickness=5
spot_diameter=55
spot_pitch=100
randomize_alignment=true
distance_breaks="0,20,40,60,80,100,150,200,Inf"
recent_mrca_hours="6,12"
permutations=20
max_pairs_per_bin=100000
max_exact_pairs=2000000
max_candidate_pairs=2000000
resume=true
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
    --replicates) replicates="$2"; shift 2 ;;
    --first-seed) first_seed="$2"; shift 2 ;;
    --target-cells) target_cells="$2"; shift 2 ;;
    --jobs) jobs="$2"; shift 2 ;;
    --modalities) modalities="$2"; shift 2 ;;
    --slice-offsets) slice_offsets="$2"; shift 2 ;;
    --slice-offsets=*) slice_offsets="${1#*=}"; shift ;;
    --plane-normal) plane_normal="$2"; shift 2 ;;
    --section-thickness) section_thickness="$2"; shift 2 ;;
    --spot-diameter) spot_diameter="$2"; shift 2 ;;
    --spot-pitch) spot_pitch="$2"; shift 2 ;;
    --randomize-alignment)
      randomize_alignment="$(printf '%s' "$2" | tr '[:upper:]' '[:lower:]')"
      shift 2
      ;;
    --distance-breaks) distance_breaks="$2"; shift 2 ;;
    --recent-mrca-hours) recent_mrca_hours="$2"; shift 2 ;;
    --permutations) permutations="$2"; shift 2 ;;
    --max-pairs-per-bin) max_pairs_per_bin="$2"; shift 2 ;;
    --max-exact-pairs) max_exact_pairs="$2"; shift 2 ;;
    --max-candidate-pairs) max_candidate_pairs="$2"; shift 2 ;;
    --resume)
      resume="$(printf '%s' "$2" | tr '[:upper:]' '[:lower:]')"
      shift 2
      ;;
    -h|--help) usage; exit 0 ;;
    *)
      printf 'Unknown option: %s\n\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

for integer_value in \
  "$replicates" "$first_seed" "$target_cells" "$jobs" "$permutations" \
  "$max_pairs_per_bin" "$max_exact_pairs" "$max_candidate_pairs"; do
  if [[ ! "$integer_value" =~ ^[0-9]+$ ]]; then
    printf 'Expected a non-negative integer, received: %s\n' "$integer_value" >&2
    exit 2
  fi
done
if ((replicates < 1 || target_cells < 1 || jobs < 1 ||
     max_pairs_per_bin < 1 || max_exact_pairs < 1 ||
     max_candidate_pairs < 1)); then
  printf 'Replicate, cell, worker, and pair limits must be positive.\n' >&2
  exit 2
fi
for boolean_value in "$randomize_alignment" "$resume"; do
  if [[ "$boolean_value" != "true" && "$boolean_value" != "false" ]]; then
    printf 'Boolean options must be true or false; received: %s\n' "$boolean_value" >&2
    exit 2
  fi
done
if ! command -v Rscript >/dev/null 2>&1; then
  printf 'Rscript is required.\n' >&2
  exit 1
fi

physicell_dir="$(cd "$physicell_dir" && pwd)"
params_path="$(cd "$(dirname "$params_path")" && pwd)/$(basename "$params_path")"
if [[ -z "$output_dir" ]]; then
  output_dir="$script_dir/output/physicell_visium_$(date '+%Y%m%d_%H%M%S')"
elif [[ "$output_dir" != /* ]]; then
  output_dir="$script_dir/$output_dir"
fi
if [[ -e "$output_dir" && "$resume" == "false" &&
      -n "$(find "$output_dir" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  printf 'Batch output exists and --resume is false: %s\n' "$output_dir" >&2
  exit 1
fi
mkdir -p "$output_dir"

manifest_path="$output_dir/replicate_manifest.csv"
printf 'replicate,seed,run_dir,simulation_status,analysis_status\n' > "$manifest_path"

for ((replicate = 1; replicate <= replicates; replicate++)); do
  seed=$((first_seed + replicate - 1))
  replicate_label="$(printf 'replicate_%03d_seed_%06d' "$replicate" "$seed")"
  replicate_dir="$output_dir/$replicate_label"
  simulation_status=completed
  analysis_status=completed

  printf '\n[%s/%s] PhysiCell replicate %s (seed %s)\n' \
    "$replicate" "$replicates" "$replicate_label" "$seed"
  if [[ "$resume" == "true" && -s "$replicate_dir/pipeline_manifest.csv" &&
        ( -s "$replicate_dir/lineage_recording/lineage_nodes.csv.gz" ||
          -s "$replicate_dir/lineage_recording/lineage_nodes.csv" ) &&
        ( -s "$replicate_dir/lineage_recording/terminal_cells.csv.gz" ||
          -s "$replicate_dir/lineage_recording/terminal_cells.csv" ) ]]; then
    printf 'Reusing completed PhysiCell lineage: %s\n' "$replicate_dir"
    simulation_status=reused
  else
    if [[ -e "$replicate_dir" &&
          -n "$(find "$replicate_dir" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
      printf 'Incomplete non-empty replicate cannot be overwritten safely: %s\n' \
        "$replicate_dir" >&2
      exit 1
    fi
    bash "$script_dir/run_physicell_10000_pipeline.sh" \
      --physicell-dir "$physicell_dir" \
      --params "$params_path" \
      --output-dir "$replicate_dir" \
      --target-cells "$target_cells" \
      --modalities "$modalities" \
      --seed "$seed" \
      --jobs "$jobs"
  fi

  analysis_dir="$replicate_dir/visium_spatial"
  if [[ "$resume" == "true" && -s "$analysis_dir/analysis_complete.txt" &&
        ( -s "$analysis_dir/cell_distance_summary.csv.gz" ||
          -s "$analysis_dir/cell_distance_summary.csv" ) ]]; then
    printf 'Reusing completed Visium analysis: %s\n' "$analysis_dir"
    analysis_status=reused
  else
    Rscript "$script_dir/analysis/cli/analyze_physicell_visium.R" \
      --run-dir "$replicate_dir" \
      --output-dir "$analysis_dir" \
      --slice-offsets "$slice_offsets" \
      --plane-normal "$plane_normal" \
      --section-thickness "$section_thickness" \
      --spot-diameter "$spot_diameter" \
      --spot-pitch "$spot_pitch" \
      --randomize-alignment "$randomize_alignment" \
      --distance-breaks "$distance_breaks" \
      --recent-mrca-hours "$recent_mrca_hours" \
      --permutations "$permutations" \
      --max-pairs-per-bin "$max_pairs_per_bin" \
      --max-exact-pairs "$max_exact_pairs" \
      --max-candidate-pairs "$max_candidate_pairs" \
      --seed "$seed"
  fi

  printf '%s,%s,%s,%s,%s\n' \
    "$replicate" "$seed" "$replicate_dir" \
    "$simulation_status" "$analysis_status" >> "$manifest_path"
done

Rscript "$script_dir/analysis/cli/aggregate_physicell_visium.R" \
  --batch-dir "$output_dir" \
  --output-dir "$output_dir/aggregate"

printf '\nVisium replicate simulation and analysis complete.\n'
printf 'Replicate manifest: %s\n' "$manifest_path"
printf 'Aggregate output: %s\n' "$output_dir/aggregate"
printf 'BATCH_OUTPUT_DIR=%s\n' "$output_dir"
