# Shared settings for the large-scale tree-building runs.
#
# Outputs land under this repository's output/ directory, which is ignored, so a
# run writes gigabytes without touching what is tracked. The harness scripts and
# the cliqueR package live in a separate checkout; point CLIQUE_DIR at it if it
# is not beside this one.
#
# Every script here is resumable: reconstruction is invoked with --resume=true,
# so re-running after an interruption picks up the conditions already built
# rather than redoing them.

set -euo pipefail

RUNS_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "$RUNS_DIR/../.." && pwd)
CLIQUE_DIR=${CLIQUE_DIR:-"$(cd -- "$REPO_ROOT/../clique_2025_12_10" 2>/dev/null && pwd || true)"}
if [ -z "${CLIQUE_DIR:-}" ] || [ ! -d "$CLIQUE_DIR" ]; then
  echo "Set CLIQUE_DIR to the clique checkout holding r_pkg/ and rust_cmd/." >&2
  exit 1
fi
HARNESS_DIR="$CLIQUE_DIR/rust_cmd/test_harness/tree_building"
CLIQUE_R_DIR="$CLIQUE_DIR/r_pkg/cliqueR"
SOURCE_BENCHMARK=${SOURCE_BENCHMARK:-"$REPO_ROOT/output/lineage_benchmark_20260819_173705"}
OUTPUT_BASE=${OUTPUT_BASE:-"$REPO_ROOT/output/tree_building"}

SEEDS=${SEEDS:-1,2,3,4,5}
WORKERS=${WORKERS:-4}
TREE_WORKERS=${TREE_WORKERS:-8}
# Cassiopeia's Python backend is loaded through reticulate, and forking an R
# worker after Python is initialised crashes on macOS. Methods that need it run
# in a single process; the rest run in parallel.
FORKING_METHODS=${FORKING_METHODS:-nj,parsimony,iqtree,mix,vine}
ALL_METHODS=${ALL_METHODS:-nj,parsimony,iqtree,mix,vine,cassiopeia}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/cliqueR-mpl-cache}

for required in "$HARNESS_DIR" "$CLIQUE_R_DIR" "$SOURCE_BENCHMARK"; do
  if [ ! -d "$required" ]; then
    echo "Missing required path: $required" >&2
    exit 1
  fi
done
mkdir -p "$OUTPUT_BASE"

# Reconstruct a benchmark directory with every method, splitting the run so
# Cassiopeia never executes inside a forked worker.
reconstruct_all_methods () {
  local benchmark_dir=$1
  local output_dir=$2
  shift 2
  Rscript "$CLIQUE_R_DIR/inst/examples/process-lineage-benchmark.R" \
    --benchmark-dir="$benchmark_dir" --package-dir="$CLIQUE_R_DIR" \
    --output-dir="$output_dir" --methods="$FORKING_METHODS" \
    --matrix=character_matrix --seeds="$SEEDS" --workers="$TREE_WORKERS" \
    --threads=1 --vine-nj-only=false --resume=true "$@"
  Rscript "$CLIQUE_R_DIR/inst/examples/process-lineage-benchmark.R" \
    --benchmark-dir="$benchmark_dir" --package-dir="$CLIQUE_R_DIR" \
    --output-dir="$output_dir" --methods="$ALL_METHODS" \
    --matrix=character_matrix --seeds="$SEEDS" --workers=1 \
    --threads=1 --vine-nj-only=false --resume=true "$@"
}

announce () {
  echo
  echo "=== $* ==="
  date "+%Y-%m-%d %H:%M:%S"
}
