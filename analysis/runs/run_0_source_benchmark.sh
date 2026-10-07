#!/usr/bin/env bash
# 0. The source benchmark that runs A, B and C all read from.
#
# A simulation-only sweep: population shapes x seeds x recorder systems x sample
# sizes x integration counts, writing the ground-truth tree and the sparse
# state and character matrices for every condition. Reconstruction happens
# later, in the runs that consume this.
#
# The settings below reproduce output/lineage_benchmark_20260819_173705, read
# from the benchmark_settings.json that run wrote. Note the systems list:
# wt_crispr is absent, because the FLARE recorder is injected by the harness at
# grid-generation time rather than being one of the sweep's native systems.
#
# This is the expensive step and the prerequisite for everything else: the
# existing copy is 7.1GB. It is resumable, and A, B and C will reuse an existing
# directory rather than rebuild it, so run this only when there is no source
# benchmark or when the sweep's definition changes.
#
#   bash analysis/runs/run_0_source_benchmark.sh
#   OUTPUT_NAME=lineage_benchmark_small TREE_SIZES=250 SEEDS=1,2 \
#     bash analysis/runs/run_0_source_benchmark.sh

source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

OUTPUT_NAME=${OUTPUT_NAME:-lineage_benchmark_$(date +%Y%m%d_%H%M%S)}
OUTPUT_DIR="$REPO_ROOT/output/$OUTPUT_NAME"
SHAPES=${SHAPES:-balanced,comb,neutral,hierarchical}
TREE_SIZES=${TREE_SIZES:-250,1000,2000,5000}
INTEGRATION_COUNTS=${INTEGRATION_COUNTS:-1,2,3,5,10,20,30,40,50}
MT_DEPTHS=${MT_DEPTHS:-1,2,3,5,10,20,30,40,50}
MT_GENOMES=${MT_GENOMES:-50}
BENCH_SYSTEMS=${BENCH_SYSTEMS:-baseline,prime,palincode,mitochondrial}
BENCH_SEEDS=${BENCH_SEEDS:-1,2,3,4,5,6,7,8,9,10}
DURATION=${DURATION:-10}

announce "0: source benchmark -> $OUTPUT_DIR"
Rscript "$REPO_ROOT/analysis/cli/run_lineage_benchmark.R" \
  --output-dir="$OUTPUT_DIR" \
  --shapes="$SHAPES" \
  --tree-sizes="$TREE_SIZES" \
  --integration-counts="$INTEGRATION_COUNTS" \
  --mt-observation-depths="$MT_DEPTHS" \
  --mt-genomes-per-cell="$MT_GENOMES" \
  --systems="$BENCH_SYSTEMS" \
  --seeds="$BENCH_SEEDS" \
  --synthetic-duration="$DURATION" \
  --write-dense-csv=false \
  --resume=true \
  --progress=true

announce "0: done -> $OUTPUT_DIR"
echo "Point the other runs at it with:"
echo "  export SOURCE_BENCHMARK=$OUTPUT_DIR"
