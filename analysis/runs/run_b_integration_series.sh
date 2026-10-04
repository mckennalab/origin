#!/usr/bin/env bash
# B. Integration series -- accuracy against the number of integrations.
#
# Every recorder at 1, 5, 10, 20, 30, 40, 50 and 100 integrations on the neutral
# shape, reconstructed with neighbour joining, five seeds. The levels are nested:
# a smaller level reads a prefix of the same realisation rather than a separately
# simulated recorder, so the curve varies how many integrations EXIST while
# holding everything else fixed.
#
# This supersedes the two earlier series that stopped at 50 or started at 30.
# Both predate PEtracer and still carry the single-outcome prime editor, whose
# nRF at 100 integrations was 0.153 against the current 0.0996.
#
#   bash analysis/runs/run_b_integration_series.sh

source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

NAME=${NAME:-neutral_best_nj_integration_series_nested_n250}
OUTPUT_ROOT="$OUTPUT_BASE/$NAME"
LEVELS=${LEVELS:-1,5,10,20,30,40,50,100}
SYSTEMS=${SYSTEMS:-baseline,wt_crispr,prime,palincode,mitochondrial}
MITO_GENOMES=${MITO_GENOMES:-100}

announce "B: integration series -> $OUTPUT_ROOT"
Rscript "$HARNESS_DIR/generate-recorder-parameter-grid.R" \
  --source-root="$REPO_ROOT" \
  --benchmark-root="$SOURCE_BENCHMARK" \
  --output-dir="$OUTPUT_ROOT/simulation" \
  --grid=best \
  --extension-doublings=0 \
  --systems="$SYSTEMS" \
  --shapes=neutral \
  --levels="$LEVELS" \
  --seeds="$SEEDS" \
  --mitochondrial-genomes="$MITO_GENOMES" \
  --workers="$WORKERS"

# Each level is read out separately: mitochondrial recovery is a sampling depth
# rather than an integration count, so it shares the axis but not the mechanism.
for level in ${LEVELS//,/ }; do
  announce "B: reconstruction at $level integrations"
  Rscript "$CLIQUE_R_DIR/inst/examples/process-lineage-benchmark.R" \
    --benchmark-dir="$OUTPUT_ROOT/simulation" \
    --package-dir="$CLIQUE_R_DIR" \
    --output-dir="$OUTPUT_ROOT/reconstruction_k_$level" \
    --methods=nj --matrix=character_matrix --shapes=neutral \
    --systems="$SYSTEMS" --integrations="$level" \
    --observation-depths="$level" --seeds="$SEEDS" \
    --workers="$TREE_WORKERS" --threads=1 --resume=true
  Rscript "$HARNESS_DIR/summarize-neutral-best-integration-level.R" \
    "$OUTPUT_ROOT" "$level"
done

announce "B: finalise and plot"
Rscript "$HARNESS_DIR/finalize-neutral-best-integration-series.R" "$OUTPUT_ROOT"
Rscript "$HARNESS_DIR/plot-neutral-best-nj-integration-series.R" "$OUTPUT_ROOT"
announce "B: done -> $OUTPUT_ROOT"
