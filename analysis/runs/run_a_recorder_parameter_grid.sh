#!/usr/bin/env bash
# A. Recorder parameter grid -- the main comparison panel.
#
# Five recorders at five editing-rate tiers, on four tree shapes, five seeds,
# reconstructed by six methods: 500 conditions. This is the run behind
# recorder_parameter_grid_by_tree_shape and the accuracy_by_recorder_* tables.
#
# FLARE is simulated with interval deletions and the allele-state character
# matrix, which is what the recorder now does by default; the published grid
# predating that change has FLARE from the point-deletion model and is not
# comparable for that one system.
#
# Roughly 500 conditions x 4 minutes of reconstruction, dominated by MIX,
# parsimony and IQ-TREE. Resumable.
#
#   bash analysis/runs/run_a_recorder_parameter_grid.sh
#   SYSTEMS=wt_crispr bash analysis/runs/run_a_recorder_parameter_grid.sh

source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

NAME=${NAME:-recorder_parameter_grid_n250_level05}
OUTPUT_ROOT="$OUTPUT_BASE/$NAME"
GRID=${GRID:-wide}
SYSTEMS=${SYSTEMS:-baseline,wt_crispr,prime,palincode,mitochondrial}
SHAPES=${SHAPES:-balanced,comb,neutral,hierarchical}
LEVEL=${LEVEL:-5}

announce "A: recorder parameter grid -> $OUTPUT_ROOT"
Rscript "$HARNESS_DIR/generate-recorder-parameter-grid.R" \
  --source-root="$REPO_ROOT" \
  --benchmark-root="$SOURCE_BENCHMARK" \
  --output-dir="$OUTPUT_ROOT/simulation" \
  --grid="$GRID" \
  --extension-doublings=0 \
  --systems="$SYSTEMS" \
  --shapes="$SHAPES" \
  --levels="$LEVEL" \
  --seeds="$SEEDS" \
  --workers="$WORKERS"

announce "A: reconstruction"
reconstruct_all_methods "$OUTPUT_ROOT/simulation" "$OUTPUT_ROOT/reconstruction"

announce "A: summarise and plot"
Rscript "$HARNESS_DIR/summarize-recorder-parameter-grid.R" "$OUTPUT_ROOT"
Rscript "$HARNESS_DIR/plot-recorder-parameter-grid.R" "$OUTPUT_ROOT"
announce "A: done -> $OUTPUT_ROOT"
