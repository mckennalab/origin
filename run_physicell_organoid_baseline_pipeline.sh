#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec bash "$script_dir/run_physicell_organoid_pipeline.sh" \
  --params "$script_dir/example_json_params/physicell_neural_organoid_baseline.json" \
  --num-integrations 5 \
  --founder-label-sites 12 \
  --mt-genomes-per-cell 32 \
  "$@"
