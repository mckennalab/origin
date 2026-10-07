#!/usr/bin/env bash
# C. Capacity-matched recorder comparisons.
#
# Matching integrations compares constructs of very different size, because the
# recorders carry different numbers of targets per integration: BASELINE 50,
# WT-CRISPR 10, PEtracer 3, PALINCODE 2. These runs instead hold a budget fixed
# and give each recorder the integration count that spends it.
#
# Three budgets:
#
#   targets      100 targets, neutral shape
#   information  133 bits of capacity, neutral shape
#   allshapes    133 bits, across all four tree shapes
#
# Every level is generated for every system so the nested-prefix construction is
# identical across recorders; only the matched level is read out per system. The
# comparison is capacity-matched but not cost-matched -- 50 integrations is a far
# harder construct to build than 2 -- and it holds targets fixed while letting
# integration count vary, so any effect of spreading the same targets over more
# independent integrations lands inside the recorder comparison rather than
# beside it.
#
#   bash analysis/runs/run_c_matched_capacity.sh            # all three
#   BUDGETS=information bash analysis/runs/run_c_matched_capacity.sh

source "$(dirname -- "${BASH_SOURCE[0]}")/common.sh"

BUDGETS=${BUDGETS:-targets information allshapes}
SYSTEMS=${SYSTEMS:-baseline,wt_crispr,prime,palincode}

# Per budget: run name, the levels to generate, each system's matched level, and
# the shapes to read out. Mitochondrial recording is excluded throughout: it has
# no integrations, and its capacity axis is recovered genome depth, so there is
# no target count to match on.
run_budget () {
  local budget=$1
  local name levels matched shapes design
  case "$budget" in
    targets)
      name=matched_targets_nj_n250
      levels=2,10,33,50
      matched="baseline=2 wt_crispr=10 prime=33 palincode=50"
      shapes=neutral ;;
    information)
      name=matched_information_nj_n250
      levels=2,13,14,30
      matched="baseline=2 wt_crispr=13 prime=14 palincode=30"
      shapes=neutral ;;
    allshapes)
      name=matched_information_nj_allshapes_n250
      levels=2,13,14,30
      matched="baseline=2 wt_crispr=13 prime=14 palincode=30"
      shapes=balanced,comb,neutral,hierarchical ;;
    *)
      echo "Unknown budget: $budget" >&2; return 1 ;;
  esac
  local output_root="$OUTPUT_BASE/$name"
  mkdir -p "$output_root"

  announce "C/$budget: simulate -> $output_root"
  Rscript "$HARNESS_DIR/generate-recorder-parameter-grid.R" \
    --source-root="$REPO_ROOT" \
    --benchmark-root="$SOURCE_BENCHMARK" \
    --output-dir="$output_root/simulation" \
    --grid=best \
    --extension-doublings=0 \
    --systems="$SYSTEMS" \
    --shapes="$shapes" \
    --levels="$levels" \
    --seeds="$SEEDS" \
    --workers="$WORKERS"

  # design.csv is what the summarising and assessment scripts read to learn
  # which level each recorder was matched at; without it they cannot tell a
  # target-matched run from an information-matched one.
  design="$output_root/design.csv"
  {
    echo "system,label,targets_per_integration,integrations"
    for pair in $matched; do
      local system=${pair%%=*}
      local level=${pair##*=}
      case "$system" in
        baseline)  echo "baseline,BASELINE (cas12a),50,$level" ;;
        wt_crispr) echo "wt_crispr,WT-CRISPR (FLARE),10,$level" ;;
        prime)     echo "prime,PEtracer (prime),3,$level" ;;
        palincode) echo "palincode,PALINCODE,2,$level" ;;
      esac
    done
  } > "$design"

  for shape in ${shapes//,/ }; do
    for pair in $matched; do
      local system=${pair%%=*}
      local level=${pair##*=}
      announce "C/$budget: $shape, $system at $level integrations"
      local destination="$output_root/reconstruction/${system}_k_${level}"
      if [ "$shapes" != "neutral" ]; then
        destination="$output_root/reconstruction/$shape/${system}_k_${level}"
      fi
      Rscript "$CLIQUE_R_DIR/inst/examples/process-lineage-benchmark.R" \
        --benchmark-dir="$output_root/simulation" \
        --package-dir="$CLIQUE_R_DIR" \
        --output-dir="$destination" \
        --methods=nj --matrix=character_matrix --shapes="$shape" \
        --systems="$system" --integrations="$level" \
        --observation-depths="$level" --seeds="$SEEDS" \
        --workers="$TREE_WORKERS" --threads=1 --resume=true
    done
  done

  announce "C/$budget: summarise, assess and plot"
  Rscript "$REPO_ROOT/analysis/cli/summarize_matched_targets_nj.R" "$output_root"
  for shape in ${shapes//,/ }; do
    Rscript "$REPO_ROOT/analysis/cli/assess_matched_targets_trees.R" \
      "$output_root" "$output_root/matched_tree_metrics_$shape" "$shape"
  done
  # The figure. A single-shape budget gets the per-recorder metric panel; the
  # all-shapes budget gets the across-shape panel instead, at k=5 clones, which
  # is the granularity that panel reads at.
  mkdir -p "$REPO_ROOT/analysis/figures"
  if [ "$budget" = "allshapes" ]; then
    Rscript "$REPO_ROOT/analysis/cli/plot_allshapes_recorder_metrics.R" \
      "$output_root" "$REPO_ROOT/analysis/figures/allshapes_matched_information" 5
  else
    Rscript "$REPO_ROOT/analysis/cli/plot_matched_recorder_metrics.R" \
      "$output_root" "$REPO_ROOT/analysis/figures/${budget}_matched_metrics" \
      "matched $budget budget, neutral trees, 250 cells" \
      "matched_tree_metrics_neutral"
  fi
  announce "C/$budget: done -> $output_root"
}

for budget in $BUDGETS; do
  run_budget "$budget"
done
