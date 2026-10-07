# Large-scale run drivers

Each script reproduces one family of simulations end to end: simulate,
reconstruct, summarise, plot. Outputs go to `output/tree_building/<run name>/`
in this repository, which is ignored, so a run writes gigabytes without touching
anything tracked.

| script | what it produces |
|--------|------------------|
| `run_0_source_benchmark.sh` | the simulated populations every other run reads from; 7.1GB, and the prerequisite for all of them |
| `run_a_recorder_parameter_grid.sh` | 5 recorders x 5 rate tiers x 4 shapes x 5 seeds x 6 methods (500 conditions) |
| `run_b_integration_series.sh` | every recorder at 1 to 100 integrations, neutral shape, neighbour joining |
| `run_c_matched_capacity.sh` | three capacity-matched comparisons: 100 targets, 133 bits, and 133 bits across all shapes |

## Running them

```bash
bash analysis/runs/run_a_recorder_parameter_grid.sh
```

The harness scripts and the cliqueR package live in a separate checkout. The
scripts look for it beside this one; set `CLIQUE_DIR` if it is somewhere else.
They also need a source benchmark to sample lineages from, defaulting to
`output/lineage_benchmark_20260819_173705`. If that directory is absent, build
one with `run_0_source_benchmark.sh` and point the others at it:

```bash
bash analysis/runs/run_0_source_benchmark.sh
export SOURCE_BENCHMARK=$PWD/output/lineage_benchmark_<timestamp>
```

Run 0 is the expensive step and is where an end-to-end rebuild starts. Note that
its systems list does not include `wt_crispr`: the FLARE recorder is injected by
the harness when the grid is generated rather than being one of the sweep's
native systems.

Everything is overridable from the environment, which is how to run a subset:

```bash
SYSTEMS=wt_crispr bash analysis/runs/run_a_recorder_parameter_grid.sh
BUDGETS=information bash analysis/runs/run_c_matched_capacity.sh
SEEDS=1 WORKERS=1 bash analysis/runs/run_b_integration_series.sh
```

## Two things worth knowing before a long run

**They are resumable.** Reconstruction is invoked with `--resume=true`, so
re-running after an interruption picks up the conditions already built. Deleting
the output directory is the way to force a clean rebuild.

**Cassiopeia cannot run in a forked worker.** Its Python backend is loaded
through reticulate, and forking an R worker after Python is initialised crashes
on macOS. `reconstruct_all_methods` therefore runs the other five methods in
parallel and then makes a second single-process pass that adds Cassiopeia.

## Scale

Reconstruction dominates. Measured per condition: VINE 0.06s, NJ 0.42s,
Cassiopeia 1.3s, IQ-TREE 41s, parsimony 64s, MIX 134s -- about four minutes for
the full six-method panel. Run A is 500 conditions, so budget hours and expect
the summariser to refuse until every method has completed, which is deliberate:
a partial grid would otherwise summarise as though it were whole.
