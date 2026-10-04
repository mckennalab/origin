# Repository layout

The code is split into three tiers. Dependencies point strictly inward: the
core depends on nothing here, analysis depends on the core, and legacy depends
on the core but nothing depends on legacy.

```
origin/          R package - the simulator core
  R/             6 simulator files, 109 exported functions
  man/           generated from roxygen; do not edit by hand
  inst/scripts/  command-line drivers for the core
analysis/        benchmarking and analysis library (sourced, not a package)
  cli/           command-line drivers for the analysis tier
legacy/          the sim5 timestep engine and its satellites
load_origin.R    sources core + analysis into the global environment
```

## origin/ - the package

The core simulator: population growth (Gillespie or imported PhysiCell division
history) and recorder replay along the resulting lineages. It imports only
Matrix, methods, stats, and utils, so installing it does not pull in
Bioconductor or plotting stacks.

```r
install.packages("origin", repos = NULL, type = "source")   # or R CMD INSTALL origin
library(origin)
```

`R CMD check` passes clean. After editing any roxygen block, regenerate the
generated files:

```r
roxygen2::roxygenise("origin", roclets = c("namespace", "rd"))
```

## analysis/ - benchmarking and analysis

`lineage_benchmark.R`, `physicell_visium.R`, `physicell_scdesign3.R`,
`scdesign3_helpers.R`, `symsim_profiles.R`, and `engine_comparison.R`, plus
their drivers in `analysis/cli/`. These are sourced rather than installed, deliberately: they
carry the heavy dependencies (Seurat, scDesign3, zellkonverter,
SingleCellExperiment, ggplot2) and they are still moving quickly. Keeping them
out of the package keeps those dependencies optional.

They are also the tier the clique test harness patches at runtime, which only
works for functions in the global environment. See below.

## legacy/ - the sim5 timestep engine

`sim5_code.R` and its satellites, unchanged apart from their paths. Still
driven by `bash_wrapper_all_combos.sh`. It depends on the core through
`load_origin.R`; nothing depends on it. `sim5_code.R` sets the working
directory to the repository root, so its `output/...` paths are unchanged.

## load_origin.R - the source-based entry point

Sourcing the simulator puts its functions in the global environment, where a
caller can replace them. A package namespace does not allow that: internal
calls resolve to the namespace binding and an override in the caller's
environment is silently ignored.

The clique test harness (`clique_2025_12_10/rust_cmd/test_harness/`) depends on
exactly that behavior - it replaces nine simulator functions to inject the
WT-CRISPR/FLARE recorder. So it loads through `load_origin.R`, not
`library(origin)`:

```r
source(file.path(source_root, "load_origin.R"))
```

Use `library(origin)` for everything that does not need to override simulator
internals. Turning those nine hooks into real extension points is what would
let the harness move to the package too; until then both entry points are
supported and load the same code.

```r
origin_include_analysis <- FALSE   # core only, skipping analysis/
source("load_origin.R")
```

## Verified after the split

- `R CMD check origin` - Status: OK, no notes or warnings
- `Rscript tests/regression_tests.R` - all regression tests pass
- all command-line drivers respond to `--help`
- `Rscript tests/symsim_tests.R` - SymSim integration tests (skips cleanly when
  SymSim is not installed)
- the clique harness loads and all 9 of its override targets remain patchable
