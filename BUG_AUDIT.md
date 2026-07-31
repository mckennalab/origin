# Bug audit

Audit date: 2026-07-25

The audit covered every tracked R and shell source file, the example JSON, the
main simulator’s function-call wiring, the PhysiCell fixed-tree adapters, and
focused executions of pure helpers. `tests/regression_tests.R` contains
dependency-light reproductions for the highest-risk fixes.

## Corrected defects

| Area | Defect and impact | Correction |
| --- | --- | --- |
| Main startup | `sim5_code.R` changed into a cluster-specific hard-coded directory, so it could not run from this checkout. | Resolve the project root from the script path. |
| Main startup | The simulator sourced `fit_plot_parameters.R`, which immediately read a calibration CSV absent from the repository even though none of its functions were used. | Treat the file as a standalone legacy utility and do not source it in the simulator. |
| Main call wiring | `setup_sim()` required `cold_startup` and `cell_population`, but `sim_arglist` never supplied them; the run failed before simulation. | Remove the unused required formals and add a regression check that `sim_arglist` exactly matches the function signature. |
| Main call wiring | The sole `multi_core_func()` call omitted score-matrix controls and the initial mt-genome count. | Supply every formal and statically verify exact argument coverage. |
| Parallel startup | `library("Matrix", "data.table")` treated the second package as `lib.loc`, leaving workers without `data.table`. | Load both packages in an explicit worker expression. |
| Barcode construction | When BE and nuclease targets were both enabled, the nuclease branch regenerated and overwrote the barcode after BE targets had been inserted. | Generate the nuclease-only backbone only when no barcode sequence already exists. |
| Barcode composition | A zero-BE-target barcode ignored requested nucleotide fractions; the rounding-adjustment signs also moved counts in the wrong direction. | Use validated largest-remainder allocation and reserve target bases before returning the shuffled backbone. |
| Target layout | Invalid configurations could leave `all_inds` undefined or allow spaced targets beyond the barcode. | Validate count, length, mode, position, and gap; fail explicitly on impossible layouts. |
| Target classes | Independent rounding of high and medium counts could make the low count negative. | Allocate all three classes with a largest-remainder calculation. |
| Editing windows | The BE close-after-edit option was read from the nuclease configuration. | Read `be_targets$editing_window$close_after_edit`. |
| Editing windows | Expanded positions were stored under target-position keys while the target itself was stored under `be_window_*`/`nuc_window_*`; reverse mappings therefore described separate, incomplete windows. | Keep each target and all expanded positions under one stable window name. |
| Editing windows | An empty eligible set at one position was discarded before intersection, reopening a window that should have been closed. | Retain empty sets in the intersection so any saturated position closes the whole window. |
| Editing windows | Direct zero-width calls returned uninitialized `growing_window_editrates`. | Initialize returned collections before branching. |
| Substitution setup | F81/HKY/GTR received raw nucleotide counts instead of fractions, and absent bases produced `NULL`. | Compute a complete named A/G/C/T fraction vector once. |
| HKY | The transition multiplier was applied to incorrect cells (mostly destination G) and multiplied an already-transition-scaled baseline a second time. | Build destination-weighted rates with transition multipliers only on A↔G and C↔T. |
| Heterogeneity | `shape_param = 0` replaced all background mutation rates with zero; a degenerate gamma draw replaced rates instead of scaling them. | Interpret zero shape as “no heterogeneity” and always apply nonzero draws multiplicatively. |
| Target-rate generation | The no-target test checked `length(total_num_targets)`, which is always one for a scalar, then attempted an invalid gamma cut. | Test `total_num_targets == 0` and return an empty list. |
| Prime editing | Integer guide translation returned the entire lookup table for every base, producing the same malformed string for every guide. | Translate each integer through its individual lookup entry. |
| Target transversions | The nonuniform transversion path omitted required destination bases and crashed; forced conversions also did not use `be_conversion_pattern`. | Build the destination vector and pass the configured BE destination base through simulator and worker calls. |
| Transversions in insertions | A list slice (`[`) was used where a numeric transversion vector (`[[`) was required, causing invalid matrix subscripting. | Extract the numeric option vector with `[[`. |
| Mitochondrial invariants | Only position 1 was eligible because the code used `length(mito_genome_length)` instead of the length value. | Use `seq_len(mito_genome_length)`. |
| Founder mt profiles | The mitochondrion map used Poisson genome counts, while the profile and heteroplasmy setup used unrelated fixed/max row counts. Map indices and profile rows could diverge. | Set the founder genome count and profile dimensions to the realized sum of per-mitochondrion genome counts. |
| Mitochondrial dropout | Dynamic list removal used the literal field `$mito_num`, so selected mitochondria remained in daughter maps after their profile rows were removed; surviving map indices were also stale after row deletion. | Remove selected dynamic list names and reindex every surviving map to the subsetted profile. |
| Mitochondrial fission | Binomial fission could create an empty daughter mitochondrion. | Condition the split so both products receive at least one genome. |
| Heteroplasmy fitness | Severity-map keys used the loop ordinal rather than the actual mt genomic position, so later mutation names did not match their seeded scores. | Key variants with `heteroplasmy_inds[ind_num]`. |
| Heteroplasmy helper | A zero-variant branch called the misspelled `iteger()`. | Use `integer()`. |
| Score matrices | Deletion adjacency had a misplaced parenthesis, so all deletions in one cell/integration became one run; boundary rows were also attributed to the wrong run. | Reimplement run grouping over sorted cell/integration/position rows. |
| Score matrices | Missing-integration masking checked the entire recovery list instead of the current cell, and skipped the case of exactly one missing feature. | Check the per-cell vector and mask whenever at least one feature is missing. |
| Result aggregation | Induced transition columns were renamed using the uninduced loop length. | Iterate over induced columns and validate both matrices against all cell-type pairs. |
| Result aggregation | Every parsed cell type received the first sampling fraction from the filename. | Parse the contiguous `cell-type-rate` block in one pass. |
| Result aggregation | `make_results_df()` returned the final `print()` value rather than its data frame. | Return the written data frame explicitly. |
| Tree comparison | The `--treefile_dir_path` option placed `type`, `default`, and `help` inside the flag vector, leaving the option unset. | Correct the `optparse::make_option()` call and validate required options. |
| Tree comparison | A caller-supplied ground-truth tree never initialized `full_gt_path`, but plotting later required it. | Normalize and retain the supplied path. |
| Tree comparison | Ambiguous snapshot/tree matches and reconstructed tips absent from ground truth failed later with unclear errors. | Require unique matches and validate the tip set before RF distance. |
| Tree plotting | Unknown cell types returned `NULL`, shortening the color vector; palettes with fewer than three or more than eight types were unsafe. | Use a gray fallback and size/interpolate the palette correctly. |
| Wrapper | Boolean checks executed `true`/`false` as commands; simulations that failed still post-processed the previous run. | Use explicit shell comparisons and skip post-processing on failure. |
| Wrapper | Selecting the newest output directory raced with concurrent simulations. | Parse the simulator’s machine-readable `UNIQUE_RUN_ID=<id>` line. |
| Wrapper | Unquoted paths, fragile `ls` loops, the wrong `.phy` suffix removal, and unconditional post-processing made edge cases fail. | Use arrays/nullglob, quote paths, remove `.fasta`, and guard optional output stages. |
| IGV utility | FASTA-to-FASTQ conversion assumed exactly two lines per record and used GNU-only `gensub`, which is unavailable in standard macOS `awk`. | Accumulate multiline records and generate qualities with POSIX `awk`. |
| Single-cell profiles | The same lineage appears at multiple stopping points, producing duplicate SCE column/row names. | Add a unique cell-timepoint `sample_id` while retaining lineage `cell_id`. |
| Single-cell profiles | Historical parents remain `alive` but are no longer terminal leaves; they were resampled at later snapshots. | Include only cells with both `alive` and `terminal` set. |
| Single-cell profiles | The scDesign3 adapter passed an unsupported `n_cores` argument and invalid `mclapply` value to `construct_data()`, so it failed against current scDesign3. | Use the 1.10 API and its supported `mcmapply` parallelization value. |
| Single-cell profiles | The complete `fit_copula()` return object was passed where `simu_new()` requires its `copula_list`; filtered-gene metadata was also omitted. | Retain the structured copula fit and pass `copula_list`, `important_feature`, and `filtered_gene` explicitly. |

## Known limitations not changed

These require a model or representation decision beyond a safe bug fix:

- Mutation-matrix insertions are encoded in floating-point decimal digits.
  Long insertions can lose digits to floating-point precision, and interactions
  between insertions and already deleted bases are difficult to represent
  unambiguously. A structured mutation representation would be safer.
- `max_mito_per_cell` is described as a saturation setting but is not enforced
  by `mito_dynamics()`. Defining whether the cap applies before replication,
  before daughter allocation, or after allocation is a biological-model choice.
- `process_results_from_bash.r::make_heatmap()` and the notebook copy target an
  older results-column schema. The maintained CLI writes the aggregation CSV
  but does not call this legacy plot helper.
- `fit_plot_parameters.R` depends on the untracked
  `imported_heatmap_plotval_dat.csv`; it remains a standalone legacy utility and
  is no longer on the simulation startup path.

## Verification

- `Rscript tests/regression_tests.R`
- A compiled end-to-end `tumor_3D_lineage` smoke run from one founder to eight
  current cells, followed by barcode and mitochondrial replay
- A live scDesign3 1.10 fit and conditional simulation using lineage
  pseudotime plus barcode edit fraction for the eight-cell PhysiCell run
- Parse every tracked `.R`/`.r` source with `Rscript::parse`
- `bash -n bash_wrapper_all_combos.sh make_igv_input.sh run_physicell_10000_pipeline.sh`
- Parse `example_json_params/short_test.json` and
  `example_json_params/physicell_10000.json` with `jsonlite`
- `git diff --check`

The full 10,000-cell simulation, BEAST workflows, MUSCLE, IQ-TREE, BWA, and
samtools were not executed in this environment. The PhysiCell/scDesign3 path
was validated at smoke-test scale.
