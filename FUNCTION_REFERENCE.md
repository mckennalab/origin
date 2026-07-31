# Function reference

This reference documents every named function in the tracked R, shell, and
notebook source. Functions defined inside another function are marked
**internal**. Anonymous callbacks passed directly to `apply`, `lapply`,
`vapply`, and related functions are implementation details and are not listed.

## Shared conventions

- Nucleotides are ordered `A, G, C, T` and encoded as integers `1, 2, 3, 4`.
- Mutation matrices use `0` for reference, `-1` for deletion, `1`–`4` for a
  substitution, and a decimal encoding for an insertion.
- A “profile” is generally a mutation matrix whose rows are mitochondrial
  genomes or barcode integrations and whose columns are sequence positions.
- Many simulator helpers use objects initialized by `sim5_code.R`; those are
  script-level helpers, not a stable package API.

## `sim5_code.R`

### Parsing, sequence construction, and rate setup

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `process_cla_string(cla_string, outputted_type)` | Semicolon-delimited value and requested output type. | Removes spaces, splits the value, and returns character, numeric, or integer values. |
| `parse_target_config(config_str)` | Target layout such as `U`, `R`, or `S:1:4`. | Returns `config`, `first_targ_pos`, and `bases_btwn`; the latter two are `NA` outside spaced mode. |
| `parse_target_count_arguments(num_targets, class_fracs)` | Target count and high/medium/low weights. | Validates and normalizes the weights, then returns integer `num_h`, `num_m`, and `num_l` counts that sum to `num_targets`. |
| `parse_be_example(example_string)` | Base-editor pattern such as `A --> G`. | Returns `from_base` and `to_base`; invalid or self-conversion patterns raise an error. |
| `classify_be_mutation_type(from_base, to_base)` | Two distinct A/C/G/T bases. | Returns `transition` or `transversion`. |
| `split_inds_into_hml(inds, num_h, num_m, num_l)` | Target indices and class counts. | Randomly assigns indices to edit-rate classes and returns a named position-to-class list. |
| `gcd(a, b)` | Two integer or finite-decimal cycle lengths. | Returns their greatest common divisor while retaining decimal precision. |
| `gcd_multiple_vals(...)` | Two or more values accepted by `gcd`. | Reduces all values to one greatest common divisor. |
| `convert_int_to_nuc(int_val)` | Integer nucleotide code. | Returns its A/G/C/T character using the script’s lookup table. |
| `convert_nuc_to_int(nuc_val)` | A/G/C/T character. | Returns its integer code as stored in the script’s lookup table. |
| `sequence_nucleotide_fractions(sequence)` | Non-empty A/G/C/T vector. | Returns named A/G/C/T fractions, including zero for absent bases. |
| `create_bc_sequence(...)` | Barcode length, target counts/classes/layouts, optional fixed sequence path, and nucleotide fractions. | Returns the barcode sequence plus named BE and nuclease position-to-class lists. A fixed sequence may be plain text or FASTA. |
| `create_prime_editing_basepos_seqs(target_inds, guide_length, num_unique_guides)` | Target positions and guide-library dimensions. | Returns position-keyed integer guides and correctly translated nucleotide guide strings. |
| `parse_sub_model_params(raw_cla_submodel, selected_sub_model, sequence_with_targets)` | Model name (`JC`, `K80`, `K81`, `F81`, `HKY`, or `GTR`), model parameters, and sequence. | Returns the parsed parameter list and the model’s 4-by-4 substitution-probability matrix. |
| `generate_transition_basepos_list(sequence_with_targets, sub_prob_mat)` | Sequence and A/G/C/T substitution matrix. | Returns one transition probability per sequence position. |
| `generate_transversion_basepos_list(sequence_with_targets, sub_prob_mat)` | Sequence and A/G/C/T substitution matrix. | Returns two named destination-base probabilities per sequence position. |
| `drop_editrate(rate, num_degrees)` | `High`, `Medium`, or `Low` and number of decay steps. | Lowers an edit-rate class; returns `FALSE` once decay falls below `Low`. |
| `get_new_be_targets(...)` | BE editing-window specification, target classes, barcode sequence, and conversion type. | Expands BE targets into their editing windows and returns expanded classes plus window-to-target mappings. |
| `get_new_nuc_targets(...)` | Nuclease editing-window specification, target classes, and barcode sequence. | Expands nuclease targets and returns expanded classes plus window-to-target mappings. |
| `estimate_prob_per_timept(prob_event_per_cell_cycle, timepoints_per_cell_cycle)` | Per-cycle probability and number of simulation steps per cycle. | Returns the per-step probability that preserves the requested per-cycle probability. |

### Heteroplasmy and population dynamics

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `draw_severity_scores(num_draws, mean1, mean2, sigma, mean2_weight)` | Mixture distribution parameters. | Draws variant severity scores from a two-component Gaussian mixture. |
| `heteroplasmy_one_site(...)` | One seeded site plus penetrance, genome, substitution, and fitness parameters. | Returns sparse-matrix triplets, variant fractions, severity scores, and counts for that site. |
| `make_dt(genomes, base)` | **Internal to heteroplasmy initialization.** Genome rows and a mutated base. | Builds an `(i, j, x)` tibble for one mutation class at the current heteroplasmy position. |
| `logistic_prob_survive_given_score(this_score, init_score, init_heteroplasmy_survive_prob, beta)` | Cell score and logistic anchor parameters. | Returns the cell’s survival probability. |
| `old_cells_at_timept(timept, cc_length)` | Timepoint and cell-cycle length. | Legacy lower-bound calculation retained for compatibility; it is not used by the current simulation. |
| `reassign_genome_inds(mito_to_genome_list)` | Mitochondrion-to-genome index list. | Renumbers flattened genome indices consecutively while preserving mitochondrial groups. |
| `mito_dynamics(...)` | Mitochondrial map/profile and fusion, fission, inheritance, and dropout settings. | Applies fusion/fission, replication, random daughter allocation, and post-mitotic dropout; returns daughter profiles and reindexed maps. Only `inheritance_pattern = "random"` is implemented. |
| `generate_downsample_cells(cell_sample_rate_vec, sim_length_stopping_points, cell_cycle_length)` | Cell sampling rates and simulation timing. | Returns a data frame of expected terminal population sizes and sampled cell indices. |
| `generate_downsample_integrations(max_ints_per_cell_vec, recovery_rate_vec)` | Integration counts and recovery rates. | Returns a data frame containing sampled integration indices for every parameter combination. |
| `setup_sim(...)` | Full initialized simulation state and parameter collections. | Creates the founder population, parallel worker cluster, worker exports, and division schedule; returns the initial cell population. |
| `get_future_div_points(sim_length, cc_length, current_timepoint)` | **Internal to `setup_sim`.** Time horizon, mean cycle length, and current time. | Draws future division times from an exponential waiting-time model, truncated at the horizon. |
| `get_initial_induction(induction_list)` | **Internal to `setup_sim`.** Induction time/count/fraction specification. | Returns whether the founder starts induced. |
| `multi_core_func(...)` | One timepoint, current population, mutation-rate collections, transition matrices, and output controls. | Divides eligible cells, applies differentiation/editing induction, death, mt/barcode mutation, and stopping-point output; returns the updated population. |

### Stopping-point output

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `all_processes_at_stopping_point(timept_savename, relative_timepoint, this_endpoint, all_recon_methods)` | Endpoint labels and requested reconstruction methods. | Coordinates timing output, population snapshots, ground-truth trees, sampled profiles, FASTA/score matrices, and reference files. |
| `describe_mutation_process_timing(...)` | **Internal.** Possible times and mt/barcode timing vectors. | Produces a one-row timing data frame for the current endpoint. |
| `save_mutation_profiles(mt_profiles, bc_profiles)` | **Internal.** Endpoint mt and barcode profile lists. | Writes non-null profile lists under `output/mut_profiles/<run_id>/`. |
| `create_lineage_strings(cell_lineage)` | **Internal.** Named lineage records. | Returns lineage strings used as tree tip/node identifiers. |
| `create_ground_truth_tree(cell_population, urid, save_path_stem)` | **Internal.** Population and output labels. | Constructs the lineage tree and writes its Newick representation. |
| `get_parent(x)` | **Internal to `create_ground_truth_tree`.** One cell identifier. | Returns that cell’s parent identifier. |
| `create_modified_profile_lists(...)` | **Internal.** Population, modality, recovery, integration, and FASTA controls. | Selects cells/integrations and returns profile lists plus recovered integration/UMI metadata. |
| `write_reference_fastas(...)` | **Internal.** Modality, reference sequences, and output names. | Writes mt/barcode reference FASTA files used by alignment and visualization tools. |
| `join_endpoint_results(unique_run_id)` | Run identifier. | Reads all per-endpoint result/specification CSVs, row-binds them, and writes one merged CSV. |
| `make_lineplot(run_id, poss_recon_modals, save_plots)` | Run identifier, modalities, and save flag. | Builds endpoint/cell-count line plots and optionally writes them under `output/lineplots/`. |

## `nonuniform_muts_heterogeneous.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `filter_elig_ints_by_edit_window(pos_to_window_inds_list, window_to_pos_inds_list, pos_to_unedited_int_list)` | Position/window mappings and unedited integration indices. | Returns position-keyed integrations that remain eligible after enforcing closed editing windows. |
| `non_uniform_editing(pos_er_list, num_integrations, eligible_ints, timepoint_savename, length1_positions)` | Position-specific rates and eligible integrations. | Samples unique mutation coordinates and returns `i_coords`/`j_coords`, or `FALSE` sentinels when no edit occurs. |
| `get_background_edit_inds(num_rows, num_cols, bg_pos_er_list, mut_type, sample_transversion, verbose)` | Matrix dimensions and background position rates. | Samples background mutation coordinates; transversion mode also returns the selected destination base per position. |
| `transition_func(...)` | Profile matrix, reference sequence, background/target transition rates, and optional editing-window controls. | Applies background and target transitions and returns the updated sparse profile. |
| `post_indices_transition_func(...)` | **Internal to `transition_func`.** Coordinates, matrix, and transition lookup. | Converts sampled coordinates into mutation values, including transitions within insertions. |
| `transversion_func(...)` | Profile matrix, reference, substitution matrix, background/target rates, optional forced destination base, and window controls. | Applies background and target transversions and returns the updated sparse profile; forced BE transversions use the configured conversion destination. |
| `post_indices_transversion_func(...)` | **Internal to `transversion_func`.** Coordinates, matrix, destination bases, and force flag. | Converts sampled transversion coordinates into sparse-matrix updates, including edits within insertions. |
| `insertion_func(...)` | Profile matrix, background/target insertion rates, prime-editing data, and window controls. | Applies ordinary or prime-editing insertions and returns the updated profile. |
| `ins_in_ins(ins_pos, current_ins, new_ins_length)` | **Internal to `insertion_func`.** Existing encoded insertion and new insertion placement. | Returns the numeric delta needed to encode an insertion inside another insertion. |
| `post_indices_insertion_func(...)` | **Internal to `insertion_func`.** Coordinates, incoming matrix, and prime-editing data. | Builds insertion values and adds them to the incoming sparse matrix. |
| `num_deletable_bases(x)` | One encoded mutation-matrix value. | Returns how many reference/inserted bases at that position can still be deleted. |
| `perform_deletion(ival, jval, del_length, mat_name, num_cols)` | Starting coordinate, deletion length, and profile matrix. | Recursively extends a deletion leftward and returns the modified matrix. |
| `all_deletions_one_mat(i, j, d, old_mat, num_cols)` | Vectors of deletion starts/lengths and a matrix. | Applies all requested deletions sequentially. |
| `deletion_func(...)` | Profile matrix, background/target deletion rates, dropout settings, and window controls. | Applies deletions and optional inter-target dropout, returning the updated profile. |
| `multi_edit_bc_dropout(...)` | **Internal to `deletion_func`.** Simultaneous deletion coordinates, radius/probability, and barcode profile. | Probabilistically deletes sequence between nearby deletion events. |
| `find_deletions_within_radius(delmat, radius)` | **Internal to `multi_edit_bc_dropout`.** Deletion coordinates and radius. | Returns pairs of deletion positions close enough to permit intervening dropout. |
| `post_indices_deletion_func(...)` | **Internal to `deletion_func`.** Deletion coordinates and incoming matrix. | Draws deletion lengths and delegates to `all_deletions_one_mat`. |
| `perform_all_mt_mutations(...)` | Mitochondrial profile, four background-rate lists, and substitution matrix. | Applies transition, transversion, insertion, and deletion processes in order. |
| `perform_all_bc_mutations(...)` | Barcode profile, background/target rates, substitution matrix, prime-editing, forced-transversion, dropout, and window controls. | Applies all four barcode mutation processes in order. |

## `physicell_lineage.R`

### Lineage import and tree construction

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `physicell_id(values, field_name)` | PhysiCell ID values and an error-label name. | Validates non-negative integer IDs and returns stable, non-scientific character IDs. |
| `validate_physicell_divisions(divisions)` | Data frame expected to contain `time`, `parent_ID`, and `daughter_ID`. | Returns a time-sorted canonical event table; rejects invalid times, self-parenting, and reused daughters. |
| `read_physicell_divisions(path)` | PhysiCell division-event CSV path. | Reads and validates the event table. |
| `read_physicell_terminal_ids(path)` | PhysiCell live-cell CSV path. | Validates the required `ID` column and returns unique live IDs. |
| `build_physicell_lineage(divisions, end_time, founder_time, founder_ids)` | Division events, final sampling time, founder start time, and optional explicit founders. | Converts persistent parent IDs into an event-resolved binary node table with exact branch durations. |
| `add_node(...)` | **Internal to `build_physicell_lineage`.** Cell/parent IDs, time, origin, and event index. | Appends one uniquely named branch segment to the growing node table. |
| `add_founder(cell_id)` | **Internal to `build_physicell_lineage`.** Previously unseen PhysiCell ID. | Creates and activates one founder segment. |
| `escape_newick_label(label)` | Arbitrary node/tip label. | Returns a Newick-safe unquoted or single-quoted label. |
| `physicell_lineage_to_newick(nodes, terminal_physicell_ids)` | Event-resolved node table and optional sampled/live IDs. | Returns a synthetic-root Newick string; prunes unrequested tips while retaining unary nodes and elapsed branch time. |
| `render_node(node_id)` | **Internal to `physicell_lineage_to_newick`.** One event-resolved node ID. | Recursively renders the requested node subtree. |

### Barcode model preparation

| Function | Inputs | Result |
| --- | --- | --- |
| `largest_remainder_counts(total, weights)` | Integer total and non-negative allocation weights. | Integer counts summing exactly to the total. |
| `physicell_target_positions(target_spec, barcode_length)` | One JSON target block and barcode length. | Named target-position-to-H/M/L-class vector for uniform, random, or spaced layouts. |
| `expand_physicell_target_windows(target_classes, window_spec, barcode_sequence, same_base_only)` | Target classes, editing-window block, barcode, and BE/nuclease mode. | Expanded position classes and named editing-window memberships. |
| `lower_class(edit_class, distance)` | **Internal to `expand_physicell_target_windows`.** H/M/L class and distance from its target. | Same or decayed class, or `NA` after decay below `Low`. |
| `physicell_barcode_reference(params, be_positions, be_from, params_dir)` | Parsed JSON, BE positions/base, and parameter-file directory. | Reads a configured barcode or generates one with the requested composition and forced BE targets. |
| `physicell_model_parameters(values, expected_length)` | Vector/list/semicolon-delimited model parameters. | Fixed-length numeric vector preserving missing values as `NA`. |
| `physicell_substitution_matrix(model_name, model_parameters, sequence)` | JC/K80/K81/F81/HKY/GTR specification and barcode reference. | Validated A/G/C/T substitution-probability matrix. |
| `derive_three(parameters, model_label)` | **Internal to `physicell_substitution_matrix`.** Ratio, transition, and transversion values. | Derives the missing member when at least two are supplied. |
| `draw_physicell_target_rates(target_classes, mean_probability)` | Named H/M/L classes and mean per-division probability. | Gamma-distributed, class-stratified position rates capped below one. |
| `physicell_rate_set(params, cell_type, editing_state, barcode_sequence, be_classes, nuc_classes, be_to)` | Parsed JSON plus one cell type/editing state and target definitions. | Per-position competing event probabilities, cycle length, substitution matrix, and invariant positions. |
| `prepare_physicell_recording_model(params, cell_type, num_integrations, params_dir, seed)` | Parsed remote_mito JSON and adapter overrides. | Complete induced/uninduced barcode model, reference, targets, windows, induction time, and integration count. |
| `physicell_probability_hazard(probability, cell_cycle_length)` | Per-division probabilities and cell-cycle duration. | Named continuous-time hazards. |

### Branch simulation and output

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `mutate_physicell_barcode_segment(profile, duration, rate_set, model)` | One inherited profile, elapsed time, state-specific rates, and model. | Applies irreversible competing barcode edits and returns the profile plus newly acquired event rows. |
| `mutate_physicell_barcode_branch(profile, start_time, end_time, model, editing_state)` | Profile, branch interval, prepared model, and `auto`/forced editing state. | Splits at global induction when needed and returns the final profile plus events with interval provenance. |
| `simulate_recording_on_physicell_lineage(nodes, model, editing_state, terminal_physicell_ids, seed)` | Event tree, prepared model, state policy, optional live IDs, and seed. | Propagates inherited profiles in topological order and returns all node profiles, sampled terminals, and mutation events. |
| `physicell_profile_to_sequence(profile_row, reference)` | One encoded integration profile and barcode reference. | Reconstructs its variable-length nucleotide sequence. |
| `write_physicell_lineage_outputs(nodes, terminal_nodes, output_dir)` | Event-resolved node table, sampled terminal rows, and destination. | Writes common node/tip CSVs and full/sampled Newick trees for either recording modality. |
| `write_physicell_recording_outputs(simulation, model, output_dir)` | Completed imported-lineage simulation and destination. | Writes common tree files plus barcode events, RDS profiles, allele/binary matrices, FASTA files, and manifest. |

## `physicell_mito.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `physicell_mito_rate_set(params, cell_type, editing_state, mitochondrial_reference)` | Parsed JSON, selected cell type/rate state, and mt reference. | Per-position substitution/indel hazards, invariant positions, cumulative sampling weights, and cell-cycle scale. |
| `prepare_physicell_mito_model(params, cell_type, genomes_per_cell, seed)` | Parsed JSON and mitochondrial replay overrides. | Generated mt reference plus induced/uninduced rates, induction time, and fixed bottleneck size. |
| `initialize_physicell_mito_profile(model, params)` | Prepared model and parsed JSON. | Founder genome profiles and baseline heteroplasmy event rows drawn from the configured site fraction, Beta penetrance, and transition probability. |
| `inherit_physicell_mito_profile(parent_profile, genomes_per_cell)` | Parent mt genome list and daughter bottleneck size. | Daughter genome list resampled with replacement from the parent. |
| `sample_physicell_mito_position(rate_set, excluded_positions)` | Weighted rate set and already mutated positions. | One mutable position sampled in proportion to its total event hazard. |
| `mutate_physicell_mito_segment(profile, duration, rate_set, model)` | Sparse genome list, elapsed time, state-specific rates, and model. | Applies low-rate irreversible substitutions/indels and returns the updated profile plus event rows. |
| `mutate_physicell_mito_branch(profile, start_time, end_time, model, editing_state)` | Profile, branch interval, model, and `auto`/forced state. | Splits a branch at induction when needed and returns its final profile and provenance-annotated events. |
| `simulate_mito_on_physicell_lineage(nodes, model, params, editing_state, terminal_physicell_ids, seed)` | Event tree, prepared model, JSON, state policy, optional current-cell IDs, and seed. | Propagates bottlenecked mt profiles in topological order and returns all profiles, sampled terminals, and events. |
| `physicell_mito_variant_fractions(simulation)` | Completed mitochondrial replay. | Long-form per-terminal-cell position/allele counts and heteroplasmy fractions. |
| `physicell_mito_haplotype_sequence(genome, reference)` | One sparse encoded mt genome and reference. | Reconstructs its nucleotide sequence, including substitutions, deletions, and one-base insertions. |
| `write_physicell_mito_outputs(simulation, model, output_dir, write_fasta)` | Completed replay, prepared model, destination, and FASTA flag. | Writes mt profiles/events/fractions/reference, optional sparse matrix and sampled haplotypes, and a manifest. |

## `simulate_physicell_lineage.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `physicell_cli_usage()` | None. | Returns the command-line usage text. |
| `parse_physicell_cli_args(arguments)` | Trailing command-line arguments. | Validates supported short/long options, including modality/mt controls, and returns a typed options list; `--help` prints usage and exits. |

The script’s top-level code loads a remote_mito JSON, imports the PhysiCell
division/current-cell files, prepares the requested barcode and/or
mitochondrial models, simulates recording on the fixed tree, writes all
artifacts, and prints `OUTPUT_DIR=<path>`.

## `scdesign3_helpers.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `split_scdesign3_columns(value, expected_length)` | Optional comma-delimited CLI value and optional required count. | Validated unique syntactic covariate names. |
| `load_scdesign3_reference(path, celltype_col, pseudotime_col, use_pseudotime, spatial_cols, other_covariates)` | SCE/Seurat RDS or h5ad path plus reference column selections. | A `SingleCellExperiment` with `counts` and standardized cell-type/pseudotime/spatial fields. |
| `scdesign3_default_mu_formula(sce, use_pseudotime, use_spatial, other_covariates)` | Standardized reference and selected predictor classes. | Default marginal-mean formula containing applicable cell-type, smooth pseudotime/spatial, and additional terms. |
| `scdesign3_fit_cache_key(...)` | Reference identity, package/model version, and covariate/formula choices. | SHA-1 cache key when `digest` is installed or a deterministic fallback. |
| `fit_or_load_scdesign3(sce, cache_path, use_pseudotime, use_spatial, other_covariates, mu_formula, family_use, ncores)` | Standardized reference and fit specification. | Loads a compatible fit or runs corrected scDesign3 data construction, marginal fitting, and Gaussian-copula fitting before caching it. |
| `scdesign3_new_covariates(fit, metadata)` | Cached fit and target-cell metadata. | Reference-compatible target covariate frame with validated factor levels/numeric values and correlation group. |
| `simulate_scdesign3_counts(fit, metadata, ncores)` | Cached fit, target metadata with sample IDs, and workers. | Feature-by-target-cell counts from `extract_para()` and `simu_new()`. |

## `physicell_scdesign3.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `read_physicell_property_manifest(path)` | Optional two-column manifest CSV. | Named character property vector, or an empty vector when absent. |
| `physicell_node_depths(nodes)` | Event-resolved node/parent table. | Named generation depths; missing parents and cycles fail explicitly. |
| `physicell_neighbor_count(values)` | Semicolon-delimited PhysiCell neighbor-ID values. | Integer neighbor counts, treating blank/missing fields as zero. |
| `build_physicell_sc_covariates(recording_dir, lineage_table_path, cell_type)` | Recording output, final PhysiCell table, and optional type override. | Terminal metadata joined across lineage timing/depth, spatial/neighborhood, barcode burden, and mitochondrial heteroplasmy. |
| `apply_physicell_celltype_map(metadata, celltype_map)` | Target metadata and named simulated-to-reference mapping. | Metadata with mapped model type and retained original `cell_type_sim`. |
| `scale_physicell_covariate_to_reference(values, reference_values)` | Simulated and reference numeric covariates. | Simulated values min/max mapped into the reference range, with degenerate-range handling. |
| `align_physicell_sc_covariates(metadata, reference_sce, use_pseudotime, use_spatial)` | Joined target metadata and standardized reference. | Adds reference-range pseudotime and/or `spatial1`/`spatial2` predictors. |
| `write_physicell_scdesign3_outputs(counts, metadata, output_dir, fit, reference_path, seed)` | Synthetic counts, linked metadata, fit/provenance, and destination. | Writes a terminal-cell SCE, count RDS, metadata/summary CSVs, and manifest. |

## `generate_physicell_sc_profiles.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `physicell_scdesign3_usage()` | None. | Returns dependency-light CLI usage text. |
| `parse_physicell_scdesign3_args(arguments)` | Trailing CLI arguments. | Validates flags and returns typed PhysiCell/scDesign3 options; help prints usage and exits. |

The top-level CLI resolves a completed PhysiCell pipeline run, joins terminal
covariates, loads/maps the reference, fits or reuses scDesign3, simulates
counts, writes `lineage_recording/sc_profiles`, and prints
`SCDESIGN3_OUTPUT_DIR=<path>`.

## `substitution_models.r`

All returned matrices use A/G/C/T rows as source bases and A/G/C/T columns as
destination bases; diagonals are zero.

| Function | Inputs | Result |
| --- | --- | --- |
| `jc_sub_rate_mat(overall_sub_rate)` | Shared non-self substitution rate. | Jukes–Cantor 4-by-4 rate matrix. |
| `k80_sub_rate_mat(transition_to_transversion_ratio, transition_rate, transversion_rate)` | Any two of the three K80 parameters. | K80 matrix with one transition and one transversion rate. |
| `k81_sub_rate_mat(transition_rate, transversion_rate_weakstrong_conserved, transversion_rate_aminoketo_conserved)` | Three K81 rates. | K81 matrix with one transition and two transversion classes. |
| `f81_sub_rate_mat(frac_a, frac_g, frac_c, frac_t, baseline_overall_sub_rate)` | Nucleotide fractions and baseline rate. | F81-style destination-frequency-weighted matrix normalized to the requested mean. |
| `hky_sub_rate_mat(...)` | Nucleotide fractions plus any two of ratio, transition rate, and transversion rate. | HKY matrix with destination-frequency weights and transition multipliers on A↔G and C↔T. |
| `gtr_sub_rate_mat(...)` | Six exchange rates and four nucleotide fractions. | General time-reversible matrix. |
| `SIMPLIFY_target_site_gamma_based_sub_rates(...)` | Sequence length, H/M/L target positions, and gamma bootstrap settings. | Named target-position rate list; returns an empty list when no targets exist. |
| `nontarget_get_invariant_inds(eligible_invariant_sites, frac_invariant)` | Eligible positions and invariant fraction. | Randomly selected invariant indices. |
| `nontarget_scale_gamma_heterogeneity(...)` | Position-rate list and gamma/discretization settings. | Copy of the rate list with independently sampled multiplicative heterogeneity; `shape_param = 0` leaves rates unchanged. |

## `add_intervening_be_targets_to_seq.r`

| Function | Inputs | Result |
| --- | --- | --- |
| `generate_non_be_target_sequence(barcode_length, nuc_fracs, target_from, be_target_count)` | Final length, A/G/C/T weights, edited base, and reserved target count. | Shuffled non-target backbone with exact length and largest-remainder nucleotide counts; target reservations take priority over composition. |
| `generate_target_indices(config, num_targets, target_pos_1, bc_length_with_targets, num_bases_btwn)` | `U`, `R`, or `S` layout and dimensions. | Validated target indices; impossible layouts raise an error. |
| `add_intervening_be_targets(...)` | Target layout/base/count and a non-target backbone. | Returns `target_inds` and the completed `seq_with_targets`. |

## `mut_to_fasta_difflen_ints.r`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `get_profiles_ints_and_umis(cell_pop, bc_or_mt, int_rec_prob, num_ints, umis)` | Population, modality, recovery probability/count, and optional UMIs. | Per-cell list of sampled mutation matrices, recovered integration indices, and optional recovered UMIs. |
| `ins_to_charvec(ins, pos_num, ref_seq, fixed_length)` | Encoded insertion, position/reference, and fixed-length mode. | Returns the number and character vector of decoded inserted/reference bases. |
| `get_one_cell_sequence(cell_int_mat, ref_seq, these_bc_int_umis)` | One cell’s integrations/genomes, reference, and optional UMIs. | List of reconstructed nucleotide strings, one per recovered row. |
| `write_all_cell_sequences(...)` | Named cell matrices, reference, output path/type, and optional UMIs. | Reconstructs sequences in parallel and writes a FASTA file with the requested type suffix. |

## `mut_to_scoremat.r`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `group_deletions(deletion_df)` | Sorted or unsorted deletion-event data frame. | One row per contiguous cell/integration deletion run; the mutation value is `d<length>`. |
| `score_mat_to_phylip(score_mat, output_phylip_path)` | Score matrix and destination path. | Writes sequential PHYLIP text. |
| `new_scoremat_to_fasta(scoremat, output_fasta_path)` | Score matrix and destination path. | Writes character-state FASTA, representing missing values as `?`. |
| `new_create_one_score_mat(...)` | Profiles, recovery metadata, modality, collapse/binarization/allelic-fraction controls, and output labels. | Main score-matrix builder; returns heteroplasmy counts in count mode or writes matrix/FASTA artifacts in output mode. |
| `write_mat(mat, suffix)` | **Internal to `new_create_one_score_mat`.** Matrix and filename suffix. | Adds row/column labels and writes RDS plus character-state FASTA artifacts. |
| `create_one_score_mat(...)` | Legacy profile-to-score-matrix parameters. | Older score-matrix implementation retained for compatibility; writes barcode or mt matrices and FASTA files. |
| `get_norm_cell_heteroplasmy_scores(...)` | Per-cell mutation counts, genome counts, variant severity map, and normalization/fitness settings. | Returns named per-cell heteroplasmy severity scores and may add new variant scores to the shared severity map. |

## `parse_cell_type_specific_args.r`

| Function | Inputs | Result |
| --- | --- | --- |
| `make_cell_type_transition_lists(cell_type_names)` | Cell-type names; transition matrices are read from script-level `input_args`. | Returns nested named uninduced and induced transition lists; incompatible dimensions raise an error. |

## `fit_plot_parameters.R`

| Function | Inputs | Result |
| --- | --- | --- |
| `get_heatmap_params(num_cells)` | Population size. | Uses splines fitted at source time to return label position, height, font size, and image dimensions. |

Sourcing this file also reads `imported_heatmap_plotval_dat.csv`, fits the
curves, and writes `diagnostic_fitvals_plots.png`.

## `generate_sc_profiles_from_bash.r`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `load_reference(path, celltype_col, pseudotime_col, use_pseudotime)` | SCE/Seurat RDS or h5ad path and metadata column choices. | Returns a `SingleCellExperiment` with standardized `counts`, `cell_type`, and optional `pseudotime`; unsupported data raises an error. |
| `list_cell_population_files(cell_pop_dir)` | Run’s population directory. | Returns population snapshot paths ordered by parsed numeric timepoint. |
| `build_covariate_df(cellpop_files, celltype_map)` | Snapshot table and optional simulator-to-reference type mapping. | Returns one metadata row per living terminal/leaf cell at every stopping point, with lineage `cell_id` and unique cell-timepoint `sample_id`. |
| `add_pseudotime(meta)` | Cell metadata with lineage depth. | Adds depth normalized by the maximum depth. |
| `maybe_downsample(meta, max_per_tp)` | Metadata and optional per-timepoint cap. | Returns all rows or a random capped subset per timepoint. |
| `ref_cache_key(ref_path, celltype_col, use_pseudotime, pseudotime_col)` | Reference path and fit-shaping options. | Delegates to the versioned shared scDesign3 cache-key helper. |
| `fit_or_load(sce, cache_path, use_pseudotime, ncores)` | Reference SCE, cache path, model choice, and cores. | Delegates to the corrected shared marginal/copula fit and cache implementation. |
| `simulate_for_meta(fit, new_meta, ncores)` | Cached fit and simulator-derived covariates. | Delegates to shared conditional count simulation; unknown factor levels fail explicitly. |

## `compare_trees_call_from_bash.r`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `save_image_of_tree(tree_path, savename_prefix, cell_pop_path, midpt)` | Tree/output labels and optional midpoint-root flag. | Reads and optionally roots a tree, then writes a PNG. `cell_pop_path` is retained for call compatibility. |
| `plot_tree_with_color(...)` | Tree, cell population, output label, tree role, and color controls. | Writes a fan-tree PDF with tips colored by cell type and optional internal-node colors/branch labels. |

The script’s top-level code also locates the matching ground-truth tree,
subsets it to reconstructed tips, computes normalized RF distance, and writes
the metric plus source parameter path.

## `process_results_from_bash.r`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `extract_subrun_details(file_name, bc_or_mt, mat_colnames, recon_method)` | RF filename, modality, requested columns, and method. | Parses timepoint, integration/recovery values, distinct cell-type sampling fractions, collapse, AF, and binarization metadata into a vector. |
| `change_transition_mat_colnames(flattened_mat)` | One-row flattened JSON parameter matrix. | Replaces numbered transition-matrix columns with source/target cell-type names. |
| `make_results_df(urid)` | Run identifier. | Combines recognized RF files with flattened JSON parameters, writes `stacked_results.csv`, and returns the data frame. |
| `split_vals(joined_val, outputted_type)` | **Internal to `make_results_df`.** Colon/semicolon-delimited string and type. | Splits and coerces legacy compound parameter values. |
| `make_heatmap(run_id)` | Run identifier. | Legacy plotting helper for an older results-column schema; it is not called by the current CLI path. |

## `make_babette_tree.r`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `extract_state_file_path(beast2_output)` | Parsed babette/BEAST output. | Extracts the state-file path from the expected output line. |
| `extract_newick_data(path_to_beast_output, newick_path)` | BEAST output and destination. | Extracts the first Newick-like tree string and writes it with a semicolon terminator. |
| `phylo_obj_from_newick(path_to_newick)` | Newick path. | Returns an `ape::phylo` object. |
| `fasta_to_phylo(...)` | FASTA, run/output labels, inference/site-model settings, and optional BEAST options. | Runs babette/BEAST, selects posterior tree output, builds a consensus tree, relabels tips, writes Newick files, and optionally returns the phylogeny. |

This module is currently optional; its `source()` call in `sim5_code.R` is
commented out.

## Shell functions

### `make_igv_input.sh`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `usage()` | None. | Prints CLI usage and exits with status 1. |
| `fasta_to_fastq(input_fasta_path, output_fastq_path)` | FASTA and FASTQ paths. | Converts single- or multi-line FASTA records to four-line FASTQ records with uniform `I` quality characters, using POSIX `awk`. |

### `run_physicell_10000_pipeline.sh`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `usage()` | None. | Prints the end-to-end PhysiCell pipeline options and defaults. |

The rest of this script is top-level orchestration: it validates inputs, stages
and patches a private PhysiCell project copy, compiles/runs it, checks lineage
outputs, launches dual-modality recording replay, and writes a pipeline
manifest.

`bash_wrapper_all_combos.sh` contains no named shell functions; it is a
top-level orchestration script.

## Notebook helpers

`parse_rf_results.ipynb` contains notebook-local versions of
`make_results_df(urid)`, its internal `split_vals(...)`, and
`make_heatmap(run_id)`. They perform exploratory RF aggregation and plotting
inside the notebook and are separate from the maintained functions in
`process_results_from_bash.r`.

## Regression-test helpers

`tests/regression_tests.R` defines:

- `load_named_function(relative_path, function_name, envir)` to load one
  top-level function without executing a CLI script;
- `expect_true(value, message)`, `expect_equal(actual, expected, message,
  tolerance)`, and `expect_error(expression, message)` as dependency-free test
  assertions;
- `collect_calls(expression, function_name)` to locate named calls in parsed R
  syntax for static argument-wiring checks; its internal `walk(node)` performs
  the recursive traversal;
- `top_level_assignment(name)` to retrieve one parsed top-level simulator
  assignment for static formal/argument comparisons.
