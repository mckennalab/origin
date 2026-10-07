# Function reference

This reference documents every named function in the tracked R, shell, and
notebook source. Functions defined inside another function are marked
**internal**. Anonymous callbacks passed directly to `apply`, `lapply`,
`vapply`, and related functions are implementation details and are not listed.

Each entry below is a one-line contract meant for scanning. The R source itself
carries the full per-function documentation as roxygen (`#'`) blocks giving the
mechanism, one `@param` per argument, the `@return` shape, and a
`@section Side effects:` wherever a function writes files, seeds the RNG, or
starts a worker cluster — read those alongside the code. The one exception is
`add_intervening_be_targets_to_seq.r`, which uses the `docstring` package and so
keeps its blocks *inside* the function bodies.

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
| `sample_induced_cells(cell_names, num_cells, frac_cells)` | Available cell identifiers and either a fixed count or fraction. | Samples one valid induction subset; fixed counts are capped at the available population. |
| `initialize_founder_population(...)` | Founder count, initial profiles/state, per-founder division schedules, and time-zero induction selections. | Returns one independently named population record per configured founder. |
| `setup_sim(...)` | Full initialized simulation state and parameter collections. | Creates the founder population, parallel worker cluster, worker exports, and division schedule; returns the initial cell population. |
| `get_future_div_points(sim_length, cc_length, current_timepoint)` | **Internal to `setup_sim`.** Time horizon, mean cycle length, and current time. | Draws future division times from an exponential waiting-time model, truncated at the horizon. |
| `multi_core_func(...)` | One timepoint, current population, mutation-rate collections, transition matrices, and output controls. | Divides eligible cells, applies differentiation/editing induction, death, mt/barcode mutation, and stopping-point output; returns the updated population. |

### Stopping-point output

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `all_processes_at_stopping_point(timept_savename, relative_timepoint, this_endpoint, all_recon_methods)` | Endpoint labels and requested reconstruction methods. | Coordinates timing output, population snapshots, ground-truth trees, sampled profiles, FASTA/score matrices, and reference files. |
| `describe_mutation_process_timing(...)` | **Internal.** Possible times and mt/barcode timing vectors. | Produces a one-row timing data frame for the current endpoint. |
| `save_mutation_profiles(mt_profiles, bc_profiles)` | **Internal.** Endpoint mt and barcode profile lists. | Writes non-null profile lists under `output/mut_profiles/<run_id>/`. |
| `create_lineage_strings(cell_lineage)` | **Internal.** Named lineage records. | Returns lineage strings used as tree tip/node identifiers. |
| `create_ground_truth_tree(cell_population, urid, save_path_stem, output_root)` | Population and output labels. | Constructs and writes the Newick lineage tree; multiple founders are joined beneath a synthetic time-zero root. |
| `create_modified_profile_lists(...)` | **Internal.** Population, modality, recovery, integration, and FASTA controls. | Selects cells/integrations and returns profile lists plus recovered integration/UMI metadata. |
| `write_reference_fastas(...)` | **Internal.** Modality, reference sequences, and output names. | Writes mt/barcode reference FASTA files used by alignment and visualization tools. |
| `join_endpoint_results(unique_run_id)` | Run identifier. | Reads all per-endpoint result/specification CSVs, row-binds them, and writes one merged CSV. |
| `make_lineplot(run_id, poss_recon_modals, save_plots)` | Run identifier, modalities, and save flag. | Builds endpoint/cell-count line plots and optionally writes them under `output/lineplots/`. |

## `nonuniform_muts_heterogeneous.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `filter_elig_ints_by_edit_window(pos_to_window_inds_list, window_to_pos_inds_list, pos_to_unedited_int_list)` | Position/window mappings and unedited integration indices. | Returns position-keyed integrations that remain eligible after enforcing closed editing windows. |
| `non_uniform_editing(pos_er_list, num_integrations, eligible_ints, timepoint_savename, length1_positions)` | Position-specific rates and eligible integrations. | Samples unique mutation coordinates and returns `i_coords`/`j_coords`, or `FALSE` sentinels when no edit occurs. |
| `get_background_edit_inds(num_rows, num_cols, bg_pos_er_list, mut_type, sample_transversion, verbose)` | Matrix dimensions and background position rates. | Independently samples each position at its configured probability without duplicate coordinates; transversion mode also returns the selected destination base per position. |
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

## `prime_editing.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `prime_editing_enabled(params)` | Parsed native parameter list. | Whether the unified backend or legacy prime-editing flag enables prime editing. |
| `prime_editing_probability_vector(value, number, name)` | Scalar/vector probability, target count, and diagnostic name. | Validated probability vector with one value per recorder target. |
| `prime_editing_dna_sequence(value, name, allow_empty)` | Candidate sequence and validation controls. | Uppercase whitespace-free A/C/G/T sequence or an explicit validation error. |
| `prime_editing_pool_frame(value)` | Data frame or JSON-style array of pegRNA objects. | Rectangular pegRNA pool with missing optional fields filled by `NA`. |
| `read_prime_editing_pool(params, params_dir)` | Parameters and parameter-file directory. | Validated known pegRNA pool loaded from embedded JSON or CSV; supports the legacy generated-pool fallback. |
| `assign_prime_editing_targets(pool, target_positions, configuration)` | Validated pool, target positions, and backend configuration. | One ordered pegRNA assignment per target using explicit, cyclic, or sampled assignment. |
| `prepare_prime_editing_backend(params, target_positions, params_dir, seed)` | Parameters, recorder positions, relative-path root, and seed. | Unified pool, target mapping, provenance, assignment mode, and categorical state encoding. |
| `prime_editing_scale_probability(base_probability, efficiency)` | Base per-cycle probability and pegRNA efficiency. | Efficiency-adjusted probability `1-(1-p)^e`, preserving multiplicative survival hazards. |
| `prime_editing_sequence_integer_map(backend)` | Prepared backend. | Legacy time-step position-to-integer-template map. |
| `prime_editing_sequence_character_map(backend)` | Prepared backend. | Exact position-to-edit-sequence map. |
| `write_prime_editing_backend_manifest(backend, path, write_csv)` | Prepared backend, output path, and optional CSV writer. | Writes the exact target/pegRNA mapping plus pool provenance and assignment mode. |

## `physicell_lineage.R`

### Lineage import and tree construction

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `physicell_id(values, field_name)` | PhysiCell ID values and an error-label name. | Validates non-negative integer IDs and returns stable, non-scientific character IDs. |
| `physicell_csv_output_path(path, compress)` | Base CSV path and compression flag. | Returns the canonical `.csv.gz` or uncompressed `.csv` output path. |
| `resolve_physicell_csv_path(path, required)` | CSV or CSV-gzip path and required flag. | Resolves compressed/uncompressed alternatives, preferring `.csv.gz`, or returns `NA` for an optional missing table. |
| `write_physicell_csv(x, path, row.names, compress, ...)` | Tabular object, base path, row-name/compression flags, and `write.csv` options. | Writes a gzip-compressed CSV by default and invisibly returns its normalized path. |
| `read_physicell_csv(path, required, ...)` | CSV or CSV-gzip path and `read.csv` options. | Transparently reads either form after resolving the available file. |
| `validate_physicell_divisions(divisions)` | Data frame expected to contain `time`, `parent_ID`, and `daughter_ID`. | Returns a time-sorted, validated canonical event table using linear-time vector matching; rejects invalid times, self-parenting, and reused daughters. |
| `read_physicell_divisions(path)` | PhysiCell division-event CSV path. | Reads and validates the event table. |
| `read_physicell_cell_ids(path, table_description, alive_only)` | PhysiCell cell-table path, error label, and alive-filter flag. | Validates the required `ID` column and optionally filters an `alive` field. |
| `read_physicell_terminal_ids(path)` | PhysiCell live-cell CSV path. | Returns unique live IDs, filtering `alive == false` when that column is present. |
| `read_physicell_founder_ids(path)` | PhysiCell founder CSV path. | Returns every explicit day-zero founder ID without alive filtering. |
| `build_physicell_lineage(divisions, end_time, founder_time, founder_ids, show_progress, progress_updates)` | Division events, final sampling time, founder start time, optional explicit founders, and progress controls. | Converts persistent parent IDs into an event-resolved binary node table in linear time using a hashed active-node map and preallocated columns. |
| `add_node(...)` | **Internal to `build_physicell_lineage`.** Cell/parent IDs, time, origin, and event index. | Fills the next preallocated branch row and returns its integer index. |
| `add_founder(cell_id)` | **Internal to `build_physicell_lineage`.** Previously unseen PhysiCell ID. | Creates and activates one founder segment in the hashed cell-to-node map. |
| `escape_newick_label(label)` | Vector of arbitrary node/tip labels. | Returns Newick-safe unquoted or single-quoted labels. |
| `prepare_physicell_newick(nodes, terminal_physicell_ids)` | Event-resolved node table and optional sampled/live IDs. | Validates parent-before-child topology, marks retained ancestors in one reverse pass, and returns integer child arrays plus preformatted Newick tokens. |
| `emit_physicell_newick(prepared, emit_token, show_progress, progress_updates, progress_label)` | Prepared Newick topology, token consumer, and progress controls. | Traverses the retained forest iteratively with an explicit event stack and emits each Newick token once. |
| `physicell_lineage_to_newick(nodes, terminal_physicell_ids)` | Event-resolved node table and optional sampled/live IDs. | Compatibility API returning one synthetic-root Newick string; uses the iterative renderer and preserves unary nodes and elapsed branch time. |
| `write_physicell_lineage_newick(nodes, path, terminal_physicell_ids, show_progress, progress_updates, buffer_bytes, progress_label)` | Event-resolved node table, destination, optional sampled/live IDs, progress controls, and output-buffer size. | Streams a synthetic-root Newick tree through a bounded buffer without recursion or repeated subtree-string copying. |

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
| `physicell_rate_set(params, cell_type, editing_state, barcode_sequence, be_classes, nuc_classes, be_to)` | Parsed JSON plus one cell type/editing state and target definitions. | Per-position probabilities plus cached event hazards, total hazards, cumulative event-selection probabilities, active positions, substitution matrix, and invariants. |
| `prepare_physicell_recording_model(params, cell_type, num_integrations, founder_label_sites, params_dir, seed)` | Parsed remote_mito JSON and adapter overrides. | Complete induced/uninduced barcode model, recorder identity, reference, targets, windows, induction time, integration count, profile-storage policy, and optional founder-code width. |
| `initialize_physicell_barcode_profile(model)` | Prepared barcode model. | Empty dense matrix or sparse per-integration profile, according to the adapter storage policy. |
| `physicell_barcode_value(profile, integration, position)` | Dense/sparse profile and one integration-position pair. | Encoded allele value, returning WT `0` for absent sparse entries. |
| `set_physicell_barcode_value(profile, integration, position, value)` | Dense/sparse profile, one integration-position pair, and encoded allele. | Updated profile without exposing its storage representation to callers. |
| `physicell_barcode_window_edited(profile, integration, window_positions)` | Dense/sparse profile, integration, and editing-window positions. | Whether any position in that editing window has already changed. |
| `initialize_physicell_founder_barcode(model, founder_index, num_founders)` | Prepared barcode model and one founder's rank/count. | Initial profile and founder-label events using stable allele-coded sites; validates coding capacity. |
| `physicell_probability_hazard(probability, cell_cycle_length)` | Per-division probabilities and cell-cycle duration. | Named continuous-time hazards. |
| `non_mendelian_selection_coefficient(params, marker, default)` | Parsed parameters, optional marker-specific key, and neutral fallback. | Validated coefficient in `[0,1]`, using the marker override first and the global `non_mendelian_selection.coefficient` second. |
| `palincode_probability_vector(value, number, name)` | Scalar/per-cBit probability, cBit count, and error label. | Validated probability vector with one value per PALINCODE cBit. |
| `palincode_positive_integer(value, name, allow_zero)` | Candidate count, diagnostic name, and zero policy. | Validated PALINCODE count. |
| `palincode_static_ids(number, identifier_length)` | Integration count and identifier length. | Unique random A/C/G/T static integration identifiers. |
| `prepare_palincode_recording_model(params, cell_type, num_integrations, founder_label_sites, seed)` | Native parameters plus lineage-overlay controls. | PALINCODE model with target-specific continuous hazards, locked left/right/both outcomes, static IDs, and sparse-output metadata. |
| `prepare_prime_editing_recording_model(params, cell_type, num_integrations, founder_label_sites, params_dir, seed)` | Native parameters, shared backend configuration, lineage-overlay controls, relative-path root, and seed. | Known-pegRNA model with target-specific effective hazards, locked edited states, exact templates, and integration static IDs. |
| `format_physicell_progress_duration(seconds)` | Non-negative elapsed/remaining seconds. | Compact seconds, minutes, or hours label for progress output. |
| `physicell_log_stage(message_text, enabled)` | Stage-transition text and enabled flag. | Emits one timestamped log message, or remains quiet when disabled. |
| `new_physicell_progress_reporter(total, label, enabled, updates, unit)` | Total work count, phase label, enabled flag, approximate update count, and unit label. | Closure that emits throttled count, percentage, elapsed time, throughput, and ETA messages. |

### Branch simulation and output

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `empty_physicell_barcode_events()` | None. | Empty typed core barcode-event table used by zero-event branches. |
| `select_physicell_position_events(positions, rate_set)` | Mutated positions and a prepared rate set. | Resolves cached single-event outcomes directly and samples among multiple competing event hazards. |
| `decode_physicell_barcode_events(selected_events)` | Selected substitution/insertion/deletion names. | Returns vectorized event labels, alternate bases, and encoded alleles. |
| `mutate_prime_editing_segment(profile, duration, rate_set, model, segment_start)` | Inherited prime-editing profile, branch interval, effective hazards, model, and absolute start. | Samples exact-time first edits, locks resolved targets, and reports exact programmed sequences. |
| `mutate_palincode_segment(profile, duration, rate_set, model, segment_start)` | Inherited cBit profile, branch interval, PALINCODE rates, model, and absolute start. | Samples exact-time first events and irreversibly assigns left, right, or simultaneous-both states. |
| `mutate_physicell_barcode_segment(profile, duration, rate_set, model, segment_start)` | One inherited profile, elapsed time, state-specific rates, model, and segment start. | Applies irreversible competing edits with vectorized coordinate draws and batch profile updates; uses sequential traversal only when close-after-edit windows require it. |
| `mutate_physicell_barcode_branch(profile, start_time, end_time, model, editing_state)` | Profile, branch interval, prepared model, and `auto`/forced editing state. | Splits at global induction when needed and returns the final profile plus events with interval provenance. |
| `simulate_recording_on_physicell_lineage(nodes, model, editing_state, terminal_physicell_ids, seed, show_progress, progress_updates)` | Event tree, prepared model, state policy, optional live IDs, seed, and progress controls. | Propagates inherited profiles in topological order and accumulates mutation events in growing typed columns before one final data-frame construction. |
| `grow_event_buffer(required_capacity)` | **Internal to `simulate_recording_on_physicell_lineage`.** Required mutation-event capacity. | Geometrically expands every typed event column while preserving accumulated values. |
| `append_barcode_events(events, node_id, physicell_id_value, parent_node_id, branch_start, branch_end)` | **Internal to `simulate_recording_on_physicell_lineage`.** One event batch and its branch metadata. | Copies a batch into the preallocated event columns without per-node row-binding. |
| `physicell_terminal_descendant_intervals(nodes, terminal_nodes)` | Event-resolved lineage and sampled terminal rows. | Computes a depth-first sampled-cell order and one contiguous terminal-descendant interval per node in linear time. |
| `physicell_event_descendant_matrix(nodes, terminal_nodes, mutation_events, event_ids, show_progress, progress_updates)` | Lineage, sampled terminals, event rows with `node_id`, optional column IDs, and progress controls. | Returns a sparse literal cell-by-unique-event matrix, descendant counts, and row-order mapping. |
| `write_physicell_event_descendant_outputs(nodes, terminal_nodes, event_tables, output_dir, show_progress, progress_updates, compress_csv)` | Lineage, sampled terminals, named modality event tables, destination, and output controls. | Writes the combined sparse event-descendant matrix and a compressed one-row-per-column event manifest. |
| `physicell_profile_to_sequence(profile_row, reference)` | One encoded integration profile and barcode reference. | Reconstructs its variable-length nucleotide sequence. |
| `physicell_sparse_recording_matrix(terminal_profiles, model, show_progress, progress_updates)` | Named terminal profiles, prepared model, and progress controls. | Counts nonzeros, preallocates exact sparse-triplet capacity, fills it in a second pass, and returns the cell-by-recording-position matrix. |
| `physicell_palincode_character_matrix(state_matrix, model)` | Sparse/dense PALINCODE state matrix and prepared model. | Sparse one-hot left/right/both character matrix, with wild-type cBits represented by three zeros. |
| `physicell_baseline_target_layout(model)` | Prepared barcode, prime-editing, or PALINCODE model. | Long target table; specialized recorders include static IDs, target assignments, calibrated rates, and outcome metadata. |
| `write_physicell_lineage_outputs(nodes, terminal_nodes, output_dir, show_progress, progress_updates, compress_csv)` | Event-resolved node table, sampled terminal rows, destination, progress controls, and compression flag. | Writes compressed node/tip CSVs by default and streams full/sampled Newick trees with substage and traversal progress logs. |
| `write_physicell_recording_outputs(simulation, model, output_dir, show_progress, write_lineage, progress_updates, compress_csv)` | Completed imported-lineage simulation, destination, logging controls, shared-lineage-output flag, and compression flag. | Writes exact-time barcode events, target layout, dense or compact sparse matrices, optional lineage/FASTA files, and manifest with substage logs. |

## `physicell_mito.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `physicell_mito_rate_set(params, cell_type, editing_state, mitochondrial_reference)` | Parsed JSON, selected cell type/rate state, and mt reference. | Per-position substitution/indel hazards, invariant positions, cumulative sampling weights, and cell-cycle scale. |
| `prepare_physicell_mito_model(params, cell_type, genomes_per_cell, seed)` | Parsed JSON and mitochondrial replay overrides. | Generated mt reference plus induced/uninduced rates, induction time, fixed bottleneck size, and mitochondrial variant selection coefficient. |
| `initialize_physicell_mito_profile(model, params)` | Prepared model and parsed JSON. | Founder genome profiles and baseline heteroplasmy event rows drawn from the configured site fraction, Beta penetrance, and transition probability. |
| `inherit_physicell_mito_profile(parent_profile, genomes_per_cell, variant_selection_coefficient)` | Parent mt genome list, daughter bottleneck size, and coefficient in `[0,1]`. | Daughter genomes resampled with replacement using relative weight `(1-s)^variant_count`; neutral when `s=0`. |
| `sample_physicell_mito_position(rate_set, excluded_positions)` | Weighted rate set and already mutated positions. | One mutable position sampled in proportion to its total event hazard. |
| `mutate_physicell_mito_segment(profile, duration, rate_set, model, segment_start)` | Sparse genome list, elapsed time, state-specific rates, model, and segment start. | Applies low-rate irreversible substitutions/indels and returns the updated profile plus exact-time event rows. |
| `mutate_physicell_mito_branch(profile, start_time, end_time, model, editing_state)` | Profile, branch interval, model, and `auto`/forced state. | Splits a branch at induction when needed and returns its final profile and provenance-annotated events. |
| `simulate_mito_on_physicell_lineage(nodes, model, params, editing_state, terminal_physicell_ids, seed, show_progress, progress_updates)` | Event tree, prepared model, JSON, state policy, optional current-cell IDs, seed, and progress controls. | Propagates bottlenecked mt profiles in topological order, optionally reports node progress/ETA, and returns all profiles, sampled terminals, and events. |
| `physicell_mito_variant_fractions(simulation)` | Completed mitochondrial replay. | Long-form per-terminal-cell position/allele counts and heteroplasmy fractions. |
| `physicell_mito_haplotype_sequence(genome, reference)` | One sparse encoded mt genome and reference. | Reconstructs its nucleotide sequence, including substitutions, deletions, and one-base insertions. |
| `write_physicell_mito_outputs(simulation, model, output_dir, write_fasta, show_progress, compress_csv)` | Completed replay, prepared model, destination, FASTA flag, logging flag, and compression flag. | Writes compressed mt event/fraction/manifest tables by default plus profiles/reference, optional sparse matrix, and sampled haplotypes. |

## `gillespie_lineage.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `gillespie_scalar(value, name, mode)` | Candidate scalar, diagnostic name, and numeric/integer mode. | Validated typed scalar or an explicit input error. |
| `gillespie_probability_hazard(probability, interval)` | Interval probability and interval length. | Exact constant hazard `-log(1-p)/interval`; positive infinity maps to zero hazard. |
| `gillespie_induction_selection(cell_indices, num_cells, frac_cells)` | Active node indices and fixed-count or binomial-fraction induction specification. | Selected active indices with the same semantics as native timestep induction. |
| `gillespie_induction_spec(params, name)` | Parsed JSON and induction-block name. | Normalized time/count/fraction specification; a time-only block induces all cells. |
| `gillespie_transition_matrices(params, cell_types)` | Parsed JSON and ordered cell types. | Validated named induced/uninduced daughter transition matrices, using identity matrices when omitted. |
| `gillespie_cell_type_rates(params, cell_types)` | Parsed JSON and cell types. | Named division, death, and total continuous-time propensities. |
| `prepare_gillespie_barcode_models(params, num_integrations, founder_label_sites, params_dir, seed)` | JSON, recorder layout, parameter directory, and seed. | Compatible prepared barcode model per cell type on one shared reference/layout. |
| `prepare_gillespie_mito_models(params, genomes_per_cell, seed)` | JSON, bottleneck size, and seed. | Compatible prepared mitochondrial model per cell type on one shared reference. |
| `gillespie_active_indices(active_nodes, active_count)` | Preallocated active-node slots and used-slot count (legacy grouped lists are also accepted). | Flat active-node index vector used at induction and final sampling. |
| `simulate_gillespie_population(params, end_time, seed, max_cells, show_progress, progress_updates)` | Native parameters plus duration, reproducibility, safety, and logging controls. | Exact Fenwick-tree-weighted birth/death lineage, division/event/checkpoint tables, live terminals, rate summary, and stop reason. |
| `gillespie_sc_cell_states(simulation)` | Completed Gillespie population result. | PhysiCell-compatible live-cell covariate table with cell type, generation, and developmental pseudotime. |
| `parse_matrix(value, label)` | **Internal to `gillespie_transition_matrices`.** Raw matrix and diagnostic label. | Validated named stochastic matrix or identity fallback. |
| `grow_nodes(required)` | **Internal to `simulate_gillespie_population`.** Required node capacity. | Geometrically grows and explicitly initializes every preallocated node column. |
| `update_propensity_tree(position, difference)` | **Internal to `simulate_gillespie_population`.** Active slot and propensity delta. | Updates the Fenwick tree in logarithmic time. |
| `select_propensity_slot(target)` | **Internal to `simulate_gillespie_population`.** Uniform draw on total propensity. | Finds the corresponding weighted active-cell slot in logarithmic time. |
| `add_active(node_index)` / `remove_active(node_index)` | **Internal to `simulate_gillespie_population`.** One node index. | Adds or swap-removes an active cell while keeping slots, reverse positions, total propensity, and Fenwick weights consistent. |
| `add_node(...)` | **Internal to `simulate_gillespie_population`.** Cell, parent, time, type, origin, generation, induction, and division fields. | Appends a parent-before-child node and activates it. |
| `append_event(...)` | **Internal to `simulate_gillespie_population`.** Population-event fields. | Appends one typed division/death/induction provenance row. |
| `select_active(specification)` / `apply_inductions(time)` | **Internal to `simulate_gillespie_population`.** Induction specification or boundary time. | Selects current cells and creates exact-time continuation nodes with heritable state changes. |
| `record_checkpoint(time)` / `report_time_progress()` | **Internal to `simulate_gillespie_population`.** Boundary/current time. | Records type counts or emits throttled time/population/event progress. |

## `gillespie_pipeline.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `gillespie_pipeline_option(configuration, name, default)` | Gillespie JSON block, key, and fallback. | Configured value or fallback. |
| `gillespie_modalities(value)` | Scalar/vector modality selection; accepts comma-separated strings, the `mito` alias, and the `both`/`all` expansions. | Validated unique modality vector drawn from `barcode`, `mitochondrial`, and `ecdna`; `both` expands to the first two and `all` to all three. |
| `gillespie_timing_recorder()` | None. | Phase-finishing and timing-summary closures. |
| `run_gillespie_lineage_pipeline(params, params_path, output_dir, overrides)` | Parsed parameters, source path, destination, and CLI overrides. | Runs population and continuous recording phases for each selected modality, including the ecDNA stage; writes lineage, feature, event-descendant, scDesign3 covariate, manifest, and timing outputs. The combined feature matrix is written only when both barcode and mitochondrial outputs exist, so an ecDNA-only run produces no combined matrix. |
| `finish(name)` / `summary()` | **Internal to `gillespie_timing_recorder`.** Phase name or no arguments. | Records a phase boundary or returns phase/total elapsed-time rows. |

## `simulate_gillespie_lineage.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `gillespie_cli_usage()` | None. | Standalone Gillespie command-line help text. |
| `parse_gillespie_cli(arguments)` | Trailing command-line arguments. | Typed standalone options with validation and help handling. |

The top-level script reads a native parameter JSON and calls
`run_gillespie_lineage_pipeline()`. `sim5_code.R` makes the same call when the
JSON contains `"simulation_engine": "gillespie"`; otherwise its historical
timestep path is unchanged.

## `engine_comparison.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `simulate_timestep_population_outcome(params, end_time, time_step, seed, max_cells)` | Native population parameters, duration, grid size, seed, and safety cap. | Vectorized discrete competing-hazard population outcome with divisions, deaths, extinction, stop reason, and final cell-type counts. |
| `assign_differentiation()` | **Internal to `simulate_timestep_population_outcome`.** No arguments. | Applies differentiation at the first grid boundary at or after its configured time. |
| `summarize_gillespie_population_outcome(simulation)` | One exact population result. | Comparison-compatible continuous-time outcome and cell-type counts. |
| `compare_population_engines(params, replicates, end_time, time_step, seed, max_cells, show_progress, progress_updates)` | Parameters and paired-benchmark controls. | Per-replicate continuous/time-step outcomes and long-form cell-type tables. |
| `summarize_engine_comparison(comparison)` | Completed paired comparison. | Means, standard deviations, absolute differences, and relative differences for population/event/runtime metrics. |
| `summarize_engine_cell_types(comparison)` | Completed paired comparison. | Per-engine cell-type count and fraction summaries. |
| `write_engine_comparison_plot(comparison, path)` | Completed comparison and PNG destination. | Writes final-population boxplots and mean division/death bars. |

## `compare_simulation_engines.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `comparison_cli_usage()` | None. | Returns comparison command-line help. |
| `parse_comparison_cli(arguments)` | Trailing command-line arguments. | Validates and returns typed replicate, duration, timestep, cap, seed, progress, compression, and path options. |

The top-level comparison script reads one native JSON, runs a population-only
paired benchmark, and writes compressed replicate/summary tables plus a PNG.

## `lineage_benchmark.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `lineage_benchmark_integer_vector(value, name, minimum)` | Numeric/list/comma-delimited values and validation labels. | Sorted unique integer vector at or above the requested minimum. |
| `lineage_benchmark_log(..., enabled)` | Message components and logging switch. | Emits one timestamped benchmark progress line when enabled. |
| `lineage_benchmark_editing_state(be_probability, mt_probability)` | BASELINE target and mitochondrial substitution probabilities. | Native editing-state block with the supplied mutation probabilities. |
| `lineage_benchmark_cell_type(cell_cycle_length, death_probability, be_probability)` | Population and recording rates. | One complete benchmark cell-type parameter block. |
| `lineage_benchmark_population_params(shape, maximum_cells, seed)` | Named topology/population shape, cell cap, and seed. | Native population parameters for balanced, comb, neutral, turnover, or hierarchical growth. |
| `lineage_benchmark_synthetic_divisions(shape, terminal_cells, duration)` | Balanced/comb shape, exact terminal count, and duration. | PhysiCell-format persistent-ID division table with the requested topology extreme. |
| `lineage_benchmark_annotate_nodes(nodes, cell_type)` | Synthetic event-resolved nodes and type label. | Nodes augmented with generation, induced state, fate, and live/death fields. |
| `simulate_lineage_benchmark_population(shape, terminal_cells, seed, synthetic_duration, show_progress)` | Shape, exact live-cell target, seed, duration, and logging. | Exact-size synthetic or Gillespie population bundle plus parameters. |
| `lineage_benchmark_sample_sets(population, tree_sizes, seed)` | Maximum population, requested sample sizes, and sampling seed. | Nested deterministic terminal-ID sets. |
| `lineage_benchmark_empty_target_spec()` | None. | Valid zero-target native target block. |
| `lineage_benchmark_prime_pool()` | None. | Six fixed benchmark pegRNAs with programmed edits and efficiencies 0.90–0.15. |
| `lineage_benchmark_targets_per_integration(system)` | BASELINE, prime, PALINCODE, or mitochondrial name. | Logical target count of 50, 6, or 2; mitochondrial returns `NA` because depth is not an integration count. |
| `lineage_benchmark_mito_params(population_params, substitution_probability)` | Population parameters and per-destination mitochondrial substitution probability. | Complete benchmark mitochondrial parameters with a 16,569-base reference and sparse profiles. |
| `lineage_benchmark_recording_params(population_params, system, maximum_integrations)` | Population parameters, recorder name, and integration count. | Complete sparse event-resolved recorder parameters using the benchmark chemistry preset. |
| `lineage_benchmark_column_integrations(column_names)` | Native matrix column names. | Parsed integration number for every column. |
| `lineage_benchmark_subset_matrix(matrix_value, sample_ids, integrations)` | Full sparse matrix, ordered sampled cells, and integration limit. | Paired cell/integration submatrix. |
| `lineage_benchmark_subset_rows(matrix_value, sample_ids)` | Full sparse matrix and ordered sampled cells. | Cell-subset matrix retaining every feature column. |
| `lineage_benchmark_logical_target_matrix(state_matrix, model, integrations)` | Native states, prepared model, and integration count. | Unified categorical matrix with one column per physical target; BASELINE windows become bitmask states. |
| `lineage_benchmark_logical_target_manifest(target_layout, model)` | Native target layout and model. | Layout augmented with ordered logical positions and categorical encoding. |
| `lineage_benchmark_write_json(value, path)` | Serializable value and destination. | Pretty reproducibility JSON with scalar unboxing. |
| `write_lineage_benchmark_population(population_bundle, output_dir, show_progress)` | Population/parameter bundle and destination. | Saves the reusable population, parameters, lineage tables, division events, and full/live Newick trees. |
| `write_lineage_benchmark_samples(population, sample_sets, output_dir, show_progress)` | Population, nested terminal sets, and destination. | Writes one sampled-cell table and ground-truth Newick per requested tree size. |
| `lineage_benchmark_mito_sampling_orders(terminal_profiles, maximum_depth, seed)` | Named terminal profiles, maximum observation depth, and seed. | Reproducible without-replacement genome order per cell for nested depth panels. |
| `lineage_benchmark_mito_variant_manifest(terminal_profiles, sampling_orders, maximum_depth, reference)` | Profiles, nested orders, maximum depth, and mitochondrial reference. | Sorted union of observed mitochondrial variants with reference/alternate annotations. |
| `lineage_benchmark_mito_observation_matrices(terminal_profiles, sampling_orders, depths, variant_manifest)` | Profiles, nested orders, depths, and shared variant definition. | Per-depth sparse heteroplasmy and binary-presence matrices on a common feature universe. |
| `simulate_lineage_benchmark_mitochondrial(population_bundle, observation_depths, genomes_per_cell, seed, output_dir, show_progress)` | Population, mitochondrial sampling panel, modeled copy number, seed, destination, and logging. | One event-resolved mitochondrial overlay plus nested observation matrices and exact-event output. |
| `simulate_lineage_benchmark_recorder(population_bundle, system, maximum_integrations, seed, output_dir, show_progress, observation_depths, mitochondrial_genomes_per_cell)` | Reusable population, system, maximum recording level, seed, destination, logging, and mitochondrial controls. | Simulates and writes one maximum integrated or mitochondrial recorder overlay for subsetting. |
| `load_lineage_benchmark_recorder(output_dir)` | Completed maximum-recorder directory. | Reloaded model, native matrices, and target layout for resumption. |
| `lineage_benchmark_relative_path(path, root)` | Candidate path and benchmark root. | Portable root-relative path when the candidate is inside the benchmark. |
| `write_lineage_benchmark_condition(...)` | Full recording, paired cell/integration condition, truth/sample paths, and output controls. | Writes condition-specific truth, native/logical sparse matrices, target/sample tables, manifest, and completion marker. |
| `load_lineage_benchmark_condition(output_dir)` | Condition directory. | Completed condition manifest, or `NULL` when unfinished. |
| `run_lineage_benchmark(...)` | Output, shapes, sizes, integrations, systems, seeds, duration, dense/resume/logging controls, mitochondrial observation depths, and mitochondrial copy number. | Resumable paired population/recorder sweep plus root reconstruction-job manifest and summary. |

## `run_lineage_benchmark.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `lineage_benchmark_cli_usage()` | None. | Simulation-only benchmark command-line help. |
| `lineage_benchmark_cli_boolean(value, flag)` | Text value and flag name. | Validated logical CLI value. |
| `lineage_benchmark_cli_sequence(value, name)` | Comma-delimited integers or `A:B` range. | Validated integer seed vector. |
| `parse_lineage_benchmark_cli(arguments)` | Trailing command-line arguments. | Typed output/grid/resume/logging options. |

The top-level script loads the shared prime-editing, event-resolved lineage,
Gillespie, and benchmark modules; it writes simulation truth and recording data
only and prints `BENCHMARK_OUTPUT_DIR=<path>` when complete.

## `physicell_visium.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `visium_numeric_vector(value, name, minimum_length, allow_infinite)` | Numeric/list/comma-delimited value and validation controls. | Validated numeric vector, optionally retaining infinite distance-bin bounds. |
| `visium_plane_basis(normal)` | Three-vector section normal. | Orthonormal section-normal, `u`, and `v` basis vectors. |
| `make_visium_6_5mm_array(spot_pitch, rotation_degrees, translation, array_rows, spots_per_row, slice_id)` | Conventional-array geometry and alignment. | Staggered Space Ranger-style row/column table; defaults produce exactly 4,992 spots in 78 rows. |
| `read_physicell_visium_inputs(run_dir)` | Completed PhysiCell pipeline directory. | Validated live-cell coordinates joined to terminal sample/node IDs plus event-resolved nodes. |
| `slice_physicell_cells(cells, basis, center, offset, thickness)` | Live cells and section geometry. | Cells whose centers intersect the finite-thickness plane, with projected coordinates. |
| `assign_cells_to_visium_spots(slice_cells, spots, spot_diameter)` | Section cells and circular spot grid. | Every section cell annotated with its nearest spot and capture status; spot table receives tissue flags and cell counts. |
| `build_visium_lca_index(nodes)` | Topologically ordered event-resolved lineage nodes. | Parent/depth/root arrays plus binary-lifting ancestors for repeated logarithmic-time MRCA queries. |
| `query_visium_lineage_pairs(index, node_a, node_b)` | LCA index and equal-length node-index vectors. | MRCA identity/time, same-founder flag, MRCA age, patristic distance, and within-founder division-edge distance for every pair. |
| `sample_visium_cell_pairs(number_cells, max_exact_pairs, max_candidate_pairs, seed)` | Cell count, pair limits, and seed. | Every unordered pair below the exact limit or uniform pairs with replacement above it. |
| `visium_distance_labels(breaks)` | Increasing spatial-distance edges. | Stable left-closed interval labels; its internal `format_bound(value)` prints finite and infinite endpoints. |
| `prepare_visium_spatial_pairs(cells, distance_breaks, max_pairs_per_bin, max_exact_pairs, max_candidate_pairs, seed)` | Captured cells and scalable pair controls. | Stratified pair table with 2-D, 3-D, spot-center distances and candidate/evaluated counts. |
| `summarize_visium_lineage_pairs(pair_table, distance_labels, candidate_counts, recent_threshold_hours, null_summaries)` | Lineage-annotated cell pairs and bin/null definitions. | Cell-distance correlogram with MRCA, patristic, division, founder, same-spot, recent-MRCA, and null-enrichment summaries. |
| `summarize_visium_spot_pair_lineage(pair_table, spot_pitch)` | Lineage-annotated cells pairs and Visium pitch. | Equal-spot-pair within/between-spot tables and distance-bin summaries, avoiding cell-rich spot domination. |
| `permute_visium_lineage_null(pair_table, captured_node_indices, lineage_index, distance_labels, recent_threshold_hours, permutations, time_units_per_hour, seed)` | Fixed spatial pairs, lineage tips, null controls, and time scale. | Mean MRCA/recent-ancestor expectations after independently permuting lineage locations. |
| `analyze_visium_slice_pairs(captured_cells, lineage_index, distance_breaks, spot_pitch, ...)` | One slice’s captured cells plus lineage/pair/null settings. | Annotated cell pairs, cell/spot summaries, and overall slice metrics. |
| `write_visium_slice_plot(slice_cells, spots, path, spot_diameter)` | Section cells, spots, destination, and footprint size. | PNG overlay distinguishing captured cells, gap cells, and occupied spot circles. |
| `write_visium_correlogram_plot(summary, path)` | Cell-distance summary and destination. | PNG observed-versus-permuted mean-MRCA curve. |
| `run_physicell_visium_analysis(...)` | Completed run, slice/array geometry, pair/null limits, time scale, seed, and output controls. | Analyzes all sections, writes compressed cell/spot/pair/summary tables and plots, and returns combined results. |
| `aggregate_physicell_visium_results(batch_dir, output_dir, compress_csv)` | Replicate batch and destination. | Pools slices within tumors, then reports across-tumor means/SD/SE for cell- and spot-distance summaries plus an aggregate plot. |

## `analyze_physicell_visium.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `physicell_visium_cli_usage()` | None. | Returns single-run Visium analysis help and defaults. |
| `parse_physicell_visium_cli(arguments)` | Trailing command-line arguments. | Typed and validated run, section, spot, pair, null, seed, and output options. |

The top-level CLI analyzes one completed PhysiCell run and prints its slice
occupancy/correlation summary plus `OUTPUT_DIR=<path>`.

## `aggregate_physicell_visium.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `visium_aggregate_usage()` | None. | Returns replicated-batch aggregation help. |

The top-level CLI parses batch/output/compression options, calls
`aggregate_physicell_visium_results()`, and prints `OUTPUT_DIR=<path>`.

## `notebooks/physicell_visium_lineage_distance.Rmd`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `resolve_notebook_csv(path)` | Compressed or uncompressed CSV base path. | Existing `.csv.gz`/`.csv` path, preferring compressed output. |
| `read_notebook_csv(path)` | CSV base path. | Data frame read transparently from gzip or plain CSV. |
| `discover_visium_analyses(input_path)` | One analysis directory or replicate-batch root. | Sorted `visium_spatial` directories containing distance summaries. |
| `notebook_replicate_id(analysis_dir)` | One `visium_spatial` directory. | Parent tumor/replicate directory name used in facets and summaries. |
| `sample_notebook_pairs(path, maximum, seed)` | Pair-table path, plotting cap, and seed. | Full pair table or a reproducible display-only row sample. |

The parameterized R notebook loads cell/spot summaries and saved cell pairs;
it plots physical distance against MRCA age and lineage-edge distance,
observed/permuted correlograms, recent-relative enrichment, spot mixing,
section projection, and tumor-replicate means/standard errors. It does not
rerun PhysiCell or lineage simulation.

## `render_physicell_visium_notebook.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `physicell_visium_notebook_usage()` | None. | Returns notebook-renderer command-line usage and defaults. |
| `parse_physicell_visium_notebook_cli(arguments)` | Trailing CLI arguments. | Validated input/output, plotting-cap, threshold, and seed options. |

The top-level CLI stages the Rmd to avoid absolute-output resource-path issues,
renders a self-contained HTML notebook, copies it to the requested location,
and prints `NOTEBOOK_OUTPUT=<path>`.

## `ecdna_lineage.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `ecdna_config_value(configuration, name, default)` | ecDNA configuration, key, and fallback. | Configured value or fallback. |
| `ecdna_probability(value, name)` | Candidate probability and diagnostic name. | Validated scalar in `[0,1]`. |
| `ecdna_positive_integer(value, name, allow_zero)` | Candidate count, diagnostic name, and zero policy. | Validated integer count. |
| `prepare_ecdna_model(params, seed)` | Native parameters and seed. | Static species IDs, labeled species selected by nearest requested fraction or exact count, initial copy counts, recorder layout/hazards, selection/replication/segregation settings, and copy cap. |
| `initialize_ecdna_profile(model)` | Prepared model. | Founder species/label/haplotype copy groups. |
| `collapse_ecdna_profile(profile)` | Grouped ecDNA profile. | Sorted profile with identical species/label/haplotype rows merged and zero counts removed. |
| `ecdna_truncated_event_times(number, hazard, duration, start_time)` | Number of observed edits and branch-process parameters. | Exact event times conditional on occurring within the segment. |
| `mutate_ecdna_profile_segment(profile, duration, segment_start, hazards, model)` | Inherited profile and active recorder interval. | Irreversibly edited grouped haplotypes plus copy-count-aware aggregate event rows. |
| `mutate_ecdna_profile_branch(profile, start_time, end_time, cell_type, model)` | Profile, cellular branch, type, and model. | Applies only the recorder-active portion using cell-type-specific hazards. |
| `cap_ecdna_profile(profile, maximum)` | Profile and cell-level copy cap. | Multinomially downsampled grouped profile when the cap is exceeded. |
| `partition_ecdna_profile(profile, model)` | Parent profile and selection/replication/segregation settings. | Two jointly sampled complementary daughter profiles after selection-weighted copy replication. |
| `simulate_ecdna_on_lineage(nodes, model, terminal_physicell_ids, seed, show_progress, progress_updates)` | Event-resolved cellular tree, prepared model, sampled cells, and run controls. | Propagates/edit ecDNA profiles; partitions only at binary divisions; returns terminal profiles and edit provenance. |
| `ecdna_terminal_tables(simulation, model)` | Completed ecDNA simulation. | Per-cell burden summary and long-form static-ID/haplotype copy table. |
| `ecdna_feature_matrices(simulation, model)` | Completed simulation and model. | Sparse static-ID copy-number, target edit-fraction, and binary character matrices based on actual terminal inheritance. |
| `write_ecdna_outputs(simulation, model, output_dir, show_progress, compress_csv)` | Completed simulation, model, destination, and output controls. | Writes compressed manifests/events/summaries/haplotypes, terminal profiles, and sparse matrices. |

## `simulate_physicell_lineage.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `physicell_cli_usage()` | None. | Returns the command-line usage text. |
| `parse_physicell_cli_args(arguments)` | Trailing command-line arguments. | Validates supported short/long options, including modality/mt controls, and returns a typed options list; `--help` prints usage and exits. |
| `new_physicell_timing_recorder(total_start, clock)` | Optional total-start timestamp and monotonic clock function. | Returns closures that finish sequential phases and produce phase/total wall-clock durations. |
| `finish_phase(label)` | **Internal to `new_physicell_timing_recorder`.** Unique phase label. | Records elapsed wall time since the preceding phase boundary. |
| `summary()` | **Internal to `new_physicell_timing_recorder`.** None. | Returns phase names, seconds, human-readable durations, percentages, and total R runtime. |
| `format_physicell_timing_summary(timing_summary)` | Timing summary data frame. | Returns aligned terminal/log lines for every phase and total runtime. |

The script’s top-level code loads a remote_mito JSON, imports the PhysiCell
division/current-cell files, writes ground-truth lineage output, optionally
prepares and simulates the requested barcode, mitochondrial, and/or ecDNA
models, writes `r_timing_summary.csv.gz` by default, prints a phase/total timing
summary, and prints `OUTPUT_DIR=<path>`. `--modalities lineage` stops after the
ground-truth tables and Newick trees.

## `scdesign3_helpers.R`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `scdesign3_log_stage(message_text)` | Non-empty stage-transition text. | Emits one timestamped scDesign3 phase message. |
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
| `build_physicell_sc_covariates(recording_dir, lineage_table_path, cell_type)` | Recording output, final PhysiCell table, and optional fallback type. | Terminal metadata joined across lineage, per-cell developmental state/pseudotime, founder, niche/spatial fields, barcode burden, mitochondrial heteroplasmy, and optional ecDNA copy/species/recorder burdens. |
| `apply_physicell_celltype_map(metadata, celltype_map)` | Target metadata and named simulated-to-reference mapping. | Metadata with mapped model type and retained original `cell_type_sim`. |
| `scale_physicell_covariate_to_reference(values, reference_values)` | Simulated and reference numeric covariates. | Simulated values min/max mapped into the reference range, with degenerate-range handling. |
| `align_physicell_sc_covariates(metadata, reference_sce, use_pseudotime, use_spatial)` | Joined target metadata and standardized reference. | Adds reference-range pseudotime and/or `spatial1`/`spatial2` predictors. |
| `write_physicell_scdesign3_outputs(counts, metadata, output_dir, fit, reference_path, seed, compress_csv)` | Synthetic counts, linked metadata, fit/provenance, destination, and compression flag. | Writes a terminal-cell SCE, count RDS, and compressed metadata/summary/manifest CSVs by default. |

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
| `get_profiles_ints_and_umis(cell_pop, bc_or_mt, int_rec_prob, num_ints, umis)` | Population, modality, recovery probability/count, and optional UMIs. | Per-cell list of sampled mutation matrices, recovered integration indices, and optional recovered UMIs; zero recovery produces a well-formed zero-row profile. |
| `ins_to_charvec(ins, pos_num, ref_seq, fixed_length)` | Encoded insertion, position/reference, and fixed-length mode. | Returns the number and character vector of decoded inserted/reference bases. |
| `get_one_cell_sequence(cell_int_mat, ref_seq, these_bc_int_umis)` | One cell’s integrations/genomes, reference, and optional UMIs. | List of reconstructed nucleotide strings, one per recovered row. |
| `write_all_cell_sequences(...)` | Named cell matrices, reference, output path/type, and optional UMIs. | Reconstructs sequences in parallel and writes a FASTA file with the requested type suffix, omitting cells with no recovered molecule and emitting an empty file when none remain. |

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

## `convert_scoremats_to_csvs.r`

This file defines no functions. It is an `optparse` CLI, invoked by
`bash_wrapper_all_combos.sh` after `process_results_from_bash.r` and only when
`output/score_mats/<run_id>/matrices/` exists:

```
Rscript convert_scoremats_to_csvs.r --urid <run_id> \
  --score_mat_path <base_path>/output/score_mats/<run_id>/matrices/
```

| Option | Meaning |
| --- | --- |
| `-P`, `--score_mat_path` | Directory of RDS score matrices to export. Every path the script uses is derived from this. |
| `-U`, `--urid` | Unique run id. Accepted for call compatibility with the other `*_from_bash.r` scripts; the script never reads it. |

The top-level code lists `*.rds` directly inside `--score_mat_path`
(non-recursively, so the nested `af/` allelic-fraction matrices are skipped),
creates a `csvs/` subdirectory of that path if needed, then for each file
densifies the sparse matrix with `as.matrix()` and writes
`<same-basename>.csv` via `write.csv`, preserving the cell row names and
`<mt|bc>_<mut_name>` column names set by `mut_to_scoremat.r`. Existing CSVs of
the same name are overwritten. Note that with the wrapper's invocation the
output directory is `output/score_mats/<run_id>/matrices/csvs/`, one level
deeper than the `output/score_mats/<run_id>/csvs/` quoted in `README.md`.

## `compare_celltype_probs.r`

This file defines no functions. It is a standalone diagnostic, outside the
`bash_wrapper_all_combos.sh` pipeline, run as `Rscript compare_celltype_probs.r`
from the repo root or stepped through interactively. It takes no command-line
arguments and writes no files.

It reads `celltype_prob_concordance/<urid>/*.rds` for a single hard-coded
`urid`, where the filenames encode cell type (`ct1`..`ct6`), mutation type, and
induction state (`induced_editing_params` / `uninduced_editing_params`). Those
inputs are only produced when the `saveRDS` block near the end of the
per-cell-type / per-induction parameter loop in `sim5_code.R` is re-enabled — it
is commented out in the committed source — so the script cannot run against the
repository as checked in.

The top-level code does three things:

- Loads the induced barcode-target **transition** parameter list for each cell
  type into a long `probs` / `pos` / `celltype` data frame. `mut_types` names
  the other three mutation classes but the loop is restricted to `'transition'`.
- Reads `basepos_erc_be_list<ct>_induced_editing_params.rds` for all six cell
  types and prints an `all.equal` result for every ordered pair, testing
  `sim5_code.R`'s assumption that a target's base-editing edit-rate class
  (`High` / `Medium` / `Low`) is shared across cell types even though the
  underlying numeric rates are not.
- Pivots the stacked frame to one row per base position and one column per cell
  type (`compare_pos_df`), and draws a dodged `ggplot2` barplot of transition
  probability by position, filled by cell type, to the active graphics device.
  `pos` is a character column, so both the table and the plot order positions
  lexicographically.

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

### `run_physicell_visium_replicates.sh`

| Function | Inputs | Result / side effects |
| --- | --- | --- |
| `usage()` | None. | Prints replicate, tumor, section, Visium-grid, pair-sampling, null, and resume controls. |

The top-level wrapper runs or resumes independent PhysiCell tumors using
lineage-only replay by default, analyzes every requested virtual section,
aggregates results at the independent-tumor level, and writes a replicate
manifest plus `BATCH_OUTPUT_DIR=<path>`.

### `run_physicell_organoid_ecdna_pipeline.sh`

This top-level wrapper delegates to
`run_physicell_organoid_baseline_pipeline.sh` with `--modalities all`, enabling
BASELINE-like chromosomal recording, mitochondrial tracing, and the configured
non-Mendelian ecDNA barcode overlay.

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
