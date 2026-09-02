# Covariate construction and output helpers for PhysiCell-linked scDesign3
# simulation. Source physicell_lineage.R and scdesign3_helpers.R first.

#' Read a two-column PhysiCell property manifest into a named vector
#'
#' Resolves `path` through `resolve_physicell_csv_path()` so either the plain or
#' the `.gz` form is accepted, and treats an absent file as "no properties"
#' rather than an error.
#'
#' @param path Path to a manifest CSV with `property` and `value` columns. The
#'   file is optional; a missing one yields an empty result.
#' @return A named character vector whose names are the `property` column and
#'   whose values are the `value` column, or a zero-length named character
#'   vector when the file does not exist.
read_physicell_property_manifest <- function(path){
  resolved_path <- resolve_physicell_csv_path(path, required = FALSE)
  if(is.na(resolved_path)){
    return(setNames(character(), character()))
  }
  manifest <- read_physicell_csv(
    resolved_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if(!all(c('property', 'value') %in% names(manifest))){
    stop(sprintf(
      'Manifest %s must contain property and value columns.',
      resolved_path
    ))
  }
  setNames(as.character(manifest$value), as.character(manifest$property))
}

#' Compute generation depth for every node of a lineage node table
#'
#' Roots (rows whose `parent_node_id` is missing or empty) get depth `0` and
#' every other node gets its parent's depth plus one. Resolution proceeds in
#' repeated sweeps over the still-unresolved rows, so parents may appear in any
#' row order; a sweep that resolves nothing means the table has a dangling
#' parent or a cycle and is reported as an error.
#'
#' @param nodes Data frame of lineage nodes; must contain `node_id` and
#'   `parent_node_id`. `node_id` must be unique and non-missing.
#' @return An integer vector of depths named by `node_id`.
#' @note Cost is one pass per generation, so a deep lineage costs roughly
#'   `depth * nrow(nodes)` parent lookups.
physicell_node_depths <- function(nodes){
  required_columns <- c('node_id', 'parent_node_id')
  if(!all(required_columns %in% names(nodes))){
    stop('Lineage nodes must contain node_id and parent_node_id.')
  }
  if(anyNA(nodes$node_id) || anyDuplicated(nodes$node_id)){
    stop('Lineage node IDs must be unique and non-missing.')
  }
  depths <- setNames(rep(NA_integer_, nrow(nodes)), nodes$node_id)
  unresolved <- seq_len(nrow(nodes))
  while(length(unresolved) > 0){
    progress <- FALSE
    next_unresolved <- integer()
    for(index in unresolved){
      parent <- nodes$parent_node_id[index]
      if(is.na(parent) || !nzchar(parent)){
        depths[[nodes$node_id[index]]] <- 0L
        progress <- TRUE
      } else if(parent %in% names(depths) && !is.na(depths[[parent]])){
        depths[[nodes$node_id[index]]] <- depths[[parent]] + 1L
        progress <- TRUE
      } else{
        next_unresolved <- c(next_unresolved, index)
      }
    }
    if(!progress){
      stop('Lineage nodes contain missing parents or a parent cycle.')
    }
    unresolved <- next_unresolved
  }
  depths
}

#' Count neighbors in semicolon-delimited PhysiCell neighbor-ID fields
#'
#' @param values Vector coerced to character, each element a `;`-delimited list
#'   of neighbor cell IDs. `NA` and blank/whitespace-only entries count as zero.
#' @return An integer vector of neighbor counts, named by the input values (the
#'   caller typically drops the names with `unname()`).
physicell_neighbor_count <- function(values){
  values <- as.character(values)
  vapply(values, function(value){
    if(is.na(value) || !nzchar(trimws(value))){
      return(0L)
    }
    length(strsplit(value, ';', fixed = TRUE)[[1]])
  }, integer(1))
}

#' Assemble the per-terminal-cell covariate table used to condition scDesign3
#'
#' Joins the recorded terminal cells (`terminal_cells.csv`) to the lineage node
#' table (`lineage_nodes.csv`) and to the final PhysiCell cell table, then
#' derives the lineage, spatial, and recorder covariates that the marginal model
#' is conditioned on. One row per terminal cell, keyed by `sample_id`.
#'
#' @details
#' Lineage depth comes from `physicell_node_depths()`; `lineage_pseudotime` is
#' that depth divided by the maximum depth (all zero when the lineage has a
#' single generation). `developmental_pseudotime` is taken from the PhysiCell
#' table when that column exists and otherwise defaults to `lineage_pseudotime`;
#' it must lie in `[0, 1]` and is copied to `pseudotime`, the name scDesign3
#' predicts on. `culture_day` defaults to `end_time / 1440`, i.e. PhysiCell
#' minutes converted to days.
#'
#' Recorder burdens are added only when their tables are present, and default to
#' zero otherwise. Barcode edits are read from
#' `barcode_binary_score_matrix_sparse.rds` if present (requires `Matrix`) or
#' else from `barcode_binary_score_matrix.csv`; that matrix is cells x barcode
#' sites with cells as row names, and any non-zero score counts as an edit, so
#' `barcode_edit_count` is a row sum and `barcode_edit_fraction` a row mean.
#' Mitochondrial variants are summarised per `sample_id` from
#' `mitochondrial_variant_fractions.csv`, whose fractions must lie in `(0, 1]`.
#' ecDNA copy/species/recorder columns come from `ecdna_cell_summary.csv`.
#'
#' @param recording_dir Directory holding the lineage recording output;
#'   `terminal_cells.csv` and `lineage_nodes.csv` are required, and
#'   `run_manifest.csv`, the barcode matrix, the mitochondrial table, and the
#'   ecDNA summary are optional. Each is resolved plain or `.gz`.
#' @param lineage_table_path Path to the final PhysiCell cell table; must exist
#'   and carry `ID`, `x`, `y`, `z`, and `neighbor_IDs` columns.
#' @param cell_type Fallback cell-type label used when the PhysiCell table has
#'   no per-cell `cell_type` column. When `NULL`, it is read from the
#'   `cell_type` property of `run_manifest.csv`; if neither supplies one, the
#'   call fails.
#' @return A data frame with one row per terminal cell containing `sample_id`,
#'   `physicell_id`, `node_id`, `cell_type`, lineage fields (`lineage_depth`,
#'   `lineage_pseudotime`, `developmental_pseudotime`, `pseudotime`,
#'   `birth_time`, `sampling_time`, `branch_length`), state fields
#'   (`culture_day`, `state_start_time`, `time_in_state`, `transition_count`),
#'   micro-environment fields (`oxygen`, `nutrient`), spatial fields (`x`, `y`,
#'   `z`, `neighbor_count`, `tumor_radius`, `radial_position`), `founder_id`,
#'   and the barcode, mitochondrial, and ecDNA burden columns.
build_physicell_sc_covariates <- function(recording_dir,
                                          lineage_table_path,
                                          cell_type = NULL){
  required_recording_files <- c(
    terminal_cells = file.path(recording_dir, 'terminal_cells.csv'),
    lineage_nodes = file.path(recording_dir, 'lineage_nodes.csv')
  )
  required_recording_files <- vapply(
    required_recording_files,
    resolve_physicell_csv_path,
    character(1),
    required = FALSE
  )
  missing_files <- names(required_recording_files)[
    is.na(required_recording_files)
  ]
  if(length(missing_files) > 0){
    stop(sprintf(
      'Missing PhysiCell recording tables: %s.',
      paste(missing_files, collapse = ', ')
    ))
  }
  lineage_table_path <- normalizePath(lineage_table_path, mustWork = TRUE)

  terminals <- read_physicell_csv(
    required_recording_files[['terminal_cells']],
    colClasses = c(
      sample_id = 'character',
      physicell_id = 'character',
      node_id = 'character'
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  nodes <- read_physicell_csv(
    required_recording_files[['lineage_nodes']],
    colClasses = c(
      node_id = 'character',
      physicell_id = 'character',
      parent_node_id = 'character'
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  current_cells <- utils::read.csv(
    lineage_table_path,
    colClasses = c(
      ID = 'character',
      parent_ID = 'character',
      neighbor_IDs = 'character'
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if(!all(c('sample_id', 'physicell_id', 'node_id') %in% names(terminals)) ||
     !all(c('ID', 'x', 'y', 'z', 'neighbor_IDs') %in% names(current_cells))){
    stop('PhysiCell terminal or current-cell table has an invalid schema.')
  }
  if(anyNA(terminals$sample_id) || anyDuplicated(terminals$sample_id)){
    stop('PhysiCell terminal sample IDs must be unique and non-missing.')
  }
  if(anyDuplicated(current_cells$ID)){
    stop('PhysiCell current-cell IDs must be unique.')
  }

  terminal_node_indices <- match(terminals$node_id, nodes$node_id)
  if(anyNA(terminal_node_indices)){
    stop('One or more terminal nodes are absent from lineage_nodes.csv.')
  }
  current_cell_indices <- match(terminals$physicell_id, current_cells$ID)
  if(anyNA(current_cell_indices)){
    stop('One or more terminal cells are absent from the PhysiCell lineage table.')
  }
  depth_by_node <- physicell_node_depths(nodes)
  lineage_depth <- unname(depth_by_node[terminals$node_id])
  max_depth <- max(lineage_depth)
  lineage_pseudotime <- if(max_depth == 0){
    rep(0, length(lineage_depth))
  } else{
    lineage_depth / max_depth
  }

  has_per_cell_type <- 'cell_type' %in% names(current_cells)
  if(!has_per_cell_type && is.null(cell_type)){
    recording_manifest <- read_physicell_property_manifest(
      file.path(recording_dir, 'run_manifest.csv')
    )
    cell_type <- unname(recording_manifest[['cell_type']])
  }
  if(!has_per_cell_type &&
     (is.null(cell_type) || length(cell_type) != 1 || is.na(cell_type) ||
      !nzchar(cell_type))){
    stop('A cell type was not supplied and could not be read from run_manifest.csv.')
  }

  node_rows <- nodes[terminal_node_indices, , drop = FALSE]
  cell_rows <- current_cells[current_cell_indices, , drop = FALSE]
  cell_types <- if(has_per_cell_type){
    trimws(as.character(cell_rows$cell_type))
  } else{
    rep(as.character(cell_type), nrow(terminals))
  }
  if(anyNA(cell_types) || any(!nzchar(cell_types))){
    stop('PhysiCell per-cell cell_type values must be non-missing strings.')
  }

  #' Internal: read an optional numeric column from the PhysiCell cell rows
  #'
  #' @param column Column name looked up in `cell_rows`, the PhysiCell table
  #'   subset to the terminal cells in `metadata` order.
  #' @param default Value used when the column is absent; either length one
  #'   (recycled) or already one value per terminal cell.
  #' @return A numeric vector with one finite value per terminal cell.
  optional_numeric <- function(column, default){
    if(column %in% names(cell_rows)){
      values <- suppressWarnings(as.numeric(cell_rows[[column]]))
      if(any(!is.finite(values))){
        stop(sprintf('PhysiCell %s values must be finite numeric data.', column))
      }
      values
    } else{
      if(length(default) == nrow(cell_rows)){
        as.numeric(default)
      } else if(length(default) == 1){
        rep(as.numeric(default), nrow(cell_rows))
      } else{
        stop(sprintf(
          'Internal default for PhysiCell %s has an invalid length.',
          column
        ))
      }
    }
  }
  developmental_pseudotime <- optional_numeric(
    'developmental_pseudotime',
    lineage_pseudotime
  )
  if(any(developmental_pseudotime < 0) ||
     any(developmental_pseudotime > 1)){
    stop('PhysiCell developmental_pseudotime must lie in [0, 1].')
  }
  metadata <- data.frame(
    sample_id = terminals$sample_id,
    physicell_id = terminals$physicell_id,
    node_id = terminals$node_id,
    cell_type = cell_types,
    lineage_depth = as.integer(lineage_depth),
    lineage_pseudotime = lineage_pseudotime,
    developmental_pseudotime = developmental_pseudotime,
    pseudotime = developmental_pseudotime,
    birth_time = as.numeric(terminals$birth_time),
    sampling_time = as.numeric(terminals$end_time),
    branch_length = as.numeric(terminals$branch_length),
    x = as.numeric(cell_rows$x),
    y = as.numeric(cell_rows$y),
    z = as.numeric(cell_rows$z),
    neighbor_count = unname(physicell_neighbor_count(cell_rows$neighbor_IDs)),
    culture_day = optional_numeric(
      'culture_day',
      as.numeric(terminals$end_time) / 1440
    ),
    state_start_time = optional_numeric(
      'state_start_time',
      as.numeric(terminals$birth_time)
    ),
    time_in_state = optional_numeric(
      'time_in_state',
      as.numeric(terminals$end_time) - as.numeric(terminals$birth_time)
    ),
    transition_count = optional_numeric('transition_count', 0),
    oxygen = optional_numeric('oxygen', 0),
    nutrient = optional_numeric('nutrient', 0),
    stringsAsFactors = FALSE
  )
  metadata$founder_id <- if('founder_ID' %in% names(cell_rows)){
    as.character(cell_rows$founder_ID)
  } else{
    NA_character_
  }
  metadata$tumor_radius <- sqrt(
    metadata$x^2 + metadata$y^2 + metadata$z^2
  )
  metadata$radial_position <- metadata$tumor_radius
  if(any(!is.finite(as.matrix(metadata[, c(
    'birth_time', 'sampling_time', 'branch_length', 'x', 'y', 'z',
    'tumor_radius'
  )])))){
    stop('PhysiCell lineage/spatial covariates must be finite.')
  }

  metadata$barcode_edit_count <- 0L
  metadata$barcode_edit_fraction <- 0
  barcode_path <- resolve_physicell_csv_path(
    file.path(recording_dir, 'barcode_binary_score_matrix.csv'),
    required = FALSE
  )
  sparse_barcode_path <- file.path(
    recording_dir,
    'barcode_binary_score_matrix_sparse.rds'
  )
  if(file.exists(sparse_barcode_path)){
    if(!requireNamespace('Matrix', quietly = TRUE)){
      stop('The Matrix package is required to read sparse barcode output.')
    }
    barcode_scores <- readRDS(sparse_barcode_path)
    if(is.null(dim(barcode_scores)) || is.null(rownames(barcode_scores))){
      stop('Sparse barcode matrix must be a row-named matrix.')
    }
  } else if(!is.na(barcode_path)){
    barcode_scores <- as.matrix(read_physicell_csv(
      barcode_path,
      row.names = 1,
      check.names = FALSE
    ))
    storage.mode(barcode_scores) <- 'numeric'
  } else{
    barcode_scores <- NULL
  }
  if(!is.null(barcode_scores)){
    barcode_indices <- match(metadata$sample_id, rownames(barcode_scores))
    if(anyNA(barcode_indices)){
      stop('Barcode matrix is missing one or more terminal sample IDs.')
    }
    selected_scores <- barcode_scores[barcode_indices, , drop = FALSE]
    edited_scores <- selected_scores != 0
    score_row_sums <- if(inherits(edited_scores, 'Matrix')){
      Matrix::rowSums(edited_scores)
    } else{
      rowSums(edited_scores)
    }
    metadata$barcode_edit_count <- as.integer(score_row_sums)
    metadata$barcode_edit_fraction <- if(ncol(selected_scores) == 0){
      0
    } else if(inherits(edited_scores, 'Matrix')){
      Matrix::rowMeans(edited_scores)
    } else{
      rowMeans(edited_scores)
    }
  }

  metadata$mt_variant_count <- 0L
  metadata$mt_heteroplasmy_burden <- 0
  metadata$mt_mean_variant_fraction <- 0
  metadata$mt_max_variant_fraction <- 0
  mitochondrial_path <- resolve_physicell_csv_path(
    file.path(
      recording_dir,
      'mitochondrial_variant_fractions.csv'
    ),
    required = FALSE
  )
  if(!is.na(mitochondrial_path)){
    mitochondrial <- read_physicell_csv(
      mitochondrial_path,
      colClasses = c(sample_id = 'character'),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if(nrow(mitochondrial) > 0){
      if(!all(c('sample_id', 'variant_fraction') %in% names(mitochondrial))){
        stop('Mitochondrial variant-fraction table has an invalid schema.')
      }
      mitochondrial$variant_fraction <- as.numeric(
        mitochondrial$variant_fraction
      )
      if(any(!is.finite(mitochondrial$variant_fraction)) ||
         any(mitochondrial$variant_fraction <= 0) ||
         any(mitochondrial$variant_fraction > 1)){
        stop('Mitochondrial variant fractions must be in (0, 1].')
      }
      by_sample <- split(
        mitochondrial$variant_fraction,
        mitochondrial$sample_id
      )
      observed_samples <- intersect(names(by_sample), metadata$sample_id)
      indices <- match(observed_samples, metadata$sample_id)
      metadata$mt_variant_count[indices] <- vapply(
        by_sample[observed_samples],
        length,
        integer(1)
      )
      metadata$mt_heteroplasmy_burden[indices] <- vapply(
        by_sample[observed_samples],
        sum,
        numeric(1)
      )
      metadata$mt_mean_variant_fraction[indices] <- vapply(
        by_sample[observed_samples],
        mean,
        numeric(1)
      )
      metadata$mt_max_variant_fraction[indices] <- vapply(
        by_sample[observed_samples],
        max,
        numeric(1)
      )
    }
  }

  metadata$ecdna_copy_number <- 0L
  metadata$ecdna_labeled_copy_number <- 0L
  metadata$ecdna_species_count <- 0L
  metadata$ecdna_labeled_species_count <- 0L
  metadata$ecdna_recorder_edit_fraction <- 0
  ecdna_path <- resolve_physicell_csv_path(
    file.path(recording_dir, 'ecdna_cell_summary.csv'),
    required = FALSE
  )
  if(!is.na(ecdna_path)){
    ecdna <- read_physicell_csv(
      ecdna_path,
      colClasses = c(sample_id = 'character'),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    required_ecdna_columns <- c(
      'sample_id', 'total_ecdna_copies', 'labeled_ecdna_copies',
      'observed_ecdna_species', 'observed_labeled_species',
      'ecdna_recorder_edit_fraction'
    )
    if(!all(required_ecdna_columns %in% names(ecdna))){
      stop('ecDNA cell-summary table has an invalid schema.')
    }
    numeric_columns <- setdiff(required_ecdna_columns, 'sample_id')
    ecdna[numeric_columns] <- lapply(ecdna[numeric_columns], as.numeric)
    if(any(!is.finite(as.matrix(ecdna[numeric_columns]))) ||
       any(as.matrix(ecdna[numeric_columns]) < 0) ||
       any(ecdna$ecdna_recorder_edit_fraction > 1)){
      stop('ecDNA covariates must be finite, non-negative, and fractions <= 1.')
    }
    observed <- intersect(ecdna$sample_id, metadata$sample_id)
    source_indices <- match(observed, ecdna$sample_id)
    target_indices <- match(observed, metadata$sample_id)
    metadata$ecdna_copy_number[target_indices] <- as.integer(
      ecdna$total_ecdna_copies[source_indices]
    )
    metadata$ecdna_labeled_copy_number[target_indices] <- as.integer(
      ecdna$labeled_ecdna_copies[source_indices]
    )
    metadata$ecdna_species_count[target_indices] <- as.integer(
      ecdna$observed_ecdna_species[source_indices]
    )
    metadata$ecdna_labeled_species_count[target_indices] <- as.integer(
      ecdna$observed_labeled_species[source_indices]
    )
    metadata$ecdna_recorder_edit_fraction[target_indices] <-
      ecdna$ecdna_recorder_edit_fraction[source_indices]
  }
  metadata
}

#' Rename simulated cell types to their reference-model equivalents
#'
#' Simulated PhysiCell labels rarely match the reference's cell-type levels, and
#' `scdesign3_new_covariates()` rejects levels the reference does not know. This
#' rewrites `cell_type` in place through the supplied map while preserving the
#' original label in `cell_type_sim`. Labels absent from the map are left
#' unchanged.
#'
#' @param metadata Covariate data frame carrying a `cell_type` column.
#' @param celltype_map Named vector, or named list flattened with `unlist()`,
#'   mapping simulated label to reference label; every element must be named.
#'   `NULL` returns `metadata` untouched with no `cell_type_sim` column added.
#' @return `metadata` with `cell_type` remapped and a `cell_type_sim` column
#'   holding the pre-mapping labels for every row.
apply_physicell_celltype_map <- function(metadata, celltype_map){
  if(is.null(celltype_map)){
    return(metadata)
  }
  if(is.list(celltype_map)){
    mapped_values <- unlist(celltype_map, use.names = TRUE)
  } else{
    mapped_values <- celltype_map
  }
  if(is.null(names(mapped_values)) || any(!nzchar(names(mapped_values)))){
    stop('The cell-type map must be a named JSON object.')
  }
  source_types <- metadata$cell_type
  replacements <- unname(mapped_values[source_types])
  replace <- !is.na(replacements)
  metadata$cell_type_sim <- source_types
  metadata$cell_type[replace] <- as.character(replacements[replace])
  metadata
}

#' Linearly rescale a simulated covariate onto the reference's observed range
#'
#' Maps the min/max of `values` onto the min/max of `reference_values` so the
#' simulated predictor stays inside the support the marginal model was fitted
#' over, avoiding extrapolation of the fitted smooth terms.
#'
#' @param values Simulated covariate; must be finite after numeric coercion.
#' @param reference_values Reference covariate defining the target range; must
#'   be finite after numeric coercion.
#' @return A numeric vector as long as `values`. When the reference range is
#'   degenerate every element is that single reference value; when the input
#'   range is degenerate every element is the midpoint of the reference range.
scale_physicell_covariate_to_reference <- function(values, reference_values){
  values <- as.numeric(values)
  reference_values <- as.numeric(reference_values)
  if(any(!is.finite(values)) || any(!is.finite(reference_values))){
    stop('Covariate scaling requires finite numeric values.')
  }
  reference_range <- range(reference_values)
  value_range <- range(values)
  if(diff(reference_range) == 0){
    return(rep(reference_range[1], length(values)))
  }
  if(diff(value_range) == 0){
    return(rep(mean(reference_range), length(values)))
  }
  reference_range[1] +
    (values - value_range[1]) / diff(value_range) * diff(reference_range)
}

#' Rescale simulated pseudotime and spatial predictors to reference ranges
#'
#' Reads the reference `colData` and rewrites the requested predictors with
#' `scale_physicell_covariate_to_reference()`. The reference must already be
#' standardized by `load_scdesign3_reference()`, which is what creates the
#' `pseudotime`, `spatial1`, and `spatial2` columns this reads.
#'
#' @param metadata Covariate data frame from `build_physicell_sc_covariates()`.
#' @param reference_sce Standardized reference `SingleCellExperiment` (rows are
#'   genes, columns are cells); only its `colData` is used.
#' @param use_pseudotime When `TRUE`, rescale `metadata$pseudotime` onto the
#'   reference `pseudotime` range.
#' @param use_spatial When `TRUE`, derive `spatial1` from `metadata$x` and
#'   `spatial2` from `metadata$y`, each rescaled onto the matching reference
#'   coordinate range.
#' @return `metadata` with the rescaled `pseudotime` and/or newly added
#'   `spatial1`/`spatial2` columns. With both flags `FALSE` it is returned
#'   unchanged.
align_physicell_sc_covariates <- function(metadata,
                                          reference_sce,
                                          use_pseudotime = FALSE,
                                          use_spatial = FALSE){
  reference_metadata <- as.data.frame(
    SummarizedExperiment::colData(reference_sce)
  )
  if(isTRUE(use_pseudotime)){
    metadata$pseudotime <- scale_physicell_covariate_to_reference(
      metadata$pseudotime,
      reference_metadata$pseudotime
    )
  }
  if(isTRUE(use_spatial)){
    metadata$spatial1 <- scale_physicell_covariate_to_reference(
      metadata$x,
      reference_metadata$spatial1
    )
    metadata$spatial2 <- scale_physicell_covariate_to_reference(
      metadata$y,
      reference_metadata$spatial2
    )
  }
  metadata
}

#' Write the simulated expression matrix, metadata, and run manifest
#'
#' Packages the simulated counts and their covariates into a
#' `SingleCellExperiment` and writes the artifacts the downstream analysis
#' scripts read. Does not require `scDesign3` itself, only the Bioconductor
#' container packages.
#'
#' @param counts Simulated count matrix in genes x cells orientation, as
#'   returned by `simulate_scdesign3_counts()`: rows are genes (row names from
#'   the reference), columns are cells named by `sample_id`. It is converted to
#'   a sparse `Matrix` when the `Matrix` package is available.
#' @param metadata Covariate data frame with one row per column of `counts` and
#'   a `sample_id` column, which becomes the `colData` row names.
#' @param output_dir Destination directory, created recursively if needed.
#' @param fit The scDesign3 fit list; `scdesign3_version`, `mu_formula`,
#'   `use_pseudotime`, `use_spatial`, and `other_covariates` are recorded in the
#'   manifest.
#' @param reference_path Path to the reference dataset, normalized into the
#'   manifest as provenance; must exist.
#' @param seed Random seed recorded in the manifest for provenance only — this
#'   function does not set or consume it.
#' @param compress_csv When `TRUE` (the default) the CSV outputs are gzipped.
#' @return Invisibly, the assembled `SingleCellExperiment`.
#' @section Side effects:
#' Creates `output_dir` and writes `sim_sce_final.rds` (the
#' `SingleCellExperiment`), `simulated_counts.rds` (the count matrix alone),
#' `cell_metadata.csv`, `cell_count_summary.csv` (cells per `cell_type`), and
#' `scdesign3_manifest.csv`, where `num_cells` is `ncol(counts)` and
#' `num_features` is `nrow(counts)`.
write_physicell_scdesign3_outputs <- function(counts,
                                               metadata,
                                               output_dir,
                                               fit,
                                               reference_path,
                                               seed,
                                               compress_csv = TRUE){
  if(!requireNamespace('SingleCellExperiment', quietly = TRUE) ||
     !requireNamespace('S4Vectors', quietly = TRUE)){
    stop('SingleCellExperiment and S4Vectors are required for output.')
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_counts <- if(requireNamespace('Matrix', quietly = TRUE)){
    Matrix::Matrix(counts, sparse = TRUE)
  } else{
    counts
  }
  simulated_sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = output_counts),
    colData = S4Vectors::DataFrame(
      metadata,
      row.names = metadata$sample_id
    )
  )
  saveRDS(
    simulated_sce,
    file.path(output_dir, 'sim_sce_final.rds')
  )
  saveRDS(
    output_counts,
    file.path(output_dir, 'simulated_counts.rds')
  )
  write_physicell_csv(
    metadata,
    file.path(output_dir, 'cell_metadata.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  cell_count_summary <- stats::aggregate(
    sample_id ~ cell_type,
    data = metadata,
    FUN = length
  )
  names(cell_count_summary)[2] <- 'n_cells'
  write_physicell_csv(
    cell_count_summary,
    file.path(output_dir, 'cell_count_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  manifest <- data.frame(
    property = c(
      'num_cells',
      'num_features',
      'reference_path',
      'scdesign3_version',
      'mu_formula',
      'use_pseudotime',
      'use_spatial',
      'other_covariates',
      'random_seed'
    ),
    value = c(
      ncol(counts),
      nrow(counts),
      normalizePath(reference_path, mustWork = TRUE),
      fit$scdesign3_version,
      fit$mu_formula,
      fit$use_pseudotime,
      fit$use_spatial,
      paste(fit$other_covariates, collapse = ','),
      seed
    ),
    stringsAsFactors = FALSE
  )
  write_physicell_csv(
    manifest,
    file.path(output_dir, 'scdesign3_manifest.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  invisible(simulated_sce)
}
