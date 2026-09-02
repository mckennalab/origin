# Visium-style sampling and spatial lineage analysis for PhysiCell tumors.
#
# The functions in this file operate on the event-resolved lineage tables
# written by simulate_physicell_lineage.R. They intentionally do not require a
# reconstructed tree: repeated MRCA queries use a binary-lifting parent index.

# ---- Input validation and section-plane geometry ----

#' Coerce and validate a numeric option supplied on the command line or in R
#'
#' Accepts a length-one comma-delimited string (the form the CLI wrappers pass),
#' a list, or a numeric vector, and flattens it to a plain numeric vector.
#' Values that fail to parse are rejected rather than silently dropped.
#'
#' @param value Numeric vector, list, or a single comma-delimited character
#'   string such as `"0,20,40,Inf"`.
#' @param name Human-readable option name used in the error message.
#' @param minimum_length Smallest acceptable number of parsed values.
#' @param allow_infinite When `TRUE`, infinite values are accepted (used for the
#'   open-ended top distance-bin edge); otherwise every value must be finite.
#' @return An unnamed numeric vector of length at least `minimum_length`.
visium_numeric_vector <- function(value,
                                  name,
                                  minimum_length = 1L,
                                  allow_infinite = FALSE){
  if(is.character(value) && length(value) == 1L){
    value <- trimws(unlist(strsplit(value, ',', fixed = TRUE)))
  }
  result <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
  valid <- !is.na(result)
  if(!isTRUE(allow_infinite)){
    valid <- valid & is.finite(result)
  }
  if(length(result) < minimum_length || !all(valid)){
    stop(sprintf('%s must contain at least %d valid numeric value(s).',
                 name, minimum_length))
  }
  result
}

#' Build an orthonormal basis for a section plane from its normal
#'
#' The returned `u` and `v` axes span the plane; together with `normal` they
#' form a right-handed orthonormal frame used to project PhysiCell cell centres
#' (microns in the simulation frame) into section coordinates. The in-plane axes
#' are otherwise arbitrary: `u` comes from crossing the normal with a helper
#' axis chosen to avoid degeneracy when the normal is near vertical.
#'
#' @param normal Three numbers giving the section normal in simulation
#'   coordinates; need not be unit length but must be nonzero.
#' @return A named list of three length-3 unit vectors, `normal`, `u`, and `v`.
visium_plane_basis <- function(normal = c(0, 0, 1)){
  normal <- visium_numeric_vector(normal, 'plane normal', minimum_length = 3L)
  if(length(normal) != 3L || sqrt(sum(normal^2)) == 0){
    stop('plane normal must contain exactly three values and have nonzero length.')
  }
  normal <- normal / sqrt(sum(normal^2))
  helper <- if(abs(normal[3]) < 0.9) c(0, 0, 1) else c(0, 1, 0)
  plane_u <- c(
    helper[2] * normal[3] - helper[3] * normal[2],
    helper[3] * normal[1] - helper[1] * normal[3],
    helper[1] * normal[2] - helper[2] * normal[1]
  )
  plane_u <- plane_u / sqrt(sum(plane_u^2))
  plane_v <- c(
    normal[2] * plane_u[3] - normal[3] * plane_u[2],
    normal[3] * plane_u[1] - normal[1] * plane_u[3],
    normal[1] * plane_u[2] - normal[2] * plane_u[1]
  )
  list(normal = normal, u = plane_u, v = plane_v)
}

# ---- Visium capture-array geometry ----

#' Lay out a Space Ranger-style Visium capture array in section coordinates
#'
#' Spots sit on a hexagonal lattice in the section plane. `array_row` counts
#' rows from 0 and `array_col` is `2 * within_row + (array_row %% 2)`, so even
#' rows carry even column indices and odd rows odd ones, matching how Space
#' Ranger numbers a 6.5 mm slide. Centres are `spot_pitch` microns apart along a
#' row (two `array_col` steps) and `spot_pitch * sqrt(3) / 2` microns apart
#' between rows, which puts every nearest neighbour exactly `spot_pitch` microns
#' away. The lattice is centred on its own bounding box before being rotated and
#' translated, so `translation` displaces the array from the section origin.
#'
#' @param spot_pitch Centre-to-centre spacing in microns; one positive finite
#'   number.
#' @param rotation_degrees Rotation of the array within the section plane; one
#'   finite number.
#' @param translation Two numbers, the array offset in microns along the plane
#'   `u` and `v` axes.
#' @param array_rows Number of spot rows; a positive integer (78 on a real
#'   6.5 mm slide).
#' @param spots_per_row Spots per row; a positive integer (64 on a real slide).
#' @param slice_id String embedded in the generated `barcode` and `spot_id`
#'   values; non-alphanumeric characters are stripped from the barcode.
#' @return A data frame with one row per spot (`array_rows * spots_per_row`
#'   rows, 4,992 with the defaults) and columns `barcode`, `spot_id`,
#'   `in_tissue` (0 until spots are populated), `array_row`, `array_col`,
#'   `plane_u` and `plane_v` (spot centres in section microns), and `n_cells`
#'   (0).
make_visium_6_5mm_array <- function(spot_pitch = 100,
                                    rotation_degrees = 0,
                                    translation = c(0, 0),
                                    array_rows = 78L,
                                    spots_per_row = 64L,
                                    slice_id = 'slice_001'){
  spot_pitch <- as.numeric(spot_pitch)
  rotation_degrees <- as.numeric(rotation_degrees)
  translation <- visium_numeric_vector(
    translation,
    'Visium array translation',
    minimum_length = 2L
  )
  array_rows <- as.integer(array_rows)
  spots_per_row <- as.integer(spots_per_row)
  if(length(spot_pitch) != 1L || !is.finite(spot_pitch) || spot_pitch <= 0){
    stop('spot_pitch must be one positive finite number.')
  }
  if(length(rotation_degrees) != 1L || !is.finite(rotation_degrees)){
    stop('rotation_degrees must be one finite number.')
  }
  if(length(translation) != 2L){
    stop('Visium array translation must contain exactly two values.')
  }
  if(length(array_rows) != 1L || is.na(array_rows) || array_rows < 1L ||
     length(spots_per_row) != 1L || is.na(spots_per_row) ||
     spots_per_row < 1L){
    stop('array_rows and spots_per_row must be positive integers.')
  }

  rows <- rep(seq.int(0L, array_rows - 1L), each = spots_per_row)
  within_row <- rep(seq.int(0L, spots_per_row - 1L), times = array_rows)
  columns <- 2L * within_row + (rows %% 2L)
  raw_u <- columns * spot_pitch / 2
  raw_v <- rows * spot_pitch * sqrt(3) / 2
  raw_u <- raw_u - mean(range(raw_u))
  raw_v <- raw_v - mean(range(raw_v))
  angle <- rotation_degrees * pi / 180
  plane_u <- cos(angle) * raw_u - sin(angle) * raw_v + translation[1]
  plane_v <- sin(angle) * raw_u + cos(angle) * raw_v + translation[2]
  data.frame(
    barcode = sprintf(
      'SIM-%s-r%02d-c%03d',
      gsub('[^A-Za-z0-9]+', '', slice_id),
      rows,
      columns
    ),
    spot_id = sprintf('%s_r%02d_c%03d', slice_id, rows, columns),
    in_tissue = 0L,
    array_row = rows,
    array_col = columns,
    plane_u = plane_u,
    plane_v = plane_v,
    n_cells = 0L,
    stringsAsFactors = FALSE
  )
}

# ---- PhysiCell run input and sectioning ----

#' Load a completed PhysiCell run's live cells and event-resolved lineage
#'
#' Reads `physicell_build/output/lineage_table.csv` (PhysiCell's final live-cell
#' table) together with `lineage_recording/terminal_cells.csv` and
#' `lineage_recording/lineage_nodes.csv`, then joins every live cell to the
#' lineage tip it ended at. Requires physicell_lineage.R to be sourced for
#' `read_physicell_csv()` and `physicell_id()`; plain and `.gz` CSVs are both
#' accepted. IDs must be unique on both sides and every live cell must have a
#' terminal node, or the function stops.
#'
#' @param run_dir Directory of a completed PhysiCell pipeline run; must exist.
#' @return A named list with `run_dir` (normalised path), `cells` (one row per
#'   live cell with `sample_id`, `physicell_id`, `node_id`, and finite `x`, `y`,
#'   `z` in simulation microns, plus whichever of `parent_ID`, `founder_ID`,
#'   `cell_type`, and `alive` the position table carried), `nodes` (the lineage
#'   node table as read), `terminal_cells`, and `recording_dir`.
read_physicell_visium_inputs <- function(run_dir){
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  position_path <- file.path(
    run_dir,
    'physicell_build',
    'output',
    'lineage_table.csv'
  )
  recording_dir <- file.path(run_dir, 'lineage_recording')
  positions <- read_physicell_csv(position_path, check.names = FALSE)
  terminal_cells <- read_physicell_csv(
    file.path(recording_dir, 'terminal_cells.csv'),
    check.names = FALSE
  )
  nodes <- read_physicell_csv(
    file.path(recording_dir, 'lineage_nodes.csv'),
    check.names = FALSE
  )
  position_required <- c('ID', 'x', 'y', 'z')
  terminal_required <- c('sample_id', 'physicell_id', 'node_id')
  node_required <- c(
    'node_id', 'parent_node_id', 'birth_time', 'end_time', 'branch_length'
  )
  if(!all(position_required %in% names(positions))){
    stop('PhysiCell lineage_table.csv must contain ID, x, y, and z.')
  }
  if(!all(terminal_required %in% names(terminal_cells))){
    stop('terminal_cells.csv is missing sample_id, physicell_id, or node_id.')
  }
  if(!all(node_required %in% names(nodes))){
    stop('lineage_nodes.csv is missing required lineage columns.')
  }

  positions$physicell_id <- physicell_id(positions$ID, 'position ID')
  terminal_cells$physicell_id <- physicell_id(
    terminal_cells$physicell_id,
    'terminal physicell_id'
  )
  if(anyDuplicated(positions$physicell_id)){
    stop('PhysiCell lineage_table.csv contains duplicated live-cell IDs.')
  }
  if(anyDuplicated(terminal_cells$physicell_id)){
    stop('terminal_cells.csv contains duplicated terminal PhysiCell IDs.')
  }
  terminal_index <- match(positions$physicell_id, terminal_cells$physicell_id)
  if(anyNA(terminal_index)){
    stop(sprintf(
      '%d live PhysiCell cells do not have matching terminal lineage nodes.',
      sum(is.na(terminal_index))
    ))
  }
  coordinate_matrix <- as.matrix(positions[, c('x', 'y', 'z'), drop = FALSE])
  storage.mode(coordinate_matrix) <- 'double'
  if(any(!is.finite(coordinate_matrix))){
    stop('PhysiCell live-cell coordinates must be finite numbers.')
  }
  cells <- data.frame(
    sample_id = as.character(terminal_cells$sample_id[terminal_index]),
    physicell_id = positions$physicell_id,
    node_id = as.character(terminal_cells$node_id[terminal_index]),
    x = coordinate_matrix[, 1],
    y = coordinate_matrix[, 2],
    z = coordinate_matrix[, 3],
    stringsAsFactors = FALSE
  )
  optional_columns <- intersect(
    c('parent_ID', 'founder_ID', 'cell_type', 'alive'),
    names(positions)
  )
  if(length(optional_columns) > 0L){
    cells <- cbind(cells, positions[, optional_columns, drop = FALSE])
  }
  list(
    run_dir = run_dir,
    cells = cells,
    nodes = nodes,
    terminal_cells = terminal_cells,
    recording_dir = recording_dir
  )
}

#' Cut a finite-thickness section through the tumour and project it
#'
#' Cell centres are translated by `center`, projected onto the plane basis, and
#' kept when their signed distance along the normal lies within `thickness / 2`
#' of `offset`. All lengths are the same microns as the PhysiCell coordinates.
#' Only cell centres are tested, so a cell whose body straddles the section but
#' whose centre does not is dropped.
#'
#' @param cells Data frame with numeric `x`, `y`, and `z` columns.
#' @param basis Plane basis from `visium_plane_basis()`: a list with unit
#'   `normal`, `u`, and `v`.
#' @param center Three coordinates giving the point the section is measured
#'   from, in simulation microns (typically the tumour centroid).
#' @param offset Signed distance of this section from `center` along the normal,
#'   in microns; one finite number.
#' @param thickness Section thickness in microns; one positive finite number.
#' @return The subset of `cells` intersecting the section, row names reset, with
#'   four added columns in microns: `plane_u` and `plane_v` (in-section
#'   coordinates relative to `center`), `plane_w` (signed distance from `center`
#'   along the normal), and `distance_from_section_midplane`
#'   (`plane_w - offset`).
slice_physicell_cells <- function(cells,
                                  basis,
                                  center,
                                  offset = 0,
                                  thickness = 5){
  if(!is.data.frame(cells) || !all(c('x', 'y', 'z') %in% names(cells))){
    stop('cells must contain x, y, and z columns.')
  }
  center <- visium_numeric_vector(center, 'slice center', minimum_length = 3L)
  if(length(center) != 3L){
    stop('slice center must contain exactly three coordinates.')
  }
  offset <- as.numeric(offset)
  thickness <- as.numeric(thickness)
  if(length(offset) != 1L || !is.finite(offset)){
    stop('slice offset must be one finite number.')
  }
  if(length(thickness) != 1L || !is.finite(thickness) || thickness <= 0){
    stop('slice thickness must be one positive finite number.')
  }
  coordinates <- as.matrix(cells[, c('x', 'y', 'z'), drop = FALSE])
  centered <- sweep(coordinates, 2, center, '-')
  plane_u <- as.numeric(centered %*% basis$u)
  plane_v <- as.numeric(centered %*% basis$v)
  plane_w <- as.numeric(centered %*% basis$normal)
  keep <- abs(plane_w - offset) <= thickness / 2
  result <- cells[keep, , drop = FALSE]
  result$plane_u <- plane_u[keep]
  result$plane_v <- plane_v[keep]
  result$plane_w <- plane_w[keep]
  result$distance_from_section_midplane <- plane_w[keep] - offset
  rownames(result) <- NULL
  result
}

#' Assign section cells to their nearest capture spot
#'
#' For each section cell the nearest spot centre is found in the 2-D section
#' plane, and the cell is captured when that distance is at most
#' `spot_diameter / 2`. Cells landing in the gaps between spot footprints are
#' kept but flagged uncaptured. Both tables are returned as modified copies.
#'
#' @param slice_cells Section cells from `slice_physicell_cells()`, carrying
#'   `plane_u` and `plane_v` in microns.
#' @param spots Spot table from `make_visium_6_5mm_array()`; its `n_cells` and
#'   `in_tissue` columns are recomputed from scratch.
#' @param spot_diameter Capture-spot footprint diameter in microns; one positive
#'   finite number (55 on a real Visium slide).
#' @return A named list with `cells` (`slice_cells` plus `spot_id`,
#'   `spot_barcode`, `spot_array_row`, `spot_array_col`, `spot_center_u`,
#'   `spot_center_v`, `distance_to_spot_center`, and logical `captured`, where
#'   the spot columns are `NA` for uncaptured cells) and `spots` (the same rows
#'   with `n_cells` set to the captured count and `in_tissue` set to 1 wherever
#'   that count is positive).
#' @note The nearest-spot search is an explicit per-cell loop that scans every
#'   spot, so the cost is O(section cells x array spots); with the default
#'   4,992-spot array this is the slowest step for a large section.
assign_cells_to_visium_spots <- function(slice_cells,
                                         spots,
                                         spot_diameter = 55){
  spot_diameter <- as.numeric(spot_diameter)
  if(length(spot_diameter) != 1L || !is.finite(spot_diameter) ||
     spot_diameter <= 0){
    stop('spot_diameter must be one positive finite number.')
  }
  spots$n_cells <- integer(nrow(spots))
  spots$in_tissue <- integer(nrow(spots))
  slice_cells$spot_id <- rep(NA_character_, nrow(slice_cells))
  slice_cells$spot_barcode <- rep(NA_character_, nrow(slice_cells))
  slice_cells$spot_array_row <- rep(NA_integer_, nrow(slice_cells))
  slice_cells$spot_array_col <- rep(NA_integer_, nrow(slice_cells))
  slice_cells$spot_center_u <- rep(NA_real_, nrow(slice_cells))
  slice_cells$spot_center_v <- rep(NA_real_, nrow(slice_cells))
  slice_cells$distance_to_spot_center <- rep(NA_real_, nrow(slice_cells))
  slice_cells$captured <- rep(FALSE, nrow(slice_cells))
  if(nrow(slice_cells) == 0L || nrow(spots) == 0L){
    return(list(cells = slice_cells, spots = spots))
  }
  nearest_index <- integer(nrow(slice_cells))
  nearest_distance <- numeric(nrow(slice_cells))
  for(cell_index in seq_len(nrow(slice_cells))){
    squared_distance <-
      (spots$plane_u - slice_cells$plane_u[cell_index])^2 +
      (spots$plane_v - slice_cells$plane_v[cell_index])^2
    nearest_index[cell_index] <- which.min(squared_distance)
    nearest_distance[cell_index] <- sqrt(squared_distance[nearest_index[cell_index]])
  }
  captured <- nearest_distance <= spot_diameter / 2
  assigned <- nearest_index[captured]
  if(length(assigned) > 0L){
    counts <- tabulate(assigned, nbins = nrow(spots))
    spots$n_cells <- counts
    spots$in_tissue <- as.integer(counts > 0L)
    slice_cells$spot_id[captured] <- spots$spot_id[assigned]
    slice_cells$spot_barcode[captured] <- spots$barcode[assigned]
    slice_cells$spot_array_row[captured] <- spots$array_row[assigned]
    slice_cells$spot_array_col[captured] <- spots$array_col[assigned]
    slice_cells$spot_center_u[captured] <- spots$plane_u[assigned]
    slice_cells$spot_center_v[captured] <- spots$plane_v[assigned]
    slice_cells$distance_to_spot_center[captured] <- nearest_distance[captured]
    slice_cells$captured[captured] <- TRUE
  }
  list(cells = slice_cells, spots = spots)
}

# ---- Lineage MRCA index and pair queries ----

#' Build a binary-lifting ancestor index over the event-resolved lineage
#'
#' Converts the lineage node table into integer arrays so repeated MRCA queries
#' cost O(log depth) each without reconstructing a tree object. Nodes must
#' already be topologically ordered (every parent before its children) and node
#' IDs unique; roots are the rows whose `parent_node_id` is missing, empty, or
#' the literal string `NA`. Times stay in the lineage table's own units
#' (PhysiCell minutes) and are not rescaled here.
#'
#' @param nodes Data frame with `node_id`, `parent_node_id`, `birth_time`,
#'   `end_time`, and `branch_length` columns.
#' @return A named list: `node_ids` (character), `node_lookup` (named integer
#'   vector mapping node ID to row index), `parent` (integer row index, 0 for a
#'   root), `ancestors` (integer matrix of nodes x levels, where
#'   `ancestors[i, k]` is the 2^(k-1)-th ancestor of node `i` and 0 means the
#'   jump runs past a root), `depth` (0 at each root), `root` (row index of each
#'   node's founder), `birth_time`, `end_time`, and `synthetic_root_time` (the
#'   earliest root birth time, used as the notional MRCA time for cells from
#'   different founders).
build_visium_lca_index <- function(nodes){
  required <- c(
    'node_id', 'parent_node_id', 'birth_time', 'end_time', 'branch_length'
  )
  if(!is.data.frame(nodes) || !all(required %in% names(nodes))){
    stop('nodes must contain node, parent, and branch-time columns.')
  }
  node_ids <- as.character(nodes$node_id)
  if(anyNA(node_ids) || any(!nzchar(node_ids)) || anyDuplicated(node_ids)){
    stop('lineage node IDs must be unique non-empty strings.')
  }
  parent_ids <- as.character(nodes$parent_node_id)
  parent_ids[is.na(nodes$parent_node_id) | parent_ids %in% c('', 'NA')] <- NA
  parent <- match(parent_ids, node_ids)
  unresolved <- !is.na(parent_ids) & is.na(parent)
  if(any(unresolved)){
    stop('Every non-root lineage parent must identify a node in the table.')
  }
  node_number <- seq_along(node_ids)
  if(any(!is.na(parent) & parent >= node_number)){
    stop('Lineage nodes must be topologically ordered with parents before children.')
  }
  birth_time <- suppressWarnings(as.numeric(nodes$birth_time))
  end_time <- suppressWarnings(as.numeric(nodes$end_time))
  if(any(!is.finite(birth_time)) || any(!is.finite(end_time)) ||
     any(end_time < birth_time)){
    stop('Lineage birth and end times must be finite and non-decreasing.')
  }
  depth <- integer(length(node_ids))
  root <- integer(length(node_ids))
  for(index in node_number){
    if(is.na(parent[index])){
      depth[index] <- 0L
      root[index] <- index
    } else{
      depth[index] <- depth[parent[index]] + 1L
      root[index] <- root[parent[index]]
    }
  }
  parent_integer <- parent
  parent_integer[is.na(parent_integer)] <- 0L
  levels <- max(1L, floor(log2(max(depth) + 1)) + 1L)
  ancestors <- matrix(0L, nrow = length(node_ids), ncol = levels)
  ancestors[, 1] <- parent_integer
  if(levels > 1L){
    for(level in 2:levels){
      previous <- ancestors[, level - 1L]
      present <- previous > 0L
      ancestors[present, level] <- ancestors[previous[present], level - 1L]
    }
  }
  list(
    node_ids = node_ids,
    node_lookup = setNames(node_number, node_ids),
    parent = parent_integer,
    ancestors = ancestors,
    depth = depth,
    root = root,
    birth_time = birth_time,
    end_time = end_time,
    synthetic_root_time = min(birth_time[parent_integer == 0L], na.rm = TRUE)
  )
}

#' Resolve MRCA identity, age, and distance for many lineage-node pairs
#'
#' Answers all pairs at once by binary lifting: the deeper node is raised to its
#' partner's depth, then both climb together while their ancestors differ. Pairs
#' whose nodes descend from different founders have no real MRCA; for those the
#' index's `synthetic_root_time` stands in for the MRCA time, `mrca_node_id` and
#' `division_distance` are `NA`, and `patristic_distance` is measured against
#' that synthetic root. Times stay in the lineage table's units, so the caller
#' divides by `time_units_per_hour` to reach hours.
#'
#' @param index Index from `build_visium_lca_index()`.
#' @param node_a Integer row indices into `index$node_ids`.
#' @param node_b Integer row indices, the same length as `node_a`.
#' @return A data frame with one row per input pair, in input order, holding
#'   `mrca_node_id`, `same_founder`, `mrca_time`, `mrca_age`
#'   (`patristic_distance / 2`), `patristic_distance` (the two tip end times
#'   less twice the MRCA time), and integer `division_distance` (the number of
#'   division edges separating the two tips). A zero-length input returns the
#'   same columns with no rows.
query_visium_lineage_pairs <- function(index, node_a, node_b){
  node_a <- as.integer(node_a)
  node_b <- as.integer(node_b)
  number_nodes <- length(index$node_ids)
  if(length(node_a) != length(node_b) || anyNA(node_a) || anyNA(node_b) ||
     any(node_a < 1L | node_a > number_nodes) ||
     any(node_b < 1L | node_b > number_nodes)){
    stop('node_a and node_b must be equal-length valid lineage node indices.')
  }
  pair_count <- length(node_a)
  if(pair_count == 0L){
    return(data.frame(
      mrca_node_id = character(), same_founder = logical(),
      mrca_time = numeric(), mrca_age = numeric(),
      patristic_distance = numeric(), division_distance = integer(),
      stringsAsFactors = FALSE
    ))
  }
  original_a <- node_a
  original_b <- node_b
  same_founder <- index$root[node_a] == index$root[node_b]
  mrca <- integer(pair_count)
  active_indices <- which(same_founder)
  if(length(active_indices) > 0L){
    left <- node_a[active_indices]
    right <- node_b[active_indices]
    swap <- index$depth[left] < index$depth[right]
    if(any(swap)){
      temporary <- left[swap]
      left[swap] <- right[swap]
      right[swap] <- temporary
    }
    difference <- index$depth[left] - index$depth[right]
    for(level in seq_len(ncol(index$ancestors))){
      move <- ((difference %/% (2^(level - 1L))) %% 2) == 1
      if(any(move)){
        left[move] <- index$ancestors[left[move], level]
      }
    }
    matched <- left == right
    local_mrca <- integer(length(left))
    local_mrca[matched] <- left[matched]
    unresolved <- which(!matched)
    if(length(unresolved) > 0L){
      for(level in rev(seq_len(ncol(index$ancestors)))){
        left_ancestor <- index$ancestors[left[unresolved], level]
        right_ancestor <- index$ancestors[right[unresolved], level]
        move <- left_ancestor != right_ancestor &
          left_ancestor > 0L & right_ancestor > 0L
        if(any(move)){
          selected <- unresolved[move]
          left[selected] <- left_ancestor[move]
          right[selected] <- right_ancestor[move]
        }
      }
      local_mrca[unresolved] <- index$parent[left[unresolved]]
    }
    mrca[active_indices] <- local_mrca
  }
  mrca_time <- rep(index$synthetic_root_time, pair_count)
  mrca_time[same_founder] <- index$end_time[mrca[same_founder]]
  end_a <- index$end_time[original_a]
  end_b <- index$end_time[original_b]
  patristic_distance <- end_a + end_b - 2 * mrca_time
  division_distance <- rep(NA_integer_, pair_count)
  division_distance[same_founder] <-
    index$depth[original_a[same_founder]] +
    index$depth[original_b[same_founder]] -
    2L * index$depth[mrca[same_founder]]
  mrca_node_id <- rep(NA_character_, pair_count)
  mrca_node_id[same_founder] <- index$node_ids[mrca[same_founder]]
  data.frame(
    mrca_node_id = mrca_node_id,
    same_founder = same_founder,
    mrca_time = mrca_time,
    mrca_age = patristic_distance / 2,
    patristic_distance = patristic_distance,
    division_distance = division_distance,
    stringsAsFactors = FALSE
  )
}

# ---- Cell-pair sampling and spatial binning ----

#' Enumerate or subsample unordered cell pairs
#'
#' At or below `max_exact_pairs` every unordered pair is enumerated exactly;
#' above it, `max_candidate_pairs` pairs are drawn uniformly at random *with
#' replacement*, so a pair may appear more than once and the enumeration is no
#' longer exhaustive. The two indices of a drawn pair are always distinct and
#' are returned in ascending order.
#'
#' @param number_cells Number of cells to pair; one non-negative integer.
#' @param max_exact_pairs Largest total pair count still enumerated
#'   exhaustively; a positive integer.
#' @param max_candidate_pairs Number of pairs drawn once that limit is exceeded;
#'   a positive integer.
#' @param seed Integer seed for the draw.
#' @return A named list with `pairs` (a 2-row integer matrix whose *columns* are
#'   the pairs: row 1 is the smaller cell index, row 2 the larger; it is 2 x 0
#'   when fewer than two cells are supplied), `total_pairs` (how many unordered
#'   pairs exist, as a double), and `sampling_mode` (`'no_pairs'`,
#'   `'all_pairs'`, or `'uniform_pairs_with_replacement'`).
#' @section Side effects:
#' Calls `set.seed(seed)`, replacing the global RNG state.
sample_visium_cell_pairs <- function(number_cells,
                                     max_exact_pairs = 2000000L,
                                     max_candidate_pairs = 2000000L,
                                     seed = 1L){
  number_cells <- as.integer(number_cells)
  max_exact_pairs <- as.integer(max_exact_pairs)
  max_candidate_pairs <- as.integer(max_candidate_pairs)
  if(length(number_cells) != 1L || is.na(number_cells) || number_cells < 0L){
    stop('number_cells must be one non-negative integer.')
  }
  if(anyNA(c(max_exact_pairs, max_candidate_pairs)) ||
     any(c(max_exact_pairs, max_candidate_pairs) < 1L)){
    stop('pair limits must be positive integers.')
  }
  total_pairs <- number_cells * (number_cells - 1) / 2
  if(number_cells < 2L){
    return(list(
      pairs = matrix(integer(), nrow = 2L),
      total_pairs = total_pairs,
      sampling_mode = 'no_pairs'
    ))
  }
  set.seed(as.integer(seed))
  if(total_pairs <= max_exact_pairs){
    pairs <- utils::combn(number_cells, 2L)
    sampling_mode <- 'all_pairs'
  } else{
    sample_size <- min(max_candidate_pairs, total_pairs)
    first <- sample.int(number_cells, sample_size, replace = TRUE)
    second <- sample.int(number_cells - 1L, sample_size, replace = TRUE)
    second <- second + as.integer(second >= first)
    pairs <- rbind(pmin(first, second), pmax(first, second))
    sampling_mode <- 'uniform_pairs_with_replacement'
  }
  list(
    pairs = pairs,
    total_pairs = total_pairs,
    sampling_mode = sampling_mode
  )
}

#' Format left-closed distance-bin labels
#'
#' Produces stable `[lower,upper)` labels that are reused as `cut()` levels, as
#' the `distance_bin` column of the pair tables, and as the join key between the
#' observed and permuted-null summaries, so every caller must generate them the
#' same way.
#'
#' @param breaks Strictly increasing distance edges in microns, at least two of
#'   them; the first must be finite while the last may be infinite.
#' @return A character vector of `length(breaks) - 1` labels.
visium_distance_labels <- function(breaks){
  breaks <- visium_numeric_vector(
    breaks,
    'distance breaks',
    minimum_length = 2L,
    allow_infinite = TRUE
  )
  if(is.infinite(breaks[1]) || any(diff(breaks) <= 0)){
    stop('distance breaks must be strictly increasing with a finite lower bound.')
  }
  #' Internal: Render one bin edge for a label
  #'
  #' @param value One numeric bin edge.
  #' @return `'Inf'` for an infinite edge, otherwise the trimmed non-scientific
  #'   decimal form.
  format_bound <- function(value){
    if(is.infinite(value)) 'Inf' else format(value, trim = TRUE, scientific = FALSE)
  }
  vapply(seq_len(length(breaks) - 1L), function(index){
    sprintf('[%s,%s)', format_bound(breaks[index]), format_bound(breaks[index + 1L]))
  }, character(1))
}

#' Sample captured-cell pairs stratified by in-section distance
#'
#' Draws candidate pairs with `sample_visium_cell_pairs()`, bins them on the 2-D
#' in-section separation (`plane_u`/`plane_v`, microns), then caps each bin at
#' `max_pairs_per_bin` so the rare near bins are not swamped by the abundant far
#' ones. Pairs falling outside the outermost break are dropped.
#' `candidate_counts` records occupancy before the cap so a summary can report
#' how many pairs a bin could have contributed.
#'
#' @param cells Captured cells; must carry `sample_id`, `physicell_id`,
#'   `node_index`, `x`, `y`, `z`, `plane_u`, `plane_v`, `spot_id`,
#'   `spot_center_u`, and `spot_center_v`.
#' @param distance_breaks Increasing bin edges in microns, the last optionally
#'   infinite.
#' @param max_pairs_per_bin Cap on the pairs retained per distance bin.
#' @param max_exact_pairs Passed to `sample_visium_cell_pairs()`.
#' @param max_candidate_pairs Passed to `sample_visium_cell_pairs()`.
#' @param seed Integer seed; the candidate draw uses `seed` and the per-bin
#'   thinning uses `seed + 1`.
#' @return A named list with `pairs` (one row per retained pair holding the two
#'   cell row indices, the sample/PhysiCell/node/spot identifiers of each side,
#'   `same_spot`, `spatial_distance_2d_microns`, `spatial_distance_3d_microns`,
#'   `spot_center_distance_microns`, and the `distance_bin` label),
#'   `total_pairs`, `sampling_mode`, `distance_labels`, and `candidate_counts`
#'   (integer, one per label, counted before thinning).
#' @section Side effects:
#' Calls `set.seed()` twice, replacing the global RNG state.
#' @note Binning uses the 2-D in-section distance; the 3-D distance is carried
#'   through for reference only.
prepare_visium_spatial_pairs <- function(cells,
                                         distance_breaks,
                                         max_pairs_per_bin = 100000L,
                                         max_exact_pairs = 2000000L,
                                         max_candidate_pairs = 2000000L,
                                         seed = 1L){
  required <- c(
    'sample_id', 'physicell_id', 'node_index', 'x', 'y', 'z',
    'plane_u', 'plane_v', 'spot_id', 'spot_center_u', 'spot_center_v'
  )
  if(!is.data.frame(cells) || !all(required %in% names(cells))){
    stop('captured cells are missing spatial, spot, or lineage columns.')
  }
  breaks <- visium_numeric_vector(
    distance_breaks,
    'distance breaks',
    minimum_length = 2L,
    allow_infinite = TRUE
  )
  labels <- visium_distance_labels(breaks)
  pair_sample <- sample_visium_cell_pairs(
    nrow(cells),
    max_exact_pairs = max_exact_pairs,
    max_candidate_pairs = max_candidate_pairs,
    seed = seed
  )
  pairs <- pair_sample$pairs
  if(ncol(pairs) == 0L){
    return(list(
      pairs = data.frame(
        cell_index_a = integer(), cell_index_b = integer(),
        sample_id_a = character(), sample_id_b = character(),
        physicell_id_a = character(), physicell_id_b = character(),
        node_index_a = integer(), node_index_b = integer(),
        spot_id_a = character(), spot_id_b = character(),
        same_spot = logical(), spatial_distance_2d_microns = numeric(),
        spatial_distance_3d_microns = numeric(),
        spot_center_distance_microns = numeric(), distance_bin = character(),
        stringsAsFactors = FALSE
      ),
      total_pairs = pair_sample$total_pairs,
      sampling_mode = pair_sample$sampling_mode,
      distance_labels = labels,
      candidate_counts = integer(length(labels))
    ))
  }
  first <- pairs[1, ]
  second <- pairs[2, ]
  distance_2d <- sqrt(
    (cells$plane_u[first] - cells$plane_u[second])^2 +
      (cells$plane_v[first] - cells$plane_v[second])^2
  )
  distance_3d <- sqrt(
    (cells$x[first] - cells$x[second])^2 +
      (cells$y[first] - cells$y[second])^2 +
      (cells$z[first] - cells$z[second])^2
  )
  distance_bin <- cut(
    distance_2d,
    breaks = breaks,
    labels = labels,
    right = FALSE,
    include.lowest = TRUE
  )
  candidate_counts <- tabulate(
    as.integer(distance_bin),
    nbins = length(labels)
  )
  keep <- integer()
  set.seed(as.integer(seed) + 1L)
  for(bin_index in seq_along(labels)){
    candidates <- which(as.integer(distance_bin) == bin_index)
    if(length(candidates) > max_pairs_per_bin){
      candidates <- sample(candidates, max_pairs_per_bin)
    }
    keep <- c(keep, candidates)
  }
  keep <- sort(keep)
  first <- first[keep]
  second <- second[keep]
  same_spot <- cells$spot_id[first] == cells$spot_id[second]
  spot_distance <- sqrt(
    (cells$spot_center_u[first] - cells$spot_center_u[second])^2 +
      (cells$spot_center_v[first] - cells$spot_center_v[second])^2
  )
  pair_table <- data.frame(
    cell_index_a = first,
    cell_index_b = second,
    sample_id_a = cells$sample_id[first],
    sample_id_b = cells$sample_id[second],
    physicell_id_a = cells$physicell_id[first],
    physicell_id_b = cells$physicell_id[second],
    node_index_a = cells$node_index[first],
    node_index_b = cells$node_index[second],
    spot_id_a = cells$spot_id[first],
    spot_id_b = cells$spot_id[second],
    same_spot = same_spot,
    spatial_distance_2d_microns = distance_2d[keep],
    spatial_distance_3d_microns = distance_3d[keep],
    spot_center_distance_microns = spot_distance,
    distance_bin = as.character(distance_bin[keep]),
    stringsAsFactors = FALSE
  )
  list(
    pairs = pair_table,
    total_pairs = pair_sample$total_pairs,
    sampling_mode = pair_sample$sampling_mode,
    distance_labels = labels,
    candidate_counts = candidate_counts
  )
}

# ---- Pair summaries and permutation null ----

#' Summarise lineage relatedness per spatial-distance bin
#'
#' Builds the cell-level correlogram: one row per distance bin holding the mean
#' spatial separation alongside mean and median MRCA age, patristic and division
#' distance, and the fractions of pairs sharing a founder or a spot. A bin with
#' no evaluated pairs still gets a row, filled with `NA`. When a permuted-null
#' table is supplied its columns are joined on `distance_bin` and the
#' observed-versus-null contrasts are appended.
#'
#' @param pair_table Pairs from `prepare_visium_spatial_pairs()` after the
#'   lineage columns (`mrca_age_hours`, `patristic_distance_hours`,
#'   `division_distance`, `same_founder`) have been attached.
#' @param distance_labels Bin labels, in the order rows should appear.
#' @param candidate_counts Pre-thinning pair counts, one per label.
#' @param recent_threshold_hours One or more MRCA-age cutoffs in hours; each
#'   adds a `recent_mrca_le_<threshold>h_fraction` column whose suffix has `.`
#'   replaced by `_`.
#' @param null_summaries Optional output of `permute_visium_lineage_null()`.
#' @return A data frame with one row per label and columns `distance_bin`,
#'   `candidate_pairs`, `evaluated_pairs`, the mean spatial and lineage
#'   statistics, `same_founder_fraction`, `same_spot_fraction`, and one
#'   recent-MRCA fraction per threshold. With `null_summaries` supplied it also
#'   carries every null column, `mrca_age_reduction_vs_null_hours` (null minus
#'   observed), and a `..._enrichment` ratio per threshold that is `Inf` where
#'   the null fraction is zero but the observed one is not.
summarize_visium_lineage_pairs <- function(pair_table,
                                           distance_labels,
                                           candidate_counts,
                                           recent_threshold_hours = c(6, 12),
                                           null_summaries = NULL){
  thresholds <- visium_numeric_vector(
    recent_threshold_hours,
    'recent MRCA thresholds',
    minimum_length = 1L
  )
  summaries <- vector('list', length(distance_labels))
  for(bin_index in seq_along(distance_labels)){
    selected <- which(pair_table$distance_bin == distance_labels[bin_index])
    if(length(selected) == 0L){
      summary_row <- data.frame(
        distance_bin = distance_labels[bin_index],
        candidate_pairs = candidate_counts[bin_index],
        evaluated_pairs = 0L,
        mean_spatial_distance_2d_microns = NA_real_,
        mean_spatial_distance_3d_microns = NA_real_,
        mean_mrca_age_hours = NA_real_,
        median_mrca_age_hours = NA_real_,
        mean_patristic_distance_hours = NA_real_,
        mean_division_distance = NA_real_,
        same_founder_fraction = NA_real_,
        same_spot_fraction = NA_real_,
        stringsAsFactors = FALSE
      )
    } else{
      summary_row <- data.frame(
        distance_bin = distance_labels[bin_index],
        candidate_pairs = candidate_counts[bin_index],
        evaluated_pairs = length(selected),
        mean_spatial_distance_2d_microns = mean(
          pair_table$spatial_distance_2d_microns[selected]
        ),
        mean_spatial_distance_3d_microns = mean(
          pair_table$spatial_distance_3d_microns[selected]
        ),
        mean_mrca_age_hours = mean(pair_table$mrca_age_hours[selected]),
        median_mrca_age_hours = stats::median(pair_table$mrca_age_hours[selected]),
        mean_patristic_distance_hours = mean(
          pair_table$patristic_distance_hours[selected]
        ),
        mean_division_distance = if(all(is.na(
          pair_table$division_distance[selected]
        ))){
          NA_real_
        } else{
          mean(pair_table$division_distance[selected], na.rm = TRUE)
        },
        same_founder_fraction = mean(pair_table$same_founder[selected]),
        same_spot_fraction = mean(pair_table$same_spot[selected]),
        stringsAsFactors = FALSE
      )
    }
    for(threshold in thresholds){
      suffix <- gsub('\\.', '_', format(threshold, trim = TRUE))
      column <- paste0('recent_mrca_le_', suffix, 'h_fraction')
      summary_row[[column]] <- if(length(selected) == 0L){
        NA_real_
      } else{
        mean(pair_table$mrca_age_hours[selected] <= threshold)
      }
    }
    summaries[[bin_index]] <- summary_row
  }
  result <- do.call(rbind, summaries)
  if(!is.null(null_summaries)){
    null_index <- match(result$distance_bin, null_summaries$distance_bin)
    null_columns <- setdiff(names(null_summaries), 'distance_bin')
    for(column in null_columns){
      result[[column]] <- null_summaries[[column]][null_index]
    }
    result$mrca_age_reduction_vs_null_hours <-
      result$null_mean_mrca_age_hours - result$mean_mrca_age_hours
    for(threshold in thresholds){
      suffix <- gsub('\\.', '_', format(threshold, trim = TRUE))
      observed_column <- paste0('recent_mrca_le_', suffix, 'h_fraction')
      null_column <- paste0(observed_column, '_null')
      enrichment_column <- paste0(observed_column, '_enrichment')
      result[[enrichment_column]] <- ifelse(
        result[[null_column]] > 0,
        result[[observed_column]] / result[[null_column]],
        ifelse(result[[observed_column]] > 0, Inf, NA_real_)
      )
    }
  }
  rownames(result) <- NULL
  result
}

#' Collapse cell pairs to spot pairs and summarise by spot separation
#'
#' Cell pairs are grouped by their unordered spot-pair key so each pair of spots
#' contributes one row however many cells it captured; without this, cell-rich
#' spots would dominate the distance summary. Spot pairs are then binned on
#' centre-to-centre separation using edges derived from `spot_pitch`: `[0,1)`
#' microns is the same-spot bin (a spot is exactly 0 from itself) and the
#' remaining edges sit at 1.25, 1.75, 2.25, 3.25, and 4.25 pitches, i.e. between
#' the neighbour shells of the hexagonal array.
#'
#' @param pair_table Lineage-annotated cell pairs; needs `spot_id_a`,
#'   `spot_id_b`, `spot_center_distance_microns`, `mrca_age_hours`,
#'   `patristic_distance_hours`, `division_distance`, and `same_founder`.
#' @param spot_pitch Array pitch in microns setting the bin edges; one positive
#'   finite number.
#' @return A named list with `spot_pairs` (one row per unordered spot pair, with
#'   `contributing_cell_pairs`, the per-pair lineage means, and its
#'   `spot_distance_bin`) and `summary` (one row per bin label in bin order,
#'   averaging over spot pairs rather than cell pairs). An empty `pair_table`
#'   returns the same columns with no spot-pair rows and an all-`NA` summary.
summarize_visium_spot_pair_lineage <- function(pair_table, spot_pitch = 100){
  required <- c(
    'spot_id_a', 'spot_id_b', 'spot_center_distance_microns',
    'mrca_age_hours', 'patristic_distance_hours', 'division_distance',
    'same_founder'
  )
  if(!is.data.frame(pair_table)) stop('pair_table must be a data frame.')
  spot_pitch <- as.numeric(spot_pitch)
  if(length(spot_pitch) != 1L || !is.finite(spot_pitch) || spot_pitch <= 0){
    stop('spot_pitch must be one positive finite number.')
  }
  breaks <- c(
    0,
    1,
    1.25 * spot_pitch,
    1.75 * spot_pitch,
    2.25 * spot_pitch,
    3.25 * spot_pitch,
    4.25 * spot_pitch,
    Inf
  )
  labels <- c(
    'same_spot',
    sprintf('(0,%.0f)', 1.25 * spot_pitch),
    sprintf('[%.0f,%.0f)', 1.25 * spot_pitch, 1.75 * spot_pitch),
    sprintf('[%.0f,%.0f)', 1.75 * spot_pitch, 2.25 * spot_pitch),
    sprintf('[%.0f,%.0f)', 2.25 * spot_pitch, 3.25 * spot_pitch),
    sprintf('[%.0f,%.0f)', 3.25 * spot_pitch, 4.25 * spot_pitch),
    sprintf('[%.0f,Inf)', 4.25 * spot_pitch)
  )
  empty_pairs <- data.frame(
    spot_id_a = character(), spot_id_b = character(), same_spot = logical(),
    spot_center_distance_microns = numeric(), contributing_cell_pairs = integer(),
    mean_mrca_age_hours = numeric(), median_mrca_age_hours = numeric(),
    mean_patristic_distance_hours = numeric(), mean_division_distance = numeric(),
    same_founder_fraction = numeric(), spot_distance_bin = character(),
    stringsAsFactors = FALSE
  )
  empty_summary <- data.frame(
    spot_distance_bin = labels,
    spot_pairs = integer(length(labels)),
    contributing_cell_pairs = integer(length(labels)),
    mean_spot_center_distance_microns = rep(NA_real_, length(labels)),
    mean_mrca_age_hours = rep(NA_real_, length(labels)),
    median_mrca_age_hours = rep(NA_real_, length(labels)),
    mean_patristic_distance_hours = rep(NA_real_, length(labels)),
    mean_division_distance = rep(NA_real_, length(labels)),
    same_founder_fraction = rep(NA_real_, length(labels)),
    stringsAsFactors = FALSE
  )
  if(nrow(pair_table) == 0L){
    return(list(spot_pairs = empty_pairs, summary = empty_summary))
  }
  if(!all(required %in% names(pair_table))){
    stop('pair_table is missing spot or lineage-distance columns.')
  }
  first_spot <- pmin(pair_table$spot_id_a, pair_table$spot_id_b)
  second_spot <- pmax(pair_table$spot_id_a, pair_table$spot_id_b)
  group_key <- paste(first_spot, second_spot, sep = '\r')
  groups <- split(seq_len(nrow(pair_table)), group_key)
  spot_pair_rows <- lapply(groups, function(selected){
    data.frame(
      spot_id_a = first_spot[selected[1]],
      spot_id_b = second_spot[selected[1]],
      same_spot = first_spot[selected[1]] == second_spot[selected[1]],
      spot_center_distance_microns = mean(
        pair_table$spot_center_distance_microns[selected]
      ),
      contributing_cell_pairs = length(selected),
      mean_mrca_age_hours = mean(pair_table$mrca_age_hours[selected]),
      median_mrca_age_hours = stats::median(pair_table$mrca_age_hours[selected]),
      mean_patristic_distance_hours = mean(
        pair_table$patristic_distance_hours[selected]
      ),
      mean_division_distance = if(all(is.na(
        pair_table$division_distance[selected]
      ))){
        NA_real_
      } else{
        mean(pair_table$division_distance[selected], na.rm = TRUE)
      },
      same_founder_fraction = mean(pair_table$same_founder[selected]),
      stringsAsFactors = FALSE
    )
  })
  spot_pairs <- do.call(rbind, spot_pair_rows)
  spot_pairs$spot_distance_bin <- as.character(cut(
    spot_pairs$spot_center_distance_microns,
    breaks = breaks,
    labels = labels,
    include.lowest = TRUE,
    right = FALSE
  ))
  summary_rows <- lapply(labels, function(label){
    selected <- which(spot_pairs$spot_distance_bin == label)
    if(length(selected) == 0L){
      return(data.frame(
        spot_distance_bin = label,
        spot_pairs = 0L,
        contributing_cell_pairs = 0L,
        mean_spot_center_distance_microns = NA_real_,
        mean_mrca_age_hours = NA_real_,
        median_mrca_age_hours = NA_real_,
        mean_patristic_distance_hours = NA_real_,
        mean_division_distance = NA_real_,
        same_founder_fraction = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      spot_distance_bin = label,
      spot_pairs = length(selected),
      contributing_cell_pairs = sum(
        spot_pairs$contributing_cell_pairs[selected]
      ),
      mean_spot_center_distance_microns = mean(
        spot_pairs$spot_center_distance_microns[selected]
      ),
      mean_mrca_age_hours = mean(spot_pairs$mean_mrca_age_hours[selected]),
      median_mrca_age_hours = stats::median(
        spot_pairs$median_mrca_age_hours[selected]
      ),
      mean_patristic_distance_hours = mean(
        spot_pairs$mean_patristic_distance_hours[selected]
      ),
      mean_division_distance = if(all(is.na(
        spot_pairs$mean_division_distance[selected]
      ))){
        NA_real_
      } else{
        mean(spot_pairs$mean_division_distance[selected], na.rm = TRUE)
      },
      same_founder_fraction = mean(
        spot_pairs$same_founder_fraction[selected]
      ),
      stringsAsFactors = FALSE
    )
  })
  list(
    spot_pairs = spot_pairs,
    summary = do.call(rbind, summary_rows)
  )
}

#' Estimate the lineage expectation under randomised cell locations
#'
#' Holds the spatial pairing fixed and permutes which lineage tip sits at which
#' captured-cell position, then re-queries every pair. Averaging over
#' permutations gives the MRCA age and recent-ancestor fractions expected per
#' distance bin when relatedness carries no spatial signal, which is what the
#' observed correlogram is contrasted against.
#'
#' @param pair_table Pairs carrying `cell_index_a`, `cell_index_b`, and
#'   `distance_bin`; the cell indices must address rows of the captured-cell
#'   table in the same order as `captured_node_indices`.
#' @param captured_node_indices Integer lineage-node index per captured cell.
#' @param lineage_index Index from `build_visium_lca_index()`.
#' @param distance_labels Bin labels defining the output rows and their order.
#' @param recent_threshold_hours One or more MRCA-age cutoffs in hours.
#' @param permutations Number of permutations; one non-negative integer, and 0
#'   short-circuits to an all-`NA` table.
#' @param time_units_per_hour Divisor converting lineage times to hours (60 for
#'   PhysiCell minutes).
#' @param seed Integer seed for the permutations.
#' @return A data frame with one row per label: `distance_bin`,
#'   `null_mean_mrca_age_hours`, and one
#'   `recent_mrca_le_<threshold>h_fraction_null` column per threshold. Bins no
#'   permutation ever populated stay `NA` rather than becoming `NaN`.
#' @section Side effects:
#' Calls `set.seed(seed)`, replacing the global RNG state.
permute_visium_lineage_null <- function(pair_table,
                                        captured_node_indices,
                                        lineage_index,
                                        distance_labels,
                                        recent_threshold_hours = c(6, 12),
                                        permutations = 20L,
                                        time_units_per_hour = 60,
                                        seed = 1L){
  permutations <- as.integer(permutations)
  if(length(permutations) != 1L || is.na(permutations) || permutations < 0L){
    stop('permutations must be one non-negative integer.')
  }
  thresholds <- visium_numeric_vector(
    recent_threshold_hours,
    'recent MRCA thresholds',
    minimum_length = 1L
  )
  output <- data.frame(
    distance_bin = distance_labels,
    null_mean_mrca_age_hours = NA_real_,
    stringsAsFactors = FALSE
  )
  for(threshold in thresholds){
    suffix <- gsub('\\.', '_', format(threshold, trim = TRUE))
    output[[paste0('recent_mrca_le_', suffix, 'h_fraction_null')]] <- NA_real_
  }
  if(permutations == 0L || nrow(pair_table) == 0L){
    return(output)
  }
  bin_index <- match(pair_table$distance_bin, distance_labels)
  mean_age <- matrix(NA_real_, nrow = permutations, ncol = length(distance_labels))
  recent <- lapply(thresholds, function(threshold){
    matrix(NA_real_, nrow = permutations, ncol = length(distance_labels))
  })
  set.seed(as.integer(seed))
  for(permutation in seq_len(permutations)){
    permuted <- sample(captured_node_indices, replace = FALSE)
    queried <- query_visium_lineage_pairs(
      lineage_index,
      permuted[pair_table$cell_index_a],
      permuted[pair_table$cell_index_b]
    )
    age_hours <- queried$mrca_age / time_units_per_hour
    for(current_bin in seq_along(distance_labels)){
      selected <- which(bin_index == current_bin)
      if(length(selected) > 0L){
        mean_age[permutation, current_bin] <- mean(age_hours[selected])
        for(threshold_index in seq_along(thresholds)){
          recent[[threshold_index]][permutation, current_bin] <- mean(
            age_hours[selected] <= thresholds[threshold_index]
          )
        }
      }
    }
  }
  output$null_mean_mrca_age_hours <- colMeans(mean_age, na.rm = TRUE)
  output$null_mean_mrca_age_hours[
    colSums(!is.na(mean_age)) == 0L
  ] <- NA_real_
  for(threshold_index in seq_along(thresholds)){
    suffix <- gsub('\\.', '_', format(thresholds[threshold_index], trim = TRUE))
    values <- recent[[threshold_index]]
    column <- paste0('recent_mrca_le_', suffix, 'h_fraction_null')
    output[[column]] <- colMeans(values, na.rm = TRUE)
    output[[column]][colSums(!is.na(values)) == 0L] <- NA_real_
  }
  output
}

#' Run the full pair analysis for one section
#'
#' Maps captured cells onto lineage node indices, samples distance-stratified
#' pairs, attaches MRCA statistics converted to hours, builds the permuted null,
#' and produces both the cell-distance and spot-distance summaries for the
#' slice.
#'
#' @param captured_cells Captured cells for one section; needs the spatial and
#'   spot columns plus `node_id`, and every `node_id` must exist in the lineage
#'   index.
#' @param lineage_index Index from `build_visium_lca_index()`.
#' @param distance_breaks Cell-distance bin edges in microns.
#' @param spot_pitch Array pitch in microns, used for the spot-distance bins.
#' @param recent_threshold_hours Recent-MRCA cutoffs in hours.
#' @param permutations Permutation count for the null.
#' @param max_pairs_per_bin Per-bin cap on retained pairs.
#' @param max_exact_pairs Exhaustive-enumeration limit.
#' @param max_candidate_pairs Candidate draw size above that limit.
#' @param time_units_per_hour Divisor converting lineage times to hours.
#' @param seed Integer seed; pair sampling uses `seed` and the null uses
#'   `seed + 2`.
#' @return A named list with `captured_cells` (the input plus `node_index`),
#'   `pairs` (the annotated pair table), `summary` (the cell-distance
#'   correlogram), `spot_pairs` and `spot_summary` (from
#'   `summarize_visium_spot_pair_lineage()`), and `overall` (a one-row data
#'   frame with the captured-cell and pair counts, `pair_sampling_mode`, the
#'   Spearman correlation of 2-D separation against MRCA age -- `NA` below three
#'   pairs -- the mean cells per occupied spot, and the occupied-spot count).
#' @section Side effects:
#' Reseeds the global RNG by way of the sampling and permutation helpers.
analyze_visium_slice_pairs <- function(captured_cells,
                                       lineage_index,
                                       distance_breaks,
                                       spot_pitch = 100,
                                       recent_threshold_hours = c(6, 12),
                                       permutations = 20L,
                                       max_pairs_per_bin = 100000L,
                                       max_exact_pairs = 2000000L,
                                       max_candidate_pairs = 2000000L,
                                       time_units_per_hour = 60,
                                       seed = 1L){
  if(nrow(captured_cells) > 0L){
    captured_cells$node_index <- unname(
      lineage_index$node_lookup[as.character(captured_cells$node_id)]
    )
    if(anyNA(captured_cells$node_index)){
      stop('One or more captured cells do not map to lineage node IDs.')
    }
  } else{
    captured_cells$node_index <- integer()
  }
  prepared <- prepare_visium_spatial_pairs(
    captured_cells,
    distance_breaks = distance_breaks,
    max_pairs_per_bin = max_pairs_per_bin,
    max_exact_pairs = max_exact_pairs,
    max_candidate_pairs = max_candidate_pairs,
    seed = seed
  )
  pairs <- prepared$pairs
  if(nrow(pairs) > 0L){
    lineage <- query_visium_lineage_pairs(
      lineage_index,
      pairs$node_index_a,
      pairs$node_index_b
    )
    pairs$mrca_node_id <- lineage$mrca_node_id
    pairs$same_founder <- lineage$same_founder
    pairs$mrca_time <- lineage$mrca_time
    pairs$mrca_age_hours <- lineage$mrca_age / time_units_per_hour
    pairs$patristic_distance_hours <-
      lineage$patristic_distance / time_units_per_hour
    pairs$division_distance <- lineage$division_distance
  } else{
    pairs$mrca_node_id <- character()
    pairs$same_founder <- logical()
    pairs$mrca_time <- numeric()
    pairs$mrca_age_hours <- numeric()
    pairs$patristic_distance_hours <- numeric()
    pairs$division_distance <- integer()
  }
  null <- permute_visium_lineage_null(
    pairs,
    captured_cells$node_index,
    lineage_index,
    prepared$distance_labels,
    recent_threshold_hours = recent_threshold_hours,
    permutations = permutations,
    time_units_per_hour = time_units_per_hour,
    seed = as.integer(seed) + 2L
  )
  summary <- summarize_visium_lineage_pairs(
    pairs,
    prepared$distance_labels,
    prepared$candidate_counts,
    recent_threshold_hours = recent_threshold_hours,
    null_summaries = null
  )
  spot_analysis <- summarize_visium_spot_pair_lineage(
    pairs,
    spot_pitch = spot_pitch
  )
  overall <- data.frame(
    captured_cells = nrow(captured_cells),
    possible_cell_pairs = prepared$total_pairs,
    evaluated_cell_pairs = nrow(pairs),
    pair_sampling_mode = prepared$sampling_mode,
    spearman_spatial_vs_mrca = if(nrow(pairs) >= 3L){
      suppressWarnings(stats::cor(
        pairs$spatial_distance_2d_microns,
        pairs$mrca_age_hours,
        method = 'spearman'
      ))
    } else{
      NA_real_
    },
    mean_cells_per_occupied_spot = if(nrow(captured_cells) > 0L){
      mean(table(captured_cells$spot_id))
    } else{
      0
    },
    occupied_spots = length(unique(captured_cells$spot_id)),
    stringsAsFactors = FALSE
  )
  list(
    captured_cells = captured_cells,
    pairs = pairs,
    summary = summary,
    spot_pairs = spot_analysis$spot_pairs,
    spot_summary = spot_analysis$summary,
    overall = overall
  )
}

# ---- Section and correlogram plots ----

#' Draw a section overlay of cells and occupied capture spots
#'
#' Plots the section in `plane_u`/`plane_v` microns at aspect ratio 1, colouring
#' captured cell centres blue and cells lying in the gaps between spots grey,
#' and outlines every occupied spot at its true footprint size. A section with
#' no cells yields a placeholder plot instead of an error.
#'
#' @param slice_cells Section cells with `plane_u`, `plane_v`, and `captured`.
#' @param spots Spot table with `plane_u`, `plane_v`, and `in_tissue`.
#' @param path Destination PNG path.
#' @param spot_diameter Spot footprint diameter in microns, drawn to scale.
#' @return `path`, invisibly.
#' @section Side effects:
#' Writes a 1100 x 1000 pixel PNG at `path`.
write_visium_slice_plot <- function(slice_cells,
                                    spots,
                                    path,
                                    spot_diameter = 55){
  grDevices::png(path, width = 1100, height = 1000, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  relevant_spots <- spots$in_tissue == 1L
  if(nrow(slice_cells) == 0L){
    graphics::plot.new()
    graphics::title('No cell centers intersect this section')
    return(invisible(path))
  }
  x_range <- range(slice_cells$plane_u) + c(-1, 1) * spot_diameter
  y_range <- range(slice_cells$plane_v) + c(-1, 1) * spot_diameter
  graphics::plot(
    slice_cells$plane_u,
    slice_cells$plane_v,
    asp = 1,
    pch = 16,
    cex = 0.45,
    col = ifelse(slice_cells$captured, '#2166AC99', '#BDBDBD55'),
    xlim = x_range,
    ylim = y_range,
    xlab = 'section coordinate u (microns)',
    ylab = 'section coordinate v (microns)',
    main = 'Simulated Visium section: captured and uncaptured cells'
  )
  if(any(relevant_spots)){
    graphics::symbols(
      spots$plane_u[relevant_spots],
      spots$plane_v[relevant_spots],
      circles = rep(spot_diameter / 2, sum(relevant_spots)),
      inches = FALSE,
      add = TRUE,
      fg = '#B2182B',
      bg = NA
    )
  }
  graphics::legend(
    'topright',
    legend = c('captured cell center', 'section cell between spots', 'occupied spot'),
    col = c('#2166AC', '#BDBDBD', '#B2182B'),
    pch = c(16, 16, 1),
    bty = 'n'
  )
  invisible(path)
}

#' Plot mean time to MRCA against in-section distance
#'
#' Draws the observed curve and, when the summary carries the permuted-null
#' column, the null curve on the same axes; the gap between them is the spatial
#' lineage signal. Bins whose means are not finite are skipped, and a summary
#' with no usable bin yields a placeholder plot.
#'
#' @param summary Cell-distance summary from
#'   `summarize_visium_lineage_pairs()`, read for
#'   `mean_spatial_distance_2d_microns`, `mean_mrca_age_hours`, and optionally
#'   `null_mean_mrca_age_hours`.
#' @param path Destination PNG path.
#' @return `path`, invisibly.
#' @section Side effects:
#' Writes a 1100 x 750 pixel PNG at `path`.
write_visium_correlogram_plot <- function(summary, path){
  grDevices::png(path, width = 1100, height = 750, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  valid <- is.finite(summary$mean_spatial_distance_2d_microns) &
    is.finite(summary$mean_mrca_age_hours)
  if(!any(valid)){
    graphics::plot.new()
    graphics::title('Too few captured cell pairs for a correlogram')
    return(invisible(path))
  }
  y_values <- c(
    summary$mean_mrca_age_hours[valid],
    summary$null_mean_mrca_age_hours[valid]
  )
  graphics::plot(
    summary$mean_spatial_distance_2d_microns[valid],
    summary$mean_mrca_age_hours[valid],
    type = 'b',
    pch = 16,
    lwd = 2,
    ylim = range(y_values, finite = TRUE),
    xlab = 'cell-cell distance in section (microns)',
    ylab = 'mean time to MRCA (hours)',
    main = 'Spatial lineage correlogram'
  )
  if(any(is.finite(summary$null_mean_mrca_age_hours[valid]))){
    graphics::lines(
      summary$mean_spatial_distance_2d_microns[valid],
      summary$null_mean_mrca_age_hours[valid],
      type = 'b',
      pch = 1,
      lty = 2,
      lwd = 2,
      col = '#B2182B'
    )
  }
  graphics::legend(
    'bottomright',
    legend = c('observed spatial arrangement', 'permuted lineage locations'),
    col = c('black', '#B2182B'),
    pch = c(16, 1),
    lty = c(1, 2),
    bty = 'n'
  )
  invisible(path)
}

# ---- Run-level and batch-level entry points ----

#' Sample a PhysiCell tumour with simulated Visium sections and analyse them
#'
#' The top-level entry point behind analyze_physicell_visium.R. Loads a
#' completed run, centres the section stack on the live-cell centroid, and for
#' each requested offset cuts a section, lays a capture array over it, assigns
#' cells to spots, and runs the pair analysis. Each slice draws its own array
#' alignment so spot boundaries are not systematically aligned with tumour
#' structure across slices.
#'
#' @param run_dir Completed PhysiCell pipeline run directory.
#' @param output_dir Destination for the analysis; created if absent.
#' @param slice_offsets One or more signed offsets in microns along the plane
#'   normal, measured from the live-cell centroid; one section per value.
#' @param plane_normal Three numbers giving the section normal in simulation
#'   coordinates.
#' @param section_thickness Section thickness in microns.
#' @param spot_diameter Capture-spot diameter in microns.
#' @param spot_pitch Spot centre-to-centre pitch in microns.
#' @param randomize_array_alignment When `TRUE`, each slice draws a rotation in
#'   [0, 60) degrees and a translation spanning one lattice cell; when `FALSE`
#'   the array is unrotated and centred.
#' @param distance_breaks Cell-distance bin edges in microns.
#' @param recent_threshold_hours Recent-MRCA cutoffs in hours; must be
#'   non-negative.
#' @param permutations Permutations per slice for the null.
#' @param max_pairs_per_bin Per-bin cap on retained pairs.
#' @param max_exact_pairs Exhaustive-enumeration limit.
#' @param max_candidate_pairs Candidate draw size above that limit.
#' @param time_units_per_hour Divisor converting lineage times to hours; one
#'   positive finite number (60 for PhysiCell minutes).
#' @param seed Base integer seed; slice `i` uses `seed + 1009 * i`.
#' @param compress_csv When `TRUE`, tables are written gzipped with a `.gz`
#'   suffix.
#' @param show_progress When `TRUE`, one progress message is emitted per slice.
#' @return A named list with `distance_summary`, `spot_distance_summary`, and
#'   `slice_summary` (the per-slice tables stacked and tagged with `slice_id`
#'   and `slice_offset_microns`), `manifest`, and the normalised `output_dir`.
#' @section Side effects:
#' Creates `output_dir` and one `slice_<nnn>/` subdirectory per offset. Each
#' slice directory receives `visium_spots.csv`, `slice_cells.csv`,
#' `spot_cell_membership.csv`, `cell_pair_lineage_sample.csv`,
#' `cell_distance_summary.csv`, `spot_pair_lineage.csv`,
#' `spot_distance_summary.csv`, `slice_summary.csv`, `slice_manifest.csv`,
#' `visium_slice.png`, and `cell_lineage_correlogram.png`. The top level
#' receives the three combined summaries, `analysis_manifest.csv`, and
#' `analysis_complete.txt`. `set.seed()` runs once per slice and again inside
#' the sampling helpers, so the global RNG state is replaced.
#' @note The written `visium_spots.csv` gains `x`, `y`, and `z`: each spot
#'   centre lifted back into simulation coordinates on that section's plane.
#' @note `slice_manifest.csv` records 78 rows of 64 spots and
#'   `analysis_manifest.csv` records 4,992 array spots as literals, which holds
#'   because `make_visium_6_5mm_array()` is called here with its default
#'   geometry.
run_physicell_visium_analysis <- function(
    run_dir,
    output_dir = file.path(run_dir, 'visium_spatial'),
    slice_offsets = 0,
    plane_normal = c(0, 0, 1),
    section_thickness = 5,
    spot_diameter = 55,
    spot_pitch = 100,
    randomize_array_alignment = TRUE,
    distance_breaks = c(0, 20, 40, 60, 80, 100, 150, 200, Inf),
    recent_threshold_hours = c(6, 12),
    permutations = 20L,
    max_pairs_per_bin = 100000L,
    max_exact_pairs = 2000000L,
    max_candidate_pairs = 2000000L,
    time_units_per_hour = 60,
    seed = 1L,
    compress_csv = TRUE,
    show_progress = TRUE){
  slice_offsets <- visium_numeric_vector(
    slice_offsets,
    'slice offsets',
    minimum_length = 1L
  )
  recent_threshold_hours <- visium_numeric_vector(
    recent_threshold_hours,
    'recent MRCA thresholds',
    minimum_length = 1L
  )
  if(any(recent_threshold_hours < 0)){
    stop('recent MRCA thresholds must be non-negative.')
  }
  time_units_per_hour <- as.numeric(time_units_per_hour)
  if(length(time_units_per_hour) != 1L || !is.finite(time_units_per_hour) ||
     time_units_per_hour <= 0){
    stop('time_units_per_hour must be one positive finite number.')
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  inputs <- read_physicell_visium_inputs(run_dir)
  basis <- visium_plane_basis(plane_normal)
  center <- colMeans(inputs$cells[, c('x', 'y', 'z'), drop = FALSE])
  lineage_index <- build_visium_lca_index(inputs$nodes)
  all_summaries <- vector('list', length(slice_offsets))
  all_spot_summaries <- vector('list', length(slice_offsets))
  all_metrics <- vector('list', length(slice_offsets))

  for(slice_index in seq_along(slice_offsets)){
    slice_id <- sprintf('slice_%03d', slice_index)
    slice_seed <- as.integer(seed) + 1009L * slice_index
    slice_dir <- file.path(output_dir, slice_id)
    dir.create(slice_dir, recursive = TRUE, showWarnings = FALSE)
    if(isTRUE(show_progress)){
      message(sprintf(
        '[Visium] %s/%s: offset %.3f microns',
        slice_index,
        length(slice_offsets),
        slice_offsets[slice_index]
      ))
    }
    sliced <- slice_physicell_cells(
      inputs$cells,
      basis,
      center,
      offset = slice_offsets[slice_index],
      thickness = section_thickness
    )
    set.seed(slice_seed)
    rotation <- if(isTRUE(randomize_array_alignment)) runif(1, 0, 60) else 0
    translation <- if(isTRUE(randomize_array_alignment)){
      c(runif(1, -spot_pitch / 2, spot_pitch / 2),
        runif(1, -spot_pitch * sqrt(3) / 4, spot_pitch * sqrt(3) / 4))
    } else{
      c(0, 0)
    }
    spots <- make_visium_6_5mm_array(
      spot_pitch = spot_pitch,
      rotation_degrees = rotation,
      translation = translation,
      slice_id = slice_id
    )
    assigned <- assign_cells_to_visium_spots(
      sliced,
      spots,
      spot_diameter = spot_diameter
    )
    sliced <- assigned$cells
    spots <- assigned$spots
    plane_origin <- center + basis$normal * slice_offsets[slice_index]
    spot_xyz <-
      outer(spots$plane_u, basis$u) +
      outer(spots$plane_v, basis$v) +
      matrix(plane_origin, nrow = nrow(spots), ncol = 3L, byrow = TRUE)
    spots$x <- spot_xyz[, 1]
    spots$y <- spot_xyz[, 2]
    spots$z <- spot_xyz[, 3]
    captured <- sliced[sliced$captured, , drop = FALSE]
    pair_analysis <- analyze_visium_slice_pairs(
      captured,
      lineage_index,
      distance_breaks = distance_breaks,
      spot_pitch = spot_pitch,
      recent_threshold_hours = recent_threshold_hours,
      permutations = permutations,
      max_pairs_per_bin = max_pairs_per_bin,
      max_exact_pairs = max_exact_pairs,
      max_candidate_pairs = max_candidate_pairs,
      time_units_per_hour = time_units_per_hour,
      seed = slice_seed
    )
    captured <- pair_analysis$captured_cells
    sliced$node_index <- unname(
      lineage_index$node_lookup[as.character(sliced$node_id)]
    )
    pair_analysis$summary$slice_id <- slice_id
    pair_analysis$summary$slice_offset_microns <- slice_offsets[slice_index]
    pair_analysis$spot_summary$slice_id <- slice_id
    pair_analysis$spot_summary$slice_offset_microns <- slice_offsets[slice_index]
    pair_analysis$overall$slice_id <- slice_id
    pair_analysis$overall$slice_offset_microns <- slice_offsets[slice_index]
    pair_analysis$overall$section_cells <- nrow(sliced)
    pair_analysis$overall$total_array_spots <- nrow(spots)
    pair_analysis$overall$section_capture_fraction <- if(nrow(sliced) > 0L){
      nrow(captured) / nrow(sliced)
    } else{
      NA_real_
    }

    write_physicell_csv(
      spots,
      file.path(slice_dir, 'visium_spots.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
    write_physicell_csv(
      sliced,
      file.path(slice_dir, 'slice_cells.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
    write_physicell_csv(
      captured,
      file.path(slice_dir, 'spot_cell_membership.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
    write_physicell_csv(
      pair_analysis$pairs,
      file.path(slice_dir, 'cell_pair_lineage_sample.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
    write_physicell_csv(
      pair_analysis$summary,
      file.path(slice_dir, 'cell_distance_summary.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
    write_physicell_csv(
      pair_analysis$spot_pairs,
      file.path(slice_dir, 'spot_pair_lineage.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
    write_physicell_csv(
      pair_analysis$spot_summary,
      file.path(slice_dir, 'spot_distance_summary.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
    write_physicell_csv(
      pair_analysis$overall,
      file.path(slice_dir, 'slice_summary.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
    write_visium_slice_plot(
      sliced,
      spots,
      file.path(slice_dir, 'visium_slice.png'),
      spot_diameter = spot_diameter
    )
    write_visium_correlogram_plot(
      pair_analysis$summary,
      file.path(slice_dir, 'cell_lineage_correlogram.png')
    )

    slice_manifest <- data.frame(
      property = c(
        'slice_id', 'slice_offset_microns', 'section_thickness_microns',
        'plane_normal_x', 'plane_normal_y', 'plane_normal_z',
        'plane_center_x', 'plane_center_y', 'plane_center_z',
        'spot_diameter_microns', 'spot_pitch_microns',
        'array_rows', 'spots_per_row', 'total_array_spots',
        'array_rotation_degrees', 'array_translation_u', 'array_translation_v',
        'section_cells', 'captured_cells', 'occupied_spots',
        'permutations', 'analysis_seed'
      ),
      value = as.character(c(
        slice_id, slice_offsets[slice_index], section_thickness,
        basis$normal, center, spot_diameter, spot_pitch,
        78L, 64L, nrow(spots), rotation, translation,
        nrow(sliced), nrow(captured), sum(spots$in_tissue),
        permutations, slice_seed
      )),
      stringsAsFactors = FALSE
    )
    write_physicell_csv(
      slice_manifest,
      file.path(slice_dir, 'slice_manifest.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
    all_summaries[[slice_index]] <- pair_analysis$summary
    all_spot_summaries[[slice_index]] <- pair_analysis$spot_summary
    all_metrics[[slice_index]] <- pair_analysis$overall
  }

  combined_summary <- do.call(rbind, all_summaries)
  combined_spot_summary <- do.call(rbind, all_spot_summaries)
  combined_metrics <- do.call(rbind, all_metrics)
  write_physicell_csv(
    combined_summary,
    file.path(output_dir, 'cell_distance_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    combined_spot_summary,
    file.path(output_dir, 'spot_distance_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    combined_metrics,
    file.path(output_dir, 'slice_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  analysis_manifest <- data.frame(
    property = c(
      'run_dir', 'number_slices', 'slice_offsets_microns',
      'section_thickness_microns', 'spot_diameter_microns',
      'spot_pitch_microns', 'visium_array_spots', 'distance_breaks_microns',
      'recent_mrca_thresholds_hours', 'permutations',
      'max_pairs_per_bin', 'max_exact_pairs', 'max_candidate_pairs',
      'time_units_per_hour', 'randomize_array_alignment', 'analysis_seed'
    ),
    value = c(
      inputs$run_dir,
      length(slice_offsets),
      paste(slice_offsets, collapse = ','),
      section_thickness,
      spot_diameter,
      spot_pitch,
      4992L,
      paste(distance_breaks, collapse = ','),
      paste(recent_threshold_hours, collapse = ','),
      permutations,
      max_pairs_per_bin,
      max_exact_pairs,
      max_candidate_pairs,
      time_units_per_hour,
      randomize_array_alignment,
      seed
    ),
    stringsAsFactors = FALSE
  )
  write_physicell_csv(
    analysis_manifest,
    file.path(output_dir, 'analysis_manifest.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  writeLines(
    sprintf('Visium spatial lineage analysis completed at %s.', Sys.time()),
    file.path(output_dir, 'analysis_complete.txt')
  )
  list(
    distance_summary = combined_summary,
    spot_distance_summary = combined_spot_summary,
    slice_summary = combined_metrics,
    manifest = analysis_manifest,
    output_dir = normalizePath(output_dir, mustWork = TRUE)
  )
}

#' Pool Visium analyses across replicate tumours
#'
#' Aggregates in two stages: slices are averaged within a replicate tumour
#' first, so a tumour analysed at many offsets does not outweigh one analysed at
#' few, and only then are the replicate means combined into across-tumour mean,
#' SD, and standard error. Both the cell-distance and spot-distance summaries
#' are pooled this way. Replicates are discovered by finding every
#' `visium_spatial/cell_distance_summary.csv` under `batch_dir`, and a
#' replicate is named for the directory that contains its `visium_spatial`.
#'
#' @param batch_dir Directory holding one subdirectory per replicate run; must
#'   exist and contain at least one such summary (plain or `.gz`).
#' @param output_dir Destination for the pooled tables; created if absent.
#' @param compress_csv When `TRUE`, tables are written gzipped.
#' @return A named list with `slice_summary` (every slice row, tagged with
#'   `replicate` and `replicate_dir`), `replicate_summary` (one row per
#'   replicate and distance bin), `aggregate_summary` (one row per distance bin
#'   with `_mean`, `_sd`, and `_se` columns), the matching
#'   `spot_slice_summary`, `spot_replicate_summary`, and
#'   `spot_aggregate_summary`, and the normalised `output_dir`.
#' @section Side effects:
#' Creates `output_dir` and writes `all_slice_distance_summaries.csv`,
#' `replicate_distance_summary.csv`, `aggregate_distance_summary.csv`,
#' `all_slice_spot_distance_summaries.csv`,
#' `replicate_spot_distance_summary.csv`,
#' `aggregate_spot_distance_summary.csv`, `aggregate_lineage_correlogram.png`,
#' and `aggregate_complete.txt`.
#' @note Counts (`candidate_pairs`, `evaluated_pairs`, `spot_pairs`,
#'   `contributing_cell_pairs`) are summed rather than averaged, and
#'   `slice_offset_microns` is excluded from the averaged columns.
aggregate_physicell_visium_results <- function(batch_dir,
                                               output_dir = file.path(
                                                 batch_dir,
                                                 'aggregate'
                                               ),
                                               compress_csv = TRUE){
  batch_dir <- normalizePath(batch_dir, mustWork = TRUE)
  candidate_files <- list.files(
    batch_dir,
    pattern = '^cell_distance_summary\\.csv(\\.gz)?$',
    recursive = TRUE,
    full.names = TRUE
  )
  candidate_files <- candidate_files[
    basename(dirname(candidate_files)) == 'visium_spatial'
  ]
  if(length(candidate_files) == 0L){
    stop('No replicate visium_spatial/cell_distance_summary.csv files were found.')
  }
  tables <- vector('list', length(candidate_files))
  for(index in seq_along(candidate_files)){
    table <- read_physicell_csv(candidate_files[index], check.names = FALSE)
    replicate_dir <- dirname(dirname(candidate_files[index]))
    table$replicate <- basename(replicate_dir)
    table$replicate_dir <- replicate_dir
    tables[[index]] <- table
  }
  slice_table <- do.call(rbind, tables)
  numeric_columns <- names(slice_table)[vapply(slice_table, is.numeric, logical(1))]
  numeric_columns <- setdiff(
    numeric_columns,
    c('slice_offset_microns', 'candidate_pairs', 'evaluated_pairs')
  )
  replicate_keys <- unique(slice_table[, c('replicate', 'distance_bin')])
  replicate_rows <- vector('list', nrow(replicate_keys))
  for(row_index in seq_len(nrow(replicate_keys))){
    selected <- slice_table$replicate == replicate_keys$replicate[row_index] &
      slice_table$distance_bin == replicate_keys$distance_bin[row_index]
    row <- data.frame(
      replicate = replicate_keys$replicate[row_index],
      distance_bin = replicate_keys$distance_bin[row_index],
      slices = sum(selected),
      candidate_pairs = sum(slice_table$candidate_pairs[selected], na.rm = TRUE),
      evaluated_pairs = sum(slice_table$evaluated_pairs[selected], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    for(column in numeric_columns){
      values <- slice_table[[column]][selected]
      row[[column]] <- if(all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
    }
    replicate_rows[[row_index]] <- row
  }
  replicate_summary <- do.call(rbind, replicate_rows)
  distance_bins <- unique(slice_table$distance_bin)
  aggregate_rows <- vector('list', length(distance_bins))
  for(bin_index in seq_along(distance_bins)){
    selected <- replicate_summary$distance_bin == distance_bins[bin_index]
    row <- data.frame(
      distance_bin = distance_bins[bin_index],
      replicates = sum(selected),
      replicates_with_pairs = sum(
        replicate_summary$evaluated_pairs[selected] > 0L,
        na.rm = TRUE
      ),
      slices = sum(replicate_summary$slices[selected]),
      candidate_pairs = sum(
        replicate_summary$candidate_pairs[selected],
        na.rm = TRUE
      ),
      evaluated_pairs = sum(
        replicate_summary$evaluated_pairs[selected],
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
    for(column in numeric_columns){
      values <- replicate_summary[[column]][selected]
      values <- values[is.finite(values)]
      row[[paste0(column, '_mean')]] <- if(length(values) > 0L){
        mean(values)
      } else{
        NA_real_
      }
      row[[paste0(column, '_sd')]] <- if(length(values) > 1L){
        stats::sd(values)
      } else{
        NA_real_
      }
      row[[paste0(column, '_se')]] <- if(length(values) > 1L){
        stats::sd(values) / sqrt(length(values))
      } else{
        NA_real_
      }
    }
    aggregate_rows[[bin_index]] <- row
  }
  aggregate_summary <- do.call(rbind, aggregate_rows)

  spot_tables <- vector('list', length(candidate_files))
  for(index in seq_along(candidate_files)){
    spot_path <- resolve_physicell_csv_path(file.path(
      dirname(candidate_files[index]),
      'spot_distance_summary.csv'
    ))
    spot_table <- read_physicell_csv(spot_path, check.names = FALSE)
    replicate_dir <- dirname(dirname(candidate_files[index]))
    spot_table$replicate <- basename(replicate_dir)
    spot_table$replicate_dir <- replicate_dir
    spot_tables[[index]] <- spot_table
  }
  spot_slice_table <- do.call(rbind, spot_tables)
  spot_numeric_columns <- names(spot_slice_table)[vapply(
    spot_slice_table,
    is.numeric,
    logical(1)
  )]
  spot_numeric_columns <- setdiff(
    spot_numeric_columns,
    c(
      'slice_offset_microns', 'spot_pairs', 'contributing_cell_pairs'
    )
  )
  spot_replicate_keys <- unique(
    spot_slice_table[, c('replicate', 'spot_distance_bin')]
  )
  spot_replicate_rows <- vector('list', nrow(spot_replicate_keys))
  for(row_index in seq_len(nrow(spot_replicate_keys))){
    selected <-
      spot_slice_table$replicate == spot_replicate_keys$replicate[row_index] &
      spot_slice_table$spot_distance_bin ==
        spot_replicate_keys$spot_distance_bin[row_index]
    row <- data.frame(
      replicate = spot_replicate_keys$replicate[row_index],
      spot_distance_bin = spot_replicate_keys$spot_distance_bin[row_index],
      slices = sum(selected),
      spot_pairs = sum(spot_slice_table$spot_pairs[selected], na.rm = TRUE),
      contributing_cell_pairs = sum(
        spot_slice_table$contributing_cell_pairs[selected],
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
    for(column in spot_numeric_columns){
      values <- spot_slice_table[[column]][selected]
      row[[column]] <- if(all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
    }
    spot_replicate_rows[[row_index]] <- row
  }
  spot_replicate_summary <- do.call(rbind, spot_replicate_rows)
  spot_bins <- unique(spot_slice_table$spot_distance_bin)
  spot_aggregate_rows <- vector('list', length(spot_bins))
  for(bin_index in seq_along(spot_bins)){
    selected <- spot_replicate_summary$spot_distance_bin == spot_bins[bin_index]
    row <- data.frame(
      spot_distance_bin = spot_bins[bin_index],
      replicates = sum(selected),
      replicates_with_spot_pairs = sum(
        spot_replicate_summary$spot_pairs[selected] > 0L,
        na.rm = TRUE
      ),
      slices = sum(spot_replicate_summary$slices[selected]),
      spot_pairs = sum(spot_replicate_summary$spot_pairs[selected], na.rm = TRUE),
      contributing_cell_pairs = sum(
        spot_replicate_summary$contributing_cell_pairs[selected],
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
    for(column in spot_numeric_columns){
      values <- spot_replicate_summary[[column]][selected]
      values <- values[is.finite(values)]
      row[[paste0(column, '_mean')]] <- if(length(values) > 0L){
        mean(values)
      } else{
        NA_real_
      }
      row[[paste0(column, '_sd')]] <- if(length(values) > 1L){
        stats::sd(values)
      } else{
        NA_real_
      }
      row[[paste0(column, '_se')]] <- if(length(values) > 1L){
        stats::sd(values) / sqrt(length(values))
      } else{
        NA_real_
      }
    }
    spot_aggregate_rows[[bin_index]] <- row
  }
  spot_aggregate_summary <- do.call(rbind, spot_aggregate_rows)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_physicell_csv(
    slice_table,
    file.path(output_dir, 'all_slice_distance_summaries.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    replicate_summary,
    file.path(output_dir, 'replicate_distance_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    aggregate_summary,
    file.path(output_dir, 'aggregate_distance_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    spot_slice_table,
    file.path(output_dir, 'all_slice_spot_distance_summaries.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    spot_replicate_summary,
    file.path(output_dir, 'replicate_spot_distance_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    spot_aggregate_summary,
    file.path(output_dir, 'aggregate_spot_distance_summary.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  mean_x <- aggregate_summary$mean_spatial_distance_2d_microns_mean
  mean_y <- aggregate_summary$mean_mrca_age_hours_mean
  valid <- is.finite(mean_x) & is.finite(mean_y)
  plot_path <- file.path(output_dir, 'aggregate_lineage_correlogram.png')
  grDevices::png(plot_path, width = 1100, height = 750, res = 140)
  if(any(valid)){
    standard_error <- aggregate_summary$mean_mrca_age_hours_se
    graphics::plot(
      mean_x[valid],
      mean_y[valid],
      type = 'b',
      pch = 16,
      lwd = 2,
      xlab = 'cell-cell distance in section (microns)',
      ylab = 'replicate-mean time to MRCA (hours)',
      main = 'Visium-sampled spatial lineage relationship'
    )
    error_valid <- valid & is.finite(standard_error) & standard_error > 0
    if(any(error_valid)){
      graphics::arrows(
        mean_x[error_valid],
        mean_y[error_valid] - standard_error[error_valid],
        mean_x[error_valid],
        mean_y[error_valid] + standard_error[error_valid],
        angle = 90,
        code = 3,
        length = 0.05
      )
    }
  } else{
    graphics::plot.new()
    graphics::title('No aggregate spatial lineage pairs available')
  }
  grDevices::dev.off()
  writeLines(
    sprintf('Visium replicate aggregation completed at %s.', Sys.time()),
    file.path(output_dir, 'aggregate_complete.txt')
  )
  list(
    slice_summary = slice_table,
    replicate_summary = replicate_summary,
    aggregate_summary = aggregate_summary,
    spot_slice_summary = spot_slice_table,
    spot_replicate_summary = spot_replicate_summary,
    spot_aggregate_summary = spot_aggregate_summary,
    output_dir = normalizePath(output_dir, mustWork = TRUE)
  )
}
