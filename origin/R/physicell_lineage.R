# PhysiCell lineage import and barcode-recording simulation.
#
# This file intentionally contains functions only so it can be sourced from
# tests, notebooks, or the simulate_physicell_lineage.R command-line wrapper.

# ---- Cell IDs and CSV input/output ----

#' Normalize PhysiCell cell IDs to stable character keys
#'
#' Coerces the IDs through `as.numeric` and rejects anything that is not a
#' finite non-negative whole number, then formats them without scientific
#' notation so the same cell always yields the same character key.
#'
#' @export
#' @param values PhysiCell ID values; factors, characters, and numerics are all
#'   accepted.
#' @param field_name Label used in the error message when validation fails.
#' @return A character vector the same length as `values`.
physicell_id <- function(values, field_name = 'cell ID'){
  if(is.factor(values)){
    values <- as.character(values)
  }
  numeric_values <- suppressWarnings(as.numeric(values))
  if(anyNA(numeric_values) || any(!is.finite(numeric_values)) ||
     any(numeric_values < 0) || any(numeric_values %% 1 != 0)){
    stop(sprintf('%s values must be non-negative integer IDs.', field_name))
  }
  format(numeric_values, scientific = FALSE, trim = TRUE, digits = 22)
}

#' Build the canonical CSV output path for a table
#'
#' @export
#' @param path Base CSV path; must be one non-empty string.
#' @param compress When `TRUE`, ensure the path ends in `.gz`; when `FALSE`,
#'   strip a trailing `.gz`.
#' @return The resolved output path as one string.
physicell_csv_output_path <- function(path, compress = TRUE){
  path <- as.character(path)
  if(length(path) != 1 || is.na(path) || !nzchar(path)){
    stop('CSV output path must be one non-empty string.')
  }
  if(isTRUE(compress)){
    if(grepl('\\.gz$', path, ignore.case = TRUE)){
      path
    } else{
      paste0(path, '.gz')
    }
  } else{
    sub('\\.gz$', '', path, ignore.case = TRUE)
  }
}

#' Resolve a CSV path to whichever of its plain/gzipped forms exists
#'
#' Checks both `<path>` and `<path>.gz` (or the `.gz`-stripped form when `path`
#' already ends in `.gz`), preferring the gzipped file.
#'
#' @export
#' @param path CSV or CSV-gzip path; must be one non-empty string.
#' @param required When `TRUE`, a missing file raises an error; when `FALSE`,
#'   `NA_character_` is returned instead.
#' @return The path of the existing file, or `NA_character_`.
resolve_physicell_csv_path <- function(path, required = TRUE){
  path <- as.character(path)
  if(length(path) != 1 || is.na(path) || !nzchar(path)){
    stop('CSV input path must be one non-empty string.')
  }
  candidates <- if(grepl('\\.gz$', path, ignore.case = TRUE)){
    c(path, sub('\\.gz$', '', path, ignore.case = TRUE))
  } else{
    c(paste0(path, '.gz'), path)
  }
  existing <- candidates[file.exists(candidates)]
  if(length(existing) > 0){
    return(existing[1])
  }
  if(isTRUE(required)){
    stop(sprintf(
      'CSV file not found; checked: %s.',
      paste(candidates, collapse = ', ')
    ))
  }
  NA_character_
}

#' Write a table as CSV, gzip-compressed by default
#'
#' @export
#' @param x Object accepted by `utils::write.csv`.
#' @param path Base output path; the `.gz` suffix is added or removed by
#'   `physicell_csv_output_path()`.
#' @param row.names Forwarded to `utils::write.csv`.
#' @param compress When `TRUE`, write through a `gzfile` connection at
#'   compression level 6.
#' @param ... Further arguments forwarded to `utils::write.csv`.
#' @return Invisibly, the normalized path of the file that was written.
#' @section Side effects: Creates or overwrites the resolved CSV file.
write_physicell_csv <- function(x,
                                path,
                                row.names = FALSE,
                                compress = TRUE,
                                ...){
  output_path <- physicell_csv_output_path(path, compress)
  if(isTRUE(compress)){
    connection <- gzfile(
      output_path,
      open = 'wt',
      compression = 6
    )
    on.exit(close(connection), add = TRUE)
    utils::write.csv(
      x,
      connection,
      row.names = row.names,
      ...
    )
    close(connection)
    on.exit(NULL, add = FALSE)
  } else{
    utils::write.csv(
      x,
      output_path,
      row.names = row.names,
      ...
    )
  }
  invisible(normalizePath(output_path, mustWork = TRUE))
}

#' Read a CSV that may be stored plain or gzipped
#'
#' @export
#' @param path CSV or CSV-gzip path, resolved by `resolve_physicell_csv_path()`.
#' @param required When `FALSE` and neither form exists, `NULL` is returned.
#' @param ... Further arguments forwarded to `utils::read.csv`.
#' @return A data frame, or `NULL` for an optional missing table.
read_physicell_csv <- function(path, required = TRUE, ...){
  input_path <- resolve_physicell_csv_path(path, required = required)
  if(is.na(input_path)){
    return(NULL)
  }
  if(grepl('\\.gz$', input_path, ignore.case = TRUE)){
    connection <- gzfile(input_path, open = 'rt')
    on.exit(close(connection), add = TRUE)
    return(utils::read.csv(connection, ...))
  }
  utils::read.csv(input_path, ...)
}

# ---- Division-event ingestion and lineage reconstruction ----

#' Validate and canonicalize a PhysiCell division-event table
#'
#' Strips a byte-order mark and surrounding whitespace from the column names,
#' keeps only `time`, `parent_ID`, and `daughter_ID`, normalizes both ID columns
#' with `physicell_id()`, and sorts events by time with input order as the
#' tie-break. A daughter that already existed as a parent in an earlier row is
#' rejected with one linear `match()` rather than a per-row scan.
#'
#' @export
#' @param divisions Data frame that must contain `time`, `parent_ID`, and
#'   `daughter_ID`. Times must be finite and non-negative, a parent may not be
#'   its own daughter, and each daughter may be born only once.
#' @return The time-sorted event table with an added `.input_order` column and a
#'   `physicell_validated` attribute set to `TRUE`.
validate_physicell_divisions <- function(divisions){
  if(!is.data.frame(divisions)){
    stop('PhysiCell divisions must be supplied as a data frame.')
  }

  names(divisions) <- sub('^\ufeff', '', trimws(names(divisions)))
  required <- c('time', 'parent_ID', 'daughter_ID')
  missing_columns <- setdiff(required, names(divisions))
  if(length(missing_columns) > 0){
    stop(sprintf(
      'PhysiCell division file is missing required column(s): %s.',
      paste(missing_columns, collapse = ', ')
    ))
  }

  divisions <- divisions[, required, drop = FALSE]
  divisions$time <- suppressWarnings(as.numeric(divisions$time))
  if(anyNA(divisions$time) || any(!is.finite(divisions$time)) ||
     any(divisions$time < 0)){
    stop('PhysiCell division times must be finite, non-negative numbers.')
  }
  divisions$parent_ID <- physicell_id(divisions$parent_ID, 'parent_ID')
  divisions$daughter_ID <- physicell_id(divisions$daughter_ID, 'daughter_ID')
  if(any(divisions$parent_ID == divisions$daughter_ID)){
    stop('A PhysiCell division cannot use the same parent_ID and daughter_ID.')
  }
  if(anyDuplicated(divisions$daughter_ID)){
    duplicated_ids <- unique(divisions$daughter_ID[duplicated(divisions$daughter_ID)])
    stop(sprintf(
      'Each PhysiCell daughter_ID must be born once; duplicated ID(s): %s.',
      paste(duplicated_ids, collapse = ', ')
    ))
  }

  divisions$.input_order <- seq_len(nrow(divisions))
  divisions <- divisions[
    order(divisions$time, divisions$.input_order),
    ,
    drop = FALSE
  ]
  rownames(divisions) <- NULL

  if(nrow(divisions) > 0){
    event_indices <- seq_len(nrow(divisions))
    first_parent_indices <- match(
      divisions$daughter_ID,
      divisions$parent_ID
    )
    reused_founder <- which(
      !is.na(first_parent_indices) &
        first_parent_indices < event_indices
    )
    if(length(reused_founder) > 0){
      row_index <- reused_founder[1]
      stop(sprintf(
        'daughter_ID %s already existed before input row %d.',
        divisions$daughter_ID[row_index],
        divisions$.input_order[row_index]
      ))
    }
  }

  attr(divisions, 'physicell_validated') <- TRUE
  divisions
}

#' Read and validate a PhysiCell division-event CSV
#'
#' @export
#' @param path Path to the division CSV; the file must exist.
#' @return The validated table from `validate_physicell_divisions()`.
read_physicell_divisions <- function(path){
  if(length(path) != 1 || is.na(path) || !file.exists(path)){
    stop(sprintf('PhysiCell division file not found: %s.', path))
  }
  divisions <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  validate_physicell_divisions(divisions)
}

#' Read the unique cell IDs out of a PhysiCell cell table
#'
#' @export
#' @param path Path to a PhysiCell cell CSV; the file must exist.
#' @param table_description Label used in the error messages.
#' @param alive_only When `TRUE` and the table has an `alive` column, keep only
#'   rows flagged `true`/`1`; other spellings raise an error.
#' @return A character vector of unique normalized IDs.
read_physicell_cell_ids <- function(path,
                                    table_description = 'PhysiCell cell table',
                                    alive_only = FALSE){
  if(length(path) != 1 || is.na(path) || !file.exists(path)){
    stop(sprintf('%s not found: %s.', table_description, path))
  }
  cell_table <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(cell_table) <- sub('^\ufeff', '', trimws(names(cell_table)))
  if(!('ID' %in% names(cell_table))){
    stop(sprintf('%s must contain an ID column.', table_description))
  }
  if(isTRUE(alive_only) && 'alive' %in% names(cell_table)){
    alive <- cell_table$alive
    if(is.character(alive)){
      normalized <- tolower(trimws(alive))
      if(any(!(normalized %in% c('true', 'false', '1', '0')))){
        stop(sprintf(
          '%s alive values must be true/false or 1/0.',
          table_description
        ))
      }
      alive <- normalized %in% c('true', '1')
    } else{
      alive <- as.logical(alive)
    }
    if(anyNA(alive)){
      stop(sprintf('%s alive values cannot be missing.', table_description))
    }
    cell_table <- cell_table[alive, , drop = FALSE]
  }
  unique(physicell_id(cell_table$ID, 'ID'))
}

#' Read the live-cell IDs that terminate the sampled lineage
#'
#' @export
#' @param path Path to the PhysiCell live-cell CSV.
#' @return Unique IDs of cells whose `alive` flag is true, when that column is
#'   present.
read_physicell_terminal_ids <- function(path){
  read_physicell_cell_ids(
    path,
    table_description = 'PhysiCell live-cell table',
    alive_only = TRUE
  )
}

#' Read the explicit day-zero founder IDs
#'
#' @export
#' @param path Path to the PhysiCell founder CSV.
#' @return Unique IDs from the table, with no alive filtering.
read_physicell_founder_ids <- function(path){
  read_physicell_cell_ids(
    path,
    table_description = 'PhysiCell founder table',
    alive_only = FALSE
  )
}

#' Convert PhysiCell division events into a binary lineage node table
#'
#' PhysiCell reuses a parent's ID across every division, so one ID spans many
#' branch segments. Each event closes the parent's current segment and opens two
#' new ones (`continuing_parent` and `new_daughter`), giving a binary tree. The
#' node columns are preallocated to `founders + 2 * divisions` rows and the
#' currently active node for each PhysiCell ID is held in a hashed environment,
#' so the sweep over events is linear.
#'
#' @export
#' @param divisions Division events; revalidated unless they already carry the
#'   `physicell_validated` attribute.
#' @param end_time Final sampling time; must be finite and no earlier than the
#'   last division. Every still-open segment ends here.
#' @param founder_time Birth time given to founder segments; must be finite and
#'   no later than `end_time`.
#' @param founder_ids Optional explicit founder IDs seeded before the sweep;
#'   parents first seen in the events become founders implicitly.
#' @param show_progress Whether to emit progress messages.
#' @param progress_updates Approximate number of progress messages to emit.
#' @return A data frame with one row per branch segment, holding `node_id`,
#'   `physicell_id`, `parent_node_id`, `birth_time`, `end_time`,
#'   `branch_length`, `division_event`, `origin`, and `is_terminal`.
build_physicell_lineage <- function(divisions,
                                    end_time,
                                    founder_time = 0,
                                    founder_ids = NULL,
                                    show_progress = FALSE,
                                    progress_updates = 20L){
  if(!isTRUE(attr(divisions, 'physicell_validated'))){
    divisions <- validate_physicell_divisions(divisions)
  }
  end_time <- as.numeric(end_time)
  founder_time <- as.numeric(founder_time)
  if(length(end_time) != 1 || !is.finite(end_time) ||
     length(founder_time) != 1 || !is.finite(founder_time) ||
     end_time < founder_time){
    stop('end_time and founder_time must be finite scalars with end_time >= founder_time.')
  }
  if(nrow(divisions) > 0 && max(divisions$time) > end_time){
    stop('end_time precedes at least one PhysiCell division event.')
  }

  if(is.null(founder_ids)){
    founder_ids <- character()
  } else{
    founder_ids <- unique(physicell_id(founder_ids, 'founder_ids'))
  }

  num_divisions <- nrow(divisions)
  possible_founder_count <- length(unique(c(
    founder_ids,
    divisions$parent_ID
  )))
  node_capacity <- possible_founder_count + 2L * num_divisions
  total_work <- num_divisions + 1L
  report_progress <- new_physicell_progress_reporter(
    total = total_work,
    label = 'PhysiCell lineage reconstruction',
    enabled = show_progress,
    updates = progress_updates,
    unit = 'steps'
  )
  report_progress(0L)

  node_id <- if(node_capacity == 0){
    character()
  } else{
    sprintf('node_%06d', seq_len(node_capacity))
  }
  node_physicell_id <- character(node_capacity)
  node_parent_node_id <- rep(NA_character_, node_capacity)
  node_birth_time <- numeric(node_capacity)
  node_end_time <- rep(NA_real_, node_capacity)
  node_branch_length <- numeric(node_capacity)
  node_division_event <- rep(NA_integer_, node_capacity)
  node_origin <- character(node_capacity)
  node_is_terminal <- rep(FALSE, node_capacity)
  node_count <- 0L

  active_nodes <- new.env(
    hash = TRUE,
    parent = emptyenv(),
    size = max(29L, possible_founder_count + num_divisions)
  )

  # Internal: fill the next preallocated branch-segment row
  #
  # @param cell_id PhysiCell ID that owns the segment.
  # @param parent_node_id `node_id` of the segment this one descends from, or
  #   `NA_character_` for a founder.
  # @param birth_time Time at which the segment starts.
  # @param origin One of `founder`, `continuing_parent`, or `new_daughter`.
  # @param event_index Row index of the division that created the segment.
  # @return The integer index of the row just filled.
  add_node <- function(cell_id, parent_node_id, birth_time, origin, event_index){
    node_count <<- node_count + 1L
    if(node_count > node_capacity){
      stop('Internal lineage node capacity was exceeded.')
    }
    node_physicell_id[node_count] <<- cell_id
    node_parent_node_id[node_count] <<- parent_node_id
    node_birth_time[node_count] <<- as.numeric(birth_time)
    node_division_event[node_count] <<- as.integer(event_index)
    node_origin[node_count] <<- origin
    node_count
  }

  # Internal: create and activate one founder segment
  #
  # @param cell_id Cell ID to seed; an already active ID is left untouched.
  # @return Invisibly, the integer node index now active for `cell_id`.
  add_founder <- function(cell_id){
    if(exists(cell_id, envir = active_nodes, inherits = FALSE)){
      return(invisible(get(
        cell_id,
        envir = active_nodes,
        inherits = FALSE
      )))
    }
    new_node_index <- add_node(
      cell_id = cell_id,
      parent_node_id = NA_character_,
      birth_time = founder_time,
      origin = 'founder',
      event_index = NA_integer_
    )
    assign(cell_id, new_node_index, envir = active_nodes)
    invisible(new_node_index)
  }

  if(length(founder_ids) > 0){
    for(founder_id in founder_ids){
      add_founder(founder_id)
    }
  }

  if(num_divisions > 0){
    for(event_index in seq_len(num_divisions)){
      event_time <- divisions$time[event_index]
      parent_id <- divisions$parent_ID[event_index]
      daughter_id <- divisions$daughter_ID[event_index]

      if(!exists(parent_id, envir = active_nodes, inherits = FALSE)){
        add_founder(parent_id)
      }
      if(exists(daughter_id, envir = active_nodes, inherits = FALSE)){
        stop(sprintf(
          'PhysiCell daughter_ID %s already exists at division event %d.',
          daughter_id,
          event_index
        ))
      }

      parent_index <- get(
        parent_id,
        envir = active_nodes,
        inherits = FALSE
      )
      if(event_time < node_birth_time[parent_index]){
        stop(sprintf(
          'Division event %d occurs before parent %s was born.',
          event_index,
          parent_id
        ))
      }
      node_end_time[parent_index] <- event_time
      node_division_event[parent_index] <- event_index

      continuing_parent_index <- add_node(
        cell_id = parent_id,
        parent_node_id = node_id[parent_index],
        birth_time = event_time,
        origin = 'continuing_parent',
        event_index = event_index
      )
      new_daughter_index <- add_node(
        cell_id = daughter_id,
        parent_node_id = node_id[parent_index],
        birth_time = event_time,
        origin = 'new_daughter',
        event_index = event_index
      )
      assign(parent_id, continuing_parent_index, envir = active_nodes)
      assign(daughter_id, new_daughter_index, envir = active_nodes)
      report_progress(event_index)
    }
  }

  if(node_count == 0){
    stop(
      paste(
        'No cells could be inferred. Supply at least one division event or',
        'provide founder_ids (normally from the PhysiCell live-cell table).'
      )
    )
  }

  node_indices <- seq_len(node_count)
  terminal_indices <- as.integer(unlist(
    as.list(active_nodes, all.names = TRUE, sorted = FALSE),
    use.names = FALSE
  ))
  missing_end_indices <- node_indices[is.na(node_end_time[node_indices])]
  node_end_time[missing_end_indices] <- end_time
  node_branch_length[node_indices] <-
    node_end_time[node_indices] - node_birth_time[node_indices]
  node_is_terminal[terminal_indices] <- TRUE
  negative_branch <- which(node_branch_length[node_indices] < 0)
  if(length(negative_branch) > 0){
    stop(sprintf(
      'Node %s has a negative branch length.',
      node_id[negative_branch[1]]
    ))
  }

  report_progress(total_work)
  data.frame(
    node_id = node_id[node_indices],
    physicell_id = node_physicell_id[node_indices],
    parent_node_id = node_parent_node_id[node_indices],
    birth_time = node_birth_time[node_indices],
    end_time = node_end_time[node_indices],
    branch_length = node_branch_length[node_indices],
    division_event = node_division_event[node_indices],
    origin = node_origin[node_indices],
    is_terminal = node_is_terminal[node_indices],
    stringsAsFactors = FALSE
  )
}

# ---- Newick rendering ----

#' Quote node labels that are not Newick-safe
#'
#' @export
#' @param label Vector of arbitrary node or tip labels.
#' @return The labels unchanged where they match `[A-Za-z0-9_.-]+`, otherwise
#'   wrapped in single quotes with any embedded quote doubled.
escape_newick_label <- function(label){
  label <- as.character(label)
  safe <- !is.na(label) & grepl("^[A-Za-z0-9_.-]+$", label)
  escaped <- label
  if(any(!safe)){
    escaped[!safe] <- paste0(
      "'",
      gsub("'", "''", label[!safe], fixed = TRUE),
      "'"
    )
  }
  escaped
}

#' Precompute the topology and tokens needed to render a Newick tree
#'
#' Validates that every parent precedes its children, marks the ancestors of the
#' retained tips in a single reverse pass over the node table, and groups the
#' kept children by parent with a radix sort so each node's children occupy a
#' contiguous slice of `child_order`. Branch lengths and labels are formatted
#' once here instead of during traversal.
#'
#' @export
#' @param nodes Node table with `node_id`, `physicell_id`, `parent_node_id`,
#'   `branch_length`, and `is_terminal`; node IDs must be unique.
#' @param terminal_physicell_ids Optional PhysiCell IDs to retain as tips; each
#'   must belong to a terminal node. `NULL` retains every terminal node.
#' @return A list with `num_nodes`, `num_kept_nodes`, `root_indices`,
#'   `child_counts`, `child_starts`, `child_ends`, `child_order`, `tip_tokens`,
#'   and `close_tokens`.
prepare_physicell_newick <- function(nodes, terminal_physicell_ids = NULL){
  required <- c(
    'node_id', 'physicell_id', 'parent_node_id', 'branch_length', 'is_terminal'
  )
  if(!is.data.frame(nodes) || !all(required %in% names(nodes))){
    stop('nodes is not a valid PhysiCell lineage node table.')
  }
  num_nodes <- nrow(nodes)
  if(num_nodes == 0){
    stop('No requested terminal cells remain in the lineage tree.')
  }
  if(anyDuplicated(nodes$node_id)){
    stop('PhysiCell lineage node IDs must be unique.')
  }

  if(!is.null(terminal_physicell_ids)){
    terminal_physicell_ids <- unique(
      physicell_id(terminal_physicell_ids, 'terminal_physicell_ids')
    )
    unknown_ids <- setdiff(
      terminal_physicell_ids,
      nodes$physicell_id[nodes$is_terminal]
    )
    if(length(unknown_ids) > 0){
      stop(sprintf(
        'Requested terminal PhysiCell ID(s) are absent from the lineage: %s.',
        paste(unknown_ids, collapse = ', ')
      ))
    }
  }

  parent_indices <- match(nodes$parent_node_id, nodes$node_id)
  non_root_indices <- which(!is.na(nodes$parent_node_id))
  missing_parent_indices <- non_root_indices[
    is.na(parent_indices[non_root_indices])
  ]
  if(length(missing_parent_indices) > 0){
    stop(sprintf(
      'Parent node %s is absent from the PhysiCell lineage.',
      nodes$parent_node_id[missing_parent_indices[1]]
    ))
  }
  if(any(parent_indices[non_root_indices] >= non_root_indices)){
    stop(
      'PhysiCell lineage nodes must be ordered with every parent before its children.'
    )
  }
  root_indices <- which(is.na(nodes$parent_node_id))
  if(length(root_indices) == 0){
    stop('No requested terminal cells remain in the lineage tree.')
  }

  selected_terminal <- as.logical(nodes$is_terminal)
  if(!is.null(terminal_physicell_ids)){
    selected_terminal <- selected_terminal &
      nodes$physicell_id %in% terminal_physicell_ids
  }
  keep <- selected_terminal
  if(length(non_root_indices) > 0){
    for(node_index in rev(non_root_indices)){
      if(keep[node_index]){
        keep[parent_indices[node_index]] <- TRUE
      }
    }
  }
  kept_root_indices <- root_indices[keep[root_indices]]
  if(length(kept_root_indices) == 0){
    stop('No requested terminal cells remain in the lineage tree.')
  }

  kept_non_root_indices <- non_root_indices[keep[non_root_indices]]
  kept_child_counts <- tabulate(
    parent_indices[kept_non_root_indices],
    nbins = num_nodes
  )
  if(length(kept_non_root_indices) > 0){
    child_order <- kept_non_root_indices[
      order(
        parent_indices[kept_non_root_indices],
        method = 'radix'
      )
    ]
  } else{
    child_order <- integer()
  }
  child_ends <- cumsum(kept_child_counts)
  child_starts <- child_ends - kept_child_counts + 1L

  branch_lengths <- vapply(
    nodes$branch_length,
    format,
    character(1),
    scientific = FALSE,
    trim = TRUE,
    digits = 15
  )
  tip_tokens <- paste0(
    escape_newick_label(paste0('cell_', nodes$physicell_id)),
    ':',
    branch_lengths
  )
  close_tokens <- paste0(
    ')',
    escape_newick_label(nodes$node_id),
    ':',
    branch_lengths
  )

  list(
    num_nodes = num_nodes,
    num_kept_nodes = sum(keep),
    root_indices = kept_root_indices,
    child_counts = kept_child_counts,
    child_starts = child_starts,
    child_ends = child_ends,
    child_order = child_order,
    tip_tokens = tip_tokens,
    close_tokens = close_tokens
  )
}

#' Traverse the retained forest and emit Newick tokens in order
#'
#' Uses an explicit integer stack rather than recursion: a positive entry is a
#' node to visit, its negation closes that node, and `0` emits a separating
#' comma. The whole forest is wrapped in one synthetic `remote_mito_root`.
#'
#' @export
#' @param prepared Topology from `prepare_physicell_newick()`.
#' @param emit_token Function called once per token; it decides whether tokens
#'   are buffered in memory or streamed to a connection.
#' @param show_progress Whether to emit progress messages.
#' @param progress_updates Approximate number of progress messages to emit.
#' @param progress_label Label used in the progress messages.
#' @return Invisibly, the number of nodes rendered.
emit_physicell_newick <- function(prepared,
                                  emit_token,
                                  show_progress = FALSE,
                                  progress_updates = 20L,
                                  progress_label = 'Newick rendering'){
  if(!is.list(prepared) || !is.function(emit_token)){
    stop('prepared must be a Newick topology and emit_token must be a function.')
  }
  report_progress <- new_physicell_progress_reporter(
    total = prepared$num_kept_nodes,
    label = progress_label,
    enabled = show_progress,
    updates = progress_updates,
    unit = 'nodes'
  )
  report_progress(0L)

  stack <- integer(max(3L * prepared$num_kept_nodes, 1L))
  stack_top <- 0L
  rendered_nodes <- 0L

  # Internal: push one traversal event onto the explicit stack
  #
  # @param event Node index to visit, its negation to close that node, or `0L`
  #   for a separating comma.
  # @return Called for its effect on the enclosing `stack` and `stack_top`.
  push_event <- function(event){
    stack_top <<- stack_top + 1L
    stack[stack_top] <<- event
  }

  emit_token('(')
  for(root_number in seq_along(prepared$root_indices)){
    if(root_number > 1L){
      emit_token(',')
    }
    push_event(prepared$root_indices[root_number])
    while(stack_top > 0L){
      event <- stack[stack_top]
      stack_top <- stack_top - 1L
      if(event == 0L){
        emit_token(',')
        next
      }
      if(event < 0L){
        emit_token(prepared$close_tokens[-event])
        next
      }

      node_index <- event
      rendered_nodes <- rendered_nodes + 1L
      if(prepared$child_counts[node_index] == 0L){
        emit_token(prepared$tip_tokens[node_index])
      } else{
        emit_token('(')
        push_event(-node_index)
        child_start <- prepared$child_starts[node_index]
        child_end <- prepared$child_ends[node_index]
        for(child_position in seq.int(child_end, child_start)){
          push_event(prepared$child_order[child_position])
          if(child_position > child_start){
            push_event(0L)
          }
        }
      }
      report_progress(rendered_nodes)
    }
  }
  emit_token(')remote_mito_root;')
  invisible(rendered_nodes)
}

#' Render a lineage node table as one Newick string
#'
#' Compatibility wrapper that collects the tokens from
#' `emit_physicell_newick()` into a preallocated character vector and pastes
#' them together, preserving unary nodes and elapsed branch times.
#'
#' @export
#' @param nodes Event-resolved lineage node table.
#' @param terminal_physicell_ids Optional PhysiCell IDs to retain as tips.
#' @return One Newick string ending in `)remote_mito_root;`.
#' @note Holds the entire tree in memory; use
#'   `write_physicell_lineage_newick()` for large lineages.
physicell_lineage_to_newick <- function(nodes, terminal_physicell_ids = NULL){
  prepared <- prepare_physicell_newick(nodes, terminal_physicell_ids)
  max_tokens <- 4L * prepared$num_kept_nodes +
    2L * length(prepared$root_indices) + 2L
  tokens <- character(max_tokens)
  token_count <- 0L

  # Internal: append one token to the in-memory token vector
  #
  # @param token One Newick token.
  # @return `NULL`, invisibly.
  emit_token <- function(token){
    token_count <<- token_count + 1L
    tokens[token_count] <<- token
    invisible(NULL)
  }
  emit_physicell_newick(
    prepared,
    emit_token,
    show_progress = FALSE
  )
  paste0(tokens[seq_len(token_count)], collapse = '')
}

#' Stream a lineage node table to a Newick file
#'
#' Renders the tree with the iterative emitter and writes it through a bounded
#' token buffer, so neither the traversal nor the output ever holds a complete
#' subtree string.
#'
#' @export
#' @param nodes Event-resolved lineage node table.
#' @param path Destination file, opened in binary mode.
#' @param terminal_physicell_ids Optional PhysiCell IDs to retain as tips.
#' @param show_progress Whether to emit progress messages.
#' @param progress_updates Approximate number of progress messages to emit.
#' @param buffer_bytes Flush threshold in bytes; must be one positive integer.
#' @param progress_label Label used in the progress messages.
#' @return Invisibly, a list with the normalized `path` and `rendered_nodes`.
#' @section Side effects: Creates or overwrites `path`.
write_physicell_lineage_newick <- function(nodes,
                                           path,
                                           terminal_physicell_ids = NULL,
                                           show_progress = FALSE,
                                           progress_updates = 20L,
                                           buffer_bytes = 1048576L,
                                           progress_label = 'Newick rendering'){
  buffer_bytes <- as.integer(buffer_bytes)
  if(length(buffer_bytes) != 1 || is.na(buffer_bytes) || buffer_bytes < 1){
    stop('buffer_bytes must be one positive integer.')
  }
  prepared <- prepare_physicell_newick(nodes, terminal_physicell_ids)
  connection <- file(path, open = 'wb')
  on.exit(close(connection), add = TRUE)
  token_buffer <- character(65536L)
  token_count <- 0L
  buffered_bytes <- 0L

  # Internal: write the buffered tokens and reset the buffer
  #
  # @return `NULL`, invisibly.
  flush_tokens <- function(){
    if(token_count == 0L){
      return(invisible(NULL))
    }
    writeChar(
      paste0(token_buffer[seq_len(token_count)], collapse = ''),
      connection,
      eos = NULL,
      useBytes = TRUE
    )
    token_count <<- 0L
    buffered_bytes <<- 0L
    invisible(NULL)
  }

  # Internal: buffer one token, flushing once the byte budget is reached
  #
  # @param token One Newick token.
  # @return `NULL`, invisibly.
  emit_token <- function(token){
    token_count <<- token_count + 1L
    token_buffer[token_count] <<- token
    buffered_bytes <<- buffered_bytes + nchar(token, type = 'bytes')
    if(buffered_bytes >= buffer_bytes ||
       token_count >= length(token_buffer)){
      flush_tokens()
    }
    invisible(NULL)
  }

  rendered_nodes <- emit_physicell_newick(
    prepared,
    emit_token,
    show_progress = show_progress,
    progress_updates = progress_updates,
    progress_label = progress_label
  )
  emit_token('\n')
  flush_tokens()
  close(connection)
  on.exit(NULL, add = FALSE)
  invisible(list(
    path = normalizePath(path, mustWork = TRUE),
    rendered_nodes = rendered_nodes
  ))
}

# ---- Barcode target layout and reference construction ----

#' Split an integer total across weights by largest remainder
#'
#' @export
#' @param total Non-negative integer to allocate.
#' @param weights Non-negative finite weights with a positive sum; they are
#'   normalized internally.
#' @return Integer counts, one per weight, summing exactly to `total`.
largest_remainder_counts <- function(total, weights){
  total <- as.integer(total)
  weights <- as.numeric(weights)
  if(length(total) != 1 || is.na(total) || total < 0 ||
     length(weights) == 0 || any(!is.finite(weights)) ||
     any(weights < 0) || sum(weights) <= 0){
    stop('Invalid total or weights for largest-remainder allocation.')
  }
  raw_counts <- total * weights / sum(weights)
  counts <- floor(raw_counts)
  remainder <- total - sum(counts)
  if(remainder > 0){
    recipients <- order(raw_counts - counts, decreasing = TRUE)[seq_len(remainder)]
    counts[recipients] <- counts[recipients] + 1L
  }
  as.integer(counts)
}

#' Place recorder targets along the barcode and assign edit-rate classes
#'
#' Three layouts are read from `target_spec$config`: `U` spreads the targets
#' uniformly, `R` samples positions without replacement, and `S:first:gap`
#' places them at a fixed spacing. Class counts come from
#' `largest_remainder_counts()` on the high/medium/low fractions and are then
#' attached to a shuffled copy of the positions.
#'
#' @export
#' @param target_spec One JSON target block with `num_targets`, `config`, and
#'   `edit_rate_class_fractions`; `NULL` or a missing `num_targets` yields no
#'   targets.
#' @param barcode_length Barcode length in bases; `num_targets` may not exceed
#'   it and a spaced layout may not run past it.
#' @return A character vector of `High`/`Medium`/`Low` classes named by the
#'   one-based target position.
physicell_target_positions <- function(target_spec, barcode_length){
  if(is.null(target_spec) || is.null(target_spec$num_targets)){
    return(setNames(character(), character()))
  }
  num_targets <- as.integer(unlist(target_spec$num_targets, use.names = FALSE))
  if(length(num_targets) != 1 || is.na(num_targets) || num_targets < 0 ||
     num_targets > barcode_length){
    stop('Target num_targets must be one integer between zero and bc_length.')
  }
  if(num_targets == 0){
    return(setNames(character(), character()))
  }

  config_text <- toupper(gsub('[[:space:]]+', '', as.character(target_spec$config)))
  config_parts <- strsplit(config_text, ':', fixed = TRUE)[[1]]
  config <- substr(config_parts[1], 1, 1)
  if(config == 'U'){
    positions <- unique(round(seq(1, barcode_length, length.out = num_targets)))
  } else if(config == 'R'){
    positions <- sample(seq_len(barcode_length), num_targets, replace = FALSE)
  } else if(config == 'S'){
    if(length(config_parts) != 3){
      stop("Spaced target configuration must have the form 'S:first:gap'.")
    }
    first <- suppressWarnings(as.integer(config_parts[2]))
    gap <- suppressWarnings(as.integer(config_parts[3]))
    if(is.na(first) || is.na(gap) || first < 1 || gap < 0){
      stop('Spaced target first position and gap must be valid integers.')
    }
    positions <- first + (seq_len(num_targets) - 1L) * (gap + 1L)
    if(max(positions) > barcode_length){
      stop('Spaced target configuration extends beyond bc_length.')
    }
  } else{
    stop("Target config must begin with 'U', 'R', or 'S'.")
  }
  if(length(positions) != num_targets){
    stop('Target configuration did not produce the requested number of unique positions.')
  }

  fractions <- target_spec$edit_rate_class_fractions
  weights <- c(fractions$high, fractions$medium, fractions$low)
  class_counts <- largest_remainder_counts(num_targets, weights)
  shuffled <- sample(positions, length(positions), replace = FALSE)
  classes <- rep(c('High', 'Medium', 'Low'), class_counts)
  setNames(classes, as.character(shuffled))
}

#' Expand each target into its editing window
#'
#' Every position within `window_spec$size` bases of a target joins that
#' target's window. With a decaying window the class drops one step inside the
#' inner half of the window and two steps beyond it, and positions that decay
#' below `Low` are excluded. Where windows overlap, the highest class wins.
#'
#' @export
#' @param target_classes Named `High`/`Medium`/`Low` classes from
#'   `physicell_target_positions()`.
#' @param window_spec JSON editing-window block supplying `size` and `decaying`.
#' @param barcode_sequence Barcode reference as a character vector of bases.
#' @param same_base_only When `TRUE` (base editing), only window positions
#'   carrying the same base as the target are eligible.
#' @return A list with `classes`, the expanded named class vector, and
#'   `windows`, the accepted positions per `be_window_<i>`/`nuc_window_<i>`.
expand_physicell_target_windows <- function(target_classes,
                                            window_spec,
                                            barcode_sequence,
                                            same_base_only){
  if(length(target_classes) == 0){
    return(list(classes = target_classes, windows = list()))
  }
  window_size <- as.integer(window_spec$size)
  if(length(window_size) != 1 || is.na(window_size) || window_size < 0){
    stop('Editing-window size must be one non-negative integer.')
  }
  decaying <- isTRUE(window_spec$decaying)
  class_order <- c(Low = 1L, Medium = 2L, High = 3L)
  expanded <- target_classes
  windows <- list()

  # Internal: decay an edit-rate class by its distance from the target
  #
  # @param edit_class `High`, `Medium`, or `Low`.
  # @param distance Absolute distance in bases from the target position.
  # @return The class unchanged for a non-decaying window, the decayed class
  #   otherwise, or `NA_character_` once it would fall below `Low`.
  lower_class <- function(edit_class, distance){
    class_index <- class_order[[edit_class]]
    if(!decaying || window_size == 0){
      return(edit_class)
    }
    degrees <- if(distance / window_size < 0.5) 1L else 2L
    new_index <- class_index - degrees
    if(new_index < 1L){
      return(NA_character_)
    }
    names(class_order)[match(new_index, class_order)]
  }

  target_positions <- as.integer(names(target_classes))
  for(target_index in seq_along(target_positions)){
    target_position <- target_positions[target_index]
    candidate_positions <- seq(
      max(1L, target_position - window_size),
      min(length(barcode_sequence), target_position + window_size)
    )
    if(same_base_only){
      candidate_positions <- candidate_positions[
        barcode_sequence[candidate_positions] == barcode_sequence[target_position]
      ]
    }

    accepted <- integer()
    for(position in candidate_positions){
      edit_class <- if(position == target_position){
        target_classes[[as.character(target_position)]]
      } else{
        lower_class(
          target_classes[[as.character(target_position)]],
          abs(position - target_position)
        )
      }
      if(is.na(edit_class)){
        next
      }
      accepted <- c(accepted, position)
      old_class <- unname(expanded[as.character(position)])
      if(length(old_class) == 0 || is.na(old_class) ||
         class_order[[edit_class]] > class_order[[old_class]]){
        expanded[as.character(position)] <- edit_class
      }
    }
    windows[[paste0(if(same_base_only) 'be' else 'nuc', '_window_', target_index)]] <-
      unique(accepted)
  }

  list(classes = expanded, windows = windows)
}

#' Obtain the barcode reference sequence
#'
#' Reads `params$barcode_sequence` when one is configured (plain text or FASTA,
#' resolved against `params_dir` unless the path is absolute). Otherwise it
#' composes a random sequence whose base counts follow `bc_nuc_composition`,
#' forcing the base-editor source base onto every BE target position and
#' borrowing from the most abundant other base when the composition cannot
#' supply enough of it.
#'
#' @export
#' @param params Parsed JSON parameter list; `bc_length` fixes the length.
#' @param be_positions Integer positions that must carry `be_from`.
#' @param be_from Base-editor source base.
#' @param params_dir Directory used to resolve a relative sequence path.
#' @return A character vector of `A`/`G`/`C`/`T` of length `bc_length`.
physicell_barcode_reference <- function(params,
                                        be_positions,
                                        be_from,
                                        params_dir = '.'){
  barcode_length <- as.integer(params$bc_length)
  supplied_sequence <- params$barcode_sequence
  if(!is.null(supplied_sequence)){
    sequence_path <- as.character(supplied_sequence)
    if(!grepl('^/', sequence_path)){
      sequence_path <- file.path(params_dir, sequence_path)
    }
    if(!file.exists(sequence_path)){
      stop(sprintf('Configured barcode_sequence not found: %s.', sequence_path))
    }
    lines <- readLines(sequence_path, warn = FALSE)
    lines <- lines[!grepl('^>', lines)]
    sequence_text <- toupper(gsub('[[:space:]]+', '', paste(lines, collapse = '')))
    sequence <- strsplit(sequence_text, '', fixed = TRUE)[[1]]
    if(length(sequence) != barcode_length ||
       any(!(sequence %in% c('A', 'G', 'C', 'T')))){
      stop('Configured barcode_sequence must match bc_length and contain only A/G/C/T.')
    }
    return(sequence)
  }

  composition <- params$bc_nuc_composition
  weights <- c(
    composition$frac_a,
    composition$frac_g,
    composition$frac_c,
    composition$frac_t
  )
  base_names <- c('A', 'G', 'C', 'T')
  base_counts <- largest_remainder_counts(barcode_length, weights)
  names(base_counts) <- base_names
  be_count <- length(be_positions)
  target_index <- match(be_from, base_names)
  deficit <- max(be_count - base_counts[target_index], 0)
  if(deficit > 0){
    for(unused in seq_len(deficit)){
      donors <- base_counts
      donors[target_index] <- -Inf
      donor <- which.max(donors)
      if(base_counts[donor] <= 0){
        stop('Unable to reconcile BE targets with barcode nucleotide composition.')
      }
      base_counts[donor] <- base_counts[donor] - 1L
      base_counts[target_index] <- base_counts[target_index] + 1L
    }
  }
  base_counts[target_index] <- base_counts[target_index] - be_count

  sequence <- character(barcode_length)
  sequence[be_positions] <- be_from
  other_positions <- setdiff(seq_len(barcode_length), be_positions)
  sequence[other_positions] <- sample(
    rep(base_names, base_counts),
    length(other_positions),
    replace = FALSE
  )
  sequence
}

# ---- Substitution models and per-position rate sets ----

#' Coerce substitution-model parameters to a fixed-length numeric vector
#'
#' @export
#' @param values Numeric vector, list, or one semicolon-delimited string.
#' @param expected_length Number of parameters the model consumes.
#' @return A numeric vector of exactly `expected_length` entries; absent values
#'   stay `NA_real_` so the model can derive them.
physicell_model_parameters <- function(values, expected_length){
  if(is.character(values) && length(values) == 1 && grepl(';', values, fixed = TRUE)){
    values <- strsplit(gsub('[[:space:]]+', '', values), ';', fixed = TRUE)[[1]]
  }
  if(is.list(values)){
    values <- vapply(values, function(value){
      if(is.null(value) || length(value) == 0){
        return(NA_real_)
      }
      as.numeric(value)[1]
    }, numeric(1))
  } else{
    values <- as.numeric(values)
  }
  if(length(values) < expected_length){
    values <- c(values, rep(NA_real_, expected_length - length(values)))
  }
  values[seq_len(expected_length)]
}

#' Build the barcode substitution-probability matrix for one model
#'
#' Supports `JC`, `K80`, `K81`, `F81`, `HKY`, and `GTR`. Equilibrium base
#' frequencies are taken from the supplied sequence; `F81` and `HKY` rescale
#' their off-diagonal entries so the mean off-diagonal rate matches the
#' configured baseline.
#'
#' @export
#' @param model_name Model label, matched case-insensitively.
#' @param model_parameters Parameters in the layout
#'   `physicell_model_parameters()` expects for that model.
#' @param sequence Barcode reference used for the base frequencies.
#' @return A 4-by-4 finite non-negative matrix with `A`, `G`, `C`, `T` dimnames
#'   and a zero diagonal.
physicell_substitution_matrix <- function(model_name, model_parameters, sequence){
  bases <- c('A', 'G', 'C', 'T')
  fractions <- as.numeric(table(factor(sequence, levels = bases))) / length(sequence)
  model_name <- toupper(as.character(model_name))

  # Internal: complete a ratio/transition/transversion triple
  #
  # @param parameters Three values ordered ratio, transition, transversion; at
  #   least two must be supplied.
  # @param model_label Model name used in the error messages.
  # @return A named vector with all three values filled in.
  derive_three <- function(parameters, model_label){
    ratio <- parameters[1]
    transition <- parameters[2]
    transversion <- parameters[3]
    supplied <- !is.na(c(ratio, transition, transversion))
    if(sum(supplied) < 2){
      stop(sprintf('%s requires at least two rate parameters.', model_label))
    }
    if(is.na(ratio)){
      ratio <- transition / transversion
    }
    if(is.na(transition)){
      transition <- ratio * transversion
    }
    if(is.na(transversion)){
      if(ratio == 0){
        stop(sprintf('%s cannot derive a transversion rate from ratio zero.', model_label))
      }
      transversion <- transition / ratio
    }
    c(ratio = ratio, transition = transition, transversion = transversion)
  }

  if(model_name == 'JC'){
    parameters <- physicell_model_parameters(model_parameters, 1)
    result <- matrix(parameters[1], 4, 4)
    diag(result) <- 0
  } else if(model_name == 'K80'){
    rates <- derive_three(physicell_model_parameters(model_parameters, 3), 'K80')
    result <- rbind(
      c(0, rates['transition'], rates['transversion'], rates['transversion']),
      c(rates['transition'], 0, rates['transversion'], rates['transversion']),
      c(rates['transversion'], rates['transversion'], 0, rates['transition']),
      c(rates['transversion'], rates['transversion'], rates['transition'], 0)
    )
  } else if(model_name == 'K81'){
    parameters <- physicell_model_parameters(model_parameters, 3)
    result <- rbind(
      c(0, parameters[1], parameters[2], parameters[3]),
      c(parameters[1], 0, parameters[3], parameters[2]),
      c(parameters[2], parameters[3], 0, parameters[1]),
      c(parameters[3], parameters[2], parameters[1], 0)
    )
  } else if(model_name == 'F81'){
    parameters <- physicell_model_parameters(model_parameters, 1)
    baseline <- parameters[1]
    if(baseline == 0){
      result <- matrix(0, 4, 4)
    } else{
      destination_rates <- baseline * (1 + fractions - 0.25)
      result <- matrix(rep(destination_rates, 4), nrow = 4, byrow = TRUE)
      diag(result) <- 0
      result <- result * (baseline / (sum(result) / 12))
    }
  } else if(model_name == 'HKY'){
    rates <- derive_three(physicell_model_parameters(model_parameters, 3), 'HKY')
    destination_weights <- 1 + fractions - 0.25
    result <- rbind(
      destination_weights * c(0, rates['transition'], rates['transversion'], rates['transversion']),
      destination_weights * c(rates['transition'], 0, rates['transversion'], rates['transversion']),
      destination_weights * c(rates['transversion'], rates['transversion'], 0, rates['transition']),
      destination_weights * c(rates['transversion'], rates['transversion'], rates['transition'], 0)
    )
    if(any(result != 0)){
      harmonic_mean <- 2 / (
        1 / rates['transition'] + 1 / rates['transversion']
      )
      result <- result * (harmonic_mean / (sum(result) / 12))
    }
  } else if(model_name == 'GTR'){
    parameters <- physicell_model_parameters(model_parameters, 6)
    result <- rbind(
      c(0, parameters[1] * fractions[2], parameters[2] * fractions[3], parameters[3] * fractions[4]),
      c(parameters[1] * fractions[1], 0, parameters[4] * fractions[3], parameters[5] * fractions[4]),
      c(parameters[2] * fractions[1], parameters[4] * fractions[2], 0, parameters[6] * fractions[4]),
      c(parameters[3] * fractions[1], parameters[5] * fractions[2], parameters[6] * fractions[3], 0)
    )
  } else{
    stop(sprintf('Unsupported barcode substitution model: %s.', model_name))
  }

  if(any(!is.finite(result)) || any(result < 0)){
    stop(sprintf('%s produced invalid substitution probabilities.', model_name))
  }
  dimnames(result) <- list(bases, bases)
  result
}

#' Draw class-stratified per-target mutation probabilities
#'
#' Samples 1000 gamma variates with shape 0.5 and the requested mean, cuts them
#' at the quantiles implied by the Low/Medium/High class proportions, and draws
#' each target's rate from its own class bin, so higher classes receive the
#' heavier tail. Falls back to the mean for every target when those quantiles
#' are not distinct.
#'
#' @export
#' @param target_classes Named `High`/`Medium`/`Low` classes; an empty vector
#'   returns an empty result.
#' @param mean_probability Mean per-division probability; must be one finite
#'   non-negative value, and zero returns all-zero rates.
#' @return A named numeric vector of probabilities capped just below one.
draw_physicell_target_rates <- function(target_classes, mean_probability){
  if(length(target_classes) == 0){
    return(setNames(numeric(), character()))
  }
  mean_probability <- as.numeric(mean_probability)
  if(length(mean_probability) != 1 || !is.finite(mean_probability) ||
     mean_probability < 0){
    stop('Target mutation probability must be one finite non-negative value.')
  }
  if(mean_probability == 0){
    return(setNames(rep(0, length(target_classes)), names(target_classes)))
  }

  bootstrap <- stats::rgamma(
    1000,
    shape = 0.5,
    scale = mean_probability / 0.5
  )
  class_counts <- table(factor(
    target_classes,
    levels = c('Low', 'Medium', 'High')
  ))
  nonempty <- names(class_counts)[class_counts > 0]
  breakpoints <- c(0, cumsum(class_counts[nonempty]) / length(target_classes))
  cutpoints <- unique(stats::quantile(
    bootstrap,
    probs = breakpoints,
    names = FALSE
  ))
  if(length(cutpoints) != length(nonempty) + 1L){
    return(setNames(rep(mean_probability, length(target_classes)), names(target_classes)))
  }
  bins <- cut(
    bootstrap,
    breaks = cutpoints,
    labels = nonempty,
    include.lowest = TRUE
  )

  rates <- vapply(target_classes, function(edit_class){
    candidates <- bootstrap[bins == edit_class]
    if(length(candidates) == 0){
      mean_probability
    } else{
      sample(candidates, 1)
    }
  }, numeric(1))
  pmin(rates, 1 - .Machine$double.eps)
}

#' Build the per-position event rates for one cell type and editing state
#'
#' Every non-target position gets the background substitution probabilities of
#' its reference base plus the background insertion and deletion rates; target
#' positions start empty and instead receive their drawn base-editing
#' substitution rate and nuclease insertion/deletion rates. A random share of
#' the remaining positions is then made invariant.
#'
#' @export
#' @details Event names are `sub_<1-4>` for a substitution to `A`, `G`, `C`, or
#'   `T`, plus `insertion` and `deletion`. The per-division probabilities are
#'   converted to continuous-time hazards once here, and the totals, cumulative
#'   event-selection probabilities, and single-event shortcuts are cached so
#'   branch simulation never recomputes them.
#'
#' @param params Parsed JSON parameter list.
#' @param cell_type Key into `params$cell_type_dict$cell_type_params`; the type
#'   must carry a positive finite `cell_cycle_length`.
#' @param editing_state Either `uninduced_editing_params` or
#'   `induced_editing_params`.
#' @param barcode_sequence Barcode reference as a character vector.
#' @param be_classes Expanded base-editing target classes.
#' @param nuc_classes Expanded nuclease target classes.
#' @param be_to Base-editor destination base.
#' @return A list with `probabilities`, `event_hazards`, `total_hazards`,
#'   `event_cumulative_probabilities`, `single_event_names`,
#'   `active_positions`, `active_position_names`, `cell_cycle_length`,
#'   `substitution_matrix`, and `invariant_positions`.
physicell_rate_set <- function(params,
                               cell_type,
                               editing_state,
                               barcode_sequence,
                               be_classes,
                               nuc_classes,
                               be_to){
  cell_params <- params$cell_type_dict$cell_type_params[[cell_type]]
  if(is.null(cell_params)){
    stop(sprintf('Cell type %s is absent from the parameter JSON.', cell_type))
  }
  state_params <- cell_params[[editing_state]]
  if(is.null(state_params)){
    stop(sprintf(
      'Editing state %s is absent for cell type %s.',
      editing_state,
      cell_type
    ))
  }
  cell_cycle_length <- as.numeric(cell_params$cell_cycle_length)
  if(length(cell_cycle_length) != 1 || !is.finite(cell_cycle_length) ||
     cell_cycle_length <= 0){
    stop('Selected cell type must have a positive finite cell_cycle_length.')
  }

  substitution_matrix <- physicell_substitution_matrix(
    state_params$bc_substitution_model,
    state_params$bc_sub_model_params,
    barcode_sequence
  )
  background_insertion <- as.numeric(state_params$bc_bg_insertion_prob_per_division)
  background_deletion <- as.numeric(state_params$bc_bg_deletion_prob_per_division)
  probabilities <- list()
  bases <- c('A', 'G', 'C', 'T')
  joint_target_positions <- union(
    as.integer(names(be_classes)),
    as.integer(names(nuc_classes))
  )

  for(position in seq_along(barcode_sequence)){
    from_index <- match(barcode_sequence[position], bases)
    position_probabilities <- substitution_matrix[from_index, ]
    names(position_probabilities) <- paste0('sub_', seq_along(bases))
    position_probabilities <- position_probabilities[position_probabilities > 0]
    position_probabilities <- c(
      position_probabilities,
      insertion = background_insertion,
      deletion = background_deletion
    )
    if(position %in% joint_target_positions){
      position_probabilities <- numeric()
    }
    probabilities[[position]] <- position_probabilities
  }

  be_rates <- draw_physicell_target_rates(
    be_classes,
    state_params$be_mutations_per_target_per_division
  )
  be_destination <- match(be_to, bases)
  for(position_name in names(be_rates)){
    position <- as.integer(position_name)
    probabilities[[position]] <- c(
      probabilities[[position]],
      setNames(be_rates[[position_name]], paste0('sub_', be_destination))
    )
  }

  insertion_rates <- draw_physicell_target_rates(
    nuc_classes,
    state_params$nuc_insertions_per_target_per_division
  )
  deletion_rates <- draw_physicell_target_rates(
    nuc_classes,
    state_params$nuc_deletions_per_target_per_division
  )
  for(position_name in union(names(insertion_rates), names(deletion_rates))){
    position <- as.integer(position_name)
    probabilities[[position]] <- c(
      probabilities[[position]],
      insertion = insertion_rates[[position_name]],
      deletion = deletion_rates[[position_name]]
    )
  }

  invariant_fraction <- as.numeric(cell_params$bc_invariant_sites)
  eligible_invariants <- setdiff(seq_along(barcode_sequence), joint_target_positions)
  num_invariants <- round(invariant_fraction * length(eligible_invariants))
  invariant_positions <- if(num_invariants > 0){
    sample(eligible_invariants, num_invariants, replace = FALSE)
  } else{
    integer()
  }
  probabilities[invariant_positions] <- replicate(
    length(invariant_positions),
    numeric(),
    simplify = FALSE
  )

  event_hazards <- lapply(
    probabilities,
    physicell_probability_hazard,
    cell_cycle_length = cell_cycle_length
  )
  total_hazards <- vapply(event_hazards, sum, numeric(1))
  active_positions <- which(total_hazards > 0)
  event_cumulative_probabilities <- lapply(
    event_hazards,
    function(hazards){
      total_hazard <- sum(hazards)
      if(total_hazard <= 0){
        return(numeric())
      }
      cumulative <- cumsum(hazards / total_hazard)
      cumulative[length(cumulative)] <- 1
      cumulative
    }
  )
  single_event_names <- vapply(
    event_hazards,
    function(hazards){
      if(length(hazards) == 1){
        names(hazards)[1]
      } else{
        NA_character_
      }
    },
    character(1)
  )

  list(
    probabilities = probabilities,
    event_hazards = event_hazards,
    total_hazards = total_hazards,
    event_cumulative_probabilities = event_cumulative_probabilities,
    single_event_names = single_event_names,
    active_positions = active_positions,
    active_position_names = as.character(active_positions),
    cell_cycle_length = cell_cycle_length,
    substitution_matrix = substitution_matrix,
    invariant_positions = invariant_positions
  )
}

# ---- PALINCODE recorder model ----

#' Validate a PALINCODE probability as one value per cBit
#'
#' @export
#' @param value Scalar probability, or one probability per cBit.
#' @param number Number of cBits.
#' @param name Configuration key used in the error message.
#' @return A numeric vector of length `number` with every entry in `[0, 1]`.
palincode_probability_vector <- function(value, number, name){
  value <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
  if(length(value) == 1L){
    value <- rep(value, number)
  }
  if(length(value) != number || any(!is.finite(value)) ||
     any(value < 0) || any(value > 1)){
    stop(sprintf(
      '%s must contain one probability or one per PALINCODE cBit.',
      name
    ))
  }
  value
}

#' Validate one configured count as an integer
#'
#' @export
#' @param value Candidate count; must be a single finite whole number no larger
#'   than `.Machine$integer.max`.
#' @param name Configuration key used in the error message.
#' @param allow_zero When `TRUE` zero is accepted, otherwise the value must be
#'   positive.
#' @return The value as an integer.
palincode_positive_integer <- function(value, name, allow_zero = FALSE){
  numeric_value <- suppressWarnings(as.numeric(value))
  minimum <- if(allow_zero) 0L else 1L
  if(length(numeric_value) != 1L || !is.finite(numeric_value) ||
     numeric_value < minimum || numeric_value > .Machine$integer.max ||
     numeric_value %% 1 != 0){
    stop(sprintf(
      '%s must be one %s integer.',
      name,
      if(allow_zero) 'non-negative' else 'positive'
    ))
  }
  as.integer(numeric_value)
}

#' Draw unique static nucleotide identifiers for integrations
#'
#' Resamples until no two identifiers collide, after checking that
#' `4^identifier_length` can accommodate `number` distinct strings.
#'
#' @export
#' @param number Number of identifiers to draw; zero returns an empty vector.
#' @param identifier_length Length in bases of each identifier.
#' @return A character vector of `number` unique A/C/G/T strings.
palincode_static_ids <- function(number, identifier_length){
  if(number == 0L){
    return(character())
  }
  if(log(number) > identifier_length * log(4) + 1e-12){
    stop(paste(
      'palincode_adapter.static_id_length is too short for unique',
      'integrations.'
    ))
  }
  repeat{
    identifiers <- vapply(seq_len(number), function(index){
      paste0(
        sample(c('A', 'C', 'G', 'T'), identifier_length, TRUE),
        collapse = ''
      )
    }, character(1))
    if(!anyDuplicated(identifiers)){
      return(identifiers)
    }
  }
}

#' Prepare a PALINCODE recording model
#'
#' Each integration carries `num_cbits_per_integration` cBits occupying barcode
#' positions `1..num_cbits`, with any founder-label sites appended after them. A
#' cBit edits at most once and resolves to `left`, `right`, or `both` in the
#' configured proportions, encoded as `1`, `2`, and `3` against wild type `0`.
#' The per-cell-cycle probabilities are converted to continuous-time hazards for
#' the induced and uninduced states.
#'
#' @export
#' @param params Parsed remote_mito parameter list; `palincode_adapter` holds
#'   the recorder configuration and `physicell_adapter` the storage policy.
#' @param cell_type Cell type whose `cell_cycle_length` scales the hazards;
#'   defaults to the configured founder cell type.
#' @param num_integrations Integrations per cell; defaults to the largest
#'   `max_bc_ints_per_cell` option.
#' @param founder_label_sites Extra positions appended for founder labels.
#' @param seed Optional seed set before the static IDs are drawn.
#' @return A prepared model list carrying `is_palincode = TRUE`, the barcode
#'   geometry, `rate_sets` for both editing states, and a `palincode` block with
#'   the cBit names, outcome fractions, static IDs, and state encoding.
#' @section Side effects: Calls `set.seed()` when `seed` is supplied.
prepare_palincode_recording_model <- function(params,
                                               cell_type = NULL,
                                               num_integrations = NULL,
                                               founder_label_sites = 0,
                                               seed = NULL){
  configuration <- params$palincode_adapter
  if(is.null(configuration)){
    configuration <- list()
  }
  if(!is.list(configuration)){
    stop('palincode_adapter must be a JSON object.')
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.null(cell_type)){
    cell_type <- as.character(params$cell_type_dict$founder_cell_type)
  }
  if(length(cell_type) != 1L || is.na(cell_type) || !nzchar(cell_type) ||
     is.null(params$cell_type_dict$cell_type_params[[cell_type]])){
    stop('PALINCODE requires one configured cell type.')
  }
  if(is.null(num_integrations)){
    integration_options <- as.integer(unlist(
      params$max_bc_ints_per_cell,
      use.names = FALSE
    ))
    num_integrations <- max(integration_options)
  }
  num_integrations <- palincode_positive_integer(
    num_integrations,
    'num_integrations'
  )
  founder_label_sites <- palincode_positive_integer(
    founder_label_sites,
    'founder_label_sites',
    allow_zero = TRUE
  )
  num_cbits <- palincode_positive_integer(
    if(is.null(configuration$num_cbits_per_integration)){
      2L
    } else{
      configuration$num_cbits_per_integration
    },
    'palincode_adapter.num_cbits_per_integration'
  )
  cbit_names <- configuration$cbit_names
  if(is.null(cbit_names)){
    width <- max(3L, nchar(num_cbits))
    cbit_names <- paste0(
      'cBit_',
      formatC(seq_len(num_cbits), width = width, flag = '0')
    )
  } else{
    cbit_names <- as.character(unlist(cbit_names, use.names = FALSE))
  }
  if(length(cbit_names) != num_cbits || anyNA(cbit_names) ||
     any(!nzchar(cbit_names)) || anyDuplicated(cbit_names)){
    stop('palincode_adapter.cbit_names must uniquely name every cBit.')
  }

  induced_source <- configuration$induced_edit_probability_per_cbit_per_cell_cycle
  if(is.null(induced_source)){
    induced_source <- configuration$edit_probability_per_cbit_per_cell_cycle
  }
  if(is.null(induced_source)){
    induced_source <- 0.1
  }
  uninduced_source <-
    configuration$uninduced_edit_probability_per_cbit_per_cell_cycle
  if(is.null(uninduced_source)){
    uninduced_source <- 0
  }
  induced_probabilities <- palincode_probability_vector(
    induced_source,
    num_cbits,
    'palincode_adapter.induced_edit_probability_per_cbit_per_cell_cycle'
  )
  uninduced_probabilities <- palincode_probability_vector(
    uninduced_source,
    num_cbits,
    'palincode_adapter.uninduced_edit_probability_per_cbit_per_cell_cycle'
  )
  left_fractions <- palincode_probability_vector(
    if(is.null(configuration$left_edit_fraction)) 0.495 else
      configuration$left_edit_fraction,
    num_cbits,
    'palincode_adapter.left_edit_fraction'
  )
  right_fractions <- palincode_probability_vector(
    if(is.null(configuration$right_edit_fraction)) 0.495 else
      configuration$right_edit_fraction,
    num_cbits,
    'palincode_adapter.right_edit_fraction'
  )
  both_fractions <- palincode_probability_vector(
    if(is.null(configuration$both_edit_fraction)) 0.01 else
      configuration$both_edit_fraction,
    num_cbits,
    'palincode_adapter.both_edit_fraction'
  )
  outcome_fractions <- cbind(
    left = left_fractions,
    right = right_fractions,
    both = both_fractions
  )
  if(any(abs(rowSums(outcome_fractions) - 1) > 1e-8)){
    stop(paste(
      'PALINCODE left_edit_fraction, right_edit_fraction, and',
      'both_edit_fraction must sum to one for every cBit.'
    ))
  }
  outcome_fractions <- outcome_fractions / rowSums(outcome_fractions)

  cell_cycle_length <- as.numeric(
    params$cell_type_dict$cell_type_params[[cell_type]]$cell_cycle_length
  )
  if(length(cell_cycle_length) != 1L || is.na(cell_cycle_length) ||
     cell_cycle_length <= 0){
    stop('Every PALINCODE cell type needs a positive cell_cycle_length.')
  }
  barcode_length <- num_cbits + founder_label_sites
  cbit_positions <- seq_len(num_cbits)
  founder_label_positions <- if(founder_label_sites == 0L){
    integer()
  } else{
    num_cbits + seq_len(founder_label_sites)
  }

  # Internal: turn per-cBit probabilities into one PALINCODE rate set
  #
  # @param probabilities Per-cell-cycle edit probability for each cBit.
  # @return A list with `total_hazards`, `active_positions`,
  #   `active_position_names`, `outcome_fractions`, and `cell_cycle_length`.
  make_rate_set <- function(probabilities){
    hazards <- physicell_probability_hazard(
      probabilities,
      cell_cycle_length
    )
    names(hazards) <- as.character(cbit_positions)
    list(
      total_hazards = hazards,
      active_positions = cbit_positions[hazards > 0],
      active_position_names = as.character(cbit_positions[hazards > 0]),
      outcome_fractions = outcome_fractions,
      cell_cycle_length = cell_cycle_length
    )
  }

  editing_induction_time <- suppressWarnings(as.numeric(
    params$editing_induction$timepoint
  ))
  if(length(editing_induction_time) != 1L ||
     !is.finite(editing_induction_time)){
    editing_induction_time <- Inf
  }
  adapter <- params$physicell_adapter
  if(is.null(adapter)){
    adapter <- list()
  }
  profile_storage <- if(is.null(adapter$profile_storage)){
    'sparse'
  } else{
    tolower(as.character(adapter$profile_storage))
  }
  if(length(profile_storage) != 1L ||
     !(profile_storage %in% c('dense', 'sparse'))){
    stop('physicell_adapter.profile_storage must be dense or sparse.')
  }
  compact_output <- if(is.null(adapter$compact_output)){
    TRUE
  } else{
    isTRUE(adapter$compact_output)
  }
  retain_internal_profiles <- if(is.null(adapter$retain_internal_profiles)){
    TRUE
  } else{
    isTRUE(adapter$retain_internal_profiles)
  }
  static_id_length <- palincode_positive_integer(
    if(is.null(configuration$static_id_length)) 12L else
      configuration$static_id_length,
    'palincode_adapter.static_id_length'
  )
  integration_static_ids <- palincode_static_ids(
    num_integrations,
    static_id_length
  )
  target_classes <- setNames(rep('PALINCODE', num_cbits), cbit_positions)

  list(
    is_palincode = TRUE,
    barcode_sequence = rep('A', barcode_length),
    barcode_length = barcode_length,
    num_integrations = num_integrations,
    founder_label_sites = founder_label_sites,
    founder_label_positions = founder_label_positions,
    cell_type = cell_type,
    be_from = 'A',
    be_to = 'G',
    be_targets = target_classes,
    be_target_classes = target_classes,
    nuc_target_classes = setNames(character(), character()),
    windows = list(),
    position_windows = setNames(
      replicate(barcode_length, character(), simplify = FALSE),
      as.character(seq_len(barcode_length))
    ),
    close_be_window = FALSE,
    close_nuc_window = FALSE,
    editing_induction_time = editing_induction_time,
    recorder_system = 'PALINCODE',
    profile_storage = profile_storage,
    compact_output = compact_output,
    retain_internal_profiles = retain_internal_profiles,
    output_positions = seq_len(barcode_length),
    rate_sets = list(
      uninduced_editing_params = make_rate_set(uninduced_probabilities),
      induced_editing_params = make_rate_set(induced_probabilities)
    ),
    palincode = list(
      num_cbits_per_integration = num_cbits,
      cbit_positions = cbit_positions,
      cbit_names = cbit_names,
      induced_edit_probabilities = induced_probabilities,
      uninduced_edit_probabilities = uninduced_probabilities,
      outcome_fractions = outcome_fractions,
      static_id_length = static_id_length,
      integration_static_ids = integration_static_ids,
      state_encoding = c(wild_type = 0L, left = 1L, right = 2L, both = 3L)
    )
  )
}

# ---- Prime-editing recorder model ----

#' Prepare a known-pegRNA prime-editing recording model
#'
#' Lays out the targets with `physicell_target_positions()` on the nuclease
#' target block, hands them to `prepare_prime_editing_backend()` for pegRNA
#' assignment, and scales each target's base per-cell-cycle probability by its
#' pegRNA editing efficiency before converting to hazards. Every target edits at
#' most once, to the locked state `1`.
#'
#' @export
#' @param params Parsed remote_mito parameter list; `prime_editing_backend`
#'   supplies the pool and the optional probability overrides.
#' @param cell_type Cell type whose `cell_cycle_length` scales the hazards;
#'   defaults to the configured founder cell type.
#' @param num_integrations Integrations per cell; defaults to the largest
#'   `max_bc_ints_per_cell` option.
#' @param founder_label_sites Extra non-target positions reserved for founder
#'   labels, taken from the lowest free positions.
#' @param params_dir Directory used to resolve a relative pegRNA-pool path.
#' @param seed Optional seed set before pool assignment and static-ID drawing.
#' @return A prepared model list carrying `is_prime_editing = TRUE`, the barcode
#'   geometry, `rate_sets` for both editing states, and a `prime_editing` block
#'   holding the backend plus the base and effective probabilities.
#' @section Side effects: Calls `set.seed()` when `seed` is supplied.
#' @note Requires `prime_editing.R` to have been sourced first.
prepare_prime_editing_recording_model <- function(params,
                                                   cell_type = NULL,
                                                   num_integrations = NULL,
                                                   founder_label_sites = 0,
                                                   params_dir = '.',
                                                   seed = NULL){
  if(!exists('prepare_prime_editing_backend', mode = 'function')){
    stop('Source prime_editing.R before preparing a prime-editing model.')
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.null(cell_type)){
    cell_type <- as.character(params$cell_type_dict$founder_cell_type)
  }
  if(length(cell_type) != 1L || is.na(cell_type) || !nzchar(cell_type) ||
     is.null(params$cell_type_dict$cell_type_params[[cell_type]])){
    stop('Prime editing requires one configured cell type.')
  }
  if(is.null(num_integrations)){
    integration_options <- as.integer(unlist(
      params$max_bc_ints_per_cell,
      use.names = FALSE
    ))
    num_integrations <- max(integration_options)
  }
  num_integrations <- palincode_positive_integer(
    num_integrations,
    'num_integrations'
  )
  founder_label_sites <- palincode_positive_integer(
    founder_label_sites,
    'founder_label_sites',
    allow_zero = TRUE
  )
  barcode_length <- suppressWarnings(as.integer(params$bc_length))
  if(length(barcode_length) != 1L || is.na(barcode_length) ||
     barcode_length < 1L){
    stop('bc_length must be one positive integer.')
  }
  nuc_targets <- physicell_target_positions(
    params$nuclease_targets,
    barcode_length
  )
  target_positions <- sort(as.integer(names(nuc_targets)))
  if(length(target_positions) == 0L){
    stop('Prime editing requires at least one configured nuclease target.')
  }
  backend <- prepare_prime_editing_backend(
    params,
    target_positions = target_positions,
    params_dir = params_dir,
    seed = seed
  )
  barcode_sequence <- physicell_barcode_reference(
    params,
    integer(),
    'A',
    params_dir
  )
  founder_label_positions <- head(
    setdiff(seq_len(barcode_length), target_positions),
    founder_label_sites
  )
  if(length(founder_label_positions) != founder_label_sites){
    stop(sprintf(
      'Only %d non-prime-editing positions are available for founder labels.',
      length(setdiff(seq_len(barcode_length), target_positions))
    ))
  }

  configuration <- params$prime_editing_backend
  if(is.null(configuration)){
    configuration <- list()
  }
  type_parameters <- params$cell_type_dict$cell_type_params[[cell_type]]

  # Internal: resolve one editing state's base per-cell-cycle probability
  #
  # @param state Either `induced_editing_params` or `uninduced_editing_params`.
  # @param fallback Value used when neither the backend override nor the cell
  #   type supplies one.
  # @return A validated probability vector with one entry per target position.
  base_probability <- function(state, fallback){
    override_name <- paste0(
      if(state == 'induced_editing_params') 'induced' else 'uninduced',
      '_edit_probability_per_cell_cycle'
    )
    value <- configuration[[override_name]]
    if(is.null(value) && !is.null(type_parameters[[state]])){
      value <- type_parameters[[state]]$nuc_insertions_per_target_per_division
    }
    if(is.null(value)){
      value <- fallback
    }
    prime_editing_probability_vector(
      value,
      length(target_positions),
      paste0('prime_editing_backend.', override_name)
    )
  }
  induced_base <- base_probability('induced_editing_params', 0.1)
  uninduced_base <- base_probability('uninduced_editing_params', 0)
  induced_probabilities <- prime_editing_scale_probability(
    induced_base,
    backend$targets$editing_efficiency
  )
  uninduced_probabilities <- prime_editing_scale_probability(
    uninduced_base,
    backend$targets$editing_efficiency
  )
  cell_cycle_length <- suppressWarnings(as.numeric(
    type_parameters$cell_cycle_length
  ))
  if(length(cell_cycle_length) != 1L || is.na(cell_cycle_length) ||
     cell_cycle_length <= 0){
    stop('Every prime-editing cell type needs a positive cell_cycle_length.')
  }

  # Internal: turn per-target probabilities into one prime-editing rate set
  #
  # @param probabilities Effective per-cell-cycle probability per target.
  # @return A list with `total_hazards`, `active_positions`,
  #   `active_position_names`, and `cell_cycle_length`.
  make_rate_set <- function(probabilities){
    hazards <- physicell_probability_hazard(
      probabilities,
      cell_cycle_length
    )
    names(hazards) <- as.character(target_positions)
    list(
      total_hazards = hazards,
      active_positions = target_positions[hazards > 0],
      active_position_names = as.character(target_positions[hazards > 0]),
      cell_cycle_length = cell_cycle_length
    )
  }
  editing_induction_time <- suppressWarnings(as.numeric(
    params$editing_induction$timepoint
  ))
  if(length(editing_induction_time) != 1L ||
     !is.finite(editing_induction_time)){
    editing_induction_time <- Inf
  }
  adapter <- params$physicell_adapter
  if(is.null(adapter)){
    adapter <- list()
  }
  profile_storage <- if(is.null(adapter$profile_storage)){
    'sparse'
  } else{
    tolower(as.character(adapter$profile_storage))
  }
  if(length(profile_storage) != 1L ||
     !(profile_storage %in% c('dense', 'sparse'))){
    stop('physicell_adapter.profile_storage must be dense or sparse.')
  }
  compact_output <- if(is.null(adapter$compact_output)){
    TRUE
  } else{
    isTRUE(adapter$compact_output)
  }
  retain_internal_profiles <- if(is.null(adapter$retain_internal_profiles)){
    TRUE
  } else{
    isTRUE(adapter$retain_internal_profiles)
  }
  static_id_length <- palincode_positive_integer(
    if(is.null(configuration$static_id_length)) 12L else
      configuration$static_id_length,
    'prime_editing_backend.static_id_length'
  )
  backend$integration_static_ids <- palincode_static_ids(
    num_integrations,
    static_id_length
  )
  backend$static_id_length <- static_id_length
  backend$induced_base_probabilities <- induced_base
  backend$uninduced_base_probabilities <- uninduced_base
  backend$induced_effective_probabilities <- induced_probabilities
  backend$uninduced_effective_probabilities <- uninduced_probabilities
  target_classes <- setNames(
    rep('prime_editing', length(target_positions)),
    target_positions
  )

  list(
    is_prime_editing = TRUE,
    barcode_sequence = barcode_sequence,
    barcode_length = barcode_length,
    num_integrations = num_integrations,
    founder_label_sites = founder_label_sites,
    founder_label_positions = founder_label_positions,
    cell_type = cell_type,
    be_from = 'A',
    be_to = 'G',
    be_targets = setNames(character(), character()),
    be_target_classes = setNames(character(), character()),
    nuc_target_classes = target_classes,
    windows = list(),
    position_windows = setNames(
      replicate(barcode_length, character(), simplify = FALSE),
      as.character(seq_len(barcode_length))
    ),
    close_be_window = FALSE,
    close_nuc_window = TRUE,
    editing_induction_time = editing_induction_time,
    recorder_system = if(is.null(adapter$recorder_system)){
      'prime editing'
    } else{
      as.character(adapter$recorder_system)
    },
    profile_storage = profile_storage,
    compact_output = compact_output,
    retain_internal_profiles = retain_internal_profiles,
    output_positions = sort(unique(c(
      target_positions,
      founder_label_positions
    ))),
    rate_sets = list(
      uninduced_editing_params = make_rate_set(uninduced_probabilities),
      induced_editing_params = make_rate_set(induced_probabilities)
    ),
    prime_editing = backend
  )
}

# ---- Recording-model dispatch and barcode profile storage ----

#' Prepare the barcode recording model for a PhysiCell replay
#'
#' Dispatches to the PALINCODE or prime-editing preparer when either recorder is
#' configured, and otherwise builds the generic base-editor/nuclease barcode
#' model: it lays out the BE and nuclease targets, loads or generates the
#' reference, expands the editing windows, reserves founder-label positions
#' among the non-recording sites, and builds both rate sets.
#'
#' @export
#' @param params Parsed remote_mito JSON parameter list.
#' @param cell_type Cell type supplying the rate parameters; defaults to the
#'   configured founder cell type.
#' @param num_integrations Integrations per cell; defaults to the largest
#'   `max_bc_ints_per_cell` option.
#' @param founder_label_sites Number of non-recording positions reserved to
#'   encode the founder index; must fit within the barcode.
#' @param params_dir Directory used to resolve relative paths in `params`.
#' @param seed Optional seed set before target layout and reference generation.
#' @return A prepared model list with the barcode reference and length, target
#'   classes, editing windows and their closure flags, induction time, recorder
#'   identity, storage and output policy, `output_positions`, and `rate_sets`.
#' @section Side effects: Calls `set.seed()` when `seed` is supplied.
#' @note PALINCODE and prime editing cannot both be enabled. A configured
#'   `nuclease_targets$prime_editing_system` that does not activate the unified
#'   backend only warns: the generic path then emits ordinary one-base insertion
#'   alleles instead of replaying guide sequences.
prepare_physicell_recording_model <- function(params,
                                              cell_type = NULL,
                                              num_integrations = NULL,
                                              founder_label_sites = 0,
                                              params_dir = '.',
                                              seed = NULL){
  if(!is.list(params)){
    stop('params must be the parsed remote_mito JSON parameter list.')
  }
  recorder_system <- params$physicell_adapter$recorder_system
  palincode_configuration <- params$palincode_adapter
  if(!is.null(palincode_configuration) && !is.list(palincode_configuration)){
    stop('palincode_adapter must be a JSON object.')
  }
  palincode_enabled <-
    (!is.null(palincode_configuration) &&
       !identical(palincode_configuration$enabled, FALSE)) ||
    (!is.null(recorder_system) &&
       grepl('PALINCODE', as.character(recorder_system), ignore.case = TRUE))
  prime_configuration_present <-
    !is.null(params$prime_editing_backend) ||
    isTRUE(params$nuclease_targets$prime_editing_system)
  if(prime_configuration_present &&
     !exists('prime_editing_enabled', mode = 'function')){
    stop('Source prime_editing.R before preparing a prime-editing model.')
  }
  prime_enabled <- prime_configuration_present && prime_editing_enabled(params)
  if(isTRUE(palincode_enabled) && isTRUE(prime_enabled)){
    stop('PALINCODE and prime editing cannot share one barcode model.')
  }
  if(isTRUE(prime_enabled)){
    return(prepare_prime_editing_recording_model(
      params,
      cell_type = cell_type,
      num_integrations = num_integrations,
      founder_label_sites = founder_label_sites,
      params_dir = params_dir,
      seed = seed
    ))
  }
  if(isTRUE(palincode_enabled)){
    return(prepare_palincode_recording_model(
      params,
      cell_type = cell_type,
      num_integrations = num_integrations,
      founder_label_sites = founder_label_sites,
      seed = seed
    ))
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  barcode_length <- as.integer(params$bc_length)
  if(length(barcode_length) != 1 || is.na(barcode_length) ||
     barcode_length < 1){
    stop('bc_length must be one positive integer.')
  }
  if(is.null(cell_type)){
    cell_type <- as.character(params$cell_type_dict$founder_cell_type)
  }
  if(is.null(num_integrations)){
    integration_options <- as.integer(unlist(
      params$max_bc_ints_per_cell,
      use.names = FALSE
    ))
    num_integrations <- max(integration_options)
  }
  num_integrations <- as.integer(num_integrations)
  if(length(num_integrations) != 1 || is.na(num_integrations) ||
     num_integrations < 1){
    stop('num_integrations must be one positive integer.')
  }
  founder_label_sites <- as.integer(founder_label_sites)
  if(length(founder_label_sites) != 1 || is.na(founder_label_sites) ||
     founder_label_sites < 0 || founder_label_sites > barcode_length){
    stop(sprintf(
      'founder_label_sites must be between 0 and the barcode length (%d).',
      barcode_length
    ))
  }

  be_targets <- physicell_target_positions(params$be_targets, barcode_length)
  nuc_targets <- physicell_target_positions(params$nuclease_targets, barcode_length)
  if(length(be_targets) == 0 && is.null(params$be_conversion_pattern)){
    be_from <- 'A'
    be_to <- 'G'
  } else{
    conversion <- toupper(gsub('[[:space:]]+', '', params$be_conversion_pattern))
    conversion_match <- regexec('^([ACGT])-+>([ACGT])$', conversion)
    conversion_parts <- regmatches(conversion, conversion_match)[[1]]
    if(length(conversion_parts) != 3 || conversion_parts[2] == conversion_parts[3]){
      stop("be_conversion_pattern must look like 'A --> G' with distinct A/C/G/T bases.")
    }
    be_from <- conversion_parts[2]
    be_to <- conversion_parts[3]
  }

  barcode_sequence <- physicell_barcode_reference(
    params,
    as.integer(names(be_targets)),
    be_from,
    params_dir
  )

  be_expanded <- expand_physicell_target_windows(
    be_targets,
    params$be_targets$editing_window,
    barcode_sequence,
    same_base_only = TRUE
  )
  nuc_expanded <- expand_physicell_target_windows(
    nuc_targets,
    params$nuclease_targets$editing_window,
    barcode_sequence,
    same_base_only = FALSE
  )
  recording_positions <- unique(c(
    as.integer(names(be_expanded$classes)),
    as.integer(names(nuc_expanded$classes))
  ))
  founder_label_positions <- head(
    setdiff(seq_len(barcode_length), recording_positions),
    founder_label_sites
  )
  if(length(founder_label_positions) != founder_label_sites){
    stop(sprintf(
      paste(
        'Only %d non-recording barcode positions are available for',
        '%d requested founder-label sites.'
      ),
      length(setdiff(seq_len(barcode_length), recording_positions)),
      founder_label_sites
    ))
  }

  if(isTRUE(params$nuclease_targets$prime_editing_system)){
    warning(
      paste(
        'The PhysiCell adapter currently emits ordinary one-base insertion',
        'alleles; prime-editing guide sequences are not replayed.'
      ),
      call. = FALSE
    )
  }

  rate_sets <- list(
    uninduced_editing_params = physicell_rate_set(
      params,
      cell_type,
      'uninduced_editing_params',
      barcode_sequence,
      be_expanded$classes,
      nuc_expanded$classes,
      be_to
    ),
    induced_editing_params = physicell_rate_set(
      params,
      cell_type,
      'induced_editing_params',
      barcode_sequence,
      be_expanded$classes,
      nuc_expanded$classes,
      be_to
    )
  )

  all_windows <- c(be_expanded$windows, nuc_expanded$windows)
  position_windows <- setNames(
    replicate(barcode_length, character(), simplify = FALSE),
    as.character(seq_len(barcode_length))
  )
  for(window_name in names(all_windows)){
    for(position in all_windows[[window_name]]){
      position_windows[[as.character(position)]] <- c(
        position_windows[[as.character(position)]],
        window_name
      )
    }
  }

  editing_induction_time <- as.numeric(params$editing_induction$timepoint)
  if(length(editing_induction_time) != 1 || !is.finite(editing_induction_time)){
    editing_induction_time <- Inf
  }
  adapter <- params$physicell_adapter
  profile_storage <- if(is.null(adapter$profile_storage)){
    'dense'
  } else{
    tolower(as.character(adapter$profile_storage))
  }
  if(length(profile_storage) != 1 ||
     !(profile_storage %in% c('dense', 'sparse'))){
    stop('physicell_adapter.profile_storage must be dense or sparse.')
  }
  compact_output <- isTRUE(adapter$compact_output)
  retain_internal_profiles <- if(is.null(adapter$retain_internal_profiles)){
    TRUE
  } else{
    isTRUE(adapter$retain_internal_profiles)
  }
  recorder_system <- if(is.null(adapter$recorder_system)){
    'generic barcode'
  } else{
    as.character(adapter$recorder_system)
  }
  if(length(recorder_system) != 1 || is.na(recorder_system) ||
     !nzchar(recorder_system)){
    stop('physicell_adapter.recorder_system must be one non-empty string.')
  }

  list(
    barcode_sequence = barcode_sequence,
    barcode_length = barcode_length,
    num_integrations = num_integrations,
    founder_label_sites = founder_label_sites,
    founder_label_positions = founder_label_positions,
    cell_type = cell_type,
    be_from = be_from,
    be_to = be_to,
    be_targets = be_targets,
    be_target_classes = be_expanded$classes,
    nuc_target_classes = nuc_expanded$classes,
    windows = all_windows,
    position_windows = position_windows,
    close_be_window = isTRUE(params$be_targets$editing_window$close_after_edit),
    close_nuc_window = isTRUE(params$nuclease_targets$editing_window$close_after_edit),
    editing_induction_time = editing_induction_time,
    recorder_system = recorder_system,
    profile_storage = profile_storage,
    compact_output = compact_output,
    retain_internal_profiles = retain_internal_profiles,
    output_positions = sort(unique(c(
      recording_positions,
      founder_label_positions
    ))),
    rate_sets = rate_sets
  )
}

#' Create an empty barcode profile in the model's storage representation
#'
#' @export
#' @param model Prepared recording model; `profile_storage` selects the form.
#' @return Either a `num_integrations`-by-`barcode_length` numeric matrix of
#'   zeros, or a `physicell_sparse_barcode` list holding one position-named
#'   numeric vector of edited sites per integration.
initialize_physicell_barcode_profile <- function(model){
  if(identical(model$profile_storage, 'sparse')){
    return(structure(
      replicate(
        model$num_integrations,
        setNames(numeric(), character()),
        simplify = FALSE
      ),
      class = c('physicell_sparse_barcode', 'list')
    ))
  }
  matrix(
    0,
    nrow = model$num_integrations,
    ncol = model$barcode_length
  )
}

#' Read one encoded allele out of a barcode profile
#'
#' @export
#' @param profile Dense matrix or `physicell_sparse_barcode` list.
#' @param integration Integration index.
#' @param position Barcode position.
#' @return The encoded allele, or `0` (reference) when a sparse profile holds no
#'   entry there.
physicell_barcode_value <- function(profile, integration, position){
  if(inherits(profile, 'physicell_sparse_barcode')){
    value <- profile[[integration]][as.character(position)]
    if(length(value) == 0 || is.na(value)){
      return(0)
    }
    return(unname(value))
  }
  profile[integration, position]
}

#' Write one encoded allele into a barcode profile
#'
#' @export
#' @param profile Dense matrix or `physicell_sparse_barcode` list.
#' @param integration Integration index.
#' @param position Barcode position.
#' @param value Encoded allele.
#' @return The updated profile, in the representation it arrived in.
set_physicell_barcode_value <- function(profile,
                                        integration,
                                        position,
                                        value){
  if(inherits(profile, 'physicell_sparse_barcode')){
    profile[[integration]][as.character(position)] <- value
  } else{
    profile[integration, position] <- value
  }
  profile
}

#' Report whether an editing window already carries an edit
#'
#' @export
#' @param profile Dense matrix or `physicell_sparse_barcode` list.
#' @param integration Integration index.
#' @param positions Barcode positions making up the editing window.
#' @return `TRUE` when any of those positions differs from reference.
physicell_barcode_window_edited <- function(profile,
                                            integration,
                                            positions){
  if(inherits(profile, 'physicell_sparse_barcode')){
    return(any(
      names(profile[[integration]]) %in% as.character(positions)
    ))
  }
  any(profile[integration, positions] != 0)
}

#' Stamp a founder's identity into a fresh barcode profile
#'
#' Writes `founder_index - 1` in binary across the model's founder-label
#' positions on integration 1, with the low bit at the first label position. A
#' `0` bit and a `1` bit are stored as the first and second non-reference base
#' at that position, so the label can be read back from the allele alone.
#'
#' @export
#' @param model Prepared recording model.
#' @param founder_index One-based rank of this founder; must lie within
#'   `num_founders`.
#' @param num_founders Total number of founders; `ceiling(log2())` of it may not
#'   exceed the configured label-site count.
#' @return A list with `profile` and an `events` data frame of
#'   `founder_barcode` records, empty when no label sites are configured.
initialize_physicell_founder_barcode <- function(model,
                                                  founder_index,
                                                  num_founders){
  founder_index <- as.integer(founder_index)
  num_founders <- as.integer(num_founders)
  if(length(founder_index) != 1 || is.na(founder_index) ||
     founder_index < 1 || founder_index > num_founders ||
     length(num_founders) != 1 || is.na(num_founders) || num_founders < 1){
    stop('Founder index and count must identify one valid founder.')
  }

  profile <- initialize_physicell_barcode_profile(model)
  label_sites <- as.integer(model$founder_label_sites)
  required_sites <- if(num_founders <= 1){
    0L
  } else{
    as.integer(ceiling(log(num_founders, base = 2)))
  }
  if(label_sites > 0 && label_sites < required_sites){
    stop(sprintf(
      paste(
        '%d founder-label sites cannot distinguish %d founders;',
        'at least %d sites are required.'
      ),
      label_sites,
      num_founders,
      required_sites
    ))
  }
  if(label_sites == 0){
    return(list(
      profile = profile,
      events = data.frame(
        integration = integer(),
        position = integer(),
        event = character(),
        reference = character(),
        alternate = character(),
        allele = numeric(),
        stringsAsFactors = FALSE
      )
    ))
  }

  bases <- c('A', 'G', 'C', 'T')
  founder_code <- founder_index - 1L
  events <- lapply(seq_len(label_sites), function(label_index){
    position <- model$founder_label_positions[label_index]
    bit <- (founder_code %/% (2^(label_index - 1L))) %% 2L
    reference_index <- match(model$barcode_sequence[position], bases)
    alternate_indices <- setdiff(seq_along(bases), reference_index)
    allele <- alternate_indices[bit + 1L]
    profile <<- set_physicell_barcode_value(
      profile,
      integration = 1L,
      position = position,
      value = allele
    )
    data.frame(
      integration = 1L,
      position = position,
      event = 'founder_barcode',
      reference = model$barcode_sequence[position],
      alternate = bases[allele],
      allele = as.numeric(allele),
      stringsAsFactors = FALSE
    )
  })
  list(profile = profile, events = do.call(rbind, events))
}

# ---- Hazards, selection coefficients, and progress logging ----

#' Convert per-division probabilities to continuous-time hazards
#'
#' Uses `-log(1 - p) / cell_cycle_length`, so a branch of exactly one cell cycle
#' reproduces the configured per-division probability.
#'
#' @export
#' @param probability Probabilities, optionally named; values are clamped into
#'   `[0, 1)`.
#' @param cell_cycle_length Cell-cycle duration in simulation time units.
#' @return Hazards per unit time, keeping the input names.
physicell_probability_hazard <- function(probability, cell_cycle_length){
  probability_names <- names(probability)
  probability <- pmin(pmax(as.numeric(probability), 0), 1 - .Machine$double.eps)
  hazards <- -log1p(-probability) / cell_cycle_length
  names(hazards) <- probability_names
  hazards
}

#' Resolve the non-Mendelian selection coefficient
#'
#' @export
#' @param params Parsed remote_mito parameter list.
#' @param marker Optional key inside `non_mendelian_selection` that overrides
#'   the global coefficient; must be one non-empty string when supplied.
#' @param default Value used when neither the marker key nor
#'   `non_mendelian_selection.coefficient` is present.
#' @return One coefficient validated to lie in `[0, 1]`.
non_mendelian_selection_coefficient <- function(params,
                                                 marker = NULL,
                                                 default = 0){
  if(!is.list(params)){
    stop('params must be a parsed remote_mito parameter list.')
  }
  configuration <- params$non_mendelian_selection
  if(is.null(configuration)){
    configuration <- list()
  }
  if(!is.list(configuration)){
    stop('non_mendelian_selection must be a JSON object.')
  }
  marker <- if(is.null(marker)) NULL else as.character(marker)
  if(!is.null(marker) &&
     (length(marker) != 1L || is.na(marker) || !nzchar(marker))){
    stop('marker must be NULL or one non-empty configuration key.')
  }
  value <- if(!is.null(marker) && !is.null(configuration[[marker]])){
    configuration[[marker]]
  } else if(!is.null(configuration$coefficient)){
    configuration$coefficient
  } else{
    default
  }
  value <- suppressWarnings(as.numeric(value))
  coefficient_name <- if(is.null(marker)){
    'non_mendelian_selection.coefficient'
  } else{
    paste0('non_mendelian_selection.', marker)
  }
  if(length(value) != 1L || !is.finite(value) || value < 0 || value > 1){
    stop(sprintf('%s must be one selection coefficient in [0, 1].',
                 coefficient_name))
  }
  value
}

#' Format an elapsed or remaining duration for progress messages
#'
#' @export
#' @param seconds One non-negative finite number of seconds.
#' @return `"unknown"` for anything else, otherwise seconds, minutes and
#'   seconds, or hours and minutes.
format_physicell_progress_duration <- function(seconds){
  if(length(seconds) != 1 || !is.finite(seconds) || seconds < 0){
    return('unknown')
  }
  if(seconds < 60){
    return(sprintf('%.1fs', seconds))
  }
  if(seconds < 3600){
    return(sprintf('%dm %02ds', floor(seconds / 60), round(seconds %% 60)))
  }
  sprintf(
    '%dh %02dm',
    floor(seconds / 3600),
    floor((seconds %% 3600) / 60)
  )
}

#' Emit one timestamped stage message
#'
#' @export
#' @param message_text One non-empty string.
#' @param enabled When `FALSE`, nothing is printed.
#' @return Invisibly, whether the message was emitted.
#' @section Side effects: Writes to the message connection.
physicell_log_stage <- function(message_text, enabled = TRUE){
  message_text <- as.character(message_text)
  if(length(message_text) != 1 || is.na(message_text) ||
     !nzchar(message_text)){
    stop('Stage log message must be one non-empty string.')
  }
  if(isTRUE(enabled)){
    message(sprintf(
      '[%s] %s',
      format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z'),
      message_text
    ))
  }
  invisible(isTRUE(enabled))
}

#' Create a throttled progress reporter
#'
#' Precomputes evenly spaced completion checkpoints and returns a closure that
#' prints only when one is crossed, reporting the count, percentage, elapsed
#' time, throughput, and ETA against the wall clock captured at construction.
#'
#' @export
#' @param total Total units of work; must be one positive integer.
#' @param label Phase name shown in every message.
#' @param enabled When `FALSE`, the closure still validates its argument but
#'   prints nothing.
#' @param updates Approximate number of messages to emit, capped at `total`.
#' @param unit Noun used for the work units, such as `nodes` or `events`.
#' @return A function of the completed count, which must lie in `[0, total]`,
#'   returning invisibly whether it printed.
#' @section Side effects: The returned closure writes to the message connection.
new_physicell_progress_reporter <- function(total,
                                            label,
                                            enabled = FALSE,
                                            updates = 20L,
                                            unit = 'nodes'){
  total <- as.integer(total)
  updates <- as.integer(updates)
  label <- as.character(label)
  unit <- as.character(unit)
  if(length(total) != 1 || is.na(total) || total < 1){
    stop('Progress total must be one positive integer.')
  }
  if(length(updates) != 1 || is.na(updates) || updates < 1){
    stop('Progress updates must be one positive integer.')
  }
  if(length(label) != 1 || is.na(label) || !nzchar(label)){
    stop('Progress label must be one non-empty string.')
  }
  if(length(unit) != 1 || is.na(unit) || !nzchar(unit)){
    stop('Progress unit must be one non-empty string.')
  }
  updates <- min(updates, total)

  checkpoints <- unique(c(
    0L,
    as.integer(ceiling(seq_len(updates) * total / updates))
  ))
  checkpoint_index <- 1L
  start_time <- unname(proc.time()[['elapsed']])

  function(completed){
    completed <- as.integer(completed)
    if(length(completed) != 1 || is.na(completed) ||
       completed < 0 || completed > total){
      stop(sprintf(
        'Progress for %s must be between 0 and %d.',
        label,
        total
      ))
    }
    if(!isTRUE(enabled) || checkpoint_index > length(checkpoints) ||
       completed < checkpoints[checkpoint_index]){
      return(invisible(FALSE))
    }

    elapsed <- max(unname(proc.time()[['elapsed']]) - start_time, 0)
    percent <- 100 * completed / total
    if(completed == 0){
      message(sprintf(
        '[%s] 0/%s %s (0.0%%): starting',
        label,
        format(total, big.mark = ',', scientific = FALSE, trim = TRUE),
        unit
      ))
    } else{
      rate <- completed / max(elapsed, 0.001)
      remaining <- if(is.finite(rate) && rate > 0){
        (total - completed) / rate
      } else{
        NA_real_
      }
      rate_text <- if(is.finite(rate)){
        paste0(
          format(round(rate, 1), big.mark = ',', scientific = FALSE, trim = TRUE),
          ' ',
          unit,
          '/s'
        )
      } else{
        'calculating rate'
      }
      message(sprintf(
        '[%s] %s/%s %s (%.1f%%): elapsed %s, %s, ETA %s',
        label,
        format(completed, big.mark = ',', scientific = FALSE, trim = TRUE),
        format(total, big.mark = ',', scientific = FALSE, trim = TRUE),
        unit,
        percent,
        format_physicell_progress_duration(elapsed),
        rate_text,
        format_physicell_progress_duration(remaining)
      ))
    }

    while(checkpoint_index <= length(checkpoints) &&
          checkpoints[checkpoint_index] <= completed){
      checkpoint_index <<- checkpoint_index + 1L
    }
    invisible(TRUE)
  }
}

# ---- Branch-segment mutation simulation ----

#' Build the empty core barcode-event table
#'
#' @export
#' @return A zero-row data frame with the columns `integration`, `position`,
#'   `event`, `reference`, `alternate`, `allele`, and `event_time`.
empty_physicell_barcode_events <- function(){
  data.frame(
    integration = integer(),
    position = integer(),
    event = character(),
    reference = character(),
    alternate = character(),
    allele = numeric(),
    event_time = numeric(),
    stringsAsFactors = FALSE
  )
}

#' Choose which event fires at each mutated position
#'
#' Positions whose rate set offers a single event are resolved straight from the
#' cached `single_event_names`; the rest draw one uniform and index into that
#' position's cached cumulative hazard-share vector.
#'
#' @export
#' @param positions Barcode positions that mutated on this segment.
#' @param rate_set Prepared rate set carrying the per-position event caches.
#' @return A character vector of event names, one per position, drawn from
#'   `sub_1`-`sub_4`, `insertion`, and `deletion`.
select_physicell_position_events <- function(positions, rate_set){
  selected_events <- unname(rate_set$single_event_names[positions])
  multiple_event_indices <- which(is.na(selected_events))
  if(length(multiple_event_indices) > 0){
    selection_draws <- stats::runif(length(multiple_event_indices))
    for(draw_index in seq_along(multiple_event_indices)){
      result_index <- multiple_event_indices[draw_index]
      position <- positions[result_index]
      cumulative <-
        rate_set$event_cumulative_probabilities[[position]]
      selected_index <- which(
        selection_draws[draw_index] <= cumulative
      )[1]
      selected_events[result_index] <-
        names(rate_set$event_hazards[[position]])[selected_index]
    }
  }
  selected_events
}

#' Turn selected event names into encoded alleles
#'
#' A substitution keeps its destination base index (`1`-`4` for `A`, `G`, `C`,
#' `T`), a deletion encodes as `-1`, and an insertion draws a uniform base and
#' encodes as that base index divided by ten.
#'
#' @export
#' @param selected_events Event names from
#'   `select_physicell_position_events()`; anything else raises an error.
#' @return A list with `event` (`substitution`, `deletion`, or `insertion`),
#'   `alternate` (the destination base, `-`, or `+<base>`), and `allele`.
decode_physicell_barcode_events <- function(selected_events){
  bases <- c('A', 'G', 'C', 'T')
  substitution <- grepl('^sub_[1-4]$', selected_events)
  deletion <- selected_events == 'deletion'
  insertion <- selected_events == 'insertion'
  if(any(!(substitution | deletion | insertion))){
    stop(sprintf(
      'Unknown recording event type: %s.',
      selected_events[which(!(substitution | deletion | insertion))[1]]
    ))
  }

  event_types <- rep.int('substitution', length(selected_events))
  alternates <- character(length(selected_events))
  alleles <- numeric(length(selected_events))
  if(any(substitution)){
    destinations <- as.integer(sub(
      '^sub_',
      '',
      selected_events[substitution]
    ))
    alternates[substitution] <- bases[destinations]
    alleles[substitution] <- destinations
  }
  if(any(deletion)){
    event_types[deletion] <- 'deletion'
    alternates[deletion] <- '-'
    alleles[deletion] <- -1
  }
  if(any(insertion)){
    inserted_bases <- sample.int(
      length(bases),
      sum(insertion),
      replace = TRUE
    )
    event_types[insertion] <- 'insertion'
    alternates[insertion] <- paste0('+', bases[inserted_bases])
    alleles[insertion] <- inserted_bases / 10
  }
  list(
    event = event_types,
    alternate = alternates,
    allele = alleles
  )
}

#' Apply prime-editing hazards over one branch segment
#'
#' For each integration only the still-unedited target positions compete. One
#' uniform per position both decides whether that target fires within `duration`
#' and, by inverting the exponential CDF, fixes the exact time it fired. An
#' edited target locks at state `1` and reports its assigned pegRNA's programmed
#' sequence.
#'
#' @export
#' @param profile Inherited profile, dense or sparse.
#' @param duration Length of the segment; a non-positive value is a no-op.
#' @param rate_set Prime-editing rate set with `total_hazards` named by target
#'   position.
#' @param model Prepared prime-editing model.
#' @param segment_start Absolute time at which the segment begins, added to the
#'   sampled waiting times.
#' @return A list with the updated `profile` and an `events` table whose
#'   `alternate` column holds the exact edited sequence.
mutate_prime_editing_segment <- function(profile,
                                         duration,
                                         rate_set,
                                         model,
                                         segment_start = 0){
  if(duration <= 0 || length(rate_set$active_positions) == 0L){
    return(list(
      profile = profile,
      events = empty_physicell_barcode_events()
    ))
  }
  active_positions <- rate_set$active_positions
  active_position_names <- as.character(active_positions)
  target_rows <- match(
    active_positions,
    model$prime_editing$targets$target_position
  )
  if(anyNA(target_rows)){
    stop('A prime-editing rate position has no assigned pegRNA.')
  }
  target_by_position <- setNames(
    target_rows,
    active_position_names
  )
  event_chunks <- vector('list', model$num_integrations)
  event_chunk_count <- 0L

  for(integration in seq_len(model$num_integrations)){
    if(inherits(profile, 'physicell_sparse_barcode')){
      edited_positions <- names(profile[[integration]])
      unedited <- if(length(edited_positions) == 0L){
        rep.int(TRUE, length(active_positions))
      } else{
        is.na(match(active_position_names, edited_positions))
      }
    } else{
      unedited <- profile[integration, active_positions] == 0
    }
    candidate_positions <- active_positions[unedited]
    if(length(candidate_positions) == 0L){
      next
    }
    hazards <- rate_set$total_hazards[as.character(candidate_positions)]
    mutation_draws <- stats::runif(length(candidate_positions))
    mutated <- mutation_draws < -expm1(-hazards * duration)
    if(!any(mutated)){
      next
    }
    mutated_positions <- candidate_positions[mutated]
    mutated_hazards <- hazards[mutated]
    event_times <- segment_start -
      log1p(-mutation_draws[mutated]) / mutated_hazards
    if(inherits(profile, 'physicell_sparse_barcode')){
      new_values <- rep.int(1, length(mutated_positions))
      names(new_values) <- as.character(mutated_positions)
      profile[[integration]][names(new_values)] <- new_values
    } else{
      profile[integration, mutated_positions] <- 1
    }
    assignment_rows <- unname(
      target_by_position[as.character(mutated_positions)]
    )
    assignments <- model$prime_editing$targets[
      assignment_rows,
      ,
      drop = FALSE
    ]
    event_chunk_count <- event_chunk_count + 1L
    event_chunks[[event_chunk_count]] <- data.frame(
      integration = rep.int(integration, length(mutated_positions)),
      position = mutated_positions,
      event = rep.int('prime_edit', length(mutated_positions)),
      reference = rep.int('unedited', length(mutated_positions)),
      alternate = assignments$edit_sequence,
      allele = rep.int(1, length(mutated_positions)),
      event_time = event_times,
      stringsAsFactors = FALSE
    )
  }
  events <- if(event_chunk_count == 0L){
    empty_physicell_barcode_events()
  } else if(event_chunk_count == 1L){
    event_chunks[[1L]]
  } else{
    do.call(rbind, event_chunks[seq_len(event_chunk_count)])
  }
  list(profile = profile, events = events)
}

#' Apply PALINCODE hazards over one branch segment
#'
#' Each still-wild-type cBit competes with one uniform that both decides whether
#' it fires within `duration` and fixes the exact event time. A second uniform
#' resolves the outcome against that cBit's left/right/both fractions, and the
#' resulting state (`1`, `2`, or `3`) is locked in.
#'
#' @export
#' @param profile Inherited profile, dense or sparse.
#' @param duration Length of the segment; a non-positive value is a no-op.
#' @param rate_set PALINCODE rate set with `total_hazards` and
#'   `outcome_fractions` indexed by cBit position.
#' @param model Prepared PALINCODE model.
#' @param segment_start Absolute time at which the segment begins.
#' @return A list with the updated `profile` and an `events` table whose `event`
#'   column is `palincode_left`, `palincode_right`, or `palincode_both`.
mutate_palincode_segment <- function(profile,
                                     duration,
                                     rate_set,
                                     model,
                                     segment_start = 0){
  if(duration <= 0 || length(rate_set$active_positions) == 0L){
    return(list(
      profile = profile,
      events = empty_physicell_barcode_events()
    ))
  }
  active_positions <- rate_set$active_positions
  active_position_names <- as.character(active_positions)
  event_chunks <- vector('list', model$num_integrations)
  event_chunk_count <- 0L
  outcome_names <- c('left', 'right', 'both')

  for(integration in seq_len(model$num_integrations)){
    if(inherits(profile, 'physicell_sparse_barcode')){
      edited_positions <- names(profile[[integration]])
      unedited <- if(length(edited_positions) == 0L){
        rep.int(TRUE, length(active_positions))
      } else{
        is.na(match(active_position_names, edited_positions))
      }
    } else{
      unedited <- profile[integration, active_positions] == 0
    }
    candidate_positions <- active_positions[unedited]
    if(length(candidate_positions) == 0L){
      next
    }
    hazards <- rate_set$total_hazards[candidate_positions]
    mutation_draws <- stats::runif(length(candidate_positions))
    mutated <- mutation_draws < -expm1(-hazards * duration)
    if(!any(mutated)){
      next
    }
    mutated_positions <- candidate_positions[mutated]
    mutated_hazards <- hazards[mutated]
    event_times <- segment_start -
      log1p(-mutation_draws[mutated]) / mutated_hazards
    outcome_draws <- stats::runif(length(mutated_positions))
    outcome_rows <- rate_set$outcome_fractions[mutated_positions, , drop = FALSE]
    outcome_states <- 1L +
      as.integer(outcome_draws > outcome_rows[, 'left']) +
      as.integer(
        outcome_draws > outcome_rows[, 'left'] + outcome_rows[, 'right']
      )
    if(inherits(profile, 'physicell_sparse_barcode')){
      new_values <- as.numeric(outcome_states)
      names(new_values) <- as.character(mutated_positions)
      profile[[integration]][names(new_values)] <- new_values
    } else{
      profile[integration, mutated_positions] <- outcome_states
    }
    outcomes <- outcome_names[outcome_states]
    event_chunk_count <- event_chunk_count + 1L
    event_chunks[[event_chunk_count]] <- data.frame(
      integration = rep.int(integration, length(mutated_positions)),
      position = mutated_positions,
      event = paste0('palincode_', outcomes),
      reference = rep.int('wild_type', length(mutated_positions)),
      alternate = outcomes,
      allele = as.numeric(outcome_states),
      event_time = event_times,
      stringsAsFactors = FALSE
    )
  }

  events <- if(event_chunk_count == 0L){
    empty_physicell_barcode_events()
  } else if(event_chunk_count == 1L){
    event_chunks[[1L]]
  } else{
    do.call(rbind, event_chunks[seq_len(event_chunk_count)])
  }
  list(profile = profile, events = events)
}

#' Apply barcode mutation over one constant-rate branch segment
#'
#' Delegates to the prime-editing or PALINCODE mutator when the model is one of
#' those recorders. Otherwise already-edited positions are skipped (recording is
#' irreversible) and the remaining active positions are drawn together: one
#' uniform per position decides whether it mutates within `duration` and fixes
#' the exact event time, then the event types are selected and decoded in a
#' batch and written back in one update. Hazard caches missing from an older
#' rate set are rebuilt here from `rate_set$probabilities`.
#'
#' @export
#' @param profile Inherited profile, dense or sparse.
#' @param duration Length of the segment; a non-positive value is a no-op.
#' @param rate_set Rate set for this segment's editing state.
#' @param model Prepared recording model.
#' @param segment_start Absolute time at which the segment begins.
#' @return A list with the updated `profile` and an `events` table.
#' @note When either editing window closes after an edit, the vectorized path
#'   cannot be used: positions are then visited one at a time in random order so
#'   an earlier edit can close its window against the later ones.
mutate_physicell_barcode_segment <- function(profile,
                                             duration,
                                             rate_set,
                                             model,
                                             segment_start = 0){
  if(isTRUE(model$is_prime_editing)){
    return(mutate_prime_editing_segment(
      profile,
      duration = duration,
      rate_set = rate_set,
      model = model,
      segment_start = segment_start
    ))
  }
  if(isTRUE(model$is_palincode)){
    return(mutate_palincode_segment(
      profile,
      duration = duration,
      rate_set = rate_set,
      model = model,
      segment_start = segment_start
    ))
  }
  if(duration <= 0){
    return(list(
      profile = profile,
      events = empty_physicell_barcode_events()
    ))
  }

  active_positions <- rate_set$active_positions
  if(length(active_positions) == 0){
    return(list(
      profile = profile,
      events = empty_physicell_barcode_events()
    ))
  }

  event_hazards <- rate_set$event_hazards
  if(is.null(event_hazards)){
    event_hazards <- lapply(
      rate_set$probabilities,
      physicell_probability_hazard,
      cell_cycle_length = rate_set$cell_cycle_length
    )
  }
  total_hazards <- rate_set$total_hazards
  if(is.null(total_hazards)){
    total_hazards <- vapply(event_hazards, sum, numeric(1))
  }
  event_cumulative_probabilities <-
    rate_set$event_cumulative_probabilities
  if(is.null(event_cumulative_probabilities)){
    event_cumulative_probabilities <- lapply(
      event_hazards,
      function(hazards){
        if(length(hazards) == 0){
          return(numeric())
        }
        cumulative <- cumsum(hazards / sum(hazards))
        cumulative[length(cumulative)] <- 1
        cumulative
      }
    )
  }
  single_event_names <- rate_set$single_event_names
  if(is.null(single_event_names)){
    single_event_names <- vapply(
      event_hazards,
      function(hazards){
        if(length(hazards) == 1){
          names(hazards)[1]
        } else{
          NA_character_
        }
      },
      character(1)
    )
  }
  active_position_names <- rate_set$active_position_names
  if(is.null(active_position_names)){
    active_position_names <- as.character(active_positions)
  }
  rate_set$event_hazards <- event_hazards
  rate_set$total_hazards <- total_hazards
  rate_set$event_cumulative_probabilities <-
    event_cumulative_probabilities
  rate_set$single_event_names <- single_event_names

  event_chunks <- vector('list', model$num_integrations)
  event_chunk_count <- 0L
  closure_enabled <-
    isTRUE(model$close_be_window) || isTRUE(model$close_nuc_window)

  if(!closure_enabled){
    for(integration in seq_len(model$num_integrations)){
      if(inherits(profile, 'physicell_sparse_barcode')){
        edited_positions <- names(profile[[integration]])
        unedited <- if(length(edited_positions) == 0){
          rep.int(TRUE, length(active_positions))
        } else{
          is.na(match(active_position_names, edited_positions))
        }
      } else{
        unedited <- profile[integration, active_positions] == 0
      }
      candidate_positions <- active_positions[unedited]
      if(length(candidate_positions) == 0){
        next
      }

      mutation_draws <- stats::runif(length(candidate_positions))
      mutation_probabilities <- -expm1(
        -total_hazards[candidate_positions] * duration
      )
      mutated <- mutation_draws < mutation_probabilities
      if(!any(mutated)){
        next
      }

      mutated_positions <- candidate_positions[mutated]
      mutated_draws <- mutation_draws[mutated]
      event_times <- segment_start -
        log1p(-mutated_draws) / total_hazards[mutated_positions]
      selected_events <- select_physicell_position_events(
        mutated_positions,
        rate_set
      )
      decoded <- decode_physicell_barcode_events(selected_events)

      if(inherits(profile, 'physicell_sparse_barcode')){
        integration_values <- profile[[integration]]
        new_values <- decoded$allele
        names(new_values) <- as.character(mutated_positions)
        integration_values[names(new_values)] <- new_values
        profile[[integration]] <- integration_values
      } else{
        profile[integration, mutated_positions] <- decoded$allele
      }

      event_chunk_count <- event_chunk_count + 1L
      event_chunks[[event_chunk_count]] <- list(
        integration = rep.int(integration, length(mutated_positions)),
        position = mutated_positions,
        event = decoded$event,
        reference = model$barcode_sequence[mutated_positions],
        alternate = decoded$alternate,
        allele = decoded$allele,
        event_time = event_times
      )
    }
  } else{
    for(integration in seq_len(model$num_integrations)){
      for(position in sample(active_positions)){
        if(physicell_barcode_value(profile, integration, position) != 0){
          next
        }

        relevant_windows <- model$position_windows[[as.character(position)]]
        if(length(relevant_windows) > 0){
          window_closed <- any(vapply(
            relevant_windows,
            function(window_name){
              close_this_window <-
                (
                  isTRUE(model$close_be_window) &&
                    startsWith(window_name, 'be_window_')
                ) ||
                (
                  isTRUE(model$close_nuc_window) &&
                    startsWith(window_name, 'nuc_window_')
                )
              close_this_window &&
                physicell_barcode_window_edited(
                  profile,
                  integration,
                  model$windows[[window_name]]
                )
            },
            logical(1)
          ))
          if(window_closed){
            next
          }
        }

        total_hazard <- total_hazards[position]
        mutation_probability <- -expm1(-total_hazard * duration)
        mutation_draw <- stats::runif(1)
        if(mutation_draw >= mutation_probability){
          next
        }

        selected_event <- select_physicell_position_events(
          position,
          rate_set
        )
        decoded <- decode_physicell_barcode_events(selected_event)
        profile <- set_physicell_barcode_value(
          profile,
          integration,
          position,
          decoded$allele
        )
        event_chunk_count <- event_chunk_count + 1L
        event_chunks[[event_chunk_count]] <- list(
          integration = integration,
          position = position,
          event = decoded$event,
          reference = model$barcode_sequence[position],
          alternate = decoded$alternate,
          allele = decoded$allele,
          event_time = segment_start -
            log1p(-mutation_draw) / total_hazard
        )
      }
    }
  }

  events <- if(event_chunk_count == 0L){
    empty_physicell_barcode_events()
  } else{
    retained_chunks <- event_chunks[seq_len(event_chunk_count)]
    data.frame(
      integration = as.integer(unlist(lapply(
        retained_chunks,
        `[[`,
        'integration'
      ), use.names = FALSE)),
      position = as.integer(unlist(lapply(
        retained_chunks,
        `[[`,
        'position'
      ), use.names = FALSE)),
      event = unlist(lapply(
        retained_chunks,
        `[[`,
        'event'
      ), use.names = FALSE),
      reference = unlist(lapply(
        retained_chunks,
        `[[`,
        'reference'
      ), use.names = FALSE),
      alternate = unlist(lapply(
        retained_chunks,
        `[[`,
        'alternate'
      ), use.names = FALSE),
      allele = as.numeric(unlist(lapply(
        retained_chunks,
        `[[`,
        'allele'
      ), use.names = FALSE)),
      event_time = as.numeric(unlist(lapply(
        retained_chunks,
        `[[`,
        'event_time'
      ), use.names = FALSE)),
      stringsAsFactors = FALSE
    )
  }
  list(profile = profile, events = events)
}

#' Mutate one branch, splitting it at editing induction
#'
#' A branch that straddles the model's global induction time is simulated as an
#' uninduced segment followed by an induced one; otherwise it is a single
#' segment.
#'
#' @export
#' @param profile Inherited profile.
#' @param start_time Branch start in simulation time.
#' @param end_time Branch end in simulation time.
#' @param model Prepared recording model supplying `editing_induction_time` and
#'   both rate sets.
#' @param editing_state `auto` to split at the induction time, or `induced` /
#'   `uninduced` to force one state across the whole branch.
#' @return A list with the updated `profile` and an `events` table carrying
#'   `segment_start`, `segment_end`, and `editing_state` provenance columns.
mutate_physicell_barcode_branch <- function(profile,
                                            start_time,
                                            end_time,
                                            model,
                                            editing_state = 'auto'){
  allowed_states <- c('auto', 'induced', 'uninduced')
  if(!(editing_state %in% allowed_states)){
    stop(sprintf(
      'editing_state must be one of: %s.',
      paste(allowed_states, collapse = ', ')
    ))
  }
  if(editing_state == 'induced'){
    segment_starts <- start_time
    segment_ends <- end_time
    segment_states <- 'induced_editing_params'
  } else if(editing_state == 'uninduced'){
    segment_starts <- start_time
    segment_ends <- end_time
    segment_states <- 'uninduced_editing_params'
  } else if(end_time <= model$editing_induction_time){
    segment_starts <- start_time
    segment_ends <- end_time
    segment_states <- 'uninduced_editing_params'
  } else if(start_time >= model$editing_induction_time){
    segment_starts <- start_time
    segment_ends <- end_time
    segment_states <- 'induced_editing_params'
  } else{
    segment_starts <- c(start_time, model$editing_induction_time)
    segment_ends <- c(model$editing_induction_time, end_time)
    segment_states <- c(
      'uninduced_editing_params',
      'induced_editing_params'
    )
  }

  all_events <- list()
  for(segment_index in seq_along(segment_starts)){
    mutation <- mutate_physicell_barcode_segment(
      profile,
      duration = segment_ends[segment_index] -
        segment_starts[segment_index],
      rate_set = model$rate_sets[[segment_states[segment_index]]],
      model = model,
      segment_start = segment_starts[segment_index]
    )
    profile <- mutation$profile
    if(nrow(mutation$events) > 0){
      mutation$events$segment_start <- segment_starts[segment_index]
      mutation$events$segment_end <- segment_ends[segment_index]
      mutation$events$editing_state <- segment_states[segment_index]
      all_events[[length(all_events) + 1L]] <- mutation$events
    }
  }

  events <- if(length(all_events) == 0){
    empty_events <- empty_physicell_barcode_events()
    empty_events$segment_start <- numeric()
    empty_events$segment_end <- numeric()
    empty_events$editing_state <- character()
    empty_events
  } else if(length(all_events) == 1L){
    all_events[[1]]
  } else{
    do.call(rbind, all_events)
  }
  list(profile = profile, events = events)
}

# ---- Lineage replay ----

#' Replay barcode recording along an imported PhysiCell lineage
#'
#' Walks the node table in order - parents precede their children, so a parent's
#' profile is always available - initializing founder profiles and mutating each
#' branch from the profile it inherited. Mutation events accumulate in
#' preallocated typed vectors that grow geometrically and are turned into one
#' data frame at the end, sorted by event time.
#'
#' @export
#' @param nodes Event-resolved lineage node table; needs `node_id`,
#'   `physicell_id`, `parent_node_id`, `birth_time`, `end_time`, and
#'   `is_terminal`, and may also carry `cell_type` and `editing_state` columns.
#' @param model One prepared recording model, or a named collection of them
#'   keyed by cell type, which requires a `cell_type` column on `nodes`.
#' @param editing_state Default branch policy: `auto`, `induced`, or
#'   `uninduced`; a per-node `editing_state` column overrides it.
#' @param terminal_physicell_ids Optional live-cell IDs forming the sampled
#'   terminal set; each must be terminal in the lineage.
#' @param seed Seed set before the simulation begins.
#' @param show_progress Whether to emit progress messages.
#' @param progress_updates Approximate number of progress messages to emit.
#' @return A list with `nodes`, `profiles` keyed by node ID, `terminal_nodes`,
#'   `mutation_events`, `editing_state`, and `seed`.
#' @section Side effects: Calls `set.seed(seed)`.
#' @note When the model does not retain internal profiles, a parent's profile is
#'   released as soon as its last child has been simulated, so `profiles` then
#'   holds only the leaves and any still-pending branches.
simulate_recording_on_physicell_lineage <- function(nodes,
                                                    model,
                                                    editing_state = 'auto',
                                                    terminal_physicell_ids = NULL,
                                                    seed = 1,
                                                    show_progress = FALSE,
                                                    progress_updates = 20L){
  if(!is.data.frame(nodes) || !all(c(
    'node_id', 'physicell_id', 'parent_node_id', 'birth_time',
    'end_time', 'is_terminal'
  ) %in% names(nodes))){
    stop('nodes is not a valid PhysiCell lineage node table.')
  }
  model_is_collection <- is.null(model$barcode_length)
  if(model_is_collection){
    if(length(model) == 0 || is.null(names(model)) ||
       any(!nzchar(names(model))) ||
       any(vapply(model, function(value) is.null(value$barcode_length), logical(1)))){
      stop('A barcode model collection must contain named prepared models.')
    }
    if(!('cell_type' %in% names(nodes))){
      stop('Lineage nodes need a cell_type column when using model collections.')
    }
    missing_models <- setdiff(unique(as.character(nodes$cell_type)), names(model))
    if(length(missing_models) > 0){
      stop(sprintf(
        'No barcode recording model is available for cell type(s): %s.',
        paste(missing_models, collapse = ', ')
      ))
    }
    output_model <- model[[1]]
  } else{
    output_model <- model
  }

  # Internal: pick the recording model that applies to one node
  #
  # @param node One-row slice of the node table.
  # @return The single prepared model, or the collection entry matching that
  #   node's `cell_type`.
  node_model <- function(node){
    if(model_is_collection){
      model[[as.character(node$cell_type)]]
    } else{
      model
    }
  }
  set.seed(seed)
  profiles <- list()
  event_capacity <- 4096L
  event_count <- 0L
  event_integration <- integer(event_capacity)
  event_position <- integer(event_capacity)
  event_type <- rep(NA_character_, event_capacity)
  event_reference <- rep(NA_character_, event_capacity)
  event_alternate <- rep(NA_character_, event_capacity)
  event_allele <- numeric(event_capacity)
  event_time <- numeric(event_capacity)
  event_segment_start <- numeric(event_capacity)
  event_segment_end <- numeric(event_capacity)
  event_editing_state <- rep(NA_character_, event_capacity)
  event_node_id <- rep(NA_character_, event_capacity)
  event_physicell_id <- rep(NA_character_, event_capacity)
  event_parent_node_id <- rep(NA_character_, event_capacity)
  event_branch_start <- numeric(event_capacity)
  event_branch_end <- numeric(event_capacity)

  # Internal: grow every typed event column to hold more events
  #
  # @param required_capacity Number of event rows that must fit.
  # @return `NULL`, invisibly; each round lengthens the enclosing event vectors
  #   by at least half again.
  grow_event_buffer <- function(required_capacity){
    if(required_capacity <= event_capacity){
      return(invisible(NULL))
    }
    new_capacity <- event_capacity
    while(new_capacity < required_capacity){
      new_capacity <- as.integer(min(
        max(ceiling(as.double(new_capacity) * 1.5), required_capacity),
        .Machine$integer.max
      ))
      if(new_capacity <= event_capacity){
        stop('Barcode mutation-event buffer exceeded R integer capacity.')
      }
    }
    length(event_integration) <<- new_capacity
    length(event_position) <<- new_capacity
    length(event_type) <<- new_capacity
    length(event_reference) <<- new_capacity
    length(event_alternate) <<- new_capacity
    length(event_allele) <<- new_capacity
    length(event_time) <<- new_capacity
    length(event_segment_start) <<- new_capacity
    length(event_segment_end) <<- new_capacity
    length(event_editing_state) <<- new_capacity
    length(event_node_id) <<- new_capacity
    length(event_physicell_id) <<- new_capacity
    length(event_parent_node_id) <<- new_capacity
    length(event_branch_start) <<- new_capacity
    length(event_branch_end) <<- new_capacity
    event_capacity <<- new_capacity
    invisible(NULL)
  }

  # Internal: copy one branch's events into the preallocated columns
  #
  # @param events Event table from `mutate_physicell_barcode_branch()`.
  # @param node_id Node the events belong to.
  # @param physicell_id_value PhysiCell ID of that node.
  # @param parent_node_id Parent node ID, or `NA` for a founder.
  # @param branch_start Branch start time.
  # @param branch_end Branch end time.
  # @return `NULL`, invisibly; `event_count` advances by `nrow(events)`.
  append_barcode_events <- function(events,
                                    node_id,
                                    physicell_id_value,
                                    parent_node_id,
                                    branch_start,
                                    branch_end){
    num_events <- nrow(events)
    if(num_events == 0){
      return(invisible(NULL))
    }
    required_capacity <- event_count + num_events
    grow_event_buffer(required_capacity)
    output_indices <- seq.int(event_count + 1L, required_capacity)
    event_integration[output_indices] <<- as.integer(events$integration)
    event_position[output_indices] <<- as.integer(events$position)
    event_type[output_indices] <<- as.character(events$event)
    event_reference[output_indices] <<- as.character(events$reference)
    event_alternate[output_indices] <<- as.character(events$alternate)
    event_allele[output_indices] <<- as.numeric(events$allele)
    event_time[output_indices] <<- as.numeric(events$event_time)
    event_segment_start[output_indices] <<-
      as.numeric(events$segment_start)
    event_segment_end[output_indices] <<- as.numeric(events$segment_end)
    event_editing_state[output_indices] <<-
      as.character(events$editing_state)
    event_node_id[output_indices] <<- as.character(node_id)
    event_physicell_id[output_indices] <<-
      as.character(physicell_id_value)
    event_parent_node_id[output_indices] <<-
      as.character(parent_node_id)
    event_branch_start[output_indices] <<- as.numeric(branch_start)
    event_branch_end[output_indices] <<- as.numeric(branch_end)
    event_count <<- required_capacity
    invisible(NULL)
  }
  root_indices <- which(is.na(nodes$parent_node_id))
  founder_number <- 0L
  remaining_children <- table(
    nodes$parent_node_id[!is.na(nodes$parent_node_id)]
  )
  report_progress <- new_physicell_progress_reporter(
    total = nrow(nodes),
    label = if(is.null(output_model$recorder_system)){
      'barcode lineage recording'
    } else{
      output_model$recorder_system
    },
    enabled = show_progress,
    updates = progress_updates
  )
  report_progress(0L)

  for(node_index in seq_len(nrow(nodes))){
    node <- nodes[node_index, , drop = FALSE]
    branch_model <- node_model(node)
    if(is.na(node$parent_node_id)){
      founder_number <- founder_number + 1L
      initialized <- initialize_physicell_founder_barcode(
        branch_model,
        founder_index = founder_number,
        num_founders = length(root_indices)
      )
      start_profile <- initialized$profile
      if(nrow(initialized$events) > 0){
        initialized$events$segment_start <- node$birth_time
        initialized$events$segment_end <- node$birth_time
        initialized$events$event_time <- node$birth_time
        initialized$events$editing_state <- 'founder_label'
        append_barcode_events(
          initialized$events,
          node_id = node$node_id,
          physicell_id_value = node$physicell_id,
          parent_node_id = NA_character_,
          branch_start = node$birth_time,
          branch_end = node$birth_time
        )
      }
    } else{
      parent_node_id <- as.character(node$parent_node_id)
      start_profile <- profiles[[parent_node_id]]
      if(is.null(start_profile)){
        stop(sprintf(
          'Parent profile %s was not available before child %s.',
          node$parent_node_id,
          node$node_id
        ))
      }
      remaining_children[parent_node_id] <-
        remaining_children[parent_node_id] - 1L
      if(!isTRUE(output_model$retain_internal_profiles) &&
         remaining_children[parent_node_id] <= 0){
        profiles[[parent_node_id]] <- NULL
      }
    }

    branch_editing_state <- if('editing_state' %in% names(nodes)){
      as.character(node$editing_state)
    } else{
      editing_state
    }
    mutation <- mutate_physicell_barcode_branch(
      start_profile,
      node$birth_time,
      node$end_time,
      branch_model,
      branch_editing_state
    )
    profiles[[node$node_id]] <- mutation$profile
    if(nrow(mutation$events) > 0){
      append_barcode_events(
        mutation$events,
        node_id = node$node_id,
        physicell_id_value = node$physicell_id,
        parent_node_id = node$parent_node_id,
        branch_start = node$birth_time,
        branch_end = node$end_time
      )
    }
    report_progress(node_index)
  }

  terminal_rows <- nodes[nodes$is_terminal, , drop = FALSE]
  if(!is.null(terminal_physicell_ids)){
    terminal_physicell_ids <- unique(
      physicell_id(terminal_physicell_ids, 'terminal_physicell_ids')
    )
    unknown_ids <- setdiff(terminal_physicell_ids, terminal_rows$physicell_id)
    if(length(unknown_ids) > 0){
      stop(sprintf(
        'Live-cell ID(s) absent from the reconstructed lineage: %s.',
        paste(unknown_ids, collapse = ', ')
      ))
    }
    terminal_rows <- terminal_rows[
      terminal_rows$physicell_id %in% terminal_physicell_ids,
      ,
      drop = FALSE
    ]
  }
  if(nrow(terminal_rows) == 0){
    stop('No terminal cells remain for lineage-recording output.')
  }

  events <- if(event_count == 0L){
    data.frame(
      integration = integer(),
      position = integer(),
      event = character(),
      reference = character(),
      alternate = character(),
      allele = numeric(),
      event_time = numeric(),
      segment_start = numeric(),
      segment_end = numeric(),
      editing_state = character(),
      node_id = character(),
      physicell_id = character(),
      parent_node_id = character(),
      branch_start = numeric(),
      branch_end = numeric(),
      stringsAsFactors = FALSE
    )
  } else{
    length(event_integration) <- event_count
    length(event_position) <- event_count
    length(event_type) <- event_count
    length(event_reference) <- event_count
    length(event_alternate) <- event_count
    length(event_allele) <- event_count
    length(event_time) <- event_count
    length(event_segment_start) <- event_count
    length(event_segment_end) <- event_count
    length(event_editing_state) <- event_count
    length(event_node_id) <- event_count
    length(event_physicell_id) <- event_count
    length(event_parent_node_id) <- event_count
    length(event_branch_start) <- event_count
    length(event_branch_end) <- event_count
    data.frame(
      integration = event_integration,
      position = event_position,
      event = event_type,
      reference = event_reference,
      alternate = event_alternate,
      allele = event_allele,
      event_time = event_time,
      segment_start = event_segment_start,
      segment_end = event_segment_end,
      editing_state = event_editing_state,
      node_id = event_node_id,
      physicell_id = event_physicell_id,
      parent_node_id = event_parent_node_id,
      branch_start = event_branch_start,
      branch_end = event_branch_end,
      stringsAsFactors = FALSE
    )
  }
  if(nrow(events) > 0){
    events <- events[
      order(
        events$event_time,
        events$node_id,
        events$integration,
        events$position
      ),
      ,
      drop = FALSE
    ]
    rownames(events) <- NULL
  }

  list(
    nodes = nodes,
    profiles = profiles,
    terminal_nodes = terminal_rows,
    mutation_events = events,
    editing_state = editing_state,
    seed = seed
  )
}

# ---- Descendant intervals and event-inheritance matrices ----

#' Order the sampled cells depth-first and index each node's descendants
#'
#' One iterative depth-first pass numbers the sampled terminal cells in
#' traversal order, which makes every node's sampled descendants a contiguous
#' run of that ordering. Recording `[start, end]` per node then answers "which
#' sampled cells inherit an event on this branch" in constant time.
#'
#' @export
#' @param nodes Lineage node table with unique `node_id` and `parent_node_id`,
#'   ordered so every parent precedes its children.
#' @param terminal_nodes Sampled terminal rows with unique `node_id` and
#'   `physicell_id`; all must be reachable from a root.
#' @return A list with `descendant_start`, `descendant_end`, and
#'   `descendant_count` (one entry per node), `terminal_order` (rows of
#'   `terminal_nodes` in traversal order), and `sample_ids` (`cell_<id>` in that
#'   same order).
physicell_terminal_descendant_intervals <- function(nodes, terminal_nodes){
  required_node_columns <- c('node_id', 'parent_node_id')
  required_terminal_columns <- c('node_id', 'physicell_id')
  if(!is.data.frame(nodes) ||
     !all(required_node_columns %in% names(nodes))){
    stop('nodes must contain node_id and parent_node_id.')
  }
  if(!is.data.frame(terminal_nodes) ||
     !all(required_terminal_columns %in% names(terminal_nodes))){
    stop('terminal_nodes must contain node_id and physicell_id.')
  }
  if(anyNA(nodes$node_id) || anyDuplicated(nodes$node_id)){
    stop('Lineage node IDs must be unique and non-missing.')
  }
  if(anyNA(terminal_nodes$node_id) ||
     anyDuplicated(terminal_nodes$node_id) ||
     anyNA(terminal_nodes$physicell_id) ||
     anyDuplicated(terminal_nodes$physicell_id)){
    stop('Sampled terminal node and PhysiCell IDs must be unique and non-missing.')
  }

  num_nodes <- nrow(nodes)
  terminal_node_indices <- match(terminal_nodes$node_id, nodes$node_id)
  if(anyNA(terminal_node_indices)){
    stop('One or more sampled terminal nodes are absent from the lineage.')
  }
  parent_indices <- match(nodes$parent_node_id, nodes$node_id)
  non_root_indices <- which(!is.na(nodes$parent_node_id))
  if(anyNA(parent_indices[non_root_indices])){
    stop('One or more lineage parents are absent from the node table.')
  }
  if(any(parent_indices[non_root_indices] >= non_root_indices)){
    stop('Lineage nodes must be ordered with every parent before its children.')
  }
  root_indices <- which(is.na(nodes$parent_node_id))
  if(num_nodes > 0 && length(root_indices) == 0){
    stop('Lineage nodes do not contain a root.')
  }

  child_counts <- tabulate(parent_indices[non_root_indices], nbins = num_nodes)
  child_order <- if(length(non_root_indices) == 0){
    integer()
  } else{
    non_root_indices[
      order(parent_indices[non_root_indices], method = 'radix')
    ]
  }
  child_ends <- cumsum(child_counts)
  child_starts <- child_ends - child_counts + 1L
  terminal_row_by_node <- rep(NA_integer_, num_nodes)
  terminal_row_by_node[terminal_node_indices] <- seq_len(nrow(terminal_nodes))

  descendant_start <- integer(num_nodes)
  descendant_end <- integer(num_nodes)
  ordered_terminal_rows <- integer(nrow(terminal_nodes))
  terminal_count <- 0L
  stack <- integer(max(2L * num_nodes, 1L))
  stack_top <- 0L

  # Internal: push one traversal event onto the explicit stack
  #
  # @param node_event Node index to visit, or its negation to close that node.
  # @return Called for its effect on the enclosing `stack` and `stack_top`.
  push_node <- function(node_event){
    stack_top <<- stack_top + 1L
    stack[stack_top] <<- node_event
  }

  for(root_position in rev(seq_along(root_indices))){
    push_node(root_indices[root_position])
  }
  while(stack_top > 0L){
    node_event <- stack[stack_top]
    stack_top <- stack_top - 1L
    if(node_event < 0L){
      descendant_end[-node_event] <- terminal_count
      next
    }

    node_index <- node_event
    descendant_start[node_index] <- terminal_count + 1L
    push_node(-node_index)
    if(child_counts[node_index] == 0L){
      terminal_row <- terminal_row_by_node[node_index]
      if(!is.na(terminal_row)){
        terminal_count <- terminal_count + 1L
        ordered_terminal_rows[terminal_count] <- terminal_row
      }
      next
    }

    child_start <- child_starts[node_index]
    child_end <- child_ends[node_index]
    for(child_position in seq.int(child_end, child_start)){
      push_node(child_order[child_position])
    }
  }
  if(terminal_count != nrow(terminal_nodes)){
    stop('One or more sampled terminal nodes were not reachable from a root.')
  }

  descendant_count <- pmax.int(
    descendant_end - descendant_start + 1L,
    0L
  )
  list(
    descendant_start = descendant_start,
    descendant_end = descendant_end,
    descendant_count = descendant_count,
    terminal_order = ordered_terminal_rows,
    sample_ids = paste0(
      'cell_',
      terminal_nodes$physicell_id[ordered_terminal_rows]
    )
  )
}

#' Build the sparse cell-by-event inheritance matrix
#'
#' A mutation event is inherited by exactly the sampled cells descending from
#' the node it occurred on, which the depth-first intervals make a contiguous
#' block of rows. The `dgCMatrix` is therefore assembled directly from its
#' column pointers and row indices, without materializing triplets.
#'
#' @export
#' @param nodes Lineage node table.
#' @param terminal_nodes Sampled terminal rows; these become the matrix rows.
#' @param mutation_events Event table containing a `node_id` column; these
#'   become the matrix columns.
#' @param event_ids Optional unique column names; generated as
#'   `mutation_event_<n>` when `NULL`.
#' @param show_progress Whether to emit progress messages.
#' @param progress_updates Approximate number of progress messages to emit.
#' @return A list with `matrix` (a `dgCMatrix` of ones), `descendant_counts` per
#'   event, and `terminal_order`.
#' @note Requires the `Matrix` package, and errors when the nonzero count would
#'   exceed the sparse integer-index limit.
physicell_event_descendant_matrix <- function(nodes,
                                              terminal_nodes,
                                              mutation_events,
                                              event_ids = NULL,
                                              show_progress = FALSE,
                                              progress_updates = 20L){
  if(!requireNamespace('Matrix', quietly = TRUE)){
    stop('The Matrix package is required for event-descendant matrix output.')
  }
  if(!is.data.frame(mutation_events) ||
     !('node_id' %in% names(mutation_events))){
    stop('mutation_events must be a data frame containing node_id.')
  }
  num_events <- nrow(mutation_events)
  if(is.null(event_ids)){
    event_width <- max(8L, nchar(max(num_events, 1L)))
    event_ids <- if(num_events == 0){
      character()
    } else{
      paste0(
        'mutation_event_',
        formatC(seq_len(num_events), width = event_width, flag = '0')
      )
    }
  }
  event_ids <- as.character(event_ids)
  if(length(event_ids) != num_events || anyNA(event_ids) ||
     any(!nzchar(event_ids)) || anyDuplicated(event_ids)){
    stop('event_ids must provide one unique non-empty ID per mutation event.')
  }

  intervals <- physicell_terminal_descendant_intervals(
    nodes,
    terminal_nodes
  )
  event_node_indices <- match(mutation_events$node_id, nodes$node_id)
  if(anyNA(event_node_indices)){
    stop(sprintf(
      'Mutation event node %s is absent from the lineage.',
      mutation_events$node_id[which(is.na(event_node_indices))[1]]
    ))
  }
  descendant_counts <- intervals$descendant_count[event_node_indices]
  num_entries <- sum(as.double(descendant_counts))
  if(num_entries > .Machine$integer.max){
    stop(sprintf(
      paste(
        'The literal event-descendant matrix requires %s nonzero entries,',
        'which exceeds the sparse Matrix integer-index limit.'
      ),
      format(num_entries, big.mark = ',', scientific = FALSE, trim = TRUE)
    ))
  }
  num_entries <- as.integer(num_entries)
  column_pointers <- as.integer(c(0, cumsum(as.double(descendant_counts))))
  row_indices <- integer(num_entries)
  report_progress <- if(num_events > 0){
    new_physicell_progress_reporter(
      total = num_events,
      label = 'Mutation-event descendant matrix',
      enabled = show_progress,
      updates = progress_updates,
      unit = 'events'
    )
  } else{
    NULL
  }
  if(!is.null(report_progress)){
    report_progress(0L)
  }

  if(num_entries > 0){
    for(event_index in seq_len(num_events)){
      descendant_count <- descendant_counts[event_index]
      if(descendant_count > 0L){
        output_indices <- seq.int(
          column_pointers[event_index] + 1L,
          column_pointers[event_index + 1L]
        )
        first_row <- intervals$descendant_start[
          event_node_indices[event_index]
        ]
        row_indices[output_indices] <- seq.int(
          first_row - 1L,
          first_row + descendant_count - 2L
        )
      }
      report_progress(event_index)
    }
  } else if(!is.null(report_progress)){
    for(event_index in seq_len(num_events)){
      report_progress(event_index)
    }
  }

  event_matrix <- methods::new(
    'dgCMatrix',
    i = row_indices,
    p = column_pointers,
    Dim = as.integer(c(nrow(terminal_nodes), num_events)),
    Dimnames = list(intervals$sample_ids, event_ids),
    x = rep.int(1, num_entries)
  )
  list(
    matrix = event_matrix,
    descendant_counts = as.integer(descendant_counts),
    terminal_order = intervals$terminal_order
  )
}

#' Write the combined event-inheritance matrix and its manifest
#'
#' Prefixes each modality's events with `event_id`, `modality`, and
#' `modality_event_row`, unions the column sets across modalities so they can be
#' row-bound, then builds one matrix spanning all the events together.
#'
#' @export
#' @param nodes Lineage node table.
#' @param terminal_nodes Sampled terminal rows.
#' @param event_tables Uniquely named list of event data frames, each with a
#'   `node_id` column and none of the reserved manifest columns.
#' @param output_dir Destination directory, created if needed.
#' @param show_progress Whether to emit stage and progress messages.
#' @param progress_updates Approximate number of progress messages to emit.
#' @param compress_csv Whether the manifest CSV is gzip-compressed.
#' @return Invisibly, a list with the sparse `matrix` and the `manifest`.
#' @section Side effects: Writes
#'   `mutation_event_descendant_matrix_sparse.rds` and
#'   `mutation_event_descendant_manifest.csv` under `output_dir`.
write_physicell_event_descendant_outputs <- function(
    nodes,
    terminal_nodes,
    event_tables,
    output_dir,
    show_progress = FALSE,
    progress_updates = 20L,
    compress_csv = TRUE){
  if(!is.list(event_tables) || length(event_tables) == 0 ||
     is.null(names(event_tables)) || any(!nzchar(names(event_tables))) ||
     anyDuplicated(names(event_tables))){
    stop('event_tables must be a uniquely named list of mutation-event tables.')
  }
  reserved_columns <- c(
    'event_id',
    'modality',
    'modality_event_row',
    'descendant_terminal_cells'
  )
  event_frames <- vector('list', length(event_tables))
  for(table_index in seq_along(event_tables)){
    modality <- names(event_tables)[table_index]
    events <- event_tables[[table_index]]
    if(!is.data.frame(events) || !('node_id' %in% names(events))){
      stop(sprintf(
        'Event table %s must be a data frame containing node_id.',
        modality
      ))
    }
    conflicts <- intersect(reserved_columns, names(events))
    if(length(conflicts) > 0){
      stop(sprintf(
        'Event table %s uses reserved column(s): %s.',
        modality,
        paste(conflicts, collapse = ', ')
      ))
    }
    num_events <- nrow(events)
    event_width <- max(8L, nchar(max(num_events, 1L)))
    event_ids <- if(num_events == 0){
      character()
    } else{
      paste0(
        modality,
        '_event_',
        formatC(seq_len(num_events), width = event_width, flag = '0')
      )
    }
    event_frames[[table_index]] <- cbind(
      data.frame(
        event_id = event_ids,
        modality = rep(modality, num_events),
        modality_event_row = seq_len(num_events),
        stringsAsFactors = FALSE
      ),
      events
    )
  }

  manifest_columns <- unique(unlist(lapply(event_frames, names)))
  event_frames <- lapply(event_frames, function(event_frame){
    missing_columns <- setdiff(manifest_columns, names(event_frame))
    for(column in missing_columns){
      event_frame[[column]] <- rep(NA, nrow(event_frame))
    }
    event_frame[, manifest_columns, drop = FALSE]
  })
  combined_events <- do.call(rbind, event_frames)
  rownames(combined_events) <- NULL

  physicell_log_stage(
    sprintf(
      paste(
        'Mutation-event output: building a %s-cell by %s-event',
        'literal descendant matrix.'
      ),
      format(
        nrow(terminal_nodes),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      ),
      format(
        nrow(combined_events),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      )
    ),
    enabled = show_progress
  )
  matrix_output <- physicell_event_descendant_matrix(
    nodes,
    terminal_nodes,
    combined_events,
    event_ids = combined_events$event_id,
    show_progress = show_progress,
    progress_updates = progress_updates
  )
  event_matrix <- matrix_output$matrix
  combined_events$descendant_terminal_cells <-
    matrix_output$descendant_counts
  combined_events <- combined_events[, c(
    reserved_columns,
    setdiff(names(combined_events), reserved_columns)
  ), drop = FALSE]

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    event_matrix,
    file.path(
      output_dir,
      'mutation_event_descendant_matrix_sparse.rds'
    )
  )
  write_physicell_csv(
    combined_events,
    file.path(
      output_dir,
      'mutation_event_descendant_manifest.csv'
    ),
    row.names = FALSE,
    compress = compress_csv
  )
  physicell_log_stage(
    sprintf(
      paste(
        'Mutation-event descendant output complete with %s inherited',
        'cell-event entries.'
      ),
      format(
        length(event_matrix@x),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      )
    ),
    enabled = show_progress
  )
  invisible(list(
    matrix = event_matrix,
    manifest = combined_events
  ))
}

# ---- Output assembly and writers ----

#' Reconstruct the nucleotide sequence encoded by one integration profile
#'
#' Reads the allele encoding position by position: `0` keeps the reference base,
#' `-1` drops it, a whole number `1`-`4` substitutes that base, and a fractional
#' value appends the inserted base after the reference one.
#'
#' @export
#' @param profile_row One integration's encoded alleles, one per position.
#' @param reference Barcode reference bases.
#' @return One character string, whose length varies with the indels applied.
physicell_profile_to_sequence <- function(profile_row, reference){
  bases <- c('A', 'G', 'C', 'T')
  sequence_parts <- vapply(seq_along(profile_row), function(position){
    allele <- profile_row[position]
    if(allele == 0){
      return(reference[position])
    }
    if(allele == -1){
      return('')
    }
    if(allele %% 1 == 0){
      return(bases[as.integer(allele)])
    }
    inserted_code <- as.integer(round(abs(allele %% 1) * 10)) %% 10
    paste0(reference[position], bases[inserted_code])
  }, character(1))
  paste0(sequence_parts, collapse = '')
}

#' Assemble the sparse cell-by-recording-position allele matrix
#'
#' Makes two passes over the terminal profiles: the first counts the nonzero
#' alleles at the model's output positions so the triplet vectors can be sized
#' exactly, the second fills them. Columns run integration by integration as
#' `int_<i>_pos_<position>`.
#'
#' @export
#' @param terminal_profiles Named list of terminal-cell profiles, dense or
#'   sparse; the names become the row names.
#' @param model Prepared recording model supplying `output_positions` and
#'   `num_integrations`.
#' @param show_progress Whether to emit progress messages.
#' @param progress_updates Approximate number of progress messages to emit.
#' @return A sparse `Matrix` of encoded alleles.
#' @note Requires the `Matrix` package.
physicell_sparse_recording_matrix <- function(terminal_profiles,
                                              model,
                                              show_progress = FALSE,
                                              progress_updates = 20L){
  if(!requireNamespace('Matrix', quietly = TRUE)){
    stop('The Matrix package is required for compact barcode output.')
  }
  output_positions <- as.integer(model$output_positions)
  output_position_names <- as.character(output_positions)
  num_cells <- length(terminal_profiles)
  num_integrations <- model$num_integrations
  column_names <- unlist(lapply(
    seq_len(num_integrations),
    function(integration){
      paste0('int_', integration, '_pos_', output_positions)
    }
  ))
  if(num_cells == 0){
    return(Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      x = numeric(),
      dims = c(0L, num_integrations * length(output_positions)),
      dimnames = list(names(terminal_profiles), column_names)
    ))
  }

  report_progress <- new_physicell_progress_reporter(
    total = 2L * num_cells,
    label = 'Barcode sparse matrix assembly',
    enabled = show_progress,
    updates = progress_updates,
    unit = 'cell-passes'
  )
  report_progress(0L)
  entry_counts <- matrix(
    0L,
    nrow = num_cells,
    ncol = num_integrations
  )

  for(cell_index in seq_len(num_cells)){
    profile <- terminal_profiles[[cell_index]]
    for(integration in seq_len(num_integrations)){
      if(inherits(profile, 'physicell_sparse_barcode')){
        integration_values <- profile[[integration]]
        entry_counts[cell_index, integration] <- sum(
          match(
            names(integration_values),
            output_position_names,
            nomatch = 0L
          ) > 0L
        )
      } else{
        entry_counts[cell_index, integration] <- sum(
          profile[integration, output_positions] != 0
        )
      }
    }
    report_progress(cell_index)
  }

  num_entries <- sum(as.double(entry_counts))
  if(num_entries > .Machine$integer.max){
    stop(
      'Compact barcode output exceeds Matrix integer-index capacity.'
    )
  }
  num_entries <- as.integer(num_entries)
  row_indices <- integer(num_entries)
  column_indices <- integer(num_entries)
  values <- numeric(num_entries)
  output_index <- 0L

  for(cell_index in seq_len(num_cells)){
    profile <- terminal_profiles[[cell_index]]
    for(integration in seq_len(num_integrations)){
      num_integration_entries <- entry_counts[cell_index, integration]
      if(num_integration_entries == 0L){
        next
      }
      if(inherits(profile, 'physicell_sparse_barcode')){
        integration_values <- profile[[integration]]
        retained_indices <- match(
          names(integration_values),
          output_position_names,
          nomatch = 0L
        )
        retained <- retained_indices > 0L
        retained_positions <- retained_indices[retained]
        retained_values <- unname(integration_values[retained])
      } else{
        integration_values <- profile[integration, output_positions]
        retained <- integration_values != 0
        retained_positions <- which(retained)
        retained_values <- unname(integration_values[retained])
      }
      output_indices <- seq.int(
        output_index + 1L,
        output_index + num_integration_entries
      )
      row_indices[output_indices] <- cell_index
      column_indices[output_indices] <-
        (integration - 1L) * length(output_positions) +
        retained_positions
      values[output_indices] <- retained_values
      output_index <- output_index + num_integration_entries
    }
    report_progress(num_cells + cell_index)
  }
  if(output_index != num_entries){
    stop('Internal compact barcode entry count changed during assembly.')
  }

  Matrix::sparseMatrix(
    i = row_indices,
    j = column_indices,
    x = values,
    dims = c(
      num_cells,
      num_integrations * length(output_positions)
    ),
    dimnames = list(names(terminal_profiles), column_names)
  )
}

#' One-hot encode PALINCODE states as phylogenetic characters
#'
#' Expands every cBit into three columns (`left`, `right`, `both`); a wild-type
#' cBit is all three zeros, so an edited cBit contributes exactly one nonzero.
#' Columns of `state_matrix` that are not cBit positions are dropped.
#'
#' @export
#' @param state_matrix Cell-by-position PALINCODE state matrix in the column
#'   layout produced by `physicell_sparse_recording_matrix()`; coerced to sparse
#'   when it is not already.
#' @param model Prepared PALINCODE model.
#' @return A sparse indicator matrix with `int_<i>_<cbit>_<state>` columns.
#' @note Errors unless every stored state is `1`, `2`, or `3`.
physicell_palincode_character_matrix <- function(state_matrix, model){
  if(!isTRUE(model$is_palincode)){
    stop('A PALINCODE model is required to encode PALINCODE characters.')
  }
  if(!requireNamespace('Matrix', quietly = TRUE)){
    stop('The Matrix package is required for PALINCODE character output.')
  }
  if(!inherits(state_matrix, 'Matrix')){
    state_matrix <- Matrix::Matrix(state_matrix, sparse = TRUE)
  }
  num_cbits <- model$palincode$num_cbits_per_integration
  num_integrations <- model$num_integrations
  state_names <- c('left', 'right', 'both')
  column_names <- unlist(lapply(seq_len(num_integrations), function(integration){
    unlist(lapply(seq_len(num_cbits), function(cbit_index){
      paste0(
        'int_', integration, '_',
        model$palincode$cbit_names[cbit_index], '_',
        state_names
      )
    }))
  }))
  entries <- methods::as(state_matrix, 'TsparseMatrix')
  if(length(entries@x) == 0L){
    return(Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      x = numeric(),
      dims = c(nrow(state_matrix), length(column_names)),
      dimnames = list(rownames(state_matrix), column_names)
    ))
  }
  entry_rows <- entries@i + 1L
  entry_columns <- entries@j + 1L
  entry_values <- entries@x
  positions_per_integration <- length(model$output_positions)
  position_indices <-
    ((entry_columns - 1L) %% positions_per_integration) + 1L
  positions <- model$output_positions[position_indices]
  integrations <-
    ((entry_columns - 1L) %/% positions_per_integration) + 1L
  cbit_indices <- match(positions, model$palincode$cbit_positions)
  retained <- !is.na(cbit_indices)
  states <- as.integer(round(entry_values[retained]))
  if(any(abs(entry_values[retained] - states) > 1e-8) ||
     any(!(states %in% 1:3))){
    stop('PALINCODE cBit states must use the integer encoding 0, 1, 2, or 3.')
  }
  character_columns <- (
    (integrations[retained] - 1L) * num_cbits +
      cbit_indices[retained] - 1L
  ) * 3L + states
  Matrix::sparseMatrix(
    i = entry_rows[retained],
    j = character_columns,
    x = rep.int(1, sum(retained)),
    dims = c(nrow(state_matrix), length(column_names)),
    dimnames = list(rownames(state_matrix), column_names)
  )
}

#' Describe every integration's target sites as a long table
#'
#' @export
#' @param model Prepared prime-editing, PALINCODE, or generic barcode model.
#' @return One row per integration and target. Prime-editing models add the
#'   static ID, pegRNA assignment, efficiency, and base and effective
#'   probabilities; PALINCODE models add the cBit names, per-state
#'   probabilities, and outcome fractions; generic models report the position,
#'   edit-rate class, and the base-editor reference and destination bases.
physicell_baseline_target_layout <- function(model){
  if(isTRUE(model$is_prime_editing)){
    targets <- model$prime_editing$targets
    return(do.call(rbind, lapply(
      seq_len(model$num_integrations),
      function(integration){
        data.frame(
          integration = integration,
          static_id =
            model$prime_editing$integration_static_ids[integration],
          target_index = targets$target_index,
          position = targets$target_position,
          pegRNA_id = targets$pegRNA_id,
          edit_sequence = targets$edit_sequence,
          editing_efficiency = targets$editing_efficiency,
          induced_base_probability_per_cell_cycle =
            model$prime_editing$induced_base_probabilities,
          induced_effective_probability_per_cell_cycle =
            model$prime_editing$induced_effective_probabilities,
          uninduced_base_probability_per_cell_cycle =
            model$prime_editing$uninduced_base_probabilities,
          uninduced_effective_probability_per_cell_cycle =
            model$prime_editing$uninduced_effective_probabilities,
          spacer_sequence = targets$spacer_sequence,
          pbs_sequence = targets$pbs_sequence,
          rtt_sequence = targets$rtt_sequence,
          description = targets$description,
          stringsAsFactors = FALSE
        )
      }
    )))
  }
  if(isTRUE(model$is_palincode)){
    num_cbits <- model$palincode$num_cbits_per_integration
    return(do.call(rbind, lapply(
      seq_len(model$num_integrations),
      function(integration){
        data.frame(
          integration = integration,
          static_id = model$palincode$integration_static_ids[integration],
          target_index = seq_len(num_cbits),
          target_name = model$palincode$cbit_names,
          position = model$palincode$cbit_positions,
          edit_rate_class = 'PALINCODE',
          reference = 'wild_type',
          alternate = 'left|right|both',
          uninduced_edit_probability_per_cell_cycle =
            model$palincode$uninduced_edit_probabilities,
          induced_edit_probability_per_cell_cycle =
            model$palincode$induced_edit_probabilities,
          left_edit_fraction =
            model$palincode$outcome_fractions[, 'left'],
          right_edit_fraction =
            model$palincode$outcome_fractions[, 'right'],
          both_edit_fraction =
            model$palincode$outcome_fractions[, 'both'],
          stringsAsFactors = FALSE
        )
      }
    )))
  }
  target_positions <- sort(as.integer(names(model$be_targets)))
  if(length(target_positions) == 0){
    return(data.frame(
      integration = integer(),
      target_index = integer(),
      position = integer(),
      edit_rate_class = character(),
      reference = character(),
      alternate = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, lapply(seq_len(model$num_integrations), function(integration){
    data.frame(
      integration = integration,
      target_index = seq_along(target_positions),
      position = target_positions,
      edit_rate_class = unname(
        model$be_targets[as.character(target_positions)]
      ),
      reference = model$be_from,
      alternate = model$be_to,
      stringsAsFactors = FALSE
    )
  }))
}

#' Write the lineage node table, tip metadata, and Newick trees
#'
#' @export
#' @param nodes Event-resolved lineage node table.
#' @param terminal_nodes Sampled terminal rows.
#' @param output_dir Destination directory, created if needed.
#' @param show_progress Whether to emit stage and progress messages.
#' @param progress_updates Approximate number of progress messages to emit.
#' @param compress_csv Whether the CSVs are gzip-compressed.
#' @return Invisibly, the tip metadata table with `sample_id`, `physicell_id`,
#'   `node_id`, `birth_time`, `end_time`, and `branch_length`.
#' @section Side effects: Writes `lineage_nodes.csv`, `terminal_cells.csv`,
#'   `physicell_lineage_full.nwk`, and `physicell_lineage_sampled.nwk` under
#'   `output_dir`.
write_physicell_lineage_outputs <- function(nodes,
                                            terminal_nodes,
                                            output_dir,
                                            show_progress = FALSE,
                                            progress_updates = 20L,
                                            compress_csv = TRUE){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  physicell_log_stage(
    sprintf(
      'Lineage output: writing %s nodes and %s sampled terminal cells.',
      format(nrow(nodes), big.mark = ',', scientific = FALSE, trim = TRUE),
      format(
        nrow(terminal_nodes),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      )
    ),
    enabled = show_progress
  )
  tip_metadata <- terminal_nodes[, c(
    'node_id', 'physicell_id', 'birth_time', 'end_time', 'branch_length'
  ), drop = FALSE]
  tip_metadata$sample_id <- paste0('cell_', tip_metadata$physicell_id)
  tip_metadata <- tip_metadata[, c(
    'sample_id', 'physicell_id', 'node_id', 'birth_time',
    'end_time', 'branch_length'
  )]
  write_physicell_csv(
    nodes,
    file.path(output_dir, 'lineage_nodes.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  write_physicell_csv(
    tip_metadata,
    file.path(output_dir, 'terminal_cells.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  physicell_log_stage(
    'Lineage output: node and terminal-cell tables written; rendering full Newick tree.',
    enabled = show_progress
  )
  write_physicell_lineage_newick(
    nodes,
    file.path(output_dir, 'physicell_lineage_full.nwk'),
    show_progress = show_progress,
    progress_updates = progress_updates,
    progress_label = 'Full Newick rendering'
  )
  physicell_log_stage(
    'Lineage output: full Newick tree written; rendering sampled Newick tree.',
    enabled = show_progress
  )
  write_physicell_lineage_newick(
    nodes,
    file.path(output_dir, 'physicell_lineage_sampled.nwk'),
    terminal_physicell_ids = terminal_nodes$physicell_id,
    show_progress = show_progress,
    progress_updates = progress_updates,
    progress_label = 'Sampled Newick rendering'
  )
  physicell_log_stage(
    'Lineage output complete.',
    enabled = show_progress
  )
  invisible(tip_metadata)
}

#' Write every output of a completed lineage-recording simulation
#'
#' Selects the terminal-cell profiles, then writes either compact sparse RDS
#' matrices or dense CSV matrices according to `model$compact_output`, always
#' emitting an allele matrix and a binary score matrix. PALINCODE models use the
#' one-hot character matrix as that score matrix, and both PALINCODE and
#' prime-editing models duplicate the pair under recorder-specific names.
#'
#' @export
#' @param simulation Result of `simulate_recording_on_physicell_lineage()`.
#' @param model Prepared recording model used for that simulation.
#' @param output_dir Destination directory, created if needed.
#' @param show_progress Whether to emit stage and progress messages.
#' @param write_lineage Whether to also write the shared lineage outputs.
#' @param progress_updates Approximate number of progress messages to emit.
#' @param compress_csv Whether the CSVs are gzip-compressed.
#' @return Invisibly, a list with the normalized `output_dir`,
#'   `terminal_profiles`, `raw_alleles`, and `binary_scores`.
#' @section Side effects: Writes the allele and binary-score matrices,
#'   `mutation_events.csv`, `barcode_target_layout.csv`, `barcode_profiles.rds`,
#'   and `run_manifest.csv` under `output_dir`; a generic barcode model also
#'   writes `barcode_reference.fasta`, and its dense output additionally writes
#'   `barcode_sequences.fasta`.
write_physicell_recording_outputs <- function(simulation,
                                              model,
                                              output_dir,
                                              show_progress = FALSE,
                                              write_lineage = TRUE,
                                              progress_updates = 20L,
                                              compress_csv = TRUE){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if(!dir.exists(output_dir)){
    stop(sprintf('Could not create output directory: %s.', output_dir))
  }

  nodes <- simulation$nodes
  terminal_nodes <- simulation$terminal_nodes
  terminal_profiles <- simulation$profiles[terminal_nodes$node_id]
  names(terminal_profiles) <- paste0('cell_', terminal_nodes$physicell_id)
  physicell_log_stage(
    sprintf(
      'Barcode output started for %s terminal cells and %s mutation events.',
      format(
        nrow(terminal_nodes),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      ),
      format(
        nrow(simulation$mutation_events),
        big.mark = ',',
        scientific = FALSE,
        trim = TRUE
      )
    ),
    enabled = show_progress
  )
  if(isTRUE(write_lineage)){
    write_physicell_lineage_outputs(
      nodes,
      terminal_nodes,
      output_dir,
      show_progress = show_progress,
      progress_updates = progress_updates,
      compress_csv = compress_csv
    )
  }

  if(isTRUE(model$compact_output)){
    physicell_log_stage(
      'Barcode output: assembling compact sparse allele and binary matrices.',
      enabled = show_progress
    )
    raw_alleles <- physicell_sparse_recording_matrix(
      terminal_profiles,
      model,
      show_progress = show_progress,
      progress_updates = progress_updates
    )
    binary_scores <- if(isTRUE(model$is_palincode)){
      physicell_palincode_character_matrix(raw_alleles, model)
    } else{
      scores <- raw_alleles
      if(length(scores@x) > 0){
        scores@x[] <- 1
      }
      scores
    }
    saveRDS(
      raw_alleles,
      file.path(output_dir, 'barcode_alleles_sparse.rds')
    )
    saveRDS(
      binary_scores,
      file.path(output_dir, 'barcode_binary_score_matrix_sparse.rds')
    )
    if(isTRUE(model$is_palincode)){
      saveRDS(
        raw_alleles,
        file.path(output_dir, 'palincode_state_matrix_sparse.rds')
      )
      saveRDS(
        binary_scores,
        file.path(output_dir, 'palincode_character_matrix_sparse.rds')
      )
    } else if(isTRUE(model$is_prime_editing)){
      saveRDS(
        raw_alleles,
        file.path(output_dir, 'prime_editing_state_matrix_sparse.rds')
      )
      saveRDS(
        binary_scores,
        file.path(output_dir, 'prime_editing_character_matrix_sparse.rds')
      )
    }
    physicell_log_stage(
      'Barcode output: compact sparse matrices written.',
      enabled = show_progress
    )
  } else{
    physicell_log_stage(
      'Barcode output: assembling dense allele and binary matrices.',
      enabled = show_progress
    )
    raw_alleles <- do.call(rbind, lapply(terminal_profiles, function(profile){
      as.vector(t(profile))
    }))
    column_names <- unlist(lapply(
      seq_len(model$num_integrations),
      function(integration){
        paste0('int_', integration, '_pos_', seq_len(model$barcode_length))
      }
    ))
    colnames(raw_alleles) <- column_names
    binary_scores <- if(isTRUE(model$is_palincode)){
      as.matrix(physicell_palincode_character_matrix(raw_alleles, model))
    } else{
      (raw_alleles != 0) * 1L
    }
    write_physicell_csv(
      raw_alleles,
      file.path(output_dir, 'barcode_alleles.csv'),
      row.names = TRUE,
      compress = compress_csv
    )
    write_physicell_csv(
      binary_scores,
      file.path(output_dir, 'barcode_binary_score_matrix.csv'),
      row.names = TRUE,
      compress = compress_csv
    )
    if(isTRUE(model$is_palincode)){
      write_physicell_csv(
        raw_alleles,
        file.path(output_dir, 'palincode_state_matrix.csv'),
        row.names = TRUE,
        compress = compress_csv
      )
      write_physicell_csv(
        binary_scores,
        file.path(output_dir, 'palincode_character_matrix.csv'),
        row.names = TRUE,
        compress = compress_csv
      )
    } else if(isTRUE(model$is_prime_editing)){
      write_physicell_csv(
        raw_alleles,
        file.path(output_dir, 'prime_editing_state_matrix.csv'),
        row.names = TRUE,
        compress = compress_csv
      )
      write_physicell_csv(
        binary_scores,
        file.path(output_dir, 'prime_editing_character_matrix.csv'),
        row.names = TRUE,
        compress = compress_csv
      )
    }
    physicell_log_stage(
      'Barcode output: dense matrices written.',
      enabled = show_progress
    )
  }

  physicell_log_stage(
    'Barcode output: writing mutation events and target layout.',
    enabled = show_progress
  )
  write_physicell_csv(
    simulation$mutation_events,
    file.path(output_dir, 'mutation_events.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  target_layout <- physicell_baseline_target_layout(model)
  write_physicell_csv(
    target_layout,
    file.path(output_dir, 'barcode_target_layout.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  if(isTRUE(model$is_prime_editing)){
    write_physicell_csv(
      target_layout,
      file.path(output_dir, 'prime_editing_target_manifest.csv'),
      row.names = FALSE,
      compress = compress_csv
    )
  }
  physicell_log_stage(
    if(isTRUE(model$is_palincode)){
      'Barcode output: saving terminal PALINCODE profiles.'
    } else if(isTRUE(model$is_prime_editing)){
      'Barcode output: saving terminal prime-editing profiles.'
    } else{
      'Barcode output: saving terminal profiles and barcode reference.'
    },
    enabled = show_progress
  )
  saveRDS(
    terminal_profiles,
    file.path(output_dir, 'barcode_profiles.rds')
  )

  if(!isTRUE(model$is_palincode) && !isTRUE(model$is_prime_editing)){
    writeLines(
      c(
        '>remote_mito_barcode_reference',
        paste0(model$barcode_sequence, collapse = '')
      ),
      file.path(output_dir, 'barcode_reference.fasta')
    )
  }

  if(!isTRUE(model$compact_output) &&
     !isTRUE(model$is_palincode) &&
     !isTRUE(model$is_prime_editing)){
    physicell_log_stage(
      'Barcode output: writing reconstructed barcode FASTA.',
      enabled = show_progress
    )
    fasta_lines <- character()
    for(cell_name in names(terminal_profiles)){
      profile <- terminal_profiles[[cell_name]]
      for(integration in seq_len(nrow(profile))){
        fasta_lines <- c(
          fasta_lines,
          paste0('>', cell_name, '_int_', integration),
          physicell_profile_to_sequence(
            profile[integration, ],
            model$barcode_sequence
          )
        )
      }
    }
    writeLines(fasta_lines, file.path(output_dir, 'barcode_sequences.fasta'))
  }

  physicell_log_stage(
    'Barcode output: writing run manifest.',
    enabled = show_progress
  )
  manifest <- data.frame(
    property = c(
      'num_division_events',
      'num_lineage_nodes',
      'num_full_terminal_cells',
      'num_sampled_terminal_cells',
      'barcode_length',
      'num_integrations',
      'num_targets_per_integration',
      'total_integrated_target_sites',
      'num_recording_events',
      'recorder_system',
      'profile_storage',
      'output_format',
      'founder_label_sites',
      'cell_type',
      'editing_state',
      'random_seed'
    ),
    value = c(
      sum(nodes$origin == 'continuing_parent'),
      nrow(nodes),
      sum(nodes$is_terminal),
      nrow(terminal_nodes),
      model$barcode_length,
      model$num_integrations,
      if(isTRUE(model$is_palincode)){
        model$palincode$num_cbits_per_integration
      } else if(isTRUE(model$is_prime_editing)){
        nrow(model$prime_editing$targets)
      } else{
        length(model$be_targets)
      },
      model$num_integrations * if(isTRUE(model$is_palincode)){
        model$palincode$num_cbits_per_integration
      } else if(isTRUE(model$is_prime_editing)){
        nrow(model$prime_editing$targets)
      } else{
        length(model$be_targets)
      },
      nrow(simulation$mutation_events),
      model$recorder_system,
      model$profile_storage,
      if(isTRUE(model$compact_output)){
        'sparse_rds'
      } else if(isTRUE(compress_csv)){
        'dense_csv_gz'
      } else{
        'dense_csv'
      },
      model$founder_label_sites,
      model$cell_type,
      simulation$editing_state,
      simulation$seed
    ),
    stringsAsFactors = FALSE
  )
  if(isTRUE(model$is_palincode)){
    manifest <- rbind(
      manifest,
      data.frame(
        property = c(
          'palincode_static_id_length',
          'palincode_state_encoding',
          'palincode_character_encoding'
        ),
        value = c(
          model$palincode$static_id_length,
          '0=wild_type;1=left;2=right;3=both',
          'one-hot edited outcomes: left,right,both; wild type is all zero'
        ),
        stringsAsFactors = FALSE
      )
    )
  } else if(isTRUE(model$is_prime_editing)){
    manifest <- rbind(
      manifest,
      data.frame(
        property = c(
          'prime_editing_pool_source',
          'prime_editing_assignment',
          'prime_editing_pool_size',
          'prime_editing_static_id_length',
          'prime_editing_state_encoding'
        ),
        value = c(
          model$prime_editing$pool_source,
          model$prime_editing$assignment,
          nrow(model$prime_editing$pool),
          model$prime_editing$static_id_length,
          '0=unedited;1=assigned pegRNA edit'
        ),
        stringsAsFactors = FALSE
      )
    )
  }
  write_physicell_csv(
    manifest,
    file.path(output_dir, 'run_manifest.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  physicell_log_stage(
    'Barcode output complete.',
    enabled = show_progress
  )

  invisible(list(
    output_dir = normalizePath(output_dir, mustWork = TRUE),
    terminal_profiles = terminal_profiles,
    raw_alleles = raw_alleles,
    binary_scores = binary_scores
  ))
}
