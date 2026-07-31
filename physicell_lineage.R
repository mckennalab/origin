# PhysiCell lineage import and barcode-recording simulation.
#
# This file intentionally contains functions only so it can be sourced from
# tests, notebooks, or the simulate_physicell_lineage.R command-line wrapper.

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

read_physicell_terminal_ids <- function(path){
  read_physicell_cell_ids(
    path,
    table_description = 'PhysiCell live-cell table',
    alive_only = TRUE
  )
}

read_physicell_founder_ids <- function(path){
  read_physicell_cell_ids(
    path,
    table_description = 'PhysiCell founder table',
    alive_only = FALSE
  )
}

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

physicell_lineage_to_newick <- function(nodes, terminal_physicell_ids = NULL){
  prepared <- prepare_physicell_newick(nodes, terminal_physicell_ids)
  max_tokens <- 4L * prepared$num_kept_nodes +
    2L * length(prepared$root_indices) + 2L
  tokens <- character(max_tokens)
  token_count <- 0L
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

physicell_substitution_matrix <- function(model_name, model_parameters, sequence){
  bases <- c('A', 'G', 'C', 'T')
  fractions <- as.numeric(table(factor(sequence, levels = bases))) / length(sequence)
  model_name <- toupper(as.character(model_name))

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

prepare_physicell_recording_model <- function(params,
                                              cell_type = NULL,
                                              num_integrations = NULL,
                                              founder_label_sites = 0,
                                              params_dir = '.',
                                              seed = NULL){
  if(!is.list(params)){
    stop('params must be the parsed remote_mito JSON parameter list.')
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

physicell_probability_hazard <- function(probability, cell_cycle_length){
  probability_names <- names(probability)
  probability <- pmin(pmax(as.numeric(probability), 0), 1 - .Machine$double.eps)
  hazards <- -log1p(-probability) / cell_cycle_length
  names(hazards) <- probability_names
  hazards
}

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

mutate_physicell_barcode_segment <- function(profile,
                                             duration,
                                             rate_set,
                                             model,
                                             segment_start = 0){
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
    label = if(is.null(model$recorder_system)){
      'barcode lineage recording'
    } else{
      model$recorder_system
    },
    enabled = show_progress,
    updates = progress_updates
  )
  report_progress(0L)

  for(node_index in seq_len(nrow(nodes))){
    node <- nodes[node_index, , drop = FALSE]
    if(is.na(node$parent_node_id)){
      founder_number <- founder_number + 1L
      initialized <- initialize_physicell_founder_barcode(
        model,
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
      if(!isTRUE(model$retain_internal_profiles) &&
         remaining_children[parent_node_id] <= 0){
        profiles[[parent_node_id]] <- NULL
      }
    }

    mutation <- mutate_physicell_barcode_branch(
      start_profile,
      node$birth_time,
      node$end_time,
      model,
      editing_state
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
      event_frame[[column]] <- NA
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

physicell_baseline_target_layout <- function(model){
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
    binary_scores <- raw_alleles
    if(length(binary_scores@x) > 0){
      binary_scores@x[] <- 1
    }
    saveRDS(
      raw_alleles,
      file.path(output_dir, 'barcode_alleles_sparse.rds')
    )
    saveRDS(
      binary_scores,
      file.path(output_dir, 'barcode_binary_score_matrix_sparse.rds')
    )
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
    binary_scores <- (raw_alleles != 0) * 1L
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
  write_physicell_csv(
    physicell_baseline_target_layout(model),
    file.path(output_dir, 'barcode_target_layout.csv'),
    row.names = FALSE,
    compress = compress_csv
  )
  physicell_log_stage(
    'Barcode output: saving terminal profiles and barcode reference.',
    enabled = show_progress
  )
  saveRDS(
    terminal_profiles,
    file.path(output_dir, 'barcode_profiles.rds')
  )

  writeLines(
    c('>remote_mito_barcode_reference', paste0(model$barcode_sequence, collapse = '')),
    file.path(output_dir, 'barcode_reference.fasta')
  )

  if(!isTRUE(model$compact_output)){
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
      length(model$be_targets),
      model$num_integrations * length(model$be_targets),
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
