# Shared prime-editing recorder backend.
#
# A prepared backend maps each recorder target to one known pegRNA and keeps
# the intended edit sequence in exact character form. Simulation engines store
# only the categorical edited/unedited target state, avoiding lossy numeric
# encodings of long prime-editing templates.

#' Report whether prime editing is switched on for a parameter set
#'
#' The unified `prime_editing_backend` block wins when it is present: prime
#' editing is on unless its `enabled` field is exactly `FALSE`. With no such
#' block the legacy `nuclease_targets$prime_editing_system` flag decides.
#'
#' @export
#' @param params Parsed remote_mito JSON parameter list. A
#'   `prime_editing_backend` entry that is not a JSON object is an error.
#' @return A single logical.
prime_editing_enabled <- function(params){
  configuration <- params$prime_editing_backend
  if(!is.null(configuration) && !is.list(configuration)){
    stop('prime_editing_backend must be a JSON object.')
  }
  if(!is.null(configuration)){
    return(!identical(configuration$enabled, FALSE))
  }
  isTRUE(params$nuclease_targets$prime_editing_system)
}

#' Validate and recycle a per-target prime-editing probability vector
#'
#' Values are per-cell-cycle edit probabilities, not continuous-time hazards. A
#' single scalar is recycled to one value per recorder target; any other length
#' must already equal `number`.
#'
#' @export
#' @param value Scalar or vector of probabilities; flattened with `unlist()`
#'   and coerced with `as.numeric`, so a JSON array is accepted.
#' @param number Number of prime-editing targets the vector must cover.
#' @param name Diagnostic name used in the error message.
#' @return A numeric vector of length `number` whose entries are all finite and
#'   in `[0, 1]`; stops otherwise.
prime_editing_probability_vector <- function(value, number, name){
  value <- suppressWarnings(as.numeric(unlist(value, use.names = FALSE)))
  if(length(value) == 1L){
    value <- rep(value, number)
  }
  if(length(value) != number || any(!is.finite(value)) ||
     any(value < 0) || any(value > 1)){
    stop(sprintf(
      '%s must contain one probability or one per prime-editing target.',
      name
    ))
  }
  value
}

#' Validate one A/C/G/T sequence string
#'
#' Upper-cases the value and strips every whitespace character before checking
#' it, so sequences may be entered with spaces or line breaks.
#'
#' @export
#' @param value Candidate sequence; coerced with `as.character`.
#' @param name Diagnostic name used in the error message.
#' @param allow_empty When `TRUE` an empty string also passes validation.
#' @return The cleaned uppercase sequence; stops unless it is exactly one
#'   non-`NA` string over `A`, `C`, `G`, `T` (or empty when `allow_empty`).
prime_editing_dna_sequence <- function(value, name, allow_empty = FALSE){
  value <- toupper(gsub('[[:space:]]+', '', as.character(value)))
  if(length(value) != 1L || is.na(value) ||
     (!allow_empty && !nzchar(value)) ||
     (nzchar(value) && !grepl('^[ACGT]+$', value))){
    stop(sprintf('%s must be one %sA/C/G/T sequence.',
                 name, if(allow_empty) 'possibly empty ' else 'non-empty '))
  }
  value
}

#' Rectangularize a JSON-style pegRNA array into a data frame
#'
#' Takes the union of the field names seen across all pegRNA objects and fills
#' the fields an individual object omits with `NA`, so objects with different
#' key sets can be row-bound. A value that is already a data frame is returned
#' unchanged.
#'
#' @export
#' @param value A data frame, or a non-empty list whose elements are all lists,
#'   one per pegRNA.
#' @return A data frame with one row per pegRNA and one column per field name
#'   observed anywhere in the array, in first-seen order.
prime_editing_pool_frame <- function(value){
  if(is.data.frame(value)){
    return(value)
  }
  if(!is.list(value) || length(value) == 0L ||
     !all(vapply(value, is.list, logical(1)))){
    stop('prime_editing_backend.pegRNAs must be a non-empty array of objects.')
  }
  column_names <- unique(unlist(lapply(value, names), use.names = FALSE))
  rows <- lapply(value, function(row){
    missing <- setdiff(column_names, names(row))
    for(name in missing){
      row[[name]] <- NA
    }
    as.data.frame(row[column_names], stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Load and validate the known pegRNA pool
#'
#' Resolves the pool from the first source that is configured:
#' `prime_editing_backend$pegRNAs` (embedded JSON objects),
#' `prime_editing_backend$pegRNA_pool_path` (a CSV, resolved relative to
#' `params_dir` unless the path is absolute), or the legacy
#' `nuclease_targets` guide-count and guide-length settings, which synthesize a
#' pool of uniformly random sequences with efficiency `1`.
#'
#' @export
#' @param params Parsed remote_mito JSON parameter list.
#' @param params_dir Directory that a relative `pegRNA_pool_path` is resolved
#'   against.
#' @return A data frame with one row per pegRNA and columns `pool_index`
#'   (row order, `1..n`), `pegRNA_id` (unique, non-empty, trimmed),
#'   `edit_sequence` (validated A/C/G/T), `editing_efficiency` (finite, in
#'   `[0, 1]`), `spacer_sequence`, `pbs_sequence`, `rtt_sequence` (each
#'   validated A/C/G/T, or `NA_character_` when absent or blank), and
#'   `description`. Its `source` attribute records provenance:
#'   `'embedded_json'`, the normalized CSV path, or `'generated_legacy_pool'`.
#' @section Side effects: The legacy fallback draws random sequences, so it
#'   consumes the R random-number stream.
read_prime_editing_pool <- function(params, params_dir = '.'){
  configuration <- params$prime_editing_backend
  if(is.null(configuration)){
    configuration <- list()
  }
  pool <- NULL
  source_label <- NULL
  if(!is.null(configuration$pegRNAs)){
    pool <- prime_editing_pool_frame(configuration$pegRNAs)
    source_label <- 'embedded_json'
  } else if(!is.null(configuration$pegRNA_pool_path)){
    path <- as.character(configuration$pegRNA_pool_path)
    if(length(path) != 1L || is.na(path) || !nzchar(path)){
      stop('prime_editing_backend.pegRNA_pool_path must be one path.')
    }
    if(!grepl('^(/|[A-Za-z]:[/\\\\])', path)){
      path <- file.path(params_dir, path)
    }
    path <- normalizePath(path, mustWork = TRUE)
    pool <- utils::read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    source_label <- path
  } else{
    legacy_count <- suppressWarnings(as.integer(
      params$nuclease_targets$num_unique_prime_editing_guides
    ))
    legacy_length <- suppressWarnings(as.integer(
      params$nuclease_targets$prime_editing_guide_length
    ))
    if(length(legacy_count) != 1L || is.na(legacy_count) || legacy_count < 1L ||
       length(legacy_length) != 1L || is.na(legacy_length) || legacy_length < 1L){
      stop(paste(
        'Prime editing requires prime_editing_backend.pegRNAs,',
        'pegRNA_pool_path, or valid legacy guide-count/length settings.'
      ))
    }
    edit_sequences <- vapply(seq_len(legacy_count), function(index){
      paste0(sample(c('A', 'G', 'C', 'T'), legacy_length, TRUE), collapse = '')
    }, character(1))
    pool <- data.frame(
      pegRNA_id = paste0('legacy_pegRNA_', seq_len(legacy_count)),
      edit_sequence = edit_sequences,
      editing_efficiency = 1,
      stringsAsFactors = FALSE
    )
    source_label <- 'generated_legacy_pool'
  }

  required <- c('pegRNA_id', 'edit_sequence', 'editing_efficiency')
  missing <- setdiff(required, names(pool))
  if(length(missing) > 0L){
    stop(sprintf(
      'The pegRNA pool is missing required column(s): %s.',
      paste(missing, collapse = ', ')
    ))
  }
  if(nrow(pool) == 0L){
    stop('The pegRNA pool must contain at least one row.')
  }
  pool$pegRNA_id <- trimws(as.character(pool$pegRNA_id))
  if(anyNA(pool$pegRNA_id) || any(!nzchar(pool$pegRNA_id)) ||
     anyDuplicated(pool$pegRNA_id)){
    stop('pegRNA_id values must be unique, non-empty strings.')
  }
  pool$edit_sequence <- vapply(seq_len(nrow(pool)), function(index){
    prime_editing_dna_sequence(
      pool$edit_sequence[index],
      sprintf('edit_sequence for pegRNA %s', pool$pegRNA_id[index])
    )
  }, character(1))
  pool$editing_efficiency <- suppressWarnings(as.numeric(
    pool$editing_efficiency
  ))
  if(any(!is.finite(pool$editing_efficiency)) ||
     any(pool$editing_efficiency < 0) ||
     any(pool$editing_efficiency > 1)){
    stop('editing_efficiency values must be probabilities in [0, 1].')
  }
  optional_sequences <- c(
    'spacer_sequence', 'pbs_sequence', 'rtt_sequence'
  )
  for(column in optional_sequences){
    if(!(column %in% names(pool))){
      pool[[column]] <- NA_character_
      next
    }
    pool[[column]] <- vapply(seq_len(nrow(pool)), function(index){
      value <- pool[[column]][index]
      if(is.na(value) || !nzchar(trimws(as.character(value)))){
        return(NA_character_)
      }
      prime_editing_dna_sequence(
        value,
        sprintf('%s for pegRNA %s', column, pool$pegRNA_id[index])
      )
    }, character(1))
  }
  if(!('description' %in% names(pool))){
    pool$description <- NA_character_
  } else{
    pool$description <- as.character(pool$description)
  }
  pool$pool_index <- seq_len(nrow(pool))
  pool <- pool[, c(
    'pool_index', 'pegRNA_id', 'edit_sequence', 'editing_efficiency',
    optional_sequences, 'description'
  ), drop = FALSE]
  attr(pool, 'source') <- source_label
  pool
}

#' Assign one pegRNA from the pool to each recorder target
#'
#' `configuration$target_pegRNA_ids` assigns explicitly and must name exactly
#' one existing `pegRNA_id` per target, in target order. Otherwise
#' `configuration$assignment` selects `cycle` (recycle the pool rows in order,
#' the default), `sample_with_replacement`, or `sample_without_replacement`
#' (which requires at least as many pegRNAs as targets).
#'
#' @export
#' @param pool Validated pegRNA pool from `read_prime_editing_pool()`.
#' @param target_positions Recorder coordinates; must be unique positive
#'   integers, and their order fixes `target_index`.
#' @param configuration The `prime_editing_backend` JSON object, or an empty
#'   list to take the defaults.
#' @return A data frame with one row per target and columns `target_index`
#'   (`1..k`), `target_position`, `pool_index`, `pegRNA_id`, `edit_sequence`,
#'   `editing_efficiency`, `spacer_sequence`, `pbs_sequence`, `rtt_sequence`,
#'   `description`.
#' @section Side effects: The two sampling modes consume the random-number
#'   stream.
assign_prime_editing_targets <- function(pool,
                                         target_positions,
                                         configuration = list()){
  target_positions <- suppressWarnings(as.integer(target_positions))
  if(length(target_positions) == 0L || anyNA(target_positions) ||
     any(target_positions < 1L) || anyDuplicated(target_positions)){
    stop('Prime-editing target positions must be unique positive integers.')
  }
  configured_ids <- configuration$target_pegRNA_ids
  if(!is.null(configured_ids)){
    assigned_ids <- as.character(unlist(configured_ids, use.names = FALSE))
    if(length(assigned_ids) != length(target_positions)){
      stop(paste(
        'prime_editing_backend.target_pegRNA_ids must contain one ID per',
        'prime-editing target.'
      ))
    }
    pool_indices <- match(assigned_ids, pool$pegRNA_id)
    if(anyNA(pool_indices)){
      stop(sprintf(
        'Unknown target_pegRNA_ids value(s): %s.',
        paste(unique(assigned_ids[is.na(pool_indices)]), collapse = ', ')
      ))
    }
  } else{
    assignment <- if(is.null(configuration$assignment)){
      'cycle'
    } else{
      tolower(as.character(configuration$assignment))
    }
    if(length(assignment) != 1L || !(assignment %in% c(
      'cycle', 'sample_with_replacement', 'sample_without_replacement'
    ))){
      stop(paste(
        'prime_editing_backend.assignment must be cycle,',
        'sample_with_replacement, or sample_without_replacement.'
      ))
    }
    pool_indices <- switch(
      assignment,
      cycle = rep(seq_len(nrow(pool)), length.out = length(target_positions)),
      sample_with_replacement = sample(
        seq_len(nrow(pool)), length(target_positions), replace = TRUE
      ),
      sample_without_replacement = {
        if(length(target_positions) > nrow(pool)){
          stop('Not enough pegRNAs for sample_without_replacement assignment.')
        }
        sample(seq_len(nrow(pool)), length(target_positions), replace = FALSE)
      }
    )
  }
  assigned <- pool[pool_indices, , drop = FALSE]
  assigned$target_position <- target_positions
  assigned$target_index <- seq_along(target_positions)
  assigned[, c(
    'target_index', 'target_position', 'pool_index', 'pegRNA_id',
    'edit_sequence', 'editing_efficiency', 'spacer_sequence',
    'pbs_sequence', 'rtt_sequence', 'description'
  ), drop = FALSE]
}

#' Expand each recorder target into its alphabet of installable marks
#'
#' PEtracer-style recorders place several pegRNAs on the same edit site, each
#' installing a different predefined mark, so an edit writes one of `N` outcomes
#' rather than a single fixed one. That distinction matters for reconstruction:
#' with one outcome per site, two cells that edit the same site independently
#' are indistinguishable from two cells sharing an ancestor, and the resulting
#' homoplasy misleads tree building. With `N` marks they differ `(N - 1) / N` of
#' the time.
#'
#' Marks are taken from the pool in blocks, so target one receives the first
#' `marks_per_target` entries, target two the next block, and so on, cycling
#' when the pool is exhausted.
#'
#' @export
#' @param pool Validated pegRNA pool.
#' @param target_positions Recorder coordinates.
#' @param marks_per_target Marks installable at each site; `1` reproduces the
#'   single-outcome behaviour exactly.
#' @return A data frame with one row per target and mark, carrying
#'   `target_index`, `target_position`, `mark_index`, `allele`, `pegRNA_id`,
#'   `edit_sequence`, `editing_efficiency`, and `mark_probability`, the
#'   efficiency-weighted chance that this mark is the one installed when the
#'   site edits.
expand_prime_editing_marks <- function(pool, target_positions,
                                       marks_per_target = 1L){
  marks_per_target <- suppressWarnings(as.integer(marks_per_target))
  if(length(marks_per_target) != 1L || is.na(marks_per_target) ||
     marks_per_target < 1L){
    stop('prime_editing_backend.marks_per_target must be a positive integer.')
  }
  target_positions <- as.integer(target_positions)
  if(marks_per_target > nrow(pool)){
    # Cycling would hand the same pegRNA to one site twice, which silently
    # halves that site's real alphabet rather than failing.
    stop(sprintf(paste(
      'prime_editing_backend.marks_per_target is %d but the pegRNA pool has',
      'only %d entr%s; a site cannot install the same mark twice.'
    ), marks_per_target, nrow(pool), if(nrow(pool) == 1L) 'y' else 'ies'))
  }
  blocks <- lapply(seq_along(target_positions), function(index){
    offset <- (index - 1L) * marks_per_target
    rows <- ((offset + seq_len(marks_per_target) - 1L) %% nrow(pool)) + 1L
    block <- pool[rows, , drop = FALSE]
    efficiency <- as.numeric(block$editing_efficiency)
    # Efficiency decides which mark wins once the site fires; it does not also
    # decide whether it fires, which the site's own rate governs. An all-zero
    # block would otherwise divide by zero, so it falls back to uniform.
    weights <- if(sum(efficiency) > 0) efficiency / sum(efficiency) else
      rep(1 / length(efficiency), length(efficiency))
    data.frame(
      target_index = index,
      target_position = target_positions[index],
      mark_index = seq_len(marks_per_target),
      # Allele 0 is "unedited", so marks start at 1.
      allele = seq_len(marks_per_target),
      pegRNA_id = as.character(block$pegRNA_id),
      edit_sequence = as.character(block$edit_sequence),
      editing_efficiency = efficiency,
      mark_probability = weights,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, blocks)
}

#' Prepare the shared prime-editing backend
#'
#' Loads the pegRNA pool, assigns one pegRNA per recorder target, and records
#' the provenance an engine needs to reproduce the mapping. Returns `NULL` when
#' prime editing is disabled, which callers treat as the "not prime editing"
#' signal.
#'
#' @export
#' @param params Parsed remote_mito JSON parameter list.
#' @param target_positions Recorder coordinates to assign pegRNAs to.
#' @param params_dir Directory that a relative `pegRNA_pool_path` is resolved
#'   against.
#' @param seed Optional seed set before the pool is read and the targets are
#'   assigned, making the legacy pool and the sampled assignment modes
#'   reproducible.
#' @return `NULL` when prime editing is off, otherwise a named list with `pool`
#'   (the validated pool), `targets` (the per-target assignment table),
#'   `pool_source` (the pool's `source` attribute), `assignment` (`'explicit'`,
#'   `'cycle'`, `'sample_with_replacement'`, or
#'   `'sample_without_replacement'`), and `state_encoding`, the fixed
#'   categorical target state `c(unedited = 0L, edited = 1L)` that engines
#'   store in place of a numeric copy of the edit template.
#' @section Side effects: Calls `set.seed()` when `seed` is supplied.
prepare_prime_editing_backend <- function(params,
                                          target_positions,
                                          params_dir = '.',
                                          seed = NULL){
  if(!prime_editing_enabled(params)){
    return(NULL)
  }
  if(!is.null(seed)){
    set.seed(seed)
  }
  configuration <- params$prime_editing_backend
  if(is.null(configuration)){
    configuration <- list()
  }
  pool <- read_prime_editing_pool(params, params_dir = params_dir)
  targets <- assign_prime_editing_targets(
    pool,
    target_positions,
    configuration
  )
  marks_per_target <- if(is.null(configuration$marks_per_target)) 1L else
    configuration$marks_per_target
  marks <- expand_prime_editing_marks(
    pool, targets$target_position, marks_per_target
  )
  marks_per_target <- max(marks$mark_index)
  if(marks_per_target > 1L){
    # assign_prime_editing_targets() walks the pool one entry per site, which
    # for a blocked mark pool hands every site the FIRST block's pegRNA and so
    # the wrong efficiency. A site's rate has to come from its own alphabet:
    # the mean over the marks it can actually install.
    block_means <- vapply(
      split(marks$editing_efficiency, marks$target_position),
      mean, numeric(1)
    )
    block_first <- marks[!duplicated(marks$target_position), , drop = FALSE]
    order_by_target <- match(targets$target_position, block_first$target_position)
    targets$editing_efficiency <- unname(
      block_means[as.character(targets$target_position)]
    )
    targets$pegRNA_id <- block_first$pegRNA_id[order_by_target]
    targets$edit_sequence <- block_first$edit_sequence[order_by_target]
  }
  # With one mark the recorder is binary, as before. With several, an edited
  # site carries which mark it received, so the encoding has to name them.
  state_encoding <- if(marks_per_target == 1L){
    c(unedited = 0L, edited = 1L)
  } else{
    setNames(
      c(0L, seq_len(marks_per_target)),
      c('unedited', sprintf('mark_%d', seq_len(marks_per_target)))
    )
  }
  list(
    pool = pool,
    targets = targets,
    marks = marks,
    marks_per_target = marks_per_target,
    pool_source = attr(pool, 'source'),
    assignment = if(is.null(configuration$target_pegRNA_ids)){
      if(is.null(configuration$assignment)) 'cycle' else
        tolower(as.character(configuration$assignment))
    } else{
      'explicit'
    },
    state_encoding = state_encoding
  )
}

#' Scale a per-cell-cycle edit probability by pegRNA efficiency
#'
#' Both the input and the result are per-cell-cycle probabilities, never
#' continuous-time hazards; callers convert the result with
#' `physicell_probability_hazard()` when they need a rate. Efficiency `e` is
#' applied as `1 - (1 - p)^e`, so unedited survival probabilities stay
#' multiplicative: efficiency `0` leaves a target unedited and efficiency `1`
#' returns `p` unchanged. Evaluated through `log1p`/`expm1` for accuracy at
#' small `p`, after clamping `p` just below one.
#'
#' @export
#' @param base_probability Per-cell-cycle edit probability, or one per target;
#'   every entry must be finite and in `[0, 1]`.
#' @param efficiency pegRNA `editing_efficiency`, recycled against
#'   `base_probability` by ordinary vector arithmetic; every entry must be
#'   finite and in `[0, 1]`.
#' @return A numeric vector of efficiency-adjusted per-cell-cycle
#'   probabilities.
prime_editing_scale_probability <- function(base_probability, efficiency){
  base_probability <- suppressWarnings(as.numeric(base_probability))
  efficiency <- suppressWarnings(as.numeric(efficiency))
  if(any(!is.finite(base_probability)) || any(base_probability < 0) ||
     any(base_probability > 1)){
    stop('Prime-editing base probabilities must be in [0, 1].')
  }
  if(any(!is.finite(efficiency)) || any(efficiency < 0) ||
     any(efficiency > 1)){
    stop('Prime-editing efficiencies must be in [0, 1].')
  }
  base_probability <- pmin(base_probability, 1 - .Machine$double.eps)
  -expm1(log1p(-base_probability) * efficiency)
}

#' Build the legacy position-to-integer edit-template map
#'
#' Encodes each assigned `edit_sequence` base by base using the repo's
#' nucleotide codes `A = 1`, `G = 2`, `C = 3`, `T = 4`.
#'
#' @export
#' @param backend Prepared backend from `prepare_prime_editing_backend()`.
#' @return A list with one integer vector per target, named by the target's
#'   `target_position` rendered as a character string.
#' @note Kept for the legacy time-step engine, which stores integer templates;
#'   `prime_editing_sequence_character_map()` is the exact, non-lossy form.
prime_editing_sequence_integer_map <- function(backend){
  base_codes <- c(A = 1L, G = 2L, C = 3L, T = 4L)
  sequences <- lapply(backend$targets$edit_sequence, function(sequence){
    unname(base_codes[strsplit(sequence, '', fixed = TRUE)[[1]]])
  })
  names(sequences) <- as.character(backend$targets$target_position)
  sequences
}

#' Build the exact position-to-edit-sequence map
#'
#' Keeps the intended edit in character form, which is what engines use to
#' write the realized edit for a target.
#'
#' @export
#' @param backend Prepared backend from `prepare_prime_editing_backend()`.
#' @return A list with one `edit_sequence` string per target, named by the
#'   target's `target_position` rendered as a character string.
prime_editing_sequence_character_map <- function(backend){
  sequences <- as.list(backend$targets$edit_sequence)
  names(sequences) <- as.character(backend$targets$target_position)
  sequences
}

#' Write the target/pegRNA assignment manifest
#'
#' Records the exact mapping a run used, so a replay can be tied back to the
#' pool it came from.
#'
#' @export
#' @param backend Prepared backend from `prepare_prime_editing_backend()`.
#' @param path Destination passed straight through to `write_csv`.
#' @param write_csv Writer invoked as `write_csv(manifest, path, row.names =
#'   FALSE)`; defaults to `utils::write.csv` and can be swapped for a
#'   compressing writer.
#' @return Invisibly, the manifest data frame: the `targets` table with
#'   `pool_source` and `assignment` appended as constant columns.
#' @section Side effects: Writes the manifest through `write_csv` to `path`.
write_prime_editing_backend_manifest <- function(backend,
                                                 path,
                                                 write_csv = utils::write.csv){
  manifest <- backend$targets
  manifest$pool_source <- backend$pool_source
  manifest$assignment <- backend$assignment
  write_csv(manifest, path, row.names = FALSE)
  invisible(manifest)
}
