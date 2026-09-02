# Turn per-cell mutation matrices ("profiles") into score matrices for
# phylogenetic reconstruction, and write them out as RDS, FASTA and PHYLIP.
# Rows of a score matrix are cells keyed by lineage string; columns are unique
# mutation events (position and mutation value, plus integration for barcodes).
# Sourced by sim5_code.R after data.table, Matrix and dplyr are attached;
# get_norm_cell_heteroplasmy_scores() additionally needs draw_severity_scores().

#' Collapse runs of consecutive deleted positions into one row each
#'
#' Sorts the deletion events by lineage string, integration and position, then
#' walks them once and closes a run whenever the next event belongs to another
#' cell or integration or does not sit at the following position. Each run is
#' reported once, at its last position, with the run length as its mutation
#' value.
#'
#' @param deletion_df Data frame of deletion events with columns `linstring`,
#'   `ints_mutated`, `positions_mutated` and a mutation-value column; only
#'   deletion rows (mutation value `-1`) are expected. Need not be sorted.
#' @return A character matrix carrying the column names of `deletion_df`, one
#'   row per contiguous run: lineage string, integration, the last position of
#'   the run, and `d<run length>` in place of the mutation value. A zero-row
#'   character matrix with the same column names when `deletion_df` is empty.
#' @note Run continuation is tested with `positions_mutated[i + 1] ==
#'   positions_mutated[i] + 1`, so the position column must be numeric; a
#'   character position column errors on that arithmetic.
group_deletions <- function(deletion_df){
  if(nrow(deletion_df) == 0){
    return(matrix(character(), nrow = 0, ncol = ncol(deletion_df),
                  dimnames = list(NULL, colnames(deletion_df))))
  }

  deletion_df <- deletion_df[
    order(deletion_df$linstring,
          deletion_df$ints_mutated,
          deletion_df$positions_mutated),
    ,
    drop = FALSE
  ]

  grouped_rows <- list()
  streak_start <- 1L
  output_index <- 1L

  for(row_index in seq_len(nrow(deletion_df))){
    is_last <- row_index == nrow(deletion_df)
    continues <- !is_last &&
      deletion_df$linstring[row_index + 1L] == deletion_df$linstring[row_index] &&
      deletion_df$ints_mutated[row_index + 1L] == deletion_df$ints_mutated[row_index] &&
      deletion_df$positions_mutated[row_index + 1L] ==
        deletion_df$positions_mutated[row_index] + 1

    if(!continues){
      streak_length <- row_index - streak_start + 1L
      grouped_rows[[output_index]] <- c(
        as.character(deletion_df$linstring[row_index]),
        as.character(deletion_df$ints_mutated[row_index]),
        as.character(deletion_df$positions_mutated[row_index]),
        paste0('d', streak_length)
      )
      output_index <- output_index + 1L
      streak_start <- row_index + 1L
    }
  }

  grouped <- do.call(rbind, grouped_rows)
  colnames(grouped) <- colnames(deletion_df)
  grouped
}

#' Write a score matrix as sequential PHYLIP
#'
#' Emits a `<taxa> <characters>` header line followed by one line per row: the
#' row name, a single space, and the row's scores concatenated with no
#' separator.
#'
#' @param score_mat Matrix whose rows are taxa (cells, named by `rownames`) and
#'   whose columns are characters (unique mutations). Each score is pasted as
#'   printed, so one score occupies one column of the alignment only when every
#'   score is a single digit.
#' @param output_phylip_path Path of the file to write.
#' @return `NULL`; called for its side effect.
#' @section Side effects: Creates or overwrites `output_phylip_path`.
#' @note Taxon names are separated from the data by one space and are not
#'   padded to the ten-character field of strict PHYLIP, so the output is the
#'   relaxed sequential form.
score_mat_to_phylip <- function(score_mat, output_phylip_path) {
  num_taxa <- nrow(score_mat)
  num_chars <- ncol(score_mat)
  
  phylip_f <- file(output_phylip_path, "w")
  
  writeLines(sprintf("%d %d", num_taxa, num_chars), phylip_f)
  
  for(i in 1:num_taxa) {
    line <- paste(rownames(score_mat)[i], paste(score_mat[i, ], collapse = ""), sep = " ")
    writeLines(line, phylip_f)
  }
  
  close(phylip_f)
}

#' Write a score matrix as a character-state FASTA
#'
#' One record per row: `>` plus the row name, then the row's states with no
#' separator. A state is `as.integer()` of the score (truncated toward zero),
#' or `?` where the score is `NA`, which is how an integration a cell never
#' recovered is represented.
#'
#' @param scoremat Matrix (dense or sparse) whose rows are cells named by
#'   `rownames` and whose columns are unique mutations. Intended for 0/1 or
#'   `NA` scores.
#' @param output_fasta_path Path of the FASTA to write.
#' @return `NULL`; called for its side effect.
#' @section Side effects: Creates or overwrites `output_fasta_path`.
#' @note Because states go through `as.integer()`, an un-binarized allelic
#'   fraction in (0, 1) is written as `0` and a raw count of 10 or more takes
#'   more than one character, which shifts the rest of that record. Binarized
#'   or otherwise 0/1 matrices are the intended input.
new_scoremat_to_fasta <- function(scoremat, output_fasta_path){
  seq_vec <- apply(scoremat, 1, function(x){
    
    # write NAs as ? (for missing integrations, e.g.)
    x_filt <- ifelse(is.na(x), '?', as.character(as.integer(x)))
    return(paste0(x_filt, collapse = ''))
  })
  names(seq_vec) <- rownames(scoremat)
  cat(paste0('>', names(seq_vec), '\n', seq_vec),
      file = output_fasta_path,
      sep = '\n')
}

#' Build the score matrix for one set of cell profiles
#'
#' Gathers every non-zero entry of every cell's mutation matrix into a
#' cell-by-mutation sparse matrix. One column is one unique mutation event: for
#' `'bc'` the integration, position and mutation value together; for `'mt'` the
#' position and mutation value only, so the same variant on several genomes of
#' a cell is one character scored by its allelic fraction. Reference bases
#' (`0`) never get a column, so a cell that lacks a mutation simply scores 0 in
#' that column.
#'
#' @details Column names are `<mt|bc>_<mutation name>`, that is
#'   `bc_<integration>_<position>_<mutation value>` or
#'   `mt_<position>_<mutation value>`, where the mutation value is the profile
#'   value itself: `-1` for a deletion, `1`-`4` for a substitution to A/G/C/T,
#'   and the decimal insertion code for an insertion, so two different
#'   insertions at one position are two different characters. A condensed
#'   deletion run replaces the single position with `<first>_<last>` and keeps
#'   the value `-1`. Scores are the number of recovered integrations carrying
#'   the mutation for `'bc'`, and that count divided by the number of genomes
#'   recovered for the cell (0 when none were) for `'mt'`.
#'
#' @param profiles Named list of mutation matrices, one per cell, named by
#'   lineage string. Rows are the recovered integrations or mitochondrial
#'   genomes, in the same order as `recovered_ints[[cell]]`; columns are
#'   sequence positions.
#' @param condense Logical. `TRUE` collapses deletions at consecutive positions
#'   on one integration into a single column; `FALSE` keeps one column per
#'   deleted position.
#' @param urid Character run id naming the `output/score_mats/<urid>` tree.
#' @param savename_prefix Character. Base file name of the written artifacts.
#' @param mt_or_bc Character, `'mt'` or `'bc'`. Decides whether the integration
#'   is part of a mutation's identity and whether scores are allelic fractions
#'   or raw counts.
#' @param recovered_ints Named list, one integer vector per cell, of the
#'   integration or genome indices backing the rows of `profiles[[cell]]`.
#'   Their lengths are the denominators of the `'mt'` allelic fractions.
#' @param binarize_score Logical vector. `'mt'` only: one artifact per element,
#'   with every non-zero score set to 1 when the element is `TRUE`.
#' @param allelic_fraction_thresh Numeric vector in [0, 1]. `'mt'` only: one
#'   artifact per element, with scores below that threshold set to 0.
#' @param return_af_fracs Logical. `TRUE` returns the per-cell scores instead
#'   of writing any file.
#' @return In count mode (`return_af_fracs = TRUE`) a list ordered like
#'   `profiles`: a cell carrying mutations maps to a list of named scores
#'   covering only the mutations it carries, while a cell carrying none maps to
#'   a named numeric vector of zeros over every mutation name. In output mode,
#'   `NULL` after the artifacts are written. In either mode, if no cell carries
#'   a single mutation the return is instead a `length(profiles)` x 1 all-zero
#'   sparse matrix whose one column is named `control` and whose rownames are
#'   the lineage strings, so a caller in count mode has to check that it got a
#'   list back.
#' @section Side effects: Creates `output/score_mats/<urid>/matrices/af` and
#'   `output/score_mats/<urid>/phylips/af`, then writes
#'   `matrices/<savename_prefix><suffix>.rds` (the sparse matrix with row and
#'   column names attached) plus the matching
#'   `phylips/<savename_prefix><suffix>.fasta`. The suffix is empty for `'bc'`
#'   and `_AF_<threshold>_B_<T|F>` for each threshold/binarize pair for `'mt'`;
#'   both are written directly under `matrices/` and `phylips/`, leaving the
#'   `af/` subdirectories empty.
#' @note With `condense = TRUE` and at least one deletion present, the deletion
#'   branch calls `releid()`, which is defined neither here nor in data.table,
#'   so that combination errors.
#' @note Barcode scores are not masked: a mutation on an integration the cell
#'   never recovered scores 0 rather than `NA`, unlike `create_one_score_mat()`.
new_create_one_score_mat <- function(profiles,
                                     condense,
                                     urid,
                                     savename_prefix,
                                     mt_or_bc,
                                     recovered_ints = NULL,
                                     binarize_score = FALSE,
                                     allelic_fraction_thresh = 0,
                                     return_af_fracs = FALSE){
  
  if(!dir.exists(file.path('output', 'score_mats', urid))){
    dir.create(file.path('output', 'score_mats', urid, 'matrices', 'af'), recursive = TRUE)
    dir.create(file.path('output', 'score_mats', urid, 'phylips', 'af'), recursive = TRUE)
  }
  
  # print('prior to creating mut_combos_mat')
  mut_combos_mat <- rbindlist(lapply(seq_along(profiles), function(cell_num){
    
    linstring <- names(profiles)[cell_num]
    
    mut_mat <- profiles[[linstring]]
    
    these_recovered_ints <- recovered_ints[[linstring]]
    
    mut_coords <- which(mut_mat != 0, arr.ind = TRUE) # new 3/11
   
    if(nrow(mut_coords) > 0){
      dt <- data.table(
        linstring = linstring,
        ints_mutated = these_recovered_ints[mut_coords[, 1]],
        positions_mutated = mut_coords[, 2],
        mut_vals = mut_mat[mut_coords]
        
      )
      return(dt)
    }
    
  }), use.names = TRUE, fill = TRUE)
  
  # if mut_combos_mat is NULL, it means no mutations happened and we can/should exit early
  if(is.null(nrow(mut_combos_mat)) || nrow(mut_combos_mat) == 0){
    
    # return an empty sparse matrix (one zero value hard-coded in)
    af_mat <- sparseMatrix(i = 1,
                           j = 1,
                           x = 0,
                           dims = c(length(profiles), 1))
    
    colnames(af_mat) <- 'control'
    rownames(af_mat) <- names(profiles) 
    return(af_mat) 
    
  }
  
  # if we want to encode deletions as a single mutation rather than a different deletion event at each deleted position:
  if(condense){
    
    setkey(mut_combos_mat, linstring, ints_mutated, positions_mutated)
    dels <- mut_combos_mat[mut_vals == -1]
    nondels <- mut_combos_mat[mut_vals != -1]

    if(nrow(dels) > 0){
      dels[, run_id := releid(positions_mutated),
           by = .(linstring, ints_mutated)]
      grouped_dels <- dels[, .(ints_mutated = ints_mutated[1],
                               positions_mutated = paste0(min(positions_mutated), '_', max(positions_mutated)),
                               mut_vals = -1), by = .(linstring, run_id)]
      mut_combos_mat <- rbindlist(list(grouped_dels, nondels), use.names = TRUE)
    }
  }
    
  
  if(mt_or_bc == 'bc'){
    unique_pos_muts <- unique(mut_combos_mat[, .(ints_mutated, positions_mutated, mut_vals)])
    unique_pos_muts[, mut_idx := .I]
    # name mutations right away:
    unique_pos_muts[, mut_name := paste(ints_mutated, positions_mutated, mut_vals, sep = '_')]
    setkey(mut_combos_mat, ints_mutated, positions_mutated, mut_vals)
    join_dt <- mut_combos_mat[unique_pos_muts, .(linstring, mut_idx, mut_name), on = .(ints_mutated,
                                                                    positions_mutated,
                                                                    mut_vals)]
  } else if(mt_or_bc == 'mt'){
    unique_pos_muts <- unique(mut_combos_mat[, .(positions_mutated, mut_vals)])
    unique_pos_muts[, mut_idx := .I]
    # name mutations right away:
    unique_pos_muts[, mut_name := paste(positions_mutated, mut_vals, sep = '_')]
    setkey(mut_combos_mat, positions_mutated, mut_vals)
    
    join_dt <- mut_combos_mat[unique_pos_muts, .(linstring, mut_idx, mut_name), on = .(positions_mutated,
                                                                             mut_vals)]
  }
  


  counts_dt <- join_dt[, .N, by = .(linstring, mut_idx, mut_name)]
  setnames(counts_dt, c('linstring', 'mut_idx', 'mut_name', 'raw_count'))
  lin_to_profnum <- setNames(seq_along(names(profiles)), names(profiles))
  counts_dt[, profnum := lin_to_profnum[linstring]]
  
  if(mt_or_bc == 'mt'){
    
    rec_ints_counts <- vapply(recovered_ints, function(x){
      if(is.null(x) || length(x) == 0){
        return(0)
      } else{
        return(length(x))
      }},
      FUN.VALUE = integer(1))
    
    counts_dt[,  tot_rec := rec_ints_counts[profnum]]
    counts_dt[, af := ifelse(tot_rec > 0, raw_count/tot_rec, 0)]
    counts_dt[, x := af]
    
    
  } else{
    counts_dt[, x := raw_count]
    
  }
  
  cellmut_mat <- sparseMatrix(i = counts_dt$profnum,
                              j = counts_dt$mut_idx,
                              x = counts_dt$x,
                              dims = c(length(profiles), nrow(unique_pos_muts)))
  
  if(return_af_fracs){
    
    
    res_list <- lapply(split(counts_dt, by = 'linstring'), function(linstring_data){
      return(as.list(setNames(linstring_data$x, linstring_data$mut_name)))
    })
    
    zero_mut_cells <- setdiff(names(profiles), unique(mut_combos_mat$linstring))
    
    if(length(zero_mut_cells) > 0){
      
      null_vec <- setNames(rep(0, length(unique_pos_muts$mut_name)), unique_pos_muts$mut_name)
      zero_mut_list <- setNames(rep(list(null_vec), length(zero_mut_cells)), zero_mut_cells)
      res_list <- c(res_list, zero_mut_list)[names(profiles)]
    } else{
      res_list <- res_list[names(profiles)]
    }
    return(res_list)
    
  }
  
  #' Internal: label one score matrix and write its RDS and FASTA
  #'
  #' @param mat Sparse score matrix with rows in `names(profiles)` order and
  #'   columns in `unique_pos_muts` order.
  #' @param suffix Character appended to `savename_prefix` in both file names.
  #' @return `NULL`. Attaches `<mt_or_bc>_<mutation name>` colnames and
  #'   lineage-string rownames to `mat`, saves it to
  #'   `output/score_mats/<urid>/matrices/<savename_prefix><suffix>.rds` and
  #'   writes the character-state FASTA alongside it in `phylips/`.
  write_mat <- function(mat, suffix = ''){

    colnames(mat) <- paste(mt_or_bc, unique_pos_muts$mut_name, sep = '_')
    rownames(mat) <- names(profiles)
    
    rds_path <- file.path('output', 'score_mats', urid, 'matrices', paste0(savename_prefix, suffix, '.rds'))
    saveRDS(mat, rds_path)
    
    fasta_path <- file.path('output', 'score_mats', urid, 'phylips', paste0(savename_prefix, suffix, '.fasta'))
    new_scoremat_to_fasta(scoremat = mat, output_fasta_path = fasta_path)
  }
  
  if(mt_or_bc == 'bc'){
    write_mat(mat = cellmut_mat, suffix = '')
    
  } else{
    for(af_thresh in allelic_fraction_thresh){
      mat_af <- cellmut_mat
      if(af_thresh > 0){
        mat_af@x[mat_af@x < af_thresh] <- 0
      }
      for(binarize in binarize_score){
        temp_mat_af <- mat_af
        if(binarize){
          temp_mat_af@x[temp_mat_af@x > 0] <- 1
        }
        suffix <- paste0('_AF_', af_thresh, '_B_', substr(binarize, 1, 1))
        write_mat(temp_mat_af, suffix)
      }
    }
  }

}

# additional functionality to facilitate heteroplasmy calculations:
#' Build the score matrix for one set of cell profiles (legacy version)
#'
#' Takes the same arguments as `new_create_one_score_mat()` and produces the
#' same mutation naming and the same `output/score_mats/<urid>` layout, but is
#' built on base matrices and dplyr instead of data.table, and differs in what
#' it writes for barcodes and in where un-binarized mitochondrial artifacts
#' land.
#'
#' @details Barcode output is densified and binarized -- every non-zero score
#'   becomes 1 -- and then masked: for each cell, every column whose
#'   integration is absent from `recovered_ints[[cell]]` is set to `NA` so that
#'   an unrecovered integration is written as `?` rather than scored as
#'   reference. It is saved as `matrices/<savename_prefix>.rds` with its FASTA
#'   in `phylips/`. Mitochondrial output loops over thresholds and binarize
#'   flags as in the new version, but only the binarized matrices go to
#'   `matrices/` and `phylips/`; un-binarized ones go to the `matrices/af/` and
#'   `phylips/af/` subdirectories. Condensed deletion runs are named by
#'   `group_deletions()`, so a run appears at its last position with the
#'   mutation value `d<run length>` rather than `-1`.
#'
#' @param profiles Named list of per-cell mutation matrices, keyed by lineage
#'   string; rows are recovered integrations or genomes, columns positions.
#' @param condense Logical. `TRUE` collapses contiguous deletion runs.
#' @param urid Character run id naming the `output/score_mats/<urid>` tree.
#' @param savename_prefix Character. Base file name of the written artifacts.
#' @param mt_or_bc Character, `'mt'` or `'bc'`.
#' @param recovered_ints Named list of the integration or genome indices
#'   backing each cell's profile rows; also the barcode mask and the `'mt'`
#'   allelic-fraction denominators.
#' @param binarize_score Logical vector; `'mt'` only, one artifact per element.
#' @param allelic_fraction_thresh Numeric vector; `'mt'` only, one artifact per
#'   element, scores below the threshold set to 0.
#' @param return_af_fracs Logical. `TRUE` returns per-cell raw counts instead
#'   of writing files.
#' @return In count mode, a list ordered like `profiles` whose entries are
#'   named numeric vectors of raw counts across every mutation name, with
#'   `NULL` for cells carrying no mutation. If no cell carries a mutation, the
#'   `length(profiles)` x 1 all-zero sparse matrix with the single column
#'   `control`. Otherwise `NULL`, after the artifacts are written.
#' @section Side effects: Creates `output/score_mats/<urid>/matrices/af` and
#'   `output/score_mats/<urid>/phylips/af` and writes the RDS and FASTA
#'   artifacts described in the details.
#' @note Retained for compatibility; sim5_code.R calls
#'   `new_create_one_score_mat()`. The mutation table is built by `cbind()`
#'   against the lineage string, so its columns are character; with
#'   `condense = TRUE` and two or more deletions those character positions
#'   reach `group_deletions()`, whose run arithmetic then errors.
create_one_score_mat <- function(profiles,condense, urid, savename_prefix, mt_or_bc, 
                                 recovered_ints = NULL, binarize_score = FALSE, allelic_fraction_thresh = 0,
                                 return_af_fracs = FALSE){
  
  # create directories that will store score matrices and phy files
  if(!dir.exists(file.path('output', 'score_mats', urid))){
    dir.create(file.path('output', 'score_mats', urid, 'matrices', 'af'), recursive = TRUE)
    dir.create(file.path('output', 'score_mats', urid, 'phylips', 'af'), recursive = TRUE)
  }
  
  # get all combinations of cell x int x position x mutation
  all_mut_combos <- lapply(seq_along(profiles), function(cell_num){
    
    
    
    
    
    linstring <- names(profiles)[cell_num]
    
    mut_mat <- profiles[[linstring]]
    these_recovered_ints <- recovered_ints[[linstring]]

    mut_coords <- which(mut_mat != 0, arr.ind = TRUE) # new 3/11
    
   
    # map the respective row num to the int that was captured
    # works because which_ints_recovered has already been sorted
    ints_mutated <- sapply(mut_coords[,1], function(resp_int){
      these_recovered_ints[resp_int]
    })
    
    positions_mutated <- mut_coords[, 2]
    
    if(nrow(mut_coords) > 0){ 
      
      # iterate through the mut_coords and get the associated mutation values
      mut_vals <- apply(mut_coords, MARGIN = 1, FUN = function(row){      
        return(mut_mat[row[1], row[2]])
      })
      
      
      # final_mat will store the cell number, mutated integration and corresponding mutation positions, and the respective mutations themselves
      # final_mat <- cbind(cell_num, ints_mutated, positions_mutated, mut_vals)
      final_mat <- cbind(linstring, ints_mutated, positions_mutated, mut_vals)
      return(final_mat)
    }
  })
  
  # stack all list entries on top of one another to create matrix with same info
  mut_combos_mat <- do.call(rbind, all_mut_combos)
  

  # if mut_combos_mat is NULL, it means no mutations happened and we can/should exit early
  if(is.null(nrow(mut_combos_mat))){
    
    # return an empty sparse matrix (one zero value hard-coded in)
    af_mat <- sparseMatrix(i = 1,
                           j = 1,
                           x = 0,
                           dims = c(length(profiles), 1))
    
    colnames(af_mat) <- 'control'
    rownames(af_mat) <- names(profiles)
    return(af_mat)
    
  }
  
  # if we want to encode deletions as a single mutation rather than a different deletion event at each deleted position:
  if(condense){
    mut_combos_df <- as.data.frame(mut_combos_mat)
    old_colnames <- colnames(mut_combos_df)
    
    # deletions-only matrix that will be used for changing how deletions are labeled
    dels <- mut_combos_df %>%
      # arrange(cell_num, ints_mutated, positions_mutated) %>%
      arrange(linstring, ints_mutated, positions_mutated) %>%
      filter(mut_vals == -1)
    
    # non-deletions-only matrix that will be concatted to the newly-formatted del matrix
    nondels <- mut_combos_df %>%
      filter(mut_vals != -1)
    
    # null rownames
    
    # if no deletions
    if(nrow(dels) == 0){
      mut_combos_mat <- matrix(sapply(nondels, as.character), ncol = length(old_colnames),
                               nrow = nrow(nondels), dimnames = list(NULL, old_colnames))
    } else if(nrow(dels) == 1){ 
      # if there's only one deletion, looking for runs as in group_deletions() won't work
      # don't have to group the deletion at all 
      mut_combos_mat <- matrix(rbind(sapply(dels, as.character), sapply(nondels, as.character)), ncol = length(old_colnames),
                               nrow = nrow(dels) + nrow(nondels), dimnames = list(NULL, old_colnames))
    }
    else{
      grouped_deletion_mat <- matrix(group_deletions(dels), ncol = length(old_colnames), 
                                     dimnames = list(NULL, old_colnames))
      
      # rewrite mut_combos_mat (convert nondels to char type to enable stacking)
      mut_combos_mat <- matrix(rbind(grouped_deletion_mat, sapply(nondels, as.character)), ncol = length(old_colnames),
                               nrow = nrow(grouped_deletion_mat) + nrow(nondels), dimnames = list(NULL, old_colnames))  
    }
  }
  
  if(mt_or_bc == 'bc'){
    
    # get unique combinations of integrations x positions x mutations
    unique_pos_muts <- unique(mut_combos_mat[, c('ints_mutated', 'positions_mutated', 'mut_vals')])
    
    if(!is.matrix(unique_pos_muts)){
      unique_pos_muts <- matrix(unique_pos_muts, nrow = 1)
      colnames(unique_pos_muts) <- c('ints_mutated', 'positions_mutated', 'mut_vals')
    }
    
    
    
    mut_combos_mat <- data.table(mut_combos_mat)
    
    setkey(mut_combos_mat, ints_mutated, positions_mutated, mut_vals)
    
  } else if(mt_or_bc == 'mt'){
    # if we want to filter by allelic fraction > 0, don't include integration number in determining unique combos 
    unique_pos_muts <- unique(mut_combos_mat[, c('positions_mutated', 'mut_vals')])
    
    # 1-row matrix gets coerced to a vector
    if(!is.matrix(unique_pos_muts)){
      unique_pos_muts <- matrix(unique_pos_muts, nrow = 1)
      colnames(unique_pos_muts) <- c('positions_mutated', 'mut_vals')
    }
    
    mut_combos_mat <- data.table(mut_combos_mat)
    
    
    setkey(mut_combos_mat, positions_mutated, mut_vals)
  }
  
  cell_nums_with_mut <- lapply(seq(1, nrow(unique_pos_muts)), function(rowvals_ind){ 
    
    # rowvals will have c(ints_mutated	positions_mutated	mut_vals)
    rowvals <- unique_pos_muts[rowvals_ind, ]
    
    if(mt_or_bc == 'bc'){
      # match the entire cell number x position x mutation matrix to this specific unique int x position x mutation, keep cell number
      cell_nums_with_mut <- mut_combos_mat[.(rowvals[1], rowvals[2], rowvals[3])]$linstring
    } else if(mt_or_bc == 'mt'){ # don't match on integration
      cell_nums_with_mut <- mut_combos_mat[.(rowvals[1], rowvals[2])]$linstring
    }
    
    return(cell_nums_with_mut)
    
  })
  
  # rows are cells
  # columns are unique mutations 
  # values are 1 if cell has that mutation
  mutnames <- apply(unique_pos_muts, MARGIN = 1, function(rowvals){return(paste(mt_or_bc, paste(rowvals, collapse = '_'), sep = '_'))})
  
  if(return_af_fracs){
    ##################### return a list that maps each cell to the RAW COUNTS of each heteroplasmy mutation in that cell across all genomes

    names(cell_nums_with_mut) <- mutnames
    
    counts_df <- data.frame(
      mut_id = rep(names(cell_nums_with_mut), lengths(cell_nums_with_mut)),
      counts = unlist(cell_nums_with_mut)
    )
    
    # Create a contingency table
    counts_table <- table(counts_df$counts, counts_df$mut_id)
    
    # Convert back to list
    result_list <- lapply(
      seq_len(nrow(counts_table)),
      function(i) {
        as.numeric(counts_table[i, ])
      }
    )
    
    names(result_list) <- rownames(counts_table)
    
    # Add mut_id names to each element
    result_list <- lapply(result_list, function(x) {
      names(x) <- colnames(counts_table)
      x
    })
    
    # get linstrings that have no mutations
    linstrings_with_no_muts <- setdiff(names(profiles), names(result_list))
    
   
    # create empty list of nulls with length == length(# missing linstrings)
    null_list <- vector('list', length(linstrings_with_no_muts))
    names(null_list) <- linstrings_with_no_muts
    # combine nulls with mutation counts and preserve order
    result_list <- c(result_list, null_list)[names(profiles)]
    
    return(result_list)
  }
  
  linstring_to_rownum <- setNames(seq_along(names(profiles)), names(profiles))
  ijx_triples <- lapply(seq_along(cell_nums_with_mut), function(col_idx){
    cells_this_mut <- cell_nums_with_mut[[col_idx]]
    cell_counts <- table(cells_this_mut)
    
    
    ivals <- linstring_to_rownum[names(cell_counts)]
    jvals <- rep(col_idx, length(cell_counts))
    
    raw_counts <- as.integer(cell_counts)
    if(mt_or_bc == 'mt'){
      
      # if there are no ints associated with a given cell, return a zero
      num_recovered_ints <- vapply(recovered_ints[names(cell_counts)], FUN = function(ints){
        if(is.null(ints) || length(ints) == 0L){
          return(0L)
        }
        return(length(ints))
      },
      FUN.VALUE = integer(1))
      afs_raw <- raw_counts / num_recovered_ints
      afs <- ifelse(is.finite(afs_raw), afs_raw, 0)
      
    } else{
      afs <- raw_counts
    }
    
    
    
    return(list('ivals' = ivals,
                'jvals' = jvals,
                'xvals' = afs))
  })
  
  i_vals <- unlist(lapply(ijx_triples, `[[`, 'ivals'))
  j_vals <- unlist(lapply(ijx_triples, `[[`, 'jvals'))
  af_vals <- unlist(lapply(ijx_triples, `[[`, 'xvals'))
  cellmut_mat <- sparseMatrix(i = i_vals,
                              j = j_vals,
                              x = af_vals,
                              dims = c(length(profiles), length(mutnames)))
  
  
  if(mt_or_bc == 'bc'){
    
    cellmut_mat_cp1 <- as(cellmut_mat, class(cellmut_mat)[1]) # introducing in case further processing steps are introduced
    
    # for barcodes, we want to binarize (don't care about allelic fractions at positions)
    cellmut_mat_cp1@x[cellmut_mat_cp1@x > 0] <- 1
    
    cellmut_mat_cp1 <- matrix(cellmut_mat_cp1, nrow = dim(cellmut_mat_cp1)[1],
                              ncol = dim(cellmut_mat_cp1)[2])
    colnames(cellmut_mat_cp1) <- mutnames
    rownames(cellmut_mat_cp1) <- names(profiles)
    
    # CREATE A NA MASK FOR MUTATIONS OF CELLS LACKING RESP RECOVERED INTEGRATIONS IN BARCODE STEPS
    # missing integrations should not be informative (which would be the case if ints encoded as missing were assigned 0s)
    int_nums_in_muts <- as.integer(sub('^bc_([0-9]+)_.*$', '\\1', colnames(cellmut_mat_cp1)))
    
    for(cell in rownames(cellmut_mat_cp1)){
      rec <- recovered_ints[[cell]]
      if(is.null(rec)){
        rec <- integer(0)
      }
      
      missing_features <- which(!(int_nums_in_muts %in% rec))
      if(length(missing_features) > 0){
        cellmut_mat_cp1[cell, missing_features] <- NA
      }
    }
    
    
    saveRDS(cellmut_mat_cp1, file.path('output', 'score_mats', urid, 'matrices', paste0(savename_prefix, '.rds')))

    new_scoremat_to_fasta(scoremat = cellmut_mat_cp1, 
                          output_fasta_path = file.path('output', 'score_mats', urid, 'phylips', paste0(savename_prefix, '.fasta')))
    
  } else if(mt_or_bc == 'mt'){
    for(af_thresh in allelic_fraction_thresh){
      
      # copies so that filtering steps operate on original, not previously-filtered mats
      cellmut_mat_cp1 <- as(cellmut_mat, class(cellmut_mat)[1])
      
      # if there's an allelic fraction threshold, any fracs below the threshold are set to 0
      cellmut_mat_cp1@x[cellmut_mat_cp1@x < af_thresh] <- 0
      
      for(binarize in binarize_score){
        
        cellmut_mat_cp2 <- as(cellmut_mat_cp1, class(cellmut_mat_cp1)[1])
        
        if(binarize){
          # if binarizing scores, set any non-zero element to 1
          cellmut_mat_cp2@x[cellmut_mat_cp2@x > 0] <- 1
          colnames(cellmut_mat_cp2) <- mutnames
          rownames(cellmut_mat_cp2) <- names(profiles)
          
          saveRDS(cellmut_mat_cp2, file.path('output', 'score_mats', urid, 'matrices', paste0(savename_prefix, 
                                                                                              '_AF_', af_thresh,
                                                                                              '_B_', substr(binarize, 1, 1),
                                                                                              '.rds')))
  
          new_scoremat_to_fasta(scoremat = cellmut_mat_cp2, 
                                output_fasta_path = file.path('output', 'score_mats', urid, 'phylips', paste0(savename_prefix, '_AF_', af_thresh,
                                                                                                              '_B_', substr(binarize, 1, 1), '.fasta')))
          
        } else{
          colnames(cellmut_mat_cp2) <- mutnames
          rownames(cellmut_mat_cp2) <- names(profiles)
          
          saveRDS(cellmut_mat_cp2, file.path('output', 'score_mats', urid, 'matrices', 'af', paste0(savename_prefix, 
                                                                                                    '_AF_', af_thresh,
                                                                                                    '_B_', substr(binarize, 1, 1),
                                                                                                    '.rds')))
         
          new_scoremat_to_fasta(scoremat = cellmut_mat_cp2, 
                                output_fasta_path = file.path('output', 'score_mats', urid, 'phylips', 'af', paste0(savename_prefix, '_AF_', af_thresh,
                                                                                                                    '_B_', substr(binarize, 1, 1), '.fasta')))
        }
        
        
        
      }
    }
  }
}




#' Score each cell's heteroplasmy burden
#'
#' Draws a severity score for every variant that is not yet in the severity
#' map, then returns, per cell, the dot product of that cell's per-mutation
#' values with the severity scores of the mutations the two share. Scores can
#' be negative or positive, since severities are drawn from a two-component
#' mixture centred at -1 and +1.
#'
#' @param cell_mut_counts Named list, one entry per cell, of named per-mutation
#'   values -- what `new_create_one_score_mat(return_af_fracs = TRUE)` returns,
#'   allelic fractions in the `'mt'` case.
#' @param heteroplasmy_severity_score_list Named list mapping mutation name to
#'   a numeric severity score. Variants missing from it are scored by
#'   `draw_severity_scores()`, which sim5_code.R defines.
#' @param positive_score_weight Numeric. Passed as `mean2_weight`: the weight
#'   of the +1 component of the severity mixture.
#' @param hetero_sd Numeric. Passed as `sigma`, the spread of both components.
#' @param cell_population List of cell objects; part of the signature but not
#'   used by the body.
#' @param cell_to_num_mito_genomes_list Named list of per-cell genome counts,
#'   read only when `normalize_cell_mut_counts` is `TRUE`.
#' @param normalize_cell_mut_counts Logical. `TRUE` divides each cell's values
#'   element-wise by its genome count into a local `norm_af`, which the scoring
#'   step below does not read, so the returned scores are unaffected either
#'   way.
#' @return A named list with one numeric score per entry of `cell_mut_counts`,
#'   in the same order.
#' @section Side effects: Consumes draws from the R random number stream, one
#'   per newly seen variant.
#' @note The variants eligible for a new severity score are read from the first
#'   entry of `cell_mut_counts` alone. A variant carried only by later cells
#'   therefore never enters the severity map, and the name intersection used
#'   for the dot product then drops it from those cells' scores.
#' @note The `<<-` update binds to this function's own
#'   `heteroplasmy_severity_score_list` argument, which shadows any outer
#'   binding, so newly drawn severities live only for the duration of the call
#'   and the caller's map comes back unchanged.
get_norm_cell_heteroplasmy_scores <- function(cell_mut_counts,
                                          heteroplasmy_severity_score_list,
                                          positive_score_weight,
                                          hetero_sd,
                                          cell_population,
                                          cell_to_num_mito_genomes_list = NULL,
                                          normalize_cell_mut_counts = FALSE){
  
  # get all unique mutation names by accessing the first cell's data
  # since counts for all unique mutations are calculated for each cell
  mutnames <- names(cell_mut_counts[[names(cell_mut_counts[1])]])
  
  
  # if cell_mut_counts 
  if(normalize_cell_mut_counts){
    
    common_names <- intersect(names(cell_mut_counts), names(cell_to_num_mito_genomes_list))
    
    
    norm_af <- setNames(Map('/', cell_mut_counts[common_names],
                            cell_to_num_mito_genomes_list[common_names]),
                        common_names)
  }
    
  # generate heteroplasmy severity score for each newly observed variant, and globally update the heteroplasmy serverity score list
  new_muts <- setdiff(mutnames, names(heteroplasmy_severity_score_list))
  invisible(
    sapply(new_muts, function(mut){
      heteroplasmy_severity_score_list[[mut]] <<- draw_severity_scores(num_draws = 1, 
                                                                       mean2_weight = positive_score_weight, 
                                                                       mean1 = -1, 
                                                                       mean2 = 1, 
                                                                       sigma = hetero_sd)
    })
  )
  
  # for each cell, find dot product between that cell's AFs and the AF severity scores
  # ensure compatible order between normalized allelic fractions and severity scores before computing dot product
  cell_hetero_scores <- lapply(cell_mut_counts, function(cell){
    
    # ensure consistent order between normalized allelic fractions and severity scores
    ordered_muts <- intersect(names(cell), names(heteroplasmy_severity_score_list))
    af_scores <- as.numeric(cell[ordered_muts])
    severity_scores <- as.numeric(heteroplasmy_severity_score_list[ordered_muts])
    weighted_score <- as.numeric(af_scores %*% severity_scores)
    return(weighted_score)
    
  })
  
  return(cell_hetero_scores)
  
  
}
