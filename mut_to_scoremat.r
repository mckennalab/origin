group_deletions <- function(deletion_df){
  
  streak <- 1
  
  new_del_mat_list <- lapply(seq(2, nrow(deletion_df)), function(deletion_num){
    
    # if part of the same cell and integration AND
    # if the current mutation position is one greater than the previous
    # BUG: misplaced parenthesis. The intended check is
    #   positions_mutated[deletion_num] == positions_mutated[deletion_num - 1] + 1
    # but the parens make it (positions == positions_prev) + 1, i.e. logical->int yielding
    # 1 (when equal) or 2 (when not — TRUE? no, FALSE+1=1, TRUE+1=2). Either way the value
    # is non-zero and `&` treats it as TRUE, so the adjacency check is dead. As written,
    # any pair of deletions sharing linstring + integration is treated as a contiguous
    # streak regardless of whether the positions actually neighbour each other.
    if((deletion_df$linstring[deletion_num] == deletion_df$linstring[deletion_num - 1]) &
       (deletion_df$ints_mutated[deletion_num] == deletion_df$ints_mutated[deletion_num - 1]) &
       (deletion_df$positions_mutated[deletion_num] == deletion_df$positions_mutated[deletion_num - 1]) + 1){
      streak <<- streak + 1
      
      # check if this deletion is the last one in the dataset
      if(deletion_num == nrow(deletion_df)){
        
        return_vec <- c(deletion_df$linstring[deletion_num],
                        deletion_df$ints_mutated[deletion_num],
                        deletion_df$positions_mutated[deletion_num], 
                        paste0('d', streak))
        
        return(return_vec)
      }
      
      return(NULL)
    } else{
      
      return_vec <- c(deletion_df$linstring[deletion_num],
                      deletion_df$ints_mutated[deletion_num],
                      deletion_df$positions_mutated[deletion_num], 
                      paste0('d', streak))
      streak <<- 1
      return(return_vec)
    }
  }) 
  
  new_del_mat <- do.call(rbind, new_del_mat_list)                                      
  
  return(new_del_mat)
}

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
create_one_score_mat <- function(profiles,condense, urid, savename_prefix, mt_or_bc, 
                                 recovered_ints = NULL, binarize_score = FALSE, allelic_fraction_thresh = 0,
                                 return_af_fracs = FALSE){
  
  # create directories that will store score matrices and phy files
  if(!dir.exists(file.path('output', 'score_mats', urid))){
    dir.create(file.path('output', 'score_mats', urid, 'matrices', 'af'), recursive = TRUE)
    dir.create(file.path('output', 'score_mats', urid, 'phylips', 'af'), recursive = TRUE)
  }
  
  # get all combinations of cell x int x position x mutation
  all_mut_combos <- lapply(seq(1, length(profiles)), function(cell_num){
    
    
    
    
    
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
      # BUG: this checks the entire `recovered_ints` list (which is never NULL at this
      # point — it's the function argument we just indexed into) instead of the per-cell
      # lookup `rec`. For cells that aren't keys in `recovered_ints`, `rec` is NULL but
      # this guard never fires, and the `which(!(int_nums_in_muts %in% rec))` below ends
      # up flagging EVERY column as missing (since `x %in% NULL` is all FALSE), so those
      # rows get fully NA-masked rather than partially. Should be `if(is.null(rec))`.
      if(is.null(recovered_ints)){
        rec <- integer(0)
      }
      
      missing_features <- which(!(int_nums_in_muts %in% rec))
      if(length(missing_features) > 1){
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


