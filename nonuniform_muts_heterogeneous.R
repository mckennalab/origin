# allow for heterogeneous rates within HMLB classes:
suppressPackageStartupMessages({
  library(parallel)
  library(Matrix)
  library(zeallot)  
})

# helper function that takes in three lists as input:
# map from edit window to genomic positions
# map from genomic positions to edit window
# map from genomic position to elig ints that HAVE NOT YET BEEN FILTERED BY WINDOW IF NECESSARY
# and returns an updated position --> elig ints mapping that conforms to the restriction
# that if ANY position in an editing window has already been edited, no further 
# non-uniform mutations can occur at positions within this same editing window
filter_elig_ints_by_edit_window <- function(pos_to_window_inds_list,
                                            window_to_pos_inds_list,
                                            pos_to_unedited_int_list){
  filt_elig_ints <- lapply(window_to_pos_inds_list, function(all_pos_in_window){
    
    growing_elig_ints <- lapply(all_pos_in_window, function(pos){
      if(as.character(pos) %in% names(pos_to_unedited_int_list)){
        return(pos_to_unedited_int_list[[as.character(pos)]])
      } else{
        return(NULL)
      }
    })
    
    non_null <- growing_elig_ints[!vapply(growing_elig_ints, is.null, logical(1))]
    if(length(non_null) == 0){
      return(integer(0))
    }
    # Keep empty vectors in the intersection. An empty eligible set at any
    # edited position closes the entire window, as required by the contract.
    return(Reduce(intersect, non_null))
  })
  
  updated_elig_ints <- lapply(pos_to_window_inds_list, function(windows){
    new_poss_eligs <- sapply(windows, function(window){
      return(filt_elig_ints[[window]])  
    })
    return(as.integer(unique(unlist(new_poss_eligs))))
  })
  return(updated_elig_ints)
}
 
# returns c(i_coords, j_coords) which can be used to build new mutation matrix
# this function only operates on High/Medium/Low bases; the background is encoded with uniform!
# have to make sure we pass in the correct bg edit rate on non-uniform
non_uniform_editing <- function(pos_er_list, num_integrations, eligible_ints, timepoint_savename = '', length1_positions = NULL){
  
  

  target_positions <- names(pos_er_list)
  pos_int_list <- setNames(lapply(target_positions, function(x){ # iterate through target base positions
    char_x <- as.character(x)
    
   
    er <- pos_er_list[[char_x]]
    
    if(!is.null(eligible_ints)){ # if we pass in eligible integrations, max number that can be edited is num unedited
      
     
      if(length(eligible_ints[[char_x]]) == 0){ # if we are fully saturated, return no edits
        return(0)
      } 
      

     
      # if there are still eligible integrations that can be edited:
      num_ints_edited <- rbinom(n = 1, size = length(eligible_ints[[char_x]]), prob = min(c(er, 1)))
      
      
      # Draw a binomial count, then choose that many distinct eligible integrations.
      if(num_ints_edited > 0){
        
        # because R is very dumb and cannot sample() an element from a length 1 integer vector, need to take an extra step if there's only one eligible integration to edit
        if(length(eligible_ints[[char_x]]) == 1){
          
          if(num_ints_edited == 1){ # this should always fire if we've entered the first if
            which_ints_edited <- as.numeric(eligible_ints[[char_x]])
          } else{ # should never fire
            cat('\nERROR: MORE THAN 1 INT EDITED WHEN ONLY 1 IS ELIGIBLE\n', 
                file = file.path('output', 'run_logs', unique_run_id, paste0('runlog_', unique_run_id, '.txt')),
                append = TRUE)
            quit(save = 'no', status = 1)
          }
          
        } else{
          
          which_ints_edited <- sample(
            x = eligible_ints[[char_x]],
            size = num_ints_edited,
            replace = FALSE
          )
          
          return(which_ints_edited)
          
        }
        
        
        
      } else{ # if num_ints_edited == 0
        return(0) # no edits made
      } 
    }
    else{ # if no list of eligible ints was passed in 
      
      # should not be firing under new logic ... 
      cat('inelig. ints list passed in\n', 
          file = file.path('output', 'run_logs', unique_run_id, paste0('runlog_', unique_run_id, '.txt')), 
          append = TRUE)
      
    }
    
  }), target_positions)
   
  # which integrations were edited for respective base positions, correspond to row values in mutation matrix
  temp_i_coords <- unname(unlist(pos_int_list))

  temp_j_coords <- rep(
    as.integer(names(pos_int_list)),
    lengths(pos_int_list)
  )
  
  
  # remove duplicates
  coords <- unique(mapply(list, temp_i_coords, temp_j_coords, SIMPLIFY=F))
  i_coords <- sapply(coords, function(x){return(x[[1]])})
  j_coords <- sapply(coords, function(x){return(x[[2]])})
  
  # i_coords was set to 0 in pos_int_list as an indicator when no edits were made 
  zero_inds <- which(i_coords == 0)
  
  if(length(zero_inds) == length(i_coords)){ # if none of the positions was edited
    
    return_list <- list('i_coords' = FALSE, 'j_coords' = FALSE)
    
  } else{ # if some of the positions were edited
    if(length(zero_inds) > 0){ # at least one zero, but not all zeros, # remove elements that were set to 0 by default (no edits)
      i_coords <- i_coords[-zero_inds]
      j_coords <- j_coords[-zero_inds]
    }
    # no need for else, but it would be no zeros at all (all edited)
    
    return_list <- list('i_coords' = i_coords, 'j_coords' = j_coords)
  }
  
  return(return_list)
  
}




get_background_edit_inds <- function(num_rows, num_cols, bg_pos_er_list, mut_type,
                                     sample_transversion = FALSE, verbose = FALSE){
  # Accepts sparse matrix as input, and adds to it the transitions that occur
  
  
  # if sample_transversion, we need to select which mutation is occurring at each site each time
  # if we are working with transversions, we have two different substitutions that can occur at each position
  # each time we call this function, we'll generate a different combination of transversion rates across positions
  # as long as force_transversion == FALSE
  # we also recover the base to which the outgoing transversion base converts
  if(sample_transversion){
   
  
    chosen_rel_base <- as.integer(sapply(bg_pos_er_list, FUN = function(x){
      
      # if no transversions can occur, return arbitrary selection. mut won't occur downstream.
      if(all(x == 0)){
        return(1)
      }
      return(sample(c(1,2), size = 1, prob = as.numeric(x)))
    }))

    # extract the mutation probability at each position according to which of the 2 bases was picked as transversion
    selected_probs <- sapply(seq(1:length(chosen_rel_base)), function(x){
      return(bg_pos_er_list[[x]][chosen_rel_base[x]])
    })
    
    # extract the to-base identity
    selected_bases_to <- sapply(names(selected_probs), function(x){
      return(which(c('A', 'G', 'C', 'T') == x))
    })

    # make probs a numeric vector (ie remove base names)
    selected_probs <- as.numeric(selected_probs)
    
    # rewrite the pos:er list
    bg_pos_er_list <- as.list(selected_probs)
  }   
  
  position_edit_probs <- as.numeric(unlist(bg_pos_er_list, use.names = FALSE))
  if(length(position_edit_probs) != num_cols){
    stop('bg_pos_er_list must contain one edit probability per matrix column.')
  }
  if(any(!is.finite(position_edit_probs)) ||
     any(position_edit_probs < 0 | position_edit_probs > 1)){
    stop('Background edit probabilities must be finite values between zero and one.')
  }

  # Drawing a binomial count independently for each position is equivalent to
  # one Bernoulli draw per matrix coordinate. Sampling rows without replacement
  # then realizes that count without creating duplicate coordinates.
  edits_per_position <- rbinom(
    n = num_cols,
    size = num_rows,
    prob = position_edit_probs
  )
  edited_positions <- which(edits_per_position > 0)
  num_edits <- sum(edits_per_position)
  if(verbose){
    cat(paste0('\nnum_edits for ', mut_type, ' == ', num_edits, '\n'), 
        file = file.path('output', 'run_logs', unique_run_id, paste0('runlog_', unique_run_id, '.txt')), 
        append = TRUE)  
  }
  
  if(num_edits == 0){
    return_list <- list('num_edits' = 0, 'i_coords' = c(), 'j_coords' = c())
    return(return_list)
    
  }

  i_coords <- unlist(
    lapply(edits_per_position[edited_positions], function(position_count){
      sample.int(num_rows, size = position_count, replace = FALSE)
    }),
    use.names = FALSE
  )
  j_coords <- rep(edited_positions, edits_per_position[edited_positions])
  

  return_list <- list('num_edits' = num_edits, 'i_coords' = i_coords, 'j_coords' = j_coords)
  
  if(sample_transversion){
    
    return_list[['selected_bases_to']] <- selected_bases_to
  }
  
  return(return_list)
  
  
}



transition_func <- function(mut_mat, num_rows, num_cols, baseline_ints, 
                            bg_transition_pos_er_list,
                            target_transition_pos_er_list = NULL,
                            timepoint_filename = '',
                            verbose = FALSE,
                            target_to_window_ind_list = NULL, 
                            window_to_target_ind_list = NULL, 
                            close_window = FALSE){
  # we only accept the BE pos er list in transitions because there shouldn't be elevated rates of indels with BE
  
  transition_matches <- c(2,1,4,3)
  
  post_indices_transition_func <- function(i_coords, j_coords, incoming_mat, match_transition_bases){
    
    
    # will look something like: also have to create a vector of length i_coords == length j_coords 
    # with the existing incoming_mat values at those i, j pairs
    
    # get the current values at that position in the incoming_mat. these will influence mutation outcome
    existing_mat_vals <- sapply(seq(1, length(i_coords)), function(x){
      return(incoming_mat[i_coords[x], j_coords[x]])
    })
    
    x_vals <- sapply(seq(1, length(i_coords)), function(x){

      
      if(existing_mat_vals[x] == 0){ # if no mutation already exists at this position
        
        # can refer to the unedited baseline sequence for the base at this position
        # this used to be: 
        # return(bases[match(baseline_ints[i_coords[x]], match_transition_bases)])
        
        this_base <- match_transition_bases[as.integer(baseline_ints[j_coords[x]])]
        if(length(this_base) > 1){
          cat(paste0('\nTRANSITIONNEWBASE HAS LENGTH > 1 in unedited position == ', this_base), 
              file = file.path('output', 'run_logs', unique_run_id, paste0('runlog_', unique_run_id, '.txt')), 
              append = TRUE)
        }
        return(this_base)
      }
      else if(existing_mat_vals[x] %% 1 != 0){ # if an insertion exists at this position
        # induce a certain base in the insertion to mutate
        
        if(existing_mat_vals[x] > 0){ # this matters when determining how many characters to remove (e.g. here, 0.23123)
          num_editable_bases <- nchar(existing_mat_vals[x]) - 1 # subtract 1 for the .
        }
        else if(existing_mat_vals[x] < 0){ # e.g. -1.23123
          num_editable_bases <- nchar(existing_mat_vals[x]) - 2 # subtract 1 for the - sign, 1 for the 1, and 1 for the .
        }
        # convoluted way of sampling a digit from a float without converting it to a string:
        rand_exp <- sample(seq(1, num_editable_bases), size = 1)
        
        base_to_mutate <- abs(round(existing_mat_vals[x] * 10**(rand_exp-1))) %% 10
        
        # if we have chosen to mutate the 0 to the left of the decimal in an insertion
        if(base_to_mutate == 0){
          base_to_mutate <- as.integer(baseline_ints[as.integer(j_coords[x])]) # rewrite base to be edited from 0 to int representation of original sequence
        }
        
        newbase <- match_transition_bases[base_to_mutate]
        
        
        # we recover the base above and match it here
        # this difference will give the sum of what has to be added to go from old base to new base
        this_base <- 10**(-1*(rand_exp - 1))*(newbase - base_to_mutate)
        if(length(this_base) > 1){
          cat(paste0('\nTRANSITION NEWBASE HAS LENGTH > 1 in insertion mutation == ', this_base), 
              file = file.path('output', 'run_logs', unique_run_id, paste0('runlog_', unique_run_id, '.txt')), 
              append = TRUE)
        }
        return(this_base)
      }
      
      else if(existing_mat_vals[x] == -1){ # if a deletion has already occurred here
        return(0) # we won't add anything to the mutation matrix
      }
      else{ # if the base is a point mutation
      
        this_base <- match_transition_bases[existing_mat_vals[x]] - existing_mat_vals[x]

        return(this_base) 
      }
    })
    
    new_muts <- Matrix::sparseMatrix(i = i_coords, j = j_coords, 
                             x = x_vals, dims = c(num_rows, num_cols))
    
    incoming_mat <- incoming_mat + new_muts
    return(incoming_mat)
  }
  
  uniform_res <- get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
                                          bg_pos_er_list = bg_transition_pos_er_list,
                                          mut_type = 'transition',
                                          verbose = verbose)

  num_transitions <- uniform_res[['num_edits']]
  transition_i_coords <- uniform_res[['i_coords']]
  transition_j_coords <- uniform_res[['j_coords']]
  
  # get_background_edit_inds() will return FALSE if no edits have occurred
  if(num_transitions != 0){
    if(verbose){
      cat(paste0('\nin num_transitions\n'), file = 'no_strings.txt', append = TRUE)  
    }
    
    mut_mat <- post_indices_transition_func(i_coords = transition_i_coords, j_coords = transition_j_coords, 
                                            incoming_mat = mut_mat, match_transition_bases = transition_matches)
    if(verbose){
      cat(paste0('\nsum(mut_mat) == ', sum(mut_mat), '\n'), 
          file = file.path('output', 'run_logs', unique_run_id, paste0('runlog_', unique_run_id, '.txt')), 
          append = TRUE)
    }
  }
  

  
  # perform uniform edits

  # if the transition pos er list is null or has length zero, there is no non-uniform editing
  if(is.null(target_transition_pos_er_list)){ # if uniform, we are done after this one step
    if(verbose){
      cat(paste0('\nin is.null(target_transition_pos_er_list)\n'), 
          file = file.path('output', 'run_logs', unique_run_id, paste0('runlog_', unique_run_id, '.txt')), 
          append = TRUE)  
    }
    
    
    return(mut_mat)
  }
  if(length(target_transition_pos_er_list) == 0){
    if(verbose){
      cat(paste0('\nin length(target_transition_pos_er_list) == 0\n'), 
          file = file.path('output', 'run_logs', unique_run_id, paste0('runlog_', unique_run_id, '.txt')), 
          append = TRUE)  
    }
    
    return(mut_mat)
  }

  
  elig_ints_list <<- list()
  
  for(pos in names(target_transition_pos_er_list)){
    
    # iterate through each position that is High/Medium/Low, and see which integrations haven't been edited at that position yet
    unedited_rowvals <- which(mut_mat[, as.integer(pos)] == 0)
    elig_ints_list[[as.character(pos)]] <- unedited_rowvals
  }
  # for non-uniform editing, we can optionally perform further refinement to determine eligible ints
  # this is based on editing windows. an edit in any position in a window disqualifies future edits 
  # from occurring in bases belonging to this same window
  
  if(close_window){
    elig_ints_list <- filter_elig_ints_by_edit_window(pos_to_window_inds_list = target_to_window_ind_list,
                                                window_to_pos_inds_list = window_to_target_ind_list,
                                                pos_to_unedited_int_list = elig_ints_list)
  }
  
  nu_res <- non_uniform_editing(pos_er_list = target_transition_pos_er_list,
                                eligible_ints = elig_ints_list,
                                num_integrations = NULL,
                                timepoint_savename = timepoint_filename
  )

  nu_transition_i_coords <- nu_res[['i_coords']]
  nu_transition_j_coords <- nu_res[['j_coords']]
  
  if(nu_transition_i_coords[1] != FALSE){
    mut_mat <- post_indices_transition_func(i_coords = nu_transition_i_coords,
                                            j_coords = nu_transition_j_coords,
                                            incoming_mat = mut_mat,
                                            match_transition_bases = transition_matches)
  }

  return(mut_mat)
  

  
}


# the reason this is working is because the targets that are being generated are all the same base. 
# so force_transversions isn’t absolutely essential because the selection of which positions to edit in the target editing workflow does that already. 
# the only difference would be in an editing window context, where bases can differ. 
# i’d argue it’s actually better to NOT force transversions in these cases.

transversion_func <- function(mut_mat, num_rows, num_cols, bg_transversion_pos_er_list, baseline_ints, bg_sub_prob_mat,
                              target_transversion_pos_er_list = NULL, force_target_transversions = FALSE, verbose = FALSE,
                              target_transversion_to_base = NULL,
                              target_to_window_ind_list = NULL, 
                              window_to_target_ind_list = NULL, 
                              close_window = FALSE){
  # we only accept the BE pos er list in transversions because there shouldn't be elevated rates of indels with BE
  
  # Accepts sparse matrix as input, and adds to it the transversions that occur
  

  transversion_matches <- list(c(3,4), c(3,4), c(1,2), c(1,2))
  forced_transversion_matches <- list(3, 4, 1, 2)
  
  post_indices_transversion_func <- function(i_coords, j_coords, incoming_mat, bases_going_to, force_transversions){
    
    # the new version of this expression closely parallels the transition one
    # get the current values at that position in the incoming_mat. these will influence mutation outcome
    existing_mat_vals <- sapply(seq(1, length(i_coords)), function(x){
      return(incoming_mat[i_coords[x], j_coords[x]])
    })
    
    
    x_vals <- sapply(seq(1, length(i_coords)), function(x){
      
      if(existing_mat_vals[x] == 0){ # if no mutation already exists at this position, can mutate

        this_base <- bases_going_to[j_coords[x]]
        
        return(this_base)
       
      }
      else if(existing_mat_vals[x] %% 1 != 0){ # if an insertion exists at this position
        
      
        # find a certain base in the insertion to mutate
        if(existing_mat_vals[x] > 0){ # this matters when determining how many characters to remove (e.g. here, 0.23123)
          
          num_editable_bases <- nchar(existing_mat_vals[x]) - 1 # subtract 1 for the .
        }
        else if(existing_mat_vals[x] < 0){ # e.g. -1.23123
          # i think this should be -2
          num_editable_bases <- nchar(existing_mat_vals[x]) - 2 # subtract 1 for the - sign, 1 for the 1, and 1 for the .
        }
        # convoluted way of sampling a digit from a float without converting it to a string:
        rand_exp <- sample(seq(1, num_editable_bases), size = 1)
        

        base_to_mutate <- abs(round(existing_mat_vals[x] * 10**(rand_exp-1))) %% 10
        
        # if we have chosen to mutate the 0 to the left of the decimal in an insertion
        if(base_to_mutate == 0){
          base_to_mutate <- baseline_ints[j_coords[x]] # rewrite base to be edited from 0 to int representation of original sequence
          
          # will look something like existing_mat_vals[x] + new base
          # take difference between new base and original base
        }
        
        # get the two options that the base can undergo a transversion into
        transversion_options <- transversion_matches[[base_to_mutate]]
        
        # sample from these two bases according to bg substitution probs
        newbase <- sample(transversion_options, size = 1, prob = bg_sub_prob_mat[base_to_mutate, transversion_options])
        
        # now have to replace the existing value with the new value that has the transversion-in-insertion
        
        this_base <- 10**(-1*(rand_exp - 1))*(newbase - base_to_mutate)
        
        # this difference will give the sum of what has to be added to go from old base to new base
        return(this_base)
        
      }
      
      else if(existing_mat_vals[x] == -1){ # if a deletion has already occurred here
        return(0) # we won't add anything to the mutation matrix (zero here means that when we sum, it'll be ok)
      }
      else{ # if the base is a point mutation
       
        this_base <- bases_going_to[existing_mat_vals[x]] - existing_mat_vals[x]
        
        return(this_base) 
       
      }
      
    })
    
    
    new_muts <- sparseMatrix(i = i_coords, j = j_coords, 
                             x = x_vals, dims = c(num_rows, num_cols))
    
    incoming_mat <- incoming_mat + new_muts
    return(incoming_mat)
  }
  
  
  uniform_res <- get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
                                          bg_pos_er_list = bg_transversion_pos_er_list,
                                          sample_transversion = TRUE,
                                          mut_type = 'transversion',
                                          verbose = verbose)
  
  num_transversions <- uniform_res[['num_edits']]

  transversion_i_coords <- uniform_res[['i_coords']]
  transversion_j_coords <- uniform_res[['j_coords']]
  
  if(length(uniform_res) == 4){
    going_to_bases <- uniform_res[['selected_bases_to']]  
  }
  
  
  # if some edits occurred, perform necessary transversions
  # we allow the user to force transversions in the non-uniform editing but not the uniform
  # time_num_tranversions == TRUE when some mutation coordinates were actually generated
  if(num_transversions != 0){
    mut_mat <- post_indices_transversion_func(i_coords = transversion_i_coords,
                                              j_coords = transversion_j_coords,
                                              incoming_mat = mut_mat,
                                              force_transversions = FALSE,
                                              bases_going_to = going_to_bases)
  }
  
  # if the transversion pos er list is null or has length zero, there is no non-uniform editing
  if(is.null(target_transversion_pos_er_list)){ # if uniform, we are done after this one step
    return(mut_mat)
  }
  if(length(target_transversion_pos_er_list) == 0){
    return(mut_mat)
  }
    
  elig_ints_list <<- list()
  
  for(pos in names(target_transversion_pos_er_list)){
    # pos_er_list has form list(pos_num1 = rate1, pos_num2 = rate2, ...)
    # iterate through each position that is High/Medium/Low, and see which integrations haven't been edited at that position yet
    unedited_rowvals <- which(mut_mat[, as.integer(pos)] == 0)
    
    # create a new list with names == column (genomic) position, and values == non-edited integrations
    elig_ints_list[[as.character(pos)]] <- unedited_rowvals
  }
  
  if(close_window){
    elig_ints_list <- filter_elig_ints_by_edit_window(pos_to_window_inds_list = target_to_window_ind_list,
                                                      window_to_pos_inds_list = window_to_target_ind_list,
                                                      pos_to_unedited_int_list = elig_ints_list)
  }
  
  nu_res <- non_uniform_editing(pos_er_list = target_transversion_pos_er_list, 
                                eligible_ints = elig_ints_list,
                                num_integrations = NULL)
  nu_transversion_i_coords <- nu_res[['i_coords']]
  nu_transversion_j_coords <- nu_res[['j_coords']]
  
  
  if(nu_transversion_i_coords[1] != FALSE){
    target_bases_going_to <- vapply(seq_len(num_cols), function(pos){
      from_base <- as.integer(baseline_ints[pos])
      if(force_target_transversions){
        if(!is.null(target_transversion_to_base)){
          return(as.integer(target_transversion_to_base))
        }
        return(as.integer(forced_transversion_matches[[from_base]]))
      }
      options <- transversion_matches[[from_base]]
      probs <- bg_sub_prob_mat[from_base, options]
      if(sum(probs) <= 0){
        return(as.integer(options[1]))
      }
      as.integer(sample(options, size = 1, prob = probs))
    }, integer(1))

    # we allow the user to force transversions in the non-uniform editing but not the uniform
    mut_mat <- post_indices_transversion_func(i_coords = nu_transversion_i_coords,
                                              j_coords = nu_transversion_j_coords,
                                              incoming_mat = mut_mat,
                                              bases_going_to = target_bases_going_to,
                                              force_transversions = force_target_transversions)
  }
  
  return(mut_mat)
    
  

}
 


insertion_func <- function(mut_mat, num_rows, num_cols, bg_ins_pos_er_list,
                           target_ins_pos_er_list = NULL,
                           run_id = NULL,
                           this_cell_num = NULL,
                           verbose = FALSE,
                           prime_editing = FALSE,
                           ind_to_prime_seq_int_map = NULL,
                           target_to_window_ind_list = NULL, 
                           window_to_target_ind_list = NULL, 
                           close_window = FALSE
                           # pos_er_nuc_list = NULL
                           ){
  # we accept both the BE pos er list AND nuc pos er list in insertions because there may be elevated substitution rates at cut sites
  
  # function for generating necessary mutation values to do insertions within insertions:
  ins_in_ins <- function(ins_pos, current_ins, new_ins_length){
    
    # we can work backwards, knowing what the insertion-in-insertion should look like by the end,
    # to find the value that has to be added to the current insertion to get the correct insertion-in-insertion
    # because recall we are returning what has to be added to incoming to get new mutated
    
    
    if(current_ins > 0){
      current_ins_length <- nchar(current_ins) - 2  
    }
    else if(current_ins < 0){
      current_ins_length <- nchar(current_ins) - 3  
    }
    
    new_ins <- as.numeric(paste(sample(seq(1,4), new_ins_length, replace = TRUE), collapse = ''))
    
    # extract the integer representation of the values found after the insertion position
    bases_after_insertion <- abs(round(current_ins * 10**(current_ins_length))) %% 10**(current_ins_length - ins_pos + 1)
    
    # shift these bases accounting for the length of the current and new insertions
    shifted_bases_after_insertion <- bases_after_insertion*10**(-(current_ins_length + new_ins_length))
    
    # now find the values before the insertion position
    # this includes the value to the left of the decimal
    bases_before_insertion <- abs(round(current_ins * 10**(ins_pos-1))) %% 10**(ins_pos)
    
    # change the sign of the "before" if the incoming insertion was negative
    
    # shift bases before insertion to that they start right after decimal
    shifted_bases_before_insertion <- bases_before_insertion*10**(-(ins_pos-1))

    # shift the new insertion so that it lines up with the insertion position
    shifted_new_ins <- new_ins*10^(-(ins_pos + new_ins_length - 1))
    
    # find the insertion-in-insertion end result
    # if our current_ins is negative, have to flip to negative since we can only add, not subtract
    if(current_ins > 0){
      final_result <- shifted_bases_before_insertion + shifted_new_ins + shifted_bases_after_insertion  
    }
    else if(current_ins < 0){
      final_result <- -1*(shifted_bases_before_insertion + shifted_new_ins + shifted_bases_after_insertion)
    }
    
    # we need a value that, when added to the incoming insertion, yields the expected new insertion
    new_mut_mat_val <- final_result - current_ins
    
    return(new_mut_mat_val)
  }
  
  
  
  # Accepts sparse matrix as input, and adds to it the insertions that occur
  
  post_indices_insertion_func <- function(i_coords, j_coords, incoming_mat, elig_ints = NULL,
                                          prime_editing = prime_editing,
                                          ind_to_prime_seq_int_map = ind_to_prime_seq_int_map){
    
    # if we have a prime editing system, we already have of interest for a given j coord:
    # note length(icoords) == length(jcoords)
    if(prime_editing){
      
      new_insertion_x_vals <- sapply(X = seq(length(i_coords)), function(x){
        return(as.numeric(paste(c(0, '.', ind_to_prime_seq_int_map[[as.character(j_coords[x])]]), collapse = '')))
      })
    } else{
    
      insertion_lengths <- sapply(rgamma(n = length(i_coords), shape = 1, rate = 1), ceiling)
      ins_pos <- which(incoming_mat %% 1 != 0, arr.ind = TRUE) # which positions in incoming_mat already have insertion?
      
      
      new_insertion_x_vals <- sapply(X = seq(length(i_coords)), function(x){
        
        # check if any of the new insertions occur at positions where insertions already exist
        if(any((ins_pos[,1] == i_coords[x]) & (ins_pos[,2] == j_coords[x])) == TRUE){ # if insertion already exists at this location
          old_insertion <- incoming_mat[i_coords[x], j_coords[x]]
          
          
          
          
          if(old_insertion > 0){
            curr_ins_length <- nchar(old_insertion) - 2
          }
          else if(old_insertion < 0){
            curr_ins_length <- nchar(old_insertion) - 3
          }
          new_ins_pos <- sample(seq(1, curr_ins_length+1), size = 1)
          # this enables insertions-in-insertions
          return(ins_in_ins(ins_pos = new_ins_pos, current_ins = old_insertion, new_ins_length = insertion_lengths[x]))
          
        }
        else{ # if this position doesn't yet have an insertion
          
          # only return the decimal, no need to shift
          return(as.numeric(paste(c(0, '.', sample(seq(1,4), insertion_lengths[x], replace = TRUE)), collapse = '')))
          
        }
      })
    }
    
    incoming_mat <- incoming_mat + sparseMatrix(i = i_coords, j = j_coords,
                                                x = new_insertion_x_vals, dims = c(num_rows, num_cols))
    
    return(incoming_mat)
  }
  
  uniform_res <- get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
                                          bg_pos_er_list = bg_ins_pos_er_list,
                                          mut_type = 'insertion',
                                          verbose = verbose)
  num_insertions <- uniform_res[['num_edits']]
  insertion_i_coords <- uniform_res[['i_coords']]
  insertion_j_coords <- uniform_res[['j_coords']]
  
  
  # get_background_edit_inds() will return FALSE if no edits have occurred
  if(num_insertions != 0){
    mut_mat <- post_indices_insertion_func(i_coords = insertion_i_coords,
                                           j_coords = insertion_j_coords,
                                           incoming_mat = mut_mat,
                                           prime_editing = FALSE,
                                           ind_to_prime_seq_int_map = NULL)
  }
  
  
  # if the insertion pos er list is null or has length zero, there is no non-uniform editing
  if(is.null(target_ins_pos_er_list)){ # if uniform, we are done after this one step
    return(mut_mat)
  }
  if(length(target_ins_pos_er_list) == 0){
    return(mut_mat)
  }
  
  elig_ints_list <<- list()
  
  len1_pos <- c()
  
  for(pos in names(target_ins_pos_er_list)){
    # iterate through each position that is High/Medium/Low, and see which integrations haven't been edited at that position yet
    unedited_rowvals <- which(mut_mat[, as.integer(pos)] == 0)
    elig_ints_list[[as.character(pos)]] <- unedited_rowvals
    
    if(length(unedited_rowvals) == 1){
      len1_pos <- append(len1_pos, as.character(pos))
    }
  }
  
  if(close_window){
    elig_ints_list <- filter_elig_ints_by_edit_window(pos_to_window_inds_list = target_to_window_ind_list,
                                                      window_to_pos_inds_list = window_to_target_ind_list,
                                                      pos_to_unedited_int_list = elig_ints_list)
  }
  
  nu_res <- non_uniform_editing(pos_er_list = target_ins_pos_er_list,
                                # mutation_type = 'Insertion', 
                                num_integrations = NULL,
                                eligible_ints = elig_ints_list,
                                length1_positions = len1_pos)
  nu_insertion_i_coords <- nu_res[['i_coords']]
  nu_insertion_j_coords <- nu_res[['j_coords']]
 
  if(nu_insertion_i_coords[1] != FALSE){
    mut_mat <- post_indices_insertion_func(i_coords = nu_insertion_i_coords,
                                           j_coords = nu_insertion_j_coords,
                                           incoming_mat = mut_mat,
                                           elig_ints = elig_ints_list,
                                           prime_editing = prime_editing,
                                           ind_to_prime_seq_int_map = ind_to_prime_seq_int_map)
  }
 
  return(mut_mat)
  
}




num_deletable_bases <- function(x){
  
  post <- x %% 1
  
  if(x >= 0){
    if(post != 0){ # if there is an insertion here
      del_bases <- nchar(post) - 1 # subtract 1 to account for decimal
    }
    else{ # if there is no deletion or insertion here
      del_bases <- 1
    }}
  else{ # if the base to the left of the decimal was deleted already
    if(post != 0){ # if there is an insertion here still
      del_bases <- nchar(post) - 2 # now have to subtract 1 for decimal and 1 for ineligible already-deleted base
    }
    else{ # if the base has already been deleted and there's no insertion, no more deletions can occur
      del_bases <- 0
    }
  }
  return(del_bases)
}



perform_deletion <- function(ival, jval, del_length, mat_name, num_cols){
  
  # if the deletion has taken us out of bounds
  if((jval > num_cols) | (jval <= 0)){
    return(mat_name)
  }
  old_val <- mat_name[ival, jval]
  deletable_here <- num_deletable_bases(old_val)
  
  if(del_length < deletable_here){ # if we can't delete every base at this position, delete as many as del_length allows
    mat_name[ival, jval] <- round(old_val, digits = nchar(old_val) - del_length - 2) 
  } else if(del_length == deletable_here){ # if we have an exact match, it's easy because we just convert to -1
    mat_name[ival, jval] <- -1L
  } else if(del_length > deletable_here){ # if the deletion has length longer than number of bases we can delete at this position
    
    # this looks right, since by setting mat_name[ival, jval] <- -1, we are deleting deletable_here bases  
    mat_name[ival, jval] <- -1L
    jval <- jval - 1 # we extend the deletion to the left for simplicity
    return(perform_deletion(ival, jval, del_length - deletable_here, mat_name, num_cols))
  }
  return(mat_name)
  
}

all_deletions_one_mat <- function(i, j, d, old_mat, num_cols){
  
  for(elem_num in seq_along(i)){
    old_mat <- perform_deletion(i[elem_num], j[elem_num], d[elem_num], old_mat, num_cols)
  }
  
  return(old_mat)
  
}

deletion_func <- function(mut_mat, num_rows, num_cols, bg_del_pos_er_list, uniform = TRUE,
                          target_del_pos_er_list = NULL, verbose = FALSE, interdeletion_dropout_prob = 0,
                          interdeletion_dropout_radius = 0,
                          target_to_window_ind_list = NULL, 
                          window_to_target_ind_list = NULL, 
                          close_window = FALSE){
  
  # we accept both the BE pos er list AND nuc pos er list in deletions because there may be elevated substitution rates at cut sites

  # Accepts sparse matrix as input, and adds to it the deletions that occur
  multi_edit_bc_dropout <- function(deletion_mut_mat, deletion_radius, dropout_prob, bc_profile){
    
    find_deletions_within_radius <- function(delmat, radius) {
      deletion_positions <- delmat[, 2]
      if(length(deletion_positions) < 2){ # no inter-target dropout if fewer than 2 deletions
        return(matrix(NA, nrow = 0, ncol = 2))
      }
      pairwise_deletion_locs <- t(combn(deletion_positions, 2))  # pairwise deletion location positions

      inter_target_dropout_pairs <- matrix(pairwise_deletion_locs[abs(pairwise_deletion_locs[, 1] - pairwise_deletion_locs[, 2]) <= radius, ],
                                           ncol = 2)
      return(inter_target_dropout_pairs)
    }

    # get integrations nums present
    unique_integrations <- unique(deletion_mut_mat[, 1])

    dropped_out_intervening_positions <- lapply(unique_integrations, function(int_num) {

      # get deletion positions for this integration
      int_delmat <- matrix(deletion_mut_mat[deletion_mut_mat[, 1] == int_num, ], ncol = 2)

      # find if/which deletion events occurred within deletion_radius of one another
      dropout_pairs <- find_deletions_within_radius(int_delmat, deletion_radius)

      # write these pairs to a matrix
      if(nrow(dropout_pairs) > 0){
        cbind(int_num, dropout_pairs)
      } else{
        NULL
      }
    })

    num_null <- sum(sapply(dropped_out_intervening_positions, is.null))

    # if not all elements of this list are null, i.e. there is at least one dropout event
    if(num_null != length(dropped_out_intervening_positions)){
      
      for(i in 1:length(dropped_out_intervening_positions)){
        sorted_cols <- sort(c(dropped_out_intervening_positions[[i]][2], dropped_out_intervening_positions[[i]][3]))
      
      }
      # stack results across integration numbers
      result_matrix <- do.call(rbind, dropped_out_intervening_positions)




      res <- apply(result_matrix, MARGIN = 1, function(row){
        # for each pair of muts that could have dropped out, probabilistically determine if dropout occurred:
        dropout_occurs <- rbinom(n = 1, size = 1, prob = dropout_prob)
        sorted_row_inds <- sort(c(row[2], row[3]))

        # only set intervening seqs to -1 if dropout occurs
        # first position is integration, second is intertarget dropout start, third is intertarget dropout end
        if(dropout_occurs){
          
          bc_profile[as.integer(row[1]), sorted_row_inds[1]:sorted_row_inds[2]] <<- -1L # global update to bc mutmat
        }

      })
      
    }
    
    return(bc_profile)
  }

  
  post_indices_deletion_func <- function(i_coords, j_coords, incoming_mat){


    deletion_lengths <- sapply(rgamma(n = length(i_coords), shape = 1, rate = 1), ceiling)

    return(all_deletions_one_mat(i = i_coords, j = j_coords, d = deletion_lengths,
                                 old_mat = incoming_mat, num_cols = num_cols))
  }
  
  uniform_res <- get_background_edit_inds(num_rows = num_rows, num_cols = num_cols,
                                          bg_pos_er_list = bg_del_pos_er_list,
                                          mut_type = 'deletion',
                                          verbose = verbose)
  num_deletions <- uniform_res[['num_edits']]
  start_i_coords <- uniform_res[['i_coords']]
  start_j_coords <- uniform_res[['j_coords']]
  

  if(num_deletions != 0){
    mut_mat <- post_indices_deletion_func(i_coords = start_i_coords,
                                          j_coords = start_j_coords,
                                          incoming_mat = mut_mat)

  
    # if we are allowing dropout of intervening barcode seqs due to >=2 simultaneous deletions:
    if((interdeletion_dropout_prob > 0) & (interdeletion_dropout_radius > 0)){
      # make a matrix of the starting i and j coords of the new deletions:
      deletion_ijs_this_timepoint <- cbind(start_i_coords, start_j_coords)
      
      mut_mat <- multi_edit_bc_dropout(deletion_mut_mat = deletion_ijs_this_timepoint,
                            bc_profile = mut_mat,
                            deletion_radius = interdeletion_dropout_radius,
                            dropout_prob = interdeletion_dropout_prob)
    }
  }
  

  # if the deletion pos er list is null or has length zero, there is no non-uniform editing
  if(is.null(target_del_pos_er_list)){ # if uniform, we are done after this one step
    return(mut_mat)
  }
  if(length(target_del_pos_er_list) == 0){
    return(mut_mat)
  }


  elig_ints_list <<- list()
  
  for(pos in names(target_del_pos_er_list)){
    # iterate through each position that is High/Medium/Low, and see which integrations haven't been edited at that position yet
    
    unedited_rowvals <- which(mut_mat[, as.integer(pos)] == 0)
    elig_ints_list[[as.character(pos)]] <- unedited_rowvals
  }
  
  if(close_window){
    elig_ints_list <- filter_elig_ints_by_edit_window(pos_to_window_inds_list = target_to_window_ind_list,
                                                      window_to_pos_inds_list = window_to_target_ind_list,
                                                      pos_to_unedited_int_list = elig_ints_list)
  }
  

  nu_res <- non_uniform_editing(pos_er_list = target_del_pos_er_list,
                                num_integrations = NULL,
                                eligible_ints = elig_ints_list)
  nu_deletion_i_coords <- nu_res[['i_coords']]
  nu_deletion_j_coords <- nu_res[['j_coords']]
  
 
  if(nu_deletion_i_coords[1] != FALSE){
    mut_mat <- post_indices_deletion_func(i_coords = nu_deletion_i_coords,
                                          j_coords = nu_deletion_j_coords,
                                          incoming_mat = mut_mat)
    # if we are allowing dropout of intervening barcode seqs due to >=2 simultaneous deletions:
    if((interdeletion_dropout_prob > 0) & (interdeletion_dropout_radius > 0)){

      # make a matrix of the starting i and j coords of the new deletions:
      deletion_ijs_this_timepoint <- cbind(nu_deletion_i_coords, nu_deletion_j_coords)

      mut_mat <- multi_edit_bc_dropout(deletion_mut_mat = deletion_ijs_this_timepoint,
                            bc_profile = mut_mat,
                            deletion_radius = interdeletion_dropout_radius,
                            dropout_prob = interdeletion_dropout_prob)
    }
  }

  
  return(mut_mat)




}


perform_all_mt_mutations <- function(incoming_mut_mat,
                                     bg_transition_list,
                                     bg_transversion_list,
                                     bg_insertion_list,
                                     bg_deletion_list,
                                     prob_sub_mat
                                     ){
  

  nonzero <- which(incoming_mut_mat != 0, arr.ind = TRUE)
  
  incoming_mut_mat <- transition_func(mut_mat = incoming_mut_mat, 
                                      # num_rows = num_rows_mt, 
                                      num_rows = dim(incoming_mut_mat)[1],
                                      num_cols = num_cols_mt, 
                                      # uniform_transition_prob = transition_prob_mt,
                                      bg_transition_pos_er_list = bg_transition_list,
                                      baseline_ints = baseline_seq_ints_mt,
                                      verbose = FALSE)
  
  
  incoming_mut_mat <- transversion_func(mut_mat = incoming_mut_mat, 
                                        # num_rows = num_rows_mt, 
                                        num_rows = dim(incoming_mut_mat)[1],
                                        num_cols = num_cols_mt,
                                        # uniform_transversion_prob = transversion_prob_mt,
                                        bg_transversion_pos_er_list = bg_transversion_list,
                                        baseline_ints = baseline_seq_ints_mt,
                                        bg_sub_prob_mat = prob_sub_mat,
                                        verbose = FALSE)
  
  incoming_mut_mat <- insertion_func(mut_mat = incoming_mut_mat, 
                                     # num_rows = num_rows_mt, 
                                     num_rows = dim(incoming_mut_mat)[1],
                                     num_cols = num_cols_mt, 
                                     bg_ins_pos_er_list = bg_insertion_list,
                                     verbose = FALSE,
                                     prime_editing = FALSE,
                                     ind_to_prime_seq_int_map = NULL)
  
  incoming_mut_mat <- deletion_func(mut_mat = incoming_mut_mat, 
                                    # num_rows = num_rows_mt, 
                                    num_rows = dim(incoming_mut_mat)[1],
                                    num_cols = num_cols_mt, 
                                    bg_del_pos_er_list = bg_deletion_list,
                                    verbose = FALSE)
  
  return(incoming_mut_mat)
  
}

perform_all_bc_mutations <- function(incoming_mut_mat, 
                                     bg_transition_list,
                                     bg_transversion_list,
                                     bg_insertion_list,
                                     bg_deletion_list,
                                     target_transition_list,
                                     target_transversion_list,
                                     target_insertion_list,
                                     target_deletion_list,
                                     prob_sub_mat,
                                     timepoint_for_label = '',
                                     urid = NULL,
                                     cell_num = NULL,
                                     interdel_dropout_radius,
                                     interdel_dropout_prob,
                                     prime_editing_system,
                                     ind_to_prime_seq_int_map,
                                     force_target_transversions,
                                     target_transversion_to_base,
                                     close_nuc_window_after_edit,
                                     close_transition_window_after_edit,
                                     close_transversion_window_after_edit,
                                     be_target_to_window_ind_list,
                                     nuc_target_to_window_ind_list,
                                     be_window_to_target_ind_list,
                                     nuc_window_to_target_ind_list
                                     ){
  
  # if we are working with uniform mutation rates for all barcode sequences
  incoming_mut_mat <- transition_func(mut_mat = incoming_mut_mat, 
                                      num_rows = num_rows_bc, 
                                      num_cols = num_cols_bc, 
                                      baseline_ints = baseline_seq_ints_bc,
                                      bg_transition_pos_er_list = bg_transition_list, 
                                      target_transition_pos_er_list = target_transition_list,
                                      timepoint_filename = timepoint_for_label,
                                      target_to_window_ind_list = be_target_to_window_ind_list, 
                                      window_to_target_ind_list = be_window_to_target_ind_list, 
                                      close_window = close_transition_window_after_edit)

  incoming_mut_mat <- transversion_func(mut_mat = incoming_mut_mat, 
                                        num_rows = num_rows_bc, 
                                        num_cols = num_cols_bc,
                                        baseline_ints = baseline_seq_ints_bc,
                                        bg_transversion_pos_er_list = bg_transversion_list, 
                                        target_transversion_pos_er_list = target_transversion_list,
                                        bg_sub_prob_mat = prob_sub_mat,
                                        force_target_transversions = force_target_transversions,
                                        target_transversion_to_base = target_transversion_to_base,
                                        target_to_window_ind_list = be_target_to_window_ind_list, 
                                        window_to_target_ind_list = be_window_to_target_ind_list, 
                                        close_window = close_transversion_window_after_edit)


  incoming_mut_mat <- insertion_func(mut_mat = incoming_mut_mat, 
                                     num_rows = num_rows_bc, 
                                     num_cols = num_cols_bc,
                                     bg_ins_pos_er_list = bg_insertion_list, 
                                     target_ins_pos_er_list = target_insertion_list,
                                     run_id = urid,
                                     this_cell_num = cell_num,
                                     prime_editing = prime_editing_system,
                                     ind_to_prime_seq_int_map = ind_to_prime_seq_int_map,
                                     target_to_window_ind_list = nuc_target_to_window_ind_list, 
                                     window_to_target_ind_list = nuc_window_to_target_ind_list, 
                                     close_window = close_nuc_window_after_edit)
  

  incoming_mut_mat <- deletion_func(mut_mat = incoming_mut_mat, 
                                    num_rows = num_rows_bc, 
                                    num_cols = num_cols_bc, 
                                    bg_del_pos_er_list = bg_deletion_list, 
                                    target_del_pos_er_list = target_deletion_list,
                                    interdeletion_dropout_prob = interdel_dropout_prob,
                                    interdeletion_dropout_radius = interdel_dropout_radius,
                                    target_to_window_ind_list = nuc_target_to_window_ind_list, 
                                    window_to_target_ind_list = nuc_window_to_target_ind_list, 
                                    close_window = close_nuc_window_after_edit)  
  
  
  return(incoming_mut_mat)
  
}
