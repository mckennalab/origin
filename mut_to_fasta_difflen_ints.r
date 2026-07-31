get_profiles_ints_and_umis <- function(cell_pop,
                                       bc_or_mt,
                                       int_rec_prob = NULL,
                                       num_ints = NULL,
                                       umis = NULL){
  
  profiles_ints_and_umis <- lapply(cell_pop, function(cell){
    
    # for mt, num_ints will be null. so we assign to the total number of rows in this matrix
    if(is.null(num_ints)){
      num_ints <- dim(cell$incoming_mt_profiles)[1]
      # print(paste0('num_ints was null, now its ', num_ints))
    }
    
    num_ints_recovered <- rbinom(
      n = 1,
      size = num_ints,
      prob = int_rec_prob
    )
    which_ints_recovered <- if(num_ints_recovered > 0){
      sort(sample.int(num_ints, size = num_ints_recovered, replace = FALSE))
    } else{
      integer(0)
    }
    
    if(bc_or_mt == 'bc'){
      
      mut_mat <- cell$incoming_bc_profiles[which_ints_recovered, , drop = FALSE]
    } else if(bc_or_mt == 'mt'){
      mut_mat <- cell$incoming_mt_profiles[which_ints_recovered, , drop = FALSE]
    } else{
      stop("bc_or_mt must be either 'bc' or 'mt'.")
    }
    
    return_list <- list()
    return_list[['mut_mat']] <- mut_mat
    return_list[['which_ints_recovered']] <- which_ints_recovered
    
    if(!is.null(umis)){
      # subset the bc integration umis to only include those corresponding to selected-for integrations
      these_bc_int_umis <- umis[which_ints_recovered]  
      return_list[['recovered_umis']] <- these_bc_int_umis
    }
    
    
    return(return_list)
  })
  
  return(profiles_ints_and_umis)
  
  
}



# might just get rid of fix_length

ins_to_charvec <- function(ins, pos_num, ref_seq, fixed_length){
  nuc_bases <- c('A', 'G', 'C', 'T')
  if(ins > 0){ 
    num_bases <- nchar(ins) - 1 # adjust for decimal point
    
    # differs from sapply in else{} by what is returned if base_int == 0
    ins_bases <- sapply(seq(1, num_bases), function(ins_basenum){
      base_int <- abs(round(ins * 10**(ins_basenum-1))) %% 10
      if(base_int == 0){ # if we have 0.xx, get the base corresponding to 0
        return(ref_seq[pos_num])
      } else{
        return(nuc_bases[base_int])
      }
    })
    
  } else{ # have to adjust for the negative sign too if less than 0
    num_bases <- nchar(ins) - 2
    
    ins_bases <- sapply(seq(1, num_bases), function(ins_basenum){
      base_int <- abs(round(ins * 10**(ins_basenum-1))) %% 10
      
      if(base_int == 0){ # if we have -0.xx, the base corresponding to "-0" is ''
        if(fixed_length){
          
          return('?')
        } else if(!fixed_length){
          return('')
        }
      }
      
      else{
        return(nuc_bases[base_int])
      }
    })
  }
  
  return_list <- list('num_bases_returned' = length(ins_bases),
                      'bases_returned' = ins_bases)
  return(return_list)
}

get_one_cell_sequence <- function(cell_int_mat, ref_seq, these_bc_int_umis = NULL){
  
  # cell_int_mat will be an n_s x l matrix where n_s is the number of downsampled integrations in the cell and l is the barcode length
  
  list_of_seqs <- lapply(seq_len(nrow(cell_int_mat)), function(int_num){
    
    this_int <- cell_int_mat[int_num, ]
    
    bases <- sapply(seq_along(this_int), function(pos_num){
      if(this_int[pos_num] %% 1){ # insertion
        res <- ins_to_charvec(ins = this_int[pos_num],
                              pos_num = pos_num,
                              ref_seq = ref_seq,
                              fixed_length = FALSE)
        
        bases <- res$bases_returned
        collapsed_bases <- paste(bases, collapse = '')
        return(collapsed_bases)
      } else if(this_int[pos_num] == 0){ # no mut
        newbase <- ref_seq[pos_num]
      } else if(this_int[pos_num] == -1){ # deletion
        newbase <- ''
      } else{ # substitution
        nuc_bases <- c('A', 'G', 'C', 'T')
        newbase <- nuc_bases[this_int[pos_num]]
      }
      return(newbase)
    })
    
    int_seq <- paste(bases, collapse = '')
    
    if(!is.null(these_bc_int_umis)){
      this_int_umi <- these_bc_int_umis[int_num]  
      int_seq <- paste0(this_int_umi, int_seq)
    }
    
    
    return(int_seq)
  })
  
  return(list_of_seqs)
}


write_all_cell_sequences <- function(cell_mutmats, reference, output_fasta_name, fasta_type, bc_integration_umis = NULL){
  
  if(!is.null(bc_integration_umis)){

    
    all_cell_seqs <- parLapply(cl = one_cluster, seq_along(cell_mutmats), function(cellnum){
      get_one_cell_sequence(cell_int_mat = cell_mutmats[[cellnum]],
                            these_bc_int_umis = bc_integration_umis[[cellnum]],
                            ref_seq = reference)
    })  
  } else{
    all_cell_seqs <- parLapply(cl = one_cluster, seq_along(cell_mutmats), function(cellnum){
      get_one_cell_sequence(cell_int_mat = cell_mutmats[[cellnum]],
                            ref_seq = reference)
    })
  }
  
  # adjust fasta name to account for terminal cells only (deprecated reason)
  updated_output_fasta_path <- gsub(pattern = '(.*)(\\.fasta)$', replacement = paste0('\\1_', fasta_type, '\\2'), x = output_fasta_name)
  cells_with_sequences <- lengths(all_cell_seqs) > 0
  if(any(!cells_with_sequences)){
    warning(
      sprintf(
        'Omitting %d cells with no recovered sequences from %s.',
        sum(!cells_with_sequences),
        basename(updated_output_fasta_path)
      ),
      call. = FALSE
    )
  }
  all_cell_seqs <- all_cell_seqs[cells_with_sequences]
  sequence_names <- names(cell_mutmats)[cells_with_sequences]
  if(length(all_cell_seqs) == 0){
    writeLines(character(), updated_output_fasta_path)
    return(invisible(updated_output_fasta_path))
  }
  write.fasta(
    all_cell_seqs,
    names = sequence_names,
    file.out = updated_output_fasta_path
  )
  invisible(updated_output_fasta_path)
}
