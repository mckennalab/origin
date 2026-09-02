# sourcing this for use in sim5_code.R

suppressPackageStartupMessages({
  library(docstring)  
})

 
 
generate_non_be_target_sequence <- function(barcode_length, nuc_fracs, target_from, be_target_count){
  #' @title Generate barcode sequence not including BE targets
  #' @description This function returns a sequence of non-BE-target nucleotides into which
  #' intervening BE target nucleotides are later added.
  #' @return Character vector of length equal to number of non-BE-targets in barcode
  #' @param barcode_length integer. The desired end length of the crispr barcode
  #' @param nuc_fracs numeric. A length 4 numeric vector with relative fractions of c(A, G, C, T) in the barcode
  #' @param target_from character. The nucleotide that is targeted by the base editor (string, 'A', 'G', 'C', or 'T')
  #' @param be_target_count integer. The number of base editing targets in the barcode
  #' @note Specified targets and their respective counts take priority over nucleotide ratios
  
  all_bases <- c('A', 'G', 'C', 'T')

  if(length(barcode_length) != 1 || is.na(barcode_length) ||
     barcode_length < 0 || barcode_length %% 1 != 0){
    stop('barcode_length must be one non-negative integer.')
  }
  if(length(be_target_count) != 1 || is.na(be_target_count) ||
     be_target_count < 0 || be_target_count %% 1 != 0 ||
     be_target_count > barcode_length){
    stop('be_target_count must be an integer between zero and barcode_length.')
  }
  if(length(nuc_fracs) != 4 || any(!is.finite(nuc_fracs)) ||
     any(nuc_fracs < 0) || sum(nuc_fracs) <= 0){
    stop('nuc_fracs must contain four non-negative A/G/C/T weights with a positive sum.')
  }
  if(length(target_from) != 1 || !(target_from %in% all_bases)){
    stop("target_from must be one of 'A', 'G', 'C', or 'T'.")
  }

  # Convert fractional composition into integer counts with the largest-remainder
  # method so the counts always sum exactly to barcode_length.
  normalized_fracs <- nuc_fracs / sum(nuc_fracs)
  raw_counts <- barcode_length * normalized_fracs
  base_counts <- floor(raw_counts)
  names(base_counts) <- all_bases
  remainder <- barcode_length - sum(base_counts)
  if(remainder > 0){
    add_to <- order(raw_counts - base_counts, decreasing = TRUE)[seq_len(remainder)]
    base_counts[add_to] <- base_counts[add_to] + 1
  }

  # Target sites take priority over the requested composition. If the rounded
  # composition contains too few target bases, transfer counts from the most
  # abundant non-target bases without allowing negative counts.
  target_index <- match(target_from, all_bases)
  target_deficit <- max(be_target_count - base_counts[target_index], 0)
  if(target_deficit > 0){
    for(unused in seq_len(target_deficit)){
      donor_counts <- base_counts
      donor_counts[target_index] <- -Inf
      donor_index <- which.max(donor_counts)
      if(base_counts[donor_index] <= 0){
        stop('Unable to reconcile target count with nucleotide composition.')
      }
      base_counts[donor_index] <- base_counts[donor_index] - 1
      base_counts[target_index] <- base_counts[target_index] + 1
    }
  }

  base_counts[target_index] <- base_counts[target_index] - be_target_count
  non_target_sequence <- rep(all_bases, times = base_counts)
  sample(non_target_sequence, size = length(non_target_sequence), replace = FALSE)
}


# create helper function that is used to find indices of either nuc or BE targets in a sequence
# according to their specified configs
generate_target_indices <- function(config, num_targets, target_pos_1, bc_length_with_targets, num_bases_btwn = NULL){
  #' @title Compute the barcode positions that will hold base-editor targets
  #' @description Three layouts are supported, chosen by config (upper-cased
  #' before matching): 'U' spreads num_targets positions evenly from 1 to
  #' bc_length_with_targets, 'R' samples them uniformly at random without
  #' replacement, and 'S' starts at target_pos_1 and steps by num_bases_btwn + 1,
  #' leaving exactly num_bases_btwn non-target bases between consecutive targets.
  #' @return Integer vector of 1-based positions into the finished barcode. 'U'
  #' and 'S' come back in increasing order; 'R' is unsorted.
  #' @param config character. Layout code: 'U' (uniform), 'R' (random), or 'S'
  #' (spaced). Any other value is an error.
  #' @param num_targets integer. Single non-negative integer, no greater than
  #' bc_length_with_targets. Zero short-circuits to integer(0).
  #' @param target_pos_1 integer. Position of the first target. Used only by 'S',
  #' where it must be a positive integer.
  #' @param bc_length_with_targets integer. Length of the finished barcode, i.e.
  #' the range the returned positions index into.
  #' @param num_bases_btwn integer. Number of bases between consecutive targets.
  #' Used only by 'S', where it must be a non-negative integer.
  #' @note 'S' errors when the implied last target, target_pos_1 +
  #' (num_bases_btwn + 1) * (num_targets - 1), would fall past
  #' bc_length_with_targets. Duplicate rounded positions are dropped under 'U',
  #' so the result can in principle be shorter than num_targets.

  # if config is Uniform, we want uniformly-spaced target indices

  config <- toupper(config)
  if(length(num_targets) != 1 || is.na(num_targets) || num_targets < 0 ||
     num_targets %% 1 != 0){
    stop('num_targets must be one non-negative integer.')
  }
  if(length(bc_length_with_targets) != 1 || is.na(bc_length_with_targets) ||
     bc_length_with_targets < 0 || bc_length_with_targets %% 1 != 0){
    stop('bc_length_with_targets must be one non-negative integer.')
  }
  if(num_targets > bc_length_with_targets){
    stop('num_targets cannot exceed bc_length_with_targets.')
  }
  if(num_targets == 0){
    return(integer(0))
  }

  if(config == 'U'){
    
    # maximize the space between successive BE targets, beginning at the first base of the barcode
    all_inds <- unique(sapply(seq(from = 1, to = bc_length_with_targets, length.out = num_targets), round))  
    
  } else if(config == 'R'){ # if config is Random, we want random target indices
    all_inds <- sample(seq(1, bc_length_with_targets), size = num_targets, replace = FALSE)
    
   
  } else if(config == 'S'){
    # if config is Spaced, we have a position of the first BE target as well as an increment
    # such that each subsequent target is increment bases after the first BE target
    if(length(target_pos_1) != 1 || is.na(target_pos_1) ||
       target_pos_1 < 1 || target_pos_1 %% 1 != 0 ||
       length(num_bases_btwn) != 1 || is.na(num_bases_btwn) ||
       num_bases_btwn < 0 || num_bases_btwn %% 1 != 0){
      stop('Spaced targets require a positive integer first position and a non-negative integer gap.')
    }
    final_target_pos <- target_pos_1 + (num_bases_btwn + 1)*(num_targets - 1)
    if(final_target_pos > bc_length_with_targets){
      stop('Incompatible barcode target configuration: final target exceeds barcode length.')
    }
    
    
    all_inds <- unique(sapply(seq(from = target_pos_1, to = final_target_pos, 
                                  by = (num_bases_btwn+1)), round))
  } else{
    stop("config must be one of 'U' (uniform), 'R' (random), or 'S' (spaced).")
  }
  
  return(all_inds)
}

add_intervening_be_targets <- function(target_pos_config, target_from, be_target_count, non_target_sequence, 
                                       first_targ_pos = 1, bases_btwn_targets = 1){
  #' @title Add in base editing targets post-hoc to previously generated non-target sequence
  #' @description BE targets are added according to configuration specified by target_pos_config
  #' @return Character vector of length equal to number of non-BE-targets + num_targets (i.e. final length of barcode)
  #' @param target_pos_config character. See CLAs for details. BE targets can be uniformly (U) or randomly (R) 
  #' distributed throughout the barcode, or spaced (S) with a fixed number of intervening non-bases throughout the barcode
  #' @param target_from character. The nucleotide that is targeted by the base editor ('A', 'G', 'C', or 'T')
  #' @param be_target_count integer. The number of base editing targets in the barcode
  #' @param non_target_sequence character. Character vector generated from generate_non_be_target_sequence() of length 
  #' num_non_be_targets (specified above)
  #' @param first_targ_pos integer. If target_pos_config == 'S', the position (on [1:length(sequence)]) of the first BE target
  #' @param bases_btwn_targets integer. If target_pos_config == 'S', the number of bases between each BE target (e.g. XyyX for bases_btwn_targets = y)
  
  # find number of non-BE-target bases using provided non-target sequence generated above
  num_non_be_targets <- length(non_target_sequence)
  bc_length <- num_non_be_targets + be_target_count
  
  
  all_inds <- generate_target_indices(config = target_pos_config, 
                                      num_targets = be_target_count, 
                                      target_pos_1 = first_targ_pos, 
                                      num_bases_btwn = bases_btwn_targets,
                                      bc_length_with_targets = bc_length)
  
  
  # initialize empty character vector that will store barcode sequence WITH targets
  seq_with_targets <- character(length = length(non_target_sequence) + length(all_inds))
  
  # force target_from at each index in all_inds
  seq_with_targets[all_inds] <- target_from
  

  
  # fill in the remaining non-target positions with the existing sequence
  non_target_inds <- setdiff(seq_along(seq_with_targets), all_inds)
  seq_with_targets[non_target_inds] <- non_target_sequence
  
  
  # we want to return the indices of the targets, as well as the sequence with the targets
  return_list <- list()
  return_list[['target_inds']] <- all_inds
  return_list[['seq_with_targets']] <- seq_with_targets
  
  frac_a <- length(which(seq_with_targets == 'A'))/length(seq_with_targets)
  frac_g <- length(which(seq_with_targets == 'G'))/length(seq_with_targets)
  frac_c <- length(which(seq_with_targets == 'C'))/length(seq_with_targets)
  frac_t <- length(which(seq_with_targets == 'T'))/length(seq_with_targets)
  
  return(return_list)
  
  
}
