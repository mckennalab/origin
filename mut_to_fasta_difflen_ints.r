# Turn integer-encoded mutation profiles into per-cell nucleotide sequences and
# write them out as FASTA. Profiles encode 0 as the reference base, -1 as a
# deletion, 1-4 as a substitution to A/G/C/T, and a decimal as an insertion.
# Deletions drop characters and insertions add them, so records come out at
# different lengths (the "difflen" of the file name) and are not mutually
# aligned. Sourced by sim5_code.R, which supplies seqinr and `one_cluster`.

#' Downsample each cell's integrations or mitochondrial genomes
#'
#' Models sequencing recovery. For every cell an independent
#' `rbinom(1, num_ints, int_rec_prob)` draw decides how many molecules are seen,
#' then that many row indices are sampled without replacement and kept in
#' increasing order. Rows are subset with `drop = FALSE`, so a cell that
#' recovers nothing still yields a well-formed zero-row matrix.
#'
#' @param cell_pop List of cell objects, each carrying `incoming_bc_profiles`
#'   and/or `incoming_mt_profiles` - profiles whose rows are barcode
#'   integrations or mitochondrial genomes and whose columns are sequence
#'   positions.
#' @param bc_or_mt Either `'bc'` to subset `incoming_bc_profiles` or `'mt'` to
#'   subset `incoming_mt_profiles`; any other value is an error.
#' @param int_rec_prob Per-molecule recovery probability, used as the binomial
#'   success probability.
#' @param num_ints Number of molecules available per cell. When `NULL` it is
#'   resolved per cell as `nrow(cell$incoming_mt_profiles)`, which is the
#'   mitochondrial path where each cell carries its own genome count.
#' @param umis Optional vector of UMIs indexed by integration number. When
#'   supplied it is subset by the recovered indices.
#' @return List parallel to `cell_pop`. Each element is a named list with
#'   `mut_mat` (the recovered rows of the profile), `which_ints_recovered` (a
#'   sorted integer vector, possibly `integer(0)`), and, only when `umis` was
#'   given, `recovered_umis`.
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

#' Decode one encoded insertion into its nucleotide characters
#'
#' An insertion is stored as a decimal whose digits after the point are the
#' inserted bases in `1 = A, 2 = G, 3 = C, 4 = T` order, so `0.132` is an
#' insertion of A, C, G. The digit before the point stands for the position's
#' own base: `0` means the reference base is still present (`ins > 0`) or has
#' been deleted (`ins < 0`), while `1`-`4` means it was substituted. Digits are
#' read left to right by shifting `ins` one decimal place at a time and taking
#' the last digit of the rounded result, so the returned vector always leads
#' with the position's own base and continues with the inserted ones.
#'
#' @param ins The encoded value at this position; expected to have a non-zero
#'   fractional part. Its sign says whether the position's own base survives,
#'   and `nchar(ins)` (minus the decimal point, and minus the sign when
#'   negative) fixes how many characters are decoded.
#' @param pos_num Column index of this position, used to look up the reference
#'   base when the leading digit is `0` and `ins` is positive.
#' @param ref_seq Character vector of reference bases, one per column.
#' @param fixed_length When `TRUE`, a base deleted at this position (`ins < 0`
#'   with a leading `0`) is emitted as `'?'` so the column still occupies one
#'   character; when `FALSE` it is emitted as `''` and the column shrinks.
#' @return Named list with `num_bases_returned` (how many elements were
#'   produced, counting an empty-string placeholder) and `bases_returned` (the
#'   character vector itself).
#' @note `get_one_cell_sequence()` always calls this with `fixed_length =
#'   FALSE`, so the `'?'` branch is unused on the current path; the comment
#'   above marks the argument as a candidate for removal.
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

#' Rebuild one cell's recovered molecules as nucleotide strings
#'
#' Walks each row of the cell's profile, converts every column into zero or more
#' characters, and pastes the result into one string per row.
#'
#' @details
#' Per-column length bookkeeping, which is what makes the records
#' variable-length: `0` emits the reference base (one character); a
#' substitution `1`-`4` emits the matching base (one character); `-1` emits
#' `''`, so a deleted base shortens the molecule by one; and an insertion (value
#' with a non-zero fractional part) is expanded by `ins_to_charvec()` with
#' `fixed_length = FALSE`, emitting the position's own base - or nothing, when
#' that base was itself deleted - followed by one character per inserted
#' nucleotide. Nothing is padded, so a molecule's length is the reference length
#' minus its deleted bases plus its inserted bases.
#'
#' @param cell_int_mat Numeric matrix for one cell: `n_s` recovered rows by `l`
#'   sequence positions, holding the 0 / -1 / 1-4 / decimal encoding.
#' @param ref_seq Character vector of the `l` reference bases.
#' @param these_bc_int_umis Optional character vector of UMIs, one per row of
#'   `cell_int_mat`. When supplied, each row's UMI is prepended to its sequence.
#' @return List of length `nrow(cell_int_mat)`, each element a single sequence
#'   string.
#' @note Because deletions remove characters and insertions add them, the
#'   returned strings differ in length both within and between cells, so the
#'   FASTA built from them is unaligned rather than a fixed-width character
#'   matrix.
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


#' Write every cell's reconstructed sequences to a FASTA file
#'
#' Runs `get_one_cell_sequence()` across the cells in parallel, then writes one
#' record per recovered molecule under its cell's name. Cells that recovered
#' nothing contribute no sequence and are dropped with a warning; if that leaves
#' no cells at all, an empty file is written instead of calling `write.fasta()`.
#' Records are emitted at whatever length `get_one_cell_sequence()` produced, so
#' the file is unaligned.
#'
#' @param cell_mutmats Named list of per-cell profile matrices; the names become
#'   the FASTA record names.
#' @param reference Character vector of reference bases, forwarded as `ref_seq`.
#' @param output_fasta_name Destination path ending in `.fasta`; `fasta_type` is
#'   spliced in ahead of the extension.
#' @param fasta_type Short label inserted into the file name, so that
#'   `foo.fasta` with `fasta_type = 'bc'` is written as `foo_bc.fasta`.
#' @param bc_integration_umis Optional list parallel to `cell_mutmats` holding
#'   each cell's UMI vector; when supplied, every sequence is prefixed with its
#'   UMI.
#' @return The path actually written, returned invisibly.
#' @section Side effects:
#'   Writes `<output_fasta_name stem>_<fasta_type>.fasta`, either through
#'   `seqinr::write.fasta` or as an empty file via `writeLines(character(),
#'   ...)` when no cell has a sequence, and raises a warning naming how many
#'   cells were omitted.
#' @note Parallelism uses the cluster object `one_cluster` from the enclosing
#'   script's global environment; it must already exist and have
#'   `get_one_cell_sequence` and `ins_to_charvec` exported to its workers.
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
