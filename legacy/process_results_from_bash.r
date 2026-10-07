# process_results_from_bash.r
#
# Final aggregation step of the sim5_code.R pipeline. Collects the per-tree
# Robinson-Foulds result files written by compare_trees_call_from_bash.r for one
# simulation run, recovers each subrun's parameters from the result filename,
# joins them to the flattened JSON parameters of the run, and writes one wide
# results table.
#
# Invoked by bash_wrapper_all_combos.sh from the repository root as:
#   Rscript process_results_from_bash.r -I <run_id>
# Command-line arguments (optparse; both have NULL defaults):
#   -I / --run_id     run id (numeric string) naming the output/ subdirectories
#   -J / --json_path  fallback parameter-file path, used only for RF files whose
#                     recorded parameter path is the literal '.json'. The
#                     wrapper does not pass this flag.
#
# Inputs (paths are relative to the working directory, which must be the repo
# root so that 'output' resolves):
#   output/rf_dist_files/<run_id>/*   two-line files: normalized RF distance on
#                                     line 1, parameter JSON path on line 2
#   the JSON parameter file named on line 2 of each such file
#
# Outputs:
#   output/param_results_files/<run_id>/stacked_results.csv
#   output/param_results_files/<run_id>/rf_heatmap.png (only if the unused
#                                     make_heatmap() helper is called)
#
# Requires stringr, ggplot2, rjson, and optparse.

suppressPackageStartupMessages({
  library(stringr)
  library(ggplot2)
  library(rjson)
  library(optparse)
})

entering_dir <- getwd()
output_dir_stem <- file.path('output')


option_list <- list(
  make_option(c('-I', '--run_id'), type = 'character', default = NULL,
              help = 'run id (numeric string)'), 
  make_option(c('-J', '--json_path'), type = 'character', default = NULL,
              help = 'path to json param file')
)

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)

#' Parse per-subrun parameters out of an RF result filename
#'
#' Result filenames inherit the naming scheme sim5_code.R uses for processed
#' profile lists: `proc_<mt|bc>_list_<n>_ints_RP_<recovery prob>_samp_<type-rate
#' pairs>_<savename>_time_<t>`, with optional `_CD_<T|F>`, `_AF_<threshold>` and
#' `_B_<T|F>` suffixes added by the score-matrix writer. Each field is pulled
#' into a local variable of the same name; the requested column names are then
#' resolved against those locals by `get()`, so a name with no matching field
#' yields NA rather than an error. A column name starting with the modality
#' prefix has its first three characters (`mt_` / `bc_`) removed before lookup.
#'
#' @param file_name Character scalar: the RF result file name, or the mt/bc half
#'   of a joint `J_` file name.
#' @param bc_or_mt Either `'mt'` or `'bc'`; the modality prefix stripped from
#'   `mat_colnames` before matching.
#' @param mat_colnames Character vector of column names to fill, in the order
#'   they should appear in the returned vector.
#' @param recon_method Reconstruction method (`'fasta'` or `'score'`); returned
#'   verbatim for the `recon_method` column.
#' @return A vector the same length as `mat_colnames`, holding the parsed value
#'   for each name and NA where the filename carries no such field. Sampling
#'   fractions are collapsed into one `cell_rec_fracs` string of
#'   `type: rate` entries joined by commas.
extract_subrun_details <- function(file_name, bc_or_mt, mat_colnames, recon_method){
  
  timept <- str_extract(file_name, '(?<=time_)\\d+(\\.\\d+)?')
  
  # max_num_ints immediately preceds _ints_
  max_ints <- as.integer(str_extract(string = file_name, pattern = '\\d+(?=_ints_)'))
  
  # integration recovery prob immediately follows _recprob_
  int_recprob <- as.numeric(str_extract(string = file_name, pattern = '(?<=_RP_).+(?=_samp)'))
  
  # Cell-type/sample-rate pairs form the contiguous name-number block after
  # "_samp_". Parse the block in one pass so each type receives its own rate.
  sampling_suffix <- str_extract(file_name, '(?<=_samp_).*')
  sampling_block <- str_match(
    sampling_suffix,
    '^((?:[^_]+-[0-9.]+)(?:_[^_]+-[0-9.]+)*)'
  )[, 2]
  sampling_matches <- str_match_all(
    ifelse(is.na(sampling_block), '', sampling_block),
    '(?:^|_)([^_]+)-([0-9.]+)'
  )[[1]]

  if(nrow(sampling_matches) == 0){
    growing_cell_type_names <- character()
    growing_cell_type_samps <- numeric()
  } else{
    growing_cell_type_names <- sampling_matches[, 2]
    growing_cell_type_samps <- as.numeric(sampling_matches[, 3])
  }
  
  cell_rec_fracs <- paste(growing_cell_type_names, growing_cell_type_samps, sep = ': ')
  cell_rec_fracs <- paste(cell_rec_fracs, collapse = ', ')
  
  
  # collapsed deletions immediately follows CD_
  # even if this doesn't exist (fasta case), will have NA returned
  condense_deletions <- str_extract(string = file_name, pattern = '(?<=_CD_)[A-Z]')
  
  af_threshold <- str_extract(string = file_name, pattern = '(?<=_AF_)\\d+(\\.\\d+)?')
  binarize <- str_extract(string = file_name, pattern = '(?<=_B_)[A-Z]')
  
  
  return_vec <- c()
  for(i in seq_along(mat_colnames)){
    
    colname <- mat_colnames[i]
    
    # if the matrix column name starts with mt or bc, remove those prefixes for now, since our vars don't have them
    # note they will not always begin with mt or bc (e.g. timept, which is constant across mt/bc)
    if(startsWith(x = colname, prefix = bc_or_mt)){
      colname <- substr(colname, 4, nchar(colname))
    }
    
    # have to be sure that this resets outside the scope of the function (i.e. not pulling previous iterations' vals at current iteration)
    if(exists(colname, inherits = FALSE)){ # inherits = FALSE restricts to local scope
      return_vec[i] <- get(colname)
    } else{
      return_vec[i] <- NA
    }
  }
  
  return(return_vec)
  
}

#' Rename numbered transition-matrix columns to source/target cell-type pairs
#'
#' `unlist()` on the parsed JSON parameters flattens the uninduced and induced
#' cell-type transition matrices to columns whose only distinguishing feature is
#' a trailing index. The cell-type names are recovered from the
#' `cell_type_dict.cell_type_params.<type>.cell_cycle_length` columns, and the
#' trailing index of each transition-matrix column is replaced by the
#' corresponding concatenated source/target pair. The pairs are generated with
#' `t(outer(...))` so their order matches the row-by-row flattening of the JSON
#' matrices.
#'
#' @param flattened_mat One-row matrix whose column names come from `unlist()`
#'   of a parsed parameter list.
#' @return The same matrix with the uninduced and induced transition-matrix
#'   columns renamed. Errors if either block's column count differs from the
#'   number of source/target cell-type combinations.
change_transition_mat_colnames <- function(flattened_mat){
  
  # clean up the uninduced and induced transition matrix labels (everything else is good enough)
  uninduced_colnums <- which(startsWith(x = colnames(flattened_mat), prefix = 'cell_type_dict.uninduced_transition_matrix'))
  induced_colnums <- which(startsWith(x = colnames(flattened_mat), prefix = 'cell_type_dict.induced_transition_matrix'))
  
  cell_types <- sapply(colnames(flattened_mat), function(colname){
    str_extract(string = colname, pattern = '(?<=cell_type_dict\\.cell_type_params\\.).+?(?=\\.cell_cycle_length)')
  })
  cell_types <- unique(cell_types[!is.na(cell_types)])
  
  # correct row-by-row order is preserved using the t()
  cell_type_combos <- as.character(t(outer(cell_types, cell_types, paste0)))

  if(length(uninduced_colnums) != length(cell_type_combos) ||
     length(induced_colnums) != length(cell_type_combos)){
    stop(sprintf(
      paste(
        'Transition-matrix columns do not match the %d expected source/target',
        'cell-type combinations (uninduced=%d, induced=%d).'
      ),
      length(cell_type_combos),
      length(uninduced_colnums),
      length(induced_colnums)
    ))
  }
  
  for(i in seq_along(uninduced_colnums)){
    stripped_colname <- sub('[0-9]+$', '', colnames(flattened_mat)[uninduced_colnums[i]])
    mod_colname <- paste0(stripped_colname, cell_type_combos[i])
    colnames(flattened_mat)[uninduced_colnums[i]] <- mod_colname
  }
  for(i in seq_along(induced_colnums)){
    stripped_colname <- sub('[0-9]+$', '', colnames(flattened_mat)[induced_colnums[i]])
    mod_colname <- paste0(stripped_colname, cell_type_combos[i])
    colnames(flattened_mat)[induced_colnums[i]] <- mod_colname
  }
  
  return(flattened_mat)
  
}

#' Assemble the stacked results table for one simulation run
#'
#' Walks every file in `output/rf_dist_files/<urid>/`. Two filename families are
#' recognized: single-modality results (`fasta_proc*` / `score_proc*`), whose
#' modality is read from the `proc_mt` / `proc_bc` token, and joint mt+bc
#' results (`fasta_J*` / `score_J*`), which are split at `proc_bc` and parsed as
#' two halves that are then merged element-wise, preferring the non-NA mt value.
#' Each result file supplies the normalized RF distance (line 1) and the path of
#' the JSON parameter file that produced it (line 2); accuracy is stored as
#' `1 - rf_dist`. The flattened parameters of that JSON file, with transition
#' matrix columns renamed, are bound onto every row.
#'
#' @param urid Run id (character) naming the `output/` subdirectory to read.
#' @return A data frame with one row per recognized RF result file: the parsed
#'   subrun columns (timepoint, reconstruction method, mt/bc integration,
#'   recovery, sampling, deletion-collapse, allele-fraction and binarization
#'   settings), `json_path`, `rf_dist`, `accuracy`, and one column per flattened
#'   JSON parameter.
#' @section Side effects: Creates `output/param_results_files/<urid>/` when
#'   absent and writes `stacked_results.csv` into it. Unrecognized filenames
#'   raise a warning and are skipped; if nothing is recognized the function
#'   stops.
make_results_df <- function(urid){
  
  # helper func for downstream processing of this dataframe
  # split vals by semicolon (and colon, if necessary for target specs)
  #' Internal: split a compound parameter string into a typed vector
  #'
  #' @param joined_val Character scalar whose entries are separated by `;` or
  #'   `:`; spaces are removed before splitting.
  #' @param outputted_type `'numeric'` (default), `'integer'`/`'int'`, or any
  #'   other value to leave the pieces as character.
  #' @return A vector of the requested type.
  #' @note Retained helper; nothing in the current `make_results_df` body calls
  #'   it.
  split_vals <- function(joined_val, outputted_type = 'numeric'){
    no_spaces <- str_replace_all(joined_val, ' ', '')
    
    # replace colons with semicolons (for target specification, if necessary)
    no_spaces <- str_replace_all(no_spaces, ':', ';')
    fracs <- unlist(str_split(no_spaces, pattern = ';'))
    
    if(outputted_type == 'numeric'){
      fracs <- as.numeric(fracs)
    }
    else if((outputted_type == 'integer') | (outputted_type == 'int')){
      fracs <- as.integer(fracs)
    }
    
    return(fracs)
  }
  
  
  
  if(!dir.exists(file.path(output_dir_stem, 'param_results_files', urid))){
    dir.create(file.path(output_dir_stem, 'param_results_files', urid), recursive = TRUE)
  }
  

  all_rf_files <- list.files(file.path(output_dir_stem, 'rf_dist_files', urid))

  const_subrun_params <- c('timept',
                           'recon_method')
  mt_only_subrun_params <- c('af_threshold',
                             'binarize')
  bc_only_subrun_params <- c()
  shared_subrun_params <- c('max_ints',
                            'int_recprob',
                            'cell_rec_fracs',
                            'condense_deletions')
  full_colnames <- c(const_subrun_params,
                     mt_only_subrun_params,
                     bc_only_subrun_params,
                     paste('mt', shared_subrun_params, sep = '_'),
                     paste('bc', shared_subrun_params, sep = '_'))
  

  # subrun_details_mat <- matrix(NA, nrow = length(all_rf_files), ncol = length(full_colnames)+3)
  subrun_details_list <- list()
  
  for(i in seq_len(length(all_rf_files))){
    
    
    file_name <- all_rf_files[i]
    
    path_to_file <- file.path(output_dir_stem, 'rf_dist_files', urid, file_name)
  
    subrun_details <- NULL

    if((startsWith(x = file_name, prefix = 'fasta_proc')) |
       (startsWith(x = file_name, prefix = 'score_proc'))){
      bc_or_mt <- str_extract(string = file_name,
                                          pattern = '(?<=proc_?)(mt|bc)(?=_list)')
      recon_method <- str_extract(string = file_name,
                              pattern = '(fasta|score)(?=_proc)')
      subrun_details <-  extract_subrun_details(file_name = file_name, 
                                                bc_or_mt = bc_or_mt, 
                                                mat_colnames = full_colnames,
                                                recon_method = recon_method)
      
    } else if(startsWith(x = file_name, prefix = 'fasta_J') |
              (startsWith(x = file_name, prefix = 'score_J'))){
      # split the file name into the mt part (first half) and bc part (second) 
      splits <- str_split(string = file_name, pattern = 'proc_bc')[[1]]
      mt_filename_info <- splits[1]
      bc_filename_info <- splits[2]
      
      recon_method <- str_extract(string = file_name,
                                  pattern = '(fasta|score)(?=_J)')
      
      mt_subrun_details <-  extract_subrun_details(file_name = mt_filename_info, 
                                                bc_or_mt = 'mt', 
                                                mat_colnames = full_colnames,
                                                recon_method = recon_method)
      
      bc_subrun_details <-  extract_subrun_details(file_name = bc_filename_info, 
                                                bc_or_mt = 'bc', 
                                                mat_colnames = full_colnames,
                                                recon_method = recon_method)
      
      # merge/stack these mt and bc subrun details into a single vector
      subrun_details <- ifelse(!is.na(mt_subrun_details), mt_subrun_details, bc_subrun_details)
    }

    if(is.null(subrun_details)){
      warning(sprintf('Skipping unrecognized RF result filename: %s', file_name))
      next
    }
    
    
    first_two_lines <- readLines(path_to_file, n = 2)
    rf_dist <- as.numeric(first_two_lines[1])
    accuracy <- 1-rf_dist
    json_path <- first_two_lines[2]
    if(json_path == '.json'){
      json_path <- input_args$json_path
    }
    

    subrun_details <- append(subrun_details, c(json_path, rf_dist, accuracy))
    subrun_details <- t(as.data.frame(subrun_details))
    colnames(subrun_details) <- c(full_colnames, 'json_path', 'rf_dist', 'accuracy')
    
    # concat these subrun details to the run-specific details
    params_list <- rjson::fromJSON(file = json_path)
    flattened_params_list <- unlist(params_list)
    flattened_params_mat <- matrix(data = flattened_params_list, nrow = 1, dimnames = list(NULL, names(flattened_params_list)))
    flattened_params_df <- as.data.frame(change_transition_mat_colnames(flattened_params_mat))
    
    subrun_and_const_params <- cbind(subrun_details, flattened_params_df)
    
    subrun_details_list[[i]] <- subrun_and_const_params

    
    
    
  }

  if(length(subrun_details_list) == 0){
    stop(sprintf('No recognized RF result files found for run %s.', urid))
  }
  
  subrun_details_mat <- do.call(rbind, subrun_details_list)
  subrun_details_df <- as.data.frame(subrun_details_mat, row.names = NULL)
  
  
  full_output_path <- file.path(output_dir_stem, 'param_results_files', urid, 'stacked_results.csv')
  
  write.csv(subrun_details_df, full_output_path, row.names = FALSE)
  
  print(paste0('wrote results to ', full_output_path))

  return(subrun_details_df)
}
  


make_results_df(input_args$run_id)


  
#' Plot an accuracy heatmap for one run
#'
#' Rebuilds the results table for `run_id` and tiles accuracy over integration
#' count and cell recovery rate, faceted by integration recovery probability and
#' timepoint.
#'
#' @param run_id Run id (character) passed straight to `make_results_df()`.
#' @return The value returned by `ggsave()`; called for its side effect.
#' @section Side effects: Writes
#'   `output/param_results_files/<run_id>/rf_heatmap.png`. Because it calls
#'   `make_results_df()`, it also rewrites that run's `stacked_results.csv`.
#' @note Legacy helper: it is not called from this script's CLI path, and the
#'   columns it references (`num_ints`, `cell_rec_rate`, `int_recovery_prob`,
#'   `timepoint`) are not the ones `make_results_df()` currently emits
#'   (`max_ints`, `cell_rec_fracs`, `int_recprob`, `timept`). `ggsave()` is also
#'   called without `plot = p`, so it saves the last plot rather than `p`.
make_heatmap <- function(run_id){
  res <- make_results_df(run_id)
  
  int_breaks <- unique(res$num_ints)
  p <- ggplot(res, aes(x = num_ints, y = cell_rec_rate, fill = accuracy)) +
    geom_tile() +
    scale_x_continuous(breaks = int_breaks) + 
    geom_text(aes(label = round(accuracy, 2)), color = "black", size = 3) +
    scale_fill_gradient(low = 'white', high = 'darkgreen') +
    theme_bw() + 
    facet_grid(int_recovery_prob ~ timepoint)  # Facet by Group1 (rows) and Group2 (columns)
  
  ggsave(file.path(output_dir_stem, 'param_results_files', run_id, 'rf_heatmap.png'))
}
