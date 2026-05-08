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

extract_subrun_details <- function(file_name, bc_or_mt, mat_colnames, recon_method){
  
  timept <- str_extract(file_name, '(?<=time_)\\d+(\\.\\d+)?')
  
  # max_num_ints immediately preceds _ints_
  max_ints <- as.integer(str_extract(string = file_name, pattern = '\\d+(?=_ints_)'))
  
  # integration recovery prob immediately follows _recprob_
  int_recprob <- as.numeric(str_extract(string = file_name, pattern = '(?<=_RP_).+(?=_samp)'))
  
  # cell type-specific names and sampling rates are found immediately following _samp_
  growing_cell_type_names <- character()
  growing_cell_type_samps <- numeric()
  
  # iteratively find next cell type name by using the previous as a search pattern
  cell_type_pattern <- '(?<=_samp_).+?(?=-)'
  while(TRUE){
    
    # extract this cell type's name
    this_cell_type_name <- str_extract(string = file_name, pattern = cell_type_pattern)
    if(is.na(this_cell_type_name)){
      break
    }
    
    # extract this cell type's sampling fraction
    this_cell_type_samp <- as.numeric(str_extract(string = file_name, pattern = '(?<=-)\\d+(\\.\\d+)?'))
    
    growing_cell_type_names <- append(growing_cell_type_names, this_cell_type_name)
    growing_cell_type_samps <- append(growing_cell_type_samps, this_cell_type_samp)
    
    cell_type_pattern <- paste0('(?<=', this_cell_type_name, '-', this_cell_type_samp, '_).+?(?=-)')
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
  
  for(i in seq_along(uninduced_colnums)){
    stripped_colname <- sub('[0-9]+$', '', colnames(flattened_mat)[uninduced_colnums[i]])
    mod_colname <- paste0(stripped_colname, cell_type_combos[i])
    colnames(flattened_mat)[uninduced_colnums[i]] <- mod_colname
  }
  # BUG: this loop renames `induced_colnums` columns but iterates over
  # `seq_along(uninduced_colnums)`. If the induced and uninduced transition matrices
  # ever produce a different number of flattened columns (e.g. one matrix is sparser
  # in the JSON), `induced_colnums[i]` will go out of bounds (NA index) for trailing
  # columns and the corresponding induced columns silently keep their original
  # `cell_type_dict.induced_transition_matrix.<n>` names. Should be
  # `seq_along(induced_colnums)`.
  for(i in seq_along(uninduced_colnums)){
    stripped_colname <- sub('[0-9]+$', '', colnames(flattened_mat)[induced_colnums[i]])
    mod_colname <- paste0(stripped_colname, cell_type_combos[i])
    colnames(flattened_mat)[induced_colnums[i]] <- mod_colname
  }
  
  return(flattened_mat)
  
}

make_results_df <- function(urid){
  
  # helper func for downstream processing of this dataframe
  # split vals by semicolon (and colon, if necessary for target specs)
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
  
  subrun_details_mat <- do.call(rbind, subrun_details_list)
  subrun_details_df <- as.data.frame(subrun_details_mat, row.names = NULL)
  
  
  full_output_path <- file.path(output_dir_stem, 'param_results_files', urid, 'stacked_results.csv')
  
  write.csv(subrun_details_df, full_output_path, row.names = FALSE)
  
  print(paste0('wrote results to ', full_output_path))
  
}
  


make_results_df(input_args$run_id)


  
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

