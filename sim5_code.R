# sim5_code.R -- main driver for the mitochondrial / barcode lineage simulator.
# Parses the run configuration (command-line flags or a parameters JSON), builds
# the reference sequences and per-cell-type mutation probability tables, runs
# the timestep simulation loop, and writes profiles, trees, FASTAs/score
# matrices and plots under output/. Sources nonuniform_muts_heterogeneous.R and
# the other helper scripts at startup; a "gillespie" simulation_engine
# short-circuits to gillespie_pipeline.R before the plotting stack loads.
print(strrep('#', 60))
options(warn = 1) # print each warning as soon as it is generated
options(error = traceback)  

# ---- Run identifier ----
# generate unique run name: 
# set.seed(42)
unique_run_id <- as.character(sample(1:10000000000000, size = 1))
cat(sprintf('UNIQUE_RUN_ID=%s\n', unique_run_id))

# ---- Early Gillespie engine dispatch ----
# Gillespie runs intentionally dispatch before the historical simulator loads
# its plotting/UI dependency stack. This keeps the exact command-line engine
# usable in lightweight R environments and leaves timestep startup unchanged.
early_arguments <- commandArgs(trailingOnly = TRUE)
early_params_path <- NULL
for(argument_index in seq_along(early_arguments)){
  argument <- early_arguments[argument_index]
  if(argument %in% c('-P', '--params_json_path')){
    if(argument_index < length(early_arguments)){
      early_params_path <- early_arguments[argument_index + 1L]
    }
    break
  }
  if(grepl('^(-P|--params_json_path)=', argument)){
    early_params_path <- sub('^[^=]*=', '', argument)
    break
  }
}
if(!is.null(early_params_path)){
  early_params_path <- normalizePath(early_params_path, mustWork = TRUE)
  early_params <- if(requireNamespace('jsonlite', quietly = TRUE)){
    jsonlite::fromJSON(early_params_path, simplifyVector = FALSE)
  } else if(requireNamespace('rjson', quietly = TRUE)){
    rjson::fromJSON(file = early_params_path)
  } else{
    stop('jsonlite or rjson is required to inspect simulation_engine.')
  }
  early_engine <- early_params$simulation_engine
  if(is.list(early_engine)){
    early_engine <- if(!is.null(early_engine$method)){
      early_engine$method
    } else{
      early_engine$name
    }
  }
  early_engine <- if(is.null(early_engine)){
    'timestep'
  } else{
    tolower(as.character(early_engine))
  }
  if(length(early_engine) != 1L ||
     !(early_engine %in% c('timestep', 'time_step', 'gillespie'))){
    stop('simulation_engine must be timestep or gillespie.')
  }
  if(identical(early_engine, 'gillespie')){
    early_script_argument <- grep(
      '^--file=',
      commandArgs(trailingOnly = FALSE),
      value = TRUE
    )[1]
    early_repo_root <- dirname(normalizePath(
      sub('^--file=', '', early_script_argument),
      mustWork = TRUE
    ))
    setwd(early_repo_root)
    run_spec_dir <- file.path('output', 'run_specs', unique_run_id)
    dir.create(run_spec_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(early_params_path, file.path(
      run_spec_dir,
      basename(early_params_path)
    ))
    source('./prime_editing.R')
    source('./physicell_lineage.R')
    source('./physicell_mito.R')
    source('./ecdna_lineage.R')
    source('./gillespie_lineage.R')
    source('./gillespie_pipeline.R')
    early_configuration <- early_params$gillespie
    if(is.null(early_configuration)){
      early_configuration <- list()
    }
    early_output_dir <- if(is.null(early_configuration$output_dir)){
      file.path('output', 'gillespie', unique_run_id)
    } else{
      as.character(early_configuration$output_dir)
    }
    cat('SIMULATION_ENGINE=gillespie\n')
    tryCatch(
      run_gillespie_lineage_pipeline(
        early_params,
        params_path = early_params_path,
        output_dir = early_output_dir
      ),
      error = function(error){
        message('Gillespie simulation failed: ', conditionMessage(error))
        quit(save = 'no', status = 1L)
      }
    )
    quit(save = 'no', status = 0)
  }
}


# ---- Library loading ----
print('Loading libraries ... ')
suppressPackageStartupMessages({
  library(shiny)
  library(shinyWidgets)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(stringr)
  library(shinyjs)
  library(DT)
  library(entropy)
  library(heatmaply)
  library(plotly)
  library(visNetwork)
  library(ggdendro)
  library(grid)
  library(fresh)
  library(ape)
  library(ggmuller)
  library(phylogram)
  library(data.table)
  library(stringr)
  library(optparse)
  library(rjson)
  library(janitor)
  library(parallel)
  library(kableExtra)
  library(cowplot)
  library(gridExtra)
  library(grid)
  library(ggplotify)
  library(ggpubr)
  library(RColorBrewer)
  library(scales)
  library(randomcoloR)
  library(docstring)
  library(seqinr)
  library(zeallot)
  library(msa)
  library(Biostrings)
  library(Matrix)
  library(profvis)
})


# ---- Command-line option parsing ----
option_list <- list(
    make_option(c('-P', '--params_json_path'), type = 'character', default = NULL,
                help = 'alternative parameter input method: path to simulation parameters json (overwrites any CLAs)')
)

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)


# ---- Working directory and sourced helper scripts ----
# Resolve sourced files and output paths relative to this script, regardless of
# the caller's current working directory.
script_path_arg <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
if(length(script_path_arg) > 0){
  script_path <- sub('^--file=', '', script_path_arg[1])
  setwd(dirname(normalizePath(script_path)))
}

print('Sourcing files ... ')
# fit_plot_parameters.R is a standalone legacy plotting utility. It is not
# sourced here because its calibration CSV is not part of the simulation.
source('./nonuniform_muts_heterogeneous.R')
source('./substitution_models.r')
source('./add_intervening_be_targets_to_seq.r')
source('./prime_editing.R')
# source('./make_babette_tree.r')
source('./mut_to_fasta_difflen_ints.r')
source('./parse_cell_type_specific_args.r')
source('./mut_to_scoremat.r')
print('Files sourced ... ')



# ---- Run log files and sub-run id alphabet ----
# create runlog file
if(!dir.exists(file.path('output', 'run_logs', unique_run_id))){
  dir.create(file.path('output', 'run_logs', unique_run_id), recursive = TRUE)
}
runlog_filename <- paste0('runlog_', unique_run_id, '.txt')
runlog_path <<- file.path('output', 'run_logs', unique_run_id, runlog_filename)
close(file(runlog_path, open = 'w'))
close(file(paste0('./output/run_logs/', unique_run_id, '/no_strings_', unique_run_id, '.txt'), open = 'w'))

# generate letters grid for subrun id generation:
diletters_grid <- expand.grid(LETTERS, LETTERS)
diletters_vec <- apply(diletters_grid, MARGIN = 1, function(x){return(paste0(x[1], x[2]))})
letters_diletters <<- append(LETTERS, diletters_vec)

# ---- Parameter ingestion (JSON or CLAs) and run-spec archive ----
# if the user supplied parameters through a json rather than CLAs, re-write input_args
# create run_specs directory if it doesn't exist:
if(!dir.exists(file.path('output', 'run_specs', unique_run_id))){
  dir.create(file.path('output', 'run_specs', unique_run_id), recursive = TRUE)
}

params_json_dir <- '.'
if(!is.null(input_args$params_json_path)){
  params_json_path <- normalizePath(input_args$params_json_path, mustWork = TRUE)
  params_json_dir <- dirname(params_json_path)
  
  old_json_name_splits <- str_split(string = params_json_path, pattern = '/')[[1]]
  old_json_name <- old_json_name_splits[length(old_json_name_splits)]
  
  # copy this entire json over to the run_specs dir, preserving the json name:
  file.copy(from = params_json_path,
            to = file.path('output', 'run_specs', unique_run_id, old_json_name))
  
  # then read in the args from that json 
  input_args <- fromJSON(file = params_json_path)
  
} else{
  
  
  input_args_mat <- do.call(rbind, input_args)
  param_names <- rownames(input_args_mat)
  rownames(input_args_mat) <- NULL
  input_args_mat <- cbind(param_names, input_args_mat)
  input_param_colnames <- c('input_param', 'val')
  input_args_mat <- rbind(input_param_colnames, input_args_mat)
  input_args_df <- as.data.frame(input_args_mat)
  colnames(input_args_df) <- input_param_colnames
  format_input_args_df <- format.data.frame(input_args_df, justify = 'left')
  
  write.table(format_input_args_df, paste0('./output/run_specs/', unique_run_id, 
                                           '/input_args_', unique_run_id, '.txt'), 
              quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE, sep = '\t')
}

# ---- Random seed and nucleotide encoding tables ----
# extract random seed info
set.seed(input_args$random_seed)


bases <- c(1,2,3,4)
transition_matches <- c(2,1,4,3)
transversion_matches <- list(c(3,4), c(3,4), c(1,2), c(1,2))



# define global force_transversions indicator
force_transversions <<- input_args$force_transversions

# parse through the arguments that can be applied to the same mutational run:
# fractions:
#' Split a semicolon-delimited command-line value into a vector
#'
#' Strips every space from the string, splits on `;`, and optionally coerces the
#' pieces. Any `outputted_type` other than the three recognised names leaves the
#' pieces as characters.
#'
#' @param cla_string One character string of values separated by `;`.
#' @param outputted_type Requested element type: `'numeric'` (the default),
#'   `'integer'`/`'int'`, or anything else to keep characters.
#' @return A vector of the split values in the requested type.
process_cla_string <- function(cla_string, outputted_type = 'numeric'){
  no_spaces <- str_replace_all(cla_string, ' ', '')
  fracs <- unlist(str_split(no_spaces, pattern = ';'))
  
  if(outputted_type == 'numeric'){
    fracs <- as.numeric(fracs)
  }
  else if((outputted_type == 'integer') | (outputted_type == 'int')){
    fracs <- as.integer(fracs)
  }
  
  return(fracs)
}

#' Interpret a target-dispersal configuration string
#'
#' The configuration is identified by its first letter after whitespace removal
#' and upper-casing (`U` uniform, `R` random, `S` spaced). Spaced layouts are
#' written `S:<first_target_position>:<bases_between>`, and those two fields are
#' returned as the raw (uncoerced) character splits.
#'
#' @param config_str Target layout string such as `'U'`, `'R'`, or `'S:1:4'`.
#' @return A list with `config` (the single upper-case letter), `first_targ_pos`
#'   and `bases_btwn`; the latter two are `NA` for any configuration but `S`.
parse_target_config <- function(config_str){
  # function that interprets user's target-dispersal throughout barcode
  
  # empty list which will return values that need to be returned
  return_list <- list()
  
  # remove all whitespace, then capitalize the first letter of the string
  config <- toupper(str_sub(str_replace_all(config_str, ' ', ''), 1, 1))
  
  return_list[['config']] <- config
  
  # downstream processing (first target position and bases between) only necessary for S config
  if(config == 'S'){
    
    splits <- str_split(str_replace_all(config_str, ' ', ''), pattern = ':')[[1]]
    return_list[['first_targ_pos']] <- splits[2]
    return_list[['bases_btwn']] <- splits[3]
  } else{
    
    return_list[['first_targ_pos']] <- NA
    return_list[['bases_btwn']] <- NA  
  }

  return(return_list)
}


#' Split a target count into high, medium, and low edit-rate classes
#'
#' The high/medium/low weights are normalized to sum to one and then turned into
#' integer counts with the largest-remainder method, so no class can go negative
#' and the three counts always add back up to `num_targets`.
#'
#' @param num_targets One non-negative integer; validated.
#' @param class_fracs List with `high`, `medium`, and `low` weights. They must
#'   be finite, non-negative, and sum to a positive value; they need not sum to
#'   one.
#' @return A named list with integer `num_h`, `num_m`, and `num_l`.
parse_target_count_arguments <- function(num_targets, class_fracs){
    
  hml_rates <- as.numeric(c(class_fracs$high,
                            class_fracs$medium,
                            class_fracs$low))

  if(length(num_targets) != 1 || is.na(num_targets) || num_targets < 0 ||
     num_targets %% 1 != 0){
    stop('num_targets must be one non-negative integer.')
  }
  if(length(hml_rates) != 3 || any(!is.finite(hml_rates)) ||
     any(hml_rates < 0) || sum(hml_rates) <= 0){
    stop('High/medium/low target-class weights must be non-negative and sum to a positive value.')
  }

  # Normalize the rates in case they do not sum to 1.
  norm_rates <- hml_rates/sum(hml_rates)

  # Allocate integer class counts with the largest-remainder method so no class
  # can become negative and the counts always sum to num_targets.
  raw_counts <- norm_rates * num_targets
  class_counts <- floor(raw_counts)
  remainder <- num_targets - sum(class_counts)
  if(remainder > 0){
    add_to <- order(raw_counts - class_counts, decreasing = TRUE)[seq_len(remainder)]
    class_counts[add_to] <- class_counts[add_to] + 1
  }
  num_h <- class_counts[1]
  num_m <- class_counts[2]
  num_l <- class_counts[3]

  return_list <- list('num_h' = num_h,
                      'num_m' = num_m,
                      'num_l' = num_l)

  # if the user doesn't specify positions, the returned list will have length 3
  # and will designate the number of H, M, and L targets
  return(return_list)
}

#' Parse a base-editor conversion pattern into its two bases
#'
#' Accepts a pattern such as `'A --> G'`: the string is upper-cased and the
#' single A/C/G/T characters flanking an arrow (one or more dashes then `>`) are
#' extracted.
#'
#' @param example_string Conversion pattern holding two different A/C/G/T bases
#'   separated by an arrow. Malformed patterns and self-conversions raise an
#'   error.
#' @return A list with `from_base` and `to_base` character elements.
parse_be_example <- function(example_string){
  
  # allow the user to specify the mutation that the BE induces (e.g. C --> A)
  # separate the two bases with an arrow (at least one dash)
  
  # convert to caps
  example_string <- toupper(example_string)
  
  # extract the single nucleotides to the left and right of the arrow
  bases <- str_match(string = example_string, pattern = '([ACGT]).*-+>.*([ACGT])')[2:3]
  
  if(length(bases) != 2 || anyNA(bases) || bases[1] == bases[2]){
    stop("be_conversion_pattern must look like 'A --> G' and contain two different A/C/G/T bases.")
  }
  
  return_list <- list()
  return_list[['from_base']] <- bases[1]
  return_list[['to_base']] <- bases[2]
  
  return(return_list)
  
}

 
# note that we can determine if transition or transversion rates should be elevated in targets
# by classifying the type of mutation the BE uses
#' Classify a base-editor conversion as a transition or a transversion
#'
#' @param from_base Single A/C/G/T character the editor converts from.
#' @param to_base Single A/C/G/T character the editor converts to; must differ
#'   from `from_base`.
#' @return `'transition'` for the A/G and C/T pairs, `'transversion'` otherwise.
classify_be_mutation_type <- function(from_base, to_base){
  
  # based on the user-provided BE conversions, classify as transition or transversion
  
  transition_list <- list('C' = 'T',
                          'T' = 'C',
                          'A' = 'G',
                          'G' = 'A')
  
  if(length(from_base) != 1 || length(to_base) != 1 ||
     !(from_base %in% names(transition_list)) ||
     !(to_base %in% names(transition_list)) ||
     from_base == to_base){
    stop('from_base and to_base must be different single A/C/G/T bases.')
  }

  if(transition_list[[from_base]] == to_base){
    return('transition')
  }
  return('transversion') 
  
}
# ---- Base-editor conversion pattern ----
if(!is.null(input_args$be_conversion_pattern)){
  be_target_fromto <- parse_be_example(input_args$be_conversion_pattern)
  be_target_from <- be_target_fromto[['from_base']]
  be_target_to <- be_target_fromto[['to_base']]
  be_mutation_type <- classify_be_mutation_type(from_base = be_target_from,
                                                to_base = be_target_to)  
  be_target_to_int <- as.integer(match(be_target_to, c('A', 'G', 'C', 'T')))
} else{
  be_target_to_int <- NULL
}


# helper function that is used in create_bc_sequence()
#' Randomly assign target indices to high/medium/low edit-rate classes
#'
#' The indices are shuffled first, so taking the three class blocks in order
#' still gives a random assignment.
#'
#' @param inds Integer vector of target positions along the barcode.
#' @param num_h Number of positions to label `'High'`.
#' @param num_m Number of positions to label `'Medium'`.
#' @param num_l Number of positions to label `'Low'`.
#' @return A list whose names are the target positions and whose values are the
#'   class labels `'High'`, `'Medium'`, or `'Low'`.
split_inds_into_hml <- function(inds, num_h, num_m, num_l){
  
  # shuffle target inds in place
  inds <- sample(inds, size = length(inds), replace = FALSE)
  
  # assign HML inds sequentially since target inds are now shuffled 
  ends <- cumsum(c(num_h, num_m, num_l))
  starts <- c(1, head(ends, -1) + 1)
  
  if(num_h > 0){
    high_inds <- inds[starts[1]:ends[1]]  
  } else{
    high_inds <- c()
  }
  
  if(num_m > 0){
    med_inds <- inds[starts[2]:ends[2]]  
  } else{
    med_inds <- c()
  }
  
  if(num_l > 0){
    low_inds <- inds[starts[3]:ends[3]]  
  } else{
    low_inds <- c()
  }
  
  # create a named list
  hml_pos_list <- as.list(c(rep('High', num_h), 
                            rep('Medium', num_m), 
                            rep('Low', num_l)))
  hml_inds <- c(high_inds, med_inds, low_inds)
  
  names(hml_pos_list) <- hml_inds
  
  return(hml_pos_list)
}

# if necessary (i.e. input_args$time_inc == 'auto'), create a greatest common divisor for the cell cycle lengths 
# of all cell types in the cell population. this can be a decimal. if only one cell type is supplied, 
# or if all cell types have the same cell cycle length, 'auto' will increment the simulation by the cell cycle length

#' Greatest common divisor of two possibly non-integer values
#'
#' Both values are made positive and scaled up by the power of ten that clears
#' the decimal places of either argument, reduced by the Euclidean algorithm,
#' then scaled back down. This keeps decimal cell-cycle lengths usable as a
#' simulation time increment.
#'
#' @param a One finite numeric scalar; validated.
#' @param b One finite numeric scalar; validated.
#' @return The greatest common divisor, in the same decimal units as the inputs.
gcd <- function(a, b) {
  
  # since we might be working with decimal vals, we first scale up our vals to ints
  # then we scale back down at the end
  if(length(a) != 1 || length(b) != 1 || !is.finite(a) || !is.finite(b)){
    stop('gcd expects two finite numeric scalars.')
  }

  a_text <- sub('0+$', '', format(abs(a), scientific = FALSE, trim = TRUE, digits = 15))
  b_text <- sub('0+$', '', format(abs(b), scientific = FALSE, trim = TRUE, digits = 15))
  num_decimal_places_a <- if(grepl('.', a_text, fixed = TRUE)){
    nchar(sub('^[^.]*\\.', '', a_text))
  } else 0
  num_decimal_places_b <- if(grepl('.', b_text, fixed = TRUE)){
    nchar(sub('^[^.]*\\.', '', b_text))
  } else 0

  max_decimal_places <- max(c(num_decimal_places_a, num_decimal_places_b))
  scale_factor <- 10**max_decimal_places
  
  a <- round(abs(a)*scale_factor)
  b <- round(abs(b)*scale_factor)
  
  while (b != 0) {
    temp <- b
    b <- a %% b
    a <- temp
  }
  return(a/scale_factor)
}

#' Greatest common divisor across two or more values
#'
#' @param ... Numeric values, or one numeric vector, each accepted by `gcd`.
#' @return One numeric greatest common divisor, obtained by folding `gcd` over
#'   the concatenated values.
gcd_multiple_vals <- function(...){
  vals <- c(...)
  return(Reduce(gcd, vals))
}


  
# ---- Simulation time increment ----
# generate the time_inc of the simulation
if(input_args$time_inc == 'auto'){
  cell_cycle_lengths <- sapply(input_args$cell_type_dict$cell_type_params, function(celltype){
    celltype$cell_cycle_length})
  time_inc <- gcd_multiple_vals(cell_cycle_lengths)
} else{
  time_inc <- as.numeric(input_args$time_inc)
}



# ---- Editing-window closure flags and barcode base composition ----
close_nuc_window_after_edit <- input_args$nuclease_targets$editing_window$close_after_edit
close_be_window_after_edit <- input_args$be_targets$editing_window$close_after_edit

# initialize these to FALSE. if close_be_window_after_edit is TRUE, will rewrite the appropriate one to TRUE
close_transition_window_after_edit <- FALSE
close_transversion_window_after_edit <- FALSE


# convert base fractions from a character string to a numeric vector
barcode_base_fracs <- c(input_args$bc_nuc_composition$frac_a,
                   input_args$bc_nuc_composition$frac_g,
                   input_args$bc_nuc_composition$frac_c,
                   input_args$bc_nuc_composition$frac_t)


# ---- Reference sequence construction ----
# first, create barcode and mt sequences
# create mito and bc sequences in chars and ints
int_to_nuc_list <- list('1' = 'A', '2' = 'G', '3' = 'C', '4' = 'T')
nuc_to_int_list <- setNames(names(int_to_nuc_list), int_to_nuc_list)
#' Convert an integer nucleotide code to its character
#'
#' @param int_val Integer code, where `1`-`4` are `A`, `G`, `C`, `T`.
#' @return The matching single-character nucleotide from `int_to_nuc_list`.
convert_int_to_nuc <- function(int_val){
  return(int_to_nuc_list[[as.character(int_val)]])
}
#' Convert a nucleotide character to its integer code
#'
#' @param nuc_val Single `A`, `G`, `C`, or `T` character.
#' @return The matching entry of `nuc_to_int_list`, which stores the codes as
#'   character strings (`'1'`-`'4'`) rather than integers.
convert_nuc_to_int <- function(nuc_val){
  return(nuc_to_int_list[[nuc_val]])
}

#' Nucleotide composition of a sequence
#'
#' @param sequence Non-empty nucleotide vector. It is upper-cased before
#'   counting and rejected if it holds anything but `A`, `G`, `C`, or `T`.
#' @return A numeric vector of length four, named `A`, `G`, `C`, `T`, giving
#'   each base's fraction of the sequence, including zeros for absent bases.
sequence_nucleotide_fractions <- function(sequence){
  if(length(sequence) == 0){
    stop('Cannot calculate nucleotide fractions for an empty sequence.')
  }
  sequence <- toupper(as.character(sequence))
  if(any(!(sequence %in% c('A', 'G', 'C', 'T')))){
    stop('Sequence contains a value other than A, G, C, or T.')
  }
  counts <- table(factor(sequence, levels = c('A', 'G', 'C', 'T')))
  setNames(as.numeric(counts) / length(sequence), c('A', 'G', 'C', 'T'))
}




# updated way to construct a barcode sequence with targets at the correct positions
#' Build the barcode sequence and its target edit-rate class maps
#'
#' Base-editor targets need the underlying sequence manipulated: a sequence that
#' avoids the editor's source base is generated first, then the targets are
#' written back in at the configured positions. Nuclease targets need only
#' positions. When the caller supplies target positions directly instead of
#' high/medium/low fractions, that position-to-class list is passed through
#' unchanged. When `path_to_bc_seq` is given, the sequence is read from disk and
#' only target positions are generated.
#'
#' @param be_target_origin Single character the base editor converts from;
#'   defaults to the script-level `be_target_from`.
#' @param bc_length Integer barcode length. A supplied sequence must match it.
#' @param be_targets_counts Base-editor target count, or a position-to-class
#'   list; `NULL` when there are no base-editor targets.
#' @param nuc_targets_counts Nuclease target count, or a position-to-class list;
#'   `NULL` when there are no nuclease targets.
#' @param be_targets_classfracs High/medium/low weights for base-editor targets.
#' @param nuc_targets_classfracs High/medium/low weights for nuclease targets.
#' @param be_targets_configs Base-editor target layout string, as understood by
#'   `parse_target_config`.
#' @param nuc_targets_configs Nuclease target layout string.
#' @param path_to_bc_seq Optional path to a fixed barcode sequence, as plain
#'   text or FASTA. Header lines and whitespace are stripped and the remainder
#'   must be A/G/C/T only and of length `bc_length`.
#' @param bc_base_fracs Numeric A, G, C, T fractions used when generating a
#'   sequence.
#' @return A list with `bc_seq` (character vector of nucleotides),
#'   `be_basepos_editrate_classes`, and `nuc_basepos_editrate_classes` (named
#'   position-to-class lists, `NULL` when that target type is absent).
create_bc_sequence <- function(be_target_origin = be_target_from,
                               bc_length = input_args$bc_length, 
                               be_targets_counts = input_args$be_targets$num_targets,
                               nuc_targets_counts = input_args$nuclease_targets$num_targets,
                               be_targets_classfracs = input_args$be_targets$edit_rate_class_fractions,
                               nuc_targets_classfracs = input_args$nuclease_targets$edit_rate_class_fractions,
                               be_targets_configs = input_args$be_targets$config,
                               nuc_targets_configs = input_args$nuclease_targets$config,
                               path_to_bc_seq = input_args$barcode_sequence,
                               bc_base_fracs = barcode_base_fracs){
  
  # be_targets_counts refers to the argument that specifies how many targets there are
  # be_targets_configs refers to how targets are dispersed throughout the barcode, ie Random, Uniform, Spaced
  
  # initialize an empty list that will be returned at the end of this function
  # this list will contain:
  # - the complete bc sequence
  # - a BE list with structure pos:{HML}
  # - a nuc list with structure pos:{HML}
  
  return_list <- list()

  # if there are no barcode targets at all (BE or nuc)
  if(is.null(nuc_targets_counts) & is.null(be_targets_counts)){
    bc_seq <- generate_non_be_target_sequence(
      barcode_length = bc_length,
      nuc_fracs = bc_base_fracs,
      target_from = 'A',
      be_target_count = 0
    )
    return_list[['bc_seq']] <- bc_seq
    
    # there are no edit rate lists to return here, so return NULL
    return_list[['be_basepos_editrate_classes']] <- NULL
    return_list[['nuc_basepos_editrate_classes']] <- NULL
    return(return_list)
  }

  # parse the nuc and BE target info
  if(!is.null(be_targets_counts)){
    be_target_setup <- parse_target_count_arguments(num_targets = be_targets_counts,
                                                    class_fracs = be_targets_classfracs)  
  }
  if(!is.null(nuc_targets_counts)){
    nuc_target_setup <- parse_target_count_arguments(num_targets = nuc_targets_counts,
                                                    class_fracs = nuc_targets_classfracs)  


  }
  
  # if there is a provided be target config
  if(!is.null(be_targets_configs)){
    # now parse the inputted target configurations
    parsed_be_target_config <- parse_target_config(be_targets_configs)
    be_target_config_pattern <- parsed_be_target_config[['config']]
    be_first_target_pos <- as.integer(parsed_be_target_config[['first_targ_pos']])
    be_target_num_bases_btwn <- as.integer(parsed_be_target_config[['bases_btwn']])
    
  }
  
  # if there is a provided nuc target config
  if(!is.null(nuc_targets_configs)){
    parsed_nuc_target_config <- parse_target_config(nuc_targets_configs)
    nuc_target_config_pattern <- parsed_nuc_target_config[['config']]
    nuc_first_target_pos <- as.integer(parsed_nuc_target_config[['first_targ_pos']])
    nuc_target_num_bases_btwn <- as.integer(parsed_nuc_target_config[['bases_btwn']])
  }
  
  # if a barcode sequence file is provided, we don't need to create a sequence
  # we just need to assign indices to targets
  if(!is.null(path_to_bc_seq)){
    # read the sequence straight into the return_list since it already contains all targets
    sequence_lines <- readLines(path_to_bc_seq, warn = FALSE)
    sequence_lines <- sequence_lines[!grepl('^>', sequence_lines)]
    sequence_text <- toupper(gsub('[[:space:]]+', '', paste(sequence_lines, collapse = '')))
    return_list[['bc_seq']] <- str_split(sequence_text, '')[[1]]
    if(length(return_list[['bc_seq']]) != bc_length){
      stop('The supplied barcode sequence length does not match bc_length.')
    }
    if(any(!(return_list[['bc_seq']] %in% c('A', 'G', 'C', 'T')))){
      stop('The supplied barcode sequence must contain only A, G, C, and T.')
    }
    
    # since no manipulation of the underlying sequence is required FOR BE AND FOR NUC, we only need to generate the indices
    
    if(!is.null(be_targets_counts)){
      if('num_h' %in% names(be_target_setup)){
        # print('in num h be target setup == ')
        # print(be_target_setup)
        num_be_targets <- sum(as.numeric(be_target_setup))
        bc_be_target_inds <- generate_target_indices(config = be_target_config_pattern, 
                                                     num_targets = num_be_targets, 
                                                     target_pos_1 = be_first_target_pos, 
                                                     bc_length_with_targets = bc_length, 
                                                     num_bases_btwn = be_target_num_bases_btwn)
        # print('past bc be target inds')
        # print(bc_be_target_inds)
        
        
        return_list[['be_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_be_target_inds, 
                                                                            num_h = be_target_setup[['num_h']], 
                                                                            num_m = be_target_setup[['num_m']], 
                                                                            num_l = be_target_setup[['num_l']])
      } else{
        
        
        return_list[['be_basepos_editrate_classes']] <- be_target_setup
      }  
    }
    
    if(!is.null(nuc_targets_counts)){
      if('num_h' %in% names(nuc_target_setup)){
        num_nuc_targets <- sum(as.numeric(nuc_target_setup))
        
        
        
        bc_nuc_target_inds <- generate_target_indices(config = nuc_target_config_pattern, 
                                                      num_targets = num_nuc_targets, 
                                                      target_pos_1 = nuc_first_target_pos, 
                                                      bc_length_with_targets = bc_length, 
                                                      num_bases_btwn = nuc_target_num_bases_btwn)
        return_list[['nuc_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_nuc_target_inds, 
                                                                             num_h = nuc_target_setup[['num_h']], 
                                                                             num_m = nuc_target_setup[['num_m']], 
                                                                             num_l = nuc_target_setup[['num_l']])
      } else{

        return_list[['nuc_basepos_editrate_classes']] <- nuc_target_setup
      } 
      
      
    }
     

    return(return_list)
    
  }
  
  # (if a barcode path is not provided)
  # if 'num_h' is in the names in the returned list, it means the user
  # did not input position-specific HML edit rates and instead inputted a HML ratio. 
  # this means we still need to manually specify where the targets are located
  # which under the current approach is done by first generating a sequence with no targets
  # then going back in and adding the targets
  if(!is.null(be_targets_counts)){
    if('num_h' %in% names(be_target_setup)){ 
      
      
      
      num_be_targets <- sum(as.numeric(be_target_setup)) 
      
      # the total number of targets is computed by summing the number of HML targets in be_target_setup

      bc_sequence_no_targets <- generate_non_be_target_sequence(barcode_length = bc_length, 
                                                                nuc_fracs = bc_base_fracs,
                                                                target_from = be_target_origin,
                                                                be_target_count = num_be_targets)
      full_seq_return_list <- add_intervening_be_targets(target_pos_config = be_target_config_pattern,
                                                         target_from = be_target_origin,
                                                         be_target_count = num_be_targets,
                                                         first_targ_pos = be_first_target_pos,
                                                         non_target_sequence = bc_sequence_no_targets,
                                                         bases_btwn_targets = be_target_num_bases_btwn)
      bc_sequence_with_targets <- full_seq_return_list[['seq_with_targets']]
      
      bc_be_target_inds <- full_seq_return_list[['target_inds']]
      
      
      return_list[['be_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_be_target_inds, 
                                                                          num_h = be_target_setup[['num_h']], 
                                                                          num_m = be_target_setup[['num_m']], 
                                                                          num_l = be_target_setup[['num_l']])

      
      return_list[['bc_seq']] <- bc_sequence_with_targets

    } else{ # if the user specified where the targets are, we don't need to go through the process of 
      # generating a sequence without BE targets then adding them on
      # can still use generate_non_be_target_sequence() to get the sequence since it considers nuc fractions
      # note that we specify num_be_targets = 0 so that we don't return a truncated sequence here
      bc_sequence_with_targets <- generate_non_be_target_sequence(
        barcode_length = bc_length, 
        nuc_fracs = bc_base_fracs,
        target_from = be_target_origin,
        be_target_count = 0
      )
      
      

      
      return_list[['bc_seq']] <- bc_sequence_with_targets
      return_list[['be_basepos_editrate_classes']] <- be_target_setup
      
    }  
  }
  
  
  # now identify the nuclease targets. this one is simpler since it won't require manipulating the barcode sequence
  # if 'num_h' is in the names in the returned list, it means the user
  # did not input position-specific HML edit rates and instead inputted a HML ratio.
  # so we need to generate the indices of each of the targets according to the specified nuc target config
  if(!is.null(nuc_targets_counts)){
    if('num_h' %in% names(nuc_target_setup)){
      

      num_nuc_targets <- sum(as.numeric(nuc_target_setup))
      
      # since no manipulation of the underlying sequence is required, we only need to generate the indices
      bc_nuc_target_inds <- generate_target_indices(config = nuc_target_config_pattern, 
                                                    num_targets = num_nuc_targets, 
                                                    target_pos_1 = nuc_first_target_pos, 
                                                    bc_length_with_targets = bc_length, 
                                                    num_bases_btwn = nuc_target_num_bases_btwn)

      
      return_list[['nuc_basepos_editrate_classes']] <- split_inds_into_hml(inds = bc_nuc_target_inds, 
                                                                          num_h = nuc_target_setup[['num_h']], 
                                                                          num_m = nuc_target_setup[['num_m']], 
                                                                          num_l = nuc_target_setup[['num_l']])

    } else{ # if actual indices are supplied along with HML edit rate classes

      
      
      
      return_list[['nuc_basepos_editrate_classes']] <- nuc_target_setup
    }
    if(is.null(return_list[['bc_seq']])){
      return_list[['bc_seq']] <- generate_non_be_target_sequence(
        barcode_length = bc_length,
        nuc_fracs = bc_base_fracs,
        target_from = 'A',
        be_target_count = 0
      )
    }
  }

  return(return_list)
  
  
}

# ---- Barcode sequence and target edit-rate classes ----
bc_generation_return_list <- create_bc_sequence()
baseline_seq_nucs_bc <<- bc_generation_return_list[['bc_seq']]
# ERC == edit rate class

# not a crazy assumption that basepos_erc_be_list and basepos_erc_nuc_list would each be the same across cell types
# for example, a target that is high edit rate in one cell type would be high edit rate in another, regardless of whether the underlying numerical rates themselves are different
# thus only have one basepos_erc list for be and for nuc across all cell types
# that said, heterogeneity is introduced within each erc

baseline_seq_ints_bc <<- sapply(baseline_seq_nucs_bc, convert_nuc_to_int) 


# helper func for generating target-specific prime editing sequences
#' Assign prime-editing insertion sequences to target positions
#'
#' A library of `num_unique_guides` random guides is drawn first, then one guide
#' is assigned per target; sampling is with replacement only when there are more
#' targets than unique guides.
#'
#' @param target_inds Integer target positions along the barcode; they become
#'   the names of both returned maps.
#' @param guide_length Number of bases in each guide sequence.
#' @param num_unique_guides Size of the guide library to draw from.
#' @return A list with `ind_to_prime_seq_int_map` (position-keyed integer
#'   nucleotide vectors) and `ind_to_prime_seq_nuc_map` (the same guides as
#'   collapsed A/G/C/T strings).
#' @note Not called by the simulation itself; the run uses the backend built by
#'   `prepare_prime_editing_backend`. It is exercised by the regression tests.
create_prime_editing_basepos_seqs <- function(target_inds, guide_length, num_unique_guides){
  
  # generate an insertion sequence library where each seq has length guide_length
  prime_guide_library <- lapply(seq(num_unique_guides), function(guide_num){
    sample(seq(1, 4), size = guide_length, replace = TRUE)
  })  
  
  # if there are more target indices than unique guides, we sample from the guides with replacement
  # otherwise, we do not sample with replacement
  sample_guides_with_replacement <- length(target_inds) > num_unique_guides
  
  # generate integer representations of insertion seqs for each target
  ind_to_prime_seq_int_map <- sample(prime_guide_library, size = length(target_inds), 
                                     replace = sample_guides_with_replacement)
  names(ind_to_prime_seq_int_map) <- target_inds
  
  # translate integer representations to nucleotide string representations
  ind_to_prime_seq_nuc_map <- lapply(ind_to_prime_seq_int_map,
                                     function(ints){
                                       nuc_vec <- vapply(
                                         ints,
                                         function(int) int_to_nuc_list[[as.character(int)]],
                                         character(1)
                                       )
                                       return(paste(nuc_vec, collapse = ''))
                                     })
  
  return_list <- list('ind_to_prime_seq_int_map' = ind_to_prime_seq_int_map,
                      'ind_to_prime_seq_nuc_map' = ind_to_prime_seq_nuc_map)
  
  return(return_list)
}

# ---- Prime editing state and mitochondrial reference sequence ----
prime_editing_system <- prime_editing_enabled(input_args)
prime_editing_backend <- NULL
ind_to_prime_seq_int_map <- NULL
ind_to_prime_seq_nuc_map <- NULL
prime_editing_efficiency_by_position <- NULL

# always going to randomly generate the mt sequence
baseline_seq_ints_mt <<- sample(seq(1,4), size = input_args$mito_genome_length, replace = TRUE)
baseline_seq_nucs_mt <<- sapply(baseline_seq_ints_mt, convert_int_to_nuc)




##### cell type-specific work ... 
# generate cell type transition matrix lists for induced and uninduced conditions:
# convert transition matrix values into a list of lists, where outer key transitions into inner key
cell_type_names <<- names(input_args$cell_type_dict$cell_type_params)

tm_lists_res <- make_cell_type_transition_lists(cell_type_names = cell_type_names)
uninduced_tm_list <- tm_lists_res[['uninduced_tm_list']]
induced_tm_list <- tm_lists_res[['induced_tm_list']]








# depending on the user's chosen nucleotide substitution model, extract the relevant parameters
# this will have to be done for both the barcode and mt mutational processes
# return substitution probability matrix
#' Build a nucleotide substitution probability matrix from model parameters
#'
#' Dispatches on the model name and reads that model's positional parameters out
#' of a semicolon-delimited string. The equilibrium base frequencies used by
#' F81, HKY, and GTR are measured from `sequence_with_targets` rather than
#' supplied.
#'
#' @param raw_cla_submodel Model name: `'JC'`, `'K80'`, `'K81'`, `'F81'`,
#'   `'HKY'`, or `'GTR'`.
#' @param selected_sub_model Semicolon-delimited model parameters, in the order
#'   the chosen model expects.
#' @param sequence_with_targets Character nucleotide vector whose composition
#'   supplies the base frequencies.
#' @return A list with `sub_model_params_list` (the named parsed parameters) and
#'   `sub_prob_mat` (the 4-by-4 A/G/C/T substitution probability matrix).
parse_sub_model_params <- function(raw_cla_submodel, selected_sub_model, sequence_with_targets){
  split_mod_params <- process_cla_string(selected_sub_model)
  nuc_fracs <- sequence_nucleotide_fractions(sequence_with_targets)

  sub_model_params_list <- list()
  
  if(raw_cla_submodel == 'JC'){
    sub_model_params_list[['model_overall_subrate']] <- split_mod_params[1]
    sub_prob_mat <- jc_sub_rate_mat(overall_sub_rate = sub_model_params_list[['model_overall_subrate']])
  }
  
  if(raw_cla_submodel == 'K80'){
    sub_model_params_list[['transition_to_transversion_ratio']] <- split_mod_params[1]
    sub_model_params_list[['transition_rate']] <- split_mod_params[2]
    sub_model_params_list[['transversion_rate']] <- split_mod_params[3]
    
    sub_prob_mat <- k80_sub_rate_mat(transition_to_transversion_ratio = sub_model_params_list[['transition_to_transversion_ratio']],
                                     transition_rate = sub_model_params_list[['transition_rate']],
                                     transversion_rate = sub_model_params_list[['transversion_rate']])
  }
  
  if(raw_cla_submodel == 'K81'){
    # If K81: 'transition_rate; transversion_rate_weakstrong_conserved; transversion_rate_aminoketo_conserved'
    # If F81: 'baseline_overall_subrate'
    # If HKY: 'transition_to_transversion_ratio; baseline_transition_rate; baseline_transversion_rate'
    # If GTR: 'AG_rate; AC_rate; AT_rate; GC_rate; GT_rate; CT_rate'")
    sub_model_params_list[['transition_rate']] <- split_mod_params[1]
    sub_model_params_list[['transverison_rate_weakstrong_conserved']] <- split_mod_params[2]
    sub_model_params_list[['transversion_rate_aminoketo_conserved']] <- split_mod_params[3]
    
    sub_prob_mat <- k81_sub_rate_mat(transition_rate = sub_model_params_list[['transition_rate']],
                                     transversion_rate_weakstrong_conserved = sub_model_params_list[['transverison_rate_weakstrong_conserved']],
                                     transversion_rate_aminoketo_conserved = sub_model_params_list[['transversion_rate_aminoketo_conserved']])
  }
  
  if(raw_cla_submodel == 'F81'){
    # calculate nucleotide fractions based on provided sequence with targets
    sub_model_params_list[['baseline_overall_subrate']] <- split_mod_params[1]
    
    frac_a <- nuc_fracs[['A']]
    frac_g <- nuc_fracs[['G']]
    frac_c <- nuc_fracs[['C']]
    frac_t <- nuc_fracs[['T']]
    
    sub_model_params_list[['frac_a']] <- frac_a
    sub_model_params_list[['frac_g']] <- frac_g
    sub_model_params_list[['frac_c']] <- frac_c
    sub_model_params_list[['frac_t']] <- frac_t
    
    sub_prob_mat <- f81_sub_rate_mat(frac_a = sub_model_params_list[['frac_a']],
                                     frac_g = sub_model_params_list[['frac_g']],
                                     frac_c = sub_model_params_list[['frac_c']],
                                     frac_t = sub_model_params_list[['frac_t']],
                                     baseline_overall_sub_rate = sub_model_params_list[['baseline_overall_subrate']])
  }
  
  if(raw_cla_submodel == 'HKY'){
    sub_model_params_list[['transition_to_transversion_ratio']] <- split_mod_params[1]
    sub_model_params_list[['baseline_transition_rate']] <- split_mod_params[2]
    sub_model_params_list[['baseline_transversion_rate']] <- split_mod_params[3]
    
    frac_a <- nuc_fracs[['A']]
    frac_g <- nuc_fracs[['G']]
    frac_c <- nuc_fracs[['C']]
    frac_t <- nuc_fracs[['T']]
    
    sub_model_params_list[['frac_a']] <- frac_a
    sub_model_params_list[['frac_g']] <- frac_g
    sub_model_params_list[['frac_c']] <- frac_c
    sub_model_params_list[['frac_t']] <- frac_t
    
    sub_prob_mat <- hky_sub_rate_mat(frac_a = sub_model_params_list[['frac_a']],
                                     frac_g = sub_model_params_list[['frac_g']],
                                     frac_c = sub_model_params_list[['frac_c']],
                                     frac_t = sub_model_params_list[['frac_t']],
                                     transition_to_transversion_ratio = sub_model_params_list[['transition_to_transversion_ratio']],
                                     baseline_transition_rate = sub_model_params_list[['baseline_transition_rate']],
                                     baseline_transversion_rate = sub_model_params_list[['baseline_transversion_rate']])
  }
  
  if(raw_cla_submodel == 'GTR'){
    sub_model_params_list[['AG_rate']] <- split_mod_params[1]
    sub_model_params_list[['AC_rate']] <- split_mod_params[2]
    sub_model_params_list[['AT_rate']] <- split_mod_params[3]
    sub_model_params_list[['GC_rate']] <- split_mod_params[4]
    sub_model_params_list[['GT_rate']] <- split_mod_params[5]
    sub_model_params_list[['CT_rate']] <- split_mod_params[6]
    
    frac_a <- nuc_fracs[['A']]
    frac_g <- nuc_fracs[['G']]
    frac_c <- nuc_fracs[['C']]
    frac_t <- nuc_fracs[['T']]
    
    sub_model_params_list[['frac_a']] <- frac_a
    sub_model_params_list[['frac_g']] <- frac_g
    sub_model_params_list[['frac_c']] <- frac_c
    sub_model_params_list[['frac_t']] <- frac_t
    
    sub_prob_mat <- gtr_sub_rate_mat(frac_a = sub_model_params_list[['frac_a']],
                                     frac_g = sub_model_params_list[['frac_g']],
                                     frac_c = sub_model_params_list[['frac_c']],
                                     frac_t = sub_model_params_list[['frac_t']],
                                     ag_rate = sub_model_params_list[['AG_rate']],
                                     ac_rate = sub_model_params_list[['AC_rate']],
                                     at_rate = sub_model_params_list[['AT_rate']],
                                     gc_rate = sub_model_params_list[['GC_rate']],
                                     gt_rate = sub_model_params_list[['GT_rate']],
                                     ct_rate = sub_model_params_list[['CT_rate']])
  }
  
  return_list <- list('sub_model_params_list' = sub_model_params_list,
                      'sub_prob_mat' = sub_prob_mat)
  
  return(return_list)
}


# accepts a substitution model and a sequence as input, 
# and returns a list of basepos:transition_prob for transitions
# and a list of basepos:[transversion_base1:transversion_prob1, transversion_base2:transversion_prob2] for transversions
#' Per-position transition probabilities for a sequence
#'
#' @param sequence_with_targets Character nucleotide vector, one element per
#'   sequence position.
#' @param sub_prob_mat 4-by-4 substitution probability matrix whose rows and
#'   columns are ordered `A, G, C, T`.
#' @return A list with one element per sequence position, holding the
#'   probability of that base's transition (A/G or C/T).
generate_transition_basepos_list <- function(sequence_with_targets, sub_prob_mat){
  # sequence with_tarets is a vector of characters of length == length(sequence) 
  
  # to generate transition list: 
  # for each element in the integer sequence, access the sub mat at the position corresponding to 
  # row == base in sequence, column == base corresponding to transition base
  # note that (as indicated by the names), transition_matches are ordered AGCT
  transition_matches <- as.list(c('G', 'A', 'T', 'C'))
  names(transition_matches) <- c('A', 'G', 'C', 'T')
  
  transition_list <- lapply(sequence_with_targets, FUN = function(base){
    base_int_from <- which(names(transition_matches) == base)
    base_int_to <- which(names(transition_matches) == transition_matches[[base]])
    return(sub_prob_mat[base_int_from, base_int_to])
  })
  
  return(transition_list)
}

#' Per-position transversion probabilities for a sequence
#'
#' @param sequence_with_targets Character nucleotide vector, one element per
#'   sequence position.
#' @param sub_prob_mat 4-by-4 substitution probability matrix whose rows and
#'   columns are ordered `A, G, C, T`.
#' @return A list with one element per sequence position; each element is a
#'   two-element list named by the two possible destination bases and holding
#'   their transversion probabilities.
generate_transversion_basepos_list <- function(sequence_with_targets, sub_prob_mat){
  # transversion_list will have list structure
  # where outer names are base position number
  # and each element of this list is another list
  # whose names are the two possible bases that the nucleotide at that position in the sequence can undergo transversion to
  # and inner list values are transversion probabilities
  transversion_matches <- list('A' = c('C', 'T'),
                               'G' = c('C', 'T'),
                               'C' = c('A', 'G'),
                               'T' = c('A', 'G'))
  transversion_list <- lapply(sequence_with_targets, FUN = function(base){
    
    base_int_from <- which(names(transversion_matches) == base)
    bases_to <- transversion_matches[[base]]
    base_ints_to <- which(names(transversion_matches) == bases_to)
    
    tv_probs <- as.list(sub_prob_mat[base_int_from, base_ints_to])
    names(tv_probs) <- bases_to
    return(tv_probs)
    
  })

  return(transversion_list)
  
}

# now create target probs lists
# first determine if there is an editing window in which targets are more of a loose positional concept
# where certain bases within a context can have identical or decaying edit rates compared to target

# helper functions for expanding targets to include bases in editing window:
#############################################################

#' Lower an edit-rate class by a number of degrees
#'
#' Recursively steps `'High'` to `'Medium'` to `'Low'`; a step below `'Low'`
#' means the position has decayed to background and is no longer a target.
#'
#' @param rate Edit-rate class: `'High'`, `'Medium'`, or `'Low'`.
#' @param num_degrees Number of decay steps to apply; `0` returns `rate`
#'   unchanged.
#' @return The lowered class label, or `FALSE` once the rate falls below
#'   `'Low'`.
drop_editrate <- function(rate, num_degrees){
  # drop an edit rate down the number of degrees
  # will return FALSE if dropped out of non-uniform range
  while(num_degrees > 0){
    num_degrees <- num_degrees - 1
    if(rate == 'High'){
      return(drop_editrate('Medium', num_degrees = num_degrees))
    }
    if(rate == 'Medium'){
      return(drop_editrate('Low', num_degrees = num_degrees))
    }
    if(rate == 'Low'){
      return(FALSE)
    }
  }
  
  return(rate)
}

#' Expand base-editor targets across their editing windows
#'
#' Only positions holding the same base as the target itself are editable by the
#' base editor, so each window is restricted to matching bases. With
#' `decaying_editing` the rate drops one class within the inner half of the
#' window and two classes beyond it, and positions that decay past `'Low'` are
#' dropped. Windows are clipped to the ends of the barcode.
#'
#' @param be_editing_window One-sided window half-width in bases; `0` disables
#'   expansion and yields empty lists.
#' @param basepos_erc_be_list Named list mapping each target position to its
#'   edit-rate class.
#' @param decaying_editing Logical; whether the rate decays with distance from
#'   the target.
#' @param baseline_seq_ints_bc Integer-encoded barcode sequence, used to find
#'   matching bases and to bound the windows.
#' @return A list with `growing_window_editrates` (the newly editable positions
#'   mapped to their classes) and `be_window_to_target_ind_list` (each
#'   `be_window_<n>` mapped to the positions it covers, target included).
get_new_be_targets <- function(be_editing_window, basepos_erc_be_list, 
                               decaying_editing, baseline_seq_ints_bc){
  
  be_window_to_target_ind_list <- list()
  growing_window_editrates <- list()
  window_num <- 0
  
  # if we have an editing window
  if(be_editing_window > 0){

    # iterate through the positions (format is position:rate)
    for(target_basepos in names(basepos_erc_be_list)){
      
      window_num <- window_num + 1
      window_name <- paste0('be_window_', window_num)
      be_window_to_target_ind_list[[window_name]] <- as.integer(target_basepos)
      
      # get the edit rate associated with this target itself
      target_editrate <- unname(unlist(basepos_erc_be_list[as.character(target_basepos)]))
      
      # get the integer representation of that base
      base_int <- baseline_seq_ints_bc[as.integer(target_basepos)]
      
      # don't let lower window == 0
      lower_window <- max(1, as.integer(target_basepos) - be_editing_window)
      
      # don't let upper window exceed length of sequence
      upper_window <- min(length(baseline_seq_ints_bc), as.integer(target_basepos) + be_editing_window)
      
      # find which other bases in the window are the same as the identified target
      # shift these indices so they are with respect to the entire sequence, rather than the window
      other_same_base_inds <- which(baseline_seq_ints_bc[lower_window:upper_window] == base_int) + (lower_window-1)
      
      # remove the target itself from the window
      other_same_base_inds <- other_same_base_inds[other_same_base_inds != target_basepos]
      
      
      # if there are no other bases in this window that are identical to target
      if(length(other_same_base_inds) == 0){
        next
      }
      
      else if(length(other_same_base_inds) > 0){
        
        # if we want the editing rate to decay as we move away from the target in the window
        if(decaying_editing){
          # create new vector that will store edit rates
          # ultimately this will be combined with other_same_base_inds 
          # into a named list and appended to the existing edit rate list
          edit_window_rates <- character()
          edit_window_pos <- integer()
          
          # iterate through each identical base in the window
          for(pos in other_same_base_inds){
            
            # find number of bases away from original base this identical base is
            bases_away <- abs(pos - as.integer(target_basepos))
            
            # find how far away from the target this is, relative to window size
            rel_dist_from_target <- bases_away/be_editing_window  
            
            # if the base is within half of the one-sided editing window, drop only one degree
            if(rel_dist_from_target < 0.5){
              adjusted_editrate <- drop_editrate(target_editrate, num_degrees = 1)
            }
            
            # if it's on the opposite half, drop two degrees
            else if(rel_dist_from_target >= 0.5){
              adjusted_editrate <- drop_editrate(target_editrate, num_degrees = 2)
            }
            
            if(adjusted_editrate != FALSE){
              # only add position and rate if the rate wasn't driven down to background
              edit_window_rates <- append(edit_window_rates, adjusted_editrate)
              edit_window_pos <- append(edit_window_pos, pos)
              be_window_to_target_ind_list[[window_name]] <- append(
                be_window_to_target_ind_list[[window_name]],
                pos
              )
              
              
            }
            
          }
          
        }
        
        # if we don't want the editing rate to decay as we move away form the target
        # then we just take the rate that was specified for the target 
        else if(!decaying_editing){
          edit_window_rates <- rep(target_editrate, length(other_same_base_inds))
          edit_window_pos <- other_same_base_inds
          be_window_to_target_ind_list[[window_name]] <- c(
            as.integer(target_basepos),
            other_same_base_inds
          )
        }
        
      }
      
      if(length(edit_window_rates) > 0){
        # only append to growing list if there were some bases that were still Low or above
        new_bases <- as.list(edit_window_rates)
        names(new_bases) <- edit_window_pos
        growing_window_editrates <- append(growing_window_editrates, new_bases)
      }
    }
  }
  
  return_list <- list('growing_window_editrates' = growing_window_editrates,
                      'be_window_to_target_ind_list' = be_window_to_target_ind_list)
  return(return_list)
  
}

#' Expand nuclease targets across their editing windows
#'
#' Every position inside the window is editable, unlike the base-editor case
#' which is restricted to matching bases. With `decaying_editing` the rate drops
#' one class within the inner half of the window and two classes beyond it, and
#' positions that decay past `'Low'` are dropped. Windows are clipped to the
#' ends of the barcode.
#'
#' @param nuc_editing_window One-sided window half-width in bases; `0` disables
#'   expansion and yields empty lists.
#' @param basepos_erc_nuc_list Named list mapping each target position to its
#'   edit-rate class.
#' @param decaying_editing Logical; whether the rate decays with distance from
#'   the target.
#' @param baseline_seq_ints_bc Integer-encoded barcode sequence, used here only
#'   for its length when bounding the windows.
#' @return A list with `growing_window_editrates` (the newly editable positions
#'   mapped to their classes) and `nuc_window_to_target_ind_list` (each
#'   `nuc_window_<n>` mapped to the positions it covers, target included).
get_new_nuc_targets <- function(nuc_editing_window, basepos_erc_nuc_list, 
                                decaying_editing, baseline_seq_ints_bc){
  
  # map og targets to new window pos 
  nuc_window_to_target_ind_list <- list()
  growing_window_editrates <- list()
  
  window_num <- 0
  
  # if we have an editing window
  if(nuc_editing_window > 0){

    # iterate through the positions (format is position:rate)
    for(target_basepos in names(basepos_erc_nuc_list)){
      
      window_num <- window_num + 1
      window_name <- paste0('nuc_window_', window_num)
      # include the target itself when keeping track of which inds are in which window
      nuc_window_to_target_ind_list[[window_name]] <- as.integer(target_basepos)
      
      # get the edit rate associated with this target itself
      target_editrate <- unname(unlist(basepos_erc_nuc_list[as.character(target_basepos)]))
      
      # don't let lower window == 0
      lower_window <- max(1, as.integer(target_basepos) - nuc_editing_window)
      
      # don't let upper window exceed length of sequence
      upper_window <- min(length(baseline_seq_ints_bc), as.integer(target_basepos) + nuc_editing_window)
      
      # create window      
      other_window_inds <- seq(lower_window, upper_window)
      
      # remove the target itself from the window
      other_window_inds <- other_window_inds[other_window_inds != target_basepos]
      
      # if we want the editing rate to decay as we move away from the target in the window
      if(decaying_editing){
        # create new vector that will store edit rates
        # ultimately this will be combined with other_window_inds 
        # into a named list and appended to the existing edit rate list
        edit_window_rates <- character()
        edit_window_pos <- integer()
        
        # iterate through each identical base in the window
        for(pos in other_window_inds){
          
          # find number of bases away from original base this identical base is
          bases_away <- abs(pos - as.integer(target_basepos))
          
          # find how far away from the target this is, relative to window size
          rel_dist_from_target <- bases_away/nuc_editing_window  
          
          # if the base is within half of the one-sided editing window, drop only one degree
          if(rel_dist_from_target < 0.5){
            adjusted_editrate <- drop_editrate(target_editrate, num_degrees = 1)
          }
          
          # if it's on the opposite half, drop two degrees
          else if(rel_dist_from_target >= 0.5){
            adjusted_editrate <- drop_editrate(target_editrate, num_degrees = 2)
          }
          
          if(adjusted_editrate != FALSE){
            # only add position and rate if the rate wasn't driven down to background
            edit_window_rates <- append(edit_window_rates, adjusted_editrate)
            edit_window_pos <- append(edit_window_pos, pos)
            nuc_window_to_target_ind_list[[window_name]] <- append(
              nuc_window_to_target_ind_list[[window_name]],
              pos
            )
          }
          
        }
        
      }
      
      # if we don't want the editing rate to decay as we move away form the target
      # then we just take the rate that was specified for the target 
      else if(!decaying_editing){
        edit_window_rates <- rep(target_editrate, length(other_window_inds))
        edit_window_pos <- other_window_inds
        nuc_window_to_target_ind_list[[window_name]] <- c(
          as.integer(target_basepos),
          other_window_inds
        )
      }
      
      if(length(edit_window_rates) > 0){
        # only append to growing list if there were some bases that were still Low or above
        new_bases <- as.list(edit_window_rates)
        names(new_bases) <- edit_window_pos
        growing_window_editrates <- append(growing_window_editrates, new_bases)
        
      }
    }
  }
  
  return_list <- list('growing_window_editrates' = growing_window_editrates,
                      'nuc_window_to_target_ind_list' = nuc_window_to_target_ind_list)
  return(return_list)
}

###################################################
# end of the helper functions for editing window

# first a quick way of estimating probabilities ...
#' Convert a per-cell-cycle event probability to a per-timepoint probability
#'
#' Assumes at most one event per site per timepoint and independence between
#' timepoints, so the per-timepoint no-event probability is the per-cycle
#' no-event probability raised to `1 / timepoints_per_cell_cycle`.
#'
#' @param prob_event_per_cell_cycle Probability of the event over one full cell
#'   cycle, in [0, 1).
#' @param timepoints_per_cell_cycle Number of simulation timepoints per cell
#'   cycle.
#' @return The per-timepoint event probability that reproduces the requested
#'   per-cycle probability.
estimate_prob_per_timept <- function(prob_event_per_cell_cycle, timepoints_per_cell_cycle){
  # the intuition here is that the user can estimate the probability of an event occurring at a 
  # timepoint whose resolution is finer than its cell cycle length resolution
  # this is applicable to estimating mutation probs per edit timepoint or death probs per edit timepoint,
  # given an overarching probability of these events per cell cycle.
  # consider the mutation logic below:
  # assume that, at most, a target can undergo 1 mutation per site per timepoint.
  # then the rate of NO MUTATIONS per target site per cell division is 1 - muts_per_site_per_division
  # assume independence between editing outcomes at each edit point within each round of division
  # the probability of observing this NO MUTATION RATE is given by (no_mut_per_edit_pt)^n,
  # where no_mut_per_edit_pt is the probability of a mutation occurring at a site at an edit timepoint
  # and n is the number of edit timepoints in each cell cycle
  # by solving for 1-no_mut_per_edit_pt, we get a heuristic estimate of mutation prob per target per edit pt
  
  prob_no_event_per_cell_cycle <- 1 - prob_event_per_cell_cycle
  n <- timepoints_per_cell_cycle
  prob_no_event_at_timepoint <- exp(log(prob_no_event_per_cell_cycle)/n)
  return(1 - prob_no_event_at_timepoint)
}



# ---- Per-cell-type parameter containers ----
# cell type-specific substitution probability matrices
cell_type_mt_sub_prob_mat <- list()
cell_type_bc_sub_prob_mat <- list()

# cell type-specific nontarget bc mutation lists
cell_type_basepos_bc_nontarget_transition_probs <- list()
cell_type_basepos_bc_nontarget_transversion_probs <- list()
cell_type_basepos_bc_nontarget_insertion_probs <- list()
cell_type_basepos_bc_nontarget_deletion_probs <- list()

# cell type-specific nontarget mt mutation lists
cell_type_basepos_mt_nontarget_transition_probs <- list()
cell_type_basepos_mt_nontarget_transversion_probs <- list()
cell_type_basepos_mt_nontarget_insertion_probs <- list()
cell_type_basepos_mt_nontarget_deletion_probs <- list()

# cell type-specific target bc mutation lists
cell_type_basepos_bc_target_transition_probs <- list()
cell_type_basepos_bc_target_transversion_probs <- list()
cell_type_basepos_bc_target_insertion_probs <- list()
cell_type_basepos_bc_target_deletion_probs <- list()

# cell type-specific death probabilities
cell_type_death_probs <- list()

# cell type-specific sampling fractions
cell_type_poss_sampling_fracs <- list()

# cell type-specific cell cycle lengths
cell_type_cell_cycle_length <- list()


# assign basepos_erc_be_list and basepos_erc_nuc_list once for all cell types

basepos_erc_be_list <- bc_generation_return_list[['be_basepos_editrate_classes']]
basepos_erc_nuc_list <- bc_generation_return_list[['nuc_basepos_editrate_classes']]

# ---- Prime editing backend ----
# Prepare prime editing only after the target positions have been constructed.
# The programmed edit belongs to the target itself, not to a surrounding
# nuclease editing window.
if(prime_editing_system){
  prime_editing_backend <- prepare_prime_editing_backend(
    input_args,
    target_positions = as.integer(names(basepos_erc_nuc_list)),
    params_dir = params_json_dir,
    seed = input_args$random_seed
  )
  ind_to_prime_seq_int_map <- prime_editing_sequence_integer_map(
    prime_editing_backend
  )
  ind_to_prime_seq_nuc_map <- prime_editing_sequence_character_map(
    prime_editing_backend
  )
  prime_editing_efficiency_by_position <- setNames(
    prime_editing_backend$targets$editing_efficiency,
    prime_editing_backend$targets$target_position
  )
  write_prime_editing_backend_manifest(
    prime_editing_backend,
    file.path(
      'output', 'run_specs', unique_run_id,
      'prime_editing_target_manifest.csv'
    )
  )
}

# ---- Editing-window target expansion ----
# if we have an editing window, add the relevant bases' positions to the editable bases
if(input_args$be_targets$editing_window$size > 0){
  new_be_targets_res <- get_new_be_targets(be_editing_window = input_args$be_targets$editing_window$size, 
                                           basepos_erc_be_list = basepos_erc_be_list,
                                           decaying_editing = input_args$be_targets$editing_window$decaying,
                                           baseline_seq_ints_bc = baseline_seq_ints_bc)
  new_be_targets <- new_be_targets_res[['growing_window_editrates']]
  be_window_to_target_ind_list <- new_be_targets_res[['be_window_to_target_ind_list']]
  
  # also flip the mapping, now from target to window name:
  be_target_to_window_ind_list <- split(
    rep(names(be_window_to_target_ind_list), vapply(be_window_to_target_ind_list, length, integer(1))),  
    unlist(be_window_to_target_ind_list)
  )
  
  basepos_erc_be_list <- append(basepos_erc_be_list, new_be_targets)
} else{
  be_window_to_target_ind_list <- NULL
  be_target_to_window_ind_list <- NULL
}


if(input_args$nuclease_targets$editing_window$size > 0 &&
   !prime_editing_system){
  new_nuc_targets_res <- get_new_nuc_targets(nuc_editing_window = input_args$nuclease_targets$editing_window$size, 
                                             basepos_erc_nuc_list = basepos_erc_nuc_list,
                                             decaying_editing = input_args$nuclease_targets$editing_window$decaying,
                                             baseline_seq_ints_bc = baseline_seq_ints_bc)
  new_nuc_targets <- new_nuc_targets_res[['growing_window_editrates']]
  nuc_window_to_target_ind_list <- new_nuc_targets_res[['nuc_window_to_target_ind_list']]
  
  # also flip the mapping, now from target to window name:
  nuc_target_to_window_ind_list <- split(
    rep(names(nuc_window_to_target_ind_list), vapply(nuc_window_to_target_ind_list, length, integer(1))),  
    unlist(nuc_window_to_target_ind_list)
  )
  
  
  basepos_erc_nuc_list <- append(basepos_erc_nuc_list, new_nuc_targets)
} else{
  nuc_window_to_target_ind_list <- NULL
  nuc_target_to_window_ind_list <- NULL
}

# ---- Per-cell-type mutation probability tables ----
for(celltype in cell_type_names){
  
  print(paste0('Assigning ', celltype, ' params'))
  
  # assign cell type-specific params that depend on editing induction
  for(induction in c('induced_editing_params', 'uninduced_editing_params')){
    
    
    
    
    temp_mt_sub_model_list <- parse_sub_model_params(raw_cla_submodel = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_substitution_model,
                                                     selected_sub_model = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_sub_model_params,
                                                     sequence_with_targets = baseline_seq_nucs_mt)
    temp_bc_sub_model_list <- parse_sub_model_params(raw_cla_submodel = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_substitution_model,
                                                     selected_sub_model = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_sub_model_params,
                                                     sequence_with_targets = baseline_seq_nucs_bc)
    
    cell_type_mt_sub_prob_mat[[celltype]][[induction]] <- temp_mt_sub_model_list[['sub_prob_mat']]
    cell_type_bc_sub_prob_mat[[celltype]][[induction]] <- temp_bc_sub_model_list[['sub_prob_mat']]
    
    basepos_bc_nontarget_transition_probs <- generate_transition_basepos_list(sequence_with_targets = baseline_seq_nucs_bc,
                                                                              sub_prob_mat = temp_bc_sub_model_list[['sub_prob_mat']])
    basepos_bc_nontarget_transversion_probs <- generate_transversion_basepos_list(sequence_with_targets = baseline_seq_nucs_bc,
                                                                                  sub_prob_mat = temp_bc_sub_model_list[['sub_prob_mat']])
    basepos_mt_nontarget_transition_probs <- generate_transition_basepos_list(sequence_with_targets = baseline_seq_nucs_mt,
                                                                              sub_prob_mat = temp_mt_sub_model_list[['sub_prob_mat']])
    basepos_mt_nontarget_transversion_probs <- generate_transversion_basepos_list(sequence_with_targets = baseline_seq_nucs_mt,
                                                                                  sub_prob_mat = temp_mt_sub_model_list[['sub_prob_mat']])
    
    # parse the input background indel rates for mt and bc WHILE CONVERTING TO PROB PER TIMEPT RATHER THAN PROB PER CELL CYCLE
    this_celltype_timepts_per_cc <- input_args$cell_type_dict$cell_type_params[[celltype]]$cell_cycle_length / time_inc
    bc_bg_insertion_prob <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_bg_insertion_prob_per_division,
                                                     timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    bc_bg_deletion_prob <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_bg_deletion_prob_per_division,
                                                    timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    mt_bg_insertion_prob <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_bg_insertion_prob_per_division,
                                                     timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    mt_bg_deletion_prob <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_bg_deletion_prob_per_division,
                                                    timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    

    basepos_bc_nontarget_insertion_probs <- as.list(rep(bc_bg_insertion_prob, input_args$bc_length))
    basepos_bc_nontarget_deletion_probs <- as.list(rep(bc_bg_deletion_prob, input_args$bc_length))
    basepos_mt_nontarget_insertion_probs <- as.list(rep(mt_bg_insertion_prob, input_args$mito_genome_length))
    basepos_mt_nontarget_deletion_probs <- as.list(rep(mt_bg_deletion_prob, input_args$mito_genome_length))
    
    
    
    
    # add gamma distribution heterogeneity to mt and bc mutation probs:
    if(!is.null(input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_nontarget_heterogeneity_gamma)){
      mt_nontarget_hetero <- input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$mt_nontarget_heterogeneity_gamma
      basepos_mt_nontarget_names_vec <- c('basepos_mt_nontarget_transition_probs',
                                          'basepos_mt_nontarget_transversion_probs',
                                          'basepos_mt_nontarget_insertion_probs',
                                          'basepos_mt_nontarget_deletion_probs')
      for (prob_list_name in basepos_mt_nontarget_names_vec){
        assign(prob_list_name, nontarget_scale_gamma_heterogeneity(position_er_list = get(prob_list_name),
                                                                   shape_param = mt_nontarget_hetero[['shape_param']],
                                                                   num_discrete_bins = mt_nontarget_hetero[['num_bins']],
                                                                   bin_agg_metric = mt_nontarget_hetero[['agg_metric']]))
        
      }
    }
    
    if(!is.null(input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_nontarget_heterogeneity_gamma)){
      bc_nontarget_hetero <- input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$bc_nontarget_heterogeneity_gamma
      basepos_bc_nontarget_names_vec <- c('basepos_bc_nontarget_transition_probs',
                                          'basepos_bc_nontarget_transversion_probs',
                                          'basepos_bc_nontarget_insertion_probs',
                                          'basepos_bc_nontarget_deletion_probs')
      for (prob_list_name in basepos_bc_nontarget_names_vec){
        assign(prob_list_name, nontarget_scale_gamma_heterogeneity(position_er_list = get(prob_list_name),
                                                                   shape_param = bc_nontarget_hetero[['shape_param']],
                                                                   num_discrete_bins = bc_nontarget_hetero[['num_bins']],
                                                                   bin_agg_metric = bc_nontarget_hetero[['agg_metric']]))
        
      }
    }
    
    
    target_insertion_prob_mean_estimate <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$nuc_insertions_per_target_per_division, 
                                                                    timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    target_deletion_prob_mean_estimate <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$nuc_deletions_per_target_per_division,
                                                                   timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    target_be_prob_mean_estimate <- estimate_prob_per_timept(prob_event_per_cell_cycle = input_args$cell_type_dict$cell_type_params[[celltype]][[induction]]$be_mutations_per_target_per_division,
                                                             timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
    
    
    
    
    if(!is.null(basepos_erc_be_list)){
      
      if(be_mutation_type == 'transition'){
        
        basepos_bc_target_transition_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                         h_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'High')]), 
                                                                                         m_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Medium')]), 
                                                                                         l_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Low')]),
                                                                                         shape_param = 0.5,
                                                                                         scale_param = target_be_prob_mean_estimate/0.5,
                                                                                         num_bootstrap_draws = 1000)
        basepos_bc_target_transversion_probs <- list()
        
        
        if(close_be_window_after_edit){
          close_transition_window_after_edit <- TRUE
        }

      } else if(be_mutation_type == 'transversion'){
        # if the BE causes transversions, the transition basepos edit rate list will be NULL
        basepos_bc_target_transversion_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                           h_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'High')]), 
                                                                                           m_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Medium')]), 
                                                                                           l_pos = as.integer(names(basepos_erc_be_list)[which(basepos_erc_be_list == 'Low')]),
                                                                                           shape_param = 0.5,
                                                                                           scale_param = target_be_prob_mean_estimate/0.5,
                                                                                           num_bootstrap_draws = 1000)
        
        
        basepos_bc_target_transition_probs <- list()
        
        if(close_be_window_after_edit){
          close_transversion_window_after_edit <- TRUE
        }
        
      }  
    } else{
      basepos_bc_target_transversion_probs <- list()
      basepos_bc_target_transition_probs <- list()
    }
    
    # insertion_HML_gamma_scale
    if(!is.null(basepos_erc_nuc_list)){
      
      basepos_bc_target_insertion_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                      h_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'High')]), 
                                                                                      m_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Medium')]), 
                                                                                      l_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Low')]),
                                                                                      shape_param = 0.5,
                                                                                      scale_param = target_insertion_prob_mean_estimate/0.5,
                                                                                      num_bootstrap_draws = 1000)
      
      basepos_bc_target_deletion_probs <- SIMPLIFY_target_site_gamma_based_sub_rates(sequence_length = input_args$bc_length, 
                                                                                     h_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'High')]), 
                                                                                     m_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Medium')]), 
                                                                                     l_pos = as.integer(names(basepos_erc_nuc_list)[which(basepos_erc_nuc_list == 'Low')]),
                                                                                     shape_param = 0.5,
                                                                                     scale_param = target_deletion_prob_mean_estimate/0.5,
                                                                                     num_bootstrap_draws = 1000)

      if(prime_editing_system){
        for(position in intersect(
          names(basepos_bc_target_insertion_probs),
          names(prime_editing_efficiency_by_position)
        )){
          basepos_bc_target_insertion_probs[[position]] <-
            prime_editing_scale_probability(
              basepos_bc_target_insertion_probs[[position]],
              prime_editing_efficiency_by_position[[position]]
            )
        }
        # Prime-editing targets resolve to their assigned template rather than
        # competing deletion alleles.
        basepos_bc_target_deletion_probs <- list()
      }
      
    } else{ # if there are no nuc targets, define empty target prob lists
      basepos_bc_target_insertion_probs <- list()
      basepos_bc_target_deletion_probs <- list()
    }
    

    # lastly, once target editing windows are finalized, 
    # force sites to be invariant as appropriate for the non-targets
    # note that this will simply be ALL of the mt inds since there are no mt targets
    mt_invariant_inds <- nontarget_get_invariant_inds(eligible_invariant_sites = seq_len(input_args$mito_genome_length),
                                                      frac_invariant = input_args$cell_type_dict$cell_type_params[[celltype]]$mt_invariant_sites)
    
    # iterate through the invariant inds and set the mutation prob at that ind to 0 for each mutation type
    for(ind in mt_invariant_inds){
      basepos_mt_nontarget_transition_probs[[as.character(ind)]] <- 0
      basepos_mt_nontarget_transversion_probs[[as.character(ind)]] <- 0
      basepos_mt_nontarget_insertion_probs[[as.character(ind)]] <- 0
      basepos_mt_nontarget_deletion_probs[[as.character(ind)]] <- 0
    }
    
    # joint target positions of targets across all four mutation types
    joint_bc_targets_vector <- as.integer(c(names(basepos_bc_target_transition_probs),
                                            names(basepos_bc_target_transversion_probs),
                                            names(basepos_bc_target_insertion_probs),
                                            names(basepos_bc_target_deletion_probs)))
    
    # eligible invariant sites are those indices of the barcode NOT in the joint target position vector                                       
    bc_eligible_invariant_sites <- setdiff(seq(1, input_args$bc_length), joint_bc_targets_vector)
    bc_invariant_inds <- nontarget_get_invariant_inds(eligible_invariant_sites = bc_eligible_invariant_sites,
                                                      frac_invariant = input_args$cell_type_dict$cell_type_params[[celltype]]$bc_invariant_sites)
    
    # now do the same iteration process for bc non-target mutation prob lists:
    # iterate through the invariant inds and set the mutation prob at that ind to 0 for each mutation type
    for(ind in bc_invariant_inds){
      basepos_bc_nontarget_transition_probs[[as.character(ind)]] <- 0
      basepos_bc_nontarget_transversion_probs[[as.character(ind)]] <- 0
      basepos_bc_nontarget_insertion_probs[[as.character(ind)]] <- 0
      basepos_bc_nontarget_deletion_probs[[as.character(ind)]] <- 0
    }
    
    
    # write all basepos nontarget params to cell-type-specific lists:
    cell_type_basepos_bc_nontarget_transition_probs[[celltype]][[induction]] <- basepos_bc_nontarget_transition_probs
    cell_type_basepos_bc_nontarget_transversion_probs[[celltype]][[induction]] <- basepos_bc_nontarget_transversion_probs
    cell_type_basepos_bc_nontarget_insertion_probs[[celltype]][[induction]] <- basepos_bc_nontarget_insertion_probs
    cell_type_basepos_bc_nontarget_deletion_probs[[celltype]][[induction]] <- basepos_bc_nontarget_deletion_probs
    cell_type_basepos_mt_nontarget_transition_probs[[celltype]][[induction]] <- basepos_mt_nontarget_transition_probs
    cell_type_basepos_mt_nontarget_transversion_probs[[celltype]][[induction]] <- basepos_mt_nontarget_transversion_probs
    cell_type_basepos_mt_nontarget_insertion_probs[[celltype]][[induction]] <- basepos_mt_nontarget_insertion_probs
    cell_type_basepos_mt_nontarget_deletion_probs[[celltype]][[induction]] <- basepos_mt_nontarget_deletion_probs
    
    # write all basepos target params to cell-type-specific lists
    cell_type_basepos_bc_target_transition_probs[[celltype]][[induction]] <- basepos_bc_target_transition_probs
    cell_type_basepos_bc_target_transversion_probs[[celltype]][[induction]] <- basepos_bc_target_transversion_probs
    cell_type_basepos_bc_target_insertion_probs[[celltype]][[induction]] <- basepos_bc_target_insertion_probs
    cell_type_basepos_bc_target_deletion_probs[[celltype]][[induction]] <- basepos_bc_target_deletion_probs
    

  #   if(!dir.exists(file.path('./celltype_prob_concordance', unique_run_id))){
  #     dir.create(file.path('./celltype_prob_concordance', unique_run_id), recursive = TRUE)
  #   }
  #   print('celltype_prob_concordance created')
  #   saveRDS(basepos_bc_target_transition_probs,
  #           file.path('./celltype_prob_concordance', unique_run_id, paste0('bc_target_transition_params_', celltype, '_', induction, '.rds')))
  #   saveRDS(basepos_bc_target_transversion_probs,
  #           file.path('./celltype_prob_concordance', unique_run_id, paste0('bc_target_transversion_params_', celltype, '_', induction, '.rds')))
  #   saveRDS(basepos_bc_target_insertion_probs,
  #           file.path('./celltype_prob_concordance', unique_run_id, paste0('bc_target_insertion_params_', celltype, '_', induction, '.rds')))
  #   saveRDS(basepos_bc_target_deletion_probs,
  #           file.path('./celltype_prob_concordance', unique_run_id, paste0('bc_target_deletion_params_', celltype, '_', induction, '.rds')))
  #   saveRDS(basepos_erc_be_list,
  #           file.path('./celltype_prob_concordance', unique_run_id, paste0('basepos_erc_be_list', celltype, '_', induction, '.rds')))
  #   
  }
  
  # now for cell type-specific params that don't depend on editing induction
  
  # write all sampling frac params to cell type specific lists
  poss_sampling_fracs <- as.numeric(input_args$cell_type_dict$cell_type_params[[celltype]]$sampling_fractions)
  cell_type_poss_sampling_fracs[[celltype]][[induction]] <- poss_sampling_fracs
  
  # cell type-specific death probabilities
  # convert prob per cell cycle to prob per timept
  cell_death_prob_per_timept <- estimate_prob_per_timept(prob_event_per_cell_cycle = as.numeric(input_args$cell_type_dict$cell_type_params[[celltype]]$death_per_cell_cycle_prob), 
                                                                  timepoints_per_cell_cycle = this_celltype_timepts_per_cc)
  cell_type_death_probs[[celltype]] <- cell_death_prob_per_timept

  # cell type-specific cell cycle lengths
  cell_type_cell_cycle_length[[celltype]] <- as.numeric(input_args$cell_type_dict$cell_type_params[[celltype]]$cell_cycle_length)

  
}

# ---- Barcode integration UMIs and reference FASTA ----
poss_num_bc_integrations <- as.integer(input_args$max_bc_ints_per_cell)

if(input_args$include_bc_umis){
  
  # generate 15 bp umis that will be prepended to barcode sequences. will be scaffold for alignment
  bc_int_umis <- as.character(lapply(seq(1:max(poss_num_bc_integrations)), function(int_num){
    paste(sample(c('A', 'G', 'C', 'T'), size = 15, replace = TRUE), collapse = '')
  }))
  
  collapsed_one_integration <- paste(baseline_seq_nucs_bc, collapse = '')
  
  # modify baseline seq to include the integration UMIs:
  bc_reference_with_int_umis <- paste(bc_int_umis, collapsed_one_integration, sep = '', collapse = '')
  
  if(!dir.exists(file.path('output', 'processed_fastas', unique_run_id, 'reference_seqs'))){
    dir.create(file.path('output', 'processed_fastas', unique_run_id, 'reference_seqs'), recursive = TRUE)
  }
  write.fasta(bc_reference_with_int_umis, 
              names = c('REFERENCE'), 
              file.out = file.path('output', 'processed_fastas', unique_run_id, 'reference_seqs', 
                                    paste0('bc_reference_including_int_umis.fasta')))
}




# ---- Output modalities and mitochondrial population settings ----
# reconstruction modalities: 
# poss_recon_modals <- process_cla_string(input_args$recon_modality, outputted_type = 'character')
poss_recon_modals <- as.character(input_args$recon_modality)


average_genomes_per_mito <- as.integer(floor(input_args$average_genomes_per_mito)) # must be a fixed integer, number of genomes per mito
starting_mito_per_cell <- as.integer(input_args$starting_mito_per_cell) # must be a fixed integer, number of mito in founding cell
max_mito_per_cell <- as.integer(input_args$max_mito_per_cell) # retained saturation setting; mito_dynamics does not currently enforce it
mito_inheritance_pattern <- tolower(as.character(input_args$mito_inheritance_pattern)) # only "random" is currently implemented
baseline_heteroplasmy_sites_frac <- as.numeric(input_args$baseline_heteroplasmy_sites_frac)
baseline_heteroplasmy_variant_frac_dist <- as.numeric(input_args$baseline_heteroplasmy_variant_frac_dist) # beta distribution params for drawing initial heteroplasmy fractions
post_mitotic_mt_deletion_frac <- as.numeric(input_args$post_mitotic_mt_deletion_frac) # what fraction of mitochondria are randomly lost after each cell division
heteroplasmy_variant_transition_prob <- as.numeric(input_args$heteroplasmy_variant_transition_prob) # given a position has a heteroplasmy variant, what fraction of the variants are of the transition base?
fusion_events_per_mito_per_division <- as.numeric(input_args$fusion_events_per_mito_per_division) # probability that a given mitochondrion undergoes a fusion event during a cell cycle
split_events_per_mito_per_division <- as.numeric(input_args$split_events_per_mito_per_division) # probability that a given mitochondrion undergoes a splitting event during a cell cycle
heteroplasmy_rel_with_fitness <- tolower(as.character(input_args$heteroplasmy_rel_with_fitness)) # is increased heteroplasmy positively, negatively, or neutrally associated with cellular fitness?
init_heteroplasmy_survive_prob <- as.numeric(input_args$init_heteroplasmy_survive_prob) # probability of cell survival for the founder cell, given its initialized heteroplasmy
consider_cell_heteroplasmy_scores <- as.logical(input_args$consider_cell_heteroplasmy_scores) # boolean for whether heteroplasmy scores should be calculated and involved in processes
hetero_sd <- as.numeric(input_args$heteroplasmy_standard_deviation)


# ---- Founder mitochondrion-to-genome map ----
# map each mitochondrion to its corresponding genome inds by first generating the number of genomes per mitochondrion (poisson dist with EV average genomes per mito)
sizes_of_mito <- rpois(n = starting_mito_per_cell, lambda = average_genomes_per_mito)

# all mito have at least one genome:
sizes_of_mito <- pmax(sizes_of_mito, 1)
init_num_mito_genomes <- sum(sizes_of_mito)

# use cumulative sum to map to starting inds
mito_size_cumsum <- cumsum(sizes_of_mito)
# prepend 0 for simplicity:
mito_size_cumsum <- c(0, mito_size_cumsum)
mito_to_genome_map = list() # will serve as a cell-specific list that maps mito number to genome row inds 
for(cumsum_ind in seq(1, length(mito_size_cumsum)-1)){
  mito_to_genome_map[[cumsum_ind]] <- seq(mito_size_cumsum[cumsum_ind] + 1, mito_size_cumsum[cumsum_ind + 1])
}


# ---- Founder heteroplasmy initialization ----
# create initial mt genome matrix according to heteroplasmy params and initial wildtype sequence:

# generate inds at which heteroplasmy will be present:
heteroplasmy_inds <- sample(seq(input_args$mito_genome_length),
                            size = baseline_heteroplasmy_sites_frac*input_args$mito_genome_length,
                            replace = FALSE)
# map each heteroplasmy ind to a fraction present across all mt genomes in a cell:
# for each position of heteroplasmy, allow for some combination of point mutations
# generate per-site heteroplasmy penetrance fraction distribution from provided beta params
penetrance_dist_vals <- rbeta(n = 10000, shape1 = baseline_heteroplasmy_variant_frac_dist[1],
                              shape2 = baseline_heteroplasmy_variant_frac_dist[2])


init_i_vals <- c()
init_j_vals <- c()
init_x_vals <- c()

# each variant will be assigned a severity score that will contribute to cellular survival
# assumes independent effects of mt heteroplasmy variants on cell survival
# init weights can associated heteroplasmy with positive, negative, or balanced fitness
# severity fractions will be generated from a bimodal mixed Gaussian distribution
# without loss of generality, variants associated with increased cellular fitness are assigned positive scores
# and variants associated with decreased cellular fitness are assigned negative scores
# user-provided weight of the fraction of deleterious heteroplasmy variants controls the severity score assignment process

# random draws based on the positive score weight: if runif draw value is above positive score weight, draw from positive dist, else draw from negative
#' Draw variant severity scores from a two-component Gaussian mixture
#'
#' Each draw comes from the `mean2` component with probability `mean2_weight`
#' and from the `mean1` component otherwise. Callers pass `mean1 = -1` for the
#' deleterious component and `mean2 = 1` for the beneficial one.
#'
#' @param num_draws Number of scores to draw.
#' @param mean1 Mean of the first mixture component.
#' @param mean2 Mean of the second mixture component.
#' @param sigma Standard deviation shared by both components.
#' @param mean2_weight Probability in [0, 1] that a draw uses the `mean2`
#'   component.
#' @return A numeric vector of length `num_draws`.
draw_severity_scores <- function(num_draws, mean1, mean2, sigma, mean2_weight){
  
  # determine if each draw will come from mean1 dist or mean2 dist:
  mean2_dist <- runif(num_draws) <= mean2_weight
  
  output_scores <- numeric(length = num_draws)
  
  # assign output scores to mean2 dists based on which were uniformly drawn to be coming from dist 2
  output_scores[mean2_dist] <- rnorm(n = sum(mean2_dist), mean = mean2, sd = sigma)
  
  # and vice versa
  output_scores[!mean2_dist] <- rnorm(n = sum(!mean2_dist), mean = mean1, sd = sigma)
  
  return(output_scores)
  
}

positive_score_weight <- 1

if(consider_cell_heteroplasmy_scores){
  positive_score_weight = 1 - as.numeric(input_args$fraction_deleterious_heteroplasmy_variants)
 
}

# structure of heteroplasmy_severity_score_list will be key = 'mt_{position_num}_{mutated_base_int}', value = severity score
heteroplasmy_severity_score_list <- list()

# store the variant fractions associated with each heteroplasmy variant 
heteroplasmy_variant_fractions <- list()

heteroplasmy_allelic_fractions <- c()
heteroplasmy_counts <- c()


#' Seed one heteroplasmic site across the founder cell's mitochondrial genomes
#'
#' Draws a penetrance fraction for the site, binomially picks which genomes
#' carry a variant, then splits the carriers between the transition base and the
#' two transversion bases (the two transversions split evenly). Returns sparse
#' `(i, j, x)` triplets rather than a matrix so many sites can be bound at once.
#'
#' @param ind_num Index into `heteroplasmy_inds` naming the site to seed.
#' @param draw_severity_scores Function used to draw severity scores, passed in
#'   rather than looked up globally.
#' @param heteroplasmy_inds Integer positions of every heteroplasmic site.
#' @param penetrance_dist_vals Pool of beta-distributed penetrance fractions to
#'   sample one value from.
#' @param init_num_mito_genomes Number of mitochondrial genome rows in the
#'   founder profile.
#' @param transition_matches Integer vector mapping each base code to its
#'   transition partner.
#' @param transversion_matches List mapping each base code to its two
#'   transversion partners.
#' @param heteroplasmy_variant_transition_prob Probability that a carrier genome
#'   takes the transition base rather than a transversion.
#' @param consider_cell_heteroplasmy_scores Logical; when `FALSE` no severity
#'   scores are drawn and `severity_scores` comes back empty.
#' @param positive_score_weight Mixture weight on the beneficial severity
#'   component.
#' @param hetero_sd Standard deviation of the severity score mixture.
#' @param baseline_seq_ints_mt Integer-encoded reference mitochondrial sequence,
#'   used to look up the wild-type base at the site.
#' @return A list with `rows` (an `(i, j, x)` tibble of genome row, position,
#'   and mutated base code), plus `variant_frac`, `severity_scores`, and
#'   `counts`, each named `transition`, `transversion_1`, `transversion_2`. All
#'   four come back empty when no genome carries the variant.
#' @note Not called by the current script; the founder-initialization loop below
#'   reimplements the same draws inline with `data.table` rows.
heteroplasmy_one_site <- function(ind_num,
                                  draw_severity_scores,
                                  heteroplasmy_inds = heteroplasmy_inds,
                                  penetrance_dist_vals = penetrance_dist_vals, 
                                  init_num_mito_genomes = init_num_mito_genomes,
                                  transition_matches = transition_matches,
                                  transversion_matches = transversion_matches,
                                  heteroplasmy_variant_transition_prob = heteroplasmy_variant_transition_prob,
                                  consider_cell_heteroplasmy_scores = consider_cell_heteroplasmy_scores,
                                  positive_score_weight = positive_score_weight,
                                  hetero_sd = hetero_sd,
                                  baseline_seq_ints_mt = baseline_seq_ints_mt
                                  ){
  
  
  penetrant_frac <- sample(penetrance_dist_vals, size = 1)
  num_genomes_with_var <- rbinom(n = 1, size = init_num_mito_genomes, prob = penetrant_frac)
  
  if(num_genomes_with_var == 0){
    return(list(rows = tibble(i = integer(),
                              j = integer(),
                              x = integer()),
                variant_frac = numeric(),
                severity_scores = numeric(),
                counts = integer()
    ))
  }
  
  mt_genomes_with_variant <- sample(seq(init_num_mito_genomes), size = num_genomes_with_var, replace = FALSE)
  wt_base <- baseline_seq_ints_mt[heteroplasmy_inds[ind_num]]
  transition_base <- transition_matches[wt_base]
  transversion_bases <- transversion_matches[[wt_base]]
  
  num_mt_genomes_with_transition <- rbinom(n = 1, size = length(mt_genomes_with_variant), 
                                           prob = heteroplasmy_variant_transition_prob)
  mt_genomes_with_transition <- sample(mt_genomes_with_variant, size = num_mt_genomes_with_transition, replace = FALSE)
  
  
  mt_genomes_with_transversion <- setdiff(mt_genomes_with_variant, mt_genomes_with_transition)
  
  # prob = 0.5 because now only dividing between two transversion options
  num_mt_genomes_with_transversion1 <- rbinom(n = 1, size = length(mt_genomes_with_transversion), 
                                              prob = 0.5)
  mt_genomes_with_transversion1 <- sample(mt_genomes_with_transversion, size = num_mt_genomes_with_transversion1, 
                                          replace = FALSE)
  
  mt_genomes_with_transversion2 <- setdiff(mt_genomes_with_transversion, mt_genomes_with_transversion1)
  num_mt_genomes_with_transversion2 <- length(mt_genomes_with_transversion2)
  
  # generate variant fractions for each of transition and transversion variants:
  transition_var_frac <- length(mt_genomes_with_transition)/length(mt_genomes_with_variant)
  transversion1_var_frac <- length(mt_genomes_with_transversion1)/length(mt_genomes_with_variant)
  transversion2_var_frac <- length(mt_genomes_with_transversion2)/length(mt_genomes_with_variant)
  
  
  
  rows_transition <- tibble(i = mt_genomes_with_transition,
                            j = heteroplasmy_inds[ind_num],
                            x = transition_base)
  
  rows_tv1 <- tibble(i = mt_genomes_with_transversion1,
                     j = heteroplasmy_inds[ind_num],
                     x = transversion_bases[1])
  
  rows_tv2 <- tibble(i = mt_genomes_with_transversion2,
                     j = heteroplasmy_inds[ind_num],
                     x = transversion_bases[2])
  
  rows_all <- bind_rows(rows_transition, rows_tv1, rows_tv2)
  
  variant_frac <- c('transition' = length(mt_genomes_with_transition)/length(mt_genomes_with_variant),
                    'transversion_1' = length(mt_genomes_with_transversion1)/length(mt_genomes_with_variant),
                    'transversion_2' = length(mt_genomes_with_transversion2)/length(mt_genomes_with_variant))
  
  severity_scores <- numeric()
  if(consider_cell_heteroplasmy_scores){
    severity_scores <- c('transition' = transition_severity <- draw_severity_scores(num_draws = 1, 
                                                                                    mean2_weight = positive_score_weight, 
                                                                                    mean1 = -1, 
                                                                                    mean2 = 1, 
                                                                                    sigma = hetero_sd),
                         'transversion_1' = transition_severity <- draw_severity_scores(num_draws = 1, 
                                                                                    mean2_weight = positive_score_weight, 
                                                                                    mean1 = -1, 
                                                                                    mean2 = 1, 
                                                                                    sigma = hetero_sd),
                         'transversion_2' = transition_severity <- draw_severity_scores(num_draws = 1, 
                                                                                    mean2_weight = positive_score_weight, 
                                                                                    mean1 = -1, 
                                                                                    mean2 = 1, 
                                                                                    sigma = hetero_sd))
  }
  
  counts <- c('transition' = length(mt_genomes_with_transition),
              'transversion_1' = length(mt_genomes_with_transversion1),
              'transversion_2' = length(mt_genomes_with_transversion2))
  return(list(rows = rows_all,
              variant_frac = variant_frac,
              severity_scores = severity_scores,
              counts = counts
  ))
  
}



# ---- Seed heteroplasmy variants across the founder genomes ----
################################################### new way
new_hetero_indexing_time_start <- Sys.time()
print('starting heteroplasmy indexing')

rows_list <- vector('list', length(heteroplasmy_inds))


variant_frac_vec <- numeric() # will store fraction transition vs tv1 vs tv2 for each variant
variant_count_vec <- integer() # will store counts of genomes for each variant
severity_score_vec <- numeric() # will map variants to severity scores

if(consider_cell_heteroplasmy_scores){
  heteroplasmy_allelic_fractions <- numeric()
  heteroplasmy_counts <- integer()
}

for(ind_num in seq_along(heteroplasmy_inds)){
  
  if(ind_num %% 50 == 0){
    print(paste0(ind_num, '/', length(heteroplasmy_inds)))
  }
  
  # what fraction of mito genomes at this ind should have a heteroplasmy variant?
  penetrant_frac <- sample(penetrance_dist_vals, size = 1)
  
  num_genomes_with_var <- rbinom(n = 1, size = init_num_mito_genomes, prob = penetrant_frac)
  
  # if no genomes have the variant, continue to the next:
  if(num_genomes_with_var == 0){
    next
  }
  
  # based on this fraction, which mito genomes have ANY variant? sample from the initial number of mito genomes in the cell, not the max
  # later will distribute these genomes into specific variants
  # mt_genomes_with_variant <- sample(seq(init_num_mito_genomes), prob = penetrant_frac)
  mt_genomes_with_variant <- sample(seq_len(init_num_mito_genomes), size = num_genomes_with_var, replace = FALSE)
  
  # assign each variant a point mutation based on the original nucleotide at that position
  wt_base <- baseline_seq_ints_mt[heteroplasmy_inds[ind_num]]
  
  # for now, split variants between transition and transversion according to user-provided heteroplasmy_variant_transition_prob
  # transversions will be equiprobable among the remaining variant but not transition bases
  transition_base <- transition_matches[wt_base]
  transversion_bases <- transversion_matches[[wt_base]]
  
  num_mt_genomes_with_transition <- rbinom(n = 1, size = length(mt_genomes_with_variant), 
                                           prob = heteroplasmy_variant_transition_prob)
  mt_genomes_with_transition <- sample(mt_genomes_with_variant, size = num_mt_genomes_with_transition, replace = FALSE)
  
  
  mt_genomes_with_transversion <- setdiff(mt_genomes_with_variant, mt_genomes_with_transition)
  
  # prob = 0.5 because now only dividing between two transversion options
  num_mt_genomes_with_transversion1 <- rbinom(n = 1, size = length(mt_genomes_with_transversion), 
                                              prob = 0.5)
  mt_genomes_with_transversion1 <- sample(mt_genomes_with_transversion, size = num_mt_genomes_with_transversion1, 
                                          replace = FALSE)
  
  mt_genomes_with_transversion2 <- setdiff(mt_genomes_with_transversion, mt_genomes_with_transversion1)
  num_mt_genomes_with_transversion2 <- length(mt_genomes_with_transversion2)
  
  #' Internal: build one mutation-class triplet table for the current site
  #'
  #' @param genomes Integer mitochondrial genome row indices carrying this
  #'   mutation class.
  #' @param base Integer code of the mutated base.
  #' @return A `data.table` of `(i, j, x)` rows at the current heteroplasmy
  #'   position, or `NULL` when `genomes` is empty.
  make_dt <- function(genomes, base){
    if(length(genomes) == 0){
      return(NULL)
    }
    return(data.table(i = genomes,
                      j = heteroplasmy_inds[ind_num],
                      x = base))
  }
  
  dt_transition <- make_dt(genomes = mt_genomes_with_transition, base = transition_base)
  dt_transversion1 <- make_dt(genomes = mt_genomes_with_transversion1, base = transversion_bases[1])
  dt_transversion2 <- make_dt(genomes = mt_genomes_with_transversion2, base = transversion_bases[2])

  
  
  # only bind non-empty tables:
  tables_to_bind <- list(dt_transition, dt_transversion1, dt_transversion2)
  # only retain the tables with at least one row (relevant if some types of muts don't have any genomes)
  tables_to_bind <- Filter(Negate(is.null), tables_to_bind)
  tables_to_bind <- tables_to_bind[vapply(tables_to_bind, nrow, integer(1)) > 0]
  rows_list[[ind_num]] <- rbindlist(tables_to_bind)
  
  frac_transition <- length(mt_genomes_with_transition)/length(mt_genomes_with_variant)
  frac_transversion1 <- length(mt_genomes_with_transversion1)/length(mt_genomes_with_variant)
  frac_transversion2 <- length(mt_genomes_with_transversion2)/length(mt_genomes_with_variant)
  
  genomic_position <- heteroplasmy_inds[ind_num]
  transition_mut_name <- paste0('mt_', genomic_position, '_', transition_base)
  transversion1_mut_name <- paste0('mt_', genomic_position, '_', transversion_bases[1])
  transversion2_mut_name <- paste0('mt_', genomic_position, '_', transversion_bases[2])
  
  variant_count_vec[transition_mut_name] <- length(mt_genomes_with_transition)
  variant_count_vec[transversion1_mut_name] <- length(mt_genomes_with_transversion1)
  variant_count_vec[transversion2_mut_name] <- length(mt_genomes_with_transversion2)
  
  variant_frac_vec[transition_mut_name] <- frac_transition
  variant_frac_vec[transversion1_mut_name] <- frac_transversion1
  variant_frac_vec[transversion2_mut_name] <- frac_transversion2
  
  if(consider_cell_heteroplasmy_scores){
    
    transition_severity <- draw_severity_scores(num_draws = 1, 
                                                mean2_weight = positive_score_weight, 
                                                mean1 = -1, 
                                                mean2 = 1, 
                                                sigma = hetero_sd)
    transversion1_severity <- draw_severity_scores(num_draws = 1, 
                                                   mean2_weight = positive_score_weight, 
                                                   mean1 = -1, 
                                                   mean2 = 1, 
                                                   sigma = hetero_sd)
    transversion2_severity <- draw_severity_scores(num_draws = 1, 
                                                   mean2_weight = positive_score_weight, 
                                                   mean1 = -1, 
                                                   mean2 = 1, 
                                                   sigma = hetero_sd)
    
    severity_score_vec[transition_mut_name] <- transition_severity
    severity_score_vec[transversion1_mut_name] <- transversion1_severity
    severity_score_vec[transversion2_mut_name] <- transversion2_severity
    
  }
}

full_init_hetero_dt <- rbindlist(rows_list)

init_i_vals <- full_init_hetero_dt$i
init_j_vals <- full_init_hetero_dt$j
init_x_vals <- full_init_hetero_dt$x

# if no heteroplasmy, have to manually define i, j, and x to avoid error
if(length(init_i_vals) == 0){
  init_i_vals <- c(1)
  init_j_vals <- c(1)
  init_x_vals <- c(0L)
}




# Initialize the founder profile with exactly the genome rows represented by
# mito_to_genome_map.
init_incoming_mt_profile <- sparseMatrix(i = init_i_vals, j = init_j_vals, x = init_x_vals,
                                         dims = c(init_num_mito_genomes, input_args$mito_genome_length))


print('finished heteroplasmy indexing')

new_hetero_indexing_time_end <- Sys.time()


initial_heteroplasmy_score <- 0

#' Survival probability for a cell given its heteroplasmy score
#'
#' A logistic curve anchored so that a cell whose score equals `init_score`
#' survives with probability `init_heteroplasmy_survive_prob`. Positive `beta`
#' makes higher scores more survivable; negative `beta` less.
#'
#' @param this_score The cell's current heteroplasmy score.
#' @param init_score Score at which the curve is anchored, normally the founder
#'   cell's score.
#' @param init_heteroplasmy_survive_prob Survival probability at `init_score`;
#'   must be strictly between 0 and 1.
#' @param beta Slope of the logit in score units; defaults to `1`.
#' @return The survival probability, a numeric scalar in (0, 1).
logistic_prob_survive_given_score <- function(this_score, 
                                              init_score,
                                              init_heteroplasmy_survive_prob, 
                                              beta = 1){
  # when beta > 0, higher score --> higher survival prob
  # when beta < 0, higher score --> lower survival prob
  
  # let p = probability of survival
  # anchor at: logit(init_heteroplasmy_survive_prob) = alpha + beta(this_score - init_score)
  # when this_score == init_score: logit(init_heteroplasmy_survive_prob) = alpha
  # logit(init_heteroplasmy_survive_prob) = log((init_heteroplasmy_survive_prob)/(1-init_heteroplasmy_survive_prob))
  # p = 1 / (1+exp(-(alpha + beta * (this_score - init_score))))
  
  alpha = log(init_heteroplasmy_survive_prob/(1-init_heteroplasmy_survive_prob))
  prob_survive <- 1 / (1+exp(-(alpha + beta*(this_score - init_score))))
  return(prob_survive)
}

if(consider_cell_heteroplasmy_scores){
 
  initial_heteroplasmy_score <- variant_frac_vec %*% severity_score_vec
  
}



# ---- Recovery, output, and endpoint settings ----
poss_mt_genome_recovery_probs <- as.numeric(input_args$mt_genome_recovery_prob)
poss_bc_integration_recovery_probs <- as.numeric(input_args$bc_integration_recovery_prob)
poss_fasta_types <- as.character(input_args$fasta_type)
include_var_pos_fasta <- input_args$include_var_pos_fasta
founder_cell_type <- as.character(input_args$cell_type_dict$founder_cell_type)






# # if sim lengths are specified using start:stop:inc, define sim lengths accordingly
# if(grepl(pattern = ':', x = input_args$sim_length)){
#   splits <- as.numeric(str_split(string = input_args$sim_length, pattern = ':')[[1]])
#   sim_length_stopping_points <- seq(splits[1], splits[2], by = splits[3])
# } else{ # else if specified using semicolons or a single time point
#   sim_length_stopping_points <- process_cla_string(input_args$sim_length, outputted_type = 'numeric')
#   # sim_length_stopping_points <- as.numeric(input_args$sim_length)
# }

sim_length_stopping_points <- as.numeric(input_args$sim_length)



scoremat_collapse_deletions <- as.logical(input_args$scoremat_collapse_deletions)
combine_mt_bc <- as.logical(input_args$combine_mt_bc)

binarize_mutation_scores <- as.logical(input_args$binarize_mutation_scores)
mt_allelic_fraction_thresholds <- as.numeric(input_args$mt_allelic_fraction_thresholds)

# rewrite savename if it was passed in as NULL
if(is.null(input_args$savename)){
  custom_savename <- paste('res_', input_args$num_init_cells, '_cells_', 
                           as.character(input_args$sim_length), '_maxsimlength', sep = '')
} else{
  custom_savename <- input_args$savename
}



# initialize empty vectors to avoid having to delay page appearance below
poss_trim_depths <- c()
poss_lin_strings <- c()

# old func, no longer necessary ... 
#' Lower bound on the number of cells existing at a timepoint (legacy)
#'
#' @param timept Simulation timepoint.
#' @param cc_length Cell-cycle length in the same units.
#' @return `0` when `timept` equals one cell cycle, otherwise the summed
#'   doubling series `2^t * init_pop_size` over the whole cycles before
#'   `timept`.
#' @note Retained for compatibility and still exported to the workers, but never
#'   called; it also reads `init_pop_size` from the enclosing scope.
old_cells_at_timept <- function(timept, cc_length){
  
  if(timept == cc_length){
    return(0)
  }
  lb <- sum(sapply(seq(0, timept-2*cc_length, cc_length), function(t){
    return(2^(t)*init_pop_size)
  }))
  return(lb)
}


# helper function for renumbering mito genomes after division
#' Renumber mitochondrial genome indices consecutively
#'
#' Flattens the mitochondrion-to-genome index list, replaces the indices with
#' `1..n` in flattened order, and re-splits them into the original per-
#' mitochondrion groups. Used after a daughter's profile is subset out of its
#' parent's, so the stored indices address rows of the new matrix.
#'
#' @param mito_to_genome_list List whose elements are the genome row indices
#'   belonging to each mitochondrion.
#' @return A list of the same element lengths, with genome indices renumbered
#'   from `1` and named by group position.
reassign_genome_inds <- function(mito_to_genome_list){
  flat <- unlist(mito_to_genome_list, use.names = FALSE)
  numbered <- seq_along(flat)
  return(split(numbered, rep(seq_along(mito_to_genome_list), lengths(mito_to_genome_list))))
  
}

#' Apply one round of mitochondrial fusion, fission, replication, and division
#'
#' Fusion and fission counts are drawn from Poisson distributions whose means
#' scale with the current mitochondrion count, then executed in a random order:
#' a fusion merges two mitochondria's genome sets, a fission splits one
#' mitochondrion's genomes into two non-empty halves. Every mitochondrion and
#' genome is then duplicated, mitochondria (not individual genomes) are
#' allocated binomially between the two daughters, and each daughter's genome
#' indices are renumbered against its own profile. Post-mitotic dropout finally
#' removes a binomial share of each daughter's mitochondria along with their
#' genome rows.
#'
#' @param mito_to_genome_map List mapping each mitochondrion to its genome row
#'   indices in `incoming_mito_mat`.
#' @param incoming_mito_mat Sparse mutation matrix whose rows are mitochondrial
#'   genomes and whose columns are sequence positions.
#' @param fusion_events_per_mito_per_division Poisson rate of fusion events per
#'   mitochondrion per division.
#' @param split_events_per_mito_per_division Poisson rate of fission events per
#'   mitochondrion per division.
#' @param inheritance_pattern Allocation rule; only `'random'` is implemented
#'   and any other value raises an error.
#' @param post_mitotic_mt_deletion_frac Per-mitochondrion probability of being
#'   dropped from a daughter after division; `0` disables dropout.
#' @return A list with `daughter1_mt_profile` and `daughter2_mt_profile` (sparse
#'   matrices) plus `daughter1_mito_to_genome_map` and
#'   `daughter2_mito_to_genome_map` (renumbered index lists).
#' @note Fusion and fission use `<<-` from inside `sapply`, which reaches the
#'   function's own frame, so the caller's `mito_to_genome_map` is left alone.
mito_dynamics <- function(mito_to_genome_map, incoming_mito_mat, fusion_events_per_mito_per_division,
                          split_events_per_mito_per_division,
                          inheritance_pattern = 'random', post_mitotic_mt_deletion_frac = 0){
  
  # first allow for fusion and splitting processes
  # then allocate mito according to pois to daughter cells
  
  ################ this function accomplishes the two above steps ... 
  # then allow for the random deletion of entire mito
  # then, outside this function, allow for cells to die according to heteroplasmy fractions
  
  # return list that will store the output of the func
  return_list <- list()

  if(inheritance_pattern != 'random'){
    stop("Only inheritance_pattern = 'random' is currently implemented.")
  }
  
  all_events <- c()
  
  # how many fusion events occur? get expected count by taking product of number of mito in this cell and fusion prob per mito
  num_fusion_events <- rpois(n = 1, lambda = length(mito_to_genome_map)*fusion_events_per_mito_per_division)
  all_events <- c(all_events, rep('fusion', num_fusion_events))
  
  # how many splitting events occur? 
  num_split_events <- rpois(n = 1, lambda = length(mito_to_genome_map)*split_events_per_mito_per_division)
  all_events <- c(all_events, rep('split', num_split_events))
  
  if(length(all_events) > 0){
    # generate random order of fusion and split events
    all_events <- sample(all_events)

    sapply(all_events, function(event){
      if(event == 'fusion'){
        if(length(mito_to_genome_map) > 1){
          # sample two mitochondrial numbers to fuse:
          fuse_mito_nums <- sample(seq(1, length(mito_to_genome_map)), size = 2, replace = FALSE)
          # will assume total fusions
          fused_genome_inds <- c(mito_to_genome_map[[fuse_mito_nums[1]]], mito_to_genome_map[[fuse_mito_nums[2]]])
          # update gained mito and remove lost mito from the map
          mito_to_genome_map[[fuse_mito_nums[1]]] <<- fused_genome_inds
          mito_to_genome_map[[fuse_mito_nums[2]]] <<- NULL
        }
        
      }
      else if(event == 'split'){
        # sample a mito number to split
        split_mito_num <- sample(seq(1, length(mito_to_genome_map)), size = 1)
        current_genomes <- mito_to_genome_map[[split_mito_num]]
        # split according to pois dist whose expected val is half the number of genomes in this mito
        # for now, assuming that all splits are binary!
        if(length(current_genomes) != 1){ # no splitting process occurs if there's only one genome in this mito to begin with 
          genomes_in_split1 <- rbinom(n = 1, size = length(current_genomes), prob = 0.5)
          while(genomes_in_split1 %in% c(0, length(current_genomes))){
            genomes_in_split1 <- rbinom(n = 1, size = length(current_genomes), prob = 0.5)
          }
          genomes_in_split2 <- length(current_genomes) - genomes_in_split1
          split1_genomes <- sample(current_genomes, size = genomes_in_split1, replace = FALSE)
          split2_genomes <- setdiff(current_genomes, split1_genomes)
          mito_to_genome_map[[split_mito_num]] <<- split1_genomes
          
          # name the next mito one more than the current max mito number
          mito_to_genome_map[[length(mito_to_genome_map)+1]] <<- split2_genomes
        }
      }
      
    })
  }
  
  # replicate all mitochondria to temporarily double the size of the mt mutation matrix
  rep_mito_mat <- rbind(incoming_mito_mat, incoming_mito_mat)
  
  # # also update the mito to genome map after creating new mito
  # # starting numbering at one higher than the previous number of mito
  old_num_mito <- length(mito_to_genome_map)
  copied_mito_nums <- seq(old_num_mito + 1, old_num_mito + length(mito_to_genome_map))

  # new mito genome inds are found by taking old ones and adding the numrows of the existing mito mat
  
  sapply(seq(length(copied_mito_nums)), function(rel_num){
    mito_to_genome_map[[copied_mito_nums[rel_num]]] <<- mito_to_genome_map[[rel_num]] + dim(incoming_mito_mat)[1]
  })

  # that concludes the parent mito to genome work. now need to format daughter lists and matrices:
  
  # for now only random dispersion of initial mito counts between daughter cells
  # make copies of all mitochondria, then assign MITO THEMSELVES, NOT MTGENOMES, to daughters
  if(inheritance_pattern == 'random'){ 
    
    # number of mtio passed to the first daughter will follow random poisson dist. with expected value == half total mito (which was just doubled, so it should ~ reach original number of mito)
    num_mito_daughter1 <- rbinom(n = 1, size = length(mito_to_genome_map), prob = 0.5)
    # get the mito numbers associated with this first daughter cell
    daughter1_mito <- sample(seq(length(mito_to_genome_map)), size = num_mito_daughter1, replace = FALSE)

    
    # create a new daughter mt matrix based on the mt genomes that these mito map to
    # first get the mt genome nums associated with the selected-for mito by defining new mito_to_genome_map
   
    daughter1_mito_to_genome_map <- sapply(daughter1_mito, function(mito_num){return(mito_to_genome_map[mito_num])})
    daughter1_genome_inds <- unname(unlist(daughter1_mito_to_genome_map))
    
   
    daughter1_mt_profile <- rep_mito_mat[daughter1_genome_inds, ]
    
    # when subsetting a sparse matrix using one row index, a numeric vector is returned 
    # further processing is required to reconstruct the sparse matrix form
    if(length(daughter1_genome_inds) == 1){
      
      # only retain non-zero elements to preserve sparsity
      nonzero_j <- which(daughter1_mt_profile != 0)
      nonzero_x <- daughter1_mt_profile[nonzero_j] 
      
      # assign all to mito genome 1
      nonzero_i <- rep(1, length(nonzero_x))
      
      daughter1_mt_profile <- sparseMatrix(i = nonzero_i,
                               j = nonzero_j,
                               x = nonzero_x,
                               dims = c(1, dim(rep_mito_mat)[2]))

    }
    
    # since we change genome numbers by subsetting the parent matrix, we have to re-index each daughter mito_to_genome list
    # to ensure genome inds are found in the daughter matrix
    reindexed_daughter1_mito_to_genome_map <- reassign_genome_inds(mito_to_genome_list = daughter1_mito_to_genome_map)
    
    # get mito genome inds that belong to daughter cell 2 by taking setdifference of all genome inds and those belonging to daughter 1
    all_mito_genome_inds <- unname(unlist(mito_to_genome_map))
    daughter2_mito <- setdiff(seq(length(mito_to_genome_map)), daughter1_mito)
    daughter2_mito_to_genome_map <- sapply(daughter2_mito, function(mito_num){return(mito_to_genome_map[mito_num])})
    daughter2_genome_inds <- unname(unlist(daughter2_mito_to_genome_map))
    
    daughter2_mt_profile <- rep_mito_mat[daughter2_genome_inds, ]
    
    # perform further processing if daughter2_genome_inds has length 1
    if(length(daughter2_genome_inds) == 1){
      
      # only retain non-zero elements to preserve sparsity
      nonzero_j <- which(daughter2_mt_profile != 0)
      nonzero_x <- daughter2_mt_profile[nonzero_j] 
      
      # assign all to mito genome 1
      nonzero_i <- rep(1, length(nonzero_x))
      
      daughter2_mt_profile <- sparseMatrix(i = nonzero_i,
                                           j = nonzero_j,
                                           x = nonzero_x,
                                           dims = c(1, dim(rep_mito_mat)[2]))
      
    }
    
    
    reindexed_daughter2_mito_to_genome_map <- reassign_genome_inds(mito_to_genome_list = daughter2_mito_to_genome_map)
    
    
    
    
  }
  
  # if random mito dropout is desired, probabilistically draw here from binomial
  if(post_mitotic_mt_deletion_frac > 0){
    
    num_mito_dropout_daughter1 <- rbinom(n = 1, 
                                         size = length(reindexed_daughter1_mito_to_genome_map),
                                         prob = post_mitotic_mt_deletion_frac)
    
    if(num_mito_dropout_daughter1 > 0){
      
      # sample the mito names to determine which ones to delete (here uniform sampling)
      delete_mito_nums <- sample(names(reindexed_daughter1_mito_to_genome_map), size = num_mito_dropout_daughter1, replace = FALSE)

      # remove mt genomes that correspond to these mt numbers from the sparse matrix
      # first find row nums by mapping removed mitos to genome numbers:
      remove_genome_nums <- unname(unlist(sapply(delete_mito_nums, function(mito_num){return(reindexed_daughter1_mito_to_genome_map[mito_num])})))
      
      daughter1_mt_profile <- daughter1_mt_profile[-remove_genome_nums, ]
      
      # setting mito to genome vals to NULL removes from mito to genome map  
      reindexed_daughter1_mito_to_genome_map[delete_mito_nums] <- NULL
      reindexed_daughter1_mito_to_genome_map <- reassign_genome_inds(
        reindexed_daughter1_mito_to_genome_map
      )
    }
    
    
    
    # repeat for daughter cell 2:
    num_mito_dropout_daughter2 <- rbinom(n = 1, 
                                         size = length(reindexed_daughter2_mito_to_genome_map),
                                         prob = post_mitotic_mt_deletion_frac)
    
    if(num_mito_dropout_daughter2 > 0){
      delete_mito_nums <- sample(names(reindexed_daughter2_mito_to_genome_map), size = num_mito_dropout_daughter2, replace = FALSE)
      remove_genome_nums <- unname(unlist(sapply(delete_mito_nums, function(mito_num){return(reindexed_daughter2_mito_to_genome_map[mito_num])})))
      
      daughter2_mt_profile <- daughter2_mt_profile[-remove_genome_nums, ]
      
      # setting mito to genome vals to NULL removes from mito to genome map  
      reindexed_daughter2_mito_to_genome_map[delete_mito_nums] <- NULL
      reindexed_daughter2_mito_to_genome_map <- reassign_genome_inds(
        reindexed_daughter2_mito_to_genome_map
      )
    }
  }
  
  return_list <- list('daughter1_mt_profile' = daughter1_mt_profile,
                      'daughter2_mt_profile' = daughter2_mt_profile,
                      'daughter1_mito_to_genome_map' = reindexed_daughter1_mito_to_genome_map,
                      'daughter2_mito_to_genome_map' = reindexed_daughter2_mito_to_genome_map)
  
  return(return_list)
}

# generate the indices of the cells that will be recovered at each timepoint according to fraction of total cells captured up front. 
# this will only work correctly on terminal cell fastas. 
#' Pre-draw which terminal cells are recovered at each stopping point
#'
#' The terminal population at a stopping point is assumed to be
#' `2^floor(timepoint / cell_cycle_length)` cells, and a `ceiling`-rounded share
#' of them is sampled without replacement for each recovery rate.
#'
#' @param cell_sample_rate_vec Numeric recovery rates in [0, 1].
#' @param sim_length_stopping_points Timepoints at which output is produced.
#' @param cell_cycle_length Cell-cycle length in simulation time units.
#' @return A data frame with one row per (stopping point, rate) pair holding
#'   `cell_recovery_rate`, `num_terminal_cells`,
#'   `num_existing_nonterminal_cells`, `num_downsampled_cells`, and the list
#'   column `which_cells_recovered`.
#' @note Only meaningful for terminal-cell output; the population-size formula
#'   assumes synchronous doubling with no death. Not called by the current run.
generate_downsample_cells <- function(cell_sample_rate_vec, sim_length_stopping_points, cell_cycle_length){
  
  cell_downsample_df <- data.frame(cell_recovery_rate = numeric(),
                                   num_terminal_cells = integer(),
                                   num_existing_nonterminal_cells = integer(),
                                   num_downsampled_cells = integer(),
                                   which_cells_recovered = I(list()))
  
  num_cells_at_timepoints <- sapply(sim_length_stopping_points, function(x){
    return(2**floor((x/cell_cycle_length)))
  })
  
  # for each number of terminal cells, generate downsample inds and add to growing cell_downsample_df
  for(num_cells_ind in seq_along(num_cells_at_timepoints)){
    
    cells_at_this_timept <- num_cells_at_timepoints[num_cells_ind]
    
    for(sample_rate in cell_sample_rate_vec){
      
      num_recovered_cells <- ceiling(cells_at_this_timept*sample_rate)
      
      which_cells_recovered <- sort(sample(seq(1, cells_at_this_timept), size = num_recovered_cells, replace = FALSE)) # sorting does not hurt here since ints are independent
      
      new_row <- data.frame(cell_recovery_rate = sample_rate,
                            num_terminal_cells = cells_at_this_timept,
                            num_existing_nonterminal_cells = cells_at_this_timept-1,
                            num_downsampled_cells = num_recovered_cells,
                            which_cells_recovered = I(list(which_cells_recovered)))
      
      cell_downsample_df <- rbind(cell_downsample_df, new_row)
    }  
  }
  
  
  return(cell_downsample_df)
  
}



# generate integrations that will be selected in downsampling approaches for mt and bc:
# this will be run before the simulation bg
#' Pre-draw which integrations are recovered for each parameter combination
#'
#' @param max_ints_per_cell_vec Integer counts of integrations (or mitochondrial
#'   genomes) present per cell.
#' @param recovery_rate_vec Numeric recovery rates in [0, 1]; the recovered
#'   count is `ceiling(max_ints_per_cell * recovery_rate)`.
#' @return A data frame with one row per combination holding
#'   `max_ints_per_cell`, `recovery_rate`, `num_recovered_ints`, and the list
#'   column `which_ints_recovered`.
#' @note Not called by the current run; integration recovery is drawn per cell
#'   by `get_profiles_ints_and_umis`.
generate_downsample_integrations <- function(max_ints_per_cell_vec, recovery_rate_vec){
  
  
  integration_downsample_df <- data.frame(max_ints_per_cell = integer(),
                                          recovery_rate = numeric(),
                                          num_recovered_ints = integer(),
                                          which_ints_recovered = I(list()))
  
  for(max_ints_per_cell in max_ints_per_cell_vec){
    
    for(recovery_rate in recovery_rate_vec){
      
      num_recovered_ints <- ceiling(max_ints_per_cell*recovery_rate)
      
      which_ints_recovered <- sort(sample(seq(1, max_ints_per_cell), size = num_recovered_ints, replace = FALSE)) # sorting does not hurt here since ints are independent
      
      new_row <- data.frame(max_ints_per_cell = max_ints_per_cell,
                            recovery_rate = recovery_rate,
                            num_recovered_ints = num_recovered_ints,
                            which_ints_recovered = I(list(which_ints_recovered)))
      
      integration_downsample_df <- rbind(integration_downsample_df, new_row)
      
      
    }
    
  }
  
  return(integration_downsample_df)
  
}

#' Sample the subset of cells that receive an induction
#'
#' Exactly one of the two size arguments drives the draw: a fixed `num_cells` is
#' capped at the number of available cells, while `frac_cells` is used as a
#' binomial success probability, so the selected count varies between runs.
#'
#' @param cell_names Cell identifiers eligible for induction; coerced to
#'   character.
#' @param num_cells Fixed number to induce; must be one non-negative integer
#'   when supplied.
#' @param frac_cells Fraction to induce; must be one value in [0, 1]. Used only
#'   when `num_cells` is `NULL`.
#' @return A character vector of selected cell names, empty when nothing is
#'   selected or no cells are available.
sample_induced_cells <- function(cell_names, num_cells = NULL, frac_cells = NULL){
  cell_names <- as.character(cell_names)
  num_available <- length(cell_names)
  if(num_available == 0){
    return(character())
  }

  if(!is.null(num_cells)){
    if(length(num_cells) != 1 || !is.finite(num_cells) || num_cells < 0 ||
       num_cells %% 1 != 0){
      stop('Induction num_cells must be one non-negative integer.')
    }
    num_selected <- min(as.integer(num_cells), num_available)
  } else{
    if(length(frac_cells) != 1 || !is.finite(frac_cells) ||
       frac_cells < 0 || frac_cells > 1){
      stop('Induction frac_cells must be one fraction between zero and one.')
    }
    num_selected <- rbinom(
      n = 1,
      size = num_available,
      prob = frac_cells
    )
  }

  if(num_selected == 0){
    return(character())
  }
  sample(cell_names, size = num_selected, replace = FALSE)
}

#' Build the time-zero founder cell records
#'
#' Every founder starts alive, terminal, and parentless, and is named by its
#' integer index. All founders share the same initial mitochondrial profile,
#' barcode profile, mitochondrion-to-genome map, and heteroplasmy score; they
#' differ only in their division schedule and induction status.
#'
#' @param init_pop_size Number of founder cells; must be one positive integer.
#' @param founder_cell_type Cell type assigned to every founder.
#' @param division_points_by_founder List of eligible division timepoints, one
#'   element per founder; its length must equal `init_pop_size`.
#' @param init_incoming_mt_profile Sparse mitochondrial mutation matrix shared
#'   by the founders.
#' @param init_incoming_bc_profile Sparse barcode mutation matrix shared by the
#'   founders.
#' @param mito_to_genome_map Mitochondrion-to-genome index list shared by the
#'   founders.
#' @param initial_heteroplasmy_score Heteroplasmy score stored on each founder.
#' @param init_heteroplasmy_survive_prob Survival probability stored on each
#'   founder.
#' @param editing_induced_founders Names of the founders whose `induced_editing`
#'   field is set to `'induced_editing_params'`; all others are uninduced.
#' @param differentiation_induced_founders Names of the founders whose
#'   `induced_differentiation` field is set to `TRUE`.
#' @return A named list of founder cell records, keyed by founder name.
initialize_founder_population <- function(init_pop_size,
                                          founder_cell_type,
                                          division_points_by_founder,
                                          init_incoming_mt_profile,
                                          init_incoming_bc_profile,
                                          mito_to_genome_map,
                                          initial_heteroplasmy_score,
                                          init_heteroplasmy_survive_prob,
                                          editing_induced_founders = character(),
                                          differentiation_induced_founders = character()){
  if(length(init_pop_size) != 1 || !is.finite(init_pop_size) ||
     init_pop_size < 1 || init_pop_size %% 1 != 0){
    stop('num_init_cells must be one positive integer.')
  }
  init_pop_size <- as.integer(init_pop_size)
  founder_names <- as.character(seq_len(init_pop_size))
  if(length(division_points_by_founder) != init_pop_size){
    stop('division_points_by_founder must contain one entry per founder cell.')
  }

  setNames(lapply(seq_along(founder_names), function(index){
    founder_name <- founder_names[index]
    list(
      linstring = founder_name,
      celltype = founder_cell_type,
      birth_time = 0,
      death_time = NA,
      parent = NULL,
      descendants = c(),
      alive = TRUE,
      terminal = TRUE,
      elig_div_points = division_points_by_founder[[index]],
      incoming_mt_profiles = init_incoming_mt_profile,
      incoming_bc_profiles = init_incoming_bc_profile,
      mito_to_genome_map = mito_to_genome_map,
      heteroplasmy_score = initial_heteroplasmy_score,
      heteroplasmy_survive_prob = init_heteroplasmy_survive_prob,
      induced_editing = if(founder_name %in% editing_induced_founders){
        'induced_editing_params'
      } else{
        'uninduced_editing_params'
      },
      induced_differentiation =
        founder_name %in% differentiation_induced_founders
    )
  }), founder_names)
}

# replace pos_er_list in here ........
#' Create the founder population and start the parallel worker cluster
#'
#' Builds the global `poss_times` grid, draws each founder's division schedule
#' from an exponential waiting-time model, applies any induction scheduled for
#' time zero, assembles the founder records, then starts a `parallel` cluster
#' and exports the mutation helpers along with every parameter collection the
#' workers read. Most arguments exist only to be exported.
#'
#' @param num_clusters Number of worker processes to start.
#' @param init_pop_size Number of founder cells; must be one positive integer.
#' @param sim_length Simulation horizon; bounds the time grid and the division
#'   schedules.
#' @param cell_type_cell_cycle_length Named list of mean cell-cycle lengths per
#'   cell type.
#' @param num_rows_mt Rows of the founder mitochondrial profile.
#' @param num_cols_mt Columns of the mitochondrial profile (genome length).
#' @param num_rows_bc Rows of the barcode profile (maximum integrations).
#' @param num_cols_bc Columns of the barcode profile (barcode length).
#' @param time_inc Simulation time increment; sets the spacing of `poss_times`.
#' @param init_incoming_mt_profile Sparse founder mitochondrial mutation matrix.
#' @param cell_type_basepos_bc_nontarget_transition_probs Per-cell-type,
#'   per-induction background barcode transition probability by position.
#' @param cell_type_basepos_bc_nontarget_transversion_probs As above, for
#'   background barcode transversions.
#' @param cell_type_basepos_bc_nontarget_insertion_probs As above, for
#'   background barcode insertions.
#' @param cell_type_basepos_bc_nontarget_deletion_probs As above, for background
#'   barcode deletions.
#' @param cell_type_basepos_mt_nontarget_transition_probs Background
#'   mitochondrial transition probability by position.
#' @param cell_type_basepos_mt_nontarget_transversion_probs Background
#'   mitochondrial transversion probability by position.
#' @param cell_type_basepos_mt_nontarget_insertion_probs Background
#'   mitochondrial insertion probability by position.
#' @param cell_type_basepos_mt_nontarget_deletion_probs Background mitochondrial
#'   deletion probability by position.
#' @param cell_type_basepos_bc_target_transition_probs Barcode target-site
#'   transition probability by position.
#' @param cell_type_basepos_bc_target_transversion_probs Barcode target-site
#'   transversion probability by position.
#' @param cell_type_basepos_bc_target_insertion_probs Barcode target-site
#'   insertion probability by position.
#' @param cell_type_basepos_bc_target_deletion_probs Barcode target-site
#'   deletion probability by position.
#' @param cell_type_mt_sub_prob_mat Per-cell-type mitochondrial substitution
#'   probability matrices.
#' @param cell_type_bc_sub_prob_mat Per-cell-type barcode substitution
#'   probability matrices.
#' @param cell_type_death_probs Per-cell-type death probability per timepoint.
#' @param uninduced_tm_list Cell-type transition probabilities without
#'   differentiation induction.
#' @param induced_tm_list Cell-type transition probabilities under
#'   differentiation induction.
#' @param differentiation_induction_timepoint Timepoint at which differentiation
#'   induction fires.
#' @param editing_induction_timepoint Timepoint at which editing induction
#'   fires.
#' @param differentiation_induction_num_cells Fixed number of cells to induce
#'   for differentiation.
#' @param editing_induction_num_cells Fixed number of cells to induce for
#'   editing.
#' @param differentiation_induction_frac_cells Fraction of cells to induce for
#'   differentiation.
#' @param editing_induction_frac_cells Fraction of cells to induce for editing.
#' @param custom_savename Output filename stem exported to the workers.
#' @param forced_transversions Whether target substitutions are forced to a
#'   fixed destination base.
#' @param be_target_to_int Integer code of the base editor's destination base.
#' @param sim_length_stopping_points Timepoints at which output is written.
#' @param founder_cell_type Cell type assigned to the founders.
#' @param poss_fasta_types Which cell sets get FASTA output (`'all_cells'`,
#'   `'terminal'`).
#' @param include_var_pos_fasta Whether variable-position FASTAs are written.
#' @param interdeletion_dropout_radius Radius used when dropping sites between
#'   paired deletions.
#' @param interdeletion_dropout_prob Probability used for that dropout.
#' @param poss_recon_modals Modalities to simulate and reconstruct (`'mt'`,
#'   `'bc'`).
#' @param mito_to_genome_map Founder mitochondrion-to-genome index list.
#' @param fusion_events_per_mito_per_division Poisson fusion rate per
#'   mitochondrion per division.
#' @param split_events_per_mito_per_division Poisson fission rate per
#'   mitochondrion per division.
#' @param post_mitotic_mt_deletion_frac Per-mitochondrion post-division dropout
#'   probability.
#' @param positive_score_weight Mixture weight on the beneficial severity
#'   component.
#' @param hetero_sd Standard deviation of the severity score mixture.
#' @param heteroplasmy_severity_score_list Variant name to severity score map.
#' @param init_heteroplasmy_survive_prob Founder survival probability anchor.
#' @param mito_inheritance_pattern Mitochondrial allocation rule; only
#'   `'random'` is implemented.
#' @param heteroplasmy_variant_fractions Variant name to variant fraction map.
#' @param initial_heteroplasmy_score Founder heteroplasmy score.
#' @param init_num_mito_genomes Number of mitochondrial genomes in the founder.
#' @param ind_to_prime_seq_int_map Target position to integer prime-editing
#'   guide map.
#' @param ind_to_prime_seq_nuc_map Target position to nucleotide prime-editing
#'   guide map.
#' @param prime_editing_system Whether prime editing is enabled.
#' @param close_nuc_window_after_edit Whether a nuclease window closes once one
#'   of its positions is edited.
#' @param close_transition_window_after_edit Whether a base-editor transition
#'   window closes after an edit.
#' @param close_transversion_window_after_edit Whether a base-editor
#'   transversion window closes after an edit.
#' @param be_target_to_window_ind_list Base-editor position to window-name map.
#' @param nuc_target_to_window_ind_list Nuclease position to window-name map.
#' @param be_window_to_target_ind_list Base-editor window-name to positions map.
#' @param nuc_window_to_target_ind_list Nuclease window-name to positions map.
#' @param consider_cell_heteroplasmy_scores Whether heteroplasmy scores drive
#'   cell survival.
#' @return The named list of founder cell records that starts the simulation.
#' @section Side effects: Assigns the globals `poss_times`, `one_cluster`,
#'   `already_assigned_editing_induction`, `already_assigned_diff_induction`,
#'   and `cluster_startup_total`. Starts a `parallel` cluster, seeds its RNG
#'   streams from `input_args$random_seed`, and restores the main process's
#'   `.Random.seed` afterwards so worker seeding does not perturb it.
setup_sim <- function(num_clusters, 
                      init_pop_size, 
                      sim_length, 
                      cell_type_cell_cycle_length,
                      num_rows_mt, 
                      num_cols_mt, 
                      num_rows_bc, 
                      num_cols_bc, 
                      time_inc,
                      init_incoming_mt_profile,
                      cell_type_basepos_bc_nontarget_transition_probs,
                      cell_type_basepos_bc_nontarget_transversion_probs,
                      cell_type_basepos_bc_nontarget_insertion_probs,
                      cell_type_basepos_bc_nontarget_deletion_probs,
                      cell_type_basepos_mt_nontarget_transition_probs,
                      cell_type_basepos_mt_nontarget_transversion_probs,
                      cell_type_basepos_mt_nontarget_insertion_probs,
                      cell_type_basepos_mt_nontarget_deletion_probs,
                      cell_type_basepos_bc_target_transition_probs,
                      cell_type_basepos_bc_target_transversion_probs,
                      cell_type_basepos_bc_target_insertion_probs,
                      cell_type_basepos_bc_target_deletion_probs,
                      cell_type_mt_sub_prob_mat,
                      cell_type_bc_sub_prob_mat,
                      cell_type_death_probs,
                      uninduced_tm_list,
                      induced_tm_list,
                      differentiation_induction_timepoint,
                      editing_induction_timepoint,
                      differentiation_induction_num_cells,
                      editing_induction_num_cells,
                      differentiation_induction_frac_cells,
                      editing_induction_frac_cells,
                      custom_savename,
                      forced_transversions,
                      be_target_to_int,
                      sim_length_stopping_points,
                      founder_cell_type,
                      poss_fasta_types,
                      include_var_pos_fasta,
                      interdeletion_dropout_radius,
                      interdeletion_dropout_prob,
                      poss_recon_modals,
                      mito_to_genome_map,
                      fusion_events_per_mito_per_division,
                      split_events_per_mito_per_division,
                      post_mitotic_mt_deletion_frac,
                      positive_score_weight,
                      hetero_sd,
                      heteroplasmy_severity_score_list,
                      init_heteroplasmy_survive_prob,
                      mito_inheritance_pattern,
                      heteroplasmy_variant_fractions,
                      initial_heteroplasmy_score,
                      init_num_mito_genomes,
                      ind_to_prime_seq_int_map,
                      ind_to_prime_seq_nuc_map,
                      prime_editing_system,
                      close_nuc_window_after_edit,
                      close_transition_window_after_edit,
                      close_transversion_window_after_edit,
                      be_target_to_window_ind_list,
                      nuc_target_to_window_ind_list,
                      be_window_to_target_ind_list,
                      nuc_window_to_target_ind_list,
                      consider_cell_heteroplasmy_scores){

  if(length(init_pop_size) != 1 || !is.finite(init_pop_size) ||
     init_pop_size < 1 || init_pop_size %% 1 != 0){
    stop('num_init_cells must be one positive integer.')
  }
  init_pop_size <- as.integer(init_pop_size)
  poss_times <<- seq(0, sim_length, time_inc)
  
  init_incoming_mt_profile <- init_incoming_mt_profile
  init_incoming_bc_profile <- sparseMatrix(i = c(1), j = c(1), x = c(0L),
                                           dims = c(num_rows_bc, num_cols_bc))
  
  
  # generate initial elig_div_points by drawing from an exponential distribution:
  #' Internal: draw a cell's future division timepoints
  #'
  #' @param sim_length Horizon past which division points are discarded.
  #' @param cc_length Mean cell-cycle length, used as the mean of the
  #'   exponential waiting time between successive divisions.
  #' @param current_timepoint Time from which to start accumulating waiting
  #'   times.
  #' @return An increasing numeric vector of division timepoints, all at or
  #'   below `sim_length`.
  get_future_div_points <- function(sim_length, cc_length, current_timepoint){
    
    # empty vector to which future div points will be appended
    elig_div_points <- numeric()
    
    while(current_timepoint < sim_length){
      next_cc_length <- rexp(n = 1, rate = 1/cc_length)
      next_timept <- current_timepoint + next_cc_length
      elig_div_points <- c(elig_div_points, next_timept)
      current_timepoint <- next_timept
    }
    
    # ensure all timepoints above sim_length are removed:
    elig_div_points <- elig_div_points[elig_div_points <= sim_length]
    
    return(elig_div_points)
    
  }
  
  founder_names <- as.character(seq_len(init_pop_size))
  founder_division_points <- lapply(founder_names, function(founder_name){
    get_future_div_points(
      sim_length = sim_length,
      cc_length = cell_type_cell_cycle_length[[founder_cell_type]],
      current_timepoint = 0
    )
  })

  editing_starts_at_zero <-
    as.numeric(input_args$editing_induction$timepoint) == 0
  differentiation_starts_at_zero <-
    as.numeric(input_args$differentiation_induction$timepoint) == 0

  editing_induced_founders <- if(editing_starts_at_zero){
    sample_induced_cells(
      founder_names,
      num_cells = input_args$editing_induction$num_cells,
      frac_cells = input_args$editing_induction$frac_cells
    )
  } else{
    character()
  }
  differentiation_induced_founders <- if(differentiation_starts_at_zero){
    sample_induced_cells(
      founder_names,
      num_cells = input_args$differentiation_induction$num_cells,
      frac_cells = input_args$differentiation_induction$frac_cells
    )
  } else{
    character()
  }

  # These flags record that assignment happened, including a valid assignment
  # of zero cells, rather than whether any cell happened to be selected.
  already_assigned_editing_induction <<- editing_starts_at_zero
  already_assigned_diff_induction <<- differentiation_starts_at_zero

  cell_population <- initialize_founder_population(
    init_pop_size = init_pop_size,
    founder_cell_type = founder_cell_type,
    division_points_by_founder = founder_division_points,
    init_incoming_mt_profile = init_incoming_mt_profile,
    init_incoming_bc_profile = init_incoming_bc_profile,
    mito_to_genome_map = mito_to_genome_map,
    initial_heteroplasmy_score = initial_heteroplasmy_score,
    init_heteroplasmy_survive_prob = init_heteroplasmy_survive_prob,
    editing_induced_founders = editing_induced_founders,
    differentiation_induced_founders = differentiation_induced_founders
  )

   
  cluster_startup_start <- Sys.time()
  one_cluster <<- makeCluster(num_clusters)
  main_rng_state <- get('.Random.seed', envir = .GlobalEnv)
  clusterSetRNGStream(
    cl = one_cluster,
    iseed = input_args$random_seed
  )
  assign('.Random.seed', main_rng_state, envir = .GlobalEnv)

  clusterEvalQ(cl = one_cluster, {
    suppressPackageStartupMessages(library(Matrix))
    suppressPackageStartupMessages(library(data.table))
  })
  clusterExport(cl = one_cluster, c('perform_all_mt_mutations', 'perform_all_bc_mutations', 'transition_func', 'transversion_func',
                                    'insertion_func', 'deletion_func', 'bases', 'transition_matches',
                                    'transversion_matches', 'baseline_seq_ints_mt', 'baseline_seq_ints_bc',
                                    'baseline_seq_nucs_mt', 'baseline_seq_nucs_bc',
                                    'num_deletable_bases', 'perform_deletion', 'all_deletions_one_mat',
                                    'num_rows_bc', 'num_cols_bc', 'num_rows_mt', 'num_cols_mt',
                                    'init_pop_size', 
                                    'sim_length',
                                    'cell_type_basepos_bc_nontarget_transition_probs',
                                    'cell_type_basepos_bc_nontarget_transversion_probs',
                                    'cell_type_basepos_bc_nontarget_insertion_probs',
                                    'cell_type_basepos_bc_nontarget_deletion_probs',
                                    'cell_type_basepos_mt_nontarget_transition_probs',
                                    'cell_type_basepos_mt_nontarget_transversion_probs',
                                    'cell_type_basepos_mt_nontarget_insertion_probs',
                                    'cell_type_basepos_mt_nontarget_deletion_probs',
                                    'cell_type_basepos_bc_target_transition_probs',
                                    'cell_type_basepos_bc_target_transversion_probs',
                                    'cell_type_basepos_bc_target_insertion_probs',
                                    'cell_type_basepos_bc_target_deletion_probs',
                                    'cell_type_mt_sub_prob_mat',
                                    'cell_type_bc_sub_prob_mat',
                                    'cell_type_cell_cycle_length', 
                                    'cell_type_death_probs',
                                    'uninduced_tm_list',
                                    'induced_tm_list',
                                    'differentiation_induction_timepoint',
                                    'editing_induction_timepoint',
                                    'differentiation_induction_num_cells',
                                    'editing_induction_num_cells',
                                    'differentiation_induction_frac_cells',
                                    'editing_induction_frac_cells',
                                    'already_assigned_editing_induction',
                                    'already_assigned_diff_induction',
                                    'old_cells_at_timept', 
                                    'forced_transversions', 'be_target_to_int',
                                    'get_background_edit_inds',
                                    'non_uniform_editing',
                                    'sim_length_stopping_points', 
                                    'mito_dynamics',
                                    'unique_run_id',
                                    'poss_fasta_types', 'include_var_pos_fasta',
                                    'founder_cell_type',
                                    'get_one_cell_sequence', 'ins_to_charvec', # for writing to fastas
                                    'bc_int_umis',
                                    'interdeletion_dropout_radius',
                                    'interdeletion_dropout_prob',
                                    'scoremat_collapse_deletions',
                                    'combine_mt_bc',
                                    'binarize_mutation_scores',
                                    'mt_allelic_fraction_thresholds',
                                    'poss_recon_modals',
                                    'fusion_events_per_mito_per_division',
                                    'split_events_per_mito_per_division',
                                    'get_future_div_points',
                                    'post_mitotic_mt_deletion_frac',
                                    'heteroplasmy_severity_score_list',
                                    'heteroplasmy_variant_fractions',
                                    'positive_score_weight',
                                    'hetero_sd',
                                    'logistic_prob_survive_given_score',
                                    'init_heteroplasmy_survive_prob',
                                    'mito_inheritance_pattern',
                                    'initial_heteroplasmy_score',
                                    'draw_severity_scores',
                                    'reassign_genome_inds',
                                    'init_num_mito_genomes',
                                    'ind_to_prime_seq_int_map',
                                    'ind_to_prime_seq_nuc_map',
                                    'prime_editing_system',
                                    'close_nuc_window_after_edit',
                                    'close_transition_window_after_edit',
                                    'close_transversion_window_after_edit',
                                    'be_target_to_window_ind_list',
                                    'nuc_target_to_window_ind_list',
                                    'be_window_to_target_ind_list',
                                    'nuc_window_to_target_ind_list',
                                    'filter_elig_ints_by_edit_window',
                                    'consider_cell_heteroplasmy_scores',
                                    'all_processes_at_stopping_point',
                                    'custom_savename',
                                    'poss_times',
                                    'sim_length_stopping_points'
                                    ),
                envir = environment())
  cluster_startup_end <- Sys.time()
  cluster_startup_total <<- difftime(cluster_startup_end, cluster_startup_start, units = 'secs')
  return(cell_population)
}

#' Advance the whole cell population by one simulation timepoint
#'
#' Runs the timepoint in a fixed order: induction assignment, division, death,
#' then mutation. Cells whose next eligible division point has been reached
#' divide on the worker cluster -- mitochondria pass through `mito_dynamics`,
#' daughter cell types are drawn from the induced or uninduced transition
#' matrix, and each daughter gets a fresh exponential division schedule; a
#' daughter that inherits no mitochondrion is born dead. Death then combines a
#' heteroplasmy score draw (only when `consider_cell_heteroplasmy_scores`), a
#' per-cell-type random draw, and having lost every mitochondrion. Surviving
#' cells finally have barcode and mitochondrial mutations applied on the
#' workers, and a timepoint that is also a stopping point triggers the full
#' output pass.
#'
#' @param timepoint The simulation time being advanced to.
#' @param sim_length Simulation horizon, used when drawing daughter division
#'   schedules.
#' @param cell_population Named list of cell records at the start of the
#'   timepoint.
#' @param cell_type_basepos_bc_nontarget_transition_probs Background barcode
#'   transition probability by position, per cell type and induction state.
#' @param cell_type_basepos_bc_nontarget_transversion_probs As above, for
#'   background barcode transversions.
#' @param cell_type_basepos_bc_nontarget_insertion_probs As above, for
#'   background barcode insertions.
#' @param cell_type_basepos_bc_nontarget_deletion_probs As above, for background
#'   barcode deletions.
#' @param cell_type_basepos_mt_nontarget_transition_probs Background
#'   mitochondrial transition probability by position.
#' @param cell_type_basepos_mt_nontarget_transversion_probs Background
#'   mitochondrial transversion probability by position.
#' @param cell_type_basepos_mt_nontarget_insertion_probs Background
#'   mitochondrial insertion probability by position.
#' @param cell_type_basepos_mt_nontarget_deletion_probs Background mitochondrial
#'   deletion probability by position.
#' @param cell_type_basepos_bc_target_transition_probs Barcode target-site
#'   transition probability by position.
#' @param cell_type_basepos_bc_target_transversion_probs Barcode target-site
#'   transversion probability by position.
#' @param cell_type_basepos_bc_target_insertion_probs Barcode target-site
#'   insertion probability by position.
#' @param cell_type_basepos_bc_target_deletion_probs Barcode target-site
#'   deletion probability by position.
#' @param cell_type_mt_sub_prob_mat Per-cell-type mitochondrial substitution
#'   probability matrices.
#' @param cell_type_bc_sub_prob_mat Per-cell-type barcode substitution
#'   probability matrices.
#' @param cell_type_death_probs Per-cell-type death probability per timepoint.
#' @param uninduced_tm_list Cell-type transition probabilities without
#'   differentiation induction.
#' @param induced_tm_list Cell-type transition probabilities under
#'   differentiation induction.
#' @param differentiation_induction_timepoint Timepoint at or after which
#'   differentiation induction is assigned.
#' @param editing_induction_timepoint Timepoint at or after which editing
#'   induction is assigned.
#' @param differentiation_induction_num_cells Fixed number of cells to induce
#'   for differentiation.
#' @param editing_induction_num_cells Fixed number of cells to induce for
#'   editing.
#' @param differentiation_induction_frac_cells Fraction of cells to induce for
#'   differentiation.
#' @param editing_induction_frac_cells Fraction of cells to induce for editing.
#' @param already_assigned_editing_induction Whether editing induction has
#'   already been assigned in an earlier timepoint.
#' @param already_assigned_diff_induction Whether differentiation induction has
#'   already been assigned in an earlier timepoint.
#' @param forced_transversions Whether target substitutions are forced to a
#'   fixed destination base.
#' @param be_target_to_int Integer code of the base editor's destination base.
#' @param unique_run_id Run identifier used in every output path.
#' @param interdeletion_dropout_radius Radius used when dropping sites between
#'   paired deletions.
#' @param interdeletion_dropout_prob Probability used for that dropout.
#' @param scoremat_collapse_deletions Whether score matrices collapse runs of
#'   deleted positions.
#' @param combine_mt_bc Whether mitochondrial and barcode score matrices are
#'   joined at stopping points.
#' @param binarize_mutation_scores Whether score matrices are binarized.
#' @param mt_allelic_fraction_thresholds Allelic fraction thresholds applied to
#'   mitochondrial scores.
#' @param poss_recon_modals Modalities to simulate (`'mt'`, `'bc'`).
#' @param fusion_events_per_mito_per_division Poisson fusion rate.
#' @param split_events_per_mito_per_division Poisson fission rate.
#' @param post_mitotic_mt_deletion_frac Per-mitochondrion post-division dropout
#'   probability.
#' @param heteroplasmy_variant_fractions Variant name to variant fraction map.
#' @param heteroplasmy_severity_score_list Variant name to severity score map.
#' @param positive_score_weight Mixture weight on the beneficial severity
#'   component.
#' @param hetero_sd Standard deviation of the severity score mixture.
#' @param mito_inheritance_pattern Mitochondrial allocation rule; only
#'   `'random'` is implemented.
#' @param init_num_mito_genomes Number of mitochondrial genomes in the founder,
#'   used to name mitochondrial output combinations.
#' @param ind_to_prime_seq_int_map Target position to integer prime-editing
#'   guide map.
#' @param ind_to_prime_seq_nuc_map Target position to nucleotide prime-editing
#'   guide map.
#' @param prime_editing_system Whether prime editing is enabled.
#' @param close_nuc_window_after_edit Whether a nuclease window closes once one
#'   of its positions is edited.
#' @param close_transition_window_after_edit Whether a base-editor transition
#'   window closes after an edit.
#' @param close_transversion_window_after_edit Whether a base-editor
#'   transversion window closes after an edit.
#' @param be_target_to_window_ind_list Base-editor position to window-name map.
#' @param nuc_target_to_window_ind_list Nuclease position to window-name map.
#' @param be_window_to_target_ind_list Base-editor window-name to positions map.
#' @param nuc_window_to_target_ind_list Nuclease window-name to positions map.
#' @param consider_cell_heteroplasmy_scores Whether heteroplasmy scores drive
#'   cell survival.
#' @param poss_times Full simulation time grid.
#' @param custom_savename Output filename stem.
#' @param sim_length_stopping_points Timepoints at which output is written.
#' @param t Index of `timepoint` within `poss_times`, forwarded as the relative
#'   timepoint of the output pass.
#' @return The updated cell population list, including the new daughter cells.
#' @section Side effects: Writes
#'   `output/induction_details/<run_id>/differentiation_induction_df.csv` and
#'   `editing_induction_df.csv` on the timepoint an induction fires; sets the
#'   globals `already_assigned_diff_induction`,
#'   `already_assigned_editing_induction`, and each cell's `heteroplasmy_score`;
#'   dispatches work to the `one_cluster` global; and calls
#'   `all_processes_at_stopping_point` at stopping points.
#' @note The whole R session is terminated with `quit(save = 'no')` if every
#'   cell is dead at this timepoint, so no tree is reconstructed for that run.
multi_core_func <- function(timepoint, 
                            sim_length,
                            cell_population,
                            cell_type_basepos_bc_nontarget_transition_probs, 
                            cell_type_basepos_bc_nontarget_transversion_probs, 
                            cell_type_basepos_bc_nontarget_insertion_probs, 
                            cell_type_basepos_bc_nontarget_deletion_probs,
                            cell_type_basepos_mt_nontarget_transition_probs, 
                            cell_type_basepos_mt_nontarget_transversion_probs, 
                            cell_type_basepos_mt_nontarget_insertion_probs, 
                            cell_type_basepos_mt_nontarget_deletion_probs, 
                            cell_type_basepos_bc_target_transition_probs, 
                            cell_type_basepos_bc_target_transversion_probs, 
                            cell_type_basepos_bc_target_insertion_probs, 
                            cell_type_basepos_bc_target_deletion_probs,
                            cell_type_mt_sub_prob_mat, 
                            cell_type_bc_sub_prob_mat, 
                            cell_type_death_probs,
                            uninduced_tm_list,
                            induced_tm_list,
                            differentiation_induction_timepoint,
                            editing_induction_timepoint,
                            differentiation_induction_num_cells,
                            editing_induction_num_cells,
                            differentiation_induction_frac_cells,
                            editing_induction_frac_cells,
                            already_assigned_editing_induction,
                            already_assigned_diff_induction,
                            forced_transversions,
                            be_target_to_int,
                            unique_run_id,
                            interdeletion_dropout_radius,
                            interdeletion_dropout_prob,
                            scoremat_collapse_deletions,
                            combine_mt_bc,
                            binarize_mutation_scores,
                            mt_allelic_fraction_thresholds,
                            poss_recon_modals,
                            fusion_events_per_mito_per_division,
                            split_events_per_mito_per_division,
                            post_mitotic_mt_deletion_frac,
                            heteroplasmy_variant_fractions,
                            heteroplasmy_severity_score_list,
                            positive_score_weight,
                            hetero_sd,
                            mito_inheritance_pattern,
                            init_num_mito_genomes,
                            ind_to_prime_seq_int_map,
                            ind_to_prime_seq_nuc_map,
                            prime_editing_system,
                            close_nuc_window_after_edit,
                            close_transition_window_after_edit,
                            close_transversion_window_after_edit,
                            be_target_to_window_ind_list,
                            nuc_target_to_window_ind_list,
                            be_window_to_target_ind_list,
                            nuc_window_to_target_ind_list,
                            consider_cell_heteroplasmy_scores,
                            poss_times,
                            custom_savename,
                            sim_length_stopping_points,
                            t){
  
  
  ######## DIVIDE, THEN DIE, THEN MUTATE

  # adding another check for alive and terminal cells that can become induced
  # at the beginning of this timepoint:
  cells_alive_here_bool_list <- sapply(cell_population,
                                       function(cell) cell$terminal & cell$alive,
                                       USE.NAMES = TRUE)
  
  print(paste0('Pop size == ', sum(cells_alive_here_bool_list)))
  
  cell_names_alive_here <- names(cells_alive_here_bool_list)[cells_alive_here_bool_list]
  
  # assign differentiation induction cells if necessary
  if(timepoint >= differentiation_induction_timepoint){
    
    # only have to assign induction statuses in the first timepoint after we cross the induction timepoint:
    if(already_assigned_diff_induction == FALSE){
      print('Inducing differentiation')
      now_induced_cells <- sample_induced_cells(
        cell_names_alive_here,
        num_cells = differentiation_induction_num_cells,
        frac_cells = differentiation_induction_frac_cells
      )
      
      for(linstring in now_induced_cells){
        cell_population[[linstring]]$induced_differentiation <- TRUE
      }
      already_assigned_diff_induction <<- TRUE
      
      # write which cells were selected to be induced (and when) to csv:
      differentiation_induction_df <- data.frame(cbind(now_induced_cells, 
                                                       rep(timepoint, length(now_induced_cells))))
      
      colnames(differentiation_induction_df) <- c('linstring', 'differentiation_induced_timepoint')
      if(!dir.exists(file.path('output', 'induction_details', unique_run_id))){
        dir.create(file.path('output', 'induction_details', unique_run_id), recursive = TRUE)
      }
      write.csv(differentiation_induction_df, 
                file.path('output', 'induction_details', unique_run_id, 'differentiation_induction_df.csv'))
    }
  }
  
  # assign editing induction cells if necessary
  if(timepoint >= editing_induction_timepoint){
    
    # only have to assign induction statuses in the first timepoint after we cross the induction timepoint:
    if(already_assigned_editing_induction == FALSE){
      print('Inducing editing')
      now_induced_cells <- sample_induced_cells(
        cell_names_alive_here,
        num_cells = editing_induction_num_cells,
        frac_cells = editing_induction_frac_cells
      )
      
      for(linstring in now_induced_cells){
        cell_population[[linstring]]$induced_editing <- 'induced_editing_params'
      }
      
      already_assigned_editing_induction <<- TRUE
      
      # write which cells were selected to be induced (and when) to csv:
      editing_induction_df <- data.frame(cbind(now_induced_cells, 
                                                       rep(timepoint, length(now_induced_cells))))
      
      colnames(editing_induction_df) <- c('linstring', 'editing_induced_timepoint')
      if(!dir.exists(file.path('output', 'induction_details', unique_run_id))){
        dir.create(file.path('output', 'induction_details', unique_run_id), recursive = TRUE)
      }
      write.csv(editing_induction_df, 
                file.path('output', 'induction_details', unique_run_id, 'editing_induction_df.csv'))
    }
  }
  
  cells_dividing_here_bool_list <- sapply(cell_population, 
                          function(cell){
                            return(length(cell$elig_div_points) > 0 & timepoint >= cell$elig_div_points[1] & cell$terminal & cell$alive)
                          },
                          USE.NAMES = TRUE)
  cell_names_dividing_here <- names(cells_dividing_here_bool_list)[cells_dividing_here_bool_list]
  
  # get lineage strings corresponding to the cells dividing here
  
  
  num_cells_dividing_here <- length(cell_names_dividing_here)
  
  new_cell_list <- parLapply(cl = one_cluster, X = cell_names_dividing_here, 
                             fun = function(cell_name){
                               
                               this_cell_type <- cell_population[[cell_name]]$celltype
                               
                               mt_dynamics_res <- mito_dynamics(mito_to_genome_map = cell_population[[cell_name]]$mito_to_genome_map, 
                                                                incoming_mito_mat = cell_population[[cell_name]]$incoming_mt_profiles, 
                                                                fusion_events_per_mito_per_division = fusion_events_per_mito_per_division,
                                                                split_events_per_mito_per_division = split_events_per_mito_per_division,
                                                                inheritance_pattern = mito_inheritance_pattern,
                                                                post_mitotic_mt_deletion_frac = post_mitotic_mt_deletion_frac)
                               
                               
                               new_mitoprofiles_1 <- mt_dynamics_res[['daughter1_mt_profile']]
                               new_mitoprofiles_2 <- mt_dynamics_res[['daughter2_mt_profile']]
                               daughter1_mito_to_genome_map <- mt_dynamics_res[['daughter1_mito_to_genome_map']]
                               daughter2_mito_to_genome_map <- mt_dynamics_res[['daughter2_mito_to_genome_map']]
                               
                               # generate daughter cell lineage strings
                               daughter_cell_linstrings <- paste(cell_name, seq(1,2), sep = '_')
                               
                               differentiation_induced <- cell_population[[cell_name]]$induced_differentiation
                               editing_induced <- cell_population[[cell_name]]$induced_editing
                               
                               # generate new daughter cell types according to whether parent has been induced
                               if(differentiation_induced){
                                 daughter_cell_types <- sample(names(induced_tm_list[[this_cell_type]]), 
                                                               size = 2, replace = TRUE, 
                                                               prob = as.numeric(induced_tm_list[[this_cell_type]]))
                               } else if(!differentiation_induced){
                                 daughter_cell_types <- sample(names(uninduced_tm_list[[this_cell_type]]), 
                                                               size = 2, replace = TRUE, 
                                                               prob = as.numeric(uninduced_tm_list[[this_cell_type]]))
                               }
                               
                               # update parent cell with descendant info and to reflect changed terminal status
                               # in order for these changes to be made, have to return this modified copy from the parallel workers, then subsequently overwrite existing vals
                               cell_population[[cell_name]]$descendants <- daughter_cell_linstrings
                               cell_population[[cell_name]]$terminal <- FALSE
                               cell_population[[cell_name]]$alive <- FALSE
                               cell_population[[cell_name]]$dwell_time <- timepoint - cell_population[[cell_name]]$birth_time
                               cell_population[[cell_name]]$divide_time <- timepoint
                               
                               # update the elig div points by removing the first timepoint so that we can continue to access first timepoint for comparisons later on
                               cell_population[[cell_name]]$elig_div_points <- cell_population[[cell_name]]$elig_div_points[2:length(cell_population[[cell_name]]$elig_div_points)]
                               
                               # new line ... 
                               cell_population[[cell_name]]$death_time <- timepoint
                               
                               
                               
                               # probabilistically draw new elig div points for each daughter based on cell type-specific cell cycle lenghts
                               daughter1_elig_div_points <- get_future_div_points(sim_length = sim_length, 
                                                                                  cc_length = cell_type_cell_cycle_length[[daughter_cell_types[1]]], 
                                                                                  current_timepoint = timepoint)
                               
                               daughter2_elig_div_points <- get_future_div_points(sim_length = sim_length, 
                                                                                  cc_length = cell_type_cell_cycle_length[[daughter_cell_types[2]]], 
                                                                                  current_timepoint = timepoint)
                               
                               daughter_cells <- list()
                               
                               # if a daughter cell is assigned zero mito, label it as dead upon creation. 
                               if(length(daughter1_mito_to_genome_map) > 0){
                                 daughter1_alive <- TRUE
                               } else{
                                 daughter1_alive <- FALSE
                     
                               }
                               if(length(daughter2_mito_to_genome_map) > 0){
                                 daughter2_alive <- TRUE
                               } else{
                                 daughter2_alive <- FALSE
               
                               }
                               daughter_cells[[daughter_cell_linstrings[1]]] <- list('linstring' = daughter_cell_linstrings[1],
                                                                                    'celltype' = daughter_cell_types[1],
                                                                                    'birth_time' = timepoint,
                                                                                    'death_time' = NA,
                                                                                    'parent' = cell_name,
                                                                                    'descendants' = c(),
                                                                                    'alive' = daughter1_alive,
                                                                                    'terminal' = TRUE,
                                                                                    'elig_div_points' = daughter1_elig_div_points,
                                                                                    'incoming_mt_profiles' = new_mitoprofiles_1,
                                                                                    'incoming_bc_profiles' = cell_population[[cell_name]]$incoming_bc_profiles,
                                                                                    'mito_to_genome_map'= daughter1_mito_to_genome_map,
                                                                                    'induced_editing' = editing_induced,
                                                                                    'induced_differentiation' = differentiation_induced)
                               daughter_cells[[daughter_cell_linstrings[2]]] <- list('linstring' = daughter_cell_linstrings[2],
                                                                                 'celltype' = daughter_cell_types[2],
                                                                                 'birth_time' = timepoint,
                                                                                 'death_time' = NA,
                                                                                 'parent' = cell_name,
                                                                                 'descendants' = c(),
                                                                                 'alive' = daughter2_alive,
                                                                                 'terminal' = TRUE,
                                                                                 'elig_div_points' = daughter2_elig_div_points,
                                                                                 'incoming_mt_profiles' = new_mitoprofiles_2,
                                                                                 'incoming_bc_profiles' = cell_population[[cell_name]]$incoming_bc_profiles,
                                                                                 'mito_to_genome_map'= daughter2_mito_to_genome_map,
                                                                                 'induced_editing' = editing_induced,
                                                                                 'induced_differentiation' = differentiation_induced)
                               
                                return_list <- list()
                                return_list[['updated_parent']] <- cell_population[[cell_name]] # modified copy of the parent that will be used to overwrite the cell in the pop
                                
                                return_list[['new_cells']] <- daughter_cells
                               return(return_list)
                               
                             })
  

  
  # split the results of this apply into updated parents and new cells:
  daughter_cells <- lapply(new_cell_list, function(cellname){
    cellname[['new_cells']]
  })
  
  updated_parents <- lapply(new_cell_list, function(cellname){
    cellname[['updated_parent']]
  })
  # give the updated parents names so that they can overwrite existing vals at these cell names
  names(updated_parents) <- cell_names_dividing_here
  
  
  
  # overwrite now-parents
  cell_population[cell_names_dividing_here] <- updated_parents[cell_names_dividing_here]
  

  # flatten the list of lists that was generated for daughter cells 
  daughter_cells <- unlist(daughter_cells, recursive = FALSE)
  
  # append the new daughter cells to the end of the growing cell pop
  for(daughter_name in names(daughter_cells)){
    cell_population[[daughter_name]] <- daughter_cells[[daughter_name]]
  }
  
  ##################################### DIE
  
  
  cells_alive_here_bool_list <- sapply(cell_population,
                                       function(cell) cell$terminal & cell$alive,
                                       USE.NAMES = TRUE)
  
  cell_names_alive_here <- names(cells_alive_here_bool_list)[cells_alive_here_bool_list]

  num_cells_alive_here <- length(cell_names_alive_here)
  # only certain cells will die at this timepoint, according to their respective cell type's death prob
  
  # only terminal (& alive) cells can die here...
  
  # death process via heteroplasmy scores:
  
  # generate a list with names == cellnames, values == named mutations and their frequency in that cell
  # we make use of the same func that is used to score allelic fractions in downstream tree reconstruction but exit early
  # recovered ints here is a map of cell name: all unique genome nums in mito to genome map
  # we do not condense to save time
  # profiles is a list of living cells for which we want heteroplasmy scores
  
  
  
  # generate heteroplasmy counts for all alive cells
  alive_mt_profiles <- lapply(cell_population[cell_names_alive_here], function(cell){
    cell$incoming_mt_profiles
  })
  
  # print('past alive_mt_profiles')
  recovered_ints_list <- lapply(cell_population[cell_names_alive_here], function(cell){
    
    prof <- cell$incoming_mt_profiles
    if(!is.null(prof) && length(dim(prof)) == 2 && nrow(prof) > 0){
      return(seq_len(nrow(prof)))
    } else{
      return(integer(0))
    }
   
  })
  
  # print('past recovered_ints_list')
  cell_to_num_mito_genomes_list <- lapply(recovered_ints_list, function(cell){
    return(length(cell))
  })
  
  if(consider_cell_heteroplasmy_scores){
    
    cell_to_heteroplasmy_counts <- new_create_one_score_mat(profiles = alive_mt_profiles,
                                                        recovered_ints = recovered_ints_list,
                                                        condense = FALSE, urid = unique_run_id, savename_prefix = 'not_used',
                                                        mt_or_bc = 'mt', binarize_score = FALSE, allelic_fraction_thresh = 0,
                                                        return_af_fracs = TRUE)
    
    if(!is.list(cell_to_heteroplasmy_counts)){
      hetero_death_draws <- rep(0, length(cell_names_alive_here))

    } else{
      norm_cell_heteroplasmy_scores <- get_norm_cell_heteroplasmy_scores(cell_mut_counts = cell_to_heteroplasmy_counts,
                                                                         cell_to_num_mito_genomes_list = cell_to_num_mito_genomes_list,
                                                                         heteroplasmy_severity_score_list = heteroplasmy_severity_score_list,
                                                                         positive_score_weight = positive_score_weight,
                                                                         hetero_sd = hetero_sd,
                                                                         cell_population = cell_population, 
                                                                         normalize_cell_mut_counts = FALSE)
      
      
     
      sapply(names(norm_cell_heteroplasmy_scores), function(name){
        cell_population[[name]]$heteroplasmy_score <<- norm_cell_heteroplasmy_scores[[name]]
      })
      
      hetero_death_probs <- vapply(norm_cell_heteroplasmy_scores[cell_names_alive_here], function(score){
        return(as.numeric(1- logistic_prob_survive_given_score(this_score = score, 
                                                               init_score = initial_heteroplasmy_score,
                                                               init_heteroplasmy_survive_prob = init_heteroplasmy_survive_prob,
                                                               beta = 1)))
      }, FUN.VALUE = numeric(1))
      
      hetero_death_draws <- rbinom(length(hetero_death_probs), 1, hetero_death_probs)
    }
    
  }
  
  # all cells undergo pruning by random selection and by zero mito
  
  cell_types_alive <- sapply(cell_population[cell_names_alive_here], `[[`, 'celltype')
  
  random_death_probs <- unname(unlist(cell_type_death_probs[cell_types_alive]))
  
  random_death_draws <- rbinom(length(cell_types_alive), 1, random_death_probs)
  
  zero_mito_draws <- sapply(cell_population[cell_names_alive_here], function(cell){
    return(length(cell$mito_to_genome_map) == 0)
  })
  
  # cell death considers heteroplasmy only if consider_cell_heteroplasmy_scores
  if(consider_cell_heteroplasmy_scores){
    death_occurs <- (hetero_death_draws == 1) | (random_death_draws == 1) | zero_mito_draws
  } else{
    death_occurs <- (random_death_draws == 1) | zero_mito_draws
  }
  
  dead_cells <- cell_names_alive_here[death_occurs]

  for(dead_cellname in dead_cells){
    cell_population[[dead_cellname]]$alive <- FALSE
    cell_population[[dead_cellname]]$death_time <- timepoint
  }

  cells_alive_here_bool_list <- sapply(cell_population,
                       function(cell) cell$terminal & cell$alive,
                       USE.NAMES = TRUE)
  cell_names_alive_here <- names(cells_alive_here_bool_list)[cells_alive_here_bool_list]
  

  # if all cells are dead at this timepoint, exit the entire simulation
  if(length(cell_names_alive_here) == 0){
    print(paste0('All cells have died at timepoint ', timepoint, '\nQuitting run without tree reconstruction'))
    quit(save = 'no')
  }
  
  # combining bc and mt processes into a single parlapply call:
  mutation_results <- parLapply(cl = one_cluster, X = cell_names_alive_here, 
                                fun = function(cell_name){
                                  
                                  this_cell <- cell_population[[cell_name]]
                                  this_cell_type <- cell_population[[cell_name]]$celltype
                                  editing_induced <- cell_population[[cell_name]]$induced_editing
                                  
                                  bc_result <- NULL
                                  
                                  if('bc' %in% poss_recon_modals){
                                    bc_result <- perform_all_bc_mutations(incoming_mut_mat = cell_population[[cell_name]]$incoming_bc_profiles,
                                                                    bg_transition_list = cell_type_basepos_bc_nontarget_transition_probs[[this_cell_type]][[editing_induced]],
                                                                    bg_transversion_list = cell_type_basepos_bc_nontarget_transversion_probs[[this_cell_type]][[editing_induced]],
                                                                    bg_insertion_list = cell_type_basepos_bc_nontarget_insertion_probs[[this_cell_type]][[editing_induced]],
                                                                    bg_deletion_list = cell_type_basepos_bc_nontarget_deletion_probs[[this_cell_type]][[editing_induced]],
                                                                    target_transition_list = cell_type_basepos_bc_target_transition_probs[[this_cell_type]][[editing_induced]],
                                                                    target_transversion_list = cell_type_basepos_bc_target_transversion_probs[[this_cell_type]][[editing_induced]],
                                                                    target_insertion_list = cell_type_basepos_bc_target_insertion_probs[[this_cell_type]][[editing_induced]],
                                                                    target_deletion_list = cell_type_basepos_bc_target_deletion_probs[[this_cell_type]][[editing_induced]],
                                                                    prob_sub_mat = cell_type_bc_sub_prob_mat[[this_cell_type]][[editing_induced]],
                                                                    timepoint_for_label = timepoint,
                                                                    urid = unique_run_id,
                                                                    interdel_dropout_radius = interdeletion_dropout_radius,
                                                                    interdel_dropout_prob = interdeletion_dropout_prob,
                                                                    prime_editing_system = prime_editing_system,
                                                                    ind_to_prime_seq_int_map = ind_to_prime_seq_int_map,
                                                                    force_target_transversions = forced_transversions,
                                                                    target_transversion_to_base = be_target_to_int,
                                                                    close_nuc_window_after_edit = close_nuc_window_after_edit,
                                                                    close_transition_window_after_edit = close_transition_window_after_edit,
                                                                    close_transversion_window_after_edit = close_transversion_window_after_edit,
                                                                    be_target_to_window_ind_list = be_target_to_window_ind_list,
                                                                    nuc_target_to_window_ind_list = nuc_target_to_window_ind_list,
                                                                    be_window_to_target_ind_list = be_window_to_target_ind_list,
                                                                    nuc_window_to_target_ind_list = nuc_window_to_target_ind_list)
                                  }
                                  
                                  
                                  
                                  mt_result <- NULL

                                  if('mt' %in% poss_recon_modals){
                                    mt_result <- perform_all_mt_mutations(incoming_mut_mat = cell_population[[cell_name]]$incoming_mt_profiles,
                                                                          bg_transition_list = cell_type_basepos_mt_nontarget_transition_probs[[this_cell_type]][[editing_induced]],
                                                                          bg_transversion_list = cell_type_basepos_mt_nontarget_transversion_probs[[this_cell_type]][[editing_induced]],
                                                                          bg_insertion_list = cell_type_basepos_mt_nontarget_insertion_probs[[this_cell_type]][[editing_induced]],
                                                                          bg_deletion_list = cell_type_basepos_mt_nontarget_deletion_probs[[this_cell_type]][[editing_induced]],
                                                                          prob_sub_mat = cell_type_mt_sub_prob_mat[[this_cell_type]][[editing_induced]])
                                  }
                                  
                                  return(list(cell_name = cell_name,
                                              bc_result = bc_result,
                                              mt_result = mt_result))
                                  
                                  
                                })
  
  for(res in mutation_results){
    cell_name <- res$cell_name
    
    if(!is.null(res$bc_result)){
      cell_population[[cell_name]]$incoming_bc_profiles <- res$bc_result
    }
    
    if(!is.null(res$mt_result)){
      cell_population[[cell_name]]$incoming_mt_profiles <- res$mt_result
    }
  }
  
  if(timepoint %in% sim_length_stopping_points){
    all_processes_at_stopping_point(timept_savename = paste0(custom_savename, '_time_', timepoint), 
                                    relative_timepoint = t, this_endpoint = timepoint)
    
    
  }
  
  return(cell_population)  
  
}

#' Build and write the ground-truth lineage tree
#'
#' Terminal cells become tips and non-terminal cells become internal nodes, with
#' one edge per parent/child pair. A population started from several independent
#' founders is a forest, so a synthetic root (renamed until it cannot collide
#' with a lineage string) is added above the founders to make a valid `phylo`;
#' that root stands for no division and carries no mutational branch.
#'
#' @param cell_population Named list of cell records, each carrying `parent` and
#'   `terminal` fields. It must contain at least one founder.
#' @param urid Unique run identifier, used as the output subdirectory name.
#' @param save_path_stem Filename stem; `.newick` is appended.
#' @param output_root Root output directory; defaults to `'output'`.
#' @return The `phylo` tree object that was written to disk.
#' @section Side effects: Creates and writes
#'   `<output_root>/processed_newicks/<urid>/<save_path_stem>.newick`. A
#'   single-cell population with no edges is written as a bare Newick tip
#'   instead of going through `ape::write.tree`.
create_ground_truth_tree <- function(cell_population,
                                     urid,
                                     save_path_stem,
                                     output_root = 'output'){
  all_lineage_strings <- names(cell_population)
  terminal_cell_booleans <- vapply(cell_population, function(cell){
    isTRUE(cell$terminal)
  }, logical(1))

  terminal_node_names <- all_lineage_strings[terminal_cell_booleans]
  internal_node_names <- all_lineage_strings[!terminal_cell_booleans]
  founder_node_names <- all_lineage_strings[vapply(cell_population, function(cell){
    is.null(cell$parent) || length(cell$parent) == 0
  }, logical(1))]
  if(length(founder_node_names) == 0){
    stop('Ground-truth population contains no founder cell.')
  }

  edge_rows <- lapply(all_lineage_strings, function(child){
    parent <- cell_population[[child]]$parent
    if(!is.null(parent) && length(parent) > 0){
      return(c(parent = as.character(parent), child = child))
    }
    NULL
  })
  edge_rows <- edge_rows[!vapply(edge_rows, is.null, logical(1))]
  edges <- if(length(edge_rows) > 0){
    do.call(rbind, edge_rows)
  } else{
    matrix(
      character(),
      nrow = 0,
      ncol = 2,
      dimnames = list(NULL, c('parent', 'child'))
    )
  }

  # Multiple independently initialized cells form a forest. A synthetic
  # time-zero root turns that forest into a valid phylo object without adding
  # a simulated division or mutational branch.
  if(length(founder_node_names) > 1){
    synthetic_root <- '__founder_root__'
    while(synthetic_root %in% all_lineage_strings){
      synthetic_root <- paste0('_', synthetic_root)
    }
    founder_edges <- cbind(
      parent = rep(synthetic_root, length(founder_node_names)),
      child = founder_node_names
    )
    edges <- rbind(founder_edges, edges)
    internal_node_names <- c(synthetic_root, internal_node_names)
  }

  all_node_names <- c(terminal_node_names, internal_node_names)
  node_map <- setNames(seq_along(all_node_names), all_node_names)
  edge_list_numeric <- matrix(
    as.integer(node_map[as.vector(t(edges))]),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c('parent', 'child'))
  )

  tree <- list(
    edge = edge_list_numeric,
    tip.label = terminal_node_names,
    Nnode = length(internal_node_names),
    node.label = internal_node_names
  )
  class(tree) <- 'phylo'

  tree_path <- file.path(
    output_root,
    'processed_newicks',
    urid,
    paste0(save_path_stem, '.newick')
  )
  if(!dir.exists(dirname(tree_path))){
    dir.create(dirname(tree_path), recursive = TRUE)
  }
  if(length(terminal_node_names) == 1 && nrow(edge_list_numeric) == 0){
    writeLines(paste0(terminal_node_names, ';'), tree_path)
  } else{
    ape::write.tree(tree, file = tree_path)
  }

  tree
}

#' Produce every output artefact for one simulation stopping point
#'
#' Saves the cell population as RDS and as a mutation-matrix-free JSON, writes
#' the ground-truth Newick tree, then, for each requested reconstruction method,
#' builds the downsampled profile lists and their FASTA or score-matrix
#' representations. The population and the run-level settings are read from the
#' enclosing scope rather than passed in.
#'
#' @param timept_savename Filename stem identifying this endpoint.
#' @param relative_timepoint Index of this timepoint within `poss_times`, used
#'   when slicing the timing vectors.
#' @param this_endpoint The simulation timepoint itself.
#' @param all_recon_methods Reconstruction methods to write output for; defaults
#'   to `input_args$reconstruction_method`. The recognised values are
#'   `'fasta_only'` and `'score'`.
#' @return `NULL`, invisibly; the function is called for its side effects.
#' @section Side effects: Writes under `output/cell_populations/`,
#'   `output/processed_newicks/`, `output/processed_lists/`,
#'   `output/processed_fastas/`, and `output/score_mats/`, all keyed by
#'   `unique_run_id`. `output/linstrings/<unique_run_id>/` is created but
#'   nothing is written into it.
all_processes_at_stopping_point <- function(timept_savename, relative_timepoint, this_endpoint, 
                                            all_recon_methods = as.character(input_args$reconstruction_method)
){
  
  #' Internal: record and plot the mutation timing for this endpoint
  #'
  #' @param poss_times Full simulation time grid.
  #' @param sim_time_vec_mt Elapsed mitochondrial mutation seconds per
  #'   timepoint.
  #' @param sim_time_vec_bc Elapsed barcode mutation seconds per timepoint.
  #' @param time_ind Number of leading timepoints to keep; defaults to the
  #'   enclosing `relative_timepoint`.
  #' @return The `ggsave` result; called for its side effects, which are a
  #'   timing CSV and a runtime scatter plot under
  #'   `output/timing_obj/<run_id>/`.
  #' @note Never called by the current stopping-point pass.
  describe_mutation_process_timing <- function(poss_times, sim_time_vec_mt, sim_time_vec_bc, time_ind = relative_timepoint){
    timing_df <- data.frame(cbind(poss_times[1:time_ind], sim_time_vec_mt, sim_time_vec_bc))
    colnames(timing_df) <- c('sim_timept', 'mt_mutation_time', 'bc_mutation_time')
    timing_df <- timing_df %>%
      mutate(tot_mutation_time = mt_mutation_time + bc_mutation_time)
    melted_timing_df <- reshape2::melt(timing_df, 
                                       measure.vars = c('mt_mutation_time', 'bc_mutation_time', 'tot_mutation_time'),
                                       variable.name = 'modality',
                                       value.name = 'seconds')
    
    # create timing_obj subdirectory if it doesn't exist
    if(!dir.exists(file.path('output', 'timing_obj', unique_run_id))){
      dir.create(file.path('output', 'timing_obj', unique_run_id), recursive = TRUE)
    }
    
    write.csv(timing_df, file.path('output', 'timing_obj', unique_run_id, paste0('mt_bc_sim_time_', timept_savename, 'NUMCORES', 
                                                                                 input_args$num_cores, '_', unique_run_id, '.csv')))  
    
    runtime_plot <- ggplot(melted_timing_df, aes(x = sim_timept, y = seconds, color = modality)) +
      geom_point() +
      theme_classic() + 
      labs(title = 'Simulation runtime', 
           x = 'Simulation timepoint',
           y = 'Elapsed seconds at timepoint')
    
    
    # create plots directory and runtime_plots subdirectory if they don't exist:
    if(!dir.exists(file.path('output', 'plots', 'runtime_plots', unique_run_id))){
      dir.create(file.path('output', 'plots', 'runtime_plots', unique_run_id), recursive = TRUE)
    }
    
    ggsave(plot = runtime_plot, filename = file.path('output', 'timing_obj', unique_run_id, 
                                                     paste0(timept_savename,
                                                            '_', unique_run_id,'.png')
    ))
    
  }
  
  
  #' Internal: write the endpoint's mutation profile lists to disk
  #'
  #' @param mt_profiles List of mitochondrial mutation matrices.
  #' @param bc_profiles List of barcode mutation matrices.
  #' @return `NULL`; writes `simresults_mt_profiles_*` and
  #'   `simresults_bc_profiles_*` RDS files under
  #'   `output/mut_profiles/<run_id>/`, one per modality present in
  #'   `poss_recon_modals`.
  #' @note Never called by the current stopping-point pass.
  save_mutation_profiles <- function(mt_profiles, bc_profiles){
    # create mut_profiles subdirectory if it doesn't exist
    if(!dir.exists(file.path('output', 'mut_profiles', unique_run_id))){
      dir.create(file.path('output', 'mut_profiles', unique_run_id), recursive = TRUE)
    }
    
    if('mt' %in% poss_recon_modals){
      saveRDS(mt_profiles, file.path('output', 'mut_profiles', unique_run_id, 
                                     paste0('simresults_mt_profiles_',
                                            timept_savename,
                                            '_', unique_run_id, '.rds')))  
    }
    
    if('bc' %in% poss_recon_modals){
      saveRDS(bc_profiles, file.path('output', 'mut_profiles', unique_run_id,
                                     paste0('simresults_bc_profiles_',
                                            timept_savename,
                                            '_', unique_run_id, '.rds')))  
    }
    
    
  }
  
  
  #' Internal: build dot-separated lineage strings from a parent index vector
  #'
  #' @param cell_lineage Integer vector whose i-th element is the position of
  #'   cell i's parent, or `0` when cell i is a founder.
  #' @return A character vector of lineage strings, each daughter appending `.1`
  #'   or `.2` to its parent's string.
  #' @note Assigns the result to the global `lineage_strings` as well. Never
  #'   called by the current stopping-point pass, which names cells by the
  #'   `linstring` already stored on each record.
  create_lineage_strings <- function(cell_lineage){
    
    edge_from <- integer(length = length(cell_lineage))
    edge_to <- integer(length = length(cell_lineage))
    
    lineage_strings <<- character(length = length(cell_lineage))
    for(i in seq_len(length(cell_lineage))){
      if(cell_lineage[i] == 0){ # if the cell has no parent, it's a founder cell
        lineage_strings[i] <- i
        
      }
      
      else{ # if the cell has a parent
        
        edge_from[i] <- as.integer(cell_lineage[i]) # relative position of parent cell
        edge_to[i] <- i # relative position of the daughter cell
        
        temp_traceback <- cell_lineage[i] # look at the position of the parent cell in the founder parents list
        
        num_occur <- length(which(lineage_strings[1:i] == paste0(lineage_strings[temp_traceback], '.1')))
        
        if(num_occur == 0){ # if this is the first daughter cell of cell_lineage[i]:
          lineage_strings[i] <- paste0(lineage_strings[temp_traceback], '.1')
        }
        else if(num_occur == 1){ # if this is the second daughter cell. because each cell will now split into two daughters
          lineage_strings[i] <- paste0(lineage_strings[temp_traceback], '.2')
        }
        
      }
      
      
    } 
    return(lineage_strings)
  }
  
  linstring_dir_path <- file.path('output', 'linstrings', unique_run_id)
  if(!dir.exists(linstring_dir_path)){
    dir.create(linstring_dir_path, recursive = TRUE)
  }
  
  
  if(!dir.exists(file.path('output', 'cell_populations', unique_run_id))){
    dir.create(file.path('output', 'cell_populations', unique_run_id), recursive = TRUE)
  }
  saveRDS(cell_population, file.path('output', 'cell_populations', unique_run_id, paste0('cell_population_', timept_savename, '.rds')))
  
  # also save a copy of the cell population as a json, which can be used for downstream analyses in python
  no_mutmats_cellpop <- lapply(cell_population, function(cell) cell[!names(cell) %in% c('incoming_mt_profiles', 'incoming_bc_profiles')])
  pop_json <- toJSON(no_mutmats_cellpop)
  json_save_path <- file.path('output', 'cell_populations', unique_run_id, paste0('cell_population_', timept_savename, '.json'))
  write(pop_json, file = json_save_path)
  
  
  # create processed_newicks dir if it doesn't already exist
  if(!dir.exists(file.path('output', 'processed_newicks', unique_run_id))){
    dir.create(file.path('output', 'processed_newicks', unique_run_id), recursive = TRUE)
  }
  
  print('Writing ground truth tree for this timepoint ... ')
  
  this_timepoint_ground_truth_tree <- create_ground_truth_tree(cell_population = cell_population, 
                                                               save_path_stem = paste0('ground_truth_tree_', timept_savename),
                                                               urid = unique_run_id)
  
  
  
  
  #' Internal: build and write the downsampled profile lists for one endpoint
  #'
  #' Handles the `'all_cells'` and `'terminal'` FASTA types. Terminal output
  #' iterates every combination of per-cell-type sampling fractions, barcode
  #' integration counts, and recovery probabilities, saving each profile list as
  #' RDS and then writing it either as a FASTA or as a score matrix. With
  #' `writeout_type = 'score'` and `combine_mt_bc`, matching mitochondrial and
  #' barcode score matrices (same collapse-deletions setting, same endpoint) are
  #' also column-bound into joint matrices and FASTAs.
  #'
  #' @param cell_population Named list of cell records for this endpoint.
  #' @param bc_integrations Integer barcode integration counts to iterate over.
  #' @param mito_recovery_probs Mitochondrial genome recovery probabilities.
  #' @param bc_recovery_probs Barcode integration recovery probabilities.
  #' @param poss_fasta_types Which cell sets to write: `'all_cells'` and/or
  #'   `'terminal'`.
  #' @param bc_umis Per-integration barcode UMI sequences.
  #' @param writeout_type `'fasta_only'` or `'score'`.
  #' @param poss_recon_modals Modalities to write (`'mt'`, `'bc'`).
  #' @param this_timept_savename Filename stem; defaults to the enclosing
  #'   `timept_savename`.
  #' @return `NULL`; writes RDS profile lists under
  #'   `output/processed_lists/<run_id>/`, FASTAs under
  #'   `output/processed_fastas/<run_id>/`, and score matrices under
  #'   `output/score_mats/<run_id>/`. A sampling combination that recovers no
  #'   terminal cell raises a warning and is skipped.
  #' @note Score output for `'all_cells'` is not implemented; that branch only
  #'   prints a message.
  create_modified_profile_lists <- function(cell_population,
                                            bc_integrations,
                                            mito_recovery_probs,
                                            bc_recovery_probs,
                                            poss_fasta_types,
                                            bc_umis,
                                            writeout_type,
                                            poss_recon_modals,
                                            this_timept_savename = timept_savename){
    
    
    # we have to iterate through number of integrations as well as recovery probs
    
    if(!dir.exists(file.path('output', 'processed_lists', unique_run_id))){
      dir.create(file.path('output', 'processed_lists', unique_run_id), recursive = TRUE)
    }
    
    if('all_cells' %in% poss_fasta_types){
      
      if('bc' %in% poss_recon_modals){
        # write all bc profiles to fasta (and save mutational profiles to list)
        all_bc_profiles <- lapply(cell_population, function(cell){
          cell$incoming_bc_profiles
        })
        
        bc_list_assign_name <- file.path('output', 'processed_lists', unique_run_id, paste0('bc_all_cells_', this_timept_savename, '.rds'))
        saveRDS(all_bc_profiles, bc_list_assign_name)
        
        if(writeout_type == 'fasta_only'){
          fasta_dir_path <- file.path('output', 'processed_fastas', unique_run_id)
          if(!dir.exists(fasta_dir_path)){
            dir.create(fasta_dir_path, recursive = TRUE)
          }
          
          fasta_savename <- file.path(fasta_dir_path, paste0('bc_all_cells_', this_timept_savename, '.fasta'))
          
          write_all_cell_sequences(cell_mutmats = all_bc_profiles, 
                                   reference = baseline_seq_nucs_bc, 
                                   bc_integration_umis = rep(list(bc_umis), length(all_bc_profiles)), # need list of all ints for each cell
                                   output_fasta_name = fasta_savename,
                                   fasta_type = 'ALL_CELLS')
        } else if(writeout_type == 'score'){
          
          print('scores for all cells not yet implemented; trees generally not built from all cells (incl internal)')
          
        }
        
      }
      
      
      
      if('mt' %in% poss_recon_modals){
        # write all mt profiles to fasta (and save mutational profiles to list)
        all_mt_profiles <- lapply(cell_population, function(cell){
          cell$incoming_mt_profiles
        })
         
        mt_list_assign_name <- file.path('output', 'processed_lists', unique_run_id, paste0('mt_all_cells_', this_timept_savename, '.rds'))
        saveRDS(all_mt_profiles, mt_list_assign_name)
        
        
        if(writeout_type == 'fasta_only'){
          fasta_dir_path <- file.path('output', 'processed_fastas', unique_run_id)
          if(!dir.exists(fasta_dir_path)){
            dir.create(fasta_dir_path, recursive = TRUE)
          }
          
          fasta_savename <- file.path(fasta_dir_path, paste0('mt_all_cells_', this_timept_savename, '.fasta'))
          
          write_all_cell_sequences(cell_mutmats = all_mt_profiles, 
                                   reference = baseline_seq_nucs_mt,  
                                   output_fasta_name = fasta_savename,
                                   fasta_type = 'ALL_CELLS')
          
          
          
          
        } else if(writeout_type == 'score'){
          
          print('scores for all cells not yet implemented; trees generally not built from all cells (incl internal)')
          
        }
        
      }
    }
    
    
    if('terminal' %in% poss_fasta_types){
      
      
      print('Creating downsampled profile lists ... ')
      
      # get all possible combinations of cell sampling fracs across cell types
      all_sampling_fracs <- expand.grid(cell_type_poss_sampling_fracs)
      colnames(all_sampling_fracs) <- names(cell_type_poss_sampling_fracs)
      
      # iterate through these cell sampling frac combos
      for(i in seq_len(nrow(all_sampling_fracs))){
        
        # generate a name that includes the sampling rate for each cell type  
        # cell type will be separated from its sampling frac by -
        # cell types will be separated from one another by _
        this_sampling_name <- paste(paste(colnames(all_sampling_fracs), 
                                          as.numeric(unlist(all_sampling_fracs[i, ])), sep = '-'), collapse = '_')
        
        # make the downsampled terminal cell population, according to terminal, alive, and 
        # cell type-specific sampling frac for this iteration
        terminal_cell_population_inds <- lapply(cell_population, function(cell){
          if((cell$terminal == FALSE) | (cell$alive == FALSE)){
            return(FALSE)
          }
          this_cell_type <- cell$celltype
          
          # use this particular iteration's combo of cell type recovery probs
          this_cell_recovery_prob <- as.numeric(all_sampling_fracs[[this_cell_type]][i])
          this_cell_recovered <- rbinom(n = 1, size = 1, prob = this_cell_recovery_prob)
          
          if(this_cell_recovered){
            return(TRUE)
          }
          return(FALSE)
        })
        
        
        terminal_cell_population_names <- names(cell_population)[which(terminal_cell_population_inds == TRUE)]
        terminal_cell_population <- cell_population[terminal_cell_population_names]

        if(length(terminal_cell_population) == 0){
          warning(sprintf(
            'No terminal cells were sampled for combination %s at %s; skipping.',
            this_sampling_name,
            this_timept_savename
          ))
          next
        }
        
        if('bc' %in% poss_recon_modals){
          
          print('Starting bc recon modal workflow ... ')
          for(num_bc_ints in bc_integrations){
            
            for(bc_int_recovery_prob in bc_recovery_probs){
              
              bc_profiles_ints_and_umis <- get_profiles_ints_and_umis(cell_pop = terminal_cell_population,
                                                                      num_ints = num_bc_ints,
                                                                      int_rec_prob = bc_int_recovery_prob,
                                                                      bc_or_mt = 'bc',
                                                                      umis = bc_umis)
              
              bc_subsetted_profiles <- lapply(bc_profiles_ints_and_umis, function(cell){
                cell[['mut_mat']]
              })
              names(bc_subsetted_profiles) <- terminal_cell_population_names
              
              bc_recovered_ints <- lapply(bc_profiles_ints_and_umis, function(cell){
                cell[['which_ints_recovered']]
              })
              
              bc_recovered_umis <- lapply(bc_profiles_ints_and_umis, function(cell){
                cell[['recovered_umis']]
              }) 
              
              bc_combo_name <- paste0('proc_bc_list_', num_bc_ints, 
                                      '_ints_RP_', bc_int_recovery_prob, 
                                      '_samp_', this_sampling_name, '_',
                                      this_timept_savename)
              bc_list_assign_name <- file.path('output', 'processed_lists', unique_run_id, paste0(bc_combo_name, '.rds'))
              
              saveRDS(bc_subsetted_profiles, bc_list_assign_name)
              
              if(writeout_type == 'fasta_only'){
                
                # immediately write the fasta (IN THE FOR LOOP)
                fasta_dir_path <- file.path('output', 'processed_fastas', unique_run_id)
                if(!dir.exists(fasta_dir_path)){
                  dir.create(fasta_dir_path, recursive = TRUE)
                }
                bc_fasta_savename <- file.path(fasta_dir_path, paste0(bc_combo_name, '.fasta'))
                
                write_all_cell_sequences(cell_mutmats = bc_subsetted_profiles, 
                                         reference = baseline_seq_nucs_bc, 
                                         bc_integration_umis = bc_recovered_umis,
                                         output_fasta_name = bc_fasta_savename,
                                         fasta_type = 'TERM')
              } else if(writeout_type == 'score'){
                
                
                for(collapse in scoremat_collapse_deletions){
                  
                  new_create_one_score_mat(profiles = bc_subsetted_profiles, 
                                       recovered_ints = bc_recovered_ints, 
                                       condense = collapse, 
                                       urid = unique_run_id, 
                                       savename_prefix = paste0(bc_combo_name, '_CD_', substr(collapse, 1, 1)),
                                       mt_or_bc = 'bc')
                  
                  
                  
                } 
              }
              
            }
          }
        }
        
        if('mt' %in% poss_recon_modals){
          
          print('Starting mt recon modal workflow ... ')
          
          for(mito_recovery_prob in mito_recovery_probs){
            
            mt_combo_name <- paste0('proc_mt_list_', init_num_mito_genomes, 
                                    '_ints_RP_', mito_recovery_prob, 
                                    '_samp_', this_sampling_name, '_',
                                    this_timept_savename)
            
            mt_list_assign_name <- file.path('output', 'processed_lists', unique_run_id, paste0(mt_combo_name, '.rds'))
            
            gpiu_start <- Sys.time()
            
            gpiu_time <- system.time(
              mt_profiles_ints <- get_profiles_ints_and_umis(cell_pop = terminal_cell_population,
                                                             num_ints = NULL,
                                                             int_rec_prob = mito_recovery_prob,
                                                             bc_or_mt = 'mt',
                                                             umis = NULL)
            )
            
            
            gpiu_end <- Sys.time()
            
            mt_subsetted_profiles <- lapply(mt_profiles_ints, function(cell){
              cell[['mut_mat']]
            })
            names(mt_subsetted_profiles) <- terminal_cell_population_names
            
            mt_recovered_ints <- lapply(mt_profiles_ints, function(cell){
              cell[['which_ints_recovered']]
            })
            
            saveRDS(mt_subsetted_profiles, mt_list_assign_name)
            
            if(writeout_type == 'fasta_only'){
              
              fasta_dir_path <- file.path('output', 'processed_fastas', unique_run_id)
              if(!dir.exists(fasta_dir_path)){
                dir.create(fasta_dir_path, recursive = TRUE)
              }
              mt_fasta_savename <- file.path(fasta_dir_path, paste0(mt_combo_name, '.fasta'))
              
              # immediately write the fasta (IN THE FOR LOOP)
              write_all_cell_sequences(cell_mutmats = mt_subsetted_profiles, 
                                       reference = baseline_seq_nucs_mt, 
                                       output_fasta_name = mt_fasta_savename,
                                       fasta_type = 'TERM')
            } else if(writeout_type == 'score'){
              
              for(collapse in scoremat_collapse_deletions){
                
                new_cosm_start <- Sys.time()
                new_cosm_time <- system.time(
                  new_create_one_score_mat(profiles = mt_subsetted_profiles, 
                                           recovered_ints = mt_recovered_ints, 
                                           condense = collapse, 
                                           urid = unique_run_id, 
                                           savename_prefix = paste0(mt_combo_name, '_CD_', substr(collapse, 1, 1)),
                                           mt_or_bc = 'mt',
                                           binarize_score = binarize_mutation_scores,
                                           allelic_fraction_thresh = mt_allelic_fraction_thresholds)
                )
                new_cosm_end <- Sys.time()
                
              }
              
              
            }
            
            
          }
          
          
          # }
        }
      }
      
      # if we're working with score matrices (binary here) and want to combine mt and bc signals into one mat
      if(writeout_type == 'score'){
        if(combine_mt_bc){
          # create pairwise score mats between 1 mt mat and 1 bc mat
          # files must have respective bc or mt label at the beginning
          # and must include the timepoint savename so that we're only comparing mats at this endpoint
          # also only want to pair CD TRUE/TRUE or CD FALSE/FALSE
          
          mt_score_mat_paths <- list.files(path = file.path('output', 'score_mats', unique_run_id, 'matrices'),
                                           pattern = '^proc_mt_',
                                           full.names = TRUE)
          mt_score_mat_paths <- mt_score_mat_paths[grepl(timept_savename, mt_score_mat_paths)]
          
          
          
          bc_score_mat_paths <- list.files(path = file.path('output', 'score_mats', unique_run_id, 'matrices'),
                                           pattern = '^proc_bc_',
                                           full.names = TRUE)
          bc_score_mat_paths <- bc_score_mat_paths[grepl(timept_savename, bc_score_mat_paths)]
          
          
          
          for(mt_score_mat_path in mt_score_mat_paths){
            mt_score_mat <- readRDS(mt_score_mat_path)
            
            # check if deletions are condensed or not (only want to pair mt & bc with same condense status)
            mt_cd <- str_extract(mt_score_mat_path, '(?<=CD_)\\w{1}(?=_)')
            
            
            for(bc_score_mat_path in bc_score_mat_paths){
              bc_cd <- str_extract(bc_score_mat_path, '(?<=CD_)\\w{1}(?=\\.)')
              
              # if both these paths have the same condensed-deletion logic, join
              if(mt_cd == bc_cd){
                trimmed_mt_name <- sub('.*\\/(.*).rds', '\\1', mt_score_mat_path)
                trimmed_bc_name <- sub('.*\\/(.*).rds', '\\1', bc_score_mat_path)
                joint_save_name <- paste0('J_', trimmed_mt_name, '_',
                                          trimmed_bc_name)
                
                bc_score_mat <- readRDS(bc_score_mat_path)
                combined_score_mat <- cbind(mt_score_mat, bc_score_mat)
                saveRDS(combined_score_mat, file.path('output', 'score_mats', unique_run_id, 'matrices', 
                                                      paste0(joint_save_name, '.rds')))
                new_scoremat_to_fasta(scoremat = combined_score_mat, 
                                      output_fasta_path = file.path('output', 'score_mats', unique_run_id, 'phylips', paste0(joint_save_name, '.fasta')))
              }  
            }
            
          }
          
          
          
        }  
      }
      
      else if(writeout_type == 'fasta_only'){
        #' Internal: write the reference FASTA for one modality
        #'
        #' @param bc_or_mt Either `'bc'` or `'mt'`; selects the filename prefix.
        #' @param reference_seq Reference nucleotide vector for a single
        #'   integration or genome.
        #' @param max_number_of_integrations Number of times the reference is
        #'   repeated so it matches the concatenated read length.
        #' @param run_id Unique run identifier naming the output subdirectory.
        #' @return `NULL`; writes `<bc_or_mt>_<n>_ints_reference.fasta` under
        #'   `output/processed_fastas/<run_id>/reference_seqs/`.
        write_reference_fastas <- function(bc_or_mt,
                                           reference_seq,
                                           max_number_of_integrations,
                                           run_id){
          # regardless of whether reference seq is supplied earlier or generated here, must make compatible with number of integrations 
          reference_seq <- rep(reference_seq, max_number_of_integrations)
          
          if(!dir.exists(file.path('output', 'processed_fastas', run_id, 'reference_seqs'))){
            dir.create(file.path('output', 'processed_fastas', run_id, 'reference_seqs'), recursive = TRUE)
          }
          if(bc_or_mt == 'bc'){
            write.fasta(reference_seq, names = c('REFERENCE'), file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', 
                                                                                    paste0(bc_or_mt, '_', max_number_of_integrations, '_ints_reference.fasta')))
          } else if(bc_or_mt == 'mt'){
            write.fasta(reference_seq, names = c('REFERENCE'), file.out = file.path('output', 'processed_fastas', run_id, 'reference_seqs', 
                                                                                    paste0(bc_or_mt, '_', max_number_of_integrations, '_ints_reference.fasta')))
          }  
        }
        
        # write the reference seqs for each provided number of max ints/genomes:
        if('bc' %in% poss_recon_modals){
          for(max_num_bc_ints in poss_num_bc_integrations){
            write_reference_fastas(bc_or_mt = 'bc',
                                   reference_seq = baseline_seq_nucs_bc,
                                   max_number_of_integrations = max_num_bc_ints,
                                   run_id = unique_run_id)  
          }  
        }
        
        if('mt' %in% poss_recon_modals){
          
          write_reference_fastas(bc_or_mt = 'mt',
                                 reference_seq = baseline_seq_nucs_mt,
                                 max_number_of_integrations = init_num_mito_genomes,
                                 run_id = unique_run_id)
          
        }
        
        
        
      }
      
      
    }
  }
  
  for(recon_method in all_recon_methods){
    create_mod_profile_lists_start_time <- Sys.time()
    cmpl_time <- system.time(
      create_modified_profile_lists(cell_population = cell_population,
                                    bc_integrations = poss_num_bc_integrations,
                                    mito_recovery_probs = poss_mt_genome_recovery_probs,
                                    bc_recovery_probs = poss_bc_integration_recovery_probs,
                                    poss_fasta_types = poss_fasta_types,
                                    bc_umis = bc_int_umis,
                                    writeout_type = recon_method,
                                    poss_recon_modals = poss_recon_modals)  
    )
    create_mod_profile_lists_end_time <- Sys.time()

  }
  
  
  
  
  
  
  
}


# ---- Assemble the simulation arguments and initialize the run ----
sim_arglist <- list(num_clusters = input_args$num_cores, 
                    init_pop_size = input_args$num_init_cells,
                    init_incoming_mt_profile = init_incoming_mt_profile,
                    sim_length = max(sim_length_stopping_points),
                    cell_type_cell_cycle_length = cell_type_cell_cycle_length,
                    num_rows_mt = init_num_mito_genomes,
                    num_cols_mt = input_args$mito_genome_length,
                    num_rows_bc = max(poss_num_bc_integrations),
                    num_cols_bc = input_args$bc_length,
                    time_inc = time_inc,
                    cell_type_basepos_bc_nontarget_transition_probs = cell_type_basepos_bc_nontarget_transition_probs,
                    cell_type_basepos_bc_nontarget_transversion_probs = cell_type_basepos_bc_nontarget_transversion_probs,
                    cell_type_basepos_bc_nontarget_insertion_probs = cell_type_basepos_bc_nontarget_insertion_probs,
                    cell_type_basepos_bc_nontarget_deletion_probs = cell_type_basepos_bc_nontarget_deletion_probs,
                    cell_type_basepos_mt_nontarget_transition_probs = cell_type_basepos_mt_nontarget_transition_probs,
                    cell_type_basepos_mt_nontarget_transversion_probs = cell_type_basepos_mt_nontarget_transversion_probs,
                    cell_type_basepos_mt_nontarget_insertion_probs = cell_type_basepos_mt_nontarget_insertion_probs,
                    cell_type_basepos_mt_nontarget_deletion_probs = cell_type_basepos_mt_nontarget_deletion_probs,
                    cell_type_basepos_bc_target_transition_probs = cell_type_basepos_bc_target_transition_probs,
                    cell_type_basepos_bc_target_transversion_probs = cell_type_basepos_bc_target_transversion_probs,
                    cell_type_basepos_bc_target_insertion_probs = cell_type_basepos_bc_target_insertion_probs,
                    cell_type_basepos_bc_target_deletion_probs = cell_type_basepos_bc_target_deletion_probs,
                    cell_type_mt_sub_prob_mat = cell_type_mt_sub_prob_mat,
                    cell_type_bc_sub_prob_mat = cell_type_bc_sub_prob_mat,
                    cell_type_death_probs = cell_type_death_probs,
                    uninduced_tm_list = uninduced_tm_list,
                    induced_tm_list = induced_tm_list,
                    differentiation_induction_timepoint = input_args$differentiation_induction$timepoint,
                    differentiation_induction_num_cells = input_args$differentiation_induction$num_cells,
                    differentiation_induction_frac_cells = input_args$differentiation_induction$frac_cells,
                    editing_induction_timepoint = input_args$editing_induction$timepoint,
                    editing_induction_num_cells = input_args$editing_induction$num_cells,
                    editing_induction_frac_cells = input_args$editing_induction$frac_cells,
                    forced_transversions = force_transversions,
                    be_target_to_int = be_target_to_int,
                    custom_savename = custom_savename, 
                    sim_length_stopping_points = sim_length_stopping_points,
                    poss_fasta_types = poss_fasta_types,
                    include_var_pos_fasta = include_var_pos_fasta,
                    founder_cell_type = founder_cell_type,
                    interdeletion_dropout_radius = input_args$nuclease_targets$interdeletion_dropout_radius,
                    interdeletion_dropout_prob = input_args$nuclease_targets$interdeletion_dropout_prob,
                    poss_recon_modals = poss_recon_modals,
                    mito_inheritance_pattern = mito_inheritance_pattern,
                    mito_to_genome_map = mito_to_genome_map,
                    fusion_events_per_mito_per_division = fusion_events_per_mito_per_division,
                    split_events_per_mito_per_division = split_events_per_mito_per_division,
                    post_mitotic_mt_deletion_frac = post_mitotic_mt_deletion_frac,
                    heteroplasmy_severity_score_list = heteroplasmy_severity_score_list,
                    heteroplasmy_variant_fractions = heteroplasmy_variant_fractions,
                    init_heteroplasmy_survive_prob = init_heteroplasmy_survive_prob,
                    initial_heteroplasmy_score = initial_heteroplasmy_score,
                    hetero_sd = hetero_sd,
                    positive_score_weight = positive_score_weight,
                    init_num_mito_genomes = init_num_mito_genomes,
                    ind_to_prime_seq_int_map = ind_to_prime_seq_int_map,
                    ind_to_prime_seq_nuc_map = ind_to_prime_seq_nuc_map,
                    prime_editing_system = prime_editing_system,
                    close_nuc_window_after_edit = close_nuc_window_after_edit,
                    close_transition_window_after_edit = close_transition_window_after_edit,
                    close_transversion_window_after_edit = close_transversion_window_after_edit,
                    be_target_to_window_ind_list = be_target_to_window_ind_list,
                    nuc_target_to_window_ind_list = nuc_target_to_window_ind_list,
                    be_window_to_target_ind_list = be_window_to_target_ind_list,
                    nuc_window_to_target_ind_list = nuc_window_to_target_ind_list,
                    consider_cell_heteroplasmy_scores = consider_cell_heteroplasmy_scores
                    
                    )
cell_population <- do.call(setup_sim, sim_arglist)

for(i in seq_along(sim_arglist)){
  assign(names(sim_arglist)[i], sim_arglist[[i]], envir = .GlobalEnv)
}

# ---- Main simulation loop ----
# actually run the simulation
for(t in seq_along(poss_times)){
  
  # has to be something like: for(t in 1:length(which.min(sim_lengths_with_breakpoints)))
  # could also make this into a while loop 
  
  if(poss_times[t] > 0){
    
    print(paste0('Now simulating timepoint t = ', poss_times[t], ' ... '))
    
    mcf_start <- Sys.time()
    
    mcf_overall <- system.time(
      
      cell_population <- multi_core_func(timepoint = poss_times[t], 
                                         sim_length = sim_length,
                                          cell_population = cell_population,
                                          cell_type_basepos_bc_nontarget_transition_probs = cell_type_basepos_bc_nontarget_transition_probs, 
                                          cell_type_basepos_bc_nontarget_transversion_probs = cell_type_basepos_bc_nontarget_transversion_probs, 
                                          cell_type_basepos_bc_nontarget_insertion_probs = cell_type_basepos_bc_nontarget_insertion_probs, 
                                          cell_type_basepos_bc_nontarget_deletion_probs = cell_type_basepos_bc_nontarget_deletion_probs,
                                          cell_type_basepos_mt_nontarget_transition_probs = cell_type_basepos_mt_nontarget_transition_probs, 
                                          cell_type_basepos_mt_nontarget_transversion_probs = cell_type_basepos_mt_nontarget_transversion_probs, 
                                          cell_type_basepos_mt_nontarget_insertion_probs = cell_type_basepos_mt_nontarget_insertion_probs, 
                                          cell_type_basepos_mt_nontarget_deletion_probs = cell_type_basepos_mt_nontarget_deletion_probs, 
                                          cell_type_basepos_bc_target_transition_probs = cell_type_basepos_bc_target_transition_probs, 
                                          cell_type_basepos_bc_target_transversion_probs = cell_type_basepos_bc_target_transversion_probs, 
                                          cell_type_basepos_bc_target_insertion_probs = cell_type_basepos_bc_target_insertion_probs, 
                                          cell_type_basepos_bc_target_deletion_probs = cell_type_basepos_bc_target_deletion_probs,
                                          cell_type_death_probs = cell_type_death_probs,
                                          cell_type_mt_sub_prob_mat = cell_type_mt_sub_prob_mat, 
                                          cell_type_bc_sub_prob_mat = cell_type_bc_sub_prob_mat, 
                                          uninduced_tm_list = uninduced_tm_list,
                                          induced_tm_list = induced_tm_list,
                                         differentiation_induction_timepoint = differentiation_induction_timepoint,
                                         editing_induction_timepoint = editing_induction_timepoint,
                                         differentiation_induction_num_cells = differentiation_induction_num_cells,
                                         editing_induction_num_cells = editing_induction_num_cells,
                                         differentiation_induction_frac_cells = differentiation_induction_frac_cells,
                                         editing_induction_frac_cells = editing_induction_frac_cells,
                                         already_assigned_editing_induction = already_assigned_editing_induction,
                                         already_assigned_diff_induction = already_assigned_diff_induction,
                                         forced_transversions = forced_transversions,
                                         be_target_to_int = be_target_to_int,
                                         unique_run_id = unique_run_id,
                                         interdeletion_dropout_prob = interdeletion_dropout_prob,
                                         interdeletion_dropout_radius = interdeletion_dropout_radius,
                                         scoremat_collapse_deletions = scoremat_collapse_deletions,
                                         combine_mt_bc = combine_mt_bc,
                                         binarize_mutation_scores = binarize_mutation_scores,
                                         mt_allelic_fraction_thresholds = mt_allelic_fraction_thresholds,
                                         poss_recon_modals = poss_recon_modals,
                                         heteroplasmy_severity_score_list = heteroplasmy_severity_score_list,
                                         fusion_events_per_mito_per_division = fusion_events_per_mito_per_division,
                                         split_events_per_mito_per_division = split_events_per_mito_per_division,
                                         post_mitotic_mt_deletion_frac = post_mitotic_mt_deletion_frac,
                                         heteroplasmy_variant_fractions = heteroplasmy_variant_fractions,
                                         positive_score_weight = positive_score_weight,
                                         hetero_sd = hetero_sd,
                                         mito_inheritance_pattern = mito_inheritance_pattern,
                                         init_num_mito_genomes = init_num_mito_genomes,
                                         prime_editing_system = prime_editing_system,
                                         ind_to_prime_seq_int_map = ind_to_prime_seq_int_map,
                                         ind_to_prime_seq_nuc_map = ind_to_prime_seq_nuc_map,
                                         close_nuc_window_after_edit = close_nuc_window_after_edit,
                                         close_transition_window_after_edit = close_transition_window_after_edit,
                                         close_transversion_window_after_edit = close_transversion_window_after_edit,
                                         be_target_to_window_ind_list = be_target_to_window_ind_list,
                                         nuc_target_to_window_ind_list = nuc_target_to_window_ind_list,
                                         be_window_to_target_ind_list = be_window_to_target_ind_list,
                                         nuc_window_to_target_ind_list = nuc_window_to_target_ind_list,
                                         consider_cell_heteroplasmy_scores = consider_cell_heteroplasmy_scores,
                                         poss_times = poss_times,
                                         custom_savename = custom_savename,
                                         sim_length_stopping_points = sim_length_stopping_points,
                                         t = t
  
                                         )
    #   }
    )
    
    
    mcf_end <- Sys.time()
    
  }
  
  if(poss_times[t] %in% sim_length_stopping_points){
    if(t == length(poss_times)){
      stopCluster(one_cluster)
    }
  }
  

  
  
}


# join together specs/results across all endpoints
#' Merge the per-endpoint result CSVs into a single file
#'
#' @param unique_run_id Run identifier naming the merged-results directory;
#'   defaults to the global of the same name.
#' @return `NULL`; reads every CSV in
#'   `output/results/merged_results/<unique_run_id>/`, row-binds them, and
#'   writes `merged_results_specs_<unique_run_id>_all_endpoints.csv` back into
#'   that same directory.
#' @note Defined here for downstream use but not called by this script.
join_endpoint_results <- function(unique_run_id = unique_run_id){
  
  merged_results_dir_path <- file.path('output', 'results', 'merged_results', unique_run_id)
  
  merged_filenames <- list.files(merged_results_dir_path, full.names = TRUE)
  list_of_merged_tables <- lapply(merged_filenames, read.csv, row.names = 1)
  all_merged_results <- do.call(rbind, list_of_merged_tables)
  
  write.csv(all_merged_results, file.path(merged_results_dir_path, 
                                          paste0('merged_results_specs_',
                                                 unique_run_id, '_all_endpoints.csv')))
  
}




#' Plot Robinson-Foulds distance against simulation length and cell count
#'
#' Reads the merged all-endpoints results, groups the parameter combinations
#' that produced identical RF distances at every endpoint into lettered colour
#' groups, and for each modality draws one line plot against simulation length,
#' one against cell count, and a colour-matched table of the parameter
#' combinations. Cell count is taken as `2^endpoint`.
#'
#' @param run_id Run identifier; names the merged-results file and the output
#'   directory.
#' @param poss_recon_modals Modalities to plot, one figure set each; defaults to
#'   the global of the same name.
#' @param save_plots Whether to write the PNGs; defaults to `TRUE`.
#' @return The combined `grid.arrange` object for the last modality plotted.
#' @section Side effects: When `save_plots` is `TRUE`, writes
#'   `<modal>_simlength_plot.png`, `<modal>_numcells_plot.png`, and
#'   `<modal>_comb_plots.png` under `output/lineplots/<run_id>/`; that directory
#'   is expected to exist already.
#' @note Defined here for downstream use but not called by this script.
make_lineplot <- function(run_id, poss_recon_modals = poss_recon_modals, save_plots = TRUE){
  
  res <- read.csv(file.path('output', 'results', 'merged_results', run_id, paste0('merged_results_specs_', run_id, '_all_endpoints.csv')),
                  row.names = 1)
  
  diletters_grid <- expand.grid(LETTERS, LETTERS)
  diletters_vec <- apply(diletters_grid, MARGIN = 1, function(x){return(paste0(x[1], x[2]))})
  letters_diletters <- append(LETTERS, diletters_vec)
  
  for(modal in poss_recon_modals){
    
    modality_df <- res %>%
      filter(modality == modal)
    
    # # determine which parameters vary 
    # making overall plots ...
    modality_df <- modality_df %>%
      filter(!is.na(rf_dist)) %>%
      group_by(sampling_frac_bc, af_thresh_bc, af_thresh_mt,
               score_type_bc, score_type_mt, num_bc_integrations_bc, num_mito_genomes_mt) %>%
      mutate(param_combo = cur_group_id(),
             subrun_id = paste(mt_sub_run_id, bc_sub_run_id, sep = '_')) %>%
      arrange(endpoint)
    
    modality_df$num_cells <- 2**modality_df$endpoint
    
    # create a group x timepoint, values = rf_dist dataframe
    wide_res <- modality_df %>%
      pivot_wider(id_cols = param_combo, names_from = endpoint, values_from = rf_dist)
    
    unique_rf_dists <- unique(wide_res[, 2:ncol(wide_res)])
    
    color_groups <- list()
    for(i in seq_len(nrow(unique_rf_dists))){
      
      # find which param combos have identical rf dists at all endpoints
      
      rf_dist_combo <- unname(unlist(as.vector(unique_rf_dists[i, ])))
      
      matching_row_indices <- which(apply(wide_res[, 2:ncol(wide_res)], 1, function(row) all(row == rf_dist_combo)))
      equal_groups <- wide_res$param_combo[matching_row_indices]
      color_groups[[i]] <- equal_groups
    }
    
    # assign each group of identical RF dists a group (defined by a letter)
    # called color_groups becuase each group will be represented by a single color in plots
    # this list maps color groups to param combos
    names(color_groups) <- letters_diletters[1:length(color_groups)]
    
    # now reverse the mapping, from param combo to color group
    param_combo_to_color_group <- list()
    for(color_group in names(color_groups)){
      for(param_combo in color_groups[[color_group]]){
        param_combo_to_color_group[[param_combo]] <- color_group
      }
    }
    
    # use each row in res's param combo to determine which color group that row belongs to
    res_color_groups <- sapply(modality_df$param_combo, function(x) param_combo_to_color_group[[x]])
    
    # and insert this color group as a new feature to res
    modality_df$group <- res_color_groups
    
    # only keep the first occurrence of each unique set of param combos
    group_table <- modality_df %>%
      select(param_combo, group, subrun_id, sampling_frac_bc, af_thresh_bc, af_thresh_mt,
             score_type_bc, score_type_mt, num_bc_integrations_bc, num_mito_genomes_mt) %>%
      group_by(param_combo) %>%
      slice(1) %>%
      arrange(group)
    
    
    
    # renumber param combos to align with sorted group numbers
    group_table$param_combo <- seq(1, nrow(group_table))
    
    group_table <- group_table %>%
      rename(`Param Combo` = param_combo,
             `Plot Group` = group,
             `Joint \nSubrun` = subrun_id,
             `Cell \nSampling \nFrac` = sampling_frac_bc,
             `BC AF \nThresh` = af_thresh_bc,
             `MT AF \nThresh` = af_thresh_mt,
             `BC Score \nType` = score_type_bc,
             `MT Score \nType` = score_type_mt,
             `BC Ints` = num_bc_integrations_bc,
             `MT \nGenomes` = num_mito_genomes_mt)
    
    num_table_rows <- nrow(group_table)
    
    
    
    rand_col_pal <- distinctColorPalette(k = 50)
    
    simlength_plot <- ggplot(modality_df, aes(x = endpoint, y = rf_dist, color = group)) + 
      geom_line() +
      theme_bw() +
      scale_color_manual(values = rand_col_pal) +
      labs(x = 'Simulation Length',
           y = 'Normalized RF Distance',
           title = 'RF Distance Over Sim Length') +
      theme(legend.position = 'bottom') +
      guides(color = guide_legend(nrow = 1))
    
    # extract legend:
    table_build <- ggplot_gtable(ggplot_build(simlength_plot))
    legend_ind <- which(sapply(table_build$grobs, function(x) x$name == 'guide-box'))
    legend <- table_build$grobs[[legend_ind]]
    
    # now remove legend from simlength_plot
    simlength_plot <- simlength_plot + theme(legend.position = 'none')
    
    numcells_plot <- ggplot(modality_df, aes(x = num_cells, y = rf_dist, color = group)) + 
      geom_line() +
      theme_bw() +
      scale_color_manual(values = rand_col_pal) +
      labs(x = 'Number of Cells',
           y = 'Normalized RF Distance',
           title = 'RF Distance Over Number of Cells') +
      theme(legend.position = 'none')
    
    unique_plot_groups <- unique(group_table$`Plot Group`)
    line_cols <- rand_col_pal[1:length(unique_plot_groups)]
    cols <- matrix(NA, nrow = nrow(group_table), ncol = ncol(group_table))
    for(i in seq_len(nrow(group_table))){
      cols[i, ] <- line_cols[which(unique_plot_groups == group_table$`Plot Group`[i])]
    }
    colored_row_theme <- ttheme_minimal(core=list(fg_params = list(col = cols),
                                                  bg_params = list(col="#FFFFFF")),
                                        rowhead=list(bg_params = list(col=NA)),
                                        colhead=list(bg_params = list(col=NA)),
                                        padding = unit(c(20,4), 'pt'))
    
    grobbed_table <- tableGrob(group_table, theme = colored_row_theme, rows = NULL)
    
    # 4, then 1, then 10
    rel_size_top_plots <- max(c(4, ceiling(num_table_rows/4)))
    rel_size_top_plots_legend <- max(c(1, ceiling(num_table_rows/10)))
    rel_size_table <- max(c(10, num_table_rows))
    
    layout_mat <- rbind(
      matrix(data = rep(c(1,2), rel_size_top_plots),
             nrow = rel_size_top_plots, byrow = TRUE),
      matrix(data = rep(c(3,3), rel_size_top_plots_legend),
             nrow = rel_size_top_plots_legend, byrow = TRUE),
      matrix(data = rep(c(4,4), rel_size_table),
             nrow = rel_size_table, byrow = TRUE)
    )
    
    
    comb_plots <- grid.arrange(simlength_plot, numcells_plot, legend, grobbed_table, nrow = 3, layout_matrix = layout_mat)
    
    if(save_plots){
      ggsave(plot = simlength_plot, filename = file.path('output', 'lineplots', run_id, paste0(modal, '_simlength_plot.png')), width = 14, 
             height = 8, units = 'in')
      ggsave(plot = numcells_plot, filename = file.path('output', 'lineplots', run_id, paste0(modal, '_numcells_plot.png')), width = 14, 
             height = 8, units = 'in')
      ggsave(plot = comb_plots, filename = file.path('output', 'lineplots', run_id, paste0(modal, '_comb_plots.png')), width = 14, 
             height = max(8, num_table_rows/3), units = 'in')
    }
    
    
  }
  
  return(comb_plots)
  
}
