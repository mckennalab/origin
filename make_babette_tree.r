# babette/BEAST 2 wrapper: turns a FASTA of reconstructed cell sequences into a
# consensus phylogeny and writes Newick files for downstream tree comparison.
# Requires the external BEAST 2 program (the conda env pins beast2=2.6.3 plus
# beagle-lib) as well as the R packages loaded below - babette shells out to it,
# so nothing here runs on a machine without BEAST 2 installed.
# Optional module: the source() line for it in sim5_code.R is commented out.

suppressPackageStartupMessages({
  library(babette)
  library(stringr)
  library(ape)
  library(ggplot2)
  library(phangorn)  
})


#' Pull the BEAST 2 state-file path out of parsed babette output
#'
#' @param beast2_output Parsed babette output; element 35 of its `output`
#'   character vector is read and everything from the first `/` onwards is taken
#'   as the path.
#' @return Single character path, or `NA` when that line holds no `/`.
#' @note The line index is hard-coded, so this breaks if BEAST 2 changes its
#'   banner. It is only reached from the commented-out parsing block inside
#'   `fasta_to_phylo()`.
extract_state_file_path <- function(beast2_output){
  text <- beast2_output$output[35]
  path <- str_extract(text, '\\/.*')
  return(path)
}

#' Extract the first Newick string from a BEAST 2 output file and save it
#'
#' Reads the whole file into one string, takes the first run that starts at `(`
#' and stops before any `<`, and appends the terminating `;` that the Newick
#' grammar requires.
#'
#' @param path_to_beast_output File to read the tree out of.
#' @param newick_path Destination file for the extracted Newick string.
#' @return Invisibly, the confirmation message; the tree itself is delivered
#'   through `newick_path`.
#' @section Side effects:
#'   Overwrites `newick_path` and prints a line naming it.
extract_newick_data <- function(path_to_beast_output, newick_path){
  
  lines <- readLines(path_to_beast_output)
  
  text <- paste(lines, collapse = '\n')
  
  res <- str_extract(text, '\\([^<]*')
  
  # newick string ends with ; 
  res <- paste0(res, ';')
  
  file_connection <-file(newick_path)
  writeLines(res, file_connection)
  close(file_connection)
  
  print(paste0('Newick data written to ', newick_path))
  # return(res)
  
  
}

#' Read a Newick file into an ape phylo object
#'
#' @param path_to_newick Path to a Newick file.
#' @return The `phylo` object produced by `ape::read.tree()`.
phylo_obj_from_newick <- function(path_to_newick){
  
  # all functions from ape
  # if want to return dendro obj
  # return(as.dendrogram(as.hclust.phylo(read.tree(path_to_newick))))
  
  # if want to return phylo obj
  return(read.tree(path_to_newick))
}

#' Infer a consensus phylogeny from a FASTA by running BEAST 2 through babette
#'
#' Assembles a babette inference model, runs BEAST 2 over `fasta_path`, picks up
#' the posterior tree file the run leaves in the working directory, and reduces
#' the posterior sample to a single consensus tree whose tips are then
#' relabelled.
#'
#' @details
#' The external BEAST 2 program must be installed and reachable by babette
#' (through `beastier`); this function only builds the model and reads the
#' results back.
#'
#' MCMC and prior settings depend on `inf_model`. With the default `'test'`, the
#' model comes from `create_test_inference_model()` with a default birth-death
#' tree prior and `create_test_mcmc(chain_length = 100000, store_every = 1000)`.
#' Any other value selects `create_inference_model()` with a birth-death prior
#' whose birth rate is pinned to 1 and death rate to 0.5 by degenerate uniform
#' distributions, after which `chain_length` (10000000) and `store_every` (1000)
#' are assigned onto the returned model object. Neither branch sets a burn-in.
#'
#' The posterior trees are read with `ape::read.nexus`, forced binary one at a
#' time with `multi2di`, and summarised with `ape::consensus(p = 0.2)`.
#'
#' @param fasta_path FASTA of per-cell sequences that BEAST 2 is run on.
#' @param this_run_id Run label. Referenced only from commented-out output-path
#'   code, so it does not affect the current path.
#' @param linstrings Lineage strings for tip relabelling. Referenced only from
#'   commented-out code.
#' @param return_phylo When `TRUE` the consensus tree is returned; when `FALSE`
#'   the function returns nothing and only the Newick files remain.
#' @param inf_model `'test'` selects the short test chain; any other value
#'   selects the long chain with fixed birth and death rates.
#' @param beast_options Forwarded to `bbt_run_from_model()`. Leave as `FALSE` so
#'   that `create_beast2_options()` fills in the defaults; the value is tested
#'   with `if(!beast_options)`, so it has to stay logical.
#' @param newick_out_path Path whose stem - everything before `.newick` - names
#'   the two Newick outputs.
#' @param site_model `'JC69'` or `'HKY'`; selects a babette site model object.
#' @return The consensus `phylo` object when `return_phylo` is `TRUE`, otherwise
#'   `NULL`.
#' @section Side effects:
#'   Runs BEAST 2, which drops `.trees`, `.log`, and `.csv` files into the
#'   working directory, then writes `<stem>_preRELABEL.newick` (tips as BEAST 2
#'   labelled them) and `<stem>_posRELABEL.newick` (each `_` in a tip label
#'   replaced by `.`). Progress messages are printed throughout.
#' @note The posterior sample is located by taking the newest `.trees` file in
#'   `./` rather than from the run's own return value, so a stale or concurrent
#'   `.trees` file in the working directory would be picked up instead. The site
#'   model object built from `site_model` is not passed into the inference
#'   model on the current code path.
fasta_to_phylo <- function(fasta_path, this_run_id, linstrings = NULL, return_phylo = TRUE, inf_model = 'test', 
                           beast_options = FALSE, newick_out_path = '',
                           site_model = 'JC69'){
  
  
  if(site_model == 'JC69'){
    beast_site_model <- create_jc69_site_model()
  }
  else if(site_model == 'HKY'){
    beast_site_model <- create_hky_site_model()
  }
  # else if(site_model == 'GTR'){
  #     beast_site_model <- create_
  # }
  # could do something like:
  # if(site_model == 'detect_sim_params'){
  # function that fully customizes background site model from simulation params
  # }
  

  
  
  #     create_bd_tree_prior(
  #   id = NA,
  #   birth_rate_distr = create_uniform_distr(),
  #   death_rate_distr = create_uniform_distr()
  # )
  
  # birth_rate_param <- create_param(id = 'birthRate', value = 1.0, lower = 1.0, upper = 1.0)
  # death_rate_param <- create_param(id = 'deathRate', value = 0.5, lower = 0.5, upper = 0.5)
  
  if(inf_model == 'test'){
    # inf_model <- create_test_inference_model(tree_prior = create_yule_tree_prior(birth_rate_distr = create_distr(name = 'uniform',
    #                                                                                                             id = 'fixed_birthrate',
    #                                                                                                             value = 1,
    #                                                                                                             lower = 1,
    #                                                                                                             upper = 1)))
    
    inf_model <- create_test_inference_model(tree_prior = create_bd_tree_prior(),
                                             mcmc = create_test_mcmc(
                                               chain_length = 100000,
                                               store_every = 1000)
                                             
    #   id = 'fixed_birth_death',
    #   birth_rate_distr = create_uniform_distr(value = 1, lower = 1, upper = 1),
    #   death_rate_distr = create_uniform_distr(value = 0.5, lower = 0.5, upper = 0.5)
    #   
    # )
    )
    # inf_model$chain_length <- 100000
    # inf_model$store_every <- 1000
    
  } else{
    # inf_model <- create_inference_model(tree_prior = create_yule_tree_prior(birth_rate_distr = create_distr(name = 'uniform',
    #                                                                                                             id = 'fixed_birthrate',
    #                                                                                                             value = 1,
    #                                                                                                             lower = 1,
    #                                                                                                             upper = 1)))
    
    inf_model <- create_inference_model(tree_prior = create_bd_tree_prior(
      id = 'fixed_birth_death',
      birth_rate_distr = create_uniform_distr(value = 1, lower = 1, upper = 1),
      death_rate_distr = create_uniform_distr(value = 0.5, lower = 0.5, upper = 0.5)
    ))
    inf_model$chain_length <- 10000000
    inf_model$store_every <- 1000
    # inf_model$site_model$name <- 'HKY'
    
  }
  
  # print('inference model created')

  if(!beast_options){
    # xml_path <- 'here_test_output.xml'
    # create_beast2_input_file(
    #   fasta_filename = fasta_path,
    #   beast2_input_filename = xml_path,
    #   inference_model = inf_model
    # )
    # beast_options <- create_beast2_options(input_filename = xml_path)
    beast_options <- create_beast2_options()
    # beast_options <- create_beast2_options(
    #   output_trees_filename = paste0("/home/Kiewit/f005c3x/.cache/beautier/", this_run_id, "_output_trees.trees"),
    #   output_log_filename = paste0("/home/Kiewit/f005c3x/.cache/beautier/", this_run_id, "_output_log.log"),
    #   output_state_filename = paste0("/home/Kiewit/f005c3x/.cache/beastier/", this_run_id, "_output_state.state")
    # )
  }
  
  # print(beast_options$input_filename)
  
  # print('created input file')
  
  
  # this_beast2_input <- create_beast2_input(
  #   input_filename = fasta_path,
  #   site_model = beast_site_model,
  #   clock_model = create_strict_clock_model(),
  #   tree_prior = create_yule_tree_prior(),
  #   mcmc = create_mcmc(chain_length = 1000000)
  # )
  
  
  
  out <- bbt_run_from_model(
    fasta_filename = fasta_path,
    inference_model = inf_model,
    beast2_options = beast_options
    # beast2_input = this_beast2_input
  )
  
  
  # output_dir_path <- '/home/Kiewit/f005c3x/.cache/beautier/' # if == 'test'
  output_dir_path <- './'
  files <- list.files(output_dir_path, full.names = TRUE)
  file_info <- file.info(files)
  
  # each babette run is associated with three output files in this directory:
  # .trees, .csv, and .log
  # we are interested in the .trees file when it comes to building a consensus tree
  
  # most_recent_files <- rownames(file_info[order(file_info$ctime, decreasing = TRUE), ])
  # trees_paths <- most_recent_files[grepl(pattern = '\\.trees$', x = most_recent_files)]
  # most_recent_trees_path <- rownames(trees_paths)[1]
  
  # extract only the .trees file name:
  
  most_recent_files <- rownames(file_info[order(file_info$ctime, decreasing = TRUE), ])
  trees_paths <- most_recent_files[grepl(pattern = '\\.trees$', x = most_recent_files)]
  most_recent_trees_path <- trees_paths[1]
  print(paste0('most_recent_trees_path == ', most_recent_trees_path))
  
  # most_recent_trees_path <- three_most_recent[grepl(pattern = '\\.trees', x = three_most_recent)]
  
  posterior_trees <- ape::read.nexus(most_recent_trees_path)
  print('in here for consensus 1')
  posterior_trees_binary <- lapply(posterior_trees, multi2di) # force into binary structure when generating consensus
  print('in here for consensus 2')
  cons_tree <- ape::consensus(posterior_trees_binary, p=0.2)
  print('finished making consensus tree')
  
  # write newick of pre-tip-change tree:
  newick_out_path_stem <- str_split(newick_out_path, '\\.newick')[[1]][1]
  pre_relabel_path <- paste0(newick_out_path_stem, '_preRELABEL.newick')
  post_relabel_path <- paste0(newick_out_path_stem, '_posRELABEL.newick')
  write.tree(cons_tree, file = pre_relabel_path)
  
  # ################## mcc tree:
  # best_tree <- phangorn::maxCladeCred(posterior_trees, tree = TRUE)
  
  # print('old tip labels of cons tree == ')
  # print(cons_tree$tip.label)
  print('made it to before new tip labels')
  new_tip_labs <- unname(unlist(sapply(cons_tree$tip.label, function(old_lab){
    return(gsub(pattern = '_', replacement = '.', x = old_lab))
  })))
  print('made it to after new tip labels')
  
  cons_tree$tip.label <- new_tip_labs
  write.tree(cons_tree, post_relabel_path)
  # print('new tip labels of cons tree == ')
  # print(cons_tree$tip.label)
  

  ##########################
  # cat(paste0('out$trees == ', out$trees, '\n'), file = 'no_strings.txt', append = TRUE)
  # 
  # # print('bbt finished')
  # 
  # parsed_output <- parse_beast2_output(out = out,
  #                                      inference_model = inf_model)
  # 
  # # print('output parsed')
  # 
  # output_path <- extract_state_file_path(parsed_output)
  # 
  # # print('path extracted')
  # 
  # if(newick_out_path == ''){
  #   newick_dir = "./output/beast_newicks"
  #   if(!dir.exists(newick_dir)){
  #     dir.create(newick_dir, recursive = TRUE)
  #   }
  #   rand_run_id <- paste(sample(seq(1, 10), 8), collapse = '')
  #   newick_out_path <- file.path(newick_dir, rand_run_id)
  # }
  # 
  # extract_newick_data(output_path, newick_path = newick_out_path)
  # ##########################
  
  # print('newick written')
  
  # can either return the dendro object, or return nothing since the Newick file has been written
  if(return_phylo){
    
    # ##########################
    # dendro_obj <- phylo_obj_from_newick(newick_out_path)
    # # print(paste0('class(dendro_obj) == ', class(dendro_obj)))
    # print(paste0('newick out path == ', newick_out_path))
    # 
    # 
    # resp_linstrings <- unname(unlist(sapply(dendro_obj$tip.label, function(old_num){
    #   return(linstrings[as.integer(old_num) + 1])
    # })))
    # 
    # dendro_obj$tip.label <- resp_linstrings
    # # print('dendro obj created')
    # return(dendro_obj)
    # 
    # ##########################
    
    return(cons_tree)
    
  }else{
    return()
  }
  
}


# NOT_TEST_barebones_16_variable_only <- fasta_to_phylo(fasta_path = 'simple_16cell_only_variable.fasta',
#                                                       return_phylo = TRUE,
#                                                       inf_model = 'TEST',
#                                                       newick_out_path = './output/beast_newicks/PRIOR_NOT_TEST_simple_16cell_only_variable.newick')
# adjusted_cell_nums <- unname(unlist(sapply(NOT_TEST_barebones_16_variable_only$tip.label, function(x){as.integer(x)+1})))
# new_tip_labels <- paste0('Cell_', adjusted_cell_nums)
# NOT_TEST_barebones_16_variable_only$tip.label <- new_tip_labels