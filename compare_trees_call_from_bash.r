suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phangorn)
  library(stringr)  
  library(RColorBrewer)
})

entering_dir <- getwd()

# this is great --v
output_dir_stem <- file.path('..', '..', '..', '..', 'output')

   
option_list <- list(
  make_option(c('-R', '--recon_tree_path'), type = 'character', default = NULL,
              help = 'path to reconstructed tree from iqtree'),
  make_option(c('-G', '--ground_truth_tree_path'), type = 'character', default = '',
              help = 'path to ground truth tree (based on cell population size)'),

  make_option(c('-I', '--run_id'), type = 'character', default = NULL,
              help = 'run id (numeric string)'),
  make_option(c('-S', '--savename_prefix'), type = 'character', default = NULL,
              help = 'savename prefix for tree image and results txt file'),
  make_option(c('-P', '--param_file'), type = 'character', default = NULL,
              help = 'name of paramter json file'),
  make_option(c('-T', '--treefile_dir_path'), type = 'character', default = NULL,
              help = 'path to the subdir where the .treefile lives')
  )

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser)

if(is.null(input_args$treefile_dir_path) || is.null(input_args$recon_tree_path) ||
   is.null(input_args$run_id) || is.null(input_args$savename_prefix)){
  stop('--treefile_dir_path, --recon_tree_path, --run_id, and --savename_prefix are required.')
}

treefile_dir_path <- normalizePath(input_args$treefile_dir_path, mustWork = TRUE)
if(!grepl('^/', input_args$recon_tree_path)){
  input_args$recon_tree_path <- file.path(treefile_dir_path, input_args$recon_tree_path)
}
input_args$recon_tree_path <- normalizePath(input_args$recon_tree_path, mustWork = TRUE)
if(!input_args$ground_truth_tree_path == ''){
  input_args$ground_truth_tree_path <- normalizePath(
    input_args$ground_truth_tree_path,
    mustWork = TRUE
  )
}

# if there is not a ground truth tree path provided, we need to find the corresponding ground truth tree
# for the provided recon tree path. to find the correct tree, we need the timepoint
# and the cell sampling fraction. we can get this info from the fasta file path that was used to build the tree

setwd(treefile_dir_path)

if(!input_args$ground_truth_tree_path == ''){ # if a ground truth tree path is provided, can just read that
  full_gt_path <- input_args$ground_truth_tree_path
  ground_truth_tree <- read.tree(full_gt_path)
} else{ # manually match to groundt truth tree at this timepoint
  
  timept <- str_extract(input_args$recon_tree_path, '(?<=time_)\\d+(\\.\\d+)?')

  ground_truth_trees <- list.files(file.path(output_dir_stem, 'processed_newicks', input_args$run_id))
  timept_pattern <- paste0('ground_truth_.*time_', timept)
  matching_gt_tree_path <- ground_truth_trees[grep(pattern = timept_pattern, x = ground_truth_trees)]
  full_gt_path <- file.path(output_dir_stem, 'processed_newicks', input_args$run_id, matching_gt_tree_path)
  
  if(length(matching_gt_tree_path) != 1){
    stop(sprintf(
      'Expected one ground-truth tree for timepoint %s; found %d.',
      timept,
      length(matching_gt_tree_path)
    ))
  }
  
  ground_truth_tree <- read.tree(full_gt_path)
  
  
}

# get cell population object from matching timepoint (to be used for coloring nodes)
timept <- str_extract(input_args$recon_tree_path, '(?<=time_)\\d+(\\.\\d+)?')
cell_pop_paths <- list.files(file.path(output_dir_stem, 'cell_populations', input_args$run_id))

timept_pattern <- paste0('.*time_', timept, '.*\\.rds')
matching_cell_pop_path <- cell_pop_paths[grep(pattern = timept_pattern, x = cell_pop_paths)]
if(length(matching_cell_pop_path) != 1){
  stop(sprintf(
    'Expected one cell-population file for timepoint %s; found %d.',
    timept,
    length(matching_cell_pop_path)
  ))
}
full_cell_pop_path <- file.path(output_dir_stem, 'cell_populations', input_args$run_id, matching_cell_pop_path)


# can build a wrapper around this if i want to allow the user to pick between ML and consensus trees
# and to specify whether midpoint root is wanted. 
save_image_of_tree <- function(tree_path, savename_prefix, cell_pop_path, midpt = FALSE){
  
  this_tree <- ape::read.tree(tree_path)
  
  if(midpt){
    this_tree <- midpoint(this_tree)
  }
  
  if(!dir.exists(file.path(output_dir_stem, 'tree_images', input_args$run_id))){
    dir.create(file.path(output_dir_stem, 'tree_images', input_args$run_id), recursive = TRUE)
  }
  
  png(file.path(output_dir_stem, 'tree_images', input_args$run_id, paste0(savename_prefix, '.png')), width = 1500, height = 4000, res = 300)
  plot(this_tree, cex = 0.4)
  dev.off()
  
  
}

plot_tree_with_color <- function(path_to_tree,
                                 path_to_cellpop,
                                 image_savename,
                                 recon_or_gt,
                                 color_internal_nodes = FALSE,
                                 color_map_list = NULL){
  
  if(!dir.exists(file.path(output_dir_stem, 'tree_images', input_args$run_id))){
    dir.create(file.path(output_dir_stem, 'tree_images', input_args$run_id), recursive = TRUE)
  }

  
  tree <- read.tree(path_to_tree)
  pop <- readRDS(path_to_cellpop)
  
  
  if(!is.null(color_map_list)){
    type_color_map_list <- color_map_list
  } else{
    cell_type_names <- unique(sapply(pop, function(cell) cell$celltype))
    if(length(cell_type_names) <= 8){
      palette_size <- max(3, length(cell_type_names))
      cell_type_colors <- brewer.pal(palette_size, 'Set2')[seq_along(cell_type_names)]
    } else{
      cell_type_colors <- colorRampPalette(brewer.pal(8, 'Set2'))(length(cell_type_names))
    }
    
    type_color_map_list <- setNames(cell_type_colors, cell_type_names)
    type_color_map_list <- type_color_map_list[!is.na(names(type_color_map_list))]
    
    
    
  }
  
  # Get node colors based on labels
  tip_colors <- sapply(tree$tip.label, function(cellname){
    celltype <- pop[[cellname]]$celltype
    if(celltype %in% names(type_color_map_list)){
      return(type_color_map_list[[celltype]])
    } else{
      return('grey70')
    }
  })

  
  pdf(file.path(output_dir_stem, 'tree_images', input_args$run_id, paste0(image_savename, '.pdf')), height = 25)
  # Plot tree with colored nodes
  # plot(tree, show.tip.label = TRUE, tip.color = tip_colors, 
  #      edge.color = "black", cex = 0.3, show.node.label = FALSE)
  plot(tree, show.tip.label = TRUE, tip.color = tip_colors, type = 'fan',
       edge.color = "black", cex = 0.3, show.node.label = FALSE)
  if(recon_or_gt == 'recon'){
    edgelabels(text = round(tree$edge.length, 3), frame = "n", cex = 0.2, col = 'orange')
  }
  
  if(recon_or_gt == 'gt'){
    if(color_internal_nodes){
      
      int_node_colors <- sapply(tree$node.label, function(cellname){
        celltype <- pop[[cellname]]$celltype
        return(type_color_map_list[[celltype]])
      })
      # Get node positions
      nodelabels(pch = 21, bg = int_node_colors, cex = 1)  # Add colored nodes
    }
  }
  
  dev.off()   
  
}

recon_tree <- read.tree(input_args$recon_tree_path)

# subset ground truth tree to only include those tip labels present in recon tree:
recon_tips <- recon_tree$tip.label
missing_recon_tips <- setdiff(recon_tips, ground_truth_tree$tip.label)
if(length(missing_recon_tips) > 0){
  stop(sprintf(
    'Reconstructed tree contains %d tips absent from the ground-truth tree: %s',
    length(missing_recon_tips),
    paste(head(missing_recon_tips, 10), collapse = ', ')
  ))
}
subset_gt_tree <- drop.tip(ground_truth_tree, setdiff(ground_truth_tree$tip.label, recon_tips))

rf_dist <-  phangorn::RF.dist(recon_tree, subset_gt_tree, normalize = TRUE)
print(paste0('rf_dist == ', rf_dist))

plot_tree_with_color(path_to_tree = input_args$recon_tree_path,
                     path_to_cellpop = full_cell_pop_path,
                     image_savename = input_args$savename_prefix,
                     color_internal_nodes = FALSE,
                     color_map_list = NULL,
                     recon_or_gt = 'recon')

plot_tree_with_color(path_to_tree = full_gt_path,
                     path_to_cellpop = full_cell_pop_path,
                     image_savename = paste0('gt_timept_', timept),
                     color_internal_nodes = TRUE,
                     color_map_list = NULL,
                     recon_or_gt = 'gt')


if(!dir.exists(file.path(output_dir_stem, 'rf_dist_files', input_args$run_id))){
  dir.create(file.path(output_dir_stem, 'rf_dist_files', input_args$run_id), recursive = TRUE)
}

res_file_path <- file.path(output_dir_stem, 'rf_dist_files', input_args$run_id, paste0(input_args$savename_prefix, '_rf.txt'))
close(file(res_file_path, open = 'w'))

cat(paste0(rf_dist, '\n'), file = res_file_path, append = TRUE)
cat(paste0(input_args$param_file, '\n'), file = res_file_path, append = TRUE)
