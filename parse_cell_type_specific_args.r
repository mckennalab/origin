make_cell_type_transition_lists <- function(cell_type_names = cell_type_names){ 
  
  return_list <- list()
  
  uninduced_TM <- input_args$cell_type_dict$uninduced_transition_matrix
  induced_TM <- input_args$cell_type_dict$induced_transition_matrix
  
  uninduced_TM_list <- list()
  induced_TM_list <- list()
  
  if(length(cell_type_names) != length(uninduced_TM)){
    stop('Uninduced transition matrix is incompatible with provided cell types.')
  }
  if(length(cell_type_names) != length(induced_TM)){
    stop('Induced transition matrix is incompatible with provided cell types.')
  }
  
  for(i in seq_along(cell_type_names)){
    for(j in seq_along(cell_type_names)){
      uninduced_TM_list[[cell_type_names[i]]][[cell_type_names[j]]] <- uninduced_TM[[i]][j]
      induced_TM_list[[cell_type_names[i]]][[cell_type_names[j]]] <- induced_TM[[i]][j]
    }
  }
  
  return_list[['uninduced_tm_list']] <- uninduced_TM_list
  return_list[['induced_tm_list']] <- induced_TM_list
  
  return(return_list)
  
  
}
