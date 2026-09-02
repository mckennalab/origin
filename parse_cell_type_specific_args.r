# Reshapes the cell-type transition matrices carried in the parsed JSON run
# parameters (`input_args$cell_type_dict`) into nested lists keyed by cell-type
# name, so that downstream code can read a transition probability as
# tm[[source_type]][[target_type]] instead of by row and column number.
# Sourced by sim5_code.R, after `input_args` and `cell_type_names` both exist.

#' Build named uninduced and induced cell-type transition lists
#'
#' Reads `input_args$cell_type_dict$uninduced_transition_matrix` and
#' `$induced_transition_matrix` from the enclosing script's global environment
#' and rewrites each as `tm[[source_type]][[target_type]]`. The matrices arrive
#' from JSON as a list of rows, so element `[[i]][j]` is the probability of
#' moving from the cell type at index `i` to the one at index `j`.
#'
#' @param cell_type_names Character vector of cell-type names, in the same order
#'   as the rows and columns of both transition matrices. Its default refers to
#'   the variable of the same name in the calling scope.
#' @return Named list with `uninduced_tm_list` and `induced_tm_list`, each a
#'   nested list of transition probabilities keyed by source name then target
#'   name.
#' @note Errors when either matrix has a different number of rows than
#'   `cell_type_names`; row widths are not checked, so a ragged matrix is
#'   carried through and surfaces later as a missing entry.
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
