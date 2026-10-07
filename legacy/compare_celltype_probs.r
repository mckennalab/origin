# compare_celltype_probs.r
#
# Standalone diagnostic script, outside the `bash_wrapper_all_combos.sh`
# pipeline. It answers two questions about the per-cell-type mutation
# parameters that `sim5_code.R` constructs for a single run:
#   1. Is the base-editing edit-rate-class list (`basepos_erc_be_list`, which
#      maps a barcode base position to 'High' / 'Medium' / 'Low') really shared
#      by every cell type, as `sim5_code.R` assumes?
#   2. How does the per-position barcode-target transition probability differ
#      between cell types?
#
# Inputs: `celltype_prob_concordance/<urid>/*.rds`, resolved relative to the
# working directory, so run this from the repo root. Each file is a named list
# keyed by barcode base position; the names of the files encode cell type
# (`ct1`..`ct6`), mutation type (transition / transversion / insertion /
# deletion) and induction state (`induced_editing_params` /
# `uninduced_editing_params`).
#
# Prerequisite: those RDS files are only written when the `saveRDS` block near
# the end of the per-cell-type / per-induction parameter loop in `sim5_code.R`
# is re-enabled — it is commented out in the committed source — and a
# simulation has been run. `urid` below is hard-coded to one run id and must be
# edited by hand to point at a different run.
#
# Invocation: `Rscript compare_celltype_probs.r`, or step through it
# interactively. It takes no command-line arguments and writes no files; the
# pairwise equality checks are printed to the console and the barplot is drawn
# to the active graphics device.

library(dplyr)
library(tidyr)
library(ggplot2)


# ---- Run selection and parameter-file discovery ----

# Hard-coded run id; change this to inspect a different simulation run.
urid <- '8348251423997'
celltype_names <- paste0('ct', seq(1, 6))
mut_types <- c('transition', 'transversion', 'insertion', 'deletion')

path_to_prob_dir <- file.path('celltype_prob_concordance', urid)
all_rds_files <- list.files(path_to_prob_dir, full.names = TRUE)

# Split the run's parameter files by induction state. 'uninduced' is the
# distinguishing substring, so everything else is an induced-editing file.
uninduced_paths <- grep(pattern = 'uninduced', x = all_rds_files, value = TRUE)
induced_paths <- setdiff(all_rds_files, uninduced_paths)


# ---- Collect per-cell-type barcode-target transition probabilities ----

# Only the transition parameters are gathered here; `mut_types` above lists the
# other three mutation classes but the loop is restricted to 'transition'.
# `transition_dfs` is assigned inside the loop but survives it, because an R
# `for` body shares the enclosing environment — the plotting section below
# relies on that.
for(mut_type in c('transition')){
  paths_with_mut <- grep(pattern = mut_type, x = induced_paths, value = TRUE)
  
  transition_dfs <- list()
  for(ct in celltype_names){
    # One file per cell type is expected to match; the cell-type tag ('ct1'..)
    # is embedded in the filename.
    ct_path <- grep(pattern = ct, x = paths_with_mut, value = TRUE)
    
    # The stored object is a named list of per-base-position probabilities, so
    # the values become `probs` and the names (base positions) become `pos`.
    params <- readRDS(ct_path)
    df <- data.frame('probs' = as.numeric(params),
                     'pos' = names(params),
                     'celltype' = ct)
    transition_dfs[[ct]] <- df
    
  }
}

# ---- Check the base-editing edit-rate-class list across cell types ----

# `basepos_erc_be_list` maps each barcode base position to its edit-rate class
# ('High' / 'Medium' / 'Low'). `sim5_code.R` builds it once and reuses it for
# every cell type, on the assumption that a target's class does not depend on
# cell type even though the underlying numeric rates do. The pairwise
# `all.equal` comparison below is the check on that assumption.
erc_ct1 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct1_induced_editing_params.rds'))
erc_ct2 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct2_induced_editing_params.rds'))
erc_ct3 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct3_induced_editing_params.rds'))
erc_ct4 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct4_induced_editing_params.rds'))
erc_ct5 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct5_induced_editing_params.rds'))
erc_ct6 <- readRDS(file.path(path_to_prob_dir, 'basepos_erc_be_listct6_induced_editing_params.rds'))

all_erc_lists <- list(erc_ct1, erc_ct2, erc_ct3,
                      erc_ct4, erc_ct5, erc_ct6)
# Compare every ordered pair of cell types, skipping the diagonal, and print
# the result of each comparison (TRUE, or the `all.equal` difference report).
for(i in seq_len(length(all_erc_lists))){
  for(j in seq_len(length(all_erc_lists))){
    if(i == j){
      next
    }
    print(paste0('all ct', i, ' equal all ct', j, ': ', all.equal(all_erc_lists[[i]], all_erc_lists[[j]])))
  }
}


# ---- Wide position-by-cell-type comparison table ----

# Stack the per-cell-type frames, then reshape to one row per base position and
# one column per cell type so the transition probabilities can be eyeballed
# side by side. Positions absent for a cell type become NA.
stacked_df <- do.call(rbind, transition_dfs)
compare_pos_df <- stacked_df %>%
  pivot_wider(id_cols = pos,
              names_from = celltype,
              values_from = probs,
              values_fill = list(probs = NA)) %>%
  arrange(pos)

# ---- Grouped barplot of transition probability by position and cell type ----

# `pos` is a character column, so it is treated as a discrete axis and the bars
# are ordered lexicographically rather than numerically. Labels are rotated
# because there is one group of bars per target base position. The plot is
# returned to the active device; nothing is written to disk.
ggplot(stacked_df, aes(x = pos, y = probs, fill = celltype)) + 
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))
