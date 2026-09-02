# convert_scoremats_to_csvs.r
#
# Post-processing utility that exports the RDS score matrices produced by
# `mut_to_scoremat.r` as human-readable CSVs. It is the last score-matrix step
# of the pipeline and is called by `bash_wrapper_all_combos.sh` as:
#
#   Rscript convert_scoremats_to_csvs.r --urid <run_id> \
#     --score_mat_path <base_path>/output/score_mats/<run_id>/matrices/
#
# and only when that `matrices/` directory exists — i.e. after `sim5_code.R`
# has run and written score matrices for the run.
#
# Input:  every `*.rds` directly inside `--score_mat_path` (the search is not
#         recursive, so the `af/` subdirectory of allelic-fraction matrices is
#         not included). Each file holds a score matrix, typically a sparse
#         `Matrix` object with cell row names and mutation column names.
# Output: one `<same-basename>.csv` per input, written to a `csvs/`
#         subdirectory of `--score_mat_path`, which is created if missing.
#         For the wrapper's invocation that is
#         `output/score_mats/<run_id>/matrices/csvs/`. Row and column names are
#         preserved by `write.csv`. Existing CSVs of the same name are
#         overwritten.
#
# `--urid` is accepted for call compatibility with the other `*_from_bash.r`
# scripts but is not used: all paths are derived from `--score_mat_path`.

suppressPackageStartupMessages({
  library(optparse)
})


# ---- Command-line interface ----

# As in the other scripts the bash wrapper calls (`sim5_code.R`,
# `compare_trees_call_from_bash.r`, `process_results_from_bash.r`), optparse's
# automatic `--help` option is disabled below via `add_help_option = FALSE`.
option_list <- list(
  make_option(c('-P', '--score_mat_path'), type = 'character', default = NULL,
              help = 'path to score mat directory (will overwrite urid input)'),
  make_option(c('-U', '--urid'), type = 'character', default = NULL,
              help = 'unique run id')
)


# ---- Resolve input and output directories ----

opt_parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
input_args <- parse_args(opt_parser) 

score_mat_path <- input_args$score_mat_path

# CSVs are written to a `csvs/` directory nested inside the score-matrix
# directory that was passed in, alongside the source RDS files.
csv_scoremat_path <- file.path(score_mat_path, 'csvs')

# Non-recursive listing: only score matrices sitting directly in the given
# directory are converted, not those in nested directories such as `af/`.
all_rds_files <- list.files(score_mat_path, pattern = '.*\\.rds',
                            full.names = TRUE)

if(!dir.exists(csv_scoremat_path)){
  dir.create(csv_scoremat_path, recursive = TRUE)
}


# ---- Convert each score matrix to CSV ----

# One CSV per RDS, keeping the source basename and swapping only the extension,
# so a matrix and its CSV export stay easy to pair up by name.
for(rds_path in all_rds_files){
  
  filename <- basename(rds_path)
  
  print(paste0('converting ', filename))
  
  csv_filename <- sub('\\.rds$', '.csv', filename, ignore.case = TRUE)
  
  csv_path <- file.path(csv_scoremat_path, csv_filename)
  
  # Score matrices are stored sparse; `as.matrix` densifies before writing, so
  # memory use here scales with cells x mutations rather than with the number
  # of non-zero entries.
  mat <- as.matrix(readRDS(rds_path))
  write.csv(mat, csv_path)
}




