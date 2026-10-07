# Exercises the code shown in vignettes/gillespie-quickstart.Rmd. The vignette
# chunks are eval = FALSE so building the package stays fast, which means these
# tests are the only thing keeping the documented examples from rotting. If a
# parameter name or a return field changes, this fails and the vignette needs
# the same edit.

quickstart_population_params <- function() {
  list(
    num_init_cells = 1,
    sim_length = list(6),
    random_seed = 1,
    editing_induction = list(timepoint = 0, num_cells = NULL, frac_cells = 1),
    differentiation_induction = list(timepoint = 0, num_cells = NULL,
                                     frac_cells = 1),
    cell_type_dict = list(
      founder_cell_type = "progenitor",
      cell_type_params = list(
        progenitor = list(cell_cycle_length = 1,
                          death_per_cell_cycle_prob = 0.01)
      ),
      uninduced_transition_matrix = list(list(1)),
      induced_transition_matrix = list(list(1))
    )
  )
}

test_that("the vignette's population example runs and returns its documented fields", {
  population <- simulate_gillespie_population(
    quickstart_population_params(),
    end_time = 6, seed = 1, max_cells = 100000L, show_progress = FALSE
  )

  expect_true(all(c("nodes", "terminal_nodes", "counts", "stop_reason") %in%
                    names(population)))
  expect_true(all(c("founders", "divisions", "deaths", "nodes",
                    "surviving_cells") %in% names(population$counts)))
  # Columns the vignette tells the reader to use.
  expect_true(all(c("generation", "cell_type", "is_terminal", "alive_at_end",
                    "died", "birth_time", "end_time") %in%
                    names(population$nodes)))
  expect_gt(population$counts[["divisions"]], 0)
  expect_equal(nrow(population$nodes), population$counts[["nodes"]])

  # The vignette warns that is_terminal includes lineages that died, so the
  # alive-at-end subset is the smaller one.
  terminal <- population$nodes$is_terminal %in% c(TRUE, "TRUE")
  alive <- population$nodes$alive_at_end %in% c(TRUE, "TRUE")
  expect_lte(sum(alive), sum(terminal))
  expect_equal(sum(alive), population$counts[["surviving_cells"]])
})

test_that("a run with no death loses no lineages", {
  params <- quickstart_population_params()
  params$cell_type_dict$cell_type_params$progenitor$death_per_cell_cycle_prob <- 0
  population <- simulate_gillespie_population(
    params, end_time = 4, seed = 2, max_cells = 100000L, show_progress = FALSE
  )
  expect_equal(population$counts[["deaths"]], 0)
  expect_equal(sum(population$nodes$died %in% c(TRUE, "TRUE")), 0)
})

test_that("the recording pipeline writes the outputs the vignette lists", {
  skip_if_not_installed("Matrix")

  editing_state <- list(
    be_mutations_per_target_per_division = 0.02,
    nuc_insertions_per_target_per_division = 0,
    nuc_deletions_per_target_per_division = 0,
    bc_bg_insertion_prob_per_division = 0,
    bc_bg_deletion_prob_per_division = 0,
    bc_substitution_model = "JC", bc_sub_model_params = list(0),
    mt_substitution_model = "JC", mt_sub_model_params = list(0),
    mt_bg_insertion_prob_per_division = 0,
    mt_bg_deletion_prob_per_division = 0
  )
  params <- quickstart_population_params()
  params$sim_length <- list(4)
  params$bc_length <- 300
  params$max_bc_ints_per_cell <- list(2)
  params$bc_nuc_composition <- list(frac_a = 0.25, frac_g = 0.25,
                                    frac_c = 0.25, frac_t = 0.25)
  params$be_conversion_pattern <- "A --> G"
  params$be_targets <- list(
    num_targets = 8, config = "S:10:24",
    edit_rate_class_fractions = list(high = 0.5, medium = 0, low = 0.5),
    editing_window = list(size = 0, decaying = FALSE, close_after_edit = FALSE)
  )
  progenitor <- params$cell_type_dict$cell_type_params$progenitor
  progenitor$bc_invariant_sites <- 0
  progenitor$mt_invariant_sites <- 0
  progenitor$induced_editing_params <- editing_state
  progenitor$uninduced_editing_params <- editing_state
  params$cell_type_dict$cell_type_params$progenitor <- progenitor

  output_dir <- file.path(tempdir(), "origin-quickstart-test")
  unlink(output_dir, recursive = TRUE)
  result <- run_gillespie_lineage_pipeline(
    params, params_path = NULL, output_dir = output_dir,
    overrides = list(progress = FALSE, seed = 1,
                     modalities = list("barcode"))
  )

  expect_true(all(c("population", "output_dir") %in% names(result)))
  for (name in c("physicell_lineage_full.nwk", "physicell_lineage_sampled.nwk",
                 "lineage_nodes.csv.gz", "terminal_cells.csv.gz",
                 "barcode_binary_score_matrix.csv.gz",
                 "barcode_target_layout.csv.gz", "run_manifest.csv.gz")) {
    expect_true(file.exists(file.path(output_dir, name)),
                info = paste("missing documented output:", name))
  }

  # The vignette states the score matrix is bc_length x integrations wide and
  # one row per cell alive at the end.
  score <- read.csv(
    gzfile(file.path(output_dir, "barcode_binary_score_matrix.csv.gz")),
    row.names = 1, check.names = FALSE
  )
  expect_equal(ncol(score), 300 * 2)
  expect_equal(nrow(score), result$population$counts[["surviving_cells"]])
  expect_true(all(grepl("^int_[0-9]+_pos_[0-9]+$", colnames(score))))
  unlink(output_dir, recursive = TRUE)
})
