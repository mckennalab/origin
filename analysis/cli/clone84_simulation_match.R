#!/usr/bin/env Rscript

# Can the simulation framework reproduce clone 84's BASELINE recording?
#
# The reference is a 10-barcode, 250-cell, 272-site subset of the real clone 84
# recording (analysis/cli/clone84_reference_subset.R). This runs the Gillespie
# pipeline configured to the same geometry and time window and compares two
# things the recorder has to get right:
#
#   edit burden    edits per cell, over observed sites
#   heterogeneity  the distribution of per-site cumulative editing rate
#
# Geometry. The real assay scores 272 sites per barcode, so the barcode is
# configured as 272 positions that are all targets, with a zero-width editing
# window, and an all-A composition under an A->G conversion so that every
# position is genuinely editable. Without the composition constraint only about
# a quarter of positions would carry the target base and the simulated barcode
# would be four times narrower than the assay it is being compared with.
#
# Population. Thirty divisions of pure birth is 2^30 cells, so the run is held
# to a workable size with per-cycle death: expected offspring is 2 * (1 - death),
# which at death = 0.40 grows 1.2-fold per cycle and reaches a few thousand cells
# after 30 cycles from a modest founder pool. Cells are then subsampled to 250.
# Realised generation depth is measured rather than assumed, and reported.
#
# Sampling floor. Rates are compared as OBSERVED in 250 cells, not as true
# probabilities, because at this depth they are not the same thing: a site
# editing at 0.2% shows zero edits in 250 cells about two thirds of the time.
# The reference subset reads 23.4% of sites at exactly zero against 0.82% in the
# full 9,349-cell recording, and that spike is the sampling floor rather than a
# property of the construct. The calibration therefore passes candidate rates
# through the same binomial sampling, at the reference's own per-site
# denominators, and the simulated matrix is given matched per-site dropout.
# Calibrating true probabilities against observed rates would have absorbed the
# floor into the fitted dispersion.
#
# Rate model. The simulator draws each target's per-division rate from
# Gamma(shape, mean/shape) in draw_physicell_target_rates(), with shape defaulting to
# 0.5 in the shipped code. Two configurations are run:
#
#   shipped   shape = 0.5, mean calibrated against the reference
#   fitted    shape and mean both calibrated against the reference
#
# The shape is set through be_targets$edit_rate_dispersion_shape, which defaults
# to 0.5. Reporting both is the point: it separates "the framework cannot
# reproduce this" from "the framework's dispersion default is wrong".
#
# A third configuration drops the family altogether:
#
#   empirical measured per-site rates, sampled directly
#
# Neither gamma fit reproduces the reference shape -- the best fitted one still
# reads 0.236 mean site rate against 0.305, with 10.4% saturated against 17.3%
# and 54.1% near-dead against 38.7% -- because one dispersion parameter cannot
# raise the top tail without also deepening the bottom one. The empirical mode
# takes the rates from the FULL recording rather than the 250-cell reference:
# at a median 7,364 cells per site the sampling error on each is under 0.006, so
# the empirical distribution is essentially the true one, whereas at 250 cells
# 23% of sites read exactly zero purely from the floor. Cumulative fractions are
# inverted to per-division rates by r = 1 - (1 - p)^(1/divisions) and passed as
# be_targets$edit_rate_empirical.
#
# The gamma fits already correct for the sampling floor by resampling candidate
# rates through the reference's own denominators, so this is a comparison of
# rate MODELS, not of one model against a better-observed target. What it does
# also change is the class stratification, which the empirical path drops; see
# physicell_target_empirical_rates() for why.
#
# Arguments:
#   --reference=<dir>    output of clone84_reference_subset.R (required)
#   --output-dir=<dir>   destination (required)
#   --empirical-input=<csv>  full clone 84 matrix. When given, a third
#                        configuration samples per-target rates from the
#                        measured per-site rates instead of a gamma.
#   --empirical-min-cells=<n>  sites needing fewer observations than this are
#                        dropped from the empirical pool. Default 50.
#   --empirical-response=<csv>  measured `rate,cumulative` response of the
#                        simulator at this geometry. When given, target
#                        fractions are inverted through it instead of through
#                        1 - (1 - p)^(1/divisions). The analytic form is wrong
#                        here: the simulator's implied exposure is about 9.5
#                        divisions at low rates rising past 20 at high ones, not
#                        the nominal 30, so inverting analytically undershoots
#                        every target and collapses the saturated tail.
#   --configurations=<list>  which of shipped_shape0.5, fitted_shape,
#                        empirical_rates to run. Default all available. Rows for
#                        configurations not run are carried over from the
#                        existing summary rather than dropped.
#   --divisions=<n>      cell cycles to simulate. Default 30.
#   --cells=<n>          cells to subsample. Default 250.
#   --integrations=<n>   Default 10.
#   --sites=<n>          sites per integration. Default 272.
#   --founders=<n>       starting cells. Default 20.
#   --death=<p>          death probability per cell cycle. Default 0.40.
#   --sim-time=<t>       simulated time units. Divisions accumulate at about two
#                        per time unit, so this defaults to divisions / 2 and the
#                        realised generation depth is measured and reported
#                        rather than assumed.
#   --seed=<n>           Default 1.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
reference_dir <- value_after("--reference")
output_dir <- value_after("--output-dir")
if (is.null(reference_dir) || is.null(output_dir)) {
  stop("--reference and --output-dir are required.", call. = FALSE)
}
reference_dir <- normalizePath(reference_dir, mustWork = TRUE)
empirical_input <- value_after("--empirical-input")
empirical_min_cells <- as.integer(value_after("--empirical-min-cells", "50"))
empirical_response <- value_after("--empirical-response")
requested_configurations <- {
  raw <- value_after("--configurations", "")
  if (nzchar(raw)) strsplit(raw, ",", fixed = TRUE)[[1L]] else character()
}
divisions <- as.integer(value_after("--divisions", "30"))
cell_count <- as.integer(value_after("--cells", "250"))
integration_count <- as.integer(value_after("--integrations", "10"))
site_count <- as.integer(value_after("--sites", "272"))
founders <- as.integer(value_after("--founders", "20"))
death <- as.numeric(value_after("--death", "0.40"))
sim_time <- as.numeric(value_after("--sim-time", as.character(divisions / 2)))
seed <- as.integer(value_after("--seed", "1"))

source_root <- normalizePath(file.path(dirname(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])), "..", ".."),
  mustWork = TRUE)
source(file.path(source_root, "load_origin.R"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

reference <- read.csv(file.path(reference_dir, "reference_sites.csv"),
                      stringsAsFactors = FALSE)
reference_cells <- read.csv(file.path(reference_dir, "reference_cells.csv"),
                            stringsAsFactors = FALSE)
observed_rate <- reference$editing_rate[!is.na(reference$editing_rate)]
quantile_grid <- seq(0.01, 0.99, by = 0.01)
observed_quantiles <- stats::quantile(observed_rate, probs = quantile_grid)

# ---- Calibrate against the reference ----------------------------------------
reference_observed <- reference$observed[!is.na(reference$editing_rate)]
cumulative_from <- function(mean_rate, shape, draws = 200000L) {
  rate <- stats::rgamma(draws, shape = shape, scale = mean_rate / shape)
  probability <- 1 - (1 - pmin(rate, 1 - .Machine$double.eps))^divisions
  # Read the rate off the same number of cells the assay did, so the comparison
  # is observed-against-observed.
  denominator <- sample(reference_observed, draws, replace = TRUE)
  stats::rbinom(draws, size = denominator, prob = probability) / denominator
}
distance_for <- function(mean_rate, shape) {
  set.seed(seed)
  mean(abs(stats::quantile(cumulative_from(mean_rate, shape),
                           probs = quantile_grid) - observed_quantiles))
}
shipped_mean <- stats::optimize(
  function(log_mean) distance_for(exp(log_mean), 0.5),
  interval = log(c(1e-4, 0.5))
)$minimum
shipped_mean <- exp(shipped_mean)
fitted <- stats::optim(
  c(log(0.02), log(0.5)),
  function(par) distance_for(exp(par[1]), exp(par[2])),
  method = "Nelder-Mead"
)
fitted_mean <- exp(fitted$par[1])
fitted_shape <- exp(fitted$par[2])
cat(sprintf("[calibrate] shipped shape 0.5 -> mean %.5f (distance %.4f)\n",
            shipped_mean, distance_for(shipped_mean, 0.5)))
cat(sprintf("[calibrate] fitted shape %.4f -> mean %.5f (distance %.4f)\n",
            fitted_shape, fitted_mean, distance_for(fitted_mean, fitted_shape)))

# ---- Empirical per-target rates ----------------------------------------------
# Read from the full recording, not the reference subset. The reference reads
# each site on ~250 cells, where a site editing at 0.5% shows zero about a
# quarter of the time; the full recording reads a median 7,364 cells per site,
# where the standard error on any one rate is under 0.006. Using the subset
# would hand the simulator the sampling floor as though it were the construct.
empirical_rates <- NULL
if (!is.null(empirical_input)) {
  cat("[empirical] reading the full recording\n")
  full <- data.table::fread(normalizePath(empirical_input, mustWork = TRUE),
                            header = TRUE, colClasses = "character",
                            showProgress = FALSE)
  target_columns <- grepl("-[0-9]+$", colnames(full))
  observed_counts <- vapply(which(target_columns), function(column) {
    sum(full[[column]] %in% c("0", "1"))
  }, numeric(1))
  cumulative <- vapply(which(target_columns), function(column) {
    values <- full[[column]]
    kept <- values %in% c("0", "1")
    if (sum(kept) < empirical_min_cells) return(NA_real_)
    mean(values[kept] == "1")
  }, numeric(1))
  rm(full)
  cumulative <- cumulative[is.finite(cumulative)]
  # Cumulative edit fraction at the end of the run, inverted to the per-division
  # rate the simulator consumes. A site read at exactly 1 gives only a lower
  # bound on its rate, so it is capped rather than sent to infinity.
  if (is.null(empirical_response)) {
    empirical_rates <- 1 - (1 - pmin(cumulative, 1 - 1e-9))^(1 / divisions)
    cat("[empirical] inverted analytically; see --empirical-response\n")
  } else {
    # Invert through what the simulator actually does rather than through the
    # nominal division count. The response is monotone, so interpolating the
    # curve backwards is well defined; ties are resolved to the lowest rate
    # that reaches the fraction.
    response <- utils::read.csv(normalizePath(empirical_response,
                                              mustWork = TRUE),
                                stringsAsFactors = FALSE)
    if (!all(c("rate", "cumulative") %in% names(response))) {
      stop("--empirical-response needs rate and cumulative columns.",
           call. = FALSE)
    }
    response <- response[order(response$cumulative, response$rate), ]
    response <- response[!duplicated(response$cumulative), ]
    empirical_rates <- stats::approx(response$cumulative, response$rate,
                                     xout = cumulative, rule = 2)$y
    cat(sprintf("[empirical] inverted through %d response points from %s\n",
                nrow(response), basename(empirical_response)))
  }
  empirical_rates <- pmin(pmax(empirical_rates, 0), 1 - .Machine$double.eps)
  cat(sprintf("[empirical] %d sites (>= %d cells, median %.0f); cumulative mean %.4f, dead %.3f, saturated %.3f\n",
              length(empirical_rates), empirical_min_cells,
              stats::median(observed_counts), mean(cumulative),
              mean(cumulative == 0), mean(cumulative == 1)))
  cat(sprintf("[empirical] per-division rate: mean %.5f, median %.5f\n",
              mean(empirical_rates), stats::median(empirical_rates)))
}

# ---- Simulation parameters ---------------------------------------------------
editing_state <- list(
  be_mutations_per_target_per_division = shipped_mean,
  nuc_insertions_per_target_per_division = 0,
  nuc_deletions_per_target_per_division = 0,
  mt_mutations_per_genome_per_division = 0,
  num_mt_genomes = 0,
  bc_bg_insertion_prob_per_division = 0,
  bc_bg_deletion_prob_per_division = 0,
  bc_substitution_model = "JC",
  bc_sub_model_params = list(0),
  mt_substitution_model = "JC",
  mt_sub_model_params = list(0),
  mt_bg_insertion_prob_per_division = 0,
  mt_bg_deletion_prob_per_division = 0
)
build_params <- function(mean_rate, shape, empirical = NULL) {
  state <- editing_state
  state$be_mutations_per_target_per_division <- mean_rate
  # The empirical path never reads the shape, but the block still has to be a
  # valid one, so it carries the default rather than the NA the summary uses to
  # mark "no gamma here".
  if (!is.null(empirical) || !is.finite(shape)) shape <- 0.5
  list(
    num_init_cells = founders,
    sim_length = list(sim_time),
    random_seed = seed,
    bc_length = site_count,
    # Required even with mitochondrial recording switched off (num_mt_genomes
    # is 0 above); the model is prepared before the genome count is consulted.
    mito_genome_length = 16600,
    max_bc_ints_per_cell = list(integration_count),
    # All-A so that every one of the 272 positions is an eligible A->G target;
    # see the header note on assay width.
    bc_nuc_composition = list(frac_a = 1, frac_g = 0, frac_c = 0, frac_t = 0),
    be_conversion_pattern = "A --> G",
    be_targets = list(
      num_targets = site_count,
      config = "S:1:0",
      edit_rate_class_fractions = list(high = 0.4, medium = 0, low = 0.6),
      edit_rate_dispersion_shape = shape,
      # Present only in the empirical configuration; when it is, the mean and
      # the shape above are ignored, and so is the class split.
      edit_rate_empirical = empirical,
      editing_window = list(size = 0, decaying = FALSE, close_after_edit = FALSE)
    ),
    editing_induction = list(timepoint = 0, num_cells = NULL, frac_cells = 1),
    differentiation_induction = list(timepoint = 0, num_cells = NULL,
                                     frac_cells = 1),
    cell_type_dict = list(
      founder_cell_type = "progenitor",
      cell_type_params = list(
        progenitor = list(
          cell_cycle_length = 1,
          death_per_cell_cycle_prob = death,
          bc_invariant_sites = 0,
          mt_invariant_sites = 0,
          induced_editing_params = state,
          uninduced_editing_params = state
        )
      ),
      uninduced_transition_matrix = list(list(1)),
      induced_transition_matrix = list(list(1))
    )
  )
}

#' Summarise one simulated run against the reference.
run_configuration <- function(label, mean_rate, shape, empirical = NULL) {
  if (is.null(empirical)) {
    cat(sprintf("\n[run] %s: shape %.4f, mean %.5f\n", label, shape, mean_rate))
  } else {
    cat(sprintf("\n[run] %s: %d measured rates, mean %.5f\n", label,
                length(empirical), mean(empirical)))
  }
  run_dir <- file.path(output_dir, label)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  result <- run_gillespie_lineage_pipeline(
    build_params(mean_rate, shape, empirical), params_path = NULL,
    output_dir = run_dir,
    # Barcode only: the mitochondrial path initialises its model regardless of
    # num_mt_genomes being zero, and this run has no mitochondrial component.
    overrides = list(progress = FALSE, seed = seed, modalities = "barcode")
  )
  score_path <- file.path(run_dir, "barcode_binary_score_matrix.csv.gz")
  if (!file.exists(score_path)) {
    stop("No score matrix written for ", label, call. = FALSE)
  }
  score <- data.table::fread(score_path, showProgress = FALSE)
  if (is.character(score[[1L]]) || !is.numeric(score[[1L]])) {
    data.table::set(score, j = 1L, value = NULL)
  }
  matrix_form <- as.matrix(score)
  storage.mode(matrix_form) <- "integer"
  set.seed(seed)
  if (nrow(matrix_form) > cell_count) {
    matrix_form <- matrix_form[sort(sample(nrow(matrix_form), cell_count)), ,
                               drop = FALSE]
  }
  nodes <- data.table::fread(file.path(run_dir, "lineage_nodes.csv.gz"),
                             showProgress = FALSE)
  alive <- nodes[nodes$alive_at_end %in% c(TRUE, "TRUE"), ]
  realised <- stats::median(as.numeric(alive$generation))
  cat(sprintf("[run] %s: %d alive, %d sampled, %d positions; generation median %.0f (IQR %.0f-%.0f), target %d\n",
              label, nrow(alive), nrow(matrix_form), ncol(matrix_form),
              realised, stats::quantile(as.numeric(alive$generation), 0.25),
              stats::quantile(as.numeric(alive$generation), 0.75), divisions))
  if (abs(realised - divisions) > 0.15 * divisions) {
    cat(sprintf("[run] WARNING: realised depth %.0f differs from the %d divisions the rate was calibrated for\n",
                realised, divisions))
  }
  # Match the reference's per-site capture: draw each simulated site a
  # denominator from the reference and score it on that many cells.
  set.seed(seed + 7L)
  site_rate <- vapply(seq_len(ncol(matrix_form)), function(column) {
    depth <- min(sample(reference_observed, 1L), nrow(matrix_form))
    rows <- sample(nrow(matrix_form), depth)
    mean(matrix_form[rows, column] == 1L)
  }, numeric(1))
  list(label = label, mean_rate = if (is.null(empirical)) mean_rate else
         mean(empirical),
       shape = if (is.null(empirical)) shape else NA_real_,
       site_rate = site_rate, edits_per_cell = rowSums(matrix_form == 1L),
       positions = ncol(matrix_form), cells = nrow(matrix_form),
       generation = stats::median(as.numeric(alive$generation)))
}

wanted <- function(label) {
  !length(requested_configurations) || label %in% requested_configurations
}
runs <- list()
if (wanted("shipped_shape0.5")) {
  runs[[length(runs) + 1L]] <- run_configuration("shipped_shape0.5",
                                                 shipped_mean, 0.5)
}
if (wanted("fitted_shape")) {
  runs[[length(runs) + 1L]] <- run_configuration("fitted_shape", fitted_mean,
                                                 fitted_shape)
}
if (!is.null(empirical_rates) && wanted("empirical_rates")) {
  runs[[length(runs) + 1L]] <- run_configuration(
    "empirical_rates", mean(empirical_rates), NA_real_,
    empirical = empirical_rates)
}
if (!length(runs)) stop("No configurations selected.", call. = FALSE)
# ---- Compare -----------------------------------------------------------------
reference_edits <- reference_cells$edited
reference_observed <- reference_cells$observed
cat("\n== Reference (clone 84 subset) ==\n")
cat(sprintf("  sites %d; mean rate %.4f; median %.4f; >=99%% %.3f; <1%% %.3f\n",
            length(observed_rate), mean(observed_rate),
            stats::median(observed_rate), mean(observed_rate >= 0.99),
            mean(observed_rate < 0.01)))
cat(sprintf("  edits per cell: mean %.1f of %.0f observed sites (%.3f)\n",
            mean(reference_edits), stats::median(reference_observed),
            mean(reference_edits / reference_observed)))

summary_rows <- lapply(runs, function(run) {
  rate <- run$site_rate
  data.frame(
    configuration = run$label, shape = run$shape, mean_rate = run$mean_rate,
    positions = run$positions, cells = run$cells, generation = run$generation,
    site_rate_mean = mean(rate), site_rate_median = stats::median(rate),
    fraction_saturated = mean(rate >= 0.99), fraction_low = mean(rate < 0.01),
    edits_per_cell_mean = mean(run$edits_per_cell),
    edited_fraction = mean(run$edits_per_cell) / run$positions,
    quantile_distance = mean(abs(
      stats::quantile(rate, probs = quantile_grid) - observed_quantiles)),
    stringsAsFactors = FALSE
  )
})
summary_table <- do.call(rbind, summary_rows)
reference_row <- data.frame(
  configuration = "reference", shape = NA_real_, mean_rate = NA_real_,
  positions = length(observed_rate), cells = nrow(reference_cells),
  generation = divisions, site_rate_mean = mean(observed_rate),
  site_rate_median = stats::median(observed_rate),
  fraction_saturated = mean(observed_rate >= 0.99),
  fraction_low = mean(observed_rate < 0.01),
  edits_per_cell_mean = mean(reference_edits),
  edited_fraction = mean(reference_edits / reference_observed),
  quantile_distance = 0, stringsAsFactors = FALSE
)
summary_table <- rbind(reference_row, summary_table)
# A partial re-run replaces only its own rows; the rest stand as last measured.
carry_over <- function(existing, ran) {
  path <- file.path(output_dir, existing)
  if (!file.exists(path)) return(NULL)
  previous <- utils::read.csv(path, stringsAsFactors = FALSE)
  previous[!(previous$configuration %in% c("reference", ran)), , drop = FALSE]
}
ran_labels <- vapply(runs, function(run) run$label, character(1))
if (length(requested_configurations)) {
  kept <- carry_over("match_summary.csv", ran_labels)
  if (!is.null(kept) && nrow(kept)) {
    cat(sprintf("[merge] carrying over: %s\n",
                paste(kept$configuration, collapse = ", ")))
    summary_table <- rbind(summary_table, kept[, names(summary_table)])
  }
}
utils::write.csv(summary_table, file.path(output_dir, "match_summary.csv"),
                 row.names = FALSE)
site_rates <- do.call(rbind, c(
  list(data.frame(configuration = "reference", editing_rate = observed_rate,
                  stringsAsFactors = FALSE)),
  lapply(runs, function(run) {
    data.frame(configuration = run$label, editing_rate = run$site_rate,
               stringsAsFactors = FALSE)
  })
))
if (length(requested_configurations)) {
  kept <- carry_over("site_rates.csv", ran_labels)
  if (!is.null(kept) && nrow(kept)) {
    site_rates <- rbind(site_rates, kept[, names(site_rates)])
  }
}
utils::write.csv(site_rates, file.path(output_dir, "site_rates.csv"),
                 row.names = FALSE)

cat("\n== Comparison ==\n")
print(summary_table[, c("configuration", "shape", "site_rate_mean",
                        "site_rate_median", "fraction_saturated",
                        "fraction_low", "edits_per_cell_mean",
                        "quantile_distance")],
      row.names = FALSE, digits = 4)
cat(sprintf("\nWrote:\n  %s\n  %s\n",
            file.path(output_dir, "match_summary.csv"),
            file.path(output_dir, "site_rates.csv")))
