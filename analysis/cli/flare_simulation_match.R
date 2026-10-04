#!/usr/bin/env Rscript

# Does the FLARE/WT-Cas9 recorder reproduce the GSM8791703 HL60 recording?
#
# The real file is molecule-resolved: every read carries a unique 28-mer
# consensus tag, so a row is one integration read out of one cell, with eight
# target calls. Cells are not resolved, which rules out anything needing a
# lineage and leaves the per-molecule structure -- which is where the
# interesting signal is, because inter-target deletions dominate it.
#
# Geometry is matched to the data at EIGHT targets rather than the recorder's
# hardcoded ten, so the distinct-call and span distributions are comparable
# rather than being offset by two extra characters.
#
# Depth is FIXED and the RATE is calibrated, which is the reverse of what was
# tried first. Calibrating depth at the grid's rate is not reachable: the two
# completed sweep points put the simulator's effective rate at 0.019-0.023 per
# division against a nominal 0.055, so 0.88 saturation needs about 109
# divisions, and at the measured growth of 1.245 per division that is 7e10
# cells. max_cells cannot help -- it ends the run rather than subsampling, so
# it truncates the depth it is meant to preserve.
#
# Depth is therefore pinned at 36 divisions, which is ~44 days at a 29h HL60
# doubling time, and the editing rate is scaled until per-target saturation
# matches. The structural features -- distinct calls per molecule, deletion
# span, target-to-target association -- are what the model either reproduces or
# does not, and they are unaffected by which knob was used to match the
# marginal.
#
# Arguments:
#   --real-dir=<dir>     directory of the real-data scan (required)
#   --output-dir=<dir>   destination (required)
#   --sweep=<list>       sim_length values to calibrate over.
#   --integrations=<n>   Default 5.
#   --founders=<n>       Default 5 for the sweep, 20 for the final run.
#   --insertion=<p>      Default 0.005.
#   --deletion=<p>       Default 0.05.
#   --seed=<n>           Default 1.

arguments <- commandArgs(trailingOnly = TRUE)
value_after <- function(prefix, default = NULL) {
  match <- arguments[startsWith(arguments, paste0(prefix, "="))]
  if (!length(match)) return(default)
  sub(paste0("^", prefix, "="), "", match[[length(match)]])
}
real_dir <- value_after("--real-dir")
output_dir <- value_after("--output-dir")
if (is.null(real_dir) || is.null(output_dir)) {
  stop("--real-dir and --output-dir are required.", call. = FALSE)
}
sweep <- as.numeric(strsplit(value_after("--sweep", "1,2,3,4,6"), ",",
                             fixed = TRUE)[[1L]])
depth <- as.numeric(value_after("--depth", "18"))
integrations <- as.integer(value_after("--integrations", "5"))
sweep_founders <- as.integer(value_after("--sweep-founders", "5"))
final_founders <- as.integer(value_after("--founders", "20"))
insertion <- as.numeric(value_after("--insertion", "0.005"))
deletion <- as.numeric(value_after("--deletion", "0.05"))
death <- as.numeric(value_after("--death", "0.45"))
# Deletion resection per side, gamma shape/rate, fitted to the recording. The
# measured lengths are two classes, not one: 33.5% of events are <= 26bp with
# mean 8.7 -- a cut repaired locally, never reaching the next target 26bp away
# -- and 66.5% are longer with mean 130.5, sd 45. Fitting one gamma across
# both gave shape 1.68, which sits between the modes and reproduces neither, so
# the resection component is fitted to the multi-target class alone and the
# short class enters as a mixture component.
extent_shape <- as.numeric(value_after("--extent-shape", "4.1475"))
extent_rate <- as.numeric(value_after("--extent-rate", "0.06407"))
local_probability <- as.numeric(value_after("--local-probability", "0.335"))
local_shape <- as.numeric(value_after("--local-shape", "0.5931"))
local_rate <- as.numeric(value_after("--local-rate", "0.15473"))
# The real array cuts at 63, 89, 115, 141 ... 246, so targets sit 26bp apart
# and span 182bp. The recorder's built-in S:15:20 puts them 21bp apart over
# 147bp, which would make every deletion cover more targets than it should.
array_gap <- as.integer(value_after("--array-gap", "25"))
seed <- as.integer(value_after("--seed", "1"))
targets <- 8L
base_insertion <- insertion
base_deletion <- deletion

source_root <- normalizePath(file.path(dirname(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])), "..", ".."),
  mustWork = TRUE)
source(file.path(source_root, "load_origin.R"))
harness <- "/Users/aaronmck/Desktop/code/clique_2025_12_10/rust_cmd/test_harness/tree_building"
source(file.path(harness, "flare-wt-crispr-recorder.R"))

# Sourcing the recorder is NOT enough to use it. That file only defines
# flare_wt_crispr_* helpers; the dispatch that routes editing into them lives in
# generate-recorder-parameter-grid.R, which wraps two simulator functions in the
# global environment. Without these wrappers model$is_wt_crispr is never set,
# mutate_physicell_barcode_segment takes the native nuclease path, and the
# inter-target deletion machinery -- the whole reason to use this recorder --
# never runs.
base_prepare_recording_model <- prepare_physicell_recording_model
base_mutate_barcode_segment <- mutate_physicell_barcode_segment
base_target_layout <- physicell_baseline_target_layout

prepare_physicell_recording_model <- function(params, cell_type = NULL,
                                              num_integrations = NULL,
                                              founder_label_sites = 0,
                                              params_dir = ".", seed = NULL) {
  model <- base_prepare_recording_model(
    params = params, cell_type = cell_type, num_integrations = num_integrations,
    founder_label_sites = founder_label_sites, params_dir = params_dir,
    seed = seed)
  attach_flare_wt_crispr_model(model, params)
}
mutate_physicell_barcode_segment <- function(profile, duration, rate_set, model,
                                             segment_start = 0) {
  if (isTRUE(model$is_wt_crispr)) {
    return(mutate_flare_wt_crispr_segment(
      profile = profile, duration = duration, rate_set = rate_set,
      model = model, segment_start = segment_start))
  }
  base_mutate_barcode_segment(profile = profile, duration = duration,
                              rate_set = rate_set, model = model,
                              segment_start = segment_start)
}
physicell_baseline_target_layout <- function(model) {
  if (isTRUE(model$is_wt_crispr)) flare_wt_crispr_target_layout(model)
  else base_target_layout(model)
}
# Marks the dispatcher itself. Checking model$is_wt_crispr does NOT prove the
# wiring: attach_flare_wt_crispr_model() sets that flag straight from the
# params whether or not this wrapper was ever installed, so a guard on the flag
# passes on the broken configuration.
attr(mutate_physicell_barcode_segment, "flare_wired") <- TRUE
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

real_rates <- read.delim(file.path(real_dir, "target_rates.tsv"))
# The marginal rate mixes two populations. 11.4% of real molecules are entirely
# unedited because those cells suppressed Cas9, and the rest are essentially
# fully wiped: the distribution of edited targets is 11.4% at zero, 0.45% at
# one to six, and 88.1% at seven or eight. Conditional on being touched the
# per-target rate is 0.9924. Calibrating to the 0.8788 marginal while the model
# has no way to produce the pristine class forces it to under-edit everything
# it does touch, which is what fills the intermediate spans.
target_rate <- as.numeric(value_after("--target-rate",
                                      as.character(mean(real_rates$rate))))
cat(sprintf("[flare] real per-target editing rate: %.4f over %d targets\n",
            target_rate, nrow(real_rates)))

editing_state <- list(
  be_mutations_per_target_per_division = 0,
  nuc_insertions_per_target_per_division = insertion,
  nuc_deletions_per_target_per_division = deletion,
  mt_mutations_per_genome_per_division = 0, num_mt_genomes = 0,
  bc_bg_insertion_prob_per_division = 0, bc_bg_deletion_prob_per_division = 0,
  bc_substitution_model = "JC", bc_sub_model_params = list(0),
  mt_substitution_model = "JC", mt_sub_model_params = list(0),
  mt_bg_insertion_prob_per_division = 0, mt_bg_deletion_prob_per_division = 0)

build_params <- function(sim_length, founders) {
  params <- list(
    num_init_cells = founders, sim_length = list(sim_length),
    random_seed = seed, bc_length = 215L, mito_genome_length = 16600,
    max_bc_ints_per_cell = list(integrations),
    bc_nuc_composition = list(frac_a = 0.25, frac_g = 0.25, frac_c = 0.25,
                              frac_t = 0.25),
    editing_induction = list(timepoint = 0, num_cells = NULL, frac_cells = 1),
    differentiation_induction = list(timepoint = 0, num_cells = NULL,
                                     frac_cells = 1),
    cell_type_dict = list(
      founder_cell_type = "progenitor",
      cell_type_params = list(progenitor = list(
        cell_cycle_length = 1, death_per_cell_cycle_prob = death,
        bc_invariant_sites = 0, mt_invariant_sites = 0,
        induced_editing_params = editing_state,
        uninduced_editing_params = editing_state)),
      uninduced_transition_matrix = list(list(1)),
      induced_transition_matrix = list(list(1))))
  params <- configure_flare_wt_crispr_params(
    params, insertion, deletion,
    deletion_extent = list(
      left = list(shape = extent_shape, rate = extent_rate),
      right = list(shape = extent_shape, rate = extent_rate),
      local = list(probability = local_probability, shape = local_shape,
                   rate = local_rate)))
  # Eight targets, not the recorder's built-in ten. The array keeps its 21bp
  # spacing from position 15, so the last target sits at 162 and the barcode is
  # trimmed to match; the 200bp dropout radius still spans the whole array, as
  # it does at ten targets, so inter-target deletion behaves the same way.
  params$nuclease_targets$num_targets <- targets
  params$nuclease_targets$config <- sprintf("S:15:%d", array_gap)
  # 15 + 7*26 = 197, so 215 leaves the last target room to resect into.
  params$bc_length <- 215L
  params
}

dispatch_checked <- FALSE
reported_spacer <- FALSE
run_once <- function(label, sim_length, founders) {
  run_dir <- file.path(output_dir, label)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  # One-time proof that the recorder is actually in the path. Getting this
  # wrong silently produces a plausible-looking run of the wrong model.
  if (!dispatch_checked) {
    wired <- isTRUE(attr(mutate_physicell_barcode_segment, "flare_wired"))
    routed <- isTRUE(attach_flare_wt_crispr_model(
      list(), build_params(sim_length, founders))$is_wt_crispr)
    if (!wired || !routed) {
      stop("The FLARE recorder is not dispatched (wrapper installed: ", wired,
           ", params routed: ", routed, "); editing would take the native ",
           "nuclease path.", call. = FALSE)
    }
    cat("[flare] dispatch confirmed: wrapper installed and params routed\n")
    dispatch_checked <<- TRUE
  }
  run_gillespie_lineage_pipeline(
    build_params(sim_length, founders), params_path = NULL,
    output_dir = run_dir,
    overrides = list(progress = FALSE, seed = seed, modalities = "barcode",
                     # Allele identity is the whole point here: the binary
                     # matrix cannot tell one deletion from another.
                     write_allele_matrix = TRUE))
  path <- file.path(run_dir, "barcode_alleles.csv.gz")
  if (!file.exists(path)) path <- file.path(run_dir, "barcode_alleles.csv")
  if (!file.exists(path)) stop("No allele matrix for ", label, call. = FALSE)
  alleles <- data.table::fread(path, showProgress = FALSE)
  if (!is.numeric(alleles[[1L]])) data.table::set(alleles, j = 1L, value = NULL)
  # The matrix spans every barcode position, of which only eight are targets;
  # the rest are spacer and always wild type. barcode_target_layout.csv is not
  # the place to look -- it covers base-editing targets and comes back empty for
  # a nuclease recorder -- so positions come from the array config, checked
  # against what actually edited.
  # The layout is populated under the FLARE path, where it is written by
  # flare_wt_crispr_target_layout(); it stays empty on the base path.
  layout_path <- file.path(run_dir, "barcode_target_layout.csv.gz")
  layout_positions <- if (file.exists(layout_path)) {
    frame <- read.csv(gzfile(layout_path), stringsAsFactors = FALSE)
    if (nrow(frame)) sort(unique(frame$position)) else integer()
  } else integer()
  list(matrix = as.matrix(alleles), positions = array_positions(),
       layout_positions = layout_positions)
}

# "S:<start>:<gap>" places the first target at <start> and steps by <gap> + 1,
# which is what puts ten targets at 15, 36, ... 204 in the documented array.
array_positions <- function() {
  15L + (array_gap + 1L) * (seq_len(targets) - 1L)
}

# Molecules are cells x integrations: one row per integration, matching the
# real file where each consensus tag is one integration from one cell.
as_molecules <- function(run) {
  matrix_form <- run$matrix
  positions <- as.integer(sub(".*_pos_", "", colnames(matrix_form)))
  units <- sub("_pos_.*$", "", colnames(matrix_form))
  is_target <- positions %in% run$positions
  if (!any(is_target)) stop("No target positions found in the allele matrix.",
                            call. = FALSE)
  # Spacer positions SHOULD change now: an interval deletion removes the
  # sequence between cuts, which is the whole point of the ported model. What
  # still has to hold is that the assumed cut sites are the real ones, checked
  # against the layout the recorder itself writes.
  edited_positions <- unique(positions[colSums(matrix_form != 0) > 0])
  spacer <- setdiff(edited_positions, run$positions)
  if (!is.null(run$layout_positions) && length(run$layout_positions)) {
    if (!identical(sort(as.integer(run$layout_positions)),
                   sort(as.integer(run$positions)))) {
      stop("Cut sites disagree: layout says ",
           paste(run$layout_positions, collapse = ", "), " but the array config ",
           "implies ", paste(run$positions, collapse = ", "), call. = FALSE)
    }
  }
  if (!reported_spacer) {
    cat(sprintf("[flare] %d of %d spacer positions edited by resection\n",
                length(spacer), max(positions) - length(run$positions)))
    reported_spacer <<- TRUE
  }
  do.call(rbind, lapply(unique(units), function(unit) {
    block <- matrix_form[, units == unit & is_target, drop = FALSE]
    block <- block[, order(positions[units == unit & is_target]), drop = FALSE]
    colnames(block) <- paste0("target", seq_len(ncol(block)))
    block
  }))
}

cat(sprintf("\n[flare] calibrating rate at a fixed depth of %.0f divisions\n",
            2 * depth))
calibration <- do.call(rbind, lapply(sweep, function(multiplier) {
  insertion <<- base_insertion * multiplier
  deletion <<- base_deletion * multiplier
  editing_state$nuc_insertions_per_target_per_division <<- insertion
  editing_state$nuc_deletions_per_target_per_division <<- deletion
  molecules <- as_molecules(run_once(sprintf("sweep_x%g", multiplier), depth,
                                     sweep_founders))
  rate <- mean(molecules != 0)
  cat(sprintf("  rate x%-4g (%.4f ins / %.4f del) -> %6d molecules, per-target rate %.4f\n",
              multiplier, insertion, deletion, nrow(molecules), rate))
  data.frame(multiplier = multiplier, insertion = insertion,
             deletion = deletion, molecules = nrow(molecules),
             target_rate = rate, stringsAsFactors = FALSE)
}))
calibration$distance <- abs(calibration$target_rate - target_rate)
utils::write.csv(calibration, file.path(output_dir, "depth_calibration.csv"),
                 row.names = FALSE)

# Interpolate rather than settling for the nearest swept point: saturation is
# monotone in depth, so the crossing is well defined between the two points
# that bracket the observed rate.
bracket <- calibration[order(calibration$target_rate), ]
chosen <- if (any(bracket$target_rate >= target_rate) &&
               any(bracket$target_rate <= target_rate)) {
  stats::approx(bracket$target_rate, bracket$multiplier, xout = target_rate)$y
} else {
  calibration$multiplier[which.min(calibration$distance)]
}
insertion <- base_insertion * chosen
deletion <- base_deletion * chosen
editing_state$nuc_insertions_per_target_per_division <- insertion
editing_state$nuc_deletions_per_target_per_division <- deletion
cat(sprintf("\n[flare] chosen rate multiplier %.2f -> %.4f insertion / %.4f deletion\n",
            chosen, insertion, deletion))
cat(sprintf("[flare] that is %.1fx the parameter grid's selected tier\n", chosen))

molecules <- as_molecules(run_once("calibrated", depth, final_founders))
cat(sprintf("[flare] calibrated run: %d molecules x %d targets, rate %.4f\n",
            nrow(molecules), ncol(molecules), mean(molecules != 0)))
saveRDS(list(molecules = molecules, sim_length = depth, multiplier = chosen,
             insertion = insertion, deletion = deletion),
        file.path(output_dir, "calibrated_molecules.rds"))
utils::write.csv(data.frame(sim_length = depth, rate_multiplier = chosen,
                            target_rate_observed = target_rate,
                            target_rate_simulated = mean(molecules != 0),
                            molecules = nrow(molecules), targets = ncol(molecules),
                            insertion = insertion, deletion = deletion),
                 file.path(output_dir, "calibrated_settings.csv"), row.names = FALSE)
cat(sprintf("\nWrote %s\n", file.path(output_dir, "calibrated_molecules.rds")))
