# Measuring what the simulator does, as opposed to what its parameters say.
#
# Two relations in this engine are not predictable from the parameter values,
# and both cost real time when assumed rather than measured.
#
# The editing rate is nominal. A target's configured per-division probability is
# not the probability it ends up edited per division: on a BASELINE geometry the
# effective rate came out at 0.019-0.023 against a nominal 0.055, so reaching an
# observed saturation needed about 109 divisions where the arithmetic said 30.
#
# Growth is extremely sensitive to the death probability near the critical
# point, and not linear in it. Measured on one geometry, death 0.45 gives 1.214
# growth per division, 0.50 gives 1.151, 0.52 gives 1.116, 0.55 gives 1.072 and
# 0.58 gives 1.006. Extrapolating 0.52 from the 0.45 point predicts 1.08 and is
# wrong by enough that, compounded over 72 divisions, it is the difference
# between 5,000 cells and 50,000.
#
# Both functions run real simulations, so they cost what a simulation costs.
# Keep the probe small: these measure a relation, not a result.

#' Measure realised population growth against the death probability
#'
#' @export
#' @param params Parameter list to probe. Its death probability is replaced for
#'   each value in `deaths`; everything else is left alone.
#' @param deaths Death probabilities per cell cycle to measure.
#' @param cell_type Cell type whose death probability is varied. Defaults to the
#'   founder cell type.
#' @param sim_length Simulated time per probe. Keep it short; growth per
#'   division is a rate and does not need depth to estimate.
#' @param founders Starting cells. More founders give a steadier estimate near
#'   the critical point, where a small population can go extinct.
#' @param seed Base seed.
#' @return A data frame with one row per death probability: `death`, `alive`,
#'   `generation` (median over surviving cells) and `growth_per_division`.
#'   Probes that go extinct report `alive = 0` and `NA` growth rather than
#'   failing, because extinction is a real answer near the critical point.
origin_growth_response <- function(params,
                                   deaths,
                                   cell_type = NULL,
                                   sim_length = 8,
                                   founders = 200L,
                                   seed = 1L){
  if(!is.numeric(deaths) || !length(deaths) || any(!is.finite(deaths)) ||
     any(deaths < 0) || any(deaths >= 1)){
    stop('deaths must be finite probabilities in [0, 1).')
  }
  if(is.null(cell_type)){
    cell_type <- params$cell_type_dict$founder_cell_type
  }
  if(is.null(params$cell_type_dict$cell_type_params[[cell_type]])){
    stop(sprintf('Cell type %s is not in the parameter list.', cell_type))
  }

  rows <- lapply(deaths, function(death){
    probe <- params
    probe$num_init_cells <- founders
    probe$sim_length <- list(sim_length)
    probe$cell_type_dict$cell_type_params[[cell_type]]$
      death_per_cell_cycle_prob <- death
    output <- file.path(tempdir(), sprintf('origin-growth-%s', death))
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(output, recursive = TRUE), add = TRUE)
    survived <- tryCatch({
      run_gillespie_lineage_pipeline(
        probe, params_path = NULL, output_dir = output,
        overrides = list(progress = FALSE, seed = seed,
                         modalities = 'barcode'))
      TRUE
    }, error = function(condition) FALSE)
    if(!survived){
      return(data.frame(death = death, alive = 0L, generation = NA_real_,
                        growth_per_division = NA_real_))
    }
    nodes <- utils::read.csv(gzfile(file.path(output, 'lineage_nodes.csv.gz')))
    alive <- nodes[nodes$alive_at_end %in% c(TRUE, 'TRUE'), , drop = FALSE]
    if(!nrow(alive)){
      return(data.frame(death = death, alive = 0L, generation = NA_real_,
                        growth_per_division = NA_real_))
    }
    generation <- stats::median(as.numeric(alive$generation))
    data.frame(death = death, alive = nrow(alive), generation = generation,
               growth_per_division = if(generation > 0){
                 (nrow(alive) / founders)^(1 / generation)
               } else NA_real_)
  })
  do.call(rbind, rows)
}

#' Measure realised editing saturation against a configured rate
#'
#' The configured rate is nominal; this reports what fraction of targets are
#' actually edited at the end of a run, and the per-division rate that implies.
#' Inverting a measured saturation through the nominal rate is what produces
#' calibrations that miss by a factor of three.
#'
#' @export
#' @param params Parameter list to probe.
#' @param rates Rates to measure.
#' @param apply_rate `function(params, rate)` returning the parameter list with
#'   that rate set. Which field carries the rate depends on the recorder, so the
#'   caller supplies this rather than the function guessing.
#' @param target_positions Positions to score. Defaults to every position that
#'   is ever edited across the probes, which for an array recorder is its
#'   targets.
#' @param seed Base seed.
#' @return A data frame with one row per rate: `rate`, `molecules`,
#'   `generation`, `target_rate` (fraction of target positions edited) and
#'   `effective_rate_per_division`, the per-division probability that
#'   saturation implies given the realised depth.
origin_rate_response <- function(params,
                                 rates,
                                 apply_rate,
                                 target_positions = NULL,
                                 seed = 1L){
  if(!is.numeric(rates) || !length(rates) || any(!is.finite(rates)) ||
     any(rates < 0)){
    stop('rates must be finite and non-negative.')
  }
  if(!is.function(apply_rate)){
    stop('apply_rate must be a function of (params, rate).')
  }

  rows <- lapply(rates, function(rate){
    probe <- apply_rate(params, rate)
    output <- file.path(tempdir(), sprintf('origin-rate-%s', rate))
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(output, recursive = TRUE), add = TRUE)
    run_gillespie_lineage_pipeline(
      probe, params_path = NULL, output_dir = output,
      overrides = list(progress = FALSE, seed = seed, modalities = 'barcode',
                       write_allele_matrix = TRUE))
    path <- file.path(output, 'barcode_alleles.csv.gz')
    if(!file.exists(path)){
      path <- file.path(output, 'barcode_alleles.csv')
    }
    alleles <- utils::read.csv(if(endsWith(path, '.gz')) gzfile(path) else path,
                               check.names = FALSE)
    if(!is.numeric(alleles[[1L]])){
      alleles <- alleles[, -1L, drop = FALSE]
    }
    matrix_form <- as.matrix(alleles)
    positions <- as.integer(sub('.*_pos_', '', colnames(matrix_form)))
    scored <- if(is.null(target_positions)){
      positions %in% unique(positions[colSums(matrix_form != 0) > 0])
    } else{
      positions %in% target_positions
    }
    if(!any(scored)){
      return(data.frame(rate = rate, molecules = nrow(matrix_form),
                        generation = NA_real_, target_rate = 0,
                        effective_rate_per_division = 0))
    }
    nodes <- utils::read.csv(gzfile(file.path(output, 'lineage_nodes.csv.gz')))
    alive <- nodes[nodes$alive_at_end %in% c(TRUE, 'TRUE'), , drop = FALSE]
    generation <- if(nrow(alive)){
      stats::median(as.numeric(alive$generation))
    } else NA_real_
    saturation <- mean(matrix_form[, scored, drop = FALSE] != 0)
    data.frame(
      rate = rate, molecules = nrow(matrix_form), generation = generation,
      target_rate = saturation,
      effective_rate_per_division = if(is.finite(generation) &&
                                       generation > 0 && saturation < 1){
        1 - (1 - saturation)^(1 / generation)
      } else NA_real_)
  })
  do.call(rbind, rows)
}
