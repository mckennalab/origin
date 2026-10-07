# Lightweight population-outcome comparison for the exact Gillespie and
# discrete time-step engines. Recording overlays are intentionally excluded so
# the benchmark isolates temporal discretization effects.
#
# This file is the LIBRARY half of the engine-comparison pair: it only defines
# functions and runs nothing when sourced. The command-line ENTRY POINT is
# `compare_simulation_engines.R`, which sources this file and writes the
# benchmark's CSV tables and PNG figure. `gillespie_lineage.R` must be sourced
# first (for `gillespie_scalar()`, `gillespie_cell_type_rates()`,
# `gillespie_transition_matrices()`, `gillespie_induction_spec()`,
# `gillespie_induction_selection()`, and `simulate_gillespie_population()`).

#' Simulate one population replicate with the discrete time-step engine
#'
#' Advances the whole population on a fixed grid of width `time_step`. Within
#' each interval every live cell independently fires an event with probability
#' `1 - exp(-total_rate * interval)` (competing hazards), and a firing cell
#' either divides into two daughters drawn from the appropriate cell-type
#' transition matrix or dies. Only aggregate outcomes are tracked: no lineage
#' tree, recording overlay, or mutation profile is produced.
#'
#' @param params Parsed native remote_mito parameter list (from the JSON
#'   file). Must supply `cell_type_dict$cell_type_params`,
#'   `cell_type_dict$founder_cell_type`, the uninduced/induced transition
#'   matrices, the optional `differentiation_induction` specification, and
#'   `num_init_cells` (validated to be positive and no larger than
#'   `max_cells`).
#' @param end_time Simulated duration; one finite non-negative number. The
#'   final interval is shortened so the grid lands exactly on `end_time`.
#' @param time_step Width of each discrete interval; one finite positive
#'   number.
#' @param seed Integer RNG seed, applied with `set.seed()` before simulating.
#' @param max_cells Integer safety cap on the live population. When a step
#'   ends with `max_cells` or more cells, any excess is removed by uniform
#'   subsampling and the run stops with `stop_reason = 'max_cells'`.
#' @return A named list with `engine` (`'time_step'`), `final_cells`,
#'   `divisions`, `deaths` (cumulative event counts), `extinct` (logical),
#'   `stop_reason` (`'end_time'` or `'max_cells'`), `simulated_end_time`,
#'   `cell_type_counts` (integer vector named by every configured cell type),
#'   and `seed`.
#' @section Side effects: Calls `set.seed(seed)`, changing the global RNG
#'   state.
simulate_timestep_population_outcome <- function(params,
                                                  end_time,
                                                  time_step,
                                                  seed,
                                                  max_cells = 1000000L){
  end_time <- gillespie_scalar(end_time, 'end_time')
  time_step <- gillespie_scalar(time_step, 'time_step')
  seed <- gillespie_scalar(seed, 'seed', 'integer')
  max_cells <- gillespie_scalar(max_cells, 'max_cells', 'integer')
  if(!is.finite(end_time) || end_time < 0){
    stop('end_time must be finite and non-negative.')
  }
  if(!is.finite(time_step) || time_step <= 0){
    stop('time_step must be finite and positive.')
  }
  if(max_cells < 1L){
    stop('max_cells must be positive.')
  }
  set.seed(seed)

  cell_types <- names(params$cell_type_dict$cell_type_params)
  founder_type <- as.character(params$cell_type_dict$founder_cell_type)
  if(length(cell_types) == 0L || !(founder_type %in% cell_types)){
    stop('The parameter file has no valid founder cell type.')
  }
  rates <- gillespie_cell_type_rates(params, cell_types)
  transitions <- gillespie_transition_matrices(params, cell_types)
  differentiation <- gillespie_induction_spec(
    params,
    'differentiation_induction'
  )
  num_founders <- gillespie_scalar(
    params$num_init_cells,
    'num_init_cells',
    'integer'
  )
  if(num_founders < 1L || num_founders > max_cells){
    stop('num_init_cells must be positive and no larger than max_cells.')
  }

  active_types <- rep(match(founder_type, cell_types), num_founders)
  differentiation_induced <- rep(FALSE, num_founders)
  differentiation_assigned <- FALSE
  divisions <- 0L
  deaths <- 0L
  current_time <- 0
  stop_reason <- 'end_time'

  #' Internal: flag cells for differentiation induction at the configured time
  #'
  #' Runs at most once: at the first grid boundary at or after the configured
  #' `differentiation_induction` timepoint it selects live cells via
  #' `gillespie_induction_selection()` and sets their entries in the enclosing
  #' `differentiation_induced` vector to `TRUE` (via `<<-`). Does nothing when
  #' the induction timepoint is not finite or has not yet been reached.
  #'
  #' @return `invisible(NULL)`; called for its effect on the enclosing
  #'   `differentiation_induced` and `differentiation_assigned` state.
  assign_differentiation <- function(){
    if(differentiation_assigned ||
       !is.finite(differentiation$timepoint) ||
       current_time + 1e-12 < differentiation$timepoint){
      return(invisible(NULL))
    }
    selected <- gillespie_induction_selection(
      seq_along(active_types),
      differentiation$num_cells,
      differentiation$frac_cells
    )
    if(length(selected) > 0L){
      differentiation_induced[selected] <<- TRUE
    }
    differentiation_assigned <<- TRUE
    invisible(NULL)
  }

  assign_differentiation()
  while(current_time < end_time - 1e-12 && length(active_types) > 0L){
    interval <- min(time_step, end_time - current_time)
    cell_total_rates <- unname(rates$total[active_types])
    event_probabilities <- -expm1(-cell_total_rates * interval)
    has_event <- stats::runif(length(active_types)) < event_probabilities
    event_indices <- which(has_event)

    dividing_indices <- dying_indices <- integer()
    if(length(event_indices) > 0L){
      division_probabilities <- ifelse(
        cell_total_rates[event_indices] > 0,
        unname(rates$division[active_types[event_indices]]) /
          cell_total_rates[event_indices],
        0
      )
      is_division <- stats::runif(length(event_indices)) <
        division_probabilities
      dividing_indices <- event_indices[is_division]
      dying_indices <- event_indices[!is_division]
    }

    unchanged_indices <- if(length(event_indices) == 0L){
      seq_along(active_types)
    } else{
      setdiff(seq_along(active_types), event_indices)
    }
    next_types <- active_types[unchanged_indices]
    next_differentiation <- differentiation_induced[unchanged_indices]

    if(length(dividing_indices) > 0L){
      daughter_types <- integer(2L * length(dividing_indices))
      daughter_differentiation <- rep(
        differentiation_induced[dividing_indices],
        each = 2L
      )
      output_offset <- 0L
      for(parent_index in dividing_indices){
        parent_type <- active_types[parent_index]
        transition_matrix <- if(differentiation_induced[parent_index]){
          transitions$induced
        } else{
          transitions$uninduced
        }
        daughter_types[output_offset + 1:2] <- sample.int(
          length(cell_types),
          2L,
          replace = TRUE,
          prob = transition_matrix[parent_type, ]
        )
        output_offset <- output_offset + 2L
      }
      next_types <- c(next_types, daughter_types)
      next_differentiation <- c(
        next_differentiation,
        daughter_differentiation
      )
    }

    divisions <- divisions + length(dividing_indices)
    deaths <- deaths + length(dying_indices)
    active_types <- next_types
    differentiation_induced <- next_differentiation
    current_time <- current_time + interval

    if(length(active_types) >= max_cells){
      if(length(active_types) > max_cells){
        retained <- sample.int(length(active_types), max_cells, replace = FALSE)
        active_types <- active_types[retained]
        differentiation_induced <- differentiation_induced[retained]
      }
      stop_reason <- 'max_cells'
      break
    }
    assign_differentiation()
  }

  type_counts <- tabulate(active_types, nbins = length(cell_types))
  names(type_counts) <- cell_types
  list(
    engine = 'time_step',
    final_cells = length(active_types),
    divisions = divisions,
    deaths = deaths,
    extinct = length(active_types) == 0L,
    stop_reason = stop_reason,
    simulated_end_time = current_time,
    cell_type_counts = type_counts,
    seed = seed
  )
}

#' Reduce an exact Gillespie population run to comparison-compatible outcomes
#'
#' Collapses a full `simulate_gillespie_population()` result into the same
#' named-list shape produced by `simulate_timestep_population_outcome()`, so
#' the two engines can be tabulated row-for-row by
#' `compare_population_engines()`.
#'
#' @param simulation Result list from `simulate_gillespie_population()`; the
#'   fields used are `terminal_nodes` (its `cell_type` column), `counts` (its
#'   `divisions` and `deaths` entries), `stop_reason`, `end_time`, `seed`, and
#'   `rates$division`, whose names define the cell-type order of the counts.
#' @return A named list with `engine` (`'continuous_time'`), `final_cells`,
#'   `divisions`, `deaths`, `extinct` (logical), `stop_reason`,
#'   `simulated_end_time`, `cell_type_counts` (named integer vector), and
#'   `seed`.
summarize_gillespie_population_outcome <- function(simulation){
  cell_types <- names(simulation$rates$division)
  terminal_types <- factor(
    simulation$terminal_nodes$cell_type,
    levels = cell_types
  )
  type_counts <- tabulate(as.integer(terminal_types), nbins = length(cell_types))
  names(type_counts) <- cell_types
  list(
    engine = 'continuous_time',
    final_cells = nrow(simulation$terminal_nodes),
    divisions = unname(simulation$counts[['divisions']]),
    deaths = unname(simulation$counts[['deaths']]),
    extinct = nrow(simulation$terminal_nodes) == 0L,
    stop_reason = simulation$stop_reason,
    simulated_end_time = simulation$end_time,
    cell_type_counts = type_counts,
    seed = simulation$seed
  )
}

#' Run a paired benchmark of the continuous-time and time-step engines
#'
#' For each replicate the same seed drives one exact Gillespie run
#' (`simulate_gillespie_population()`) and one discrete time-step run
#' (`simulate_timestep_population_outcome()`), pairing the replicates across
#' engines. Wall-clock time is recorded separately for each run.
#'
#' @param params Parsed native remote_mito parameter list, forwarded to both
#'   engines. Also supplies the defaults below.
#' @param replicates Integer number of replicate pairs; must be at least two
#'   so the downstream summaries can compute standard deviations.
#' @param end_time Simulated duration for both engines, or `NULL` to use the
#'   maximum value in `params$sim_length`.
#' @param time_step Grid width for the time-step engine, or `NULL` to use the
#'   numeric `params$time_inc` when it is one finite positive number, and
#'   otherwise one twentieth of the shortest configured cell-cycle length.
#' @param seed Integer seed for the first replicate; replicate `i` uses
#'   `seed + i - 1` for both engines.
#' @param max_cells Integer per-replicate population cap passed to both
#'   engines.
#' @param show_progress Logical; when `TRUE`, progress lines are printed to
#'   the console as replicates complete.
#' @param progress_updates Integer number of progress lines to aim for across
#'   the run; must be positive.
#' @return A named list with `replicate_summary` (data frame, one row per
#'   replicate x engine with `replicate`, `seed`, `engine`, `final_cells`,
#'   `divisions`, `deaths`, `extinct`, `reached_max_cells`,
#'   `simulated_end_time`, and `elapsed_seconds`), `cell_type_summary`
#'   (long-form data frame, one row per replicate x engine x cell type with
#'   `cell_count` and `cell_fraction`; fractions are zero when a replicate
#'   went extinct), and the resolved `end_time`, `time_step`, `replicates`,
#'   `seed`, and `max_cells`.
#' @section Side effects: Each engine call reseeds the global RNG via
#'   `set.seed()`; progress is printed with `cat()` when `show_progress` is
#'   `TRUE`.
compare_population_engines <- function(params,
                                       replicates = 50L,
                                       end_time = NULL,
                                       time_step = NULL,
                                       seed = 1L,
                                       max_cells = 1000000L,
                                       show_progress = TRUE,
                                       progress_updates = 10L){
  replicates <- gillespie_scalar(replicates, 'replicates', 'integer')
  seed <- gillespie_scalar(seed, 'seed', 'integer')
  progress_updates <- gillespie_scalar(
    progress_updates,
    'progress_updates',
    'integer'
  )
  if(replicates < 2L){
    stop('replicates must be at least two.')
  }
  if(progress_updates < 1L){
    stop('progress_updates must be positive.')
  }
  if(is.null(end_time)){
    end_time <- max(as.numeric(unlist(params$sim_length, use.names = FALSE)))
  }
  if(is.null(time_step)){
    configured_step <- suppressWarnings(as.numeric(params$time_inc))
    time_step <- if(length(configured_step) == 1L &&
                    is.finite(configured_step) && configured_step > 0){
      configured_step
    } else{
      min(vapply(
        params$cell_type_dict$cell_type_params,
        function(type_params) as.numeric(type_params$cell_cycle_length),
        numeric(1)
      ), na.rm = TRUE) / 20
    }
  }
  end_time <- gillespie_scalar(end_time, 'end_time')
  time_step <- gillespie_scalar(time_step, 'time_step')
  max_cells <- gillespie_scalar(max_cells, 'max_cells', 'integer')

  cell_types <- names(params$cell_type_dict$cell_type_params)
  replicate_rows <- vector('list', 2L * replicates)
  type_rows <- vector('list', 2L * replicates)
  next_progress <- 1L
  progress_points <- unique(pmax(
    1L,
    ceiling(seq(1, replicates, length.out = progress_updates))
  ))

  for(replicate_index in seq_len(replicates)){
    replicate_seed <- seed + replicate_index - 1L
    continuous_start <- unname(proc.time()[['elapsed']])
    continuous_simulation <- suppressWarnings(simulate_gillespie_population(
      params,
      end_time = end_time,
      seed = replicate_seed,
      max_cells = max_cells,
      show_progress = FALSE
    ))
    continuous <- summarize_gillespie_population_outcome(
      continuous_simulation
    )
    continuous_elapsed <- unname(proc.time()[['elapsed']]) - continuous_start

    timestep_start <- unname(proc.time()[['elapsed']])
    timestep <- simulate_timestep_population_outcome(
      params,
      end_time = end_time,
      time_step = time_step,
      seed = replicate_seed,
      max_cells = max_cells
    )
    timestep_elapsed <- unname(proc.time()[['elapsed']]) - timestep_start

    outcomes <- list(continuous, timestep)
    elapsed <- c(continuous_elapsed, timestep_elapsed)
    for(engine_index in 1:2){
      outcome <- outcomes[[engine_index]]
      output_index <- 2L * (replicate_index - 1L) + engine_index
      replicate_rows[[output_index]] <- data.frame(
        replicate = replicate_index,
        seed = replicate_seed,
        engine = outcome$engine,
        final_cells = outcome$final_cells,
        divisions = outcome$divisions,
        deaths = outcome$deaths,
        extinct = outcome$extinct,
        reached_max_cells = outcome$stop_reason == 'max_cells',
        simulated_end_time = outcome$simulated_end_time,
        elapsed_seconds = elapsed[engine_index],
        stringsAsFactors = FALSE
      )
      type_rows[[output_index]] <- data.frame(
        replicate = replicate_index,
        seed = replicate_seed,
        engine = outcome$engine,
        cell_type = cell_types,
        cell_count = as.integer(outcome$cell_type_counts[cell_types]),
        cell_fraction = if(outcome$final_cells > 0){
          as.numeric(outcome$cell_type_counts[cell_types]) /
            outcome$final_cells
        } else{
          rep(0, length(cell_types))
        },
        stringsAsFactors = FALSE
      )
    }

    if(isTRUE(show_progress) && next_progress <= length(progress_points) &&
       replicate_index >= progress_points[next_progress]){
      cat(sprintf(
        '[engine comparison] %d / %d replicates complete (%.0f%%)\n',
        replicate_index,
        replicates,
        100 * replicate_index / replicates
      ))
      flush.console()
      next_progress <- next_progress + 1L
    }
  }

  list(
    replicate_summary = do.call(rbind, replicate_rows),
    cell_type_summary = do.call(rbind, type_rows),
    end_time = end_time,
    time_step = time_step,
    replicates = replicates,
    seed = seed,
    max_cells = max_cells
  )
}

#' Summarize per-metric agreement between the two engines
#'
#' Averages each outcome metric across replicates within each engine and
#' reports the between-engine gap. Because the benchmark asks whether the
#' discrete engine reproduces the exact one, the quality criterion for the
#' outcome metrics is agreement: `time_step_minus_continuous` and
#' `relative_difference_percent` closer to zero are better, in either
#' direction. `elapsed_seconds` is the exception where a direction is
#' meaningful on its own: a lower mean is better (that engine ran faster).
#'
#' @details The six metrics summarized, one row each, are: `final_cells`
#'   (live cells when the run stopped), `divisions` and `deaths` (cumulative
#'   event counts per replicate), `extinct` and `reached_max_cells` (0/1
#'   indicators, so their means are the fraction of replicates that went
#'   extinct or hit the `max_cells` cap), and `elapsed_seconds` (wall-clock
#'   runtime per replicate). `relative_difference_percent` is
#'   `100 * (time_step_mean - continuous_time_mean) / continuous_time_mean`,
#'   and is `NA` when the continuous-time mean is zero.
#'
#' @param comparison Result list from `compare_population_engines()`; only
#'   its `replicate_summary` data frame is used.
#' @return A data frame with one row per metric and columns `metric`,
#'   `continuous_time_mean`, `continuous_time_sd`, `time_step_mean`,
#'   `time_step_sd`, `time_step_minus_continuous`, and
#'   `relative_difference_percent`.
summarize_engine_comparison <- function(comparison){
  replicate_summary <- comparison$replicate_summary
  metrics <- c(
    'final_cells', 'divisions', 'deaths', 'extinct',
    'reached_max_cells', 'elapsed_seconds'
  )
  rows <- lapply(metrics, function(metric){
    continuous <- as.numeric(replicate_summary[
      replicate_summary$engine == 'continuous_time',
      metric
    ])
    timestep <- as.numeric(replicate_summary[
      replicate_summary$engine == 'time_step',
      metric
    ])
    continuous_mean <- mean(continuous)
    timestep_mean <- mean(timestep)
    difference <- timestep_mean - continuous_mean
    data.frame(
      metric = metric,
      continuous_time_mean = continuous_mean,
      continuous_time_sd = stats::sd(continuous),
      time_step_mean = timestep_mean,
      time_step_sd = stats::sd(timestep),
      time_step_minus_continuous = difference,
      relative_difference_percent = if(abs(continuous_mean) > 0){
        100 * difference / continuous_mean
      } else{
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Summarize final cell-type composition per engine
#'
#' Averages the per-replicate final counts and fractions of each cell type
#' within each engine. Neither a higher nor a lower value is better on its
#' own; the desired result is close agreement between the two engines' means
#' for the same cell type.
#'
#' @param comparison Result list from `compare_population_engines()`; only
#'   its `cell_type_summary` data frame is used.
#' @return A data frame ordered by cell type then engine, with columns
#'   `engine`, `cell_type`, `mean_cell_count`, `sd_cell_count`,
#'   `mean_cell_fraction`, and `sd_cell_fraction` (means and standard
#'   deviations across replicates).
summarize_engine_cell_types <- function(comparison){
  values <- comparison$cell_type_summary
  groups <- split(values, interaction(
    values$engine,
    values$cell_type,
    drop = TRUE
  ))
  rows <- lapply(groups, function(group){
    data.frame(
      engine = group$engine[1],
      cell_type = group$cell_type[1],
      mean_cell_count = mean(group$cell_count),
      sd_cell_count = stats::sd(group$cell_count),
      mean_cell_fraction = mean(group$cell_fraction),
      sd_cell_fraction = stats::sd(group$cell_fraction),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[order(result$cell_type, result$engine), , drop = FALSE]
}

#' Write the two-panel engine-comparison figure as a PNG
#'
#' Left panel: boxplots of final live cells per replicate, one box per
#' engine. Right panel: grouped bars of the mean divisions and deaths per
#' replicate for each engine.
#'
#' @param comparison Result list from `compare_population_engines()`; its
#'   `replicate_summary`, `time_step`, and `end_time` supply the data and
#'   panel titles.
#' @param path Destination PNG file path.
#' @return `path`, invisibly.
#' @section Side effects: Opens and closes a PNG graphics device at `path`
#'   (1400 x 650 pixels at 140 dpi), overwriting any existing file there.
write_engine_comparison_plot <- function(comparison, path){
  replicate_summary <- comparison$replicate_summary
  grDevices::png(path, width = 1400, height = 650, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
  engine_labels <- c(
    continuous_time = 'Continuous time',
    time_step = sprintf('Time step (dt=%g)', comparison$time_step)
  )
  engine_order <- c('continuous_time', 'time_step')
  population_values <- split(
    replicate_summary$final_cells,
    factor(replicate_summary$engine, levels = engine_order)
  )
  graphics::boxplot(
    population_values,
    names = unname(engine_labels[engine_order]),
    ylab = 'Final live cells',
    main = sprintf('Population at time %g', comparison$end_time),
    col = c('#4477AA', '#EE9944')
  )
  event_values <- rbind(
    tapply(
      replicate_summary$divisions,
      factor(replicate_summary$engine, levels = engine_order),
      mean
    ),
    tapply(
      replicate_summary$deaths,
      factor(replicate_summary$engine, levels = engine_order),
      mean
    )
  )
  graphics::barplot(
    event_values,
    beside = TRUE,
    names.arg = unname(engine_labels[engine_order]),
    ylab = 'Mean events per replicate',
    main = 'Division and death outcomes',
    col = c('#66AA55', '#CC6677')
  )
  graphics::legend(
    'topright',
    legend = c('Divisions', 'Deaths'),
    fill = c('#66AA55', '#CC6677'),
    bty = 'n'
  )
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  invisible(path)
}
