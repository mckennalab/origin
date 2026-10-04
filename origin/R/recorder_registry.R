# Recorder registration.
#
# A recorder is a set of hooks the simulator calls at four points: after a
# recording model is prepared, when a branch segment mutates, when the target
# layout is written, and when the character matrix tree building reads is built.
#
# These used to be installed by assigning over the simulator's own functions in
# the global environment, which has one failure mode that costs days. Setting
# the parameters looks sufficient -- it even sets the model flag the dispatch
# checks -- so a run configured that way starts, finishes, and writes plausible
# output while silently taking the default path. A guard on the flag passes too,
# because the flag comes from the parameters rather than from the wiring. The
# registry makes the wiring itself the thing that is registered and the thing
# that is checked, so a recorder is either dispatched or absent, never half
# installed.

origin_recorder_registry <- new.env(parent = emptyenv())

#' Register a lineage recorder
#'
#' @export
#' @param name Recorder name, used in diagnostics and to unregister.
#' @param applies Predicate taking a prepared model and returning `TRUE` when
#'   this recorder should handle it.
#' @param prepare Optional `function(model, params)` run after the model is
#'   built; it returns the model, and is where a recorder attaches its own
#'   state. Every registered recorder's `prepare` runs, in registration order.
#' @param mutate Optional replacement for `mutate_physicell_barcode_segment()`,
#'   with the same arguments.
#' @param layout Optional replacement for `physicell_baseline_target_layout()`.
#' @param character_matrix Optional `function(raw_alleles, binary_scores, model,
#'   integrations)` building the matrix tree building reads.
#' @return The name, invisibly.
register_origin_recorder <- function(name,
                                     applies,
                                     prepare = NULL,
                                     mutate = NULL,
                                     layout = NULL,
                                     character_matrix = NULL){
  if(!is.character(name) || length(name) != 1L || !nzchar(name)){
    stop('A recorder name must be one non-empty string.')
  }
  if(!is.function(applies)){
    stop('A recorder must supply an `applies` predicate.')
  }
  for(hook in list(prepare, mutate, layout, character_matrix)){
    if(!is.null(hook) && !is.function(hook)){
      stop('Recorder hooks must be functions or NULL.')
    }
  }
  if(is.null(mutate) && is.null(layout) && is.null(character_matrix) &&
     is.null(prepare)){
    stop('A recorder must supply at least one hook.')
  }
  assign(name, list(name = name, applies = applies, prepare = prepare,
                    mutate = mutate, layout = layout,
                    character_matrix = character_matrix),
         envir = origin_recorder_registry)
  invisible(name)
}

#' Remove a registered recorder
#'
#' @export
#' @param name Recorder name. Removing one that is not registered is silent, so
#'   teardown in a test does not depend on registration having succeeded.
#' @return `TRUE` when a recorder was removed, `FALSE` otherwise, invisibly.
unregister_origin_recorder <- function(name){
  if(exists(name, envir = origin_recorder_registry, inherits = FALSE)){
    rm(list = name, envir = origin_recorder_registry)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

#' Names of the registered recorders
#'
#' @export
#' @return A character vector, in registration order.
origin_recorders <- function(){
  ls(origin_recorder_registry, sorted = FALSE)
}

#' Report which recorder a prepared model dispatches to
#'
#' The question a guard should ask. It reflects what is registered, so it
#' cannot report a recorder that was configured but never wired.
#'
#' @export
#' @param model Prepared recording model.
#' @return The recorder name, or `NULL` when none applies.
origin_active_recorder <- function(model){
  for(name in origin_recorders()){
    recorder <- get(name, envir = origin_recorder_registry, inherits = FALSE)
    applies <- tryCatch(isTRUE(recorder$applies(model)),
                        error = function(condition) FALSE)
    if(applies){
      return(name)
    }
  }
  NULL
}

#' Fetch one hook of the recorder a model dispatches to
#'
#' @export
#' @param model Prepared recording model.
#' @param hook One of `prepare`, `mutate`, `layout`, `character_matrix`.
#' @return The hook function, or `NULL` when no recorder applies or the one that
#'   does leaves that hook unset.
origin_recorder_hook <- function(model, hook){
  hook <- match.arg(hook, c('prepare', 'mutate', 'layout', 'character_matrix'))
  name <- origin_active_recorder(model)
  if(is.null(name)){
    return(NULL)
  }
  get(name, envir = origin_recorder_registry, inherits = FALSE)[[hook]]
}

#' Run every registered recorder's prepare hook over a model
#'
#' Called by `prepare_physicell_recording_model()`. All prepare hooks run, in
#' registration order, because a recorder decides from the parameters whether it
#' applies and signals that by what it attaches; `applies()` is only consulted
#' afterwards, once the model carries that state.
#'
#' @export
#' @param model Prepared recording model.
#' @param params Parameter list the model was built from.
#' @return The model, after every prepare hook has seen it.
origin_apply_recorder_prepare <- function(model, params){
  for(name in origin_recorders()){
    recorder <- get(name, envir = origin_recorder_registry, inherits = FALSE)
    if(is.null(recorder$prepare)){
      next
    }
    model <- recorder$prepare(model, params)
    if(!is.list(model)){
      stop(sprintf('Recorder %s prepare hook did not return a model.', name))
    }
  }
  model
}
