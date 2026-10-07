# Load the ORIGIN simulator by sourcing its files directly.
#
# The core simulator now also builds as an R package under origin/, so the
# normal way to use it is library(origin). This loader exists for the callers
# that cannot use the package yet: code that replaces simulator functions with
# its own versions after loading. That works only when the functions live in
# the global environment, which is what source() gives and a package namespace
# does not. The clique test harness relies on it to inject the WT-CRISPR/FLARE
# recorder, so until those hooks become real extension points this file is the
# supported entry point for it.
#
# Sourcing order does not matter: every file listed here defines functions and
# nothing else, and R resolves calls when they run rather than when they load.
#
# Usage:
#   source("/path/to/remote_mito_clean/load_origin.R")            # core + analysis
#   origin_include_analysis <- FALSE
#   source("/path/to/remote_mito_clean/load_origin.R")            # core only
#
# Set origin_root beforehand to load from a checkout other than this file's own
# directory. Definitions always land in the global environment.

if (!exists("origin_root", inherits = TRUE)) {
  origin_root <- tryCatch({
    # When this file is source()d, one of the active frames carries its path in
    # $ofile. Walk outward from the innermost frame to find it.
    located <- NULL
    for (frame in rev(sys.frames())) {
      if (!is.null(frame$ofile)) {
        located <- frame$ofile
        break
      }
    }
    if (is.null(located)) stop("not sourced from a file")
    dirname(normalizePath(located))
  }, error = function(condition) getwd())
}

origin_core_files <- c(
  "recorder_registry.R",
  "response_curves.R",
  "params_builder.R",
  "prime_editing.R",
  "physicell_lineage.R",
  "physicell_mito.R",
  "ecdna_lineage.R",
  "gillespie_lineage.R",
  "gillespie_pipeline.R"
)
origin_analysis_files <- c(
  "lineage_benchmark.R",
  "engine_comparison.R",
  "physicell_visium.R",
  "scdesign3_helpers.R",
  "physicell_scdesign3.R"
)

origin_source_all <- function(paths) {
  missing_paths <- paths[!file.exists(paths)]
  if (length(missing_paths)) {
    stop(
      "Cannot load ORIGIN; these files are missing:\n  ",
      paste(missing_paths, collapse = "\n  "),
      "\nSet origin_root to the repository checkout that contains origin/R.",
      call. = FALSE
    )
  }
  for (path in paths) source(path, local = FALSE)
  invisible(paths)
}

origin_source_all(file.path(origin_root, "origin", "R", origin_core_files))

if (!exists("origin_include_analysis", inherits = TRUE) ||
    isTRUE(origin_include_analysis)) {
  origin_source_all(file.path(origin_root, "analysis", origin_analysis_files))
}
