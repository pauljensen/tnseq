
#' Parallel processing features.
#' 
#' Tn-seq computations can be carried out in parallel on multiple cores.  
#' Parallelization is accomplished using the "parallel" package.  The following
#' functions enable and configure parallel processing.
#' 
#' @name parallel
NULL

#' @describeIn parallel Set the number of cores.  Set to 1 to turn off
#' parallelization.
#' 
#' @export
use_multiple_cores <- function(ncores) {
  options(mc.cores=ncores)
}

#' @describeIn parallel Alias for use_multiple_cores(1).
#' 
#' @export
use_single_core <- function() {
  use_multiple_cores(1)
}

#' @describeIn parallel Return the number of cores to be used.
#' 
#' @export
get_cores <- function() {
  getOption("mc.cores")
}

configuration_defaults <- list(
  QUALITY_CUTOFF = NA,
  MATCH_TN_SEQUENCE = FALSE,
  MAX_CYCLES = NA,
  DUMP_MATCH_FAILS = FALSE
)

configure <- function(...) {
  opts <- list(...)
  
}

cfg <- function(opt) {
  getOption(deparse(substitute(opt)))
}
