#' Function to apply either biocParallel or standard lapply
#' 
#' @importFrom pbapply pblapply
#' @param parallel Logical, should parallel processing be used?
#' @param x The main list object which the function will be applied to.
#' @param FUN The function to be applied.
#' @param ... Additional arguments for bplapply or lapply.
#' @return results of lapply or bplapply


smartApply <- function(parallel, x, FUN, ...) {
  if (parallel) {
    BiocParallel::bplapply(x, FUN = FUN, ... = ...)
  } else {
    # lapply(x, FUN = FUN, ... = ...)
    pbapply::pblapply(x, FUN = FUN, ... = ...)
  }
}