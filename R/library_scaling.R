#' Perform library scaling on a hicexp object
#' 
#' @param hicexp A hicexp object.
#' @details This function will perform library scaling 
#'     on a hicexp object. Scaling is performed 
#'     separately for each chromosome. This is an 
#'     alternative normalization method to the
#'     cyclic loess and fastlo methods also
#'     provided in multiHiCcompare. Use this normalization
#'     method if for some reason you do not want to remove
#'     trends in the data and only want to normalize based
#'     on library size.
#' @importFrom dplyr %>%
#' @importFrom data.table rbindlist
#' @export
#' @return A hicexp object.
#' @examples 
#' data("hicexp2")
#' hicexp2 <- hic_scale(hicexp2)
#' 
hic_scale <- function(hicexp) {
  # check if data already normalized
  if (hicexp@normalized) {
    stop("Data has already been normalized.")
  }
  # extract hic_table
  tab <- hicexp@hic_table
  # split up by chromosome
  tab <- split(tab, f = tab$chr)
  # scale each chr
  new_tab <- lapply(tab, .scale_chr)
  # recombine
  new_tab <- data.table::rbindlist(new_tab)
  hicexp@hic_table <- new_tab
  hicexp@normalized <- TRUE
  return(hicexp)
}



# internal function for scaling each chr
.scale_chr <- function(tab) {
  # extract IF matrix
  IFs <- tab[, 5:ncol(tab)] %>% as.matrix()
  # get min library size
  min.lib <- colSums(IFs) %>% min()
  # perform scaling based on min lib size
  scaled_IFs <- apply(IFs, 2, function(x) {
    scale.factor <- sum(x) / min.lib
    new.IF <- x / scale.factor
    return(new.IF)
  })
  # recombine table
  new_table <- cbind(tab[, 1:4], scaled_IFs)
  return(new_table)
}

