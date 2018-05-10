#' Perform filtering on a Hi-C experiment
#' 
#' @param hicexp A hicexp object.
#' @param zero.p The proportion of zeros in a row to filter by. 
#'     If the proportion of zeros in a row is less than zero.p
#'     the row will be filtered out.
#' @param A.min The minimum average expression value
#'     (row mean) for an interaction pair. If the 
#'     interaction pair has an average expression 
#'     value less than A.min the row will be filtered
#'     out.
#'
#' @export 

hic_filter <- function(hicexp, zero.p = 0.8, A.min = 5) {
  if (zero.p < 0 | zero.p > 1) {
    stop("zero.p must be in [0,1]")
  }
  if (A.min < 0) {
    stop("A.min must be >= 0")
  }
  # make matrix of IFs
  IF_mat <- hicexp@hic_table[, 5:ncol(hicexp@hic_table), with = FALSE] %>% as.matrix()
  # get row Avg expression
  A <- apply(IF_mat, 1, mean)
  # filter by Avg expression of interacting pair
  keep.A <- A > A.min
  # filter by proportion of zero IFs for each interaction pair
  zeros <- rowSums(IF_mat == 0) / ncol(IF_mat) # proportion of zeros by row
  keep.z <- zeros < zero.p
  # final keep vector
  keep <- keep.A & keep.z
  # filter table
  hicexp@hic_table <- hicexp@hic_table[keep,]
  return(hicexp)
}