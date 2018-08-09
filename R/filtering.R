#' Perform filtering on a Hi-C experiment
#' 
#' @param hicexp A hicexp object.
#' @param zero.p The proportion of zeros in a row to filter by. 
#'     If the proportion of zeros in a row is <= zero.p
#'     the row will be filtered out, i.e. zero.p = 1 means
#'     nothing is filtered based on zeros and zero.p = 0 
#'     will filter rows that have any zeros.
#' @param A.min The minimum average expression value
#'     (row mean) for an interaction pair. If the 
#'     interaction pair has an average expression 
#'     value less than A.min the row will be filtered
#'     out.
#' @details This function is used to filter out
#'     the interactions that have low average IFs
#'     or large numbers of 0 IF values. If you have
#'     already performed filtering when making your
#'     hicexp object do not use this again. As these
#'     interactions are not very interesting and
#'     are commonly false positives during difference
#'     detection it is better to remove them from
#'     the dataset. Additionally, filtering will
#'     help speed up the run time of multiHiCcompare. 
#'     Filtering can be performed before or after 
#'     normalization, however the best computational
#'     speed gain will occur when filtering is done
#'     before normalization.
#' @return A hicexp object.
#' @export
#' @examples 
#' data("hicexp")
#' hicexp <- hic_filter(hicexp)

hic_filter <- function(hicexp, zero.p = 0.8, A.min = 5) {
  if (zero.p < 0 | zero.p > 1) {
    stop("zero.p must be in [0,1]")
  }
  if (A.min < 0) {
    stop("A.min must be >= 0")
  }
  # make matrix of IFs
  IF_mat <- hic_table(hicexp)[, 5:ncol(hic_table(hicexp)), with = FALSE] %>% as.matrix()
  # get row Avg expression
  A <- apply(IF_mat, 1, mean)
  # filter by Avg expression of interacting pair
  keep.A <- A > A.min
  # filter by proportion of zero IFs for each interaction pair
  zeros <- rowSums(IF_mat == 0) / ncol(IF_mat) # proportion of zeros by row
  keep.z <- zeros <= zero.p
  # final keep vector
  keep <- keep.A & keep.z
  # filter table
  slot(hicexp, "hic_table") <- hic_table(hicexp)[keep,]
  return(hicexp)
}
