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
#' @param remove.regions A GenomicRanges object indicating
#'     specific regions to be filtered out. By default
#'     this is the hg19 centromeric, gvar, and stalk
#'     regions. Also included in the package is
#'     hg38_cyto. If your data is not hg19 you will 
#'     need to substitute this file. To choose not 
#'     to filter any regions set regions = NULL.
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
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @examples 
#' data("hicexp2")
#' hicexp2 <- hic_filter(hicexp2)

hic_filter <- function(hicexp, zero.p = 0.8, A.min = 5, remove.regions = hg19_cyto) {
  if (zero.p < 0 || zero.p > 1) {
    stop("zero.p must be in [0,1]")
  }
  if (A.min < 0) {
    stop("A.min must be >= 0")
  }
  # make matrix of IFs
  IF_mat <- as.matrix(hic_table(hicexp)[, -c("chr", "region1", "region2", "D"), with = FALSE])
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
  
  ## TODO: make this generalizable to any genome
  if(!is.null(remove.regions)) {
    # filter out hg19 acen, gvar, stalk regions
    tab <- hic_table(hicexp)
    # make GRanges object for regions of hic_table
    region1 <- tab[, c('chr', 'region1')]
    region1[, `:=`(end = region1 + resolution(hicexp) - 1)]
    region2 <- tab[, c('chr', 'region2')]
    region2[, `:=`(end = region2 + resolution(hicexp) - 1)]
    region1 <- GenomicRanges::makeGRangesFromDataFrame(region1, seqnames.field = 'chr', start.field = 'region1', end.field = 'end')
    region2 <- GenomicRanges::makeGRangesFromDataFrame(region2, seqnames.field = 'chr', start.field = 'region2', end.field = 'end')
    # overlap
    olaps1 <- GenomicRanges::findOverlaps(region1, remove.regions)
    olaps2 <- GenomicRanges::findOverlaps(region2, remove.regions)
    remove <- unique(c(olaps1@from, olaps2@from))
    # filter table
    if (length(remove) > 0) {
      slot(hicexp, "hic_table") <- hic_table(hicexp)[-remove,]
    }
  }
  
  return(hicexp)
}
