#' Perform a permutation test to check enrichment of a genomic
#'     feature with DIRs detected by multiHiCcompare
#'    
#' @param hicexp A Hicexp object which has been compared.  
#' @param feature A GRanges object containing locations for
#'     a genomic feature you would like to test for enrichment
#'     in the differentially interacting regions (DIRs).
#' @param p.adj_cutoff The adjusted p-value cutoff for
#'      declaring a region significant. See ?topDirs for
#'      more information. Defaults to 10^-10
#' @param logfc_cutoff The log fold change cutoff for
#'      a region to be declared significant. See ?topDirs
#'      for more information. Defaults to 1.
#' @param num.perm The number of permutations to run.
#'      Defaults to 1000. 
#' @return The permutation p-value
#' @export
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps GRanges
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @examples 
#' \dontrun{
#' data("hicexp_diff")
#' data("hg19_cyto")
#' perm_test(hicexp_diff, hg19_cyto)
#' }

perm_test <- function(hicexp, feature, p.adj_cutoff = 10^-10, logfc_cutoff = 1, num.perm = 1000) {
  # check that feature is a GRanges object
  if (!is(feature, "GRanges")) {
    stop("feature must be a GRanges object")
  }
  # convert seqnames of feature to include chr
  seqlevelsStyle(feature) <- "UCSC"
  # check that data has been compared
  if (nrow(results(hicexp)) < 1) {
    stop("Differences must be detected before making a manhattan plot.")
  }
  
  # make a GRanges object for the DIRs
  sig.regions <- topDirs(hicexp, logfc_cutoff = logfc_cutoff, p.adj_cutoff = p.adj_cutoff, return_df = 'bed') 
  sig.regions <- makeGRangesFromDataFrame(sig.regions, 
                                             seqnames.field = 'chr', 
                                             start.field = 'start', 
                                             end.field = 'end',
                                             keep.extra.columns = TRUE)
  
  # find overlaps between feature and sig.regions
  olaps <- findOverlaps(sig.regions, feature)
  message("There are ", length(olaps), " overlaps between your genomic feature and the significant regions")
  
  # set up permutation test
  sample.size <- length(sig.regions) # Number of DIRs
  n.olap      <- length(olaps) # How many of them overlap with DE genes
  num.olap    <- vector(mode = "numeric", length = num.perm) # A vector to store the number of overlaps during permutations
  
  # make background regions GRanges object
  regions <- topDirs(hicexp, logfc_cutoff = 0, logcpm_cutoff = -10,
                     D_cutoff = 0, p.adj_cutoff = 1, alpha = 2, 
                     return_df = 'bed' ) 
  # Order regions
  regions <- regions[order(chr, start, end), ]
  # Remove unnecessary columns
  regions <- regions[, c('chr', 'start', 'end')]
  
  # perform permutation test
  for(i in 1:num.perm) {
    # Select random regions from the genome of the same size as DIRs
    sampled <- sample(1:nrow(regions), sample.size) 
    sampled.regions <- regions[sampled, ]
    sampled.gr <- GRanges(sampled.regions$chr, 
                          IRanges(start = sampled.regions$start, 
                                  end = sampled.regions$end))
    # Find how many times they overlap with DE genes
    tmp.olap <- findOverlaps(sampled.gr, feature)
    # Keep the results
    num.olap[i] <- length(tmp.olap)
  }
  # Calculate permutation p-value
  p.value <- (sum(num.olap > n.olap) + 1) / (num.perm + 1)
  return(p.value)
}
