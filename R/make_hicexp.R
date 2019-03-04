#' Make Hi-C experiment object from data
#' 
#' @export
#' @importFrom dplyr full_join left_join right_join
#' @import data.table
#' 
#' @param ... Hi-C data. Data must in sparse upper triangular 
#'     format with 4 columns: chr, region1, region2, IF or
#'     in 7 column BEDPE format with columns chr, start1, 
#'     end1, chr, start2, end2, IF.
#' @param data_list Alternate way to enter data. If you have
#'     your Hi-C data in the form of a list already with each
#'     entry of the list representing a sample use this option.
#' @param groups A vector of the experimental groups 
#'     corresponding to each Hi-C data object entered.
#'     If it is not in factor form when entered it will
#'     be converted to a factor.
#' @param covariates Optional data.frame containing 
#'     covariate information for your Hi-C experiment.
#'     Some examples are enzyme used, batch number, etc.
#'     Should have the same number of rows as the number
#'     of Hi-C data objects entered and columns corresponding
#'     to covariates. 
#' @param remove_zeros Logical, should rows with 1 or more
#'     zero IF values be removed?
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
#' @param filter Logical, should filtering be performed?
#'     Defaults to TRUE. If TRUE it will filter out
#'     the interactions that have low average IFs
#'     or large numbers of 0 IF values. As these
#'     interactions are not very interesting and
#'     are commonly false positives during difference
#'     detection it is better to remove them from
#'     the dataset. Additionally, filtering will
#'     help speed up the run time of multiHiCcompare. 
#'     Filtering can be performed before or after 
#'     normalization, however the best computational
#'     speed gain will occur when filtering is done
#'     before normalization. Filtering parameters
#'     are controlled by the zero.p and A.min options.
#' @param remove.regions A GenomicRanges object indicating
#'     specific regions to be filtered out. By default
#'     this is the hg19 centromeric, gvar, and stalk
#'     regions. Also included in the package is
#'     hg38_cyto. If your data is not hg19 you will 
#'     need to substitute this file. To choose not 
#'     to filter any regions set regions = NULL. NOTE:
#'     if you set filter = FALSE these regions will NOT
#'     be removed. This occurs in conjuction with the 
#'     filtering step.
#' @details Use this function to create a hicexp object for
#'     analysis in multiHiCcompare. Filtering can also be 
#'     performed in this step if the filter option is 
#'     set to TRUE. Filtering parameters are controlled
#'     by the zero.p and A.min options.
#'     
#' @return A hicexp object.
#' @examples 
#' # load data in sparse upper triangular format
#' data("HCT116_r1", "HCT116_r2", "HCT116_r3", "HCT116_r4", 
#'     "HCT116_r5", "HCT116_r6")
#' # make groups & covariate input
#' groups <- factor(c(1, 1, 1, 2, 2, 2))
#' covariates <- data.frame(enzyme = factor(c('mobi', 'mboi', 'mboi',
#'  'dpnii', 'dpnii', 'dpnii')), batch = c(1, 2, 1, 2, 1, 2))
#' # make the hicexp object
#' hicexp <- make_hicexp(HCT116_r1, HCT116_r2, HCT116_r3, HCT116_r4,
#'      HCT116_r5, HCT116_r6, groups = groups, 
#'      covariates = covariates)


make_hicexp <- function(..., data_list = NA, groups, covariates = NULL, 
                        remove_zeros = FALSE,
                        zero.p = 0.8, A.min = 5, filter = TRUE, remove.regions = hg19_cyto) {
  if (!is.na(data_list[1])) {
    tabs <- data_list
  } else {
    tabs <- list(...)
  }
  # check number of columns of input data
  if(ncol(tabs[[1]]) != 4 &  ncol(tabs[[1]]) != 7) {
    stop("You must enter data in 4 column sparse matrix format
         or in 7 column BEDPE format.")
  }
  # if BEDPE pull out only columns we need
  if (ncol(tabs[[1]]) == 7) {
    tabs <- lapply(tabs, function(x) {
      colnames(x) <- c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'IF')
      new.x <- x[, c('chr1', 'start1', 'start2', 'IF')]
      return(new.x)
    })
  }
  # set column names of data 
  tabs <- lapply(tabs, function(x) {
    colnames(x) <- c('chr', 'region1', 'region2', 'IF') 
    return(x)
  })
  
  # initialize values
  # hic_table <- list() # initialize list for tables for each condition
  num_per_group <- table(groups)
  uniq_groups <- unique(groups)
  res <- vector(length = length(uniq_groups))
  
  # convert groups to a factor if it is not already
  if (!is.factor(groups)) {
    groups <- as.factor(groups)
  }
  
  # check for correct input
  if (length(groups) != length(tabs)) {
    stop("Length of groups must equal the number of Hi-C data objects entered")
  }
  if (sum(num_per_group < 2) > 0) {
    warning("Each experimental condition should have at least 2 samples. 
            If you have less than 2 samples per group use HiCcompare instead")
  }
  if (!is.null(covariates)) {
    # check for data.frame
    if (!is(covariates, "data.frame")) {
      stop("Enter a data.frame for covariates")
    }
    if (nrow(covariates) != length(tabs)) {
      stop("Number of rows in covariates should correspond to number 
           of Hi-C data objects entered")
    }
    }
  if (zero.p < 0 | zero.p > 1) {
    stop("zero.p must be in [0,1]")
  }
  if (A.min < 0) {
    stop("A.min must be >= 0")
  }
  
  # cycle through groups to create a table for each experimental condition
  
  tmp_table <- dplyr::full_join(tabs[[1]], tabs[[2]], 
                                by = c('chr' = 'chr', 'region1' = 'region1',
                                       'region2' = 'region2'))
  if (length(tabs) > 2) {
    for(j in 3:length(tabs)) {
      tmp_table <- dplyr::full_join(tmp_table, tabs[[j]], 
                                    by = c('chr' = 'chr', 
                                           'region1' = 'region1',
                                           'region2' = 'region2'))
    }
  }
  # rename table columns & replace NA IFs with 0
  colnames(tmp_table)[4:ncol(tmp_table)] <- paste0("IF", 1:(ncol(tmp_table) - 3))
  tmp_table[is.na(tmp_table)] <- 0
  
  
  # calculate resolution
  bins <- unique(c(tmp_table$region1, tmp_table$region2))
  bins <- bins[order(bins)]
  resolution <- min(diff(bins))
  tmp_table <- data.table::as.data.table(tmp_table)
  tmp_table[, `:=`(D, abs(region2 - region1) / resolution)]
  # rearrange columns
  tmp_table <- tmp_table[, c('chr', 'region1', 'region2', 'D',
                             paste0("IF", 1:(ncol(tmp_table) - 4))), with = FALSE]
  
  # set table in place on condition list
  hic_table <- tmp_table
  # check chr format, if "chr#" change to just the number
  if (!is.numeric(hic_table$chr)) {
    # replace any "chr"
    hic_table[, chr := sub("chr", "", chr)]
    # replace any X 
    hic_table[, chr := sub("X", "23", chr)]
    # replace any Y
    hic_table[, chr := sub("Y", "24", chr)]
    # convert to numeric
    hic_table[, chr := as.numeric(chr)]
  }
  # sort hic_table
  hic_table <- hic_table[order(chr, region1, region2),]
  
  # check that all resolutions are equal
  if (length(unique(res)) > 1) {
    stop("Resolution of all datasets must be equal.")
  }
  # remove rows with 0's
  if (remove_zeros) {
    IFs <- as.matrix(hic_table[, -c("chr", "region1", "region2", "D"), with = FALSE])
    IFs <- IFs == 0 # matrix of TRUE where zeros occur
    keep <- rowSums(IFs) == 0 # if row has no zeros it will sum to 0
    hic_table <- hic_table[keep,]
  }
  
  
  # make metadata 
  metadata <- data.frame(group = groups)
  row.names(metadata) <- paste0("Sample", 1:length(groups))
  if (!is.null(covariates)) {
    metadata <- cbind(metadata, covariates)
  }
  
  # put into hicexp object
  experiment <- new("Hicexp", hic_table = hic_table, 
                    comparison = data.table::data.table(), 
                    metadata = metadata, resolution = resolution, 
                    normalized = FALSE)
  
  # filter
  if (filter) {
    experiment <- hic_filter(experiment, zero.p = zero.p, A.min = A.min, remove.regions = remove.regions)
  }
  
  return(experiment)
  }
