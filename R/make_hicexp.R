#' Make Hi-C experiment object from data
#' 
#' @export
#' @importFrom dplyr full_join
#' @import data.table
#' 
#' @param ... Hi-C data. Data must in sparse upper triangular 
#'     format with 4 columns: chr, region1, region2, IF.
#' @param groups A vector of the experimental groups 
#'     corresponding to each Hi-C data object entered.
#' @param covariates Optional data.frame containing 
#'     covariate information for your Hi-C experiment.
#'     Some examples are enzyme used, batch number, etc.
#'     Should have the same number of rows as the number
#'     of Hi-C data objects entered and columns corresponding
#'     to a covariates. 
#' @details Use this function to create a hicexp object for
#'     analysis in HiCcompare2.
#' 
#' 
make_hicexp <- function(..., groups, covariates = NULL) {
  tabs <- list(...)
  # set column names of data 
  tabs <- lapply(tabs, function(x) {
     colnames(x) <- c('chr', 'region1', 'region2', 'IF') 
     return(x)
     })
  
  # initialize values
  hic_tables <- list() # initialize list for tables for each condition
  num_per_group <- table(groups)
  uniq_groups <- unique(groups)
  res <- vector(length = length(uniq_groups))
  
  # check for correct input
  if (length(groups) != length(tabs)) {
    stop("Length of groups must equal the number of Hi-C data objects entered")
  }
  if (sum(num_per_group < 2) > 0) {
    stop("Each experimental condition must have at least 2 samples. If you have less than 2 samples per group use HiCcompare instead")
  }
  if (!is.null(covariates)) {
    # check for data.frame
    if (!is(covariates, "data.frame")) {
      stop("Enter a data.frame for covariates")
    }
    if (nrow(covariates) != length(tabs)) {
      stop("Number of rows in covariates should correspond to number of Hi-C data objects entered")
    }
  }
  
  # cycle through groups to create a table for each experimental condition
  for(i in uniq_groups) {
    current_group <- which(groups == i)
    tmp_table <- dplyr::full_join(tabs[[current_group[1]]], tabs[[current_group[2]]], by = c('chr' = 'chr', 'region1' = 'region1', 'region2' = 'region2'))
    if (num_per_group[i] > 2) {
      for(j in 3:num_per_group[i]) {
        tmp_table <- dplyr::full_join(tmp_table, tabs[[current_group[j]]], by = c('chr' = 'chr', 'region1' = 'region1', 'region2' = 'region2'))
      }
    }
    # rename table columns & replace NA IFs with 0
    colnames(tmp_table)[4:ncol(tmp_table)] <- paste0("IF", 1:(ncol(tmp_table) - 3))
    tmp_table[is.na(tmp_table)] <- 0
    
    
    # calculate resolution
    bins <- unique(c(tmp_table$region1, tmp_table$region2))
    bins <- bins[order(bins)]
    resolution <- min(diff(bins))
    res[i] <- resolution
    # calculate distance
    tmp_table <- data.table::as.data.table(tmp_table)
    tmp_table[, `:=`(D, abs(region2 - region1) / resolution)]
    # rearrange columns
    tmp_table <- tmp_table[, c(1, 2, 3, ncol(tmp_table), 4:(ncol(tmp_table) - 1)), with = FALSE]
    
    # set table in place on condition list
    hic_tables[[i]] <- tmp_table
  }
  # check that all resolutions are equal
  if (length(unique(res)) > 1) {
    stop("Resolution of all datasets must be equal.")
  }
  
  # name items in hic_table
  names(hic_tables) <- paste0("Group", 1:length(hic_tables))
  
  # make metadata 
  metadata <- data.frame(group = groups)
  row.names(metadata) <- paste0("Sample", 1:length(groups))
  if (!is.null(covariates)) {
    metadata <- cbind(metadata, covariates)
  }
  
  # put into hicexp object
  experiment <- new("hicexp", hic_tables = hic_tables, metadata = metadata, resolution = resolution)
  return(experiment)
}


# tabs <- list(r1, r2, r3, r4, r5, r6, r7)
# tabs <- lapply(tabs, function(x) cbind('chr22', x))
# groups <- c(1, 1, 1, 2, 2, 2, 2)
# covariates <- data.frame(enzyme = c('mobi', 'mboi', 'mboi', 'dpnii', 'dpnii', 'dpnii', 'dpnii'), batch = c(1, 2, 1, 2, 1, 2, 2))
# 
# hic_exp <- make_hicexp(r1, r2, r3, r4, r5, r6, r7, groups = groups, covariates = covariates)
