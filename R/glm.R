#' Perform difference detection on a Hi-C experiment
#' 
#' 
#' @import edgeR
#' @importFrom dplyr %>%

hic_glm <- function(hicexp, parallel = FALSE) {
  # check to make sure hicexp is normalized
  if (!hicexp@normalized) {
    warning("You should normalize the data before entering it into hic_glm")
  }
  # First need to split up hic_table by chr and then by distance
  # split up data by chr
  table_list <- split(hicexp@hic_table, hicexp@hic_table$chr)
  # for each chr create list of distance matrices
  table_list <- lapply(table_list, .get_dist_tables)
  # combine list of lists into single list of tables
  table_list <- do.call(c, table_list)
  # input each of the chr-distance tables into a DGElist object for edgeR
  dge_list <- lapply(table_list, .hictable2DGEList, covariates = hicexp@metadata)
}


# function to convert hic_table to a DGEList object
# ??? Might need to add "gene" names as an ID for each row if the original hic_table in case edgeR removes rows/moves them around during processing so that they match
# back up with the original chrs, region1, region2, etc.
.hictable2DGEList <- function(hic_table, covariates) {
  # convert IFs into a matrix
  IFs <- hic_table[, 5:(ncol(hic_table)), with = FALSE] %>% as.matrix
  # create DGEList 
  dge <- DGEList(counts = IFs, samples = covariates)
  return(dge)
}
