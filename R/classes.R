#' An S4 class for working with Hi-C data 
#' 
#' @slot hic_table A data.table containing the sparse upper triangular matrix 
#'     for your Hi-C data.
#' @slot comparison The results of a multiHiCcompare comparison.
#' @slot metadata Data.frame for covariate information.
#' @slot resolution The resolution of the dataset.
#' @slot normalized Indicator for if data has been normalized.

setClass("Hicexp",
         representation(
           hic_table = "data.table",
           comparison = "data.frame",
           metadata= "data.frame",
           resolution = "numeric",
           normalized = "logical")
         # prototype =  prototype(hic_tables= NULL, metadata=NULL)
)


setMethod("show", "Hicexp", function(object) {
  # get unique groups
  uniq_groups <- length(unique(object@metadata$group))
  # get number samples per groups
  s_per_groups <- table(object@metadata$group)
  cat("Hi-C Experiment Object \n")
  cat(paste0(uniq_groups, " experimental groups \n"))
  for (i in 1:length(s_per_groups)) {
    cat(paste0("Group ", i, ' has ', s_per_groups[i], ' samples \n'))
  }
  if (object@normalized) {
    cat("Data has been normalized")
  }
})
