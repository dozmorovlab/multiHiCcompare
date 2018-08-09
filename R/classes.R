#' An S4 class for working with Hi-C data 
#' 
#' @slot hic_table A data.table containing the sparse upper triangular matrix 
#'     for your Hi-C data.
#' @slot comparison The results of a multiHiCcompare comparison.
#' @slot metadata Data.frame for covariate information.
#' @slot resolution The resolution of the dataset.
#' @slot normalized Indicator for if data has been normalized.
#' @exportClass Hicexp
#' @import methods

setClass("Hicexp",
         representation(
           hic_table = "data.table",
           comparison = "data.table",
           metadata= "data.frame",
           resolution = "numeric",
           normalized = "logical")
         # prototype =  prototype(hic_tables= NULL, metadata=NULL)
)

#' Print information about a HiCexp object
#' @rdname show
#' @aliases show,character,ANY-method
#' @exportMethod show
#' @param object A Hicexp object
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

# set generics

#' @docType methods
#' @rdname hic_table
setGeneric("hic_table", function(x) {
  standardGeneric("hic_table")
})

#' @docType methods
#' @rdname results
setGeneric("results", function(x) {
  standardGeneric("results")
})

#' @docType methods
#' @rdname meta
setGeneric("meta", function(x) {
  standardGeneric("meta")
})

#' @docType methods
#' @rdname resolution
setGeneric("resolution", function(x) {
  standardGeneric("resolution")
})

#' @docType methods
#' @rdname normalized
setGeneric("normalized", function(x) {
  standardGeneric("normalized")
})

# Accessors
#' Print the hic_table
#' @exportMethod hic_table
#' @rdname hic_table
#' @aliases hic_table,character,ANY-method
#' @param x The Hicexp object
setMethod("hic_table", "Hicexp", function(x) x@hic_table)

#' Print the results
#' @exportMethod results
#' @rdname results
#' @aliases results,character,ANY-method
#' @param x The Hicexp object
setMethod("results", "Hicexp", function(x) x@comparison)

#' Print the metadata
#' @exportMethod meta
#' @rdname meta
#' @aliases meta,character,ANY-method
#' @param x The Hicexp object
setMethod("meta", "Hicexp", function(x) x@metadata)

#' Print the resolution
#' @exportMethod resolution
#' @rdname resolution
#' @aliases resolution,character,ANY-method
#' @param x The Hicexp object
setMethod("resolution", "Hicexp", function(x) x@resolution)

#' Print the indicator for if the data is normalized
#' @exportMethod normalized
#' @rdname normalized
#' @aliases normalized,character,ANY-method
#' @param x The Hicexp object
setMethod("normalized", "Hicexp", function(x) x@normalized)
