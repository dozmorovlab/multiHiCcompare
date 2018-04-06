# Set up data classes

setClass("hicexp",
         representation(
           hic_tables = "list",
           metadata= "data.frame",
           resolution = "numeric")
         # prototype =  prototype(hic_tables= NULL, metadata=NULL)
)


setMethod("show", "hicexp", function(object) {
  # get unique groups
  uniq_groups <- length(object@hic_tables)
  # get number samples per groups
  s_per_groups <- table(object@metadata$group)
  cat("Hi-C Experiment Object \n")
  cat(paste0(uniq_groups, " experimental groups \n"))
  for (i in 1:length(s_per_groups)) {
    cat(paste0("Group ", i, ' has ', s_per_groups[i], ' samples \n'))
  }
})
