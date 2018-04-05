# Set up data classes

setClass("hicexp",
         representation(
           hic_matrices = "list",
           groups= "numeric"),#,
           # covariates = "data.frame"),
         prototype =  prototype(hic_matrices= NULL, groups=NULL)#, covariates=data.frame())
)


setMethod("show", "hicexp", function(object) {
  # get unique groups
  uniq_groups <- unique(object@groups)
  # get number samples per groups
  s_per_groups <- table(object@groups)
  cat(paste0("Hi-C Experiment Object \n", 
             ""))
})
