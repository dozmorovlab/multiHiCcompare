#' Perform exact test based difference detection on a Hi-C experiment
#' 
#' @export
#' @import edgeR
#' @importFrom dplyr %>%

hic_exactTest <- function(hicexp, parallel = FALSE, p.method = "fdr", Plot = TRUE) {
  # check to make sure hicexp is normalized
  if (!hicexp@normalized) {
    warning("You should normalize the data before entering it into hic_glm")
  }
  # check to make sure there are only 2 groups
  if (hicexp@metadata$group %>% unique() %>% length() != 2) {
    stop("If you are making a comparison where the number of groups is not 2 or you have covariates use hic_glm() instead.")
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
  # estimate dispersion for data
  # check number of samples per group
  if (hicexp@metadata$group %>% length() == 2) { # IF no replicates use edgeR's recommended method to estimate dispersion
    if (parallel) {
      dge_list <- BiocParallel::bplapply(dge_list, edgeR::estimateGLMCommonDisp, method="deviance", robust=TRUE, subset=NULL)
    } else {
      dge_list <- lapply(dge_list, edgeR::estimateGLMCommonDisp, method="deviance", robust=TRUE, subset=NULL)
    }
  } else { # If replicates for condition then use standard way for estimating dispersion
    if (parallel) {
      dge_list <- BiocParallel::bplapply(dge_list, edgeR::estimateDisp, design=model.matrix(~hicexp@metadata$group))
    } else {
      dge_list <- lapply(dge_list, edgeR::estimateDisp, design=model.matrix(~hicexp@metadata$group))
    }
  }
  
  # perform exact test when the only factor is group
  et <- lapply(dge_list, edgeR::exactTest)
  
  # reformat results back into hicexp object
  comparison <- mapply(.et_reformat, et, table_list, SIMPLIFY = FALSE, MoreArgs = list(p.method = p.method))
  comparison <- data.table::rbindlist(comparison)
  hicexp@comparison <- comparison
  
  # plot
  if (Plot) {
    MD.composite(hicexp)
  }
  
  # return results
  return(hicexp)
}


# Reformat exact test results and adjust p-values
## !!! current p-value adjustment is taking place on per distance basis
.et_reformat <- function(et_result, hic_table, p.method) {
  # create table of location info and p-value results
  result <- cbind(hic_table[, 1:4, with = FALSE], et_result$table)
  colnames(result)[7] <-"p.value"
  # adjust p-values
  result$p.adj <- p.adjust(result$p.value, method = p.method)
  return(result)
}



#' Function to perform GLM differential analysis on Hi-C experiment
#' 
#' @param hicexp A hicexp object,
#' @param design A design matrix for the GLM.
#' @param contrast A vector containing the contrast you want to test 
#' @param coef A value or vect or the coefficients to test.
#' @param method The test method to be performed. Should be one of
#'    "QLFTest", "LRTest", or "Treat".
#' @param M The log2 fold change value for a TREAT analysis.
#' @param parallel Logical, Shoudl parallel processing be used?

hic_glm <- function(hicexp, design, contrast = NA, coef = NA, method = "QLFTest", M = 1, parallel = FALSE) {
  # check to make sure hicexp is normalized
  if (!hicexp@normalized) {
    warning("You should normalize the data before entering it into hic_glm")
  }
  # contrast & coef input
  if (is.na(contrast) & is.na(coef)) {
    stop("You must enter a value for contrast or a coef, but not both")
  }
  if (!is.na(contrast) & !is.na(coef)) {
    stop("Please enter either a value for contrast or coef but NOT both.")
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
  # estimate dispersion for data
  # ??? might need to adjust this to use the design matrix
  if (parallel) {
    dge_list <- BiocParallel::bplapply(dge_list, edgeR::estimateDisp, design = design)
  } else {
    dge_list <- lapply(dge_list, edgeR::estimateDisp, design = design)
  }
  
  # fit the GLM
  if (parallel) {
    fit <- BiocParallel::bplapply(dge_list, edgeR::glmQLFit, design = design)
  } else {
    fit <- lapply(dge_list, edgeR::glmQLFit, design = design)
  }
  
  ## Perform test based on method and contrast/coef specified 
  # ?? should probably add parallel option here too
  # apply QL F-test test
  if (method == "QLFTest") {
    if (is.na(coef)) {
      result <- lapply(fit, edgeR::glmQLFTest, contrast = contrast)
    } else {
      result <- lapply(fit, edgeR::glmQLFTest, coef = coef)
    }
  }
  # Likelihood Ratio Test
  if (method == "LRTest") {
    if (is.na(coef)) {
      result <- lapply(fit, edgeR::glmLRT, contrast = contrast)
    } else {
      result <- lapply(fit, edgeR::glmLRT, coef = coef)
    }
  }
  # TREAT Analysis based on a minimum log2 fold change specified as M
  if (method == "Treat") {
    if (is.na(coef)) {
      result <- lapply(fit, edgeR::glmTreat, contrast = contrast, lfc = M)
    } else {
      result <- lapply(fit, edgeR::glmTreat, coef = coef, lfc = M)
    }
  }
  
  # format results for differentially interacting regions
  
}


# reformat results of GLM
.glm_reformat <- function(result, hic_table, p.method) {
  # create table of location info and p-value results
  result <- cbind(hic_table[, 1:4, with = FALSE], result$table)
  colnames(result)[ncol(result)] <-"p.value"
  # adjust p-values
  result$p.adj <- p.adjust(result$p.value, method = p.method)
  return(result)
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
