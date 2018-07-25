#' Manhattan plot function for results of HiCcompare2
#' 
#' @param hicexp A hicexp object that has had differences
#'     detected
#' @param method string denoting the p-value method to
#'     use for plotting. Options are "standard", "fisher",
#'     and "stouffer". "standard" plots a manhattan plot
#'     using all individual p-values. "fisher" uses 
#'     Fisher's method for combining p-values to combine
#'     the p-values for each region which are then plotted.
#'     "stouffer" uses the Stouffer-Liptak method for 
#'     combining p-values for each region which are then 
#'     plotted. 
#' @param return_df Logical, should the data.frame used to
#'     generate the plot be returned?
#' @details This function is used to create a manhattan
#'     plot for the significance of all genomic regions 
#'     in the dataset. These correspond to the rows (or columns)
#'     of the upper triangle of the full Hi-C matrix. Each genomic
#'     region of the Hi-C dataset has multiple interactions it
#'     is involved in and the significance of all of these can 
#'     be visualized with \code{method = "standard"}. 
#'     Alternatively the p-values for all these interactions
#'     can be combined using either Fisher's method or the
#'     Stouffer-Liptak method of combining p-values. The 
#'     manhattan plot can be used to identify "hotspot"
#'     regions of the genome where major differences
#'     seem to be located based on the results of a HiCcompare2
#'     analysis.
#' @return A manhattan plot and optionally the data.frame used
#'     to generate the manhattan plot.
#' @importFrom qqman manhattan
#' @export
#' @examples
#' data("hicexp_diff")
#' manhattan_hicexp(hicexp_diff)


manhattan_hicexp <- function(hicexp, method = "standard", return_df = FALSE) {
  # check input
  method <- match.arg(method, c("standard", "fisher", "stouffer"), several.ok = FALSE)
  # check that data has been compared
  if (nrow(hicexp@comparison) < 1) {
    stop("Differences must be detected before making a manhattan plot.")
  }
  
  if (method == "standard") {
    # make data.frame for plotting
    man.df <- data.frame(BP = c(hicexp@comparison$region1, hicexp@comparison$region2), CHR = as.numeric(c(hicexp@comparison$chr, hicexp@comparison$chr)), 
                         P = c(hicexp@comparison$p.adj, hicexp@comparison$p.adj))
    # plot
    suppressWarnings(qqman::manhattan(man.df))
  }
  
  if (method == "fisher") {
    # make aggregate p-value for regions
    regions <- c(paste0(hicexp@comparison$chr, ':', hicexp@comparison$region1), paste0(hicexp@comparison$chr, ':', hicexp@comparison$region2))
    p.values <- c(hicexp@comparison$p.adj, hicexp@comparison$p.adj)
    
    ## Fisher method
    fisher_aggregate <- aggregate(p.values, by = list(regions), FUN = function(p) {
      combined <- -2 * sum(log(p))
      c.pval <- pchisq(combined, df = 2 * length(p))
      return(c.pval)
    })
    
    fisher_aggregate <- cbind(read.table(text = fisher_aggregate$Group.1, sep = ":"), fisher_aggregate$x)
    colnames(fisher_aggregate) <- c("CHR", "BP", "P")
    
    # plot combined p-value manahttan plots
    suppressWarnings(qqman::manhattan(fisher_aggregate))
    
    man.df <- fisher_aggregate
  }
  
  if (method == "stouffer") {
    # make aggregate p-value for regions
    regions <- c(paste0(hicexp@comparison$chr, ':', hicexp@comparison$region1), paste0(hicexp@comparison$chr, ':', hicexp@comparison$region2))
    p.values <- c(hicexp@comparison$p.adj, hicexp@comparison$p.adj)
    
    ## Stouffer-Liptak method
    stouffer_liptak_aggregate <- aggregate(p.values, by = list(regions), FUN = function(p) {
      zi <- qnorm(p/2)
      Z <- sum(zi) / sqrt(length(p))
      c.pval <- pnorm(Z)
      return(c.pval)
    })
    
    stouffer_liptak_aggregate <- cbind(read.table(text = stouffer_liptak_aggregate$Group.1, sep = ":"), stouffer_liptak_aggregate$x)
    colnames(stouffer_liptak_aggregate) <- c("CHR", "BP", "P")
    # make sure there are no zero p-values
    stouffer_liptak_aggregate$P[stouffer_liptak_aggregate$P == 0] <- .Machine$double.xmin
    
    suppressWarnings(qqman::manhattan(stouffer_liptak_aggregate))
    
    man.df <- stouffer_liptak_aggregate
  }
  
  
  # return man.df if requested
  if (return_df) {
    return(man.df)
  } 
}   



