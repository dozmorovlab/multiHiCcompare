#' Manhattan plot function for results of multiHiCcompare
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
#'     plotted. "count" produces a plot where the y-values
#'     are the -log10(1 / S) where S is the number of times
#'     the region was found significant. 
#' @param return_df Logical, should the data.frame used to
#'     generate the plot be returned?
#' @param alpha The p-value cut off to be used for 
#'     calling an interaction significant. This is
#'     only used if method = 'count'. Defaults to 
#'     0.05.
#' @details This function is used to create a manhattan
#'     plot for the significance of all genomic regions 
#'     in the dataset. These correspond to the rows (or columns)
#'     of the upper triangle of the full Hi-C matrix. Each genomic
#'     region of the Hi-C dataset has multiple interactions it
#'     is involved in and the significance of all of these can 
#'     be visualized with \code{method = "standard"}. 
#'     Alternatively the p-values for all these interactions
#'     can be combined using either Fisher's method or the
#'     Stouffer-Liptak method of combining p-values. Additionally
#'     the "count" option will plot based on the number of times
#'     each region was found to be involved in a signficantly 
#'     different interaction. The 
#'     manhattan plot can be used to identify "hotspot"
#'     regions of the genome where major differences
#'     seem to be located based on the results of a multiHiCcompare
#'     analysis.
#' @return A manhattan plot and optionally the data.frame used
#'     to generate the manhattan plot.
#' @importFrom qqman manhattan
#' @importFrom metap sumz sumlog
#' @export
#' @examples
#' data("hicexp_diff")
#' manhattan_hicexp(hicexp_diff, method = "fisher")


manhattan_hicexp <- function(hicexp, method = "standard", return_df = FALSE, alpha = 0.05) {
  # check input
  method <- match.arg(method, c("standard", "fisher", "stouffer", "count"), 
                      several.ok = FALSE)
  # check that data has been compared
  if (nrow(results(hicexp)) < 1) {
    stop("Differences must be detected before making a manhattan plot.")
  }
  
  if (method == "standard") {
    # make data.frame for plotting
    man.df <- data.frame(BP = c(results(hicexp)$region1, 
                                results(hicexp)$region2),
                         CHR = as.numeric(c(results(hicexp)$chr, 
                                            results(hicexp)$chr)), 
                         P = c(results(hicexp)$p.adj, 
                               results(hicexp)$p.adj))
    # plot
    suppressWarnings(qqman::manhattan(man.df))
  }
  
  if (method == "fisher") {
    # make aggregate p-value for regions
    regions <- c(paste0(results(hicexp)$chr, ':', results(hicexp)$region1),
                 paste0(results(hicexp)$chr, ':', results(hicexp)$region2))
    p.values <- c(results(hicexp)$p.adj, results(hicexp)$p.adj)
    
    ## Fisher method
    # fisher_aggregate <- aggregate(p.values, by = list(regions), 
    #                               FUN = function(p) {
    #   combined <- -2 * sum(log(p))
    #   c.pval <- pchisq(combined, df = 2 * length(p))
    #   return(c.pval)
    # })
    
    fisher_aggregate <- aggregate(p.values, 
                                  by = list(regions), 
                                  FUN = function(x) {
                                    if (length(x) > 1) {
                                      metap::sumlog(x)$p
                                    } else {
                                      x
                                    }
                                    })
    
    fisher_aggregate <- cbind(read.table(text = fisher_aggregate$Group.1, 
                                         sep = ":"), fisher_aggregate$x)
    colnames(fisher_aggregate) <- c("CHR", "BP", "P")
    
    # plot combined p-value manahttan plots
    suppressWarnings(qqman::manhattan(fisher_aggregate))
    
    man.df <- fisher_aggregate
  }
  
  if (method == "stouffer") {
    # make aggregate p-value for regions
    regions <- c(paste0(results(hicexp)$chr, ':', results(hicexp)$region1),
                 paste0(results(hicexp)$chr, ':', results(hicexp)$region2))
    p.values <- c(results(hicexp)$p.adj, results(hicexp)$p.adj)
    
    ## Stouffer-Liptak method
    # stouffer_liptak_aggregate <- aggregate(p.values, by = list(regions), 
    #                                        FUN = function(p) {
    #   zi <- qnorm(p, lower.tail = FALSE)
    #   Z <- sum(zi) / sqrt(length(p))
    #   c.pval <- pnorm(Z, lower.tail = FALSE)
    #   return(c.pval)
    # })
    
    stouffer_liptak_aggregate <- aggregate(p.values, 
                                           by = list(regions), 
                                           FUN = function(x) {
                                             if (length(x) > 1) {
                                               suppressWarnings(metap::sumz(x)$p) 
                                             } else {
                                               x
                                             }
                                            })
    
    stouffer_liptak_aggregate <- cbind(read.table(text = stouffer_liptak_aggregate$Group.1,
                                                  sep = ":"), stouffer_liptak_aggregate$x)
    colnames(stouffer_liptak_aggregate) <- c("CHR", "BP", "P")
    # make sure there are no zero p-values
    stouffer_liptak_aggregate$P[stouffer_liptak_aggregate$P == 0] <- .Machine$double.xmin
    
    suppressWarnings(qqman::manhattan(stouffer_liptak_aggregate))
    
    man.df <- stouffer_liptak_aggregate
  }
  
  
  if (method == 'count') {
    # make aggregate count for regions
    regions <- c(paste0(results(hicexp)$chr, ':', results(hicexp)$region1),
                 paste0(results(hicexp)$chr, ':', results(hicexp)$region2))
    p.values <- c(results(hicexp)$p.adj, results(hicexp)$p.adj)
    count <- ifelse(p.values < alpha, 1, 0)
    
    ## count method
    count_aggregate <- aggregate(count, by = list(regions), 
                                 FUN = function(cnt) {
      c.sum <- sum(cnt)
      c.pval <- 1 / sqrt(c.sum)
      return(c.pval)
    })
    
    count_aggregate$x[is.infinite(count_aggregate$x)] <- 1 # replace any regions that had counts of 0 significant with pseudo-pval of 1
    count_aggregate <- cbind(read.table(text = count_aggregate$Group.1, sep = ":"),
                             count_aggregate$x)
    colnames(count_aggregate) <- c("CHR", "BP", "P")
    
    # plot combined p-value manahttan plots
    suppressWarnings(qqman::manhattan(count_aggregate, ylab = "-log10(1/S)"))
    
    man.df <- count_aggregate
  }
  
  
  # return man.df if requested
  if (return_df) {
    return(man.df)
  } 
}