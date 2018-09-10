#' Filter results of multiHiCcompare
#' 
#' @param hicexp A hicexp object which has been compared.
#' @param logfc_cutoff The logFC value you wish to filter
#'    by. Defaults to 1.
#' @param logcpm_cutoff The logCPM cutoff you wish to
#'    filter by. Defaults to 1.
#' @param p.adj_cutoff The p-value cutoff you wish to filter
#'    by. Defaults to 0.01.
#' @param D_cutoff The distance cutoff you wish to filter
#'    by. Interactions with a D < D_cutoff will be filtered.
#'    Defaults to 1. 
#' @param return_df The format for the data.frame returned
#'    by the function. Options are "bed" and "pairedbed".
#' @param alpha The p-value cut off for determining the count
#'     of number of times a region is significant. The adjusted
#'     p-value is used for this. Defaults
#'     to 0.05.
#' @details This function is meant to filter the results of
#'     multiHiCcompare. The top differentially interacting 
#'     regions (DIRs) can be returned by using this function.
#'     When the return_df = "bed" option is set the resulting
#'     data.frame can be input into the plot_pvals or plot_counts
#'     functions to visualize the top DIRs. 
#' @return A data.table containing the filtered results.
#' @export
#' @examples 
#' data('hicexp_diff')
#' topDirs(hicexp_diff)

topDirs <- function(hicexp, logfc_cutoff = 1, logcpm_cutoff = 1, p.adj_cutoff = 0.01,
                           D_cutoff = 1, alpha = 0.05, return_df = "pairedbed") {
  # check that data has been compared
  if (nrow(results(hicexp)) < 1) {
    stop("Differences must be detected before making a manhattan plot.")
  }
  # check input
  return_df <- match.arg(return_df, c("bed", "pairedbed"), 
                      several.ok = FALSE)
  
  # make results object
  res <- results(hicexp)
  
  # if pairedbed just need to filter
  if (return_df == "pairedbed") {
    # filter
    res <- res[logCPM >= logcpm_cutoff & abs(logFC) >= logfc_cutoff & p.adj <= p.adj_cutoff & D >= D_cutoff, ]
    # reformat
    res$chr1 <- paste0('chr', res$chr)
    res$chr2 <- paste0('chr', res$chr)
    res$end1 <- res$region1 + resolution(hicexp) - 1
    res$end2 <- res$region2 + resolution(hicexp) - 1
    res <- res[, c("chr1", "region1", "end1", 'chr2', 'region2', 'end2', 'D', 'logFC', 'logCPM', 'p.value', 'p.adj')]
    colnames(res)[c(2,5)] <- c('start1', 'start2')
  }
  
  # if BED need to aggregate regions
  if (return_df == "bed") {
    # subset to only significant interactions
    res <- res[p.adj < alpha, ]
    # make vector of regions and p-values
    regions <- c(paste0(res$chr, ':', res$region1),
                 paste0(res$chr, ':', res$region2))
    p.values <- c(res$p.adj, res$p.adj)
    dist <- c(res$D, res$D)
    logfc <- c(res$logFC, res$logFC)
    logcpm <- c(res$logCPM, res$logCPM)
    
    
    # aggregate into fisher pvalue
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
    # # aggregate into stouffer pvalue
    # p.values[p.values == 1] <- 0.99999 # change pvalues from 1 so sumz works correctly
    # stouffer_liptak_aggregate <- aggregate(p.values, 
    #                                        by = list(regions), 
    #                                        FUN = function(x) {
    #                                          if (length(x) > 1) {
    #                                            suppressWarnings(metap::sumz(x)$p) 
    #                                          } else {
    #                                            x
    #                                          }
    #                                        })
    # 
    # stouffer_liptak_aggregate <- cbind(read.table(text = stouffer_liptak_aggregate$Group.1,
    #                                               sep = ":"), stouffer_liptak_aggregate$x)
    
    # aggregate counts
    count <- ifelse(p.values < alpha, 1, 0)
    count_aggregate <- aggregate(count, by = list(regions), sum)
    count_aggregate <- cbind(read.table(text = count_aggregate$Group.1, sep = ":"),
                             count_aggregate$x)
    
    # aggregate D
    dist_aggregate <- aggregate(dist, by = list(regions), mean)
    dist_aggregate <- cbind(read.table(text = dist_aggregate$Group.1, sep = ":"),
                             dist_aggregate$x)
    # aggregate logfc
    logfc_aggregate <- aggregate(logfc, by = list(regions), mean)
    logfc_aggregate <- cbind(read.table(text = logfc_aggregate$Group.1, sep = ":"),
                             logfc_aggregate$x)
    # aggregate logcpm
    logcpm_aggregate <- aggregate(logcpm, by = list(regions), mean)
    logcpm_aggregate <- cbind(read.table(text = logcpm_aggregate$Group.1, sep = ":"),
                             logcpm_aggregate$x)
    # Format results
    res <- dplyr::left_join(count_aggregate, dist_aggregate, by = c('V1' = 'V1', 'V2' = 'V2'))
    res <- dplyr::left_join(res, logfc_aggregate, by = c('V1' = 'V1', 'V2' = 'V2'))
    res <- dplyr::left_join(res, logcpm_aggregate, by = c('V1' = 'V1', 'V2' = 'V2'))
    # res <- dplyr::left_join(res, stouffer_liptak_aggregate, by = c('V1' = 'V1', 'V2' = 'V2'))
    res <- dplyr::left_join(res, fisher_aggregate, by = c('V1' = 'V1', 'V2' = 'V2'))
    res$V1 <- paste0('chr', res$V1)
    res$end <- res$V2 + resolution(hicexp) - 1
    colnames(res) <- c('chr', 'start', 'count', 'avgD', 'avgLogFC', 'avgLogCPM',  'avgP.adj', 'end')
    res <- res[, c('chr', 'start', 'end', 'count', 'avgD', 'avgLogFC', 'avgLogCPM',  'avgP.adj')]
    res <- data.table::as.data.table(res)
    res <- res[order(res$count, decreasing = TRUE), ]
    
    # filter
    res <- res[avgLogCPM >= logcpm_cutoff & abs(avgLogFC) >= logfc_cutoff & avgP.adj <= p.adj_cutoff & avgD >= D_cutoff, ]
  }
  
  return(res)
}
