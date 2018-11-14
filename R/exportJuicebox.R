#' Export multiHiCcompare results for visualization in Juicebox
#' 
#' @param hicexp A hicexp object which has been compared.
#' @param logfc_cutoff The logFC value you wish to filter
#'    by. Defaults to 1.
#' @param logcpm_cutoff The logCPM cutoff you wish to
#'    filter by. Defaults to 1.
#' @param p.adj_cutoff The adjusted p-value cutoff you wish to filter
#'    by. Defaults to 0.01.
#' @param D_cutoff The distance cutoff you wish to filter
#'    by. Interactions with a D < D_cutoff will be filtered.
#'    Defaults to 1. 
#' @param file_name The file name of the text file to be saved.
#' @param color A decimal RGB color code. Should be a character
#'     value in form of "0,0,255". Defaults to color code for
#'     blue. This will determine the color of the annotations
#'     on the Juicebox heatmap. 
#' @details This function is meant to filter the results of
#'     multiHiCcompare and export the significant 
#'     differentially interacting regions into a text file which
#'     can be imported into Juicebox as a 2D annotations file. This
#'     will allow you to visualize where your DIRs occur on the heatmap
#'     of the interactions. Please see the included vignette on using
#'     Juicebox to visualize multiHiCcompare results. This can be accessed
#'     with browseVignettes("multiHiCcompare").
#' @return A text file containing annotations for input into Juicebox.
#' @export
#' @examples 
#' data('hicexp_diff')
#' exportJuicebox(hicexp_diff, file_name = "juiceboxAnnotations.txt") 

exportJuicebox <- function(hicexp, logfc_cutoff = 1, logcpm_cutoff = 1, p.adj_cutoff = 0.01,
                           D_cutoff = 1, file_name = "juiceboxAnnotations.txt",
                           color = "0,0,255") {
  # check that data has been compared
  if (nrow(results(hicexp)) < 1) {
    stop("Differences must be detected before making a manhattan plot.")
  }
  # make results object
  res <- results(hicexp)
  # filter
  res <- res[logCPM >= logcpm_cutoff & abs(logFC) >= logfc_cutoff & p.adj <= p.adj_cutoff & D >= D_cutoff, ]
  # reformat
  res$chr1 <- paste0('chr', res$chr)
  res$chr2 <- paste0('chr', res$chr)
  res$end1 <- res$region1 + resolution(hicexp) #- 1
  res$end2 <- res$region2 + resolution(hicexp) #- 1
  # res <- res[, c("chr1", "region1", "region2", 'chr2', 'end1', 'end2')]
  res <- res[, c("chr1", "region1", "end1", 'chr2', 'region2', 'end2')]
  res <- res[order(res$chr1, res$region1, res$region2)]
  res$color <- color
  res$comment <- "significant differentially interacting region"
  colnames(res) <- c('chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'color', 'comment')
  write.table(res, file = file_name, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
}
