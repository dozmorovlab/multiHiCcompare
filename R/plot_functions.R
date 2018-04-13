#' Make MD plots for all combinations of a condition
#' 
#' 
#' @import gridGraphics
#' @importFrom gridExtra grid.arrange

# Plan: make a grid of MD plots for all the pairwise combinations of comparisons between replicates in a condition
# should use a modified bersion of MD.plot1 to be able to accpet new titles / labels for this function

MD.hicexp <- function(hicexp) {
 # get list of plots by condition
  plot_by_condition <- lapply(hicexp@hic_tables, .MD_condition)
  # arrange into grid
  p <- unlist(plot_by_condition)
  grid.arrange(p)
}


.MD_condition <- function(hic_table) {
  # get all unique pairs of samples which need to be compared
  samples <- 5:(ncol(hic_table))
  combinations <- combn(samples, 2)
  # make matrix of M values, each col corresponding to combination of samples specified in combinations object
  M_matrix <- matrix(nrow = nrow(hic_table), ncol = ncol(combinations))
  plot_list <- list()
  for (j in 1:ncol(combinations)) {
    M_matrix[,j] <- log2( (hic_table[, combinations[1,j], with = FALSE] + 1)[[1]] / (hic_table[, combinations[2,j], with = FALSE] + 1)[[1]] )
    # make plots
    HiCcompare::MD.plot2(M_matrix[,j], hic_table$D)
    grid.echo()
    plot_list[[j]] <- grid.grab()
    plot.new()
  }
  return(plot_list)
}
