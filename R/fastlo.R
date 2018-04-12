#' Perform fast loess normalization on a Hi-C experiment
#' 
#' 

fastlo <- function(hicexp, iterations = 3, span = NA, parallel = FALSE, verbose = TRUE, Plot = TRUE) {
  # input conditions to fastlo
  normalized <- lapply(hicexp@hic_tables, .fastlo_condition, iterations = iterations, span = span, parallel = parallel, verbose = verbose, Plot = Plot)
}

# bacground functions

# perform fastlo on a condition
.fastlo_condition <- function(hic_table, iterations, span, parallel, verbose, Plot) {
  # split up data by chr
  table_list <- split(hic_table, hic_table$chr)
  # for each chr create list of distance matrices
  table_list <- lapply(table_list, .get_dist_tables)
  # combine list of lists into single list of tables
  table_list <- do.call(c, table_list)
}

.get_dist_tables <- function(chr_table) {
  all_dist <- sort(unique(chr_table$D))
  dist_85 <- ceiling(0.85 * length(all_dist))
  table_by_dist <- list()
  idx <- 1
  for(i in all_dist[1:dist_85]) {
    table_by_dist[[idx]] <- chr_table[D == i,]
    idx <- idx + 1
  }
  table_by_dist[[idx]] <- chr_table[D > all_dist[dist_85],]
  return(table_by_dist)
}
