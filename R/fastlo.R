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
  table_by_chr <- split(hic_table, hic_table$chr)
  
}