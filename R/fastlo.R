#' Perform fast loess normalization on a Hi-C experiment
#' 
#' @param hicexp A hicexp object
#' @param iterations The number of iterations (cycles) 
#'     for fastlo to proceed through.
#' @param span The span of loess fitting. Defaults to 0.7. To
#'     automatically calculate a span using the GCV set 
#'     span = NA. However note that setting span = NA will
#'     significantly slow down the normalization process.
#' @param parallel Logical. Should parallel processing be
#'     used?
#' 
#' @details This function performs the fast loess (fastlo)
#'    normalization procedure on an hicexp object. 
#' 
#' @export
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr %>%
#' @importFrom data.table rbindlist

fastlo <- function(hicexp, iterations = 3, span = 0.7, parallel = FALSE, verbose = TRUE, Plot = TRUE) {
  # input conditions to fastlo
  normalized <- lapply(hicexp@hic_tables, .fastlo_condition, iterations = iterations, span = span, parallel = parallel, verbose = verbose, Plot = Plot)
  # put back into hicexp object
  hicexp@hic_tables <- normalized
  hicexp@normalized <- TRUE
  return(hicexp)
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
  # apply fastlo to list
  if (parallel) {
    normalized <- BiocParallel::bplapply(table_list, .fastlo, span = span, iterations = iterations)
  } else {
    normalized <- lapply(table_list, .fastlo, span = span, iterations = iterations)
  }
  # recombine list of tables
  hic_table <- rbindlist(normalized)
  return(hic_table)
}

# make list of tables by distance
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

# perform fastlo on table
.fastlo <- function(tab, span, iterations, loess.criterion = "gcv", degree = 1) {
  # make matrix of IFs
  IF_mat <- tab[, 5:(ncol(tab)), with = FALSE] %>% as.matrix()
  # log the IF matrix
  IF_mat <- log2(IF_mat + 1)
  n <- ncol(IF_mat)
  # outer loop for number of cycles
  for (i in 1:iterations) {
    # calculate rowmeans (A) 
    A <- rowMeans(IF_mat, na.rm = TRUE)
    for (j in 1:n) {
      M <- IF_mat[,j] - A
      # fit loess curve
      if (is.na(span)) {
        f <- .loess.as(x = A, y = M, degree = degree, 
                       criterion = loess.criterion,
                       control = loess.control(surface = "interpolate",
                                               statistics = "approximate", trace.hat = "approximate"))
      } else {
        f <- .loess.as(x = A, y = M, degree = degree, user.span = span,
                       criterion = loess.criterion,
                       control = loess.control(surface = "interpolate",
                                               statistics = "approximate", trace.hat = "approximate"))
      }
      # adjust 
      IF_mat[,j] <- IF_mat[,j] - f$fitted
    }
  }
  # anti-log the IF_mat
  IF_mat <- (2^IF_mat) - 1
  # deal with negative values after normalization
  IF_mat[IF_mat < 0] <- 0
  # rebuild table
  tab <- cbind(tab[, 1:4, with = FALSE], IF_mat)
  return(tab)
}
