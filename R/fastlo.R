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
#' @param verbose Logical, should messages about
#'     the normalization be printed?
#' @param max.pool The proportion of unit distances after
#'     which all further distances will be pooled. Distances
#'     before this value will be progressively pooled and
#'     any distances after this value will be combined
#'     into a single pool. Defaults to 0.7. Warning: do
#'     not adjust this value from the default unless you
#'     are getting errors related to the lfproc function
#'     or due to sparsity in fastlo normalization. If these
#'     errors occur it is due to either sparsity or low 
#'     variance and max.pool will need to be lowered; 
#'     typically to 0.5 or 0.6. 
#' 
#' @details This function performs the fast loess (fastlo)
#'    normalization procedure on a hicexp object. the fast linear 
#'    loess ("fastlo") method 
#'    of Ballman (2004) that is adapted to Hi-C data on a per-distance basis.
#'     To perform 
#'    "fastlo" on Hi-C data we first split the data into p pooled matrices. 
#'    The "progressive pooling" is used to split up the Hi-C matrix by unit 
#'    distance. Fastlo
#'    is then performed on the MA plots for each distance pool. 
#'    See Stansfield et al (2018)
#'    for full description. 
#'   
#' @return A hicexp object that is normalized.
#' @export
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr left_join
#' @importFrom data.table rbindlist
#' @examples 
#' data("hicexp2")
#' hicexp2 <- fastlo(hicexp2)

fastlo <- function(hicexp, iterations = 3, span = 0.7, parallel = FALSE, 
                   verbose = FALSE, max.pool = 0.7) {
  # check if data already normalized
  if (normalized(hicexp)) {
    stop("Data has already been normalized.")
  }
  # check span input
  if (!is.na(span)) {
    if (!is.numeric(span) || span <= 0 || span > 1) {
      stop("span must be set to NA or a value between 0 and 1")
    }
  }
  # check iterations input
  if (iterations != 3L) {
    warning("Typically it takes about 3 iterations for 
            cyclic loess to converge.")
  }

  # input conditions to fastlo
  normalized <- .fastlo_condition(hic_table(hicexp), 
                                  iterations = iterations, span = span, 
                                  parallel = parallel, verbose = verbose, 
                                  max.pool = max.pool)
  # sort hic_table
  normalized <- normalized[order(chr, region1, region2),]
  # put back into hicexp object
  slot(hicexp, "hic_table") <- normalized
  slot(hicexp, "normalized") <- TRUE

  
  return(hicexp)
}

# background functions

# perform fastlo on a condition
.fastlo_condition <- function(hic_table, iterations, span, 
                              parallel, verbose, max.pool) {
  # split up data by chr
  table_list <- split(hic_table, hic_table$chr)
  # for each chr create list of distance matrices
  table_list <- lapply(table_list, .get_dist_tables, max.pool = max.pool)
  # combine list of lists into single list of tables
  table_list <- do.call(c, table_list)
  # apply fastlo to list
  normalized <- smartApply(parallel, table_list, .fastlo, span = span, 
                                      iterations = iterations)
  
  # recombine list of tables
  hic_table <- rbindlist(normalized)
  return(hic_table)
}


## This is the new version where we use progressive pooling of distances
# make list of tables by distance pool
.get_dist_tables <- function(chr_table, max.pool) {
  # check for correct max.pool input
  if (!is.numeric(max.pool)) {
    stop("max.pool must be a numeric value between 0 and 1")
  }
  if (max.pool < 0 || max.pool > 1) {
    stop("max.pool must be between 0 and 1")
  }
  if (max.pool < 0.5) {
    warning("Setting max.pool < 0.5 may cause issues")
  }
  if (max.pool > 0.7) {
    warning("Setting max.pool > 0.7 may cause issues")
  }
  D <- sort(unique(chr_table$D)) 
  # use triangular number series to solve for number of pools
  x <- length(D)
  n <- (sqrt((8 * x) + 1) - 1) / 2
  n <- ceiling(n)
  pools <- rep(D[1:n], 1:n)[1:x]
  # add the pools column to the data.table corresponding to the distances
  id.df <- as.data.frame(cbind(D, pools))
  # get maximum distance
  dist_max <- ceiling(max.pool * nrow(id.df))
  # combine pools for everything at or above maximum distance
  pool_max <- id.df[dist_max,2] 
  id.df[id.df$pools >= pool_max,2] <- pool_max
  table_by_dist <- data.table::as.data.table(left_join(chr_table, id.df, by = c("D" = "D")))
  # split up tables by pool
  table_by_dist <- split(table_by_dist, table_by_dist$pools)
  # drop pools column
  table_by_dist <- lapply(table_by_dist, function(x) x[, pools := NULL])
  # check number of rows per table
  table_by_dist <- .check_tables(table_by_dist)
  return(table_by_dist)
}

# check to make sure number of rows for each table is not less than 500
.check_tables <- function(table_by_dist) {
  min.row <- 500
  # get row numbers
  row.nums <- lapply(table_by_dist, nrow)
  # find which tables have number of rows < min.row
  small_D <- which(row.nums < min.row)
  # check if any tables have less than min.row
  if (length(small_D) > 0) {
    while (length(small_D) > 0) {
      i <- small_D[1]
      # check if first distance has less than min.row
      if (table_by_dist[[i]]$D[1] == 0) {
        new_table <- rbind(table_by_dist[[i]], table_by_dist[[i+1]])
        table_by_dist[[i]] <- new_table
        table_by_dist[[i+1]] <- NA
      } else {
        # otherwise combine with previous table
        new_table <- rbind(table_by_dist[[i-1]], table_by_dist[[i]])
        table_by_dist[[i-1]] <- new_table
        table_by_dist[[i]] <- NA
      }
      # remove NAs
      table_by_dist <- table_by_dist[!is.na(table_by_dist)]
      # update row numbers
      row.nums <- lapply(table_by_dist, nrow)
      # update which tables have number of rows < min.row
      small_D <- which(row.nums < min.row)
    }
    return(table_by_dist)
  } else {
    return(table_by_dist)
  }
}



# perform fastlo on table
.fastlo <- function(tab, span, iterations, loess.criterion = "gcv", degree = 1) {
  # make matrix of IFs
  IF_mat <- as.matrix(tab[, -c("chr", "region1", "region2", "D"), with = FALSE])
  # make indicator matrix
  idx_mat <- IF_mat
  idx_mat[idx_mat != 0] <- 1
  # log the IF matrix
  IF_mat <- log2(IF_mat + 1)
  n <- ncol(IF_mat)
  # outer loop for number of cycles
  for (i in seq(from = 1, to = iterations)) {
    # calculate rowmeans (A) 
    A <- rowMeans(IF_mat, na.rm = TRUE)
    for (j in seq(from = 1, to = n)) {
      M <- IF_mat[,j] - A
      # fit loess curve
      if (is.na(span)) {
        f <- .loess.as(x = A, y = M, degree = degree, 
                       criterion = loess.criterion,
                       control = loess.control(surface = "interpolate",
                                               statistics = "approximate",
                                               trace.hat = "approximate"))
      } else {
        f <- .loess.as(x = A, y = M, degree = degree, user.span = span,
                       criterion = loess.criterion,
                       control = loess.control(surface = "interpolate",
                                               statistics = "approximate",
                                               trace.hat = "approximate"))
      }
      # adjust 
      IF_mat[,j] <- IF_mat[,j] - f$fitted
    }
  }
  # anti-log the IF_mat
  IF_mat <- (2^IF_mat) - 1
  # revert zeros
  IF_mat <- IF_mat * idx_mat
  # deal with negative values after normalization
  IF_mat[IF_mat < 0] <- 0
  # fix any potential Infs or NaN's
  IF_mat[is.nan(IF_mat)] <- 0
  IF_mat[is.infinite(IF_mat)] <- 0
  # rebuild table
  tab <- cbind(tab[, c('chr', 'region1', 'region2', 'D'), with = FALSE], IF_mat)
  return(tab)
}
