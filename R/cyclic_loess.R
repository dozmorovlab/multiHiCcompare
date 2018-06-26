#' Cyclic Loess normalization for Hi-C data
#' 
#' @param hicexp A hicexp object
#' @param iterations The number of iterations (cycles) 
#'     of loess normalization to perform. Defaults to 3.
#' @param span The span for loess normalization. Defaults to
#'     NA indicating that span will be automatically calculated
#'     using generalized cross validation.
#' @param parallel Logical. Should parallel processing be used?
#' @param verbose Logical. Should messages about loess normalization
#'     be printed to the screen.
#' @param Plot Logical. Should MD plots be printed?
#' 
#' @details This function performs cyclic loess normalization
#'    on a Hi-C experiment. HiCcompare2's cyclic loess procedure
#'    is a modified version of Ballman's (2004) cyclic loess and
#'    the joint loess normalization used in the original HiCcompare.
#'    For each unique pair of samples in the hicexp object an MD plot
#'    is generated. A loess curve is fit to the MD plot and then the 
#'    fitted values are used to adjust the data. This is performed on
#'    all unique pairs and then repeated until convergence. 
#' @return A hicexp object that has been normalized. 
#' @export
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr %>%
#' @importFrom HiCcompare MD.plot1
#' @importFrom data.table rbindlist
#' @examples 
#' data("hicexp")
#' hicexp <- cyclic_loess(hicexp)

cyclic_loess <- function(hicexp, iterations = 3, span = NA, parallel = FALSE, verbose = TRUE, Plot = TRUE) {
  # check if data already normalized
  if (hicexp@normalized) {
    stop("Data has already been normalized.")
  }
  # check span input
  if (!is.na(span)) {
    if (!is.numeric(span)) {
      stop("span must be set to NA or a value between 0 and 1")
    }
    if (span <= 0 | span > 1) {
      stop("span must be set to NA or a value between 0 and 1")
    }
  }
  # check iterations input
  if (iterations < 2) {
    warning("Typically it takes about 3 iterations for cyclic loess to converge.")
  }
  if (iterations > 4) {
    warning("Typically it takes about 3 iterations for cyclic loess to converge.")
  }
  # split up data by condition and perform cyclic loess
  normalized <- .loess_condition(hicexp@hic_table, iterations = iterations, parallel = parallel, verbose = verbose, span = span)
  # sort hic_table
  normalized <- normalized[order(chr, region1, region2),]
  # put back into hicexp object
  hicexp@hic_table <- normalized
  hicexp@normalized <- TRUE
  
  # plot
  if (Plot) {
    MD.hicexp(hicexp)
  }
  
  return(hicexp)
}




# background functions
### Perform cyclic loess for a condition
.loess_condition <- function(hic_table, iterations, parallel, verbose, span) {
  # split up data by chr
  table_list <- split(hic_table, hic_table$chr)
  # plug into parallelized loess function
  if (parallel) {
    normalized <- BiocParallel::bplapply(table_list, .cloess, iterations = iterations, verbose = verbose, span = span)
  } else {
    normalized <- lapply(table_list, .cloess, iterations = iterations, verbose = verbose, span = span)
  }
  # recombine tables
  normalized <- data.table::rbindlist(normalized)
  return(normalized)
}

# perform cyclic loess on a table 
.cloess <- function(tab, iterations, verbose, span, degree = 1, loess.criterion = "gcv") {
  # make matrix of IFs
  IF_mat <- tab[, 5:(ncol(tab)), with = FALSE] %>% as.matrix()
  # make index matrix
  idx_mat <- IF_mat
  idx_mat[idx_mat != 0] <- 1
  # log the matrix
  IF_mat <- log2(IF_mat + 1)
  n <- ncol(IF_mat)
  # begin cyclic loess
  for (i in 1:iterations) {
    for(j in 1:(n-1)) {
      for (k in (j+1):n) {
        # # get rows with zeros
        # zeros1 <-  IF_mat[,k] == 0
        # zeros2 <- IF_mat[,j] == 0
        # zeros <- zeros1 | zeros2
        # M <- IF_mat[!zeros, k] - IF_mat[!zeros, j]
        # D <- tab$D[!zeros]
        M <- IF_mat[,k] - IF_mat[,j]
        D <- tab$D
        if (is.na(span)) {
          l <- .loess.as(x = D, y = M, degree = degree, 
                         criterion = loess.criterion,
                         control = loess.control(surface = "interpolate",
                                                 statistics = "approximate", trace.hat = "approximate"))
        } else {
          l <- .loess.as(x = D, y = M, degree = degree, user.span = span,
                         criterion = loess.criterion,
                         control = loess.control(surface = "interpolate",
                                                 statistics = "approximate", trace.hat = "approximate"))
        }
        # calculate gcv and AIC
        traceL <- l$trace.hat
        sigma2 <- sum(l$residuals^2)/(l$n - 1)
        aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(l$n - traceL -2)
        gcv <- l$n * sigma2/(l$n - traceL)^2
        # print the span picked by gcv
        if (verbose) {
          message("Span for loess: ", l$pars$span)
          message("GCV for loess: ", gcv)
          message("AIC for loess: ", aicc)
        }
        # adjust IFs
        # IF_mat[!zeros,j] <- IF_mat[!zeros,j] + l$fitted/2
        # IF_mat[!zeros,k] <- IF_mat[!zeros,k] - l$fitted/2
        IF_mat[,j] <- IF_mat[,j] + l$fitted/2
        IF_mat[,k] <- IF_mat[,k] - l$fitted/2
      }
    }
  }
  # anti-log IFs
  IF_mat <- (2^IF_mat) - 1
  # reset zeros
  IF_mat <- IF_mat * idx_mat
  # set negative values to 0
  IF_mat[IF_mat < 0] <- 0
  # fix any potential Infs or NaN's
  IF_mat[is.nan(IF_mat)] <- 0
  IF_mat[is.infinite(IF_mat)] <- 0
  # recombine table
  tab <- cbind(tab[, 1:4, with = FALSE], IF_mat)
  return(tab)
}



# loess with Automatic Smoothing Parameter Selection adjusted possible
# range of smoothing originally from fANCOVA package
.loess.as <- function(x, y, degree = 1, criterion = c("aicc", "gcv"),
                      family = c("gaussian",
                                 "symmetric"), user.span = NULL, plot = FALSE, ...) {
  criterion <- match.arg(criterion)
  family <- match.arg(family)
  x <- as.matrix(x)
  
  if ((ncol(x) != 1) & (ncol(x) != 2))
    stop("The predictor 'x' should be one or two dimensional!!")
  if (!is.numeric(x))
    stop("argument 'x' must be numeric!")
  if (!is.numeric(y))
    stop("argument 'y' must be numeric!")
  if (any(is.na(x)))
    stop("'x' contains missing values!")
  if (any(is.na(y)))
    stop("'y' contains missing values!")
  if (!is.null(user.span) && (length(user.span) != 1 || !is.numeric(user.span)))
    stop("argument 'user.span' must be a numerical number!")
  if (nrow(x) != length(y))
    stop("'x' and 'y' have different lengths!")
  if (length(y) < 3)
    stop("not enough observations!")
  
  data.bind <- data.frame(x = x, y = y)
  if (ncol(x) == 1) {
    names(data.bind) <- c("x", "y")
  } else {
    names(data.bind) <- c("x1", "x2", "y")
  }
  
  opt.span <- function(model, criterion = c("aicc", "gcv"), span.range = c(0.01,
                                                                           0.9)) {
    as.crit <- function(x) {
      span <- x$pars$span
      traceL <- x$trace.hat
      sigma2 <- sum(x$residuals^2)/(x$n - 1)
      aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(x$n - traceL -
                                                          2)
      gcv <- x$n * sigma2/(x$n - traceL)^2
      result <- list(span = span, aicc = aicc, gcv = gcv)
      return(result)
    }
    criterion <- match.arg(criterion)
    fn <- function(span) {
      mod <- update(model, span = span)
      as.crit(mod)[[criterion]]
    }
    result <- optimize(fn, span.range)
    return(list(span = result$minimum, criterion = result$objective))
  }
  
  if (ncol(x) == 1) {
    if (is.null(user.span)) {
      fit0 <- loess(y ~ x, degree = degree, family = family, data = data.bind,
                    ...)
      span1 <- opt.span(fit0, criterion = criterion)$span
    } else {
      span1 <- user.span
    }
    fit <- loess(y ~ x, degree = degree, span = span1, family = family,
                 data = data.bind, ...)
  } else {
    if (is.null(user.span)) {
      fit0 <- loess(y ~ x1 + x2, degree = degree, family = family,
                    data.bind, ...)
      span1 <- opt.span(fit0, criterion = criterion)$span
    } else {
      span1 <- user.span
    }
    fit <- loess(y ~ x1 + x2, degree = degree, span = span1, family = family,
                 data = data.bind, ...)
  }
  if (plot) {
    if (ncol(x) == 1) {
      m <- 100
      x.new <- seq(min(x), max(x), length.out = m)
      fit.new <- predict(fit, data.frame(x = x.new))
      plot(x, y, col = "lightgrey", xlab = "x", ylab = "m(x)", ...)
      lines(x.new, fit.new, lwd = 1.5, ...)
    } else {
      m <- 50
      x1 <- seq(min(data.bind$x1), max(data.bind$x1), len = m)
      x2 <- seq(min(data.bind$x2), max(data.bind$x2), len = m)
      x.new <- expand.grid(x1 = x1, x2 = x2)
      fit.new <- matrix(predict(fit, x.new), m, m)
      persp(x1, x2, fit.new, theta = 40, phi = 30, ticktype = "detailed",
            xlab = "x1", ylab = "x2", zlab = "y", col = "lightblue",
            expand = 0.6)
    }
  }
  return(fit)
}





