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
#'    on a Hi-C experiment. 
#' @export
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr %>%
#' @importFrom HiCcompare MD.plot1
#' @importFrom data.table rbindlist

cyclic_loess <- function(hicexp, iterations = 3, span = NA, parallel = FALSE, verbose = TRUE, Plot = TRUE) {
  # split up data by condition and perform cyclic loess
  normalized <- .loess_condition(hicexp@hic_table, iterations = iterations, parallel = parallel, verbose = verbose, span = span, Plot = Plot)
  # put back into hicexp object
  hicexp@hic_table <- normalized
  hicexp@normalized <- TRUE
  return(hicexp)
}




# background functions
### Perform cyclic loess for a condition
.loess_condition <- function(hic_table, iterations, parallel, verbose, span, Plot) {
  # cycles <- 0
  for(iter in 1:iterations) {
    # repeat {
    # get all unique pairs of samples which need to be compared
    samples <- 5:(ncol(hic_table))
    combinations <- combn(samples, 2)
    # make matrix of M values, each col corresponding to combination of samples specified in combinations object
    M_matrix <- matrix(nrow = nrow(hic_table), ncol = ncol(combinations))
    for (j in 1:ncol(combinations)) {
      M_matrix[,j] <- log2( (hic_table[, combinations[1,j], with = FALSE] + 1)[[1]] / (hic_table[, combinations[2,j], with = FALSE] + 1)[[1]] )
    }
    # cbind M_matrix to hic table
    hic_table <- cbind(hic_table, M_matrix)
    # split up data by chr
    table_by_chr <- split(hic_table, hic_table$chr)
    # create MD list; this splits up the data into a list containing all data that will need to be loess normalized so that it can put into parallel apply statement
    MD_list <- list()
    idx <- 1
    for (k in 1:length(table_by_chr)) {
      for(l in (5 + length(samples)):(ncol(table_by_chr[[k]]))) {
        MD_list[[idx]] <- cbind(D = table_by_chr[[k]]$D, M = table_by_chr[[k]][, l, with = FALSE])
        colnames(MD_list[[idx]]) <- c("D", "M")
        idx <- idx + 1
      }
    }
    # input MD_list into loess function
    # return correction factor to correct IFs
    if (parallel) {
      mhat_list <- BiocParallel::bplapply(MD_list, .MD_loess, verbose = verbose, span = span)
    } else {
      mhat_list <- lapply(MD_list, .MD_loess, verbose = verbose, span = span)
    }

    # go through samples and apply loess correction to IFs
    idx <- 1
    for (k in 1:length(table_by_chr)) {
      for(l in 1:ncol(combinations)) {
        # adjust the IFs
        IF1 <- table_by_chr[[k]][, combinations[1,l], with = FALSE]
        IF2 <- table_by_chr[[k]][, combinations[2,l], with = FALSE]
        IF1 <- 2^(log2(IF1 + 1) - mhat_list[[idx]]) %>% as.vector()
        IF2 <- 2^(log2(IF2 + 1) + mhat_list[[idx]]) %>% as.vector()
        ## plot MD plot
        if (Plot) HiCcompare::MD.plot1((table_by_chr[[k]][, (4 + length(samples) + l), with = FALSE])[[1]], table_by_chr[[k]]$D, mhat_list[[idx]]*2)
        ##
        # update table
        table_by_chr[[k]][, combinations[1,l]] <- IF1
        table_by_chr[[k]][, combinations[2,l]] <- IF2
        idx <- idx + 1
      }
    }
    
    # recombine table_by_chr
    hic_table <- rbindlist(table_by_chr)[, 1:(ncol(hic_table) - ncol(combinations)), with = FALSE] # drop the M columns
    # check convergence criteria
    # cycles <- cycles + 1
    # print( quantile(abs(unlist(mhat_list))))
    # if ( max(abs(unlist(mhat_list))) < converge ) {
    #   message(paste0("\n Cyclic loess converged in ", cycles, " cycles. \n"))
    #   break()
    # }
    # if (cycles > 5) {
    #   warning("Cyclic loess performed over 5 cycles without convergence; ending loop")
    #   break()
    # }
  }
  return(hic_table)
}



# loess on MD_list object
.MD_loess <- function(MD_item, verbose, degree = 1, span, loess.criterion = "gcv") {
  if (is.na(span)) {
    l <- .loess.as(x = MD_item$D, y = MD_item$M, degree = degree, 
                   criterion = loess.criterion,
                   control = loess.control(surface = "interpolate",
                                           statistics = "approximate", trace.hat = "approximate"))
  } else {
    l <- .loess.as(x = MD_item$D, y = MD_item$M, degree = degree, user.span = span,
                   criterion = loess.criterion,
                   control = loess.control(surface = "interpolate",
                                           statistics = "approximate", trace.hat = "approximate"))
  }
  
  # calculate gcv and AIC
  traceL <- l$trace.hat
  sigma2 <- sum(l$residuals^2)/(l$n - 1)
  aicc <- log(sigma2) + 1 + 2 * (2 * (traceL + 1))/(l$n - traceL -
                                                      2)
  gcv <- l$n * sigma2/(l$n - traceL)^2
  # print the span picked by gcv
  if (verbose) {
    message("Span for loess: ", l$pars$span)
    message("GCV for loess: ", gcv)
    message("AIC for loess: ", aicc)
  }
  # get the correction factor
  mc <- predict(l, MD_item$D)
  mhat <- mc/2
  return(mhat)
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





