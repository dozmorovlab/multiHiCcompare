#' Make MD plots for all combinations of a condition
#' 
#' @importFrom graphics abline legend lines par persp
#'     points smoothScatter
#' @param hicexp A hicexp object.
#' @param prow The number of rows to use for the 
#'    grid of MD plots. Defaults to 3.
#' @param pcol The number of columns to use for
#'     the grid of MD plots. Defaults to 3.
#' @param plot.chr A specific chromosome or 
#'     set of chromosome which you want to plot.
#'     This should be a numeric value, i.e. to
#'     plot chromosome 1 set plot.chr = 1, to
#'     plot chromosomes 1 and 5 set plot.chr 
#'     = c(1, 5). Defaults to NA indicating that
#'     all chromosomes present in the hicexp
#'     will be plotted. 
#'
#' @return A set of MD plots.
#' 
#'
#' @export
#' @examples
#' data("hicexp2")
#' MD_hicexp(hicexp2)




MD_hicexp <- function(hicexp, prow = 3, pcol = 3, plot.chr = NA) {
  # check if more than one chromosome then split
  if (length(unique(hic_table(hicexp)$chr)) > 1) {
    chr_table <- split(hic_table(hicexp), hic_table(hicexp)$chr)
    # check if chr to plot is specified
    if (!is.na(plot.chr[1])) {
      chrs.to.plot <- which(as.numeric(names(chr_table)) %in% plot.chr)
      if (length(chrs.to.plot) < 1) {
        stop("Chr selected in plot.chr does not exist in the data")
      }
      tmp <- lapply(chr_table[chrs.to.plot], .MD.hicexp.chr, prow = prow, 
                    pcol = pcol)
    } else {
      # otherwise plot every chr
      tmp <- lapply(chr_table, .MD.hicexp.chr, prow = prow, pcol = pcol)
    }
  } else {
    # if only a single chr just plot
    .MD.hicexp.chr(hic_table(hicexp), prow = prow, pcol = pcol)
  }
  
}



.MD.hicexp.chr <- function(chr_table, prow, pcol) {
  # save par
  old_par <- par()
  # get all unique pairs
  samples <- 5:ncol(chr_table)
  combinations <- combn(samples, 2)
  # make M matrix
  M_matrix <- matrix(nrow = nrow(chr_table), ncol = ncol(combinations))
  plot_list <- list()
  par(mfrow = c(prow, pcol), mai = c(0.3, 0.3, 0.2, 0.1))
  for (j in 1:ncol(combinations)) {
    M_matrix[,j] <- log2( (chr_table[, combinations[1,j], with = FALSE] + 1)[[1]] / 
                            (chr_table[, combinations[2,j], with = FALSE] + 1)[[1]] )
    # make MD plot
    .MD.smooth(M_matrix[,j], chr_table$D, title = paste0('chr', chr_table$chr[1], 
                                                         ' ', 
                                                         'Sample ', 
                                                         combinations[1,j] - 4,
                                                         ' vs. ', combinations[2,j] - 4),
                                                          ylab = '', xlab = '')
  }
  # reset par
  suppressWarnings(par(old_par))
}


#' Plot a composite MD plot with the results of a comparison
#' @param hicexp A hicexp object which has 
#'     had a multiHiCcompare comparison step performed on it.
#' @param plot.chr A specific chromosome or 
#'     set of chromosome which you want to plot.
#'     This should be a numeric value, i.e. to
#'     plot chromosome 1 set plot.chr = 1, to
#'     plot chromosomes 1 and 5 set plot.chr 
#'     = c(1, 5). Defaults to NA indicating that
#'     all chromosomes present in the hicexp
#'     will be plotted. 
#' @return An MD plot
#' @examples 
#' data("hicexp_diff")
#' MD_composite(hicexp_diff)
#' @export


MD_composite <- function(hicexp, plot.chr = NA) {
  # check to make sure data has been compared
  if (nrow(results(hicexp)) < 1) {
    stop("You must compare the Hi-C data first before using 
         this plot function")
  }
  
  # check if more than one chromosome then split
  if (length(unique(results(hicexp)$chr)) > 1) {
    chr_table <- split(results(hicexp), results(hicexp)$chr)
    # check if chr to plot is specified
    if (!is.na(plot.chr[1])) {
      chrs.to.plot <- which(as.numeric(names(chr_table)) %in% plot.chr)
      if (length(chrs.to.plot) < 1) {
        stop("Chr selected in plot.chr does not exist in the data")
      }
      tmp <- lapply(chr_table[chrs.to.plot], function(x) {
        .MD.smooth(x$logFC, x$D, x$p.adj, title = paste0('chr', x$chr[1]))
      })
    } else {
      # otherwise plot every chr
      tmp <- lapply(chr_table, function(x) {
        .MD.smooth(x$logFC, x$D, x$p.adj, title = paste0('chr', x$chr[1]))
      })
    }
  } else {
    # if only a single chr just plot
    .MD.smooth(M = results(hicexp)$logFC, D = results(hicexp)$D, 
               p.val = results(hicexp)$p.adj, 
               title = "Composite MD Plot")
  }
}


.MD.smooth <- function(M, D, p.val = NA, title = 'MD Plot', ylab = 'M', 
                       xlab = 'Distance') {
  # smooth scatter version
  smoothScatter(D, M, xlab = xlab, ylab = ylab, main = title, cex.main = 0.85)
  abline(h = 0)
  if (!is.na(p.val[1])) {
    p0.001 <- which(p.val < 0.001)
    p0.05 <- which(p.val >= 0.001 & p.val < 0.05)
    points(D[p0.001], M[p0.001], col = "red", pch = 20)
    points(D[p0.05], M[p0.05], col = 'yellow', pch = 20)
    legend('bottomright', legend = c('P < 0.001', 'P < 0.05'), 
           fill = c('red', 'yellow'), bty = 'n', horiz = TRUE)
  }
}





#' Function to visualize p-values from multiHiCcompare results
#' 
#' @param hicexp A hicexp object that has been
#'     normalized and has had differences detected.
#' @param alpha The alpha level at which you will 
#'     call a p-value significant. If this is set to
#'     a numeric value then any p-values >= alpha will
#'     be set to 1 for the visualization in the heatmap.
#'     Defaults to NA for visualization of all p-values.
#' @param chr The numeric value for the chromosome that 
#'     you want to plot. Set to 0 to plot all chromosomes
#'     in the dataset.
#' @details The goal of this function is to visualize
#'     where in the Hi-C matrix the differences are
#'     occuring between two experimental conditions.
#'     The function will produce a heatmap of the
#'     -log10(p-values) * sign(logFC) 
#'     to visualize where the
#'     significant differences between the datasets
#'     are occuring on the genome. 
#' @return A heatmap
#' @examples 
#' data("hicexp_diff")
#' pval_heatmap(hicexp_diff, chr = 22)
#' @importFrom pheatmap pheatmap
#' @export

pval_heatmap <- function(hicexp, alpha = NA, chr = 0) {
  # check input
  if (nrow(results(hicexp)) < 1) {
    stop("You must compare the hicexp first.")
  }
  if (!is.numeric(chr)) {
    stop("chr should be a numeric value")
  }
  if (chr != 0 & sum(chr == as.numeric(unique(results(hicexp)$chr))) < 1) {
    stop("The value of chr selected does not appear in the hicexp")
  }
  
  # if chr = 0 split data up by chr
  if (chr == 0) {
    chr.list <- split(results(hicexp), results(hicexp)$chr)
  } else {
    # otherwise subset data to just the selected chr
    chr.list <- list(subset(results(hicexp)[chr == chr,]))
  }
  
  # convert to sparse matrix
  m <- lapply(chr.list, function(x) {
    new.table <- cbind(x$region1, x$region2, x$p.adj)
    return(new.table)
  })
  
  # convert to full matrix
  m <- lapply(m, HiCcompare::sparse2full) 
  
  # get fold change
  fc<- lapply(chr.list, function(x) {
    new.table <- cbind(x$region1, x$region2, x$logFC)
    return(new.table)
  })
  fc <- lapply(fc, HiCcompare::sparse2full)
  # remove non significant values from matrix if alpha is set to a value
  if (!is.na(alpha)) {
    for (i in seq_len(length(m))) {
      m[[i]][m[[i]] >= alpha] <- 1
    }
  }
  # plot heatmap
  mapply(function(m, fc) {
    pheatmap::pheatmap(-log10(m) * sign(fc), cluster_rows = FALSE, 
                       cluster_cols = FALSE, show_rownames = FALSE, 
                       show_colnames = FALSE)
    }, 
    m, fc, SIMPLIFY = FALSE)
  
  
}
