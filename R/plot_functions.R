#' Make MD plots for all combinations of a condition
#' 
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
#' @importFrom dplyr %>%
#' @export


# @import gridGraphics
# @importFrom gridExtra grid.arrange

# Plan: make a grid of MD plots for all the pairwise combinations of comparisons between replicates in a condition
# should use a modified bersion of MD.plot1 to be able to accpet new titles / labels for this function

MD.hicexp <- function(hicexp, prow = 3, pcol = 3, plot.chr = NA) {
  # check if more than one chromosome then split
  if (length(unique(hicexp@hic_table$chr)) > 1) {
    chr_table <- split(hicexp@hic_table, hicexp@hic_table$chr)
    # check if chr to plot is specified
    if (!is.na(plot.chr[1])) {
      chrs.to.plot <- which(as.numeric(names(chr_table)) %in% plot.chr)
      if (length(chrs.to.plot) < 1) {
        stop("Chr selected in plot.chr does not exist in the data")
      }
      tmp <- lapply(chr_table[chrs.to.plot], .MD.hicexp.chr, prow = prow, pcol = pcol)
    } else {
      # otherwise plot every chr
      tmp <- lapply(chr_table, .MD.hicexp.chr, prow = prow, pcol = pcol)
    }
  } else {
    # if only a single chr just plot
    .MD.hicexp.chr(hicexp@hic_table, prow = prow, pcol = pcol)
  }
  
}



.MD.hicexp.chr <- function(chr_table, prow, pcol) {
  # get all unique pairs
  samples <- 5:ncol(chr_table)
  combinations <- combn(samples, 2)
  # make M matrix
  M_matrix <- matrix(nrow = nrow(chr_table), ncol = ncol(combinations))
  plot_list <- list()
  par(mfrow = c(prow, pcol), mai = c(0.3, 0.3, 0.2, 0.1))
  for (j in 1:ncol(combinations)) {
    M_matrix[,j] <- log2( (chr_table[, combinations[1,j], with = FALSE] + 1)[[1]] / (chr_table[, combinations[2,j], with = FALSE] + 1)[[1]] )
    # make MD plot
    .MD.smooth(M_matrix[,j], chr_table$D, title = paste0('chr', chr_table$chr[1], ' ', 
                                                         'Sample ', combinations[1,j] - 4, ' vs. ', combinations[2,j] - 4), ylab = '', xlab = '')
    # grid.echo()
    # plot_list[[j]] <- grid.grab()
    # plot.new()
  }
  # cowplot::plot_grid(plotlist = plot_list[1:6], align='hv')
  .reset_par()
}


#' Plot a composite MD plot with the results of a comparison
#' @param hicexp A hicexp object.
#' @param plot.chr A specific chromosome or 
#'     set of chromosome which you want to plot.
#'     This should be a numeric value, i.e. to
#'     plot chromosome 1 set plot.chr = 1, to
#'     plot chromosomes 1 and 5 set plot.chr 
#'     = c(1, 5). Defaults to NA indicating that
#'     all chromosomes present in the hicexp
#'     will be plotted. 
#' @export

MD.composite <- function(hicexp, plot.chr = NA) {
  # check to make sure data has been compared
  if (nrow(hicexp@comparison) < 1) {
    stop("You must compare the Hi-C data first before using this plot function")
  }
  
  # check if more than one chromosome then split
  if (length(unique(hicexp@comparison$chr)) > 1) {
    chr_table <- split(hicexp@comparison, hicexp@comparison$chr)
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
    .MD.smooth(M = hicexp@comparison$logFC, D = hicexp@comparison$D, p.val = hicexp@comparison$p.adj, 
               title = "Composite MD Plot")
  }
}


.MD.smooth <- function(M, D, p.val = NA, title = 'MD Plot', ylab = 'M', xlab = 'Distance') {
  # smooth scatter version
  smoothScatter(D, M, xlab = xlab, ylab = ylab, main = title, cex.main = 0.85)
  abline(h = 0)
  if (!is.na(p.val[1])) {
    p0.01 <- which(p.val < 0.01)
    p0.05 <- which(p.val >= 0.01 & p.val < 0.05)
    points(D[p0.01], M[p0.01], col = "red", pch = 20)
    points(D[p0.05], M[p0.05], col = 'yellow', pch = 20)
    legend('bottomright', legend = c('P < 0.01', 'P < 0.05'), fill = c('red', 'yellow'), bty = 'n', horiz = TRUE)
  }
}



# function to reset par
.reset_par <- function(){
  op <- structure(list(xlog = FALSE, ylog = FALSE, adj = 0.5, ann = TRUE,
                       ask = FALSE, bg = "transparent", bty = "o", cex = 1, cex.axis = 1,
                       cex.lab = 1, cex.main = 1.2, cex.sub = 1, col = "black",
                       col.axis = "black", col.lab = "black", col.main = "black",
                       col.sub = "black", crt = 0, err = 0L, family = "", fg = "black",
                       fig = c(0, 1, 0, 1), fin = c(6.99999895833333, 6.99999895833333
                       ), font = 1L, font.axis = 1L, font.lab = 1L, font.main = 2L,
                       font.sub = 1L, lab = c(5L, 5L, 7L), las = 0L, lend = "round",
                       lheight = 1, ljoin = "round", lmitre = 10, lty = "solid",
                       lwd = 1, mai = c(1.02, 0.82, 0.82, 0.42), mar = c(5.1, 4.1,
                                                                         4.1, 2.1), mex = 1, mfcol = c(1L, 1L), mfg = c(1L, 1L, 1L,
                                                                                                                        1L), mfrow = c(1L, 1L), mgp = c(3, 1, 0), mkh = 0.001, new = FALSE,
                       oma = c(0, 0, 0, 0), omd = c(0, 1, 0, 1), omi = c(0, 0, 0,
                                                                         0), pch = 1L, pin = c(5.75999895833333, 5.15999895833333),
                       plt = c(0.117142874574832, 0.939999991071427, 0.145714307397962,
                               0.882857125425167), ps = 12L, pty = "m", smo = 1, srt = 0,
                       tck = NA_real_, tcl = -0.5, usr = c(0.568, 1.432, 0.568,
                                                           1.432), xaxp = c(0.6, 1.4, 4), xaxs = "r", xaxt = "s", xpd = FALSE,
                       yaxp = c(0.6, 1.4, 4), yaxs = "r", yaxt = "s", ylbias = 0.2), .Names = c("xlog",
                                                                                                "ylog", "adj", "ann", "ask", "bg", "bty", "cex", "cex.axis",
                                                                                                "cex.lab", "cex.main", "cex.sub", "col", "col.axis", "col.lab",
                                                                                                "col.main", "col.sub", "crt", "err", "family", "fg", "fig", "fin",
                                                                                                "font", "font.axis", "font.lab", "font.main", "font.sub", "lab",
                                                                                                "las", "lend", "lheight", "ljoin", "lmitre", "lty", "lwd", "mai",
                                                                                                "mar", "mex", "mfcol", "mfg", "mfrow", "mgp", "mkh", "new", "oma",
                                                                                                "omd", "omi", "pch", "pin", "plt", "ps", "pty", "smo", "srt",
                                                                                                "tck", "tcl", "usr", "xaxp", "xaxs", "xaxt", "xpd", "yaxp", "yaxs",
                                                                                                "yaxt", "ylbias"))
  par(op)
}


# 
# MD.hicexp <- function(hicexp) {
#  # get list of plots by condition
#   plot_by_condition <- lapply(hicexp@hic_tables, .MD_condition)
#   # arrange into grid
#   p <- unlist(plot_by_condition)
#   grid.arrange(p)
# }
# 
# 
# .MD_condition <- function(hic_table) {
#   # get all unique pairs of samples which need to be compared
#   samples <- 5:(ncol(hic_table))
#   combinations <- combn(samples, 2)
#   # make matrix of M values, each col corresponding to combination of samples specified in combinations object
#   M_matrix <- matrix(nrow = nrow(hic_table), ncol = ncol(combinations))
#   plot_list <- list()
#   for (j in 1:ncol(combinations)) {
#     M_matrix[,j] <- log2( (hic_table[, combinations[1,j], with = FALSE] + 1)[[1]] / (hic_table[, combinations[2,j], with = FALSE] + 1)[[1]] )
#     # make plots
#     HiCcompare::MD.plot2(M_matrix[,j], hic_table$D)
#     grid.echo()
#     plot_list[[j]] <- grid.grab()
#     plot.new()
#   }
#   return(plot_list)
# }
