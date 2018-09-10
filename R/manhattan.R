#' Manhattan plot function for results of multiHiCcompare
#' 
#' @param hicexp A hicexp object that has had differences
#'     detected
#' @param method string denoting the p-value method to
#'     use for plotting. Options are "standard", "fisher",
#'     "stouffer", and count. "standard" plots a manhattan plot
#'     using all individual p-values. "fisher" uses 
#'     Fisher's method for combining p-values to combine
#'     the p-values for each region which are then plotted.
#'     "stouffer" uses the Stouffer-Liptak method for 
#'     combining p-values for each region which are then 
#'     plotted. "count" produces a plot where the y-values
#'     are the -log10(1 / S) where S is the number of times
#'     the region was found significant. 
#' @param return_df Logical, should the data.frame used to
#'     generate the plot be returned?
#' @param alpha The p-value cut off to be used for 
#'     calling an interaction significant. This is
#'     only used if method = 'count'. Defaults to 
#'     0.05.
#' @param plot.chr A numeric value indicating a specific
#'     chromosome number to subset the plot to. Defaults
#'     to NA indicating that all chromosomes will be plotted.
#' @details This function is used to create a manhattan
#'     plot for the significance of all genomic regions 
#'     in the dataset. These correspond to the rows (or columns)
#'     of the upper triangle of the full Hi-C matrix. Each genomic
#'     region of the Hi-C dataset has multiple interactions it
#'     is involved in and the significance of all of these can 
#'     be visualized with \code{method = "standard"}. 
#'     Alternatively the p-values for all these interactions
#'     can be combined using either Fisher's method or the
#'     Stouffer-Liptak method of combining p-values. Additionally
#'     the "count" option will plot based on the number of times
#'     each region was found to be involved in a signficantly 
#'     different interaction. The 
#'     manhattan plot can be used to identify "hotspot"
#'     regions of the genome where major differences
#'     seem to be located based on the results of a multiHiCcompare
#'     analysis.
#' @return A manhattan plot and optionally the data.frame used
#'     to generate the manhattan plot.
#' @importFrom qqman manhattan
#' @importFrom metap sumz sumlog
#' @export
#' @examples
#' data("hicexp_diff")
#' manhattan_hicexp(hicexp_diff, method = "fisher")


manhattan_hicexp <- function(hicexp, method = "standard", return_df = FALSE, 
                             alpha = 0.05, plot.chr = NA) {
  # check input
  method <- match.arg(method, c("standard", "fisher", "stouffer", "count"), 
                      several.ok = FALSE)
  if (!is.na(plot.chr)) {
    if (!is.numeric(plot.chr)) {
      stop("plot.chr must be either NA or a numeric value.")
    }
    if (!(plot.chr %in% unique(results(hicexp)$chr))) {
      stop("The chr chosen as a subset must exist in the data object")
    }
  }
  # check that data has been compared
  if (nrow(results(hicexp)) < 1) {
    stop("Differences must be detected before making a manhattan plot.")
  }
  
  if (method == "standard") {
    # make data.frame for plotting
    man.df <- data.frame(BP = c(results(hicexp)$region1, 
                                results(hicexp)$region2),
                         CHR = as.numeric(c(results(hicexp)$chr, 
                                            results(hicexp)$chr)), 
                         P = c(results(hicexp)$p.adj, 
                               results(hicexp)$p.adj))
    # subset by chr is option is not NA
    if (!is.na(plot.chr)) {
      man.df <- man.df[man.df$CHR == plot.chr,]
    }
    # plot
    suppressWarnings(qqman::manhattan(man.df, suggestiveline = FALSE, genomewideline = FALSE))
  }
  
  if (method == "fisher") {
    # make aggregate p-value for regions
    regions <- c(paste0(results(hicexp)$chr, ':', results(hicexp)$region1),
                 paste0(results(hicexp)$chr, ':', results(hicexp)$region2))
    p.values <- c(results(hicexp)$p.adj, results(hicexp)$p.adj)
    
    ## Fisher method
    # fisher_aggregate <- aggregate(p.values, by = list(regions), 
    #                               FUN = function(p) {
    #   combined <- -2 * sum(log(p))
    #   c.pval <- pchisq(combined, df = 2 * length(p))
    #   return(c.pval)
    # })
    
    fisher_aggregate <- aggregate(p.values, 
                                  by = list(regions), 
                                  FUN = function(x) {
                                    if (length(x) > 1) {
                                      metap::sumlog(x)$p
                                    } else {
                                      x
                                    }
                                    })
    
    fisher_aggregate <- cbind(read.table(text = fisher_aggregate$Group.1, 
                                         sep = ":"), fisher_aggregate$x)
    colnames(fisher_aggregate) <- c("CHR", "BP", "P")
    fisher_aggregate$P[fisher_aggregate$P == 0] <- 10^-100
    
    # subset by chr is option is not NA
    if (!is.na(plot.chr)) {
      fisher_aggregate <- fisher_aggregate[fisher_aggregate$CHR == plot.chr,]
    }
    # plot combined p-value manahttan plots
    suppressWarnings(qqman::manhattan(fisher_aggregate, suggestiveline = FALSE, genomewideline = FALSE))
    
    man.df <- fisher_aggregate
  }
  
  if (method == "stouffer") {
    # make aggregate p-value for regions
    regions <- c(paste0(results(hicexp)$chr, ':', results(hicexp)$region1),
                 paste0(results(hicexp)$chr, ':', results(hicexp)$region2))
    p.values <- c(results(hicexp)$p.adj, results(hicexp)$p.adj)
    p.values[p.values == 1] <- 0.99999 # change pvalues from 1 so sumz works correctly
    
    
    stouffer_liptak_aggregate <- aggregate(p.values, 
                                           by = list(regions), 
                                           FUN = function(x) {
                                             if (length(x) > 1) {
                                               suppressWarnings(metap::sumz(x)$p) 
                                             } else {
                                               x
                                             }
                                            })
    
    stouffer_liptak_aggregate <- cbind(read.table(text = stouffer_liptak_aggregate$Group.1,
                                                  sep = ":"), stouffer_liptak_aggregate$x)
    colnames(stouffer_liptak_aggregate) <- c("CHR", "BP", "P")
    # make sure there are no zero p-values
    stouffer_liptak_aggregate$P[stouffer_liptak_aggregate$P == 0] <- .Machine$double.xmin
    
    # subset by chr is option is not NA
    if (!is.na(plot.chr)) {
      stouffer_liptak_aggregate <- stouffer_liptak_aggregate[stouffer_liptak_aggregate$CHR == plot.chr,]
    }
    
    suppressWarnings(qqman::manhattan(stouffer_liptak_aggregate, suggestiveline = FALSE, genomewideline = FALSE))
    
    man.df <- stouffer_liptak_aggregate
  }
  
  
  if (method == 'count') {
    # make aggregate count for regions
    regions <- c(paste0(results(hicexp)$chr, ':', results(hicexp)$region1),
                 paste0(results(hicexp)$chr, ':', results(hicexp)$region2))
    p.values <- c(results(hicexp)$p.adj, results(hicexp)$p.adj)
    count <- ifelse(p.values < alpha, 1, 0)
    
    ## count method
    count_aggregate <- aggregate(count, by = list(regions), 
                                 FUN = function(cnt) {
      c.sum <- sum(cnt)
      # c.pval <- 1 / sqrt(c.sum)
      return(c.sum)
    })
    
    # count_aggregate$x[is.infinite(count_aggregate$x)] <- 1 # replace any regions that had counts of 0 significant with pseudo-pval of 1
    count_aggregate <- cbind(read.table(text = count_aggregate$Group.1, sep = ":"),
                             count_aggregate$x)
    colnames(count_aggregate) <- c("CHR", "BP", "P")
    
    # subset by chr is option is not NA
    if (!is.na(plot.chr)) {
      count_aggregate <- count_aggregate[count_aggregate$CHR == plot.chr,]
    }
    # plot combined p-value manahttan plots
    suppressWarnings(.count_manhattan(count_aggregate, ylab = "Number of times significant"))
    
    man.df <- count_aggregate
    colnames(man.df)[3] <- 'count'
  }
  
  
  # return man.df if requested
  if (return_df) {
    return(man.df)
  } 
}



#' Plot the p-value results from topDirs
#' 
#' @export
#' @param dirs The output of the topDirs function when the return_df 
#'     option is set to "bed". 
#' @param plot.chr A numeric value indicating a specific
#'     chromosome number to subset the plot to. Defaults
#'     to NA indicating that all chromosomes will be plotted.
#' @return A plot.
#' @examples 
#' data('hicexp_diff')
#' dirs <- topDirs(hicexp_diff, return_df = 'bed')
#' plot_pvals(dirs)

plot_pvals <- function(dirs, plot.chr = NA) {
  # check input
  if (ncol(dirs) != 8) {
    stop("Use only the output of topDirs() where return_df = 'bed' for this function")
  }
  if (!is.na(plot.chr)) {
    if (!is.numeric(plot.chr)) {
      stop("plot.chr must be either NA or a numeric value.")
    }
    if (!(paste0('chr', plot.chr) %in% unique(dirs$chr))) {
      stop("The chr chosen as a subset must exist in the data object")
    }
  }
  # make dataframe for plotting
  df <- data.frame(CHR = as.numeric(sub('chr', '', dirs$chr)), BP = dirs$start, P = dirs$avgP.adj)
  # subset if plot.chr is not NA
  if(!is.na(plot.chr[1])) {
    df <- df[df$CHR == plot.chr,]
  }
  # plot combined p-value manahttan plots
  suppressWarnings(qqman::manhattan(df, suggestiveline = FALSE, genomewideline = FALSE))
}


#' Plot the count results from topDirs
#' 
#' @export
#' @param dirs The output of the topDirs function when the return_df 
#'     option is set to "bed". 
#' @param plot.chr A numeric value indicating a specific
#'     chromosome number to subset the plot to. Defaults
#'     to NA indicating that all chromosomes will be plotted.
#' @return A plot.
#' @examples 
#' data('hicexp_diff')
#' dirs <- topDirs(hicexp_diff, return_df = 'bed')
#' plot_counts(dirs)

plot_counts <- function(dirs, plot.chr = NA) {
  # check input
  if (ncol(dirs) != 8) {
    stop("Use only the output of topDirs() where return_df = 'bed' for this function")
  }
  if (!is.na(plot.chr)) {
    if (!is.numeric(plot.chr)) {
      stop("plot.chr must be either NA or a numeric value.")
    }
    if (!(paste0('chr', plot.chr) %in% unique(dirs$chr))) {
      stop("The chr chosen as a subset must exist in the data object")
    }
  }
  # make dataframe for plotting
  df <- data.frame(CHR = as.numeric(sub('chr', '', dirs$chr)), BP = dirs$start, P = dirs$count)
  # subset if plot.chr is not NA
  if(!is.na(plot.chr[1])) {
    df <- df[df$CHR == plot.chr,]
  }
  # plot
  suppressWarnings(.count_manhattan(df, ylab = "Number of times significant"))
}




# modified manhattan plot function for counts
.count_manhattan <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
                                                                                       "gray60"), chrlabs = NULL, suggestiveline = FALSE, 
                              genomewideline = FALSE, highlight = NULL, 
                              annotatePval = NULL, annotateTop = TRUE, ...) 
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) 
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  d$logp <- d$P
  
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index == 
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))+1), xlab = xlabel)
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos, 
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  # if (suggestiveline) 
  #   abline(h = suggestiveline, col = "blue")
  # if (genomewideline) 
  #   abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = "green3", pch = 20, 
                             ...))
  }
  if (!is.null(annotatePval)) {
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), textxy(pos, -log10(P), 
                                                offset = 0.625, labs = topHits$SNP, cex = 0.45), 
           ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
             labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  par(xpd = FALSE)
}
