#' Plot signal and breakpoints with segment-level signal estimates
#'
#' Plot signal and breakpoints with segment-level signal estimates
#'
#' Argument 'binCols' is mainly used to avoid calculating mean levels for
#' allelic ratios, which would not make sense as they are typically
#' multimodal.
#'
#' @param dat A \code{matrix} or data frame whose rows correspond to loci
#' sorted along the genome, or a \code{numeric} \code{vector}.
#' @param breakpoints A vector of breakpoints positions, or a \code{list} of
#' such vectors.
#' @param regNames Region labels, a vector of length
#' \code{length(breakpoints)+1} (if \code{breakpoints} is a vector) or of
#' length \code{length(breakpoints[[1]])+1} (if \code{breakpoints} is a list).
#' @param exclNames A vector of column names corresponding to columns that
#' should not be plotted.
#' @param ylabs A vector of 'y' labels (column names or indices) that should be
#' plotted.
#' @param ylims An optional \eqn{2*d} matrix with \code{ylim} values for each
#' of the \eqn{d} dimensions to be plotted.
#' @param binExclPattern A vector of column names or indices in
#' \code{colnames(dat)} for which segment-level signal estimates should *not*
#' be drawn.
#' @param col Color of plotting symbol, see \code{\link{par}}
#' @param pch Plotting symbol, see \code{\link{par}}
#' @param cex Magnification factor for plotting symbol, see \code{\link{par}}
#' @author Morgane Pierre-Jean and Pierre Neuvial
#' @examples
#'
#' affyDat <- acnr::loadCnRegionData(dataSet="GSE29172", tumorFraction=1)
#' sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)
#' dat <- sim$profile
#' res <- PSSeg(dat, method="RBS", stat=c("c", "d"), K=50)
#' bkpList <- list(true=sim$bkp, est=res$bestSeg)
#' plotSeg(dat, breakpoints=bkpList)

#' @export plotSeg
plotSeg <- function(dat, breakpoints=NULL, regNames=NULL,
                    exclNames=c("genotype", "region", "bT", "bN", "cellularity"),
                    ylabs=colnames(dat), ylims=NULL, binExclPattern="^b[N|T]*$",
                    col="#33333333", pch=19, cex=0.3){
    if (is.null(dim(dat))) {
        ## coerce to a matrix
        dat <- as.matrix(dat)
    }
    n <- nrow(dat)
    pos <- 1:n

    ## Argument 'regNames'
    if (is.null(regNames)) {
        mm <- match("region", colnames(dat))
        if (!is.na(mm)) {
            regNames <- dat[, mm]
        }
    }
    regLabs <- NULL

    ## Argument 'exclNames'
    idxsE <- stats::na.omit(match(exclNames, colnames(dat)))
    if (length(idxsE)) {
        dat <- dat[, -idxsE, drop=FALSE]
    }

    p <- ncol(dat)
    ## Argument 'ylabs'
    if (is.null(ylabs)) {
        ylabs <- paste("y", 1:p, sep="")
    } else if (length(ylabs) != p) {
        stop("Argument 'ylabs' does not match signal dimension")
    }

    ## Argument 'ylim'
    if (!is.null(ylims)) {
        stopifnot(nrow(ylims)==2)
    }

    ## Argument 'binExclPattern'
    binCols <- grep(binExclPattern, ylabs, invert=TRUE)  ## those to include !
    if(!is.null(breakpoints)){
        if (!is.list(breakpoints)) {  ## coerce to a list
            breakpoints <- list(breakpoints)
        }
        breakpoints <- lapply(breakpoints, sort)
        if (!is.null(regNames)) {
            regLabs <- regNames[c(breakpoints[[1]], n)]
        }

        meanList <- lapply(breakpoints, FUN=function(bkp) {
            xOut <- c(min(pos), bkp, max(pos))
            xOut <- sort(unique(xOut))
            binDat <- matrix(NA, nrow=length(xOut)-1, ncol=length(binCols))
            colnames(binDat) <- colnames(dat)[binCols]
            for (cc in seq(along=binCols)) {
                bc <- binCols[cc]
                means <- matrixStats::binMeans(y=dat[, bc], x=pos, bx=xOut)
                binDat[, cc] <- means
            }
            binDat
        })
    }

    xlim <- range(pos)
    nl <- is.null(ylims)
    graphics::par(mfrow = c(p, 1), mar=c(3.2, 4, 1, 0)+0.5)
    for (cc in 1:p) {
        y <- dat[, cc]
        ## xlab <- ifelse(cc==p, "position", "")
        xlab <- ""
        if (nl) {
            ylim <- stats::quantile(y, c(0.001, 0.999), na.rm=TRUE)
        } else {
            ylim <- ylims[, cc]
        }
        graphics::plot(NA , ylim=ylim, xlim=xlim, xlab=xlab, ylab=ylabs[cc])
        graphics::points(pos, y, cex=cex, col=col, pch=pch)
        if(!is.null(breakpoints)){
            for(ll in seq(along=breakpoints)){
                bkp <- breakpoints[[ll]]
                bkpStart <- c(1, bkp+1)
                bkpEnd <- c(bkp, max(pos))
                mm <- match(cc, binCols)
                if (!is.na(mm)) {
                    val <- meanList[[ll]][, mm]
                    graphics::segments(bkpStart, val, bkpEnd, val, col=ll+1, lwd=2, lty=ll)
                }
                graphics::abline(v=bkp, col=ll+1, lwd=2, lty=ll)
                if (ll==1 & !is.null(regNames)) {  ## add region labels
                    graphics::mtext(regLabs, side=3, line=0, at=(bkpStart+bkpEnd)/2)
                }
            }
        }
    }
}

############################################################################
## HISTORY:
## 2014-02-14:
## o BUG FIX: Calculation of segment means would fail when a breakpoint
##   was detected between the first and second position.
## 2014-01-22
## o BUG FIX: only means of 'c' were plotted.
## 2013-05-30
## o Added argument 'regNames' so that region labels are plotted if
##  available.
##  'plotSeg' can handle not only copy number signals.
## 2013-01-23
## o Added arguments 'exclNames', 'ylabs', and 'binExclPattern' so that
##  'plotSeg' can handle not only copy number signals.
## 2013-01-09
## o Replace all jumps by bkp
## 2012-12-27
## o Some code and doc cleanups.
## o Now using matrixStats::binMeans.
## 2012-11-30
## o Created.
############################################################################

