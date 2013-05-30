plotSeg <- structure(function(# Plot signal and breakpoints with segment-level signal estimates
### Plot signal and breakpoints with segment-level signal estimates
                              dat,
### A \code{matrix} or data frame whose rows correspond to loci sorted along the genome.
                              breakpoints=NULL,
### A vector of breakpoints positions, or a \code{list} of such vectors.
                              regNames=NULL,
### Region labels, a vector of length \code{length(breakpoints)+1} (if
### \code{breakpoints} is a vector) or of length
### \code{length(breakpoints[[1]])+1} (if \code{breakpoints} is a
### list).                              
                              exclNames=c("genotype", "region"),
### A vector of column names corresponding to columns that should not
### be plotted.
                              ylabs=colnames(dat),
### A vector of 'y' labels (among columns that should be plotted).
                              ylims=NULL,
                              ### An optional \eqn{2*d} matrix with \code{ylim} values for each of the \eqn{d} dimensions to be plotted.
                              binExclPattern="^b$"
### A vector of column indices in \code{colnames(dat)} for which
### segment-level signal estimates should be drawn.
                  ){
  ##details<<Argument 'binCols' is mainly used to avoid
  ##calculating mean levels for allelic ratios, which would not make
  ##sense are they are typically multimodal.
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
  idxsE <- na.omit(match(exclNames, colnames(dat)))
  if (length(idxsE)) {
    dat <- dat[, -idxsE]
  }
  ## if (length(idxsE)) {
  ##   mm <- match(idxsE, 1:ncol(dat))
  ##   if (any(is.na(mm))) {
  ##     ww <- paste(which(is.na(mm)), collapse=",")
  ##     stop("Could not exclude column ", ww, ". Please check argument 'exclude'")
  ##   }
  ##   dat <- dat[, -mm]
  ## }
  
  p <- ncol(dat)
  ## Argument 'ylabs'
  if (is.null(ylabs)) {
    ylabs <- rep("", p)
  } else if (length(ylabs) != p) {
    stop("Argument 'ylabs' does not match signal dimension")
  }
  
  ## Argument 'binExclPattern'
  binCols <- grep(binExclPattern, colnames(dat), invert=TRUE)  ## those to include !
  ## mm <- match(binCols, 1:ncol(dat))
  ## if (any(is.na(mm))) {
  ##   ww <- paste(which(is.na(mm)), collapse=",")
  ##   stop("Column index not found: ", ww, ". Please check argument 'binColPatterns'")
  ## }

  if(!is.null(breakpoints)){
    if (!is.list(breakpoints)) {  ## coerce to a list
      breakpoints <- list(breakpoints)
    }
    breakpoints <- lapply(breakpoints, sort)
    if (!is.null(regNames)) {
      regLabs <- regNames[c(breakpoints[[1]], n)]
    }
    
    meanList <- lapply(breakpoints, FUN=function(bkp) {
      binDat <- NULL
      for (cc in binCols) {
        xOut <- c(min(pos), bkp, max(pos))
        xOut <- sort(unique(xOut))
        means <- matrixStats::binMeans(y=dat[, "c"], x=pos, bx=xOut)
        binDat <- cbind(binDat, means)
      }
      colnames(binDat) <- colnames(dat)[binCols]
      binDat
    })
  }

  xlim <- range(pos)
  nl <- is.null(ylims)
  par(mfrow = c(p, 1), mar=c(3.2, 4, 1, 0)+0.5)
  for (cc in 1:p) {
    y <- dat[, cc]
    ## xlab <- ifelse(cc==p, "position", "")
    xlab <- ""
    if (nl) {
      ylim <- quantile(y, c(0.001, 0.999), na.rm=TRUE)
    } else {
      ylim <- ylims[, cc]
    }
    plot(NA , ylim=ylim, xlim=xlim, xlab=xlab, ylab=ylabs[cc])
    points(pos, y, cex=0.3)
    if(!is.null(breakpoints)){  
      for(ll in seq(along=breakpoints)){
        bkp <- breakpoints[[ll]]
        bkpStart <- c(1, bkp+1)
        bkpEnd <- c(bkp, max(pos))
        mm <- match(cc, binCols)
        if (!is.na(mm)) {
          val <- meanList[[ll]][, mm]
          segments(bkpStart, val, bkpEnd, val, col=ll+1, lwd=2, lty=ll)
        }
        abline(v=bkp, col=ll+1, lwd=2, lty=ll)
        if (ll==1 & !is.null(regNames)) {  ## add region labels
          mtext(regLabs, side=3, line=0, at=(bkpStart+bkpEnd)/2)
        }
      }
    }
  }
},ex = function(){	
  affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
  sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)
  dat <- sim$profile
  res <- PSSeg(dat, K=50)
  bkpList <- list(true=sim$bkp, est=res$bestSeg)
  plotSeg(dat, breakpoints=bkpList)
})

############################################################################
## HISTORY:
## 2013-05-30
## o Added argument 'regNames' so that region labels are plotted if
##  available.
## 'plotSeg' can handle not only copy number signals.
## 2013-01-23
## o Added arguments 'exclNames', 'ylabs', and 'binExclPattern' so that 
## 'plotSeg' can handle not only copy number signals.
## 2013-01-09
## o Replace all jumps by bkp
## 2012-12-27
## o Some code and doc cleanups.
## o Now using matrixStats::binMeans.
## 2012-11-30
## o Created.
############################################################################

