doRBS <- structure(function(#Recursive Binary Segmentation 
### High-level function for multivariate Recursive Binary Segmentation (RBS)
                            Y,
### A \code{vector}, \code{matrix} or \code{data.frame}) containing signals to be segmented
                            K,
### The number of change points to find
                            stat=NULL,
### A vector containing the names or indices of the columns of \code{Y} to be segmented
                            ...,
### Further aguemnts to be passed to 'segmentByRBS'.
                            verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                            ) {
  ##details<<This function is a wrapper aroung the lower-level
  ##segmentation function \code{\link{segmentByRBS}}. It can be run on
  ##p-dimensional, piecewise-constant data in order to defined a set
  ##of candidate change points. It is recommended to prune this list
  ##of candidates using dynamic programming (\code{\link{pruneByDP}}),
  ##combined with a selection of the best number of change points. The
  ##\code{\link{jointSeg}} function provides a convenient wrapper for
  ##performing segmentation, pruning and model selection.
  
  ##details<<For the specific case of DNA copy number data segmentation, see the
  ##dedicated wrapper \code{\link{PSSeg}}.

  ##seealso<<\code{\link{segmentByRBS}}, \code{\link{PSSeg}},
  ##\code{\link{pruneByDP}}
  
  ##references<<Gey, S., & Lebarbier, E. (2008). Using CART to Detect
  ##Multiple Change Points in the Mean for Large
  ##Sample. http://hal.archives-ouvertes.fr/hal-00327146/

  ## Argument 'Y'
  if (is.null(dim(Y)) || is.data.frame(Y)) {
    if (verbose) {
      print("Coercing 'Y' to a matrix")
    }
    Y <- as.matrix(Y)
  } else if (!is.matrix(Y)){
    stop("Argument 'Y' should be a matrix, vector or data.frame")
  }

  ## Assert that 'Y' is numeric
  if (!is.numeric(Y)) {
    stop("The signal to be segmented must be numeric!")
  }

  ## Argument 'stat'
  if (!is.null(stat)) {
    if (is.numeric(stat)) {
      mm <- match(stat, 1:ncol(Y)) 
    } else if (is.character(stat)) {
      mm <- match(stat, colnames(Y))
    }
    if (sum(is.na(mm))) {
      guilty <- paste("'", stat[which(is.na(mm))], "'", sep="", collapse=",")
      stop("Undefined column(s) selected in 'Y':", guilty, ". Please check argument 'stat'")
    } else {
      Y <- Y[, mm, drop=FALSE]
    }
  }
  
  ## Pass on to 'segmentByRBS'
  res <- segmentByRBS(Y, K, ..., verbose=verbose)
###  \item{bkp}{A \code{vector} of \code{K} estimated breakpoint positions, ranked
###    by order of appearance}
###  \item{gain}{The gain provided by the breakpoints in terms of difference between RSE} 
}, ex=function(){
  p <- 2
  trueK <- 10
  len <- 1e4
  sim <- randomProfile(len, trueK, 1, p)
  Y <- sim$profile
  K <- 2*trueK
  res <- doRBS(Y, K)
  getTpFp(res$bkp, sim$bkp, tol=10)   ## true and false positives
  
  cols <- rep(2, K)
  cols[1:trueK] <- 3
  par(mfrow=c(p,1))
  for (ii in 1:p) {
    plot(Y[, ii], pch=19, cex=0.2)
    abline(v=res$bkp[1:trueK], col= cols)
    abline(v=sim$bkp, col=8, lty=2)
  }
})

############################################################################
## HISTORY:
## 2014-05-15
## o Now a wrapper calling the low-level 'segmentByRBS'.
## 2013-12-09
## o Renamed to 'doRBS'.
## 2013-12-05
## o Now dropping row names of 'Y'.
## 2013-03-07
## o Add return parameter RSE to compute model selection on 'jointSeg'.
## 2013-01-23
## o BUG FIX: Empty candidate list would give an error.  Now returning
## early when 'minRegionSize' is too large for 'K'.
## o Replace all jumps by bkp
## 2013-01-09
## o Replace all jumps by bkp
## 2012-12-31
## o Now using 'anotherBkp' instead of 'oneBkp' in order for missing
## values to be handled.
## 2012-12-30
## o Now properly dealing with the special case K=0.
## 2012-12-27
## o Renamed to 'segmentByRBS'.
## o Some code and doc cleanups.
## 2012-12-23
## o SPEEDUP: removed redundant calls to 'getRSE'.
## 2012-12-07
## o BUG FIX: index shift in correspondence b/w breakpoint position and interval.
## 2012-12-05
## o Created.
############################################################################

