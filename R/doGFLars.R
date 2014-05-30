doGFLars <- structure(function(#Group fused Lars segmentation
### High-level function for multivariate fused Lars (GFLars) segmentation
                               Y,
### A \code{n*p} signal to be segmented
                               K,
### The number of change points to find
                               stat=NULL,
### A vector containing the names or indices of the columns of \code{Y} to be segmented
                               ...,
### Further arguments to be passed to 'segmentByGFLars'
                               verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                               ){
  ##details<<This function is a wrapper around the lower-level
  ##segmentation function \code{\link{segmentByGFLars}}. It can be run
  ##on p-dimensional, piecewise-constant data in order to defined a
  ##set of candidate change points. It is recommended to prune this
  ##list of candidates using dynamic programming
  ##(\code{\link{pruneByDP}}), combined with a selection of the best
  ##number of change points. The \code{\link{jointSeg}} function
  ##provides a convenient wrapper for performing segmentation, pruning
  ##and model selection.
  
  ##details<<For the specific case of DNA copy number data segmentation, see the
  ##dedicated wrapper \code{\link{PSSeg}}.

  ##details<<The default weights \eqn{\sqrt{n/(i*(n-i))}} are calibrated as
  ##suggested by Bleakley and Vert (2011).  Using this calibration,
  ##the first breakpoint maximizes the likelihood ratio test (LRT)
  ##statistic.

  ##references<< Bleakley, K., & Vert, J. P. (2011). The group fused
  ##lasso for multiple change-point detection. arXiv preprint
  ##arXiv:1106.4199.
  
  ##references<< Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
  ##change-points shared by many signals using group LARS. Advances in
  ##Neural Information Processing Systems, 23, 2343-2351.
  
  ##note<<This implementation is derived from the MATLAB code
  ##by Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.
  
  ## Argument 'Y'
  if (is.null(dim(Y)) || is.data.frame(Y)) {
    if (verbose) {
      print("Coercing 'Y' to a matrix")
    }
    Y <- as.matrix(Y)
  } else if (!is.matrix(Y)){
    stop("Argument 'Y' should be a matrix, vector or data.frame")
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
  

  ## Handle missing values by smoothing
  outPos <- 1:nrow(Y)
  if (any(is.na(Y))) {
    warning("Missing values detected. Smoothing will be performed.")
    Yb <- binMissingValues(Y)
    pos <- attr(Yb, "idxs")
    outPos <- mapPositionsBack(pos)
    Y <- Yb
  }

  n <- as.numeric(nrow(Y))
  p <- dim(Y)[2]

  ## Pass on to 'segmentByGFLars'
  res <- segmentByGFLars(Y, K, ..., verbose=verbose)

  ## map breakpoint positions back to original space (if required)
  res$bkp <- outPos[res$bkp]
  
  res
###An object of the same structure as the output of \code{\link{segmentByGFLars}}
}, ex=function(){
  p <- 2
  trueK <- 10
  sim <- randomProfile(1e4, trueK, 1, p)
  Y <- sim$profile
  K <- 2*trueK
  res <- doGFLars(Y, K)
  print(res$bkp)
  print(sim$bkp)
  plotSeg(Y, res$bkp)

  ## a toy example with missing values
  sim <- randomProfile(1e2, 1, 0.1, 2)
  Y <- sim$profile
  Y[3:50, 2] <- NA

  res <- doGFLars(Y, 10, 2, verbose=TRUE)
  print(res$bkp)
  print(sim$bkp)
  plotSeg(Y, res$bkp)
})

############################################################################
## HISTORY:
## 2014-05-15
## o Added argument 'stat'.
## o Now a wrapper calling the low-level 'segmentByGFLars'.
## 2013-12-09
## o Renamed to 'doGFLars'
## 2013-01-09
## o Replace all jumps by bkp
## 2012-12-27
## o Renamed to segmentByGFLars.
## o Some code and doc cleanups.
## 2012-12-13
## o Updated example.
## 2012-12-06
## o Replaced the function defaultWeigths by a (n-1)*1 vector of
## weigths and updated example
## 2012-09-13
## o Some code cleanups.
## o Tentative bug fix: indices in gammaTemp.
## o SPEEDUP: removed unnecessary calls to 'complex'.
## 2012-08-13
## o Created.
############################################################################

