segmentByRBS <- structure(function(#Recursive Binary Segmentation (low-level)
### Low-level function for multivariate Recursive Binary Segmentation (RBS)
                                   Y,
### A \code{n*p} signal to be segmented
                                   K,
### The number of change points to find
                                   minRegionSize=2,
### Regions with less than \code{minRegionSize} are not split
                                   verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                                   ) {
  ##details<<This function recrusively looks for the best candidate
  ##change point according to binary segmentation. This is the
  ##low-level function. It is generally advised to use the wrapper
  ##\code{\link{doRBS}} which also works on data frames and has a
  ##convenient argument \code{stat}.

  ##details<<See \code{\link{jointSeg}} for combining recursive binary
  ##segmentation with pruning by dynamic programming
  ##(\code{\link{pruneByDP}}).

  ##details<<See \code{\link{PSSeg}} for segmenting genomic signals
  ##from SNP arrays.

  ##seealso<<\code{\link{PSSeg}}, \code{\link{jointSeg}}, \code{\link{doRBS}}, \code{\link{pruneByDP}}
  
  ##references<<Gey, S., & Lebarbier, E. (2008). Using CART to Detect
  ##Multiple Change Points in the Mean for Large
  ##Sample. http://hal.archives-ouvertes.fr/hal-00327146/

  if (!is.matrix(Y) || !is.numeric(Y)){
    stop("Argument 'Y' should be a numeric matrix.\nPlease see 'doRBS' for using RBS directly on a data.frame or a numeric vector")
  }

  rownames(Y) <- NULL
  n <- as.numeric(nrow(Y))
  p <- dim(Y)[2]

  ##details<< Each dimension of the original signal is scaled before
  ##segmentation, using \code{\link{estimateSd}}.
  Y <- sweep(Y, MARGIN=2, STATS=estimateSd(Y), FUN="/")

  if (missing(K)) {
    stop("Please provide argument 'K'")
  }
  if (K>=n) {
    stop("Too many breakpoints are required")
  }

  getRSE <- function(interval) {
    idxs <- seq(from=interval[1], to=interval[2], by=1)
    ZZ <- Y[idxs,, drop=FALSE]
    ZZ <- sweep(ZZ, 2, colMeans(ZZ, na.rm=TRUE))
    rse <- sum(colSums(ZZ^2, na.rm=TRUE))
  }

  getSplit <- function(interval, allowNA=TRUE) {
    left <- interval[1]
    right <- interval[2]
    idxs <- seq(from=left, to=right, by=1)
    ZZ <- Y[idxs,, drop=FALSE]
    ZZ <- sweep(ZZ, 2, colMeans(ZZ, na.rm=TRUE))
    cand <- anotherBkp(ZZ)+left-1
    int1 <- c(left, cand)
    int2 <- c(cand+1, right)
    rseL <- getRSE(int1)
    rseR <- getRSE(int2)
    if (cand %in% interval) {
      ## Nothing to do.
    }
    c(cand=cand, left=left, right=right, rseL=rseL, rseR=rseR)
  }

  activeSet <- NULL

  for (kk in seq(length=K)) {
    if (kk==1) { ## First split
      interval <- c(1, n)
      res <- getSplit(interval)
      rse <- getRSE(interval)
      res[["gain"]] <- rse-(res[["rseL"]]+res[["rseR"]])
  
      split <- res  ## the only candidate split at this stage
      
      candSplits <- NULL
      activeSet <- rbind(activeSet, res)
    } else {  ## Next splits
      if (FALSE && verbose) {
        print(kk)
        print(split)
      }
      pos <- split[["cand"]]
      ## left interval
      leftInt <- c(split[["left"]], pos)
      if (diff(leftInt)>minRegionSize) {
        ## interval length is not too small: split !
        res <- getSplit(leftInt)
        res[["gain"]] <- split[["rseL"]]-(res[["rseL"]]+res[["rseR"]])
        candSplits <- rbind(candSplits, res)  ## add splits to candidates
      }
      
      ## right interval
      rightInt <- c(pos+1, split[["right"]])
      if (diff(rightInt)>minRegionSize) {
        ## interval length is not too small: split !
        res <- getSplit(rightInt)
        res[["gain"]] <- split[["rseR"]]-(res[["rseL"]]+res[["rseR"]])
        candSplits <- rbind(candSplits, res)  ## add splits to candidates
      }
      nCand <- length(candSplits)
      if (nCand==0) {
        warning("No more candidate splits for the desired 'minRegionSize': returning with only ", nrow(activeSet), " breakpoints")
        break;
      }
      
      wBest <- which.max(candSplits[, "gain"])
      split <- candSplits[wBest, ]
      candSplits <- candSplits[-wBest,, drop=FALSE]
      activeSet <- rbind(activeSet, split)
    } ## if (kk==1)
  } ## for (kk ...
  rownames(activeSet) <- NULL
  ##value<<A list with elements:
  list(
      bkp=activeSet[, "cand"],##<< A \code{vector} of \code{K} estimated breakpoint positions, sorted by order of appearance
      rse = rse-cumsum(c(0,activeSet[, "gain"])), ##<< the residual squared error (RSE) for the successive segmentations
      gain=activeSet[, "gain"]) ##<< The gain provided by each breakpoints in terms of difference between RSE
}, ex=function(){
  p <- 2
  trueK <- 10
  len <- 1e4
  sim <- randomProfile(len, trueK, 1, p)
  Y <- sim$profile
  K <- 2*trueK
  res <- segmentByRBS(Y, K)
  getTpFp(res$bkp, sim$bkp, tol=10, relax = -1)   ## true and false positives
  
  cols <- rep(2, K)
  cols[1:trueK] <- 3
  par(mfrow=c(p,1))
  for (ii in 1:p) {
    plot(Y[, ii], pch=19, cex=0.2)
    abline(v=res$bkp[1:trueK], col= cols)
    abline(v=sim$bkp, col=8, lty=2)
  }

  ## NA:s in one dimension at a true breakpoint
  jj <- sim$bkp[1]
  Y[jj-seq(-10, 10), p] <- NA
  res2 <- segmentByRBS(Y, K)
  getTpFp(res2$bkp, sim$bkp, tol=10, relax = -1)   ## true and false positives
  
  ## NA:s in both dimensions at a true breakpoint
  Y[jj-seq(-10, 10), ] <- NA
  res3 <- segmentByRBS(Y, K)
  getTpFp(res3$bkp, sim$bkp, tol=10, relax = -1)   ## true and false positives
})

############################################################################
## HISTORY:
## 2014-05-15
## o Renamed back to 'segmentByRBS', so that 'doRBS' is a *wrapper*
## around the core segmentation function.
## 2013-12-09
## o Renamed to 'doRBS'
## 2013-12-05
## o Now dropping row names of 'Y'.
## 2013-03-07
## o Add return parameter RSE to compute model selection on 'jointSeg'
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

