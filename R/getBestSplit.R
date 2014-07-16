getBestSplit <- structure(function(#Get best candidate change point
### Get best candidate change point according to binary segmentation
                                   left,
### \code{numeric} or \code{integer} value specifying the begining
### position of the interval to be split
                                   right,
### \code{numeric} or \code{integer} value specifying the end position
### of the interval to be split
                                   Y,
### signal
                                   weights=NULL,
### segmentation weights.  If \code{NULL}, weights are set to 1 for all non missing entries.
                                   checkLRT=FALSE,
### Should the consistency between the likelihood ratio test (LRT)
### statistic and the residual squared error (RSE) be tested ?
                                   returnAllGains=FALSE,
### Should the "gains" (that is, increases in RSE) be reported for all candidate change points in \code{interval} or only for the best one? By default, only the best one is reported. If \code{TRUE}, the gains are returned via the attribute \code{gains} of the return value.
                                   verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                                 ){
  ##keyword<<internal

  ## sanity checks
  p <- ncol(Y)
  n <- nrow(Y)
  stopifnot(left>0)
  stopifnot(left<right)
  stopifnot(right<=n)
  Y <- Y[left:right,, drop=FALSE ]
  n <- nrow(Y)

  ## Argument 'weights'
  if (is.null(weights)) {
    weights <- 0+(!is.na(Y))
  }
  stopifnot(identical(dim(weights), dim(Y)))

  ## handle missing values (part 1)
  dummy <- 0;
  Y[which(weights==0)] <- dummy;
  rm(dummy);
    
  rseL <- 0
  rseR <- 0
  rseT <- 0
  LRT2 <- 0

  meanDD <- NA
  
  for (dd in 1:p) {
    yy <- Y[, dd]
    ww <- weights[, dd]

    S0 <- cumsum(ww)
    S1 <- cumsum(yy*ww)
    S2 <- cumsum(yy*yy*ww)
    
    S0n <- S0[n]
    if (S0n==0) {  ## No non-zero weight value in this dimension (possibly because all entries are missing): break
      break
    }
    S1n <- S1[n]
    S2n <- S2[n]

    S0 <- S0[-n]
    S1 <- S1[-n]
    S2 <- S2[-n]
  
    ## RSE
    wL <- 1/S0
    wL[wL==Inf] <- 0       ## handle missing values at the beginning of interval
    rseL <- rseL +  S2      - wL*S1^2
    
    wR <- 1/(S0n-S0)
    wR[wR==Inf] <- 0       ## handle missing values at the end of interval
    rseR <- rseR + (S2n-S2) - wR*(S1n-S1)^2

    wT <- 1/S0n
    rseT <- rseT +  S2n     - wT*S1n^2

    ## LRT
    if (checkLRT) {
      score <- (wL*S1-wR*(S1n-S1))^2
      sf <- wT*S0*(S0n-S0)
      LRT2 <- LRT2 + score*sf

      stopifnot(max(abs(rseL+rseR-rseT+LRT2))<1e-10)
    }
  }

  idx <- which.min(rseL+rseR)
  cand <- idx + left-1
  
  rseLc <- rseL[idx]
  rseRc <- rseR[idx]
  gain <- rseT-(rseLc+rseRc)
  res <- c(left=left, right=right, cand=cand, rseL=rseLc, rseR=rseRc, rseT=rseT, gain=gain)
  if (returnAllGains) {
    attr(res, "gains") <-  rseT-(rseL+rseR)
  }
  res
}, ex=function() {
  n <- 1e4
  sim <- randomProfile(length=n, nBkp=1L, noiseLevel=1, dim=2L, minLength=3L)
  Y <- sim$profile
  
  YY <- sweep(Y, MARGIN=2, STATS=estimateSd(Y), FUN="/")
  YY[is.na(YY)] <- 0
  
  split <- getBestSplit(1L, n, Y)
  print(split[["cand"]]-sim$bkp)  ## generally not to far from 0

  n <- 30
  sim <- randomProfile(length=n, nBkp=1L, noiseLevel=0.5, dim=2L, minLength=3L)
  Y <- sim$profile
  split <- getBestSplit(1L, n, Y)
  print(split[["cand"]]-sim$bkp)  ## generally not to far from 0

  Y <- sim$profile
  Y[6:10, 2] <- NA
  split <- getBestSplit(1L, n, Y)
  print(split[["cand"]]-sim$bkp)  ## generally not to far from 0
  
  Y <- sim$profile
  Y[1:10, 2] <- NA
  split <- getBestSplit(1L, n, Y)
  print(split[["cand"]]-sim$bkp)  ## generally not to far from 0

  Y <- sim$profile
  Y[5:30, 2] <- NA
  split <- getBestSplit(1L, n, Y)
  print(split[["cand"]]-sim$bkp)  ## generally not to far from 0
})

############################################################################
## HISTORY:
## 2014-07-02
## o Created from 'getSplit' former local function in 'segmentByRBS'.
############################################################################
