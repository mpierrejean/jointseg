anotherBkp <- structure(function(#Get best candidate change point
### Get best candidate change point according to binary segmentation
                                 Y,
### A \code{n*p} matrix, \code{p} signals of length \code{n} to be
### segmented (centered by column)
                                 weightFUN=defaultWeights,
### A \code{function} returning a \code{(n-1)*1} vector of weights for
### the candidate change point positions. Default weights yield the
### likelihood ratio test (LRT) statistic for the identification of a
### single change point.
                                 verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                                 ){
  ##keyword<<internal

  ##details<<Contrary to \code{oneBkp}, \code{anotherBkp} handles missing values (NA:s).
  if (!is.matrix(Y)) {
    stop("Y is not a matrix, please check dimension of Y")
  }
  C <- apply(Y, 2, getUnivStat, weightFUN=weightFUN)
  cNorm <- rowSums(C^2, na.rm=TRUE)
  which.max(cNorm)
}, ex=function(){
  p <- 2
  n <- 100
  
  sim <- randomProfile(n, 1, 1, p)
  Y <- sim$profile
  bkp <- anotherBkp(Y)
  print(bkp)
  print(oneBkp(Y))
  ##  stopifnot(identical(oneBkp(Y), bkp))
  plotSeg(Y, list(sim$bkp, bkp)
  
  ## robustness to NA:s
  h <- 2
  idxs <- seq(from=max(sim$bkp[1]-h, 1), min(sim$bkp[1]+h, n))
  Y[idxs, p] <- NA
  oneBkp(Y)  ## does not work
  bkp <- anotherBkp(Y)  ## works
  bkp-sim$bkp
})

############################################################################
## HISTORY:
## 2013-01-30
## o Now an internal function (not exported anymore).
## o Added 'jointSeg:::' to example because function is not exported anymore.
## 2013-01-09
## o Replaced 'jump' by 'bkp'.
## 2012-12-30
## o Created from 'oneBkp'.
############################################################################

