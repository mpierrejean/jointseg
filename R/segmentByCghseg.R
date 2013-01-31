segmentByCghseg <- structure(function(#Run cghseg segmentation
### This function is a wrapper for convenient use of the \code{cghseg}
### segmentation method by \code{\link{PSSeg}}.  It applies the
### \code{segmeanCO} function from package \code{cghseg} and reshapes
### the results.
                                    y,
### A numeric vector, the signal to be segmented
                                   K,
### The number of change points to find
                                    verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
) {
  ##references<<Rigaill, G. (2010). Pruned dynamic programming for
  ##optimal multiple change-point detection. arXiv preprint
  ##arXiv:1004.0887.

  if (!require("cghseg")) {
    cat("Please install the 'cghseg' package to run the 'segmentByCghseg' function")
    return()
  }
  n <- length(y)
  res <- cghseg:::segmeanCO(y, K=K+1)
  bkpList <- lapply(1:K+1, FUN=function(kk) {
    res$t.est[kk, 1:(kk-1)]
  })
  dpseg <- list(bkp=bkpList, rse=res$J.est)
  res <- list(bkp=bkpList[[K]], dpseg=dpseg)
  return(res)
### \item{bkp}{A vector of \code{K} indices for candidate change points}
### \item{dpseg}{A list of two elements \describe{
###   \item{bkp}{A list of vectors of change point positions for the best model with k change points, for k=1, 2, ... K}
###   \item{rse}{A vector of K+1 residual squared errors}
###  }}
}, ex=function(){
  ## load known real copy number regions
  affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
  sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)

  ## generate a synthetic CN profile
  K <- 10
  len <- 1e5
  sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
  datS <- sim$profile

  ## run cghseg segmentation
  resDP <- segmentByCghseg(datS[["c"]], K=K)
  getTprTnr(resDP$bkp, sim$bkp, nrow(datS), 5)
  plotSeg(datS, breakpoints=list(sim$bkp, resDP$bkp))
})

############################################################################
## HISTORY:
## 2013-01-09
## o Replace all jumps by bkp
## 2013-01-04
## o BUG FIX: index shift when reshaping the results.
## 2013-01-03
## o Created.
############################################################################

