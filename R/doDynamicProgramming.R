doDynamicProgramming <- structure(function(#Run segmentation by dynamic programming
### High-level function for univariate or multivariate segmentation by dynamic programming
                                           Y,
### A numeric vector or a matrix, the signal to be segmented
                                           K,
### The number of change points to find
                                           stat=NULL,
### A vector containing the names or indices of the columns of \code{Y} to be segmented
                                           verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                                           ) {
  ##note<< This is essentially a wrapper for convenient segmentation
  ##by dynamic programming using the \code{\link{PSSeg}} function.
  
  ##references<<Rigaill, G. (2010). Pruned dynamic programming for
  ##optimal multiple change-point detection. arXiv preprint
  ##arXiv:1004.0887.
  
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
  
  if (is.null(dim(Y)) || (ncol(Y)==1)) {
    ##details<< if the signal is uni-dimensional, this function simply
    ##uses the segmentation method provided in the
    ##\code{cghseg} package reshapes the results.
    if (!require("cghseg")) {
      cat("Please install the 'cghseg' package to run the 'doDynamicProgramming' function")
      return()
    }

    n <- length(Y)
    res <- cghseg:::segmeanCO(Y, Kmax=K+1)
    ## Note: segmeanCO is a low-level segmentation method in package 'cghseg'.
    ## It is not exported.
    bkpList <- lapply(1:K+1, FUN=function(kk) {
      res$t.est[kk, 1:(kk-1)]
    })
    dpseg <- list(bkp=bkpList, rse=res$J.est)
    res <- list(bkp=bkpList[[K]], dpseg=dpseg)
    ##stop("Argument 'y' should be a numeric vector")
  } else {
    ##details<< if the signal is multi-dimensional, this function
    ##applies the \code{\link{pruneByDP}} function and reshapes the
    ##results.
    res <- pruneByDP(Y, K=K+1)
    dpseg <- list(bkp=res$bkpList, rse=res$rse)
    res <- list(bkp=res$bkpList[[K]], dpseg=dpseg)
  }
  
  return(res)
### \item{bkp}{A vector of \code{K} indices for candidate change points}
### \item{dpseg}{A list of two elements \describe{
###   \item{bkp}{A list of vectors of change point positions for the best model with k change points, for k=1, 2, ... K}
###   \item{rse}{A vector of K+1 residual squared errors}
###  }}
}, ex=function(){
  ## load known real copy number regions
  affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=1)

  ## generate a synthetic CN profile
  K <- 10
  len <- 1e4
  sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
  datS <- sim$profile

  ## run cghseg segmentation
  resDP <- doDynamicProgramming(datS[["c"]], K=K)
  getTpFp(resDP$bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
  plotSeg(datS, breakpoints=list(sim$bkp, resDP$bkp))
  
  ## run 2d dynamic programming segmentation
  K <- 2
  len <- 1e3
  sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
  datS <- sim$profile
  datS$d <- 2*abs(datS$b-1/2)
  datS[which(datS$genotype!=0.5),"d"] <- NA
  Y = cbind(datS$c,datS$d)
  resDP2d <- doDynamicProgramming(Y, K = K)
})

############################################################################
## HISTORY:
## 2014-05-30
## o Added argument 'stat'.
## o Updated doc.
## 2013-12-09
## o Renamed to 'doDynamicProgramming'
## o Added 2d dynamic programming
## 2013-01-09
## o Replace all jumps by bkp
## 2013-01-04
## o BUG FIX: index shift when reshaping the results.
## 2013-01-03
## o Created.
############################################################################

