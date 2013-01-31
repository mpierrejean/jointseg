segmentByPairedPSCBS <- structure(function(#Run Paired PSCBS segmentation
### This function is a wrapper for convenient use of the \code{PSCBS}
### segmentation method by \code{\link{PSSeg}}.  It applies the
### \code{\link[PSCBS]{segmentByPairedPSCBS}} function and reshapes the results
                                    y,
### A numeric vector, the signal to be segmented
                                   ...,
### Arguments to be passed to \code{\link[PSCBS]{segmentByPairedPSCBS}}
                                    verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
) {
  ##seealso<<\code{\link[PSCBS]{segmentByPairedPSCBS}}
  
  if (!require("PSCBS")) {
    cat("Please install the 'DNAcopy' package to run the 'segmentByCBS' function")
    return()
  }
  n <- length(y)
  chrom <- rep(1, n)
  ## TODO: update for PSCBS  !!
  maploc <- 1:n
  genomdat <- y
  cna <- CNA(genomdat, chrom, maploc)
  res <- segment(cna)
  
  bkp <- res$output$loc.start[-1]
  ## END TODO
  res <- list(bkp=bkp) 
  return(res)
###  \item{bkp}{breakppoint indices}
}, ex=function(){
  ## load known real copy number regions
  affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
  sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)

  ## generate a synthetic CN profile
  K <- 10
  len <- 1e5
  sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
  datS <- sim$profile

  ## run Paired PSCBS segmentation
##  res <- segmentByPairedPSCBS(datS[["c"]])
##  getTprTnr(res$bkp, sim$bkp, nrow(datS), 5)
##  plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
})

############################################################################
## HISTORY:
## 2013-01-09
## o Replace all jumps by bkp
## 2013-01-04
## o Created from 'segmentByCBS'.
############################################################################

