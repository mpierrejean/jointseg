segmentByCBS <- structure(function(#Run CBS segmentation
### This function is a wrapper for convenient use of the \code{CBS}
### segmentation method by \code{\link{PSSeg}}.  It applies the
### \code{\link[DNAcopy]{segment}} function and reshapes the results
                                    y,
### A numeric vector, the signal to be segmented
                                   ...,
### Arguments to be passed to \code{\link[DNAcopy]{segment}}
                                    verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
) {
  ##seealso<<\code{\link[DNAcopy]{segment}}
  
  if (!require("DNAcopy")) {
    cat("Please install the 'DNAcopy' package to run the 'segmentByCBS' function")
    return()
  }
  n <- length(y)
  chrom <- rep(1, n)
  maploc <- 1:n
  genomdat <- y
  cna <- DNAcopy::CNA(genomdat, chrom, maploc)
  res <- DNAcopy::segment(cna)
  
  bkp <- res$output$loc.start[-1]

  res <- list(bkp=bkp) 
  return(res)
###  \item{bkp}{breakpoint positions}
}, ex=function(){
  ## load known real copy number regions
  affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
  
  ## generate a synthetic CN profile
  K <- 10
  len <- 1e5
  sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
  datS <- sim$profile

  ## run CBS segmentation
  res <- jointSeg:::segmentByCBS(datS[["c"]])
  getTprTnr(res$bkp, sim$bkp, nrow(datS), 5,relax = -1)
  plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
})

############################################################################
## HISTORY:
## 2013-01-09
## o Replace all jumps by bkp
## 2013-01-04
## o Created from 'segmentByCghseg'.
############################################################################

