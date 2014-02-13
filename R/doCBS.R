doCBS <- structure(function(#Run CBS segmentation
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
    cat("Please install the 'DNAcopy' package to run the 'doCBS' function")
    return()
  }
  if (!is.null(dim(y)) || mode(y)!="numeric") {
    stop("Argument 'y' should be a numeric vector")
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
  if (require("DNAcopy")) {
    ## load known real copy number regions
    affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
    
    ## generate a synthetic CN profile
    K <- 10
    len <- 1e5
    sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
    datS <- sim$profile
    
    ## run CBS segmentation
    res <- doCBS(datS[["c"]])
    getTpFp(res$bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
    plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
  }
})

############################################################################
## HISTORY:
## 2013-12-09
## o Renamed to 'doCBS'
## 2013-05-16
## o Example code now embedded in a 'require()' statement to avoid
##   problems in the R CMD check mechanism of R-forge.
## 2013-01-09
## o Replaced 'jump' by 'bkp'.
## 2013-01-04
## o Created from 'segmentByCghseg'.
############################################################################

