doPSCN <- structure(function(#Run PSCN segmentation
### This function is a wrapper for convenient use of the \code{\link[PSCN]{PSCN}}
### segmentation method by \code{\link{PSSeg}}.  It applies the
### \code{\link[PSCN]{smoothing}} and \code{\link[PSCN]{segmentation}}
### functions.
                                    Y,
### The signal to be segmented, a matrix containing the following columns: \describe{
### \item{c}{Total copy number (log scale)}
### \item{b}{Allele B fraction (a.k.a. BAF)}
### }
                                    alpha=0.01,
### sensitivity level in [0,1] to be passed to
### \code{\link[PSCN]{segmentation}}.
                                    platform=c("Illumina", "Affymetrix"),
### Specifies form which array platform 'Y' was generated: Illumina or GSE29172
                                    verbose=FALSE,
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                                    ...
### Further arguments to be passed to \code{\link[PSCN]{smoothing}}
) {
  ##references<<Chen, H., Xing, H., & Zhang, N. R. (2011). Estimation
  ##of parent specific DNA copy number in tumors using high-density
  ##genotyping arrays. PLoS computational biology, 7(1), e1001060.

  ##seealso<<\code{\link[PSCN]{smoothing}}
  ##seealso<<\code{\link[PSCN]{segmentation}}
  ##seealso<<\code{\link{PSSeg}}
  
  if (version$major>=3) {
    stop("The 'PSCN' package is not available for R >= 3.0.0.\nSee http://cran.r-project.org/web/packages/PSCN/index.html")
  }    

  if (!require(PSCN)) {
    cat("Please install the 'PSCN' package to run the 'doPSCN' function")
    return()
  }
  n <- nrow(Y)
  cn <- colnames(Y)
  ecn <- c("c", "b") ## expected
  mm <- match(ecn, cn)
  if (any(is.na(mm))) {
    str <- sprintf("('%s')", paste(ecn, collapse="','"))
    stop("Argument 'Y' should contain columns named ", str)
  }
  
  platform <- match.arg(platform)
  if (verbose) print(platform)
  isLogScaled <- (sum(Y[, "c"]<0)>1)
  if (!isLogScaled) {
    stop("Input copy number data does not seem to be on the log-scale (no negative values)\nPlease make sure that this is the case")
  }
  logR <- Y[, "c"]
  bfreq <- Y[, "b"]
  ## ad hoc: make sure 'bfreq' is set to "0" for CN probes !!!  This
  ## is required by PSCN because NA:s are dropped by PSCN even for
  ## total CN calculation...
  wCN <- which(is.nan(bfreq))
  bfreq[wCN] <- 0
  genotype.freq = matrix(NA, n, 4)
  genotype.freq[wCN, ] <- c(1, 0, 0, 0)

  ##details<<The data are assumed to come from a single chromosome.
  chr <- 1;
  inputdata <- data.frame(Position=1:n, Chr=rep(chr, n), logR=logR, bfreq=bfreq)
  
  tmpPath <- tempdir()
  sampleName <- "dummy"
  pathname <- file.path(tmpPath, sampleName)
  ## "Smoothing" (ie the HMM)
  resSmt <- PSCN::smoothing(pathname, inputdata, platform=platform,...)

  ## "Segmentation" (ie some post-processing)
  resSeg <- PSCN::segmentation(pathname, chr=chr, verbose=verbose, combine.alpha = alpha)
  
  suffix <- paste(".Chr", chr, ".Segment.Rdata", sep = "")
  pathnameS <- paste(pathname, suffix, sep="")
  cnvobj <- NULL; ## To please R CMD check
  load(pathnameS)
  bkp <- cnvobj$chpts.new ## bkp indices
##  bkpPos <- cnvobj$pos[bkp] ## bkp position in the original data
  rm(cnvobj);

  res <- list(bkp=bkp) 
  return(res)
###  \item{bkp}{breakpoint positions}
}, ex=function(){
  if (require("PSCN") && FALSE) { ## PSCN has been removed from CRAN
                                  ## and does not seem to work anymore
                                  ## with R >= 3.0.0
    ## load known real copy number regions
    affyDat <- loadCnRegionData(platform="GSE29172", tumorFraction=1)
    sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)
    datS <- sim$profile
    
    ## run PSCN segmentation
    Y <- cbind(c=log2(datS[, "c"])-1, datS[, "b"])  ## Convert to log ('LRR') scale
    resPSCN <- doPSCN(Y)
    getTpFp(resPSCN$bkp, sim$bkp, tol=20, relax = -1)   ## true and false positives
    plotSeg(datS, breakpoints=list(resPSCN$bkp, sim$bkp))
  }
})

############################################################################
## HISTORY:
## 2014-02-13
## o Moved to (newly created) 'defunct' directory.
## 2013-12-09
## o Renamed to 'doPSCN'
## 2013-05-16
## o Example code now embedded in a 'require()' statement to avoid
##   problems in the R CMD check mechanism of R-forge.
## 2013-02_27
## o By default alpha=0.01.
## 2013-01-15
## o Added argument 'platform'.
## 2013-01-09
## o Replaced 'jump' by 'bkp'
## 2013-01-03
## o Now only performs segmentation: dynamic programming is handled by
##   'jointSeg' or 'PSSeg'.
## o Added flavor 'PSCN' segmentation.
## 2013-01-01
## o Added ad hoc fix to run PSCN at full resolution (w/o discarding CN
##   probes).
## 2012-12-30
## o Added example.
## 2012-12-27
## o Some code and doc cleanups.
## 2012-12-15
## o Created.
############################################################################

