doPSCBS <- structure(function(#Run Paired PSCBS segmentation
### This function is a wrapper for convenient use of the \code{PSCBS}
### segmentation method by \code{\link{PSSeg}}.  It applies the
### \code{\link[PSCBS]{segmentByPairedPSCBS}} function and reshapes the results
                                    Y,
### A matrix of signals to be segmented, containing the
### following columns \describe{
### \item{c}{total copy numbers}
### \item{b}{allele B fractions (a.k.a. BAF)}
### \item{genotype}{germline genotypes}
### }
                                     ...,
### Arguments to be passed to \code{\link[PSCBS]{segmentByPairedPSCBS}}
                                     verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                                     ) {
  ##seealso<<\code{\link[PSCBS]{segmentByPairedPSCBS}}
  
  if (!require("PSCBS")) {
    cat("Please install the 'PSCBS' package to run the 'doPSCBS' function")
    return()
  }
  if (!is.matrix(Y)){
    stop("Y is not a matrix, please check dimension of Y")
  }

  cn <- colnames(Y)
  ecn <- c("c", "b", "genotype") ## expected
  mm <- match(ecn, cn)
  if (any(is.na(mm))) {
    str <- sprintf("('%s')", paste(ecn, collapse="','"))
    stop("Argument 'Y' should contain columns named ", str)
  }
  
  n <- as.numeric(nrow(Y))
  p <- dim(Y)[2]
  chrom <- rep(1, n)
  x <- 1:n
  genomdat <- cbind(CT=Y[, "c"], betaT=Y[, "b"], muN=Y[, "genotype"])
  data <- data.frame(genomdat, x=x)
  str(data)
  fit <- PSCBS::segmentByPairedPSCBS(data, tbn=FALSE) ##tbn=FALSE permits to use 'muN' and not 'betaN'
  res <- PSCBS::getSegments(fit, simplify=TRUE)
  bkp <- round(res$start[-1],0)
  res <- list(bkp=bkp) 
  return(res)
###  \item{bkp}{breakpoint positions}
}, ex=function(){
  if (require("PSCBS") && require("aroma.light")) {
    ## load known real copy number regions
    affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
    sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)
    
    ## generate a synthetic CN profile
    K <- 10
    len <- 1e5
    sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
    datS <- sim$profile
    
    ## run PSCBS segmentation
    Y <- as.matrix(subset(datS, select=c(c,b,genotype)))
    res <- doPSCBS(Y)
    getTpFp(res$bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
    plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
  }
})

############################################################################
## HISTORY:
## 2013-12-09
## o Renamed to 'doPSCBS'
## 2013-05-30
## o Now explicitly requiring "aroma.light" as well.
## 2013-05-16
## o Example code now embedded in a 'require()' statement to avoid
##   problems in the R CMD check mechanism of R-forge.
## 2013-02-18
## o Bug fix.
## 2013-01-09
## o Replaced 'jump' by 'bkp'.
## 2013-01-04
## o Created from 'segmentByCBS'.
############################################################################

