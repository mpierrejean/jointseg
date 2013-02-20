segmentByPSCBS <- structure(function(#Run Paired PSCBS segmentation
### This function is a wrapper for convenient use of the \code{PSCBS}
### segmentation method by \code{\link{PSSeg}}.  It applies the
### \code{\link[PSCBS]{segmentByPairedPSCBS}} function and reshapes the results
                                    Y,
### A n*p signal to be segmented
                                   ...,
### Arguments to be passed to \code{\link[PSCBS]{segmentByPairedPSCBS}}
                                    verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
) {
  ##seealso<<\code{\link[PSCBS]{segmentByPairedPSCBS}}
  
  if (!require("PSCBS")) {
    cat("Please install the 'DNAcopy' package to run the 'segmentByPSCBS' function")
    return()
  }
  if (!is.matrix(Y)){
    stop("Y is not a matrix, please check dimension of Y")
  }
  n <- as.numeric(nrow(Y))
  p <- dim(Y)[2]
  chrom <- rep(1, n)
  x <- 1:n
  genomdat <- cbind(CT= Y[,"c"], betaT=Y[,"b"], muN= Y[,"genotype"])
  data <- data.frame(genomdat,x = x)
  fit <- PSCBS::segmentByPairedPSCBS(data,seed = 48879,tbn=FALSE) ##tbn =FALSE permit to use muN and no betaN 
  res <- PSCBS::getSegments(fit, simplify = TRUE)
  bkp <- round(res$start[-1],0)
  res <- list(bkp=bkp) 
  return(res)
###  \item{bkp}{breakpoint positions}
}, ex=function(){
  ## load known real copy number regions
  affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
  sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)

  ## generate a synthetic CN profile
  K <- 10
  len <- 1e5
  sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
  datS <- sim$profile

  ## run PSCBS segmentation
  Y <- as.matrix(subset(datS, select = c(c,b,genotype)))
  res <- segmentByPSCBS(Y)
  getTprTnr(res$bkp, sim$bkp, nrow(datS), 5)
  plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
})

############################################################################
## HISTORY:
## 2013-02-18
## o Bug fix
## 2013-01-09
## o Replace all jumps by bkp
## 2013-01-04
## o Created from 'segmentByCBS'.
############################################################################

