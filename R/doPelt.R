doPelt <- structure(function(## Run Pelt segmentation,
### This function is a wrapper for convenient use of the \code{Pelt}
### segmentation method by \code{\link{PSSeg}}.  It applies the
### \code{PELT.mean.norm} function from package \code{changepoint} and reshapes
### the results.
                                    y,
### A numeric vector, the signal to be segmented
                                    ...,
                                    verbose = FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                                    ){
  ##references<<Hinkley, D. V. (1970).
  ##Inference About the Change-Point in a Sequence
  ##of Random Variables.Biometrika 57,1-17
  ##references<<Killick R, Fearnhead P, Eckley IA (2012)
  ##Optimal detection of changepoints with
  ##a linear computational cost, JASA 107(500),1590-1598  

  if (!require("changepoint")) {
    cat("Please install the 'changepoint' package to run the 'doPelt' function")
    return()
  }
  if (!is.null(dim(y)) || mode(y)!="numeric") {
    stop("Argument 'y' should be a numeric vector")
  }
    
  cpt <- changepoint::cpt.mean(y,method="PELT")
  res <- list(bkp=cpt@cpts[-length(cpt@cpts)])
  return(res)
},ex=function(){
  ## load known real copy number regions
  affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
  
  ## generate a synthetic CN profile
  K <- 10
  len <- 1e5
  sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
  datS <- sim$profile

  ## run Pelt segmentation
  res <- doPelt(datS[["c"]])
  getTpFp(res$bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
  plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
})
############################################################################
## HISTORY:
## 2014-02-13
## o Change segment function due to update in "changepoint" package.
## o Remove last change point return by function in "changepoint" package which is the length of vector.
## 2013-12-09
## o Renamed to 'doPelt'
## 2013-03-27
## o Created
############################################################################

