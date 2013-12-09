modelSelection <- structure(function(#Model selection
### Select the best number of breakpoints
                           rse,
### RSE as output by \code{\link{pruneByDP}}
                           n,
### Length of the profile 
                           c=2.5,
### Parameter for the model selection
                           lambdas=NULL,       
### A list of candidate values for the calibration of the penalty
                           meth=c('Birge', 'Lebarbier')
### Method to calibrate the constant in the penalty for the model selection
){
  ##details<<This function is not meant to be called directly, but
  ##implicitly through \code{\link{jointSeg}} or \code{\link{PSSeg}}.
  
  ##seealso<<\code{\link{jointSeg}}
  ##seealso<<\code{\link{PSSeg}}
  
  ##references<<Lebarbier, E. (2005). Detecting multiple change-points
  ##in the mean of Gaussian process by model selection. Signal
  ##processing, 85(4), 717-736

  ##references<<Birg\'{e}, L. (2001). Gaussian model selection. J.Eur
  ##Math. Soc, 3(3):203-268

  ##If there is only one value of rse, the best model has only one segment kbest = 0
  ##If D is too small, Bige can't be used.

  D <- length(rse)  ## KMax+1 (rse[1] corresponds to 0 breakpoints)
  if(D == 1){
      return(list(kbest=0, lambda=NA))
  }
  if (is.null(lambdas)) {
    lambdas <- seq(from=0.01, to=3, length=100)
  }
  ChooseK <- function(lambda){
    pen <- lambda*(1:D)*(c+log(n/1:D))/n
    DHat <- which.min(rse/n+pen)
    return(DHat)
  }
  ERMajustment <- function(ERM,n){
    D <- length(ERM)
    ind <- (D%/%2):D
    y <- ERM[ind]/n
    var <- ind
    model <- lm(y~var)
    lambda <- -model$coefficient["var"]*n
    return(lambda)
  }
  if (meth=="Lebarbier"){
     Dhat <- sapply(lambdas, ChooseK)
     maxBkpInDhat <- which.max(-diff(Dhat))
     bestLambda <- lambdas[maxBkpInDhat]*2
  }
  if (meth=="Birge") {
    bestLambda <- ERMajustment(rse,n)
  }
  if(is.na(bestLambda)){
    warnings("Too small D, Lebarbier's method was used")
    Dhat <- sapply(lambdas, ChooseK)
    maxBkpInDhat <- which.max(-diff(Dhat))
    bestLambda <- lambdas[maxBkpInDhat]*2
  }
  bestD <- ChooseK(bestLambda)
  bestK <- bestD-1
  if (bestK>D-1) {  ## do not select a model larger than the largest candidate
    bestK <- D-1
    bestLambda <- bestLambda; ## TODO: something else for better consistency wrt K ?
  }
  if(D<= 5){
    warnings("Model selection may be not very accurate")
  }
  return(list(kbest=bestK, lambda=bestLambda))
### \item{kbest}{the best number of breakpoints}
### \item{lambda}{result of function to select the best number of breakpoints}
}, ex = function(){
  ## load known real copy number regions
  affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
  sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)
  Y <- as.matrix(sim$profile[, "c"])
  
  ## Find candidate breakpoints
  K <- 50
  resRBS <- doRBS(Y, K=K)
  ## Prune candidate breakpoints
  resDP <- pruneByDP(Y, candCP=resRBS$bkp)
  selectedModel <- modelSelection(rse=resDP$rse, n=nrow(Y), meth='Birge')
  str(selectedModel)
  
  ## breakpoints of the best model
  print(resDP$bkp[[selectedModel$kbest]])

  ## truth
  print(sim$bkp)
})

############################################################################
## HISTORY:
## 2013-01-23
## o Updated doc and example.
## o Now using 'ERMajustment' for model selection.
## 2013-01-09
## o Replace all jumps by bkp
## o BUG FIX : return the number of breakpoints now and no the number of segments
## 2012-12-30
## o Some code and doc cleanups.
## o BUG FIX: could return a best model larger than the best candidates.
## 2012-11-27
## o Created.
############################################################################

