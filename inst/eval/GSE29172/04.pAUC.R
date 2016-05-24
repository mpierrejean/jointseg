## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## summarization of ROC curves : AUC cumputing
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ComputeAUC <- function(rocArray) {
  FPs <- 0:FPSup
  TPs <- sapply(FPs, FUN=function(FP) {
    apply(rocArray, 1, FUN=function(roc) {
      ww <- which(roc[, "FP"]<=FP)
      max(roc[ww, "TP"])
    })
  })
  TPs[is.infinite(TPs)] <- NaN
  denom <- sum(lintegrate(c(0, 0, FPSup), c(0, K, K), xint=c(0, 0, FPSup)))    
  aucs <-  apply(TPs, 1, FUN=function(tp) {
    y = c(0, tp)
    x = c(0, FPs)
    sum(lintegrate(x, y, xint=x)) 
  })
  res <- aucs/denom
  res
}

for (mm in seq(along=methTags)) {
  print(methTags[mm])
  bestMatList <- NULL
  rocArrayList <- NULL
  rocArrayInitList <- NULL
  aucData <- NULL
  aucArray <- NULL
  methTag <- methTags[mm]
  for (tol in tols) {
    ## dyn prog path
    filename <- sprintf("%s,B=%s,%s,rocArray,tol=%s,relax=%s.xdr", simNameNF, B, methTag, tol, relax)
    pathname <- file.path(epath, filename)
    print(file.exists(pathname))
    if(file.exists(pathname)){ 
      rocArray <- loadObject(file=pathname)## RocArray for tol and methTag
      aggFUN <- eval(as.name(sprintf("ComputeAUC")))
      auc <- aggFUN(rocArray)
      ## roc curves
      masomeno <-  2*sd(auc)/sqrt(B)
      aucMean <- mean(auc)
      aucData <- rbind(aucData,c(aucMean, masomeno))
      aucArray <- rbind(aucArray,auc)
    }
  }## End tol
  AucPath <- "aucData"
  fpath <- file.path(AucPath, simName)
  fpath <- Arguments$getWritablePath(fpath)
  filename <- sprintf("%s,B=%s,%s,aucArray,relax=%s.xdr", simNameNF, B, methTag, relax)
  pathname <- file.path(fpath, filename)
  if(!is.null(aucArray)){
    saveObject(aucArray, pathname)
  }
  fpath <- file.path(AucPath, simName)
  fpath <- Arguments$getWritablePath(fpath)
  filename <- sprintf("%s,B=%s,%s,aucData,relax=%s.xdr", simNameNF, B, methTag, relax)
  pathname <- file.path(fpath, filename)
  if(!is.null(aucData)){
    dimnames(aucData) <- list(NULL, c("meanAUC", "masOmenoAUC"))
    saveObject(aucData, pathname)
  }
  ## aucData contains for method methTag, a table with mean of Auc for various tolerance
}
