filenames <- sprintf("%s,b=%s.xdr", simName, 1:B)
pathnames <- file.path(spath, filenames)
simList <- lapply(pathnames, loadObject)
for (mm in seq(along=methTags)) {
  methTag <- methTags[mm]
  print(methTag)
  filenames <- sprintf("%s,b=%s,%s.xdr", simNameNF, 1:B, methTag)
  pathnames <- file.path(bpath, filenames)
  there <- sapply(pathnames, file.exists)
  idxs <- which(there)
  if (sum(there)==0) {
    cat("No result for method ", methTag, "\n")
    cat("Skipping this method\n")
    next
  }
  print(sum(!there))
  if (sum(!there)>0) {
    cat("Missing results for method ", methTag, ": b=", which(!there), "\n")
    cat("Skipping these\n")
    resList <- list()
    resList[idxs] <- lapply(pathnames[idxs], loadObject)
  } else {
    resList <- lapply(pathnames, loadObject)
  }

  rm(pathnames)

  ## Best segmentation according to model selection or to PSCN
  for (tol in tols) {
    print(tol)
    filename <- sprintf("%s,B=%s,%s,bestModel,tol=%s,relax=%s.xdr", simNameNF, B, methTag, tol,relax)
    pathname <- file.path(epath, filename)
    if (!file.exists(pathname) || evalForce) {
      bestMat <- NULL
      for (bb in idxs) {
        truth <- simList[[bb]]$bkp
        res <- resList[[bb]]
        if (length(grep("PSCN|CBS", methTag))) {
          best <- getTpFp(res$initBkp, truth, tol, relax)
        } else {
          best <- getTpFp(res$bestBkp, truth, tol, relax)
        }
        bestMat <- rbind(bestMat, best)
      }
      saveObject(bestMat, file=pathname)
    }
    ## dyn prog path
    filename <- sprintf("%s,B=%s,%s,rocArray,tol=%s,relax=%s.xdr", simNameNF, B, methTag, tol,relax)
    pathname <- file.path(epath, filename)

    if (!file.exists(pathname) || evalForce) {
      ## ad hoc: take maximum common length for easier aggregation
      lens <- sapply(resList, FUN=function(res) {length(res$dpBkpList)})
      wPos <- which(lens!=0)  ## exclude experiments with no bkp detected
      if (length(wPos)<B) {
        warning("No breakpoints detected in ", B-length(wPos), " experiments")
        
      }
      ll <- max(lens[wPos])
      rocArray <- array(dim=c(B, ll, 2),
                        dimnames=list(b=1:B, K=1:ll, stat=c("FP", "TP")))
      for (bb in wPos) {
        bkp <- resList[[bb]]$dpBkpList
        rBB <- sapply(bkp, getTpFp, simList[[bb]]$bkp, tol,relax)
        rocArray[bb, 1:length(bkp), "TP"] <- rBB["TP",]
        rocArray[bb, 1:length(bkp), "FP"] <- rBB["FP",]
      }
      rocArray <- rocArray[wPos,, ]
      saveObject(rocArray, file=pathname)
    }

    ## before DP (only if DP !)
    gg <- grep("+DP", methTag)
    if (length(gg)==0) next
    filename <- sprintf("%s,B=%s,%s,rocArray,init,tol=%s,relax=%s.xdr", simNameNF, B, methTag, tol,relax)
    pathname <- file.path(epath, filename)
    if (!file.exists(pathname) || evalForce) {
      ## ad hoc: take maximum common length for easier aggregation
      lens <- sapply(resList, FUN=function(res) {length(res$initBkp)})
      wPos <- which(lens!=0)  ## exclude experiments with no bkp detected
      if (length(wPos)<B) {
        warning("No breakpoints detected in ", B-length(wPos), " experiments")
      }
      ll <- max(lens[wPos])
      rocArray <- array(dim=c(B, ll, 2),
                        dimnames=list(b=1:B, K=1:ll, stat=c("FP", "TP")))
     
      for (bb in wPos) {
        bkps <- resList[[bb]]$initBkp
        for (kk in 1:length(bkps)) {
          rbk <- getTpFp(bkps[1:kk], simList[[bb]]$bkp, tol,relax)
          rocArray[bb, kk, "TP"] <- rbk["TP"]
          rocArray[bb, kk, "FP"] <- rbk["FP"]
        }
      }
      rocArray <- rocArray[wPos,, ]
      saveObject(rocArray, file=pathname)
    }
  }
}
