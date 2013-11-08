getTpFp <- function(#Calculate the number of true positives and false positives
### Calculate the number of true positives and false positives among candidate breakpoints
                      candidates,
### Breakpoints found by the methods
                      trueBkp,
### True breakpoints
                      tol,
### Tolerance on the position of candidate breakpoints called true
    	      	     relax=c(-1,0,1)
### Controls the way multiple breapoints within tolerance area are recorded. 
### \describe{
### \item{1}{count one true positive if there is more than one breakpoint in
###   tolerance area}
### \item{0}{count one true positive only if there is exactly one breakpoint
###   in tolerance area}
### \item{1}{count only one true positive if there is exactly one breakpoint
###   in tolerance area; other breakpoints are counted as false positives
### }
                    ){
  trueBkp <- sort(trueBkp)

  ## TODO: discard regions that should not be taken into account in the
  ## evaluation because they are too small
  minRegionSize <- 2*tol+1
  diffs <- diff(c(-tol, trueBkp, Inf))
  wH1 <- which(abs(diffs)>=minRegionSize)
  if (length(wH1)<length(trueBkp)+1) {
    warning("Some breakpoints are too close for the chosen tolerance")
  }
  goodHits <- numeric(length(trueBkp))
  badHits <- numeric(length(trueBkp)+1)

  for (cc in candidates) {
    distC <- abs(trueBkp-cc)
    minC <- min(distC)
    idx <- which(distC==minC)
    if (minC<=tol) {  ## True positive
      goodHits[idx] <- goodHits[idx]+1
    } else {          ## False positive
      idx <- idx+(trueBkp[idx]<cc)
      badHits[idx] <- badHits[idx]+1
    }
  }
  addedFP <- 0
  if(relax==1){
    TP <- sum(goodHits>0)
  }else if(relax==0){
    TP <- sum(goodHits==1)
  }else if(relax==-1){
    TP <- sum(goodHits>0)
    addedFP <- sum(goodHits)-TP
  }
  FP <- sum(badHits)+addedFP
  c(TP=TP, FP=FP)
###  \item{TP}{The number of true positives} 
###  \item{FP}{The number of false positives} 
}

############################################################################
## HISTORY:
## 2013-08-04
## o Created from "getTprTnr.R".
############################################################################

