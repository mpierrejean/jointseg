getTprTnr <- function(#Calculate the proportion of true positives and true negatives
### Calculate the proportion of true positives and true negatives among candidate breakpoints
                      candidates,
### Breakpoints found by the methods
                      trueBkp,
### True breakpoints
                      len,
### length of profile
                      tol,
### tolerance on the position of candidate breakpoints called true
    	      	     relax=c(-1,0,1)
### count one true positive if there is more than one breakpoint in tolerance area if relax = 1, count only if there is exactly one bkp in tolerance area if relax=0, count only one true positive if there is exactly one bkp in tolerance area if relax=-1 other bkp are count as false positive
                      ){
  trueBkp <- sort(trueBkp)

  ## TODO: discard regions that should not be taken into account in the
  ## evaluation because they are too small
  minRegionSize <- 2*tol+1
  diffs <- diff(c(-tol, trueBkp, len+tol))
  wH1 <- which(abs(diffs)>=minRegionSize)
  ## Note: not used yet !!
  if (length(wH1)<length(trueBkp)+1) {
    warning("Some breakpoints are too close for the chosen tolerance")
  }
  goodHits <- numeric(length(trueBkp))
  badHits <- numeric(length(trueBkp)+1)

##  distC <- abs(sapply(candidates, "-", trueBkp))
##  minC <- colMins(distC)
  
  
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
  c(TPR=TP/length(trueBkp), TNR=1-FP/(len-1-length(trueBkp)))
###  \item{TPR}{The true positive rate}
###  \item{TNR}{The true negative rate} 
}

############################################################################
## HISTORY:
## 2013-03-01
## o Guillem Rigaill's suggestion , count false positive if there is more than one bkp in tolerance area.
## 2013-02-27
## o Add parameter 'relax' count a true positive if there is exactly one bkp in tolerance area
## 2012-12-30
## o BUG FIX in the calculation of false positives:
##  '(trueBkp[idx]>cc)' -> '(trueBkp[idx]<cc)'
## 2012-12-13
## o Created from "tptn,simulations.R".
############################################################################

