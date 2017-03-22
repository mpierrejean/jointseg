
#' Extract endpoint matrix from DP result
#' 
#' @param path the path vector of the "colibri_sn_R_c C" function
retour_sn <- function(path){
  n <- ncol(path)
  res3 <- matrix(NA, nrow=nrow(path), ncol=nrow(path))
  res3[1, 1] <- 0
  for(i in 2: nrow(path)){
    res3[i, i-1] <- path[i, n]
    for(k in 1:(i-1)){
      res3[i, i-1-k] <- path[i-k, res3[i, i-k]]
    }
    
  }
  diag(res3) <- ncol(path)
  return(res3)
}

#' Call the pruned DPA algorithm, use functional pruning and segment neighborhood. 
#' 
#' @param x A vector of double : the signal to be segmented
#' @param Kmax Max number of segments
#' @param mini Min value for the mean parameter of the segment
#' @param maxi Max value for the mean parameter of the segment
#' @return A list with a vector containing the position of the change-points
#' @details The pruned DPA is implemented for the L2-loss function
#' @author Guillem Rigaill
#' @export
Fpsn <- function(x, Kmax,  mini=min(x), maxi=max(x)){
  n <- length(x)
  A <- .C("colibri_sn_R_c", signal=as.double(x), n=as.integer(n), 
		Kmax=as.integer(Kmax),   min=as.double(mini), 
		max=as.double(maxi), path=integer(Kmax*n), cost=double(Kmax), allCost=double(n*Kmax)
	, PACKAGE="jointseg")
    A$path <- matrix(A$path, nrow=Kmax, byrow=TRUE)
    A$allCost <- matrix(A$allCost, nrow=Kmax, byrow=TRUE)
    A$t.est <- retour_sn(A$path)
    A$K <- length(A$t.est)
    A$J.est <- A$cost #+ sum(x^2)
    return(A);	
} 









