getUnivJ <- function(# Get the contribution of one dimension to the RSE.
### Get the contribution of one dimension to the Residual Squared Error (RSE)
                     y,##<<Input vector of signals for the considered dimension
                     candCP
### A vector of candidate change points
                     ) {
  ##keyword<<internal

  ##details<<This function is used internally by \code{\link{pruneByDP}}.
  ##seealso<<\code{\link{pruneByDP}}
  isNotNA <- !is.na(y)
  idxsR <- c(0, cumsum(isNotNA))  ## index of original data in result
  ## idxsR <- idxsR[-length(idxsR)]

  ww <- which(isNotNA)
  y <- y[ww]
  n <- length(y)
  
  ## Compute boundaries of the smallest intervals considered
  cand <- sort(c(0, idxsR[1+candCP], n))
  k <- length(candCP) + 1 # number of intervals  

  S <- c(0, cumsum(y))
  V <- c(0, cumsum(y^2))

  JJ <- matrix(numeric(k*k), ncol=k)
  for (ii in 1:k) {
    Istart <- cand[ii] +1
    idxs <- seq(from=ii, to=k)
    Iend <- cand[idxs+1]
    JJ[ii, idxs] <- V[Iend+1] - V[Istart] - (S[Iend+1] - S[Istart])^2/(Iend-Istart+1)
  }
  return(JJ)
}

############################################################################
## HISTORY:
## 2013-01-30
## o Now an internal function (not exported anymore).
## 2012-12-XX
## o Created.
############################################################################

