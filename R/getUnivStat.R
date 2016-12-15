#' Get the binary test statistic for one dimension
#' 
#' Get the binary test statistic for one dimension
#' 
#' This function is used internally by \code{\link{anotherBkp}}.
#' 
#' @param y Input vector of signals for the considered dimension
#' @param weightFUN A \code{function} returning a \code{(n-1)*1} vector of
#' weights for the candidate change point positions. See
#' \code{\link{anotherBkp}}
#' @return A vector of length \code{length(y)-1}, the binary test statistic for one dimension
#' @author Morgane Pierre-Jean and Pierre Neuvial
#' @seealso \code{\link{anotherBkp}}
#' @keywords internal
getUnivStat <- function(y, weightFUN=defaultWeights) {
    isNotNA <- !is.na(y)
    if (!(any(isNotNA))) {  ## that is, everyone is NA !
        return(rep(NA, length(y)-1))
    }
    idxsR <- cumsum(isNotNA)  ## index of original data in result
    ## ad hoc: handle beginning and end
    ## idxsR[idxsR==0] <- idxsR[length(idxsR)]  ## circular ;)
    idxsR[idxsR==0] <- NA
    idxsR <- idxsR[-length(idxsR)]

    ww <- which(isNotNA)
    y <- y[ww]
    n <- length(y)
    weights <- weightFUN(n)
    S <- cumsum(y)
    
    idxs <- seq(length=n-1)
    stat <- weights*(idxs/n*S[n]-S[idxs])
    
    stat[idxsR]
}


############################################################################
## HISTORY:
## 2013-01-30
## o Now an internal function (not exported anymore).
## 2012-12-XX
## o Created.
############################################################################

