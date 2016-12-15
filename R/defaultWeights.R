#' Compute default weights for the weighted group fused Lasso
#' 
#' Compute default weights for the weighted group fused Lasso
#' 
#' 
#' @param n Number of observations
#' @return Vector of default weights in the reference article.
#' @note This implementation is derived from the MATLAB code by Vert and
#' Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.
#' @author Morgane Pierre-Jean and Pierre Neuvial
#' @references Bleakley, K., & Vert, J. P. (2011). The group fused lasso for
#' multiple change-point detection. arXiv preprint arXiv:1106.4199.
#' @keywords internal
#' @examples
#' 
#' defaultWeights(10)
#' 
#' @export defaultWeights
defaultWeights <- function(n){
    a <- seq(length=n-1)/n
    b <- a*(1-a)
    1/sqrt(b*n)
}

############################################################################
## HISTORY:
## 2013-01-30
## o Now an internal function (not exported anymore).
## o Added 'jointSeg:::' to example because function is not exported anymore.
## 2012-12-30
## o Avoid integer underflow by changing the order of operations.
## 2012-12-27
## o Some code and doc cleanups.
## 2012-09-13
## o Some code cleanups.
## 2012-08-13
## o Created.
############################################################################

