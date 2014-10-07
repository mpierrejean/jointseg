leftMultiplyByXt <-structure(function(
### Compute X'*Y where X is the n*(n-1) design matrix for the weighted group fused Lasso, with weights defined by the vector w, and Y is any n*p matrix. The computation is done in O(np).
    Y,
### A n*p matrix
    w=defaultWeights(nrow(Y)),
### (n-1)*1 vector of weights
    verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
    ){
    ##keyword<<internal

    ##details<<This implementation is derived from the MATLAB code
    ##of Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.\\

    ##details<<Contrary to \code{\link{getUnivStat}} it does not handle missing values.
    
    ##references<< Bleakley, K., & Vert, J. P. (2011). The group fused
    ##lasso for multiple change-point detection. arXiv preprint
    ##arXiv:1106.4199.\\
    ##Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
    ##change-points shared by many signals using group LARS. Advances in
    ##Neural Information Processing Systems, 23, 2343-2351.
    n <- as.numeric(dim(Y)[1])
    p <- ncol(Y)
    u <- apply(Y, 2, cumsum)
    if(length(w)!= (n-1)) {
        stop("w needs to be of length nrow(Y)-1")
    }
    C <- apply(u, 2, function(x){
        w*((1:(n-1))*x[n]/n-x[1:(n-1)])
    })
    dim(C) <- c(n-1,p)  ## so that the code also works with n=2...
    return(C)
### \item{C}{The (n-1)*p matrix equal to X'*Y}
}, ex = function(){
    Y <- matrix(rnorm(20), ncol=2)
    C <- leftMultiplyByXt(Y)
})

############################################################################
## HISTORY:
## 2013-01-30
## o Now an internal function (not exported anymore).
## o Added 'jointSeg:::' to example because function is not exported anymore.
## 2012-12-27
## o Some code and doc cleanups.
## 2012-12-06
## o replaced the function defaultWeigths by a (n-1)*1 vector of weigths 
## o BUG FIX: when n=2 the result would loose its 'dim' attribute.
## 2012-09-13
## o Some code cleanups.
## 2012-08-13
## o Created.
############################################################################


