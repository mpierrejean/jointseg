multiplyXtXBySparse <- structure(function(
### Compute \code{C = t(X)*X*val} , where \code{val} is a row-sparse
### \code{(n-1)*p} matrix and \code{X} is the \code{n*(n-1)} design
### matrix for the weighted group fused lasso.
    n,
### Size of the problem
    ind,
### a*1 vector of indices of the non-zero rows of b (each in [1,n-1])
    val,
### a*p matrix whose rows are the non-zero rows of b (same order as ind)
    w=defaultWeights(n),
### (n-1)*1 vector of weights
    verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}
    ){
    ##keyword<<internal

    ##details<<This implementation is derived from the MATLAB code
    ##of Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.

    ##references<<Bleakley, K., & Vert, J. P. (2011). The group fused
    ##lasso for multiple change-point detection. arXiv preprint
    ##arXiv:1106.4199. \url{arxiv.org/arXiv:1106.4199}
    
    ##references<<Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
    ##change-points shared by many signals using group LARS. Advances in
    ##Neural Information Processing Systems, 23, 2343-2351.
    if(length(w)!= (n-1)) {
        stop("Argument 'w' has to be of length nrow(Y)-1")
    }
    a <- nrow(val)
    p <- ncol(val)
    indrev <- (n-1):1
    
    if (a!=0){
        ## Sort ind and val
        o <- order(ind)
        ind <- ind[o]
        val <- val[o, , drop = FALSE]
        ## First multiply val by the weights
        r <- val*w[ind]
        ## compute the sum S
        s <- ind*r
        s <- colSums(s)/n
        ## compute the matrix M : the cumsum of the reverse of r
        M <- matrix(numeric((n-1)*p), nrow = (n-1), ncol = p)
        M[indrev[ind],] <- r[,, drop=FALSE]
        M <- apply(M, 2, cumsum)
        M <- M[indrev,,drop=FALSE]
        ## compute the matrix U
        u <- sweep(M, 2, s)
        ## u <- M-matrix(s, n-1, p, byrow=TRUE)
        U <- apply(u, 2, cumsum)
        ## compute C
        C <- U*w
    }
    return(C)
### \item{C}{The \code{(n-1)*p} matrix equal to \code{t(X)*X*val}}
}, ex = function(){
    val <- matrix(c(1.56, 1.35, 1.26, 1.15), ncol=2)
    ind <- c(5,6)
    n <- 10
    res <- multiplyXtXBySparse(n=n, ind=ind, val=val)
    res
    ##           [,1]      [,2]
    ## [1,] 0.8874235 0.7329904
    ## [2,] 1.3311352 1.0994855
    ## [3,] 1.7428651 1.4395645
    ## [4,] 2.1737347 1.7954524
    ## [5,] 2.6622704 2.1989711
    ## [6,] 2.6237347 2.1787857
    ## [7,] 2.1036678 1.7469149
    ## [8,] 1.6067028 1.3342283
    ## [9,] 1.0711352 0.8894855
})

############################################################################
## HISTORY:
## 2013-01-30
## o Now an internal function (not exported anymore).
## o Added 'jointSeg:::' to example because function is not exported anymore.
## 2012-12-27
## o Some code and doc cleanups.
## 2012-12-06
## o Replaced the function defaultWeigths by a (n-1)*1 vector of weigths 
## 2012-09-13
## o Some code cleanups.
## 2012-08-13
## o Created.
############################################################################

