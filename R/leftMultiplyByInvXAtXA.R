leftMultiplyByInvXAtXA <- structure(function(
### Compute r = inv(X(:,ind)'*X(:,ind))*b , where X is the n*(n-1) design matrix for the weighted group fused lasso.
                                             n,
### The size of X is n*(n-1)
                                             ind,
### a*1 vector of indices between 1 and n-1, sorted in increasing order
                                             val,
### a*p matrix
                                             w=defaultWeights(n),
### (n-1)*1 vector of weights
                                             verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                                             ){
  ##keyword<<internal

  ##details<<This implementation is derived from the MATLAB code
  ##of Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.\\

  ##references<< Bleakley, K., & Vert, J. P. (2011). The group fused
  ##lasso for multiple change-point detection. arXiv preprint
  ##arXiv:1106.4199.

  ##references<<Vert, J. P., & Bleakley, K. (2010). Fast detection of
  ##multiple change-points shared by many signals using group
  ##LARS. Advances in Neural Information Processing Systems, 23,
  ##2343-2351.
  a <- dim(val)[1]
  p <- dim(val)[2]
  o <- order(ind)
  ind <- ind[o]
  val <- val[o,, drop=FALSE]
  r <- matrix(numeric(a*p), nrow= a, ncol=p)
  if(length(w)!= (n-1)) {
    stop("'w' needs to be of length n-1")
  }
  if (a!=0){
    ## see paper for explanation of this formula
    v <- diff(c(0,ind, n))
    d <- w[ind]
    R <- matrix(numeric((a+2)*p), ncol=p)
    val <- apply(val, 2, function(x) {x/d})
    R[1,] <- numeric(p)
    R[2:(a+1), ] <- val
    R[(a+2),] <- numeric(p)
    gamma <- apply(R, 2, diff)
    delta <- gamma/v
    r <- -diff(delta)
    r <- r/d
  }
  return(r)
### \item{r}{the (n-1)*p matrix equal to X'*Y}
}, ex=function(){
  val <- matrix(c(1.56, 1.35, 1.26, 1.15), ncol=2)
  ind <- c(5,6)
  n <- 10
  res <- leftMultiplyByInvXAtXA(n=n, ind=ind, val=val)
  res
  ##         [,1]      [,2]
  ## [1,] 1.373189 0.9630868
  ## [2,] 0.228796 0.3636429
})

############################################################################
## HISTORY:
## 2013-01-30
## o Now an internal function (not exported anymore).
## o Added 'jointSeg:::' to example because function is not exported anymore.
## 2012-12-27
## o Some code and doc cleanups.
## 2012-12-06
## o Replaced the function defaultweigths by a (n-1)*1 vector of weigths 
## 2012-09-13
## o Some code cleanups.
## 2012-08-13
## o Created.
############################################################################

