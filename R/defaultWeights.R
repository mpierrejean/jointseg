defaultWeights <- structure(function(#Compute default weights for the weighted group fused Lasso
### Compute default weights for the weighted group fused Lasso
                                     n
### Number of observations
                                     ){
    ##keyword<<internal
    
    ##references<< Bleakley, K., & Vert, J. P. (2011). The group fused
    ##lasso for multiple change-point detection. arXiv preprint
    ##arXiv:1106.4199.
    
    ##note<<This implementation is derived from the MATLAB code
    ##by Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.
    a <- seq(length=n-1)/n
    b <- a*(1-a)
    1/sqrt(b*n)
### Vector of default weights in the reference article.
}, ex=function() {
    defaultWeights(10)
})

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

