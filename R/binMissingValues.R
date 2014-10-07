binMissingValues <- structure(function(
###Perform binning in order to remove missing values
    Y,
### A numeric matrix
    verbose=FALSE
### A \code{logical} value: should extra information be output ?
### Defaults to \code{FALSE}.
    ) {
    ##details<<Some segmentation methods (in particular, GFLars) do not
    ##natively handle the situation when some observations have missing
    ##values in one or more dimensions. In order to avoid dropping the
    ##corresponding observations entirely, \code{binMissingValues} bins
    ##the signal values of the last complete observation before a (range
    ##of) observations with missing entries using the
    ##\code{\link[matrixStats]{binMeans}} function.
    
    ##references<< Bleakley, K., & Vert, J. P. (2011). The group fused
    ##lasso for multiple change-point detection. arXiv preprint
    ##arXiv:1106.4199.
    
    ##references<< Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
    ##change-points shared by many signals using group LARS. Advances in
    ##Neural Information Processing Systems, 23, 2343-2351.

    ##note<<Currently this function is only used by
    ##\code{\link{doGFLars}} in order to make it possible to run GFLars
    ##segmentation on SNP array data where most markers (on the order of
    ##2/3 to 5/6) have missing values, because of uninformative or
    ##missing allelic ratio signals.

    ##note<<The \code{binMissingValues} function may be used for other
    ##segmentation methods suffering from the same limitation.  However,
    ##we emphasize that handling missing values natively in the
    ##segmentation method would be a better solution.
    n <- nrow(Y)
    idxs <- 1:n
    nbNA <- rowSums(is.na(Y))
    wNA <- which(nbNA>0)

    if (wNA[1]==1) {  
        ##details<<In the specific case when the first row has NA values,
        ##the first non-missing entry is replicated in order to make
        ##smoothing possible.  This choice is arbitrary but some arbitrary
        ##choice is needed in that case.
        if (verbose) {
            print("First row has NA values")
        }
        idx <- min(which(nbNA==0))
        y <- colMeans(Y[1:(idx-1),, drop=FALSE], na.rm=TRUE)
        col <- which(is.nan(y))
        y[col] <- Y[idx, col]
        Y[1:(idx-1), ] <-  matrix(y, nrow=idx-1, ncol=ncol(Y), byrow=TRUE)
        Y[2:(idx-1), col] <- NA
    }

    idxsOK <- union(1, idxs[-wNA]) ## indices of non-missing entries (+ '1')
    xOut <- union(idxsOK, n+1)
    
    Ys <- apply(Y, 2, binMeans, x=idxs, bx=xOut)
    attr(Ys, "idxs") <- idxsOK
    Ys  ## A matrix with the same columns as the input argument \code{Y}
    ## and as many rows as non-missing entries in \code{Y}.
}, ex=function(){
    sim <- randomProfile(10, 1, 0.1, 3)
    Y <- sim$profile
    Y[c(4, 8), 2] <- NA
    Y[c(7, 8), 3] <- NA

    res <- binMissingValues(Y)

    Y <- sim$profile
    Y[1:5, 2] <- NA
    Yb <- binMissingValues(Y)

    Y <- sim$profile
    Y[3:5, 2] <- NA
    Yb <- binMissingValues(Y)
})


############################################################################
## HISTORY:
## 2014-05-15
## o Created. 
############################################################################

