segmentByGFLars <- structure(function(#Group fused Lars segmentation (low-level)
### Low-level function for multivariate fused Lars segmentation (GFLars)
                                      Y,
### A \code{n*p} matrix of signals to be segmented
                                      K,
### The number of change points to find
                                      weights=defaultWeights(nrow(Y)),
### A \code{(n-1)*1} vector of weights for the weigthed group
###   fused Lasso penalty. See Details.
                                      epsilon=1e-9,
### Values smaller than epsilon are considered null. Defaults to \code{1e-9}.
                                      verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                                      ){
  ##details<<This function recrusively looks for the best candidate
  ##change point according to group-fused LARS. This is a low-level
  ##function. It is generally advised to use the wrapper
  ##\code{\link{doGFLars}} which also works on data frames, has a
  ##convenient argument \code{stat}, and includes a basic workaround
  ##for handling missing values.

  ##details<<See also \code{\link{jointSeg}} for combining group fused
  ##LARS segmentation with pruning by dynamic programming
  ##(\code{\link{pruneByDP}}).

  ##details<<See \code{\link{PSSeg}} for segmenting genomic signals
  ##from SNP arrays.

  ##seealso<<\code{\link{PSSeg}}, \code{\link{jointSeg}}, \code{\link{doGFLars}}, \code{\link{pruneByDP}}
  
  ##details<<The default weights \eqn{\sqrt{n/(i*(n-i))}} are calibrated as
  ##suggested by Bleakley and Vert (2011).  Using this calibration,
  ##the first breakpoint maximizes the likelihood ratio test (LRT)
  ##statistic.

  
  ##references<< Bleakley, K., & Vert, J. P. (2011). The group fused
  ##lasso for multiple change-point detection. arXiv preprint
  ##arXiv:1106.4199.
  
  ##references<< Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
  ##change-points shared by many signals using group LARS. Advances in
  ##Neural Information Processing Systems, 23, 2343-2351.
  
  ##note<<This implementation is derived from the MATLAB code
  ##by Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.
  
  if (!is.matrix(Y) || !is.numeric(Y)){
    stop("Argument 'Y' should be a numeric matrix.\nPlease see 'doGFLars' for using GFLars directly on a data.frame or a numeric vector")
  }
  
  if (any(is.na(Y))) {
    stop("Missing values are not handled by the current implementation of group fused LARS")
  }

  rownames(Y) <- NULL
  n <- as.numeric(nrow(Y))
  p <- dim(Y)[2]

  ## Initialization
  if (missing(K)) {
    stop("Please provide argument 'K'")
  }
  if (K>=n) {
    stop("Too many breakpoints are required")
  }
  if (is.null(weights)) {
    weights <- defaultWeights(n)
  } 
  res.meansignal <- colMeans(Y);
  res.lambda <- numeric(K);
  res.bkp <- numeric(K)
  res.value <- list()
  res.c <- matrix(NA, n-1, K);
  Y <- sweep(Y, 2, colMeans(Y))  ## [PN:2013-01-02] Not needed ?? (implicitly done within 'leftMultiplyByXt')
  AS <- numeric(0)
  c <- leftMultiplyByXt(Y=Y, w=weights, verbose=verbose)   
  for (ii in 1:K){
    cNorm <- rowSums(c^2)
    res.c[, ii] <- cNorm
    bigcHat <- max(cNorm)
    
    if (verbose) {
      print(paste('optimize LARS : ', ii))
    }
    ## First breakpoint
    if (ii==1) {
      AS <- which.max(cNorm)
      res.bkp[ii] <- AS
    }
    I <- order(AS)
    AS <- AS[I]
    w <- leftMultiplyByInvXAtXA(n, AS, matrix(c[AS,], ncol=p), weights, verbose=verbose)
    a <- multiplyXtXBySparse(n=n, ind=AS, val=w, w=weights, verbose=verbose)
    a1 <- bigcHat - rowSums(a^2)
    u <- a*c
    a2 <- bigcHat - rowSums(u)
    a3 <- bigcHat - cNorm

    ## We solve it
    gammaTemp = matrix(NA, n-1, 2); # to get the indices right

    ## First those where we really have a second-order polynomial
    subset <- which(a1 > epsilon)
    delta <- a2[subset]^2 - a1[subset]*a3[subset]
    delta.neg <- subset[which(delta<0)]
    delta.pos <- subset[which(delta>=0)]
    gammaTemp[delta.neg, 1] <- NA;
    gammaTemp[delta.pos, 1] <- (a2[delta.pos]+sqrt(delta[which(delta>=0)]))/a1[delta.pos];
    gammaTemp[delta.neg, 2] <- NA;
    gammaTemp[delta.pos, 2] <- (a2[delta.pos]-sqrt(delta[which(delta>=0)]))/a1[delta.pos];

    ## Then those where the quadratic term vanishes and we have a
    ## first-order polynomial
    subset <- which((a1 <= epsilon) & (a2 > epsilon))
    gammaTemp[subset, ] = a3[subset] / (2*a2[subset])

    ## Finally the active set should not be taken into account, as well as
    ## those for which the computation gives dummy solutions
    maxg <- max(gammaTemp, na.rm=TRUE)+1;
    subset <- which((a1 <= epsilon) & (a2 <= epsilon));
    gammaTemp[subset, ] <- maxg;
    gammaTemp[AS, ] <- maxg;
    gamma <- min(gammaTemp, na.rm=TRUE)
    idx <- which.min(gammaTemp)
    nexttoadd <- 1 + (idx-1) %% (n-1);

    ## Update beta
    res.lambda[ii] <- sqrt(bigcHat);
    res.value[[ii]] <- matrix(numeric(ii*p), ncol = p)
    res.value[[ii]][I,] <- gamma*w;
    if (ii>1){
      res.value[[ii]][1:(ii-1),] <- res.value[[ii]][1:(ii-1),] + res.value[[ii-1]];
    }
    ## Update active indexes
    if (ii<K){
      AS <- c(AS,nexttoadd)
      res.bkp[ii+1] <- nexttoadd;
      c <- c-gamma*a;
    }
    ## print(AS)
  }
  ##value<< A list with elements:
  list(
      bkp=res.bkp, ##<<A vector of \code{k} candidate change-point positions
      lambda=res.lambda, ##<<The estimated lambda values for each change-point
      mean=res.meansignal, ##<<A vector of length \code{p}, the mean signal per column
      value=res.value, ##<<A \code{i x p} matrix of change-point values for the first i change-points
      c=res.c) ##<<\eqn{\hat{c}}, a \code{n-1 x K} matrix
}, ex=function(){
  p <- 2
  trueK <- 10
  sim <- randomProfile(1e4, trueK, 1, p)
  Y <- sim$profile
  K <- 2*trueK
  res <- segmentByGFLars(Y, K)
  print(res$bkp)
  print(sim$bkp)
  plotSeg(Y, res$bkp)
})

############################################################################
## HISTORY:
## 2014-05-15
## o Renamed back to 'segmentByGFLars', so that 'doGFLars' is a *wrapper*
## around the core segmentation function.
## 2013-12-09
## o Renamed to 'doGFLars'
## 2013-01-09
## o Replace all jumps by bkp
## 2012-12-27
## o Renamed to segmentByGFLars.
## o Some code and doc cleanups.
## 2012-12-13
## o Updated example.
## 2012-12-06
## o Replaced the function defaultWeigths by a (n-1)*1 vector of
## weigths and updated example
## 2012-09-13
## o Some code cleanups.
## o Tentative bug fix: indices in gammaTemp.
## o SPEEDUP: removed unnecessary calls to 'complex'.
## 2012-08-13
## o Created.
############################################################################

