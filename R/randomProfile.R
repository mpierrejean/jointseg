randomProfile <- structure(function(# Generate a random multi-dimensional profile with breakpoints and noise
### Generate a random multi-dimensional profile with breakpoints and noise
                             length,
### length of the profile
                             nBkp,
### number of breakpoints
                             noiseLevel,
### variance of the signal between two breakpoints
                             dim,
### dimension of the profile
                             minLength = 0
### minimum length of region between breakpoints by default minLength = 0
                             ){
  ##details<<Generate a random profile (vector) of length \code{length},
  ##with \code{nBkp} breakpoints randomly chosen. Between two breakpoints, the
  ##profile is constant, uniformly chosen between 0 and 1, and a
  ##Gaussian noice of variance \code{noiseLevel} is added.

  ##note<<This implementation is derived from the MATLAB code
  ##by Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.

  ## validate arguments
  maxMinLength <- length/nBkp/2
  if (minLength >= maxMinLength) {
    stop("Cannot squeeze ", nBkp, " breakpoints in a profile of length ", length,  " with minimal segment length of 2*", minLength)
  }

  ## First make the noise
  profile <- matrix(rnorm(length*dim, mean=0, sd=noiseLevel), ncol=dim);

  ## Choose the breakpoints
  interval <-1:(length-1)
  u <- numeric(0)
  for(i in seq(length=nBkp)){
    j <- sample(x=interval, size=1 , replace=FALSE)
    u <- c(u, j)
    b.inf <- max(1, j-minLength)
    b.sup <- min(length, j+minLength)
    v <- b.inf:b.sup
    interval <- setdiff(interval, v)
  }
  ## u <- sample(length-1, nBkp, replace=FALSE);
  bkp <- sort(u);
  bkpB <- c(0, bkp, length);

  ## Add the random piecewise linear profile
  for (ii in 1:(nBkp+1)) {
    idxsII <- seq(from=bkpB[ii]+1, to=bkpB[ii+1])
    meanII <- matrix(rnorm(dim), length(idxsII), dim, byrow=TRUE)
    profile[idxsII, ] <- profile[idxsII, ] + meanII
  }
  return(list(bkp=bkp, profile=profile))
###a \code{list} that contains two fields: \describe{
###\item{profile}{the profile (a \code{length} by \code{dim} matrix)}
###\item{bkp}{the list of breakpoints positions (the last position at the
###left of a breakpoint)}}
}, ex=function() {
  len <- 1e4
  nBkp <- 1e2
  noiseLevel <- 1
  dim <- 2

  sim <- randomProfile(len, nBkp, noiseLevel, dim)
  res <- doGFLars(sim$profile, K=5*nBkp)
  str(res)
})
############################################################################
## HISTORY:
## 2013-01-09
## o Replace all jumps by bkp
## 2012-12-27
## o Some code and doc cleanups.
## 2012-10-19
## o Added some sanity checks for arguments.
## o BUG FIX: 'minLength' would not work as expected.
## 2012-09-18
## o Added example code.
## 2012-09-13
## o Some code cleanups.
## 2012-08-13
## o Created.
############################################################################

