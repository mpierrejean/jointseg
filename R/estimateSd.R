estimateSd <- structure(function(#Robust standard deviation estimator
### Estimate standard deviation of an unimodal signal with possible
### changes in mean
                                 y,
### A numeric vector
                                 method='HALL'
### Method used to estimate standard deviation c("VonNeuman","HALL")
                                 ) {
  ##details<< VonNeuman's estimator is proportional to the mean absolute
  ##deviation (\code{mad}) of the first-order differences of the
  ##original signals: \code{mad(diff(y)}.  By construction this
  ##estimator is robust to 1) changes in the mean of the signal
  ##(through the use of differences) and 2) outliers (through the use
  ##of \code{mad} instead of \code{mean}).

  ##details<< The proportionality constant \eqn{1/\sqrt 2 \times
  ##1/\Phi^{-1}(3/4)} ensures that the resulting estimator is
  ##consistent for the estimation of the standard deviation in the
  ##case of Gaussian signals.

  ##references<<Von Neumann, J., Kent, R. H., Bellinson, H. R., &
  ##Hart, B. T. (1941). The mean square successive difference. The
  ##Annals of Mathematical Statistics, 153-162.

  ##details<<HALL's estimator is based on square weighted sum of y. Let m = 3.
  ##\eqn{sigma^2 = (\sum_{k=1}^{n-m}\sum_{j=1}^{m+1}(\code{wei[i]}\code{y}[i+k])^2)/(n-m)}
  
  ##references<<Peter Hall, J. W. Kay and D. M. Titterington (1990).
  ##Asymptotically Optimal Difference-Based Estimation of Variance in
  ##Nonparametric Regression
  ##Biometrika,521-528
  
  y <- na.omit(y)
  if(!(method%in%c("VonNeuman","HALL"))){
    stop("Argument method should be VonNeuman or HALL")
  }
  if(method == "VonNeumann"){ 
    dy <- diff(y)
    Sd <- mad(dy)/sqrt(2)  ## Note: 'mad' is already scaled to be consistent for Gaussian signals.c("VonNeuman","HALL")
### The estimated standard deviation
  }
  if(method == "HALL"){
    n = length(y)
    wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
    mat <- wei %*% t(y)
    mat[2, -n] <- mat[2, -1]
    mat[3, -c(n-1, n)] <- mat[3, -c(1, 2)]
    mat[4, -c(n-2, n-1, n)] <- mat[4, -c(1, 2, 3)]
    Sd <- sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2) / (n-3))
  }
  return (Sd)
}, ex=function() {
  n <- 1e4
  y <- rnorm(n)  ## a signal with no change in mean
  estimateSd(y)
  sd(y)
  mad(y)
  
  z <- y + rep(c(0,2), each=n/2)  ## a signal with *a single* change in mean
  estimateSd(z)
  sd(z)
  mad(z)

  z <- y + rep(c(0,2), each=100)  ## a signal with many changes in mean
  estimateSd(z)
  sd(z)
  mad(z)  
})

############################################################################
## HISTORY:
## 2013-03-19
## o Added method 'HALL'
## 2013-01-02
## o Created.
############################################################################

