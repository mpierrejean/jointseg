\name{estimateSd}
\alias{estimateSd}
\title{Robust standard deviation estimator}
\description{Estimate standard deviation of an unimodal signal with possible
changes in mean}
\usage{estimateSd(y, method = c("Hall", "von Neumann"))}\arguments{
  \item{y}{A numeric vector}
  \item{method}{Method used to estimate standard deviation}
}
\details{von Neumann's estimator is proportional to the mean absolute
deviation (\code{mad}) of the first-order differences of the
original signals: \code{mad(diff(y)}.  By construction this
estimator is robust to 1) changes in the mean of the signal
(through the use of differences) and 2) outliers (through the use
of \code{mad} instead of \code{mean}).

The proportionality constant \eqn{1/\sqrt 2 \times
1/\Phi^{-1}(3/4)} ensures that the resulting estimator is
consistent for the estimation of the standard deviation in the
case of Gaussian signals.

Hall's estimator is a weigthed sum of squared elements of y. Let m=3.
\eqn{sigma^2 = (\sum_{k=1}^{n-m}\sum_{j=1}^{m+1}(\code{wei[i]}\code{y}[i+k])^2)/(n-m)}}

\references{Von Neumann, J., Kent, R. H., Bellinson, H. R., &
Hart, B. T. (1941). The mean square successive difference. The
Annals of Mathematical Statistics, 153-162.

Peter Hall, J. W. Kay and D. M. Titterington (1990).
Asymptotically Optimal Difference-Based Estimation of Variance in
Nonparametric Regression
Biometrika,521-528}
\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
n <- 1e4
y <- rnorm(n)  ## a signal with no change in mean
estimateSd(y)
estimateSd(y, method="von Neumann")
sd(y)
mad(y)

z <- y + rep(c(0,2), each=n/2)  ## a signal with *a single* change in mean
estimateSd(z)
estimateSd(z, method="von Neumann")
sd(z)
mad(z)

z <- y + rep(c(0,2), each=100)  ## a signal with many changes in mean
estimateSd(z)
estimateSd(z, method="von Neumann")
sd(z)
mad(z)  
}
