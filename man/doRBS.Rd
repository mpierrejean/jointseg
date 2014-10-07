\name{doRBS}
\alias{doRBS}
\title{Run RBS segmentation}
\description{High-level function for multivariate recursive binary (RBS) segmentation}
\usage{doRBS(Y, K, stat = NULL, ..., verbose = FALSE)}\arguments{
  \item{Y}{A \code{n*p} signal to be segmented}
  \item{K}{The number of change points to find}
  \item{stat}{A vector containing the names or indices of the columns of \code{Y} to be segmented}
  \item{\dots}{Further arguments to be passed to 'segmentByRBS'}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{This function is a wrapper around the lower-level
segmentation function \code{\link{segmentByRBS}}. It can be run
on p-dimensional, piecewise-constant data in order to defined a
set of candidate change points. It is recommended to prune this
list of candidates using dynamic programming
(\code{\link{pruneByDP}}), combined with a selection of the best
number of change points. The \code{\link{jointSeg}} function
provides a convenient wrapper for performing segmentation, pruning
and model selection.}
\value{An object of the same structure as the output of \code{\link{segmentByRBS}}}
\references{Gey, S., & Lebarbier, E. (2008). Using CART to Detect
Multiple Change Points in the Mean for Large
Sample. http://hal.archives-ouvertes.fr/hal-00327146/}
\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{PSSeg}}, \code{\link{pruneByDP}}}
\examples{
p <- 2
trueK <- 10
len <- 5e4
sim <- randomProfile(len, trueK, 1, p)
Y <- sim$profile
K <- 2*trueK
res <- doRBS(Y, K)
getTpFp(res$bkp, sim$bkp, tol=10, relax = -1)   ## true and false positives

cols <- rep(2, K)
cols[1:trueK] <- 3
par(mfrow=c(p,1))
for (ii in 1:p) {
    plot(Y[, ii], pch=19, cex=0.2)
    abline(v=res$bkp[1:trueK], col= cols)
    abline(v=sim$bkp, col=8, lty=2)
}

## NA:s in one dimension at a true breakpoint
jj <- sim$bkp[1]
Y[jj-seq(-10, 10), p] <- NA
res2 <- doRBS(Y, K)
getTpFp(res2$bkp, sim$bkp, tol=10, relax = -1)   ## true and false positives

## NA:s in both dimensions at a true breakpoint
Y[jj-seq(-10, 10), ] <- NA
res3 <- doRBS(Y, K)
getTpFp(res3$bkp, sim$bkp, tol=10, relax = -1)   ## true and false positives
}
