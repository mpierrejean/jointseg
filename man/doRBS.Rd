\name{doRBS}
\alias{doRBS}
\title{Recursive Binary Segmentation }
\description{High-level function for multivariate Recursive Binary Segmentation (RBS)}
\usage{doRBS(Y, K, stat = NULL, ..., verbose = FALSE)}
\arguments{
  \item{Y}{A \code{vector}, \code{matrix} or \code{data.frame}) containing signals to be segmented}
  \item{K}{The number of change points to find}
  \item{stat}{A vector containing the names or indices of the columns of \code{Y} to be segmented}
  \item{\dots}{Further aguemnts to be passed to 'segmentByRBS'.}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{This function is a wrapper aroung the lower-level
segmentation function \code{\link{segmentByRBS}}. It can be run on
p-dimensional, piecewise-constant data in order to defined a set
of candidate change points. It is recommended to prune this list
of candidates using dynamic programming (\code{\link{pruneByDP}}),
combined with a selection of the best number of change points. The
\code{\link{jointSeg}} function provides a convenient wrapper for
performing segmentation, pruning and model selection.

For the specific case of DNA copy number data segmentation, see the
dedicated wrapper \code{\link{PSSeg}}.}
\value{\item{bkp}{A \code{vector} of \code{K} estimated breakpoint positions, ranked
by order of appearance}
\item{gain}{The gain provided by the breakpoints in terms of difference between RSE} }
\references{Gey, S., & Lebarbier, E. (2008). Using CART to Detect
Multiple Change Points in the Mean for Large
Sample. http://hal.archives-ouvertes.fr/hal-00327146/}
\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{segmentByRBS}}, \code{\link{PSSeg}},
\code{\link{pruneByDP}}}
\examples{
p <- 2
trueK <- 10
len <- 1e4
sim <- randomProfile(len, trueK, 1, p)
Y <- sim$profile
K <- 2*trueK
res <- doRBS(Y, K)
getTpFp(res$bkp, sim$bkp, tol=10)   ## true and false positives

cols <- rep(2, K)
cols[1:trueK] <- 3
par(mfrow=c(p,1))
for (ii in 1:p) {
  plot(Y[, ii], pch=19, cex=0.2)
  abline(v=res$bkp[1:trueK], col= cols)
  abline(v=sim$bkp, col=8, lty=2)
}
}
