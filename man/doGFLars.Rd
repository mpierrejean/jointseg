\name{doGFLars}
\alias{doGFLars}
\title{Group fused Lars segmentation}
\description{High-level function for multivariate fused Lars (GFLars) segmentation}
\usage{doGFLars(Y, K, stat = NULL, ..., verbose = FALSE)}\arguments{
\item{Y}{A \code{n*p} signal to be segmented}
\item{K}{The number of change points to find}
\item{stat}{A vector containing the names or indices of the columns of \code{Y} to be segmented}
\item{\dots}{Further arguments to be passed to 'segmentByGFLars'}
\item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{This function is a wrapper around the lower-level
segmentation function \code{\link{segmentByGFLars}}. It can be run
on p-dimensional, piecewise-constant data in order to defined a
set of candidate change points. It is recommended to prune this
list of candidates using dynamic programming
(\code{\link{pruneByDP}}), combined with a selection of the best
number of change points. The \code{\link{jointSeg}} function
provides a convenient wrapper for performing segmentation, pruning
and model selection.

For the specific case of DNA copy number data segmentation, see the
dedicated wrapper \code{\link{PSSeg}}.

The default weights \eqn{\sqrt{n/(i*(n-i))}} are calibrated as
suggested by Bleakley and Vert (2011).  Using this calibration,
the first breakpoint maximizes the likelihood ratio test (LRT)
statistic.}
\value{An object of the same structure as the output of \code{\link{segmentByGFLars}}}
\references{Bleakley, K., & Vert, J. P. (2011). The group fused
lasso for multiple change-point detection. arXiv preprint
arXiv:1106.4199.

Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
change-points shared by many signals using group LARS. Advances in
Neural Information Processing Systems, 23, 2343-2351.}
\author{Morgane Pierre-Jean and Pierre Neuvial}
\note{This implementation is derived from the MATLAB code
by Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.}



\examples{
p <- 2
trueK <- 10
sim <- randomProfile(1e4, trueK, 1, p)
Y <- sim$profile
K <- 2*trueK
res <- doGFLars(Y, K)
print(res$bkp)
print(sim$bkp)
plotSeg(Y, res$bkp)

## a toy example with missing values
sim <- randomProfile(1e2, 1, 0.1, 2)
Y <- sim$profile
Y[3:50, 2] <- NA

res <- doGFLars(Y, 10, 2, verbose=TRUE)
print(res$bkp)
print(sim$bkp)
plotSeg(Y, res$bkp)
}
