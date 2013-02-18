\name{segmentByGFLars}
\alias{segmentByGFLars}
\title{Solve the group fused Lasso optimization problem using group LARS}
\description{Solve the group fused Lasso optimization problem using group LARS}
\usage{segmentByGFLars(Y, K, weightFUN = defaultWeights, epsilon = 1e-09, 
    verbose = FALSE)}
\arguments{
  \item{Y}{A \code{n*p} signal to be segmented}
  \item{K}{The number of change points to find}
  \item{weightFUN}{A \code{function} returning a \code{(n-1)*1} vector of weights for
the candidate change point positions. Default weights yield the
likelihood ratio test (LRT) statistic for the identification of a
single change point.}
  \item{epsilon}{Values smaller than epsilon are considered null. Defaults to \code{1e-9}.}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{The default weights \eqn{\sqrt{n/(i*(n-i))}} are calibrated as
suggested by Bleakley and Vert (2011).  Using this calibration,
the first breakpoint maximizes the likelihood ratio test (LRT)
statistic.

Each dimension of the original signal is scaled before
segmentation, using \code{\link{estimateSd}}.}
\value{\item{lambda}{The estimated lambda values for each change-point}
\item{bkp}{The successive change-point positions (1*k)}
\item{value[[i]]}{A i*p matrix of change-point values for the first i change-points}
\item{c}{A n}}
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
sim <- randomProfile(1e5, trueK, 1, p)
Y <- sim$profile
K <- 2*trueK
res <- segmentByGFLars(Y, K, verbose=TRUE)
print(res$bkp)
print(sim$bkp)
plotSeg(Y, res$bkp)
}
