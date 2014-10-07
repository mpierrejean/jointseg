\name{pruneByDP}
\alias{pruneByDP}
\title{Exact segmentation of a multivariate signal using dynamic programming.}
\description{Exact segmentation of a multivariate signal using dynamic programming.}
\usage{pruneByDP(Y, candCP = 1:(nrow(Y) - 1), K = length(candCP), allowNA = TRUE, 
    verbose = FALSE)}\arguments{
  \item{Y}{A \code{n*p} signal to be segmented}
  \item{candCP}{A vector of candidate change point positions (defaults to 1:(n-1))}
  \item{K}{The maximum number of change points to find}
  \item{allowNA}{A boolean value specifying whether missing values should be allowed or not.}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{This function retrieves the maximum likelihood solution
of the gaussian homoscedastic change model into segments, for
\eqn{K \in {1 \dots length(candCP)}}. The dynamic programming algorithm used
is quadratic in time. For signals containing more than 1000
points, we recommend using a first pass segmentation (see
\code{\link{segmentByRBS}}) to find a smaller number of
candidates, and to run \code{pruneByDP} on these candidates only,
as initially suggested by Gey and Lebarbier (2008). These two
steps can be performed using \code{\link{jointSeg}} for generic
multivariate signals, and using \code{\link{PSSeg}} for copy
number signals from SNP array data.

if \code{allowNA}, the calulation of the cost of removing
candidate breakpoints between i and j for i<j tolerates missing
values. Method \code{!allowNA} is maintained in order to check
consistency with the original dynamic programming in the absence
of NA:s.}
\value{A list with elements:
\item{bkpList}{A list of vectors of change point positions for the best model with k change points, for k=1, 2, ... K}
\item{rse}{A vector of K+1 residual squared errors}
\item{V}{V[i,j] is the best RSE for segmenting intervals 1 to j
}}
\references{Bellman, R. (1961). On the approximation of curves by
line segments using dynamic programming. Communications of the
ACM, 4(6), 284.

Gey, S., & Lebarbier, E. (2008). Using CART to Detect
Multiple Change Points in the Mean for Large
Sample. http://hal.archives-ouvertes.fr/hal-00327146/}
\author{Morgane Pierre-Jean and Pierre Neuvial}
\note{This implementation is derived from the MATLAB code
by Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.}


\seealso{\code{\link{jointSeg}}, \code{\link{PSSeg}}}
\examples{
p <- 2
trueK <- 10
sim <- randomProfile(1e4, trueK, 1, p)
Y <- sim$profile
K <- 2*trueK
res <- segmentByRBS(Y, K)
resP <- pruneByDP(Y, res$bkp)

##   Note that all of the above can be dmethod=="other"one directly using 'jointSeg'
resJ <- jointSeg(sim$profile, method="RBS", K=K)
stopifnot(identical(resP$bkpList, resJ$dpBkp))

## check consistency when no NA
resP2 <- pruneByDP(Y, res$bkp, allowNA=FALSE)
max(abs(resP$rse-resP2$rse))

plotSeg(Y, list(resP$bkp[[trueK]], sim$bkp), col=1)
}
