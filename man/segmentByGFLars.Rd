\name{segmentByGFLars}
\alias{segmentByGFLars}
\title{Group fused Lars segmentation (low-level)}
\description{Low-level function for multivariate fused Lars segmentation (GFLars)}
\usage{segmentByGFLars(Y, K, weights = defaultWeights(nrow(Y)), epsilon = 1e-09, 
    verbose = FALSE)}\arguments{
  \item{Y}{A \code{n*p} matrix of signals to be segmented}
  \item{K}{The number of change points to find}
  \item{weights}{A \code{(n-1)*1} vector of weights for the weigthed group
fused Lasso penalty. See Details.}
  \item{epsilon}{Values smaller than epsilon are considered null. Defaults to \code{1e-9}.}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{This function recrusively looks for the best candidate
change point according to group-fused LARS. This is a low-level
function. It is generally advised to use the wrapper
\code{\link{doGFLars}} which also works on data frames, has a
convenient argument \code{stat}, and includes a basic workaround
for handling missing values.

See also \code{\link{jointSeg}} for combining group fused
LARS segmentation with pruning by dynamic programming
(\code{\link{pruneByDP}}).

See \code{\link{PSSeg}} for segmenting genomic signals
from SNP arrays.

The default weights \eqn{\sqrt{n/(i*(n-i))}} are calibrated as
suggested by Bleakley and Vert (2011).  Using this calibration,
the first breakpoint maximizes the likelihood ratio test (LRT)
statistic.}
\value{A list with elements:
\item{bkp}{A vector of \code{k} candidate change-point positions}
\item{lambda}{The estimated lambda values for each change-point}
\item{mean}{A vector of length \code{p}, the mean signal per column}
\item{value}{A \code{i x p} matrix of change-point values for the first i change-points}
\item{c}{\eqn{\hat{c}}, a \code{n-1 x K} matrix
}}
\references{Bleakley, K., & Vert, J. P. (2011). The group fused
lasso for multiple change-point detection. arXiv preprint
arXiv:1106.4199.

Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
change-points shared by many signals using group LARS. Advances in
Neural Information Processing Systems, 23, 2343-2351.}
\author{Morgane Pierre-Jean and Pierre Neuvial}
\note{This implementation is derived from the MATLAB code
by Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.}


\seealso{\code{\link{PSSeg}}, \code{\link{jointSeg}}, \code{\link{doGFLars}}, \code{\link{pruneByDP}}}
\examples{
p <- 2
trueK <- 10
sim <- randomProfile(1e4, trueK, 1, p)
Y <- sim$profile
K <- 2*trueK
res <- segmentByGFLars(Y, K)
print(res$bkp)
print(sim$bkp)
plotSeg(Y, res$bkp)
}
