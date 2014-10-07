\name{leftMultiplyByXt}
\alias{leftMultiplyByXt}
\title{leftMultiplyByXt}
\description{Compute X'*Y where X is the n*(n-1) design matrix for the weighted group fused Lasso, with weights defined by the vector w, and Y is any n*p matrix. The computation is done in O(np).}
\usage{leftMultiplyByXt(Y, w = defaultWeights(nrow(Y)), verbose = FALSE)}\arguments{
\item{Y}{A n*p matrix}
\item{w}{(n-1)*1 vector of weights}
\item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{This implementation is derived from the MATLAB code
of Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.\\

Contrary to \code{\link{getUnivStat}} it does not handle missing values.}
\value{\item{C}{The (n-1)*p matrix equal to X'*Y}}
\references{Bleakley, K., & Vert, J. P. (2011). The group fused
lasso for multiple change-point detection. arXiv preprint
arXiv:1106.4199.\\
Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
change-points shared by many signals using group LARS. Advances in
Neural Information Processing Systems, 23, 2343-2351.}
\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
Y <- matrix(rnorm(20), ncol=2)
C <- leftMultiplyByXt(Y)
}

\keyword{internal}
