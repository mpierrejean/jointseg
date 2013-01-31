\name{defaultWeights}
\alias{defaultWeights}
\title{Compute default weights for the weighted group fused Lasso}
\description{Compute default weights for the weighted group fused Lasso}
\usage{defaultWeights(n)}
\arguments{
  \item{n}{Number of observations}
}

\value{Vector of default weights in the reference article.}
\references{Bleakley, K., & Vert, J. P. (2011). The group fused
lasso for multiple change-point detection. arXiv preprint
arXiv:1106.4199.}
\author{Morgane Pierre-Jean and Pierre Neuvial}
\note{This implementation is derived from the MATLAB code
by Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.}



\examples{
jointSeg:::defaultWeights(10)
}

\keyword{internal}
