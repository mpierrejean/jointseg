\name{leftMultiplyByInvXAtXA}
\alias{leftMultiplyByInvXAtXA}
\title{leftMultiplyByInvXAtXA}
\description{Compute r = inv(X(:,ind)'*X(:,ind))*b , where X is the n*(n-1) design matrix for the weighted group fused lasso.}
\usage{leftMultiplyByInvXAtXA(n, ind, val, w, verbose = FALSE)}
\arguments{
  \item{n}{The size of X is n*(n-1)}
  \item{ind}{a*1 vector of indices between 1 and n-1, sorted in increasing order}
  \item{val}{a*p matrix}
  \item{w}{(n-1)*1 vector of weights}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{This implementation is derived from the MATLAB code
of Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.\\}
\value{\item{r}{the (n-1)*p matrix equal to X'*Y}}
\references{Bleakley, K., & Vert, J. P. (2011). The group fused
lasso for multiple change-point detection. arXiv preprint
arXiv:1106.4199.

Vert, J. P., & Bleakley, K. (2010). Fast detection of
multiple change-points shared by many signals using group
LARS. Advances in Neural Information Processing Systems, 23,
2343-2351.}
\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
val <- matrix(c(1.56, 1.35, 1.26, 1.15), ncol=2)
ind <- c(5,6)
n <- 10
w <- jointSeg:::defaultWeights(n)
res <- jointSeg:::leftMultiplyByInvXAtXA(n=n, ind=ind, val=val, w=w)
res
##         [,1]      [,2]
## [1,] 1.373189 0.9630868
## [2,] 0.228796 0.3636429
}

\keyword{internal}
