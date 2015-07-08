\name{multiplyXtXBySparse}
\alias{multiplyXtXBySparse}
\title{multiplyXtXBySparse}
\description{Compute \code{C = t(X)*X*val} , where \code{val} is a row-sparse
\code{(n-1)*p} matrix and \code{X} is the \code{n*(n-1)} design
matrix for the weighted group fused lasso.}
\usage{multiplyXtXBySparse(n, ind, val, w = defaultWeights(n), verbose = FALSE)}
\arguments{
  \item{n}{Size of the problem}
  \item{ind}{a*1 vector of indices of the non-zero rows of b (each in [1,n-1])}
  \item{val}{a*p matrix whose rows are the non-zero rows of b (same order as ind)}
  \item{w}{(n-1)*1 vector of weights}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}}
}
\details{This implementation is derived from the MATLAB code
of Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.}
\value{\item{C}{The \code{(n-1)*p} matrix equal to \code{t(X)*X*val}}}
\references{Bleakley, K., & Vert, J. P. (2011). The group fused
lasso for multiple change-point detection. arXiv preprint
arXiv:1106.4199. \url{arxiv.org/arXiv:1106.4199}

Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
change-points shared by many signals using group LARS. Advances in
Neural Information Processing Systems, 23, 2343-2351.}
\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
val <- matrix(c(1.56, 1.35, 1.26, 1.15), ncol=2)
ind <- c(5,6)
n <- 10
res <- multiplyXtXBySparse(n=n, ind=ind, val=val)
res
##           [,1]      [,2]
## [1,] 0.8874235 0.7329904
## [2,] 1.3311352 1.0994855
## [3,] 1.7428651 1.4395645
## [4,] 2.1737347 1.7954524
## [5,] 2.6622704 2.1989711
## [6,] 2.6237347 2.1787857
## [7,] 2.1036678 1.7469149
## [8,] 1.6067028 1.3342283
## [9,] 1.0711352 0.8894855
}

\keyword{internal}
