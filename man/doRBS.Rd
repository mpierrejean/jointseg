\name{doRBS}
\alias{doRBS}
\title{Run RBS segmentation}
\description{Segment a multivariate signal using recursive binary segmentation (RBS)}
\usage{doRBS(Y, K, minRegionSize = 2, verbose = FALSE)}
\arguments{
  \item{Y}{A \code{n*p} signal to be segmented}
  \item{K}{The number of change points to find}
  \item{minRegionSize}{Regions with less than \code{minRegionSize} are not split}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{This function recrusively looks for the best candidate
change point according to binary segmentation (as given by
\code{\link{anotherBkp}}).

For segmenting genomics signals from SNP arrays, see the
convenient wrapper \code{\link{PSSeg}} for combining recursive
binary segmentation with pruning by dynamic programming
(\code{\link{pruneByDP}}).

Each dimension of the original signal is scaled before
segmentation, using \code{\link{estimateSd}}.}
\value{\item{bkp}{A \code{vector} of \code{K} estimated breakpoint positions, ranked
by order of appearance}
\item{gain}{The gain provided by the breakpoints in terms of difference between RSE} }

\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{PSSeg}}, \code{\link{pruneByDP}}}
\examples{
p <- 2
trueK <- 10
len <- 1e5
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
