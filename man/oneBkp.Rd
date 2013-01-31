\name{oneBkp}
\alias{oneBkp}
\title{Get best candidate change point}
\description{Get best candidate change point according to binary segmentation}
\usage{oneBkp(Y, weights = NULL, verbose = FALSE)}
\arguments{
  \item{Y}{A \code{n*p} matrix, \code{p} signals of length \code{n} to be
segmented (centered by column)}
  \item{weights}{a \code{(n-1)*1} vector of weights for the candidate change point
positions. Default weights yield the likelihood ratio test (LRT)
statistic for the identification of a single change point.}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}



\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
p <- 2
sim <- randomProfile(1e5, 1, 1, p)
Y <- sim$profile
bkp <- jointSeg:::oneBkp(Y)
par(mfrow=c(p,1))
for (ii in 1:p) {
  plot(Y[, ii], pch=19, cex=0.2)
  abline(v=bkp, col=3)
  abline(v=sim$bkp, col=8, lty=2)
}
}

\keyword{internal}
