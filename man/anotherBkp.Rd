\name{anotherBkp}
\alias{anotherBkp}
\title{Get best candidate change point}
\description{Get best candidate change point according to binary segmentation}
\usage{anotherBkp(Y, weightFUN = defaultWeights, verbose = FALSE)}
\arguments{
  \item{Y}{A \code{n*p} matrix, \code{p} signals of length \code{n} to be
segmented (centered by column)}
  \item{weightFUN}{A \code{function} returning a \code{(n-1)*1} vector of weights for
the candidate change point positions. Default weights yield the
likelihood ratio test (LRT) statistic for the identification of a
single change point.}
  \item{verbose}{A \code{logical} value: should extra information be output ?
Defaults to \code{FALSE}.}
}
\details{Contrary to \code{oneBkp}, \code{anotherBkp} handles
missing values (NA:s).}


\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
p <- 2
n <- 100

sim <- randomProfile(n, 1, 1, p)
Y <- sim$profile
bkp <- anotherBkp(Y)
print(bkp)
print(oneBkp(Y))
##  stopifnot(identical(oneBkp(Y), bkp))
plotSeg(Y, list(sim$bkp, bkp))

## robustness to NA:s
h <- 2
idxs <- seq(from=max(sim$bkp[1]-h, 1), min(sim$bkp[1]+h, n))
Y[idxs, p] <- NA
oneBkp(Y)  ## does not work
bkp <- anotherBkp(Y)  ## works
bkp-sim$bkp
}

\keyword{internal}
