\name{doPelt}
\alias{doPelt}
\title{# Run Pelt segmentation,}
\description{This function is a wrapper for convenient use of the \code{Pelt}
segmentation method by \code{\link{PSSeg}}.  It applies the
\code{PELT.mean.norm} function from package \code{changepoint} and reshapes
the results.}
\usage{doPelt(y, ..., verbose = FALSE)}
\arguments{
  \item{y}{A numeric vector, the signal to be segmented}
  \item{\dots}{
}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}


\references{Hinkley, D. V. (1970).
Inference About the Change-Point in a Sequence
of Random Variables.Biometrika 57,1-17

Killick R, Fearnhead P, Eckley IA (2012)
Optimal detection of changepoints with
a linear computational cost, JASA 107(500),1590-1598  }
\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
## load known real copy number regions
affyDat <- loadCnRegionData(platform="GSE29172", tumorFraction=1)

## generate a synthetic CN profile
K <- 10
len <- 1e4
sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
datS <- sim$profile

## run Pelt segmentation
res <- doPelt(datS[["c"]])
getTpFp(res$bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
}
