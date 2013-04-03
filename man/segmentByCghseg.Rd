\name{segmentByCghseg}
\alias{segmentByCghseg}
\title{Run cghseg segmentation}
\description{This function is a wrapper for convenient use of the \code{cghseg}
segmentation method by \code{\link{PSSeg}}.  It applies the
\code{segmeanCO} function from package \code{cghseg} and reshapes
the results.}
\usage{segmentByCghseg(y, K, verbose = FALSE)}
\arguments{
  \item{y}{A numeric vector, the signal to be segmented}
  \item{K}{The number of change points to find}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}

\value{\item{bkp}{A vector of \code{K} indices for candidate change points}
\item{dpseg}{A list of two elements \describe{
\item{bkp}{A list of vectors of change point positions for the best model with k change points, for k=1, 2, ... K}
\item{rse}{A vector of K+1 residual squared errors}
}}}
\references{Rigaill, G. (2010). Pruned dynamic programming for
optimal multiple change-point detection. arXiv preprint
arXiv:1004.0887.}
\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
## load known real copy number regions
affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)

## generate a synthetic CN profile
K <- 10
len <- 1e5
sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
datS <- sim$profile

## run cghseg segmentation
resDP <- segmentByCghseg(datS[["c"]], K=K)
getTprTnr(resDP$bkp, sim$bkp, nrow(datS), 5, relax = -1)
plotSeg(datS, breakpoints=list(sim$bkp, resDP$bkp))
}
