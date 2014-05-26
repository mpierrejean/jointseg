\name{doDynamicProgramming}
\alias{doDynamicProgramming}
\title{Run cghseg segmentation}
\description{This function is a wrapper for convenient use of the \code{cghseg}
segmentation method by \code{\link{PSSeg}}.  It applies the
\code{segmeanCO} function from package \code{cghseg} and reshapes
the results.}
\usage{doDynamicProgramming(Y, K, verbose = FALSE)}
\arguments{
  \item{Y}{A numeric vector or a matrix, the signal to be segmented}
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

## generate a synthetic CN profile
K <- 10
len <- 1e4
sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
datS <- sim$profile

## run cghseg segmentation
resDP <- doDynamicProgramming(datS[["c"]], K=K)
getTpFp(resDP$bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
plotSeg(datS, breakpoints=list(sim$bkp, resDP$bkp))

## run 2d dynamic programming segmentation
K <- 2
len <- 1e3
sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
datS <- sim$profile
datS$d <- 2*abs(datS$b-1/2)
datS[which(datS$genotype!=0.5),"d"]=NA
Y = cbind(datS$c,datS$d)
resDP2d <- doDynamicProgramming(Y, K = K)

}
