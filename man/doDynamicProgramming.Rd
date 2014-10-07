\name{doDynamicProgramming}
\alias{doDynamicProgramming}
\title{Run segmentation by dynamic programming}
\description{High-level function for univariate or multivariate segmentation by dynamic programming}
\usage{doDynamicProgramming(Y, K, stat = NULL, verbose = FALSE)}\arguments{
\item{Y}{A numeric vector or a matrix, the signal to be segmented}
\item{K}{The number of change points to find}
\item{stat}{A vector containing the names or indices of the columns of \code{Y} to be segmented}
\item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{if the signal is uni-dimensional, this function simply
uses the segmentation method provided in the
\code{cghseg} package reshapes the results.

if the signal is multi-dimensional, this function
applies the \code{\link{pruneByDP}} function and reshapes the
results.}
\value{\item{bkp}{A vector of \code{K} indices for candidate change points}
\item{dpseg}{A list of two elements \describe{
\item{bkp}{A list of vectors of change point positions for the best model with k change points, for k=1, 2, ... K}
\item{rse}{A vector of K+1 residual squared errors}
}}}
\references{Rigaill, G. (2010). Pruned dynamic programming for
optimal multiple change-point detection. arXiv preprint
arXiv:1004.0887.}
\author{Morgane Pierre-Jean and Pierre Neuvial}
\note{This is essentially a wrapper for convenient segmentation
by dynamic programming using the \code{\link{PSSeg}} function.}



\examples{
## load known real copy number regions
affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=1)

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
datS[which(datS$genotype!=0.5),"d"] <- NA
Y = cbind(datS$c,datS$d)
resDP2d <- doDynamicProgramming(Y, K = K)
}
