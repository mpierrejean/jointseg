\name{getTpFp}
\alias{getTpFp}
\title{Calculate the number of true positives and false positives}
\description{Calculate the number of true positives and false positives among candidate breakpoints}
\usage{getTpFp(candidates, trueBkp, tol, relax = -1)}\arguments{
\item{candidates}{Breakpoints found by the methods}
\item{trueBkp}{True breakpoints}
\item{tol}{Tolerance on the position of candidate breakpoints called true}
\item{relax}{Controls the way multiple breapoints within tolerance area are recorded. 
\describe{
\item{1}{count one true positive if there is at least one breakpoint within
tolerance area}
\item{0}{count one true positive only if there is exactly one breakpoint
within tolerance area}
\item{-1}{count only one true positive if there is exactly one breakpoint
within tolerance area; other breakpoints are counted as false positives
}}}
}

\value{A list with elements:
\item{TP}{The number of true positives}
\item{FP}{The number of false positives
}}

\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
## load known real copy number regions
affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=0.7)

## generate a synthetic CN profile
K <- 10
len <- 2e4
sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
datS <- sim$profile

## (group-)fused Lasso segmentation
res <- PSSeg(data=datS, K=2*K, method="GFLars", stat="c", profile=TRUE)

## results of the initial (group-)fused lasso segmentation
getTpFp(res$initBkp, sim$bkp, tol=10, relax=-1)
getTpFp(res$initBkp, sim$bkp, tol=10, relax=0)
getTpFp(res$initBkp, sim$bkp, tol=10, relax=1)
plotSeg(datS, breakpoints=list(sim$bkp, res$initBkp))

## results after pruning (group-)fused Lasso candidates by dynamic programming)
getTpFp(res$bestBkp, sim$bkp, tol=10, relax=-1)
getTpFp(res$bestBkp, sim$bkp, tol=10, relax=0)
getTpFp(res$bestBkp, sim$bkp, tol=10, relax=1)
plotSeg(datS, breakpoints=list(sim$bkp, res$bestBkp))
}
