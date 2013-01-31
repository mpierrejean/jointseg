\name{segmentByPSCBS}
\alias{segmentByPSCBS}
\title{Run Paired PSCBS segmentation}
\description{This function is a wrapper for convenient use of the \code{PSCBS}
segmentation method by \code{\link{PSSeg}}.  It applies the
\code{\link[PSCBS]{segmentByPairedPSCBS}} function and reshapes the results}
\usage{segmentByPSCBS(Y, ..., verbose = FALSE)}
\arguments{
  \item{Y}{A n*p signal to be segmented}
  \item{\dots}{Arguments to be passed to \code{\link[PSCBS]{segmentByPairedPSCBS}}}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}

\value{\item{bkp}{breakpoint positions}}

\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link[PSCBS]{segmentByPairedPSCBS}}}
\examples{
## load known real copy number regions
affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)

## generate a synthetic CN profile
K <- 10
len <- 1e5
sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
datS <- sim$profile

## run PSCBS segmentation
## Y <- as.matrix(subset(datS, select = c(c,b,genotype)))
## res <- segmentByPSCBS(Y)
## getTprTnr(res$bkp, sim$bkp, nrow(datS), 5)
## plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
}
