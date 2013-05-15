\name{segmentByCBS}
\alias{segmentByCBS}
\title{Run CBS segmentation}
\description{This function is a wrapper for convenient use of the \code{CBS}
segmentation method by \code{\link{PSSeg}}.  It applies the
\code{\link[DNAcopy]{segment}} function and reshapes the results}
\usage{segmentByCBS(y, ..., verbose = FALSE)}
\arguments{
  \item{y}{A numeric vector, the signal to be segmented}
  \item{\dots}{Arguments to be passed to \code{\link[DNAcopy]{segment}}}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}

\value{\item{bkp}{breakpoint positions}}

\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link[DNAcopy]{segment}}}
\examples{
if (require("DNAcopy")) {
  ## load known real copy number regions
  affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
  
  ## generate a synthetic CN profile
  K <- 10
  len <- 1e5
  sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
  datS <- sim$profile
  
  ## run CBS segmentation
  res <- jointSeg:::segmentByCBS(datS[["c"]])
  getTprTnr(res$bkp, sim$bkp, nrow(datS), 5,relax = -1)
  plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
}
}
