\name{doPSCBS}
\alias{doPSCBS}
\title{Run Paired PSCBS segmentation}
\description{This function is a wrapper for convenient use of the \code{PSCBS}
segmentation method by \code{\link{PSSeg}}.  It applies the
\code{\link[PSCBS]{segmentByPairedPSCBS}} function and reshapes the results}
\usage{doPSCBS(Y, ..., verbose = FALSE)}\arguments{
  \item{Y}{A matrix of signals to be segmented, containing the
following columns \describe{
\item{c}{total copy numbers}
\item{b}{allele B fractions (a.k.a. BAF)}
\item{genotype}{germline genotypes}
}}
  \item{\dots}{Arguments to be passed to \code{\link[PSCBS]{segmentByPairedPSCBS}}}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}

\value{\item{bkp}{breakpoint positions}}

\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link[PSCBS]{segmentByPairedPSCBS}}}
\examples{
if (require("PSCBS") && require("aroma.light")) {
    ## load known real copy number regions
    affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=1)
    
    ## generate a synthetic CN profile
    K <- 10
    len <- 1e4
    sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
    datS <- sim$profile
    
    ## run PSCBS segmentation
    Y <- as.matrix(subset(datS, select=c(c, b, genotype)))
    res <- doPSCBS(Y)
    getTpFp(res$bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
    plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
}
}
