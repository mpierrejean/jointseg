\name{segmentByPSCN}
\alias{segmentByPSCN}
\title{Run PSCN segmentation}
\description{This function is a wrapper for convenient use of the \code{\link[PSCN]{PSCN}}
segmentation method by \code{\link{PSSeg}}.  It applies the
\code{\link[PSCN]{smoothing}} and \code{\link[PSCN]{segmentation}}
functions.}
\usage{segmentByPSCN(Y, alpha = 0.01, platform = c("Illumina", "Affymetrix"), 
    verbose = FALSE, ...)}
\arguments{
  \item{Y}{The signal to be segmented, a matrix containing the following columns: \describe{
\item{c}{Total copy number}
\item{b}{Allele B fraction (a.k.a. BAF)}
}}
  \item{alpha}{sensitivity level in [0,1] to be passed to
\code{\link[PSCN]{segmentation}}.}
  \item{platform}{Specifies form which array platform 'Y' was generated: Illumina or Affymetrix}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
  \item{\dots}{Further arguments to be passed to \code{\link[PSCN]{smoothing}}}
}
\details{The data are assumed to come from a single chromosome.}
\value{\item{bkp}{breakpoint positions}}
\references{Chen, H., Xing, H., & Zhang, N. R. (2011). Estimation
of parent specific DNA copy number in tumors using high-density
genotyping arrays. PLoS computational biology, 7(1), e1001060.}
\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link[PSCN]{smoothing}}

\code{\link[PSCN]{segmentation}}

\code{\link{PSSeg}}}
\examples{
## load known real copy number regions
affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)
datS <- sim$profile

## run PSCN segmentation
resPSCN <- segmentByPSCN(datS)
getTprTnr(resPSCN$bkp, sim$bkp, nrow(datS), 20, relax = -1)
plotSeg(datS, breakpoints=list(resPSCN$bkp, sim$bkp))
}
