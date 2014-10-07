\name{doPSCN}
\alias{doPSCN}
\title{Run PSCN segmentation (defunct)}
\description{The 'PSCN' package is not maintained anymore and it is not available for
R >= 3.0.0.  The original 'doPSCN' function has been moved to the directory
'zzz.defunct'.  The skeleton of that function is kept for backward
compatibility.}
\usage{doPSCN(Y, alpha = 0.01, platform = c("Illumina", "Affymetrix"), 
verbose = FALSE, ...)}\arguments{
\item{Y}{The signal to be segmented, a matrix containing the following columns: \describe{
\item{c}{Total copy number (log scale)}
\item{b}{Allele B fraction (a.k.a. BAF)}
}}
\item{alpha}{sensitivity level in [0,1] to be passed to
\code{PSCN::segmentation}.}
\item{platform}{Specifies form which array platform 'Y' was generated: Illumina or Affymetrix}
\item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
\item{\dots}{Further arguments to be passed to \code{PSCN::smoothing}}
}


\references{Chen, H., Xing, H., & Zhang, N. R. (2011). Estimation
of parent specific DNA copy number in tumors using high-density
genotyping arrays. PLoS computational biology, 7(1), e1001060.}
\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{PSSeg}}}
\examples{
print("The 'PSCN' package is not available for R >= 3.0.0.")
print("See http://cran.r-project.org/web/packages/PSCN/index.html")
}
