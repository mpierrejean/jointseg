\name{plotSeg}
\alias{plotSeg}
\title{Plot signal and breakpoints with segment-level signal estimates}
\description{Plot signal and breakpoints with segment-level signal estimates}
\usage{plotSeg(dat, breakpoints = NULL, regNames = NULL, exclNames = c("genotype", 
    "region"), ylabs = colnames(dat), ylims = NULL, binExclPattern = "^b$", 
    col = "#33333333", pch = 19, cex = 0.3)}
\arguments{
  \item{dat}{A  \code{matrix} or data frame whose rows correspond to loci sorted along the genome, or a \code{numeric} \code{vector}.}
  \item{breakpoints}{A vector of breakpoints positions, or a \code{list} of such vectors.}
  \item{regNames}{Region labels, a vector of length \code{length(breakpoints)+1} (if
\code{breakpoints} is a vector) or of length
\code{length(breakpoints[[1]])+1} (if \code{breakpoints} is a
list).                              }
  \item{exclNames}{A vector of column names corresponding to columns that should not
be plotted.}
  \item{ylabs}{A vector of 'y' labels (column names or indices) that should be plotted.}
  \item{ylims}{An optional \eqn{2*d} matrix with \code{ylim} values for each of the \eqn{d} dimensions to be plotted.}
  \item{binExclPattern}{A vector of column names or indices in \code{colnames(dat)} for which
segment-level signal estimates should *not* be drawn.}
  \item{col}{Color of plotting symbol, see \code{\link{par}}}
  \item{pch}{Plotting symbol, see \code{\link{par}}}
  \item{cex}{Magnification factor for plotting symbol, see \code{\link{par}}}
}
\details{Argument 'binCols' is mainly used to avoid
calculating mean levels for allelic ratios, which would not make
sense are they are typically multimodal.}


\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=1)
sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)
dat <- sim$profile
res <- PSSeg(dat, method="RBS", stat=c("c", "d"), K=50)
bkpList <- list(true=sim$bkp, est=res$bestSeg)
plotSeg(dat, breakpoints=bkpList)
}
