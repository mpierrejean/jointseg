\name{segmentByCnaStruct}
\alias{segmentByCnaStruct}
\title{# Run cnaStruct segmentation,}
\description{This function is a wrapper for convenient use of the \code{CnaStruct}
segmentation method by \code{\link{PSSeg}}.  It applies the
\code{breakpoints} function from package \code{CnaStruct} and reshapes
the results.}
\usage{segmentByCnaStruct(Y, K, p = 0.5, homthr = 0.95, beta = 1, s = 1, 
    maxk = 1000, ..., verbose = FALSE)}
\arguments{
  \item{Y}{A n*p matrix of signals to be segmented}
  \item{K}{The maximum number of segments to find,}
  \item{p}{Order of Minkowski distances}
  \item{homthr}{homthr sets the BAF limit to consider a SNP homozygous}
  \item{beta}{beta controls the relative importance of BAF with respect to LRR}
  \item{s}{BIC criterion for the best model if s=1}
  \item{maxk}{maxk is the maximum segment length allowed}
  \item{\dots}{
}
  \item{verbose}{
}
}


\references{Mosen-Ansorena, D. and Aransay, A.-M.(2013)
Bivariate segmentation of SNP-array data for allele-specific
copy number analysis in tumour samples
BMC Bioinformatics, 14:84  }
\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
if (require("CnaStruct")) {
  ## load known real copy number regions
  IlluDat <- loadCnRegionData(platform="Illumina", tumorFraction=0.5)
  
  ## generate a synthetic CN profile
  K <- 3
  len <- 1000
  sim <- getCopyNumberDataByResampling(len, K, minLength=10, regData=IlluDat)
  datS <- sim$profile
  
  ## run CnaStruct segmentation
  Y <- as.matrix(subset(datS, select=c(c,b)))
  Y <- cbind(Y, d=2*abs(datS$b-1/2))
  res <- jointSeg:::segmentByCnaStruct(Y, K = 3*10, maxk=500)  
  getTprTnr(res$bkp, sim$bkp, nrow(datS), 5,relax = -1)
  plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
}
}
