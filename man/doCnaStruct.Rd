\name{doCnaStruct}
\alias{doCnaStruct}
\title{# Run cnaStruct segmentation}
\description{This function is a wrapper for convenient use of the \code{CnaStruct}
segmentation method by \code{\link{PSSeg}}.  It applies the
\code{breakpoints} function from package \code{CnaStruct} and reshapes
the results.}
\usage{doCnaStruct(Y, K, p = 0.5, homthr = 0.95, beta = 1, s = 1, maxk = NULL, 
    ..., verbose = FALSE)}
\arguments{
  \item{Y}{The signal to be segmented, a matrix containing the following columns: \describe{
\item{c}{Total copy number (log scale)}
\item{b}{Allele B fraction (a.k.a. BAF)}
}}
  \item{K}{Maximum number of segments to find,}
  \item{p}{Order of Minkowski distances}
  \item{homthr}{Sets the BAF limit to consider a SNP homozygous}
  \item{beta}{Controls the relative importance of BAF with respect to LRR}
  \item{s}{BIC criterion for the best model if s=1}
  \item{maxk}{Maximum segment length allowed}
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
  Y <- cbind(c=log2(datS[, "c"])-1, b=datS[, "b"])  ## Convert to log ('LRR') scale
  res <- jointSeg:::doCnaStruct(Y, K=K*10, maxk=500)  
  getTpFp(res$bkp, sim$bkp, tol=5, relax=-1)   ## true and false positives
  plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
}
}
