\name{modelSelection}
\alias{modelSelection}
\title{Model selection}
\description{Select the best number of breakpoints}
\usage{modelSelection(rse, n, c = 2.5, lambdas = NULL, meth = c("Birge", 
    "Lebarbier"))}
\arguments{
  \item{rse}{RSE as output by \code{\link{pruneByDP}}}
  \item{n}{Length of the profile }
  \item{c}{Parameter for the model selection}
  \item{lambdas}{A list of candidate values for the calibration of the penalty}
  \item{meth}{Method to calibrate the constant in the penalty for the model selection}
}
\details{This function is not meant to be called directly, but
implicitly through \code{\link{jointSeg}} or \code{\link{PSSeg}}.}
\value{\item{kbest}{the best number of breakpoints}
\item{lambda}{result of function to select the best number of breakpoints}}
\references{Lebarbier, E. (2005). Detecting multiple change-points
in the mean of Gaussian process by model selection. Signal
processing, 85(4), 717-736

Birg\'{e}, L. (2001). Gaussian model selection. J.Eur
Math. Soc, 3(3):203-268}
\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{jointSeg}}

\code{\link{PSSeg}}}
\examples{
## load known real copy number regions
affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)
Y <- as.matrix(sim$profile[, "c"])

## Find candidate breakpoints
K <- 50
resRBS <- segmentByRBS(Y, K=K)
## Prune candidate breakpoints
resDP <- pruneByDP(Y, candCP=resRBS$bkp)
selectedModel <- modelSelection(rse=resDP$rse, n=nrow(Y), meth='Birge')
str(selectedModel)

## breakpoints of the best model
print(resDP$bkp[[selectedModel$kbest]])

## truth
print(sim$bkp)
}
