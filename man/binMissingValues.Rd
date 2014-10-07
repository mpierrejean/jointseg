\name{binMissingValues}
\alias{binMissingValues}
\title{binMissingValues}
\description{Perform binning in order to remove missing values}
\usage{binMissingValues(Y, verbose = FALSE)}\arguments{
  \item{Y}{A numeric matrix}
  \item{verbose}{A \code{logical} value: should extra information be output ?
Defaults to \code{FALSE}.}
}
\details{Some segmentation methods (in particular, GFLars) do not
natively handle the situation when some observations have missing
values in one or more dimensions. In order to avoid dropping the
corresponding observations entirely, \code{binMissingValues} bins
the signal values of the last complete observation before a (range
of) observations with missing entries using the
\code{\link[matrixStats]{binMeans}} function.

In the specific case when the first row has NA values,
the first non-missing entry is replicated in order to make
smoothing possible.  This choice is arbitrary but some arbitrary
choice is needed in that case.}

\references{Bleakley, K., & Vert, J. P. (2011). The group fused
lasso for multiple change-point detection. arXiv preprint
arXiv:1106.4199.

Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
change-points shared by many signals using group LARS. Advances in
Neural Information Processing Systems, 23, 2343-2351.}
\author{Morgane Pierre-Jean and Pierre Neuvial}
\note{Currently this function is only used by
\code{\link{doGFLars}} in order to make it possible to run GFLars
segmentation on SNP array data where most markers (on the order of
2/3 to 5/6) have missing values, because of uninformative or
missing allelic ratio signals.

The \code{binMissingValues} function may be used for other
segmentation methods suffering from the same limitation.  However,
we emphasize that handling missing values natively in the
segmentation method would be a better solution.}



\examples{
sim <- randomProfile(10, 1, 0.1, 3)
Y <- sim$profile
Y[c(4, 8), 2] <- NA
Y[c(7, 8), 3] <- NA

res <- binMissingValues(Y)

Y <- sim$profile
Y[1:5, 2] <- NA
Yb <- binMissingValues(Y)

Y <- sim$profile
Y[3:5, 2] <- NA
Yb <- binMissingValues(Y)
}
