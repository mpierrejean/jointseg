\name{randomProfile}
\alias{randomProfile}
\title{Generate a random multi-dimensional profile with breakpoints and noise}
\description{Generate a random multi-dimensional profile with breakpoints and noise}
\usage{randomProfile(length, nBkp, noiseLevel, dim, minLength = 0)}\arguments{
  \item{length}{length of the profile}
  \item{nBkp}{number of breakpoints}
  \item{noiseLevel}{variance of the signal between two breakpoints}
  \item{dim}{dimension of the profile}
  \item{minLength}{minimum length of region between breakpoints by default minLength = 0}
}
\details{Generate a random profile (vector) of length \code{length},
with \code{nBkp} breakpoints randomly chosen. Between two breakpoints, the
profile is constant, uniformly chosen between 0 and 1, and a
Gaussian noice of variance \code{noiseLevel} is added.}
\value{a \code{list} that contains two fields: \describe{
\item{profile}{the profile (a \code{length} by \code{dim} matrix)}
\item{bkp}{the list of breakpoints positions (the last position at the
left of a breakpoint)}}}

\author{Morgane Pierre-Jean and Pierre Neuvial}
\note{This implementation is derived from the MATLAB code
by Vert and Bleakley: \url{http://cbio.ensmp.fr/GFLseg}.}



\examples{
len <- 1e4
nBkp <- 1e2
noiseLevel <- 1
dim <- 2

sim <- randomProfile(len, nBkp, noiseLevel, dim)
res <- doGFLars(sim$profile, K=5*nBkp)
str(res)
}
