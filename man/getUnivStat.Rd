\name{getUnivStat}
\alias{getUnivStat}
\title{Get the binary test statistic for one dimension}
\description{Get the binary test statistic for one dimension}
\usage{getUnivStat(y, weightFUN = defaultWeights)}\arguments{
\item{y}{Input vector of signals for the considered dimension}
\item{weightFUN}{A \code{function} returning a \code{(n-1)*1} vector of weights for
the candidate change point positions. See \code{\link{anotherBkp}}}
}
\details{This function is used internally by \code{\link{anotherBkp}}.}


\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{anotherBkp}}}


\keyword{internal}
