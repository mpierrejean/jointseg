\name{getUnivJ}
\alias{getUnivJ}
\title{Get the contribution of one dimension to the RSE.}
\description{Get the contribution of one dimension to the Residual Squared Error (RSE)}
\usage{getUnivJ(y, candCP)}\arguments{
  \item{y}{Input vector of signals for the considered dimension}
  \item{candCP}{A vector of candidate change points}
}
\details{This function is used internally by \code{\link{pruneByDP}}.}


\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{pruneByDP}}}


\keyword{internal}
