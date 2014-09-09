\name{prof}
\alias{prof}
\title{profile time and memory usage of a given R expression}
\description{profile time and memory usage of a given R expression}
\usage{prof(expr, doit = TRUE)}\arguments{
  \item{expr}{An \code{R} expression to be evaluated}
  \item{doit}{A boolean variable specifying whether profiling should be
performed or not (intended for internal use).}
}
\details{Profiling is performed using \code{summaryRprof(memory="both")$by.self}.

Memory profiling is not satistfactory yet.}


\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{Rprof}}

\code{\link{summaryRprof}}}


\keyword{internal}
