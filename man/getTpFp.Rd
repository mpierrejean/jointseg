\name{getTpFp}
\alias{getTpFp}
\title{Calculate the number of true positives and false positives}
\description{Calculate the number of true positives and false positives among candidate breakpoints}
\usage{getTpFp(candidates, trueBkp, tol, relax = c(-1, 0, 1))}
\arguments{
  \item{candidates}{Breakpoints found by the methods}
  \item{trueBkp}{True breakpoints}
  \item{tol}{Tolerance on the position of candidate breakpoints called true}
  \item{relax}{Controls the way multiple breapoints within tolerance area are recorded. 
\describe{
\item{1}{count one true positive if there is more than one breakpoint in
tolerance area}
\item{0}{count one true positive only if there is exactly one breakpoint
in tolerance area}
\item{1}{count only one true positive if there is exactly one breakpoint
in tolerance area; other breakpoints are counted as false positives
}}}
}

\value{\item{TP}{The number of true positives} 
\item{FP}{The number of false positives} }

\author{Morgane Pierre-Jean and Pierre Neuvial}





