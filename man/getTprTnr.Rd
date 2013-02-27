\name{getTprTnr}
\alias{getTprTnr}
\title{Calculate the proportion of true positives and true negatives}
\description{Calculate the proportion of true positives and true negatives among candidate breakpoints}
\usage{getTprTnr(candidates, trueBkp, len, tol, relax = TRUE)}
\arguments{
  \item{candidates}{Breakpoints found by the methods}
  \item{trueBkp}{True breakpoints}
  \item{len}{length of profile}
  \item{tol}{tolerance on the position of candidate breakpoints called true}
  \item{relax}{count one true positive if there is more than one breakpoint in tolerance area if relax = TRUE, count only if there is exactly one bkp in tolerance area if relax=FALSE  }
}

\value{\item{TPR}{The true positive rate}
\item{TNR}{The true negative rate} }

\author{Morgane Pierre-Jean and Pierre Neuvial}





