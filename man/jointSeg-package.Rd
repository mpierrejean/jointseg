\name{jointSeg-package}
\alias{jointSeg-package}
\docType{package}
\title{Joint segmentation of multivariate (copy number) signals}
\description{This package implements functions to quickly segment multivariate signals into piecewise constant profiles. A typical application is the joint segmentation of total DNA copy numbers and allelic ratios obtained from Single Nucleotide Polymorphism (SNP) microarrays in cancer studies.  We use a two step approach where the first step (over-)segments the signal using recursive binary segmentation, and the second step prunes candiate change points obtained at the first step using dynamic programming.}
\details{
\tabular{ll}{Package: \tab jointSeg\cr
Type: \tab Package\cr
Title: \tab Joint segmentation of multivariate (copy number) signals\cr
Version: \tab 0.3.6\cr
Date: \tab 2012-02-27\cr
Author: \tab Morgane Pierre-Jean and Pierre Neuvial\cr
Maintainer: \tab Morgane Pierre-Jean <morgane.pierrejean@genopole.cnrs.fr>\cr
License: \tab LGPL (>= 2.1)\cr
Depends: \tab R (>= 2.10), R.utils, matrixStats (>= 0.6.0), acnr (>= 0.1.2)\cr
Suggests: \tab PSCN, cghseg, DNAcopy, PSCBS, aroma.light, R.cache, digest\cr}
}
\author{Morgane Pierre-Jean and Pierre Neuvial}

\keyword{ package }


