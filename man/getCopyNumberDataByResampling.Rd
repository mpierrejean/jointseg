\name{getCopyNumberDataByResampling}
\alias{getCopyNumberDataByResampling}
\title{Generate a copy number profile by resampling input data}
\description{Generate a copy number profile by resampling input data}
\usage{getCopyNumberDataByResampling(length, nBkp = NA, bkp = NULL, 
    regData = NULL, regions = NULL, regAnnot = NULL, minLength = 0, 
    regionSize = 0, connex = TRUE)}
\arguments{
  \item{length}{length of the profile}
  \item{nBkp}{number of breakpoints.  If \code{NULL}, then argument \code{bkp} is expected to be provided. }
  \item{bkp}{a numeric vector of breakpoint positions that may be used to
bypass the breakpoint generation step.  Defaults to \code{NULL}.}
  \item{regData}{a data.frame containing copy number data for different types of copy number regions.  Columns:\describe{
\item{c}{Total copy number}
\item{b}{Allele B fraction (a.k.a. BAF)}
\item{region}{a character value, annotation label for the
region. See Details.}
\item{genotype}{the (germline) genotype of SNPs. By definition, rows with
missing genotypes are interpreted as non-polymorphic loci (a.k.a. copy
number probes).}
}}
  \item{regions}{a character vector of region labels that may be used to bypass the
region label generation step.  Defaults to \code{NULL}.}
  \item{regAnnot}{a data.frame containing annotation data for each copy number region.  Columns:
\describe{
\item{region}{label of the form (must match \code{regData[["region"]]}).}
\item{freq}{frequency (in [0,1]) of this type of region in the genome.}
}
If \code{NULL} (the default), frequencies of regions (0,1), (0,2), (1,1) and (1,2) (the most common alterations) are set to represent 90\% of the regions.
\code{sum(regAnnot[["freq"]])} should be 1.}
  \item{minLength}{minimum length of region between breakpoints.  Defaults to 0.}
  \item{regionSize}{If \code{regionSize>0}, breakpoints are included by pairs, where the
distance within pair is set to \code{regionSize}.  \code{nBkp}
is then required to be an even number.}
  \item{connex}{If \code{TRUE}, any two successive regions are constrained to be
connex in the (minor CN, major CN) space.  See 'Details'.}
}
\details{This function generates a random copy number profile of
length 'length', with 'nBkp' breakpoints randomly chosen. Between two
breakpoints, the profile is constant and taken among the different
types of regions in \code{regData}.

Elements of \code{regData[["region"]]} must be of the form
\code{"(C1,C2)"}, where \code{C1} denotes the minor copy number
and \code{C2} denotes the major copy number.  For example,
\describe{
\item{(1,1)}{Normal}
\item{(0,1)}{Hemizygous deletion}
\item{(0,0)}{Homozygous deletion}
\item{(1,2)}{Single copy gain}
\item{(0,2)}{Copy-neutral LOH}
\item{(2,2)}{Balanced two-copy gain}
\item{(1,3)}{Unbalanced two-copy gain}
\item{(0,3)}{Single-copy gain with LOH}
}

If 'connex' is set to TRUE (the default), transitions
between copy number regions are constrained in such a way that for
any breakpoint, one of the minor and the major copy number does
not change.  Equivalently, this means that all breakpoints can be
seen in both total copy numbers and allelic ratios.}
\value{\item{profile}{the profile (a \code{length} by \code{2} data.frame
containing the same fields as the input data \code{regData}.}
\item{bkp}{a vector of bkp positions (the last row index before a
breakpoint)}
\item{regions}{a character vector of region labels}}

\author{Morgane Pierre-Jean and Pierre Neuvial}




\examples{
affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=1)
sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)
plotSeg(sim$profile, sim$bkp)

## another run with identical parameters
bkp <- sim$bkp
regions <- sim$regions
sim2 <- getCopyNumberDataByResampling(1e4, regData=affyDat, bkp=bkp, regions=regions)
plotSeg(sim2$profile, bkp)

## change tumor fraction but keep same "truth"
affyDatC <- loadCnRegionData(dataSet="GSE29172", tumorFraction=0.5)
simC <- getCopyNumberDataByResampling(1e4, regData=affyDatC, bkp=bkp, regions=regions)
plotSeg(simC$profile, bkp)

## restrict to only normal, single copy gain, and copy-neutral LOH
## with the same bkp
affyDatR <- subset(affyDat, region \%in\% c("(1,1)", "(0,2)", "(1,2)"))
simR <- getCopyNumberDataByResampling(1e4, regData=affyDatR, bkp=bkp)
plotSeg(simR$profile, bkp)

## Same 'truth', on another dataSet
regions <- simR$regions
illuDat <- loadCnRegionData(dataSet="GSE11976", tumorFraction=1)
sim <- getCopyNumberDataByResampling(1e4, regData=illuDat, bkp=bkp, regions=regions)
plotSeg(sim$profile, sim$bkp)
}
