\name{PSSeg}
\alias{PSSeg}
\title{Parent-Specific copy number segmentation}
\description{This function splits (bivariate) copy number signals into
parent-specific (PS) segments using recursive binary segmentation}
\usage{PSSeg(data, flavor = c("RBS", "GFLars", "PSCN", "cghseg", "CBS", 
    "PSCBS", "CnaStruct", "PELT"), statistic = c("c,d|het", "(c,d)|het", 
    "c", "d|het"), jitter = NULL, methModelSelection = "Birge", 
    DP = TRUE, ..., profile = FALSE, verbose = FALSE)}
\arguments{
  \item{data}{Data frame containing the following columns: \describe{
\item{c}{Total copy number}
\item{b}{Allele B fraction (a.k.a. BAF)}
\item{genotype}{(germline) genotype of the SNP, coded as 0 for AA, 1/2 for AB, 1 for BB}
}
These data are assumed to be ordered by genome position.}
  \item{flavor}{A \code{character} value, the type of segmentation method used:
\describe{
\item{"RBS"}{Recursive Binary Segmentation (the default), see
\code{\link{segmentByRBS}}}
\item{"GFLars"}{Group fused LARS as described in Bleakley and
Vert (2011).}
\item{"PSCN"}{Hidden Markov Model proposed by Chen et al (2011)}
\item{"cghseg"}{Univariate pruned dynamic programming Rigail et al (2010)}
\item{"PSCBS"}{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation by Olshen A. et al
(2011)}
\item{"CnaStruct"}{Bivariate segmentation of SNP-array data for allele-specific copy number analysis in tumour samples by Mosen-Ansorena D. et al
(2013)}
\item{"PELT"}{Optimal detection of changepoints with a linear computational cost by  Killick R. et al
(2012)}}}
  \item{statistic}{Statistic to be segmented                            }
  \item{jitter}{Uncertainty on breakpoint position after initial segmentation.  See \code{\link{jointSeg}} for details.}
  \item{methModelSelection}{Which method is used to do the model selection}
  \item{DP}{If DP =False, model selection is done on initial segmentation, else model selection is done on segmentation after dynamic programming for flavor RBS}
  \item{\dots}{Further arguments to be passed to \code{jointSeg}}
  \item{profile}{Trace time and memory usage ?}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{Before segmentation, the input copy number data are
converted to a matrix whose contents depends on the value of
argument \code{statistic}.  By default, i.e. when
\code{(statistic=="c,d|het")}, this matrix contains two columns:
\describe{
\item{c}{total copy numbers, at the "full" resolution}
\item{d}{the decrease in heterozygosity DH defined in Bengtsson et al, 2010,
\code{2|b-1/2|} in order to take advantage of the bimodality of
allele B fractions for heterozygous SNPs.  \code{d} is only
defined for heterozygous SNPs, that is, SNPs for which
\code{data$genotype==1/2}. The rationale for the above "mirror"
transformation is that allele B fractions (\code{b}) are only
informative for heterozygous SNPs (see e.g. Staaf et al, 2008).}
}

The resulting data are then segmented using the
\code{\link{jointSeg}} function for flavors \code{RBS} and
\code{GFLars}, and using the \code{\link[PSCN]{PSCN}} package for
flavor \code{PSCN}.}

\references{Bengtsson, H., Neuvial, P., & Speed,
T. P. (2010). TumorBoost: Normalization of allele-specific tumor
copy numbers from a single pair of tumor-normal genotyping
microarrays. BMC bioinformatics, 11(1), 245.

Staaf, J., Lindgren, D., Vallon-Christersson, J.,
Isaksson, A., Goransson, H., Juliusson, G., ... & Ringn\'er,
M. (2008). Segmentation-based detection of allelic imbalance and
loss-of-heterozygosity in cancer cells using whole genome SNP
arrays. Genome Biol, 9(9), R136.

Chen, H., Xing, H., & Zhang, N. R. (2011). Estimation
of parent specific DNA copy number in tumors using high-density
genotyping arrays. PLoS computational biology, 7(1), e1001060.

Bleakley, K., & Vert, J. P. (2011). The group fused
lasso for multiple change-point detection. arXiv preprint
arXiv:1106.4199.

Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
change-points shared by many signals using group LARS. Advances in
Neural Information Processing Systems, 23, 2343-2351.

Rigaill, G. (2010). Pruned dynamic programming for
optimal multiple change-point detection. arXiv preprint
arXiv:1004.0887.}
\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{jointSeg}}}
\examples{
## load known real copy number regions
affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=0.5)

## generate a synthetic CN profile
K <- 10
len <- 1e5
sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
datS <- sim$profile

## run binary segmentation (+ dynamic programming)
resRBS <- PSSeg(data=datS, K=2*K, profile=TRUE)
resRBS$prof
##getTprTnr(resRBS$initSeg$bkp, sim$bkp, nrow(datS), 10)
getTprTnr(resRBS$bestBkp, sim$bkp, nrow(datS), 10, relax = -1)
plotSeg(datS, breakpoints=list(sim$bkp, resRBS$bestBkp))
}
