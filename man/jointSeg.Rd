\name{jointSeg}
\alias{jointSeg}
\title{Joint segmentation of multivariate signals}
\description{Joint segmentation of multivariate signals in two steps:
\enumerate{
\item{first-pass segmentation.  By default, a fast, greedy
approach is used (see \code{flavor}).}
\item{pruning of the candidate change points obtained by
dynamic programming}
}}
\usage{jointSeg(Y, flavor = c("RBS", "GFLars", "PSCN", "cghseg", "CBS", 
    "PSCBS", "CnaStruct", "PELT"), jitter = NULL, methModelSelection = "Birge", 
    DP = TRUE, ..., profile = FALSE, verbose = FALSE)}
\arguments{
  \item{Y}{The signal to be segmented (must be a matrix)}
  \item{flavor}{A \code{character} value, the type of segmentation method used:
\describe{
\item{"RBS"}{Recursive Binary Segmentation (the default), see
\code{\link{segmentByRBS}}}
\item{"GFLars"}{Group fused LARS as described in Bleakley and
Vert (2011).}
\item{"PSCN"}{Hidden Markov Model proposed by Chen et al (2011).
Can only be used for copy number signals from SNP arrays}
\item{"cghseg"}{Univariate pruned dynamic programming Rigail et al
(2010)}
\item{"PSCBS"}{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation by Olshen A. et al
(2011)}
\item{"CnaStruct"}{Bivariate segmentation of SNP-array data for allele-specific copy number analysis in tumour samples by Mosen-Ansorena D. et al
(2013)}
\item{"PELT"}{Optimal detection of changepoints with a linear computational cost by  Killick R. et al
(2012)}}   }
  \item{jitter}{Uncertainty on breakpoint position after initial segmentation.  Defaults to NULL.  See Details.}
  \item{methModelSelection}{Which method is used to perform model selection}
  \item{DP}{If DP =False, model selection is done on initial segmentation, else model selection is done on segmentation after dynamic programming for flavor RBS}
  \item{\dots}{Further arguments to be passed to the lower-level segmentation
method determined by argument \code{flavor}.}
  \item{profile}{Trace time and memory usage ?}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{If \code{flavor=="PSCN"} or  \code{flavor=="CnaStruct"}, \code{Y} should contain the
following columns \describe{
\item{c}{total copy numbers}
\item{b}{allele B fractions (a.k.a. BAF)}
\item{d}{decrease in heterozygosity (\code{2|b-1/2|}) for
heterozygous SNPs}
}

If \code{flavor=="PSCBS"}, \code{Y} should contain the
following columns \describe{
\item{c}{total copy numbers}
\item{b}{allele B fractions (a.k.a. BAF)}
\item{d}{decrease in heterozygosity (\code{2|b-1/2|}) for
heterozygous SNPs}
\item{genotype}{genotype}
}

If \code{jitter} is not NULL, it should be a vector of
integer indices. The set of candidate breakpoints passed on to
dynamic programming is augmented by all indices distant from an
element of \code{jitter} from one of the candidates. For
example, if \code{jitter==c(-1, 0, 1)} and the initial set of
breakpoints is \code{c(1,5)} then dynamic programming is run on
\code{c(1,2,4,5,6)}.}
\value{list with elements:
\item{bestBkp}{Best set of breakpoints after dynamic programming}
\item{initBkp}{Results of the initial segmentation, using 'segmentByNnn', where 'Nnn' corresponds to argument \code{flavor}}
\item{dpBkpList}{Results of dynamic programming, a list of vectors of breakpoint positions for the best model with k breakpoints for k=1, 2, ... K where \code{K=length(initBkp)}}
\item{prof}{a \code{matrix} providing time usage (in seconds) and memory usage (in Mb) for the main steps of the program.  Only defined if argument \code{profile} is set to \code{TRUE}
}}
\references{Bleakley, K., & Vert, J. P. (2011). The group fused
lasso for multiple change-point detection. arXiv preprint
arXiv:1106.4199.

Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
change-points shared by many signals using group LARS. Advances in
Neural Information Processing Systems, 23, 2343-2351.

Chen, H., Xing, H., & Zhang, N. R. (2011). Estimation
of parent specific DNA copy number in tumors using high-density
genotyping arrays. PLoS computational biology, 7(1), e1001060.

Rigaill, G. (2010). Pruned dynamic programming for
optimal multiple change-point detection. arXiv preprint
arXiv:1004.0887.

Bellman, R. (1956). Dynamic programming and Lagrange multipliers.
Proceedings of the National Academy of Sciences of the United States
of America, 42(10), 767.

Mosen-Ansorena, D. and Aransay, A.-M.(2013)
Bivariate segmentation of SNP-array data for allele-specific
copy number analysis in tumour samples
BMC Bioinformatics, 14:84 

Killick R, Fearnhead P, Eckley IA (2012)
Optimal detection of changepoints with
a linear computational cost, JASA 107(500),1590-1598  }
\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{segmentByRBS}}

\code{\link{pruneByDP}}

\code{\link{segmentByGFLars}}

\code{\link{segmentByPSCN}}

\code{\link{segmentByPSCBS}}

\code{\link{segmentByCBS}}

\code{\link{segmentByPelt}}

\code{\link{segmentByCnaStruct}}}
\examples{
p <- 2
trueK <- 10
len <- 1e4
sim <- randomProfile(len, trueK, 1, p)
Y <- sim$profile
K <- 2*trueK
res <- jointSeg(Y, K=K)
bkp <- res$bestBkp
getTpFp(bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives

par(mfrow=c(p,1))
for (ii in 1:p) {
  plot(Y[, ii], pch=19, cex=0.2)
  abline(v=bkp+0.5, col=2, lwd=2)
  abline(v=sim$bkp+0.5, col=8, lty=2, lwd=3)
}

## Now we add some NA:s in one dimension
jj <- sim$bkp[1]
Y[jj-seq(-10,10), p] <- NA
res2 <- jointSeg(Y, K=K, verbose=TRUE)
bkp <- res2$bestBkp
getTpFp(res2$bestBkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
}
