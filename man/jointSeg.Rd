\name{jointSeg}
\alias{jointSeg}
\title{Joint segmentation of multivariate signals}
\description{Joint segmentation of multivariate signals in two steps:
\enumerate{
\item{first-pass segmentation.  By default, a fast, greedy
approach is used (see \code{method}).}
\item{pruning of the candidate change points obtained by
dynamic programming}
}}
\usage{jointSeg(Y, method, stat = 1:ncol(Y), dpStat = stat, segFUN = NULL, 
    jitter = NULL, modelSelectionMethod = "Birge", modelSelectionOnDP = TRUE, 
    ..., profile = FALSE, verbose = FALSE)}
\arguments{
  \item{Y}{The signal to be segmented (must be a matrix)}
  \item{method}{A \code{character} value, the type of segmentation method used.  May be one of:
\describe{
\item{"RBS"}{Recursive Binary Segmentation (the default), see
\code{\link{doRBS}} as described in Gey and Lebarbier (2005)}
\item{"GFLars"}{Group fused LARS as described in Bleakley and
Vert (2011).}
\item{"DP"}{Dynamic Programming (Bellman, 1956). For univariate signals the pruned DP of  Rigaill et al (2010) is used.}
\item{"other"}{The segmentation method is passed as a function using argument \code{segFUN}}
}}
  \item{stat}{A vector containing the names or indices of the columns of \code{Y} to be segmented}
  \item{dpStat}{A vector containing the names or indices of the columns of \code{Y} to be segmented at the second round of segmentation. Defaults to the value of argument \code{stat}.}
  \item{segFUN}{The segmentation function to be used when \code{method} is set to \code{other}. Not used otherwise.}
  \item{jitter}{Uncertainty on breakpoint position after initial segmentation.  Defaults to \code{NULL}.  See Details.}
  \item{modelSelectionMethod}{Which method is used to perform model selection.}
  \item{modelSelectionOnDP}{If \code{TRUE} (the default), model selection is performed on
segmentation after dynamic programming; else model selection is
performed on initial segmentation.  Only applies to methods "DP",
"RBS" and "GFLars".}
  \item{\dots}{Further arguments to be passed to the lower-level segmentation
method determined by argument \code{method}.}
  \item{profile}{Trace time and memory usage ?}
  \item{verbose}{A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.}
}
\details{If the return value of the initial segmentation has an
element named \code{dpseg}, then initial segmentation results
are not pruned by dynamic programming.

If \code{jitter} is not \code{NULL}, it should be a
vector of integer indices. The set of candidate breakpoints
passed on to dynamic programming is augmented by all indices
distant from an element of \code{jitter} from one of the
candidates. For example, if \code{jitter==c(-1, 0, 1)} and the
initial set of breakpoints is \code{c(1,5)} then dynamic
programming is run on \code{c(1,2,4,5,6)}.

For methods "DP", "RBS", "GFLars", if
\code{modelSelectionOnDP} is set to \code{FALSE}, then model
selection is run on the sets of the form \code{bkp[1:k]} for
\eqn{1 \leq k \leq length(bkp)}, where \code{bkp} is the set of
breakpoints identified by the initial segmentation.  In
particular, this implies that the candidate breakpoints in
\code{bkp} are sorted by order of appearance and not by
position.}
\value{A list with elements:
\item{bestBkp}{Best set of breakpoints after dynamic programming}
\item{initBkp}{Results of the initial segmentation, using
'doNnn', where 'Nnn' corresponds to argument
\code{method}}
\item{dpBkpList}{Results of dynamic programming, a list
of vectors of breakpoint positions for
the best model with k breakpoints for
k=1, 2, ... K where
\code{K=length(initBkp)}}
\item{prof}{a \code{matrix} providing time usage (in
seconds) and memory usage (in Mb) for the main steps
of the program. Only defined if argument
\code{profile} is set to \code{TRUE}
}}
\references{Bleakley, K., & Vert, J. P. (2011). The group fused
lasso for multiple change-point detection. arXiv preprint
arXiv:1106.4199.

Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
change-points shared by many signals using group LARS. Advances in
Neural Information Processing Systems, 23, 2343-2351.

Gey, S., & Lebarbier, E. (2008). Using CART to Detect
Multiple Change Points in the Mean for Large
Sample. http://hal.archives-ouvertes.fr/hal-00327146/

Rigaill, G. (2010). Pruned dynamic programming for
optimal multiple change-point detection. arXiv preprint
arXiv:1004.0887.

Bellman, R. (1956). Dynamic programming and Lagrange multipliers.
Proceedings of the National Academy of Sciences of the United States
of America, 42(10), 767.}
\author{Morgane Pierre-Jean and Pierre Neuvial}



\seealso{\code{\link{pruneByDP}}}
\examples{
## A two-dimensional signal
p <- 2
trueK <- 10
len <- 1e4
sim <- randomProfile(len, trueK, 1, p)
Y <- sim$profile
K <- 2*trueK
res <- jointSeg(Y, method="RBS", K=K)
bkp <- res$bestBkp
getTpFp(bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives

plotSeg(Y, list(sim$bkp, res$bestBkp), col=1)
      
## Now we add some NA:s in one dimension
jj <- sim$bkp[1]
Y[jj-seq(-10,10), p] <- NA
res2 <- jointSeg(Y, method="RBS", K=K, verbose=TRUE)
bkp <- res2$bestBkp
getTpFp(res2$bestBkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
}
