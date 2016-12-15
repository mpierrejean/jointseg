#' Run PSCN segmentation (defunct)
#' 
#' The 'PSCN' package is not maintained anymore and it is not available for R
#' >= 3.0.0.  The original 'doPSCN' function has been moved to the directory
#' 'zzz.defunct'.  The skeleton of that function is kept for backward
#' compatibility.
#' 
#' 
#' @param Y The signal to be segmented, a matrix containing the following
#' columns: \describe{ \item{c}{Total copy number (log scale)} \item{b}{Allele
#' B fraction (a.k.a. BAF)} }
#' @param alpha sensitivity level in [0,1] to be passed to
#' \code{PSCN::segmentation}.
#' @param platform Specifies form which array platform 'Y' was generated:
#' Illumina or Affymetrix
#' @param verbose A \code{logical} value: should extra information be output ?
#' Defaults to \code{FALSE}.
#' @param \dots Further arguments to be passed to \code{PSCN::smoothing}
#' @author Morgane Pierre-Jean and Pierre Neuvial
#' @seealso \code{\link{PSSeg}}
#' @references Chen, H., Xing, H., & Zhang, N. R. (2011). Estimation of parent
#' specific DNA copy number in tumors using high-density genotyping arrays.
#' PLoS computational biology, 7(1), e1001060.
#' @examples
#' 
#' print("The 'PSCN' package is not available for R >= 3.0.0.")
#' print("See http://cran.r-project.org/web/packages/PSCN/index.html")
#' 
#' @export doPSCN
doPSCN <- function(Y, alpha=0.01, platform=c("Illumina", "Affymetrix"), verbose=FALSE, ...) {
    stop("The 'PSCN' package is not available for R >= 3.0.0.\nSee http://cran.r-project.org/web/packages/PSCN/index.html")
}
############################################################################
## HISTORY:
## 2014-02-13
## o Created from zzz.defunct/doPSCN.R.
############################################################################

