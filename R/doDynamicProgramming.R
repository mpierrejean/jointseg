#' Run segmentation by dynamic programming
#'
#' High-level function for univariate or multivariate segmentation by dynamic
#' programming
#'
#' If the signal is uni-dimensional, this function simply uses the segmentation
#' method provided in the \code{cghseg} package reshapes the results.
#'
#' If the signal is multi-dimensional, this function applies the
#' \code{\link{pruneByDP}} function and reshapes the results.
#'
#' @param Y A numeric vector or a matrix, the signal to be segmented
#' @param K The number of change points to find
#' @param stat A vector containing the names or indices of the columns of
#'   \code{Y} to be segmented
#' @param verbose A \code{logical} value: should extra information be output ?
#'   Defaults to \code{FALSE}.
#' @return \item{bkp}{A vector of \code{K} indices for candidate change points}
#'   \item{dpseg}{A list of two elements \describe{ \item{bkp}{A list of vectors
#'   of change point positions for the best model with k change points, for k=1,
#'   2, ... K} \item{rse}{A vector of K+1 residual squared errors} }}
#' @note This is essentially a wrapper for convenient segmentation by dynamic
#'   programming using the \code{\link{PSSeg}} function.
#' @author Morgane Pierre-Jean and Pierre Neuvial
#' @encoding utf-8
#' @references Rigaill, G. (2015). A pruned dynamic programming algorithm to recover the best segmentations with 1 to K_max change-points. Journal de la Societe Francaise de Statistique, 156(4), 180-205.
#' @examples
#'
#' ## load known real copy number regions
#' affyDat <- acnr::loadCnRegionData(dataSet="GSE29172", tumorFraction=1)
#'
#' ## generate a synthetic CN profile
#' K <- 10
#' len <- 1e4
#' sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
#' datS <- sim$profile
#'
#' ## run pruned DPA segmentation
#' resDP <- doDynamicProgramming(datS[["c"]], K=K)
#' getTpFp(resDP$bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
#' plotSeg(datS, breakpoints=list(sim$bkp, resDP$bkp))
#'
#' ## run 2d dynamic programming segmentation
#' K <- 2
#' len <- 1e3
#' sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
#' datS <- sim$profile
#' datS$d <- 2*abs(datS$b-1/2)
#' datS[which(datS$genotype!=0.5),"d"] <- NA
#' Y = cbind(datS$c,datS$d)
#' resDP2d <- doDynamicProgramming(Y, K = K)
#'
#' @export doDynamicProgramming
doDynamicProgramming <- function(Y, K, stat=NULL, verbose=FALSE){
    ## Argument 'Y'
    if (is.null(dim(Y)) || is.data.frame(Y)) {
        if (verbose) {
            print("Coercing 'Y' to a matrix")
        }
        Y <- as.matrix(Y)
    } else if (!is.matrix(Y)) {
        stop("Argument 'Y' should be a matrix, vector or data.frame")
    }

    ## Argument 'stat'
    if (!is.null(stat)) {
        if (is.numeric(stat)) {
            mm <- match(stat, 1:ncol(Y))
        } else if (is.character(stat)) {
            mm <- match(stat, colnames(Y))
        }
        if (sum(is.na(mm))) {
            guilty <- paste("'", stat[which(is.na(mm))], "'", sep="", collapse=",")
            stop("Undefined column(s) selected in 'Y':", guilty, ". Please check argument 'stat'")
        } else {
            Y <- Y[, mm, drop=FALSE]
        }
    }

    if (is.null(dim(Y)) || (ncol(Y)==1)) {
        res <- Fpsn(Y, Kmax=K+1)
        ## convert matrix of breakpoints to list of breakpoints
        bkpMat <- res$t.est
        bkpList <- lapply(1:K+1, FUN=function(kk) {
            bkpMat[kk, 1:(kk-1)]
        })
        dpseg <- list(bkp=bkpList, rse=res$J.est, V=res$allCost)
        res <- list(bkp=bkpList[[K]], dpseg=dpseg)
    } else {
        res <- pruneByDP(Y, K=K+1)
        dpseg <- list(bkp=res$bkpList, rse=res$rse, V=res$V)
        res <- list(bkp=res$bkpList[[K]], dpseg=dpseg)
    }
    return(res)
}

############################################################################
## HISTORY:
## 2014-05-30
## o Added argument 'stat'.
## o Updated doc.
## 2013-12-09
## o Renamed to 'doDynamicProgramming'
## o Added 2d dynamic programming
## 2013-01-09
## o Replace all jumps by bkp
## 2013-01-04
## o BUG FIX: index shift when reshaping the results.
## 2013-01-03
## o Created.
############################################################################

