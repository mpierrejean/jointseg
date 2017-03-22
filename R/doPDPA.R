#' High-level function for segmentation by pruned dynamic programming algorithm (pDPA)
#' 
#' @param Y A numeric vector the signal to be segmented
#' @param K The number of change points to find
#' @param \dots Further parameters to be passed to Fpsn
#' @author Morgane Pierre-Jean, Guillem Rigaill and Pierre Neuvial
#' @references Rigaill, G. (2015). A pruned dynamic programming algorithm to recover the best segmentations with 1 to K_max change-points. Journal de la Société Française de Statistique, 156(4), 180-205.
#' @examples
#' ## load known real copy number regions
#' affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=1)
#' 
#' ## generate a synthetic CN profile
#' K <- 4
#' len <- 1e3
#' sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
#' datS <- sim$profile
#' 
#' ## run pruned DPA segmentation
#' res <- doPDPA(datS[["c"]], K=2*K)
#' getTpFp(res$bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
#' plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
#' @export
doPDPA <- function(Y, K, ...) {
    res <- Fpsn(Y, K+1, ...)
    
    ## convert matrix of breakpoints to list of breakpoints
    bkpMat <- res$t.est
    bkpList <- lapply(1:K+1, FUN=function(kk) {
        bkpMat[kk, 1:(kk-1)]
    })
    
    rse <- res$J.est ## residual squared error
    V <- res$allCost ## cost matrix
    
    list(bkpList=bkpList, rse=rse, V=V)
}
