#' profile time and memory usage of a given R expression
#' 
#' profile time and memory usage of a given R expression
#' 
#' Profiling is performed using \code{summaryRprof(memory="both")$by.self}.
#' 
#' @note Our memory profiling is not satistfactory yet!
#' 
#' @param expr An \code{R} expression to be evaluated
#' @param doit A boolean variable specifying whether profiling should be
#' performed or not (intended for internal use).
#' @author Morgane Pierre-Jean and Pierre Neuvial
#' @seealso \code{\link{Rprof}}
#' 
#' \code{\link{summaryRprof}}
#' @keywords internal
prof <- function(expr, doit=TRUE) {
    prf <- NULL
    if (!doit) {
        res <- eval(expr)
    } else {
        tf <- tempfile()
        utils::Rprof(tf, memory.profiling=TRUE)
        res <- eval(expr)
        utils::Rprof(NULL)
        
        ## check that something has been reported
        rl <- readLines(tf)
        if (length(rl)>3) {
            sp <- utils::summaryRprof(tf, memory="both")
            file.remove(tf)
            prf <- colSums(sp$by.self[, c("self.time", "mem.total")])
            names(prf) <- c("time", "memory")
        } else {
            prf <- c(0, NA)
        }
        names(prf) <- c("time", "memory")
    }
    list(res=res, prof=prf)
}
############################################################################
## HISTORY:
## 2013-01-30
## o Now an internal function (not exported anymore).
## 2012-12-27
## o Some code and doc cleanups.
## 2012-12-15
## o Created.
############################################################################

