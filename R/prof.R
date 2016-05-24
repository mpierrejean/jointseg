prof <- function(# profile time and memory usage of a given R expression
### profile time and memory usage of a given R expression
                 expr,
### An \code{R} expression to be evaluated
                 doit=TRUE
### A boolean variable specifying whether profiling should be
### performed or not (intended for internal use).
                 ) {
    ##keyword<<internal

    ##details<<Profiling is performed using \code{summaryRprof(memory="both")$by.self}.
    ##details<<Memory profiling is not satistfactory yet.
    ##seealso<<\code{\link{Rprof}}
    ##seealso<<\code{\link{summaryRprof}}
    prof <- NULL
    if (!doit) {
        res <- eval(expr)
    } else {
        tf <- tempfile()
        Rprof(tf, memory.profiling=TRUE)
        res <- eval(expr)
        Rprof(NULL)
        
        ## check that something has been reported
        rl <- readLines(tf)
        if (length(rl)>3) {
            sp <- summaryRprof(tf, memory="both")
            file.remove(tf)
            prof <- colSums(sp$by.self[, c("self.time", "mem.total")])
            names(prof) <- c("time", "memory")
        } else {
            prof <- c(0, NA)
        }
        names(prof) <- c("time", "memory")
    }
    list(res=res, prof=prof)
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

