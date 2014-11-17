jointSeg <- structure(function(# Joint segmentation of multivariate signals
### Joint segmentation of multivariate signals in two steps:
### \enumerate{
###   \item{first-pass segmentation.  By default, a fast, greedy
###     approach is used (see \code{method}).}
###   \item{pruning of the candidate change points obtained by
###     dynamic programming}
### }
                               Y,
### The signal to be segmented (a matrix or a numeric vector)
                               method="RBS",
### A \code{character} value, the type of segmentation method used.  May be one of:
### \describe{
###   \item{"RBS"}{Recursive Binary Segmentation (the default), see
### \code{\link{segmentByRBS}} as described in Gey and Lebarbier (2005)}
###   \item{"GFLars"}{Group fused LARS as described in Bleakley and
###   Vert (2011).}
###   \item{"DP"}{Dynamic Programming (Bellman, 1956). For univariate signals the pruned DP of  Rigaill et al (2010) is used.}
###   \item{"other"}{The segmentation method is passed as a function using argument \code{segFUN} (see examples in directory \code{otherMethods} of the \code{jointseg} package).}
###}
                               stat=NULL,
### A vector containing the names or indices of the columns of \code{Y} to be segmented
                               dpStat=stat,
### A vector containing the names or indices of the columns of \code{Y} to be segmented at the second round of segmentation. Defaults to the value of argument \code{stat}.
                               segFUN=NULL,
### The segmentation function to be used when \code{method} is set to \code{other}. Not used otherwise.
                               jitter=NULL,
### Uncertainty on breakpoint position after initial segmentation.  Defaults to \code{NULL}.  See Details.
                               modelSelectionMethod=ifelse(method %in% c("DynamicProgramming", "RBS", "GFLars"), "Lebarbier", "none"),
### Which method is used to perform model selection.
                               modelSelectionOnDP=(method %in% c("DynamicProgramming", "RBS", "GFLars")),
### If \code{TRUE} (the default), model selection is performed on
### segmentation after dynamic programming; else model selection is
### performed on initial segmentation.  Only applies to methods "DP",
### "RBS" and "GFLars".
                               ...,
### Further arguments to be passed to the lower-level segmentation
### method determined by argument \code{method}.
                               profile=FALSE,
### Trace time and memory usage ?
                               verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                               ){
    ##references<<Bleakley, K., & Vert, J. P. (2011). The group fused
    ##lasso for multiple change-point detection. arXiv preprint
    ##arXiv:1106.4199.

    ##references<<Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
    ##change-points shared by many signals using group LARS. Advances in
    ##Neural Information Processing Systems, 23, 2343-2351.

    ##references<<Gey, S., & Lebarbier, E. (2008). Using CART to Detect
    ##Multiple Change Points in the Mean for Large
    ##Sample. http://hal.archives-ouvertes.fr/hal-00327146/

    ##references<<Rigaill, G. (2010). Pruned dynamic programming for
    ##optimal multiple change-point detection. arXiv preprint
    ##arXiv:1004.0887.

    ##references<<Bellman, R. (1956). Dynamic programming and Lagrange multipliers.
    ##Proceedings of the National Academy of Sciences of the United States
    ##of America, 42(10), 767.

    ##seealso<<\code{\link{pruneByDP}}

    if (is.null(dim(Y))) {
        ## coerce to a matrix
        Y <- as.matrix(Y)
    }
    
    ## drop row names
    rownames(Y) <- NULL
    
    ## Arguments 'stat' and 'dpStat'
    if (is.null(stat)) {
        stat <- 1:ncol(Y)
        if (is.null(dpStat)) {
            dpStat <- stat
        }
    }
    if (mode(stat) != mode(dpStat)) {
        stop("Arguments 'stat' and 'dpStat' should be of the same mode ('numeric' or 'character')")
    }
    arg <- union(stat, dpStat)
    if (!is.null(arg)) {
        if (is.numeric(arg)) {
            mm <- match(arg, 1:ncol(Y)) 
        } else if (is.character(arg)) {
            mm <- match(arg, colnames(Y))
        } 
        if (sum(is.na(mm))) {
            guilty <- paste("'", arg[which(is.na(mm))], "'", sep="", collapse=",")
            stop("Undefined column(s) selected in 'Y':", guilty, ". Please check arguments 'stat' and 'dpStat'")
        } 
    }

    ## Case of rows with all entries missing (typically occurs when 'stat' is 'd')
    n <- nrow(Y)
    allNA <- apply(Y[, mm, drop=FALSE], 1, FUN=function(x) all(is.na(x)))
    ww <- which(!allNA)

    outPos <- 1:n  ## in order to be able to map back to 'position' indices from the input data
    if (length(ww)<n) { 
        outPos <- mapPositionsBack(outPos[ww])
    }
    
    Yseg <- Y[ww, stat]        ## for the initial segmentation
    nSeg <- ifelse(!is.null(dim(Yseg)), nrow(Yseg), length(Yseg))
    Ydp <- Y[ww, dpStat]       ## for pruning by DP
    nDp <- ifelse(!is.null(dim(Ydp)), nrow(Ydp), length(Ydp))
    Y <- NULL; rm(Y);                      ## not needed anymore
    
    ## Segmentation function
    if (method=="other") {
        if (is.null(segFUN) || mode(segFUN)!="function") {
            stop("A segmentation function should be provided using argument 'segFUN' for method 'other'")
        }
        nms <- names(formals(segFUN))
        inter <- intersect(nms, c("jitter", "methModelSelection", "DP"))
        if (length(inter)) {
            warning("Argument(s) ", paste("'", inter, "'", collapse=",", sep=""), " will not be passed to argument 'segFUN' as they match existing arguments of the 'jointSeg' function")
        }
    } else {
        if (!is.null(segFUN)) {
            warning("Argument 'segFUN' is only used when 'method' is set to 'other'")
        }
        segName <- paste("do", method, sep="")

        ## retrieve segmentation function and assert that it does exist
        if (!exists(segName, mode="function")) {
            stop("Function '", segName, "' should be provided for method '", method, "'")
        } else {
            segFUN <- eval(as.name(segName))    
        }
    }
    
    if (verbose) {
        cat("Start initial segmentation by", method, "\n")
        cat("Statistic:", stat, "\n")
        str(Yseg)
    }
    prof <- NULL
    resSeg <- prof(segFUN(Yseg, ...), doit=profile)
    initSeg <- resSeg$res
    prof <- rbind(prof, segmentation=resSeg$prof)
    if (verbose) {
        cat("End initial segmentation by", method, "\n")
        if (profile) {
            print(resSeg$prof)
        }
    }

    dpseg <- initSeg$dpseg
    if (!is.null(dpseg)) {
        ##details<<If the return value of the initial segmentation has an
        ##element named \code{dpseg}, then initial segmentation results
        ##are not pruned by dynamic programming.
    } else {
        ## Prune candidate breakpoints
        if (verbose) {
            cat("Start pruning by dynamic programming\n")
            cat("Statistic:", dpStat, "\n")
            str(Ydp)
        }
        bkp <- initSeg$bkp
        ##details<<If \code{jitter} is not \code{NULL}, it should be a
        ##vector of integer indices. The set of candidate breakpoints
        ##passed on to dynamic programming is augmented by all indices
        ##distant from an element of \code{jitter} from one of the
        ##candidates. For example, if \code{jitter==c(-1, 0, 1)} and the
        ##initial set of breakpoints is \code{c(1,5)} then dynamic
        ##programming is run on \code{c(1,2,4,5,6)}.
        if (!is.null(jitter)) {
            jitter <- as.integer(jitter)
            bkpJ <- sapply(bkp, FUN=function(x) {
                x+jitter
            })
            bkpJ <- unique(bkpJ)  ## remove duplicates
            bkpJ <- bkpJ[bkpJ>=1]
            bkpJ <- bkpJ[bkpJ<n]
            bkp <- bkpJ
        }

        resDP <- prof(pruneByDP(Ydp, candCP=bkp, verbose=verbose), doit=profile)
        dpseg <- resDP$res
        prof <- rbind(prof, dpseg=resDP$prof)
        if (verbose) {
            cat("End pruning by dynamic programming\n")
            if (profile) {
                print(resDP$prof)
            }
        }
    }

    ## Find the best segmentation
    if (modelSelectionMethod == "none") {
        bestSeg <- initSeg$bkp
    } else {
        if (modelSelectionOnDP) {         ## run model selection on results of dynamic programming
            mS <- modelSelection(dpseg$rse, n=nDp, method=modelSelectionMethod)
            if (verbose) {
                str(mS)
            }
            bestSeg <- integer(0L)
            if (mS$kbest!=0) {
                bestSeg <- dpseg$bkp[[mS$kbest]]
            } else {                        ## run model selection on initial segmentation
                ##details<<If
                ##\code{modelSelectionOnDP} is set to \code{FALSE}, then model
                ##selection is run on the sets of the form \code{bkp[1:k]} for
                ##\eqn{1 \leq k \leq length(bkp)}, where \code{bkp} is the set of
                ##breakpoints identified by the initial segmentation.  In
                ##particular, this implies that the candidate breakpoints in
                ##\code{bkp} are sorted by order of appearance and not by
                ##position.
                mS <- modelSelection(initSeg$rse, n=nSeg, method=modelSelectionMethod)
                bestSeg <- integer(0L)
                if (mS$kbest!=0) {
                    bestSeg <- sort(initSeg$bkp[1:mS$kbest])
                }
            }
        }
    }

    ## map breakpoint positions back to original space (if required)
    bestBkp <- outPos[bestSeg]
    initBkp <- outPos[initSeg$bkp]
    dpBkpList <- lapply(dpseg$bkp,function(bkp) outPos[bkp])
    ##value<< A list with elements:
    list(
        bestBkp=bestBkp, ##<< Best set of breakpoints after dynamic programming
        initBkp=initBkp, ##<< Results of the initial segmentation, using
        ##'doNnn', where 'Nnn' corresponds to argument
        ##\code{method}
        dpBkpList=dpBkpList, ##<< Results of dynamic programming, a list
        ##of vectors of breakpoint positions for
        ##the best model with k breakpoints for
        ##k=1, 2, ... K where
        ##\code{K=length(initBkp)}
        prof=prof) ##<< a \code{matrix} providing time usage (in
    ##seconds) and memory usage (in Mb) for the main steps
    ##of the program. Only defined if argument
    ##\code{profile} is set to \code{TRUE}
}, ex=function(){
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
})

############################################################################
## HISTORY:
## 2014-07-16
## o It is now possible to run 'modelSelection' on initial
## segmentation.  Default behavior of the function is unchanged.
## 2014-05-14
## o Argument 'flavor' renamed to 'method'.
## o Added argument 'segFUN' and method "other" to enable a
## user-defined segmentation method to be used.
## o Removed default argument values for 'method', in order to make
## the choice of a segmentation method explicit for the user.
## o Moved all copy-number-specific code to downstream methods (e.g. 'doCBS').
## 2013-12-09
## o Replaced 'cghseg' by 'DP'
## o Renamed all 'segmentByNnn' to 'doNnn'
## o Force Ydp to be a matrix
## 2013-12-05
## o Now dropping row names of 'Yseg' and 'Ydp'.
## 2013-11-29
## Added flavor : 'DP'.
## 2013-03-28
## Added flavors : 'CnaStruct' and 'Pelt'
## 2013-03-07
## Added option 'DP' for flavor "RBS" to do selection on initial segmentation
## 2013-02-26
## o Added option 'jitter' to allow more precise breakpoint identification by DP.
## 2013-02-18
## o Add flavor 'PSCBS'
## 2013-01-25
## o Cleanups in doc and return values.
## o Now returning 'bestBkp'. 
## o Removed 'position' from required fields in input data.
## 2013-01-09
## o Replaced all jumps by bkp.
## 2013-01-03
## o Added flavors 'PSCN' and 'cghseg'.
## 2012-12-27
## o Default flavor is now recursive binary segmentation.
## o Some code and doc cleanups.
## 2012-12-15
## o Added argument 'profile' for optional reporting of CPU and memory
## usage.
## 2012-12-06
## o Added argument flavor.
## 2012-12-04
## o Created.
############################################################################

