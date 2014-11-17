PSSeg <- structure(function(#Parent-Specific copy number segmentation
### This function splits (bivariate) copy number signals into
### parent-specific (PS) segments using recursive binary segmentation
                            data,
### Data frame containing the following columns: \describe{
### \item{c:}{Total copy number (logged or non-logged)}
### \item{b:}{Allele B fraction}
### \item{genotype:}{(germline) genotype of the SNP, coded as 0 for AA, 1/2 for AB, 1 for BB}
### }
### These data are assumed to be ordered by genome position.
                            method,## A \code{character} value, the type of segmentation method used. May be one of 
### \describe{
###   \item{"RBS"}{Recursive Binary Segmentation, see
### \code{\link{doRBS}}}
###   \item{"GFLars"}{Group fused LARS as described in Bleakley and
###   Vert (2011).}
###   \item{"DP"}{Univariate pruned dynamic programming Rigaill et al (2010) or bivariate dynamic programming}
###   \item{"PSCBS"}{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation by Olshen A. et al
###     (2011)}
###   \item{"other"}{The segmentation method is passed as a function using argument \code{segFUN} (see examples in directory \code{otherMethods}).}
###}
                            stat=NULL,
### A vector containing the names or indices of the columns of \code{Y} to be segmented
                            dropOutliers=TRUE,
### If TRUE, outliers are droped by using DNAcopy package
                            ...,
### Further arguments to be passed to \code{jointSeg}
                            profile=FALSE,
### Trace time and memory usage ?
                            verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                            ){
    ## Argument
    cn <- colnames(data)
    ecn <- c("c", "b", "genotype") ## expected
    mm <- match(ecn, cn)
    if (any(is.na(mm))) {
        str <- sprintf("('%s')", paste(ecn, collapse="','"))
        stop("Argument 'data' should contain columns named ", str)
    }

    prof <- NULL

    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Pre-processing
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    n <- nrow(data)
    data$idx <- 1:n
    data$d <- 2*abs(data$b-1/2)
    isHet <- (data[["genotype"]]==0.5)
    data[which(!isHet), "d"] <- NA  
    ##details<<Before segmentation, the decrease in heterozygosity
    ##\code{d=2|b-1/2|} defined in Bengtsson et al, 2010 is calculated
    ##from the input data.  \code{d} is only defined for heterozygous
    ##SNPs, that is, SNPs for which \code{data$genotype==1/2}. \code{d}
    ##may be seen as a "mirrored" version of allelic ratios (\code{b}):
    ##it converts them to a piecewise-constant signals by taking
    ##advantage of the bimodality of \code{b} for heterozygous SNPs.
    ##The rationale for this transformation is that allelic ratios
    ##(\code{b}) are only informative for heterozygous SNPs (see
    ##e.g. Staaf et al, 2008).

    if(dropOutliers){
        CNA.object <- CNA(data$c, rep(1, n), 1:n)
        smoothed.CNA.obj <- smooth.CNA(CNA.object)
        data$c <- smoothed.CNA.obj$Sample.1
    }
    ##details<<Before segmentation, the outliers in the copy number signal are
    ##droped according the method explained by
    ##Venkatraman, E. S. and Olshen, A. B., 2007.

    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ## Segmentation followed by pruning using dynamic programming
    ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ##details<<The resulting data are then segmented using the
    ##\code{\link{jointSeg}} function, which combines an initial
    ##segmentation according to argument \code{method} and pruning of
    ##candidate change points by dynamic programming (skipped when the
    ##initial segmentation *is* dynamic programming).
    if (is.null(stat)) {
        ##details<<If argument \code{stat} is not provided, then dynamic
        ##programming is run on the two dimensional statistic
        ##\code{"(c,d)"}.
        mm <- match(c("c", "d"), colnames(data))
        res <- jointSeg(data, method=method, dpStat=mm, ..., profile=profile, verbose=verbose)
    } else {
        ##details<<If argument \code{stat} is provided, then dynamic
        ##programming is run on \code{stat}; in this case we implicitly
        ##assume that \code{stat} is a piecewise-constant signal.
        res <- jointSeg(data, method=method, stat=stat, dpStat=stat, ..., profile=profile, verbose=verbose)    
    }
    prof <- rbind(prof, res$prof)

    

    ##references<<Bengtsson, H., Neuvial, P., & Speed,
    ##T. P. (2010). TumorBoost: Normalization of allele-specific tumor
    ##copy numbers from a single pair of tumor-normal genotyping
    ##microarrays. BMC bioinformatics, 11(1), 245.

    ##references<<Staaf, J., Lindgren, D., Vallon-Christersson, J.,
    ##Isaksson, A., Goransson, H., Juliusson, G., ... & Ringn\'er,
    ##M. (2008). Segmentation-based detection of allelic imbalance and
    ##loss-of-heterozygosity in cancer cells using whole genome SNP
    ##arrays. Genome Biol, 9(9), R136.

    ##seealso<<\code{\link{jointSeg}}
    list(
        bestBkp=res$bestBkp, ##<< Best set of breakpoints after
        ##dynamic programming
        initBkp=res$initBkp, ##<< Results of the initial segmentation,
        ##using 'doNnn', where 'Nnn' corresponds
        ##to argument \code{method}
        dpBkpList=res$dpBkpList, ##<< Results of dynamic programming,
        ##a list of vectors of breakpoint
        ##positions for the best model with k
        ##breakpoints for k=1, 2, ... K where
        ##\code{K=length(initBkp)}
        prof=prof ##<< a \code{matrix} providing time usage (in
        ##seconds) and memory usage (in Mb) for the main
        ##steps of the program.  Only defined if argument
        ##\code{profile} is set to \code{TRUE}
        )
}, ex=function(){
    ## load known real copy number regions
    affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=0.5)

    ## generate a synthetic CN profile
    K <- 10
    len <- 1e4
    sim <- getCopyNumberDataByResampling(len, K, regData=affyDat)
    datS <- sim$profile

    ## run binary segmentation (+ dynamic programming) resRBS <-
    resRBS <- PSSeg(data=datS, method="RBS", stat=c("c", "d"), K=2*K, profile=TRUE)
    resRBS$prof

    getTpFp(resRBS$bestBkp, sim$bkp, tol=5)
    plotSeg(datS, breakpoints=list(sim$bkp, resRBS$bestBkp))
})
############################################################################
## HISTORY:
## 2014-05-20
## o Argument 'flavor' renamed to 'method'.
## o Argument 'statistic' renamed to 'stat'.
## o Removed arguments 'jitter', 'DP' and 'methModelSelection'. These
## arguments may still be used through '...'.
## 2014-05-06
## Changed the mapping: when all probes are not used : final
## breakpoints are the median between two successive used probes and
## not the used probes before the breakpoint.
## 2013-12-09
## Replaced flavor 'cghseg' by 'DP'. 
## 2013-12-06
## Removed 'log(c)' statistic (left up to the user).
## For flavors 'CnaStruct' and 'PSCN', force conversion of total CNs to log
## scale if input data are not logged.
## 2013-11-29
## Added flavor : 'DP'.
## Cleanups in default arguments.
## 2013-03-28
## Added flavors : 'CnaStruct' and 'Pelt'.
## 2013-03-07
## Added option 'DP' for flavor "RBS" to do selection on initial segmentation.
## 2013-02-27
## o Bug fixed : flavor "GFLars" could not be run at full resolution.
## o Added statistic 'd|het'.
## 2013-02-26
## o Added option 'jitter' to allow more precise breakpoint identification by DP.
## 2013-02-15
## o Added flavor 'PSCBS'.
## 2013-01-28
## o Bug fix returned bkp from dp.
## 2013-01-25
## o Cleanups in doc and return values.
## 2013-01-09
## o Replaced all jumps by bkp.
## 2013-01-03
## o Added argument 'flavor' (passed to 'jointSeg').
## o Added flavors 'PSCN' and 'cghseg'.
## 2012-12-31
## o Added Argument 'statistic'.
## o Each dimension is now scaled to unit variance using 'estimateSd'.
## 2012-12-27
## o Some code and doc cleanups.
## 2012-12-15
## o Added argument 'profile' for optional reporting of CPU and memory usage.
## 2012-12-03
## o Added smoothing function from matrixStats.
## 2012-11-30
## o Created.
############################################################################

