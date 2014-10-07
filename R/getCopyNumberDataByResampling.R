getCopyNumberDataByResampling <- structure(function(# Generate a copy number profile by resampling input data
### Generate a copy number profile by resampling input data
                                                    length,
### length of the profile
                                                    nBkp=NA,
### number of breakpoints.  If \code{NULL}, then argument \code{bkp} is expected to be provided. 
                                                    bkp=NULL,
### a numeric vector of breakpoint positions that may be used to
### bypass the breakpoint generation step.  Defaults to \code{NULL}.
                                                    regData=NULL,
### a data.frame containing copy number data for different types of copy number regions.  Columns:\describe{
### \item{c}{Total copy number}
### \item{b}{Allele B fraction (a.k.a. BAF)}
### \item{region}{a character value, annotation label for the
### region. See Details.}
### \item{genotype}{the (germline) genotype of SNPs. By definition, rows with
###  missing genotypes are interpreted as non-polymorphic loci (a.k.a. copy
###  number probes).}
###}
                                                    regions=NULL,
### a character vector of region labels that may be used to bypass the
### region label generation step.  Defaults to \code{NULL}.
                                                    regAnnot=NULL,
### a data.frame containing annotation data for each copy number region.  Columns:
### \describe{
### \item{region}{label of the form (must match \code{regData[["region"]]}).}
### \item{freq}{frequency (in [0,1]) of this type of region in the genome.}
### }
### If \code{NULL} (the default), frequencies of regions (0,1), (0,2), (1,1) and (1,2) (the most common alterations) are set to represent 90% of the regions.
### \code{sum(regAnnot[["freq"]])} should be 1.
                                                    minLength=0,
### minimum length of region between breakpoints.  Defaults to 0.
                                                    regionSize=0,
### If \code{regionSize>0}, breakpoints are included by pairs, where the
### distance within pair is set to \code{regionSize}.  \code{nBkp}
### is then required to be an even number.
                                                    connex=TRUE
### If \code{TRUE}, any two successive regions are constrained to be
### connex in the (minor CN, major CN) space.  See 'Details'.
                                                    ){
    ##details<<This function generates a random copy number profile of
    ## length 'length', with 'nBkp' breakpoints randomly chosen. Between two
    ## breakpoints, the profile is constant and taken among the different
    ## types of regions in \code{regData}.

    ## Sanity check: lengths
    lb <- length(bkp)
    if ((lb==0) & (is.na(nBkp))) {
        stop("Please provide one of 'nBkp' and 'bkp'")
    }
    ls <- length(regions)
    if (ls>0) {  ## 'regions' is not NULL
        if (lb>0) { ## 'bkp' is not NULL
            if (lb+1!=ls) {
                stop("length(regions) should be equal to length(bkp)+1")
            }
        } else { ## 'bkp' is NULL
            if (nBkp+1!=ls) {
                stop("length(regions) should be equal to length(bkp)+1")
            }
        }
    }
    
    ## Sanity check: consistent region names
    regNames <- unique(regData[["region"]])
    if (!is.null(regions)) {
        mm <- match(regions, regNames)
        if (any(is.na(mm))) {
            warning("Regions in argument 'regions' not found in 'regData[[\"region\"]]'\nSetting argument 'regions' to NULL")
            regions <- NULL
        }
    }
    if (is.null(regAnnot)) {
        commonReg <- intersect(c("(0,1)", "(0,2)", "(1,1)", "(1,2)"), regNames)
        nCommonReg <- length(commonReg)
        ## check if common region are in regDat
        if (nCommonReg>0) {
            freqC <- rep(0.90/nCommonReg, nCommonReg)
            otherReg <- setdiff(regNames, commonReg)
            nOtherRef <- length(otherReg)
            freqO <- rep(0.10/nOtherRef, nOtherRef)
            regAnnot <- data.frame(region=c(commonReg, otherReg),
                                   freq=c(freqC, freqO),
                                   stringsAsFactors=FALSE)
        } else {
            regAnnot <- data.frame(region=regNames,
                                   freq=1/length(regNames),
                                   stringsAsFactors=FALSE)
        }
    } else {
        mm <- match(regNames, regAnnot[["region"]])
        if (any(is.na(mm))) {
            stop("Annotation not found for: ", regNames[which(is.na(mm))])
        }
        freq <- regAnnot[["freq"]]
        if (min(freq)<0 || max(freq)>1) {
            stop("Elements of regAnnot[[\"freq\"]] should be in [0,1]")
        }
        if (sum(freq)>1) {
            stop("Elements of regAnnot[[\"freq\"]] should sum up to 1")
        }
    }
    
    ## Sanity check: regions 
    if (!is.null(regions)) {
        mm <- match(regNames, regAnnot[["region"]])
    }
    ## Order regAnnot as regNames
    o <- order(regAnnot$region)
    regAnnot <- regAnnot[o,]
    ##details<<Elements of \code{regData[["region"]]} must be of the form
    ##\code{"(C1,C2)"}, where \code{C1} denotes the minor copy number
    ##and \code{C2} denotes the major copy number.  For example,
    ##\describe{
    ## \item{(1,1)}{Normal}
    ## \item{(0,1)}{Hemizygous deletion}
    ## \item{(0,0)}{Homozygous deletion}
    ## \item{(1,2)}{Single copy gain}
    ## \item{(0,2)}{Copy-neutral LOH}
    ## \item{(2,2)}{Balanced two-copy gain}
    ## \item{(1,3)}{Unbalanced two-copy gain}
    ## \item{(0,3)}{Single-copy gain with LOH}
    ##}
    pattern <- "\\(([0-9]),([0-9])\\)"
    regAnnot$C1 <- as.numeric(gsub(pattern, "\\1", regNames))
    regAnnot$C2 <- as.numeric(gsub(pattern, "\\2", regNames))
    ww <- which(is.na(regAnnot$C1) | is.na(regAnnot$C2))
    if (length(ww)==length(regNames)) {
        stop("Could not retrieve minor and major CNs from regData[[\"region\"]]\nSee ?getCopyNumberDataByResampling")
    } else if (length(ww)) {
        warning("Could not retrieve region for: ", regNames[ww], "\nDropping these regions'")
        regNames <- regNames[-ww]
        regAnnot <- regAnnot[-ww, ]
    }
    names(regNames) <- regNames

    regTypes <- lapply(regNames, FUN=function(reg) {
        which(regData[["region"]]==reg)
    })

    if (is.null(bkp)) {
        ## Choose the breakpoints
        interval <- 1:(length-1)
        u <- numeric(0)
        if (regionSize==0){
            while (length(u)<nBkp) {
                j <- sample(x=interval, size=1 , replace=FALSE)
                jend <- j
                u <- c(u, j)
                b.inf <- max(1, j-minLength)
                b.sup <- min(length, j+minLength)
                v <- b.inf:b.sup
                interval <- setdiff(interval,v)
            }
        } else {
            interval <-1:(length-regionSize-1)
            while (length(u)<nBkp) {
                j <- sample(x=interval, size=1, replace=FALSE)
                jend <- j+regionSize
                u <- c(u, j)
                if (length(u)<nBkp) {  ## add the other member of the pair
                    u <- c(u, jend)
                    b.inf <- max(1, j-minLength-regionSize)
                    b.sup <- min(length, jend+minLength)
                    v <- b.inf:b.sup
                    interval <- setdiff(interval, v)
                }
            }
        }
        ## u <- sample(length-1, nBkp, replace=FALSE);
        bkp <- sort(u);
    }
    bkpB <- c(0, bkp, length);

    ##details<<If 'connex' is set to TRUE (the default), transitions
    ##between copy number regions are constrained in such a way that for
    ##any breakpoint, one of the minor and the major copy number does
    ##not change.  Equivalently, this means that all breakpoints can be
    ##seen in both total copy numbers and allelic ratios.
    candidateRegions <- function(regName) {
        if (is.null(regName)) return(regAnnot);
        reg <- subset(regAnnot, region==regName)
        d1 <- regAnnot[, "C1"]-reg[, "C1"]
        d2 <- regAnnot[, "C2"]-reg[, "C2"]
        ## todo: make sure that all regions are connex...
        if (connex){
            ww <- which((d1 & !d2) | (!d1 & d2))
            if (!length(ww)) {
                stop("No candidate region found !")
            }
        } else {
            ww <- which(regAnnot$region!=regName)
        }
        regAnnot[ww, ]
    }
    
    ## Add the random piecewise linear profile
    idxs <- NULL
    regs <- NULL
    region <- NULL
    for (ii in 1:(length(bkp)+1)) {
        n <- bkpB[ii+1]-bkpB[ii]
        ##meanII <- matrix(rnorm(dim), length(idxsII), dim, byrow=TRUE)
        if (is.null(regions)) {
            candReg <- candidateRegions(region)
            probs <- candReg[["freq"]]/sum(candReg[["freq"]])
            region <- sample(x=candReg[["region"]], prob=probs, size=1)
        } else {
            region <- regions[ii]
        }
        idxsII <- sample(x=regTypes[[region]], size=n, replace=TRUE)
        idxs <- c(idxs, idxsII)
        regs <- c(regs, region)
    }
    profile <- regData[idxs, ]
    return(list(bkp=bkp, profile=profile, regions=regs))
###\item{profile}{the profile (a \code{length} by \code{2} data.frame
### containing the same fields as the input data \code{regData}.}
###\item{bkp}{a vector of bkp positions (the last row index before a
### breakpoint)}
###\item{regions}{a character vector of region labels}
}, ex=function() {
    affyDat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=1)
    sim <- getCopyNumberDataByResampling(1e4, 5, minLength=100, regData=affyDat)
    plotSeg(sim$profile, sim$bkp)

    ## another run with identical parameters
    bkp <- sim$bkp
    regions <- sim$regions
    sim2 <- getCopyNumberDataByResampling(1e4, regData=affyDat, bkp=bkp, regions=regions)
    plotSeg(sim2$profile, bkp)

    ## change tumor fraction but keep same "truth"
    affyDatC <- loadCnRegionData(dataSet="GSE29172", tumorFraction=0.5)
    simC <- getCopyNumberDataByResampling(1e4, regData=affyDatC, bkp=bkp, regions=regions)
    plotSeg(simC$profile, bkp)
    
    ## restrict to only normal, single copy gain, and copy-neutral LOH
    ## with the same bkp
    affyDatR <- subset(affyDat, region %in% c("(1,1)", "(0,2)", "(1,2)"))
    simR <- getCopyNumberDataByResampling(1e4, regData=affyDatR, bkp=bkp)
    plotSeg(simR$profile, bkp)

    ## Same 'truth', on another dataSet
    regions <- simR$regions
    illuDat <- loadCnRegionData(dataSet="GSE11976", tumorFraction=1)
    sim <- getCopyNumberDataByResampling(1e4, regData=illuDat, bkp=bkp, regions=regions)
    plotSeg(sim$profile, sim$bkp)
})

############################################################################
## HISTORY:
## 2014-07-16
## Now, if 'regAnnot' is NULL (the default), frequencies of regions
## (0,1), (0,2), (1,1) and (1,2) (the most common alterations) are set
## to represent 90% of the regions.
## 2013-02-27
## o Added parameter 'connex' that forces adjacent regions to be connex if connex = TRUE
## 2013-01-23
## o Removed field 'position' from output.
## 2013-01-22
## o Added more sanity checks.
## 2013-01-16
## o Added arguments 'bkp' and 'regions' to allow for bypassing the breakpoint generation step.
## o Added argument 'regAnnot', through which theoretical frequencies
## for each CN regions can be specified.
## o Made the constraints on CN transitions more generic.
## 2013-01-09
## o Replaced all 'jumps' by 'bkp'.
## 2012-12-20
## o Added new argument 'regionSize'.
## 2012-12-01
## o Added example data files and script based on public data set GSE19539.
## 2012-11-25
## o Now adapts to the set of candidate regions provided as input data.
## o Renamed to 'getCopyNumberDataByResampling'.
## o Updated documentation
## o Added an example (requires non-public data for now).
## o Some code cleanups.
## 2012-10-19
## o Added some sanity checks for arguments.
## 2012-10-16
## o Created from randomProfile.R.
############################################################################

