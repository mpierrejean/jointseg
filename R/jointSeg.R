jointSeg <- structure(function(# Joint segmentation of multivariate signals
### Joint segmentation of multivariate signals in two steps:
### \enumerate{
###   \item{first-pass segmentation.  By default, a fast, greedy
###     approach is used (see \code{flavor}).}
###   \item{pruning of the candidate change points obtained by
###     dynamic programming}
### }
                               Y,
### The signal to be segmented (must be a matrix)
                               flavor=c("RBS", "GFLars", "PSCN", "cghseg", "CBS", "PSCBS","CnaStruct","PELT"),
### A \code{character} value, the type of segmentation method used:
### \describe{
###   \item{"RBS"}{Recursive Binary Segmentation (the default), see
### \code{\link{segmentByRBS}}}
###   \item{"GFLars"}{Group fused LARS as described in Bleakley and
###   Vert (2011).}
###   \item{"PSCN"}{Hidden Markov Model proposed by Chen et al (2011).
###     Can only be used for copy number signals from SNP arrays}
###   \item{"cghseg"}{Univariate pruned dynamic programming Rigail et al
###     (2010)}
###   \item{"PSCBS"}{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation by Olshen A. et al
###     (2011)}
###   \item{"CnaStruct"}{Bivariate segmentation of SNP-array data for allele-specific copy number analysis in tumour samples by Mosen-Ansorena D. et al
###     (2013)}
###   \item{"PELT"}{Optimal detection of changepoints with a linear computational cost by  Killick R. et al
###     (2012)}}   
                               jitter=NULL,
### Uncertainty on breakpoint position after initial segmentation.  Defaults to NULL.  See Details.
                               methModelSelection='Birge',
### Which method is used to perform model selection
                               DP = TRUE,
### If DP =False, model selection is done on initial segmentation, else model selection is done on segmentation after dynamic programming for flavor RBS
                               ...,
### Further arguments to be passed to the lower-level segmentation
### method determined by argument \code{flavor}.
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

  ##references<<Chen, H., Xing, H., & Zhang, N. R. (2011). Estimation
  ##of parent specific DNA copy number in tumors using high-density
  ##genotyping arrays. PLoS computational biology, 7(1), e1001060.

  ##references<<Rigaill, G. (2010). Pruned dynamic programming for
  ##optimal multiple change-point detection. arXiv preprint
  ##arXiv:1004.0887.

  ##references<<Bellman, R. (1956). Dynamic programming and Lagrange multipliers.
  ##Proceedings of the National Academy of Sciences of the United States
  ##of America, 42(10), 767.

  ##references<<Mosen-Ansorena, D. and Aransay, A.-M.(2013)
  ##Bivariate segmentation of SNP-array data for allele-specific
  ##copy number analysis in tumour samples
  ##BMC Bioinformatics, 14:84 

  ##references<<Killick R, Fearnhead P, Eckley IA (2012)
  ##Optimal detection of changepoints with
  ##a linear computational cost, JASA 107(500),1590-1598  

  ##seealso<<\code{\link{segmentByRBS}}
  ##seealso<<\code{\link{pruneByDP}}
  ##seealso<<\code{\link{segmentByGFLars}}
  ##seealso<<\code{\link{segmentByPSCN}}
  ##seealso<<\code{\link{segmentByPSCBS}}
  ##seealso<<\code{\link{segmentByCBS}}
  ##seealso<<\code{\link{segmentByPelt}}
  ##seealso<<\code{\link{segmentByCnaStruct}}


  flavor <- match.arg(flavor)

  if (flavor %in% c("PSCN","CnaStruct")) {
    ##details<<If \code{flavor=="PSCN"} or  \code{flavor=="CnaStruct"}, \code{Y} should contain the
    ##following columns \describe{
    ## \item{c}{total copy numbers}
    ## \item{b}{allele B fractions (a.k.a. BAF)}
    ## \item{d}{decrease in heterozygosity (\code{2|b-1/2|}) for
    ## heterozygous SNPs}
    ##}
    cn <- colnames(Y)
    ecn <- c("c", "b", "d") ## expected
    mm <- match(ecn, cn)
    if (any(is.na(mm))) {
      str <- sprintf("('%s')", paste(ecn, collapse="','"))
      stop("Argument 'Y' should contain columns named ", str)
    }
    Yseg <- Y[, c("c", "b")]
    Ydp <- Y[, c("c", "d")]
  }
  if (flavor %in% c("PSCBS")) {
     ##details<<If \code{flavor=="PSCBS"}, \code{Y} should contain the
    ##following columns \describe{
    ## \item{c}{total copy numbers}
    ## \item{b}{allele B fractions (a.k.a. BAF)}
    ## \item{d}{decrease in heterozygosity (\code{2|b-1/2|}) for
    ## heterozygous SNPs}
    ## \item{genotype}{genotype}
    ##}
    cn <- colnames(Y)
    ecn <- c("c", "b", "d","genotype") ## expected
    mm <- match(ecn, cn)
    if (any(is.na(mm))) {
      str <- sprintf("('%s')", paste(ecn, collapse="','"))
      stop("Argument 'Y' should contain columns named ", str)
    }
    Yseg <- Y[, c("c", "b","genotype")]
    Ydp <- Y[, c("c", "d")]
  } else {
    Yseg <- Y
    Ydp <- Y
  }
  if (verbose) {
    cat("Start ", flavor, "\n")
  }
  prof <- NULL
  resSeg <- prof(switch(flavor,
                        "RBS"= segmentByRBS(Yseg, ..., verbose=verbose),
                        "GFLars"=segmentByGFLars(Yseg, ..., verbose=verbose),
                        "PSCN"=segmentByPSCN(Yseg, ..., verbose=verbose),
                        "cghseg"=segmentByCghseg(Yseg[, "c"], ..., verbose=verbose),
                        "CBS"=segmentByCBS(Yseg[, "c"], ..., verbose=verbose),
                        "PSCBS"=segmentByPSCBS(Yseg, ..., verbose=verbose),
                        "PELT"=segmentByPelt(Yseg[,"c"],...,verbose=verbose),
                        "CnaStruct"=segmentByCnaStruct(Yseg,...,verbose=verbose)
                        ), doit=profile)
  initSeg <- resSeg$res
  prof <- rbind(prof, segmentation=resSeg$prof)
  if (verbose) {
    cat("end ", flavor, "\n")
    print(resSeg$prof)
  }

  if (flavor %in% c("cghseg","CnaStruct")) {
    ## dynamic programming already run !  Just reshape results.
    dpseg <- initSeg$dpseg
  } else {
    ## Prune candidate breakpoints
    if (verbose) {
      print("Start dynamic programming")
    }
    bkp <- initSeg$bkp
    ##details<<If \code{jitter} is not NULL, it should be a vector of
    ##integer indices. The set of candidate breakpoints passed on to
    ##dynamic programming is augmented by all indices distant from an
    ##element of \code{jitter} from one of the candidates. For
    ##example, if \code{jitter==c(-1, 0, 1)} and the initial set of
    ##breakpoints is \code{c(1,5)} then dynamic programming is run on
    ##\code{c(1,2,4,5,6)}.
    if (!is.null(jitter)) {
      jitter <- as.integer(jitter)
      bkpJ <- sapply(bkp, FUN=function(x) {
        x+jitter
      })
      bkpJ <- unique(bkpJ)  ## remove duplicates
      bkpJ <- bkpJ[bkpJ>1]
      bkpJ <- bkpJ[bkpJ<nrow(Y)]
      bkp <- bkpJ
    }

    resDP <- prof(pruneByDP(Ydp, candCP=bkp), doit=profile)
    dpseg <- resDP$res
    prof <- rbind(prof, dpseg=resDP$prof)
    if (verbose) {
      print("End dynamic programming")
      print(resDP$prof)
    }
  }

  ## Find the best segmentation
  if (flavor%in%c("cghseg","RBS","GFLars") && DP) {
    mS <- modelSelection(dpseg$rse, n=nrow(Y),meth=methModelSelection)
    bestSeg <- integer(0L)
    if (mS$kbest!=0) {
      bestSeg <- dpseg$bkp[[mS$kbest]]
    }
  }else if(flavor=="RBS" && !DP){
    mS <- modelSelection(initSeg$rse, n=nrow(Y),meth=methModelSelection)
    bestSeg <- integer(0L)
     if (mS$kbest!=0) {
      bestSeg <- sort(initSeg$bkp[1:mS$kbest])
    }
  }else{
    bestSeg <- initSeg$bkp
  }
  ##value<< list with elements:
  list(
       bestBkp=bestSeg, ##<< Best set of breakpoints after dynamic programming
       initBkp=initSeg$bkp, ##<< Results of the initial segmentation, using 'segmentByNnn', where 'Nnn' corresponds to argument \code{flavor}
       dpBkpList=dpseg$bkp, ##<< Results of dynamic programming, a list of vectors of breakpoint positions for the best model with k breakpoints for k=1, 2, ... K where \code{K=length(initBkp)}
       prof=prof ##<< a \code{matrix} providing time usage (in seconds) and memory usage (in Mb) for the main steps of the program.  Only defined if argument \code{profile} is set to \code{TRUE}
       )  
}, ex=function(){
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
})

############################################################################
## HISTORY:
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

