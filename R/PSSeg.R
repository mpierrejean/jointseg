PSSeg <- structure(function(#Parent-Specific copy number segmentation
### This function splits (bivariate) copy number signals into
### parent-specific (PS) segments using recursive binary segmentation
                            data,
### Data frame containing the following columns: \describe{
### \item{c}{Total copy number (non-logged)}
### \item{b}{Allele B fraction}
### \item{genotype}{(germline) genotype of the SNP, coded as 0 for AA, 1/2 for AB, 1 for BB}
### }
### These data are assumed to be ordered by genome position.
                            flavor=c("RBS", "GFLars", "PSCN", "cghseg", "CBS", "PSCBS", "CnaStruct", "PELT", "DP"),
### A \code{character} value, the type of segmentation method used:
### \describe{
###   \item{"RBS"}{Recursive Binary Segmentation (the default), see
### \code{\link{segmentByRBS}}}
###   \item{"GFLars"}{Group fused LARS as described in Bleakley and
###   Vert (2011).}
###   \item{"PSCN"}{Hidden Markov Model proposed by Chen et al (2011)}
###   \item{"cghseg"}{Univariate pruned dynamic programming Rigaill et al (2010)}
###   \item{"PSCBS"}{Parent-specific copy number in paired tumor-normal studies using circular binary segmentation by Olshen A. et al
###     (2011)}
###   \item{"CnaStruct"}{Bivariate segmentation of SNP-array data for allele-specific copy number analysis in tumour samples by Mosen-Ansorena D. et al
###     (2013)}
###   \item{"PELT"}{Optimal detection of changepoints with a linear computational cost by  Killick R. et al
###     (2012)}
###   \item{"DP"}{Dynamic programming}}
                            ## statistic=c("c,d|het", "sqrt(c),d|het", "log(c),d|het", "(c,d)|het", "c|het", "c,(c1,c2)|het", "c|(CN,hom,het),d|het", "c"),
                            statistic=c("c,d|het", "log(c),d|het", "(c,d)|het", "c", "log(c)", "d|het"),
### Statistic to be segmented                            
                            jitter=NULL,
### Uncertainty on breakpoint position after initial segmentation.  See \code{\link{jointSeg}} for details.
                            methModelSelection = 'Birge',
### Which method is used to do the model selection
                            DP=TRUE,
### If DP =False, model selection is done on initial segmentation, else model selection is done on segmentation after dynamic programming for flavor RBS
                            ...,
### Further arguments to be passed to \code{jointSeg}
                            profile=FALSE,
### Trace time and memory usage ?
                            verbose=FALSE
### A \code{logical} value: should extra information be output ? Defaults to \code{FALSE}.
                            ){  
  ##details<<Before segmentation, the input copy number data are
  ##converted to a matrix whose contents depends on the value of
  ##argument \code{statistic}.  By default, i.e. when
  ##\code{(statistic=="c,d|het")}, this matrix contains two columns:
  ##\describe{
  ## \item{c}{total copy numbers, at the "full" resolution}
  ## \item{d}{the decrease in heterozygosity DH defined in Bengtsson et al, 2010,
  ## \code{2|b-1/2|} in order to take advantage of the bimodality of
  ## allele B fractions for heterozygous SNPs.  \code{d} is only
  ## defined for heterozygous SNPs, that is, SNPs for which
  ## \code{data$genotype==1/2}. The rationale for the above "mirror"
  ## transformation is that allele B fractions (\code{b}) are only
  ## informative for heterozygous SNPs (see e.g. Staaf et al, 2008).}
  ##}
  flavor <- match.arg(flavor)
  if (verbose) {
    cat("Flavor: ", flavor, "\n")
  }
  statistic <- match.arg(statistic)
  if (flavor%in% c("PSCN", "PSCBS", "CnaStruct")) {  ## See last round of 'if' statements in the function
    print(paste("Setting 'statistic' to (c,b) for flavor", flavor))  ## PSCN will convert 'c' to log scale
    statistic <- "(c,b)"
  } else if (flavor %in% c("cghseg", "CBS", "PELT")) {
    if (!(statistic %in% c("c", "log(c)", "d|het"))) {
      stop("Argument 'statistic' should be either 'c' or 'log(c)' for flavor ", flavor)
    }
  } else if (flavor=="GFL") {
    if (!(statistic %in% c("c,d|het", "log(c),d|het"))) {
      stop("Missing values are not handled by the current
implementation of group-fused LARS\nPlease choose another statistic
than", statistic)
    }
  }  
  ##details<<The resulting data are then segmented using the
  ##\code{\link{jointSeg}} function for flavors \code{RBS} and
  ##\code{GFLars}, and using the \code{\link[PSCN]{PSCN}} package for
  ##flavor \code{PSCN}.

  ##references<<Bengtsson, H., Neuvial, P., & Speed,
  ##T. P. (2010). TumorBoost: Normalization of allele-specific tumor
  ##copy numbers from a single pair of tumor-normal genotyping
  ##microarrays. BMC bioinformatics, 11(1), 245.

  ##references<<Staaf, J., Lindgren, D., Vallon-Christersson, J.,
  ##Isaksson, A., Goransson, H., Juliusson, G., ... & Ringn\'er,
  ##M. (2008). Segmentation-based detection of allelic imbalance and
  ##loss-of-heterozygosity in cancer cells using whole genome SNP
  ##arrays. Genome Biol, 9(9), R136.

  ##references<<Chen, H., Xing, H., & Zhang, N. R. (2011). Estimation
  ##of parent specific DNA copy number in tumors using high-density
  ##genotyping arrays. PLoS computational biology, 7(1), e1001060.

  ##references<< Bleakley, K., & Vert, J. P. (2011). The group fused
  ##lasso for multiple change-point detection. arXiv preprint
  ##arXiv:1106.4199.
  
  ##references<<Vert, J. P., & Bleakley, K. (2010). Fast detection of multiple
  ##change-points shared by many signals using group LARS. Advances in
  ##Neural Information Processing Systems, 23, 2343-2351.

  ##references<<Rigaill, G. (2010). Pruned dynamic programming for
  ##optimal multiple change-point detection. arXiv preprint
  ##arXiv:1004.0887.

  ##seealso<<\code{\link{jointSeg}}

  prof <- NULL

  n <- nrow(data)
  data$idx <- 1:n
  data$d <- 2*abs(data$b-1/2)
  isHet <- (data[["genotype"]]==0.5)
  data[which(!isHet), "d"] <- NA
  
  if (statistic %in% c("(c,d)|het", "c|het","d|het")) {
    library(matrixStats)  ## for 'binMeans'
    datHet <- data[which(isHet), ]
    CN <- matrix(data$c, ncol=1)
    xOut <-c(1, datHet$idx)
    
    ## Initial smoothing
    if (verbose) {
      print("Start smoothing")
    }
    resS <- prof(binMeans(y=CN, x=data$idx, bx=xOut), doit=profile)
    datHet$cnSmooth <- resS$res
    profS <- resS$prof
    if (verbose) {
      print("End smoothing")
    }
    prof <- rbind(prof, smoothing=profS)

    if (flavor=="PSCN") {
      ## ad hoc: update 'bSeg' so that it is at the heterozygous SNP-level
      bSeg <- bSeg[which(isHet)]  ## Not used ??? /PN,2013-11-29
    }
  } ## if (statistic ...

  Y <- switch(statistic,
              "(c,d)|het"=cbind(c=datHet[, "cnSmooth"], b=datHet[, "d"]),
              "log(c),d|het"=cbind(c=log2(data[, "c"])-1, b=data[, "d"]), 
              "c,d|het"=cbind(c=data[, "c"], b=data[, "d"]), 
              "c"=cbind(c=data[, "c"]),
              "log(c)"=cbind(c=log2(data[, "c"])-1),
              "d|het"=cbind(b=datHet[, "d"])
              )
  pos <- switch(statistic,
                "(c,d)|het"=datHet[, "idx"],
                "c,d|het"=data[, "idx"],
                "log(c),d|het"=data[, "idx"],
                "c"=data[, "idx"],
                "d|het"=cbind(b=datHet[, "idx"])
                )
  if (flavor == "PSCN") {
    Y <- as.matrix(data[, c("c", "b", "d")])             ## 'd' is used by DP 
    pos <- data[, "idx"]
  } else if (flavor == "CnaStruct") {
    Y <- as.matrix(data[, c("c", "b")])
    pos <- data[, "idx"]
  } else if (flavor %in% c("PSCBS")) {
    Y <- as.matrix(data[, c("c", "b", "d","genotype")])  ## 'd' is used by DP
    pos <- data[, "idx"]
  }
  if (verbose) {
    str(Y)
  }
  
  ## Segmentation followed by pruning using dynamic programming
  res <- jointSeg(Y, flavor=flavor, profile=profile, jitter=jitter, DP=DP,verbose=verbose, ...)
  prof <- rbind(prof, res$prof)
  ## back to original positions (in case of smoothing)
  list(
       bestBkp=pos[res$bestBkp], ##<< Best set of breakpoints after dynamic programming
       initBkp=pos[res$initBkp], ##<< Results of the initial segmentation, using 'segmentByNnn', where 'Nnn' corresponds to argument \code{flavor}
       dpBkpList=lapply(res$dpBkpList,function(bkp) pos[bkp]), ##<< Results of dynamic programming, a list of vectors of breakpoint positions for the best model with k breakpoints for k=1, 2, ... K where \code{K=length(initBkp)}
       prof=prof ##<< a \code{matrix} providing time usage (in seconds) and memory usage (in Mb) for the main steps of the program.  Only defined if argument \code{profile} is set to \code{TRUE}
       )
},ex=function(){	
  ## load known real copy number regions
  affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=0.5)

  ## generate a synthetic CN profile
  K <- 10
  len <- 1e5
  sim <- getCopyNumberDataByResampling(len, K, minLength=100, regData=affyDat)
  datS <- sim$profile

  ## run binary segmentation (+ dynamic programming)
  resRBS <- PSSeg(data=datS, K=2*K, profile=TRUE)
  resRBS$prof
  ##getTpFp(resRBS$initSeg$bkp, sim$bkp, 10)
  getTpFp(resRBS$bestBkp, sim$bkp, 10, relax = -1)
  plotSeg(datS, breakpoints=list(sim$bkp, resRBS$bestBkp))
})
############################################################################
## HISTORY:
## 2013-11-29
## Added flavor : 'DP'
## Cleanups in default arguments.
## 2013-03-28
## Added flavors : 'CnaStruct' and 'Pelt'
## 2013-03-07
## Added option 'DP' for flavor "RBS" to do selection on initial segmentation
## 2013-02-27
## o Bug fixed : flavor "GFLars" could not be run at full resolution
## o Added statistic 'd|het'
## 2013-02-26
## o Added option 'jitter' to allow more precise breakpoint identification by DP.
## 2013-02-18
## o Now flavor "GFLars" can be run at full resolution (as 'segmentByGFLars' handles
## missing values). 
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

