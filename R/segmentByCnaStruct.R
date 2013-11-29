segmentByCnaStruct <- structure(function(## Run cnaStruct segmentation,
### This function is a wrapper for convenient use of the \code{CnaStruct}
### segmentation method by \code{\link{PSSeg}}.  It applies the
### \code{breakpoints} function from package \code{CnaStruct} and reshapes
### the results.
                                         Y,
### The signal to be segmented, a matrix containing the following columns: \describe{
### \item{c}{Total copy number (non-logged scale)}
### \item{b}{Allele B fraction (a.k.a. BAF)}
### }
                                         K,
### The maximum number of segments to find,
                                         p=0.5,
### Order of Minkowski distances
                                         homthr=0.95,
### homthr sets the BAF limit to consider a SNP homozygous
                                         beta= 1,                                     
### beta controls the relative importance of BAF with respect to LRR
                                         s=1,
### BIC criterion for the best model if s=1
                                         maxk=1000,
### maxk is the maximum segment length allowed
                                         ...,
                                         verbose = FALSE
                                         ){
  ##references<<Mosen-Ansorena, D. and Aransay, A.-M.(2013)
  ##Bivariate segmentation of SNP-array data for allele-specific
  ##copy number analysis in tumour samples
  ##BMC Bioinformatics, 14:84  
  
  if (!require("CnaStruct")) {
    stop("Please install the 'CnaStruct' package to run the 'segmentByCnaStruct' function\nURL:http://genomics.cicbiogune.es/cnastruct")
  }
  LRR <- log2(Y[, "c"])-1
  BAF <- Y[, "b"] ## TO DO: invert BAF to beta transformation ?

  segm <- CnaStruct:::segment(LRR, BAF, maxseg=(K+1), maxk=maxk, p=p, homthr=homthr, beta = beta)
  bestll = which.max(CnaStruct:::logLik(segm) - (1:length(CnaStruct:::logLik(segm))) * log(length(Y[,1])) * s)
  bkp <- as.vector(segm@breakpoints[[bestll]])
  dpseg <-  list(bkp=lapply(segm@breakpoints[-1], as.vector))
  res <- list(bkp=bkp, dpseg= dpseg) 
  return(res)
}, ex=function(){
  if (require("CnaStruct")) {
    ## load known real copy number regions
    IlluDat <- loadCnRegionData(platform="Illumina", tumorFraction=0.5)
    
    ## generate a synthetic CN profile
    K <- 3
    len <- 1000
    sim <- getCopyNumberDataByResampling(len, K, minLength=10, regData=IlluDat)
    datS <- sim$profile
    
    ## run CnaStruct segmentation
    Y <- datS[, c("c", "b")]
    res <- jointSeg:::segmentByCnaStruct(Y, K=K*10, maxk=500)  
    getTpFp(res$bkp, sim$bkp, tol=5, relax = -1)   ## true and false positives
    plotSeg(datS, breakpoints=list(sim$bkp, res$bkp))
  }
})
############################################################################
## HISTORY:
## 2013-11-29
## o Now working on the log scale internally for total copy numbers.  Non-
##   logged values are still expected as an input.
## 2013-05-16
## o Example code now embedded in a 'require()' statement to avoid
##   problems in the R CMD check mechanism of R-forge.
## 2013-03-27
## o Created.
############################################################################




