library(jointSeg)
library(R.utils)
library(xtable)
##- - - - - - - - - - - - - - - - -
## Data generation
## - - - - - - - - - - - - - - - - 
regions <- c("(0,1)", "(0,2)", "(0,3)", "(1,1)", "(1,2)", "(1,3)", "(2,2)", "(2,3)","(0,0)") 
affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
illuDat <- loadCnRegionData(platform="Illumina", tumorFraction=1)

Affymetrix <- sapply(regions,function(rr){length(which(affyDat$region==rr))})
Illumina <- sapply(regions,function(rr){length(which(illuDat$region==rr))})

AnnotedSizeTable <- rbind(Affymetrix,Illumina)

xtable(AnnotedSizeTable, caption = "Size of annoted copy-number regions for each of the 2 data sets")
