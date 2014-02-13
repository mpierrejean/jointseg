if (!require("R.menu")) {
  source("http://aroma-project.org/hbLite.R");
  hbLite("R.menu")
  library("R.menu")
}
library(DNAcopy)
library("jointSeg")
dataSet <- "GSE29172,ASCRMAv2,H1395vsBL1395"
Chip <- "GenomeWideSNP_6/"
##pct <- c("100","70","50","30")
##pp <- textMenu(pct, value=TRUE)
pathname <- system.file(paste("extdata/", Chip,dataSet,',', pp,",cnRegions.xdr", sep=""), package="acnr")
dat <- loadObject(pathname)

## - - - - - - - - - - - - - - 
## Parameters of the experiment
## - - - - - - - - - - - - - - 
##len <- textMenu(c(1e3,2e4,5e4, 1e5), value=TRUE)
##onlySNP <- textMenu(c(TRUE, FALSE), value=TRUE)

K <- 10
B <- 10
regSize <- 0
minL <- 100
normFrac <- NA
if (onlySNP) {
  dat <- subset(dat, !is.na(b))
}
candK <- K*5
## simTag <- sprintf("n=%s,K=%s,regSize=%s,minL=%s", len, K, regSize, minL)
simTag <- sprintf("ROC,n=%s,K=%s,regSize=%s,minL=%s,pct=%s", len, K, regSize, minL,pp)
print(simTag)
##simTag <- sprintf("ROC,n=%s,K=%s,regSize=%s,minL=%s,pct=%s,region=(1,1),(0,1)", len, K, regSize, minL,pp)
##simTag <- sprintf("ROC,n=%s,K=%s,regSize=%s,minL=%s,pct=%s,region=(1,1),(0,1)", len, K, regSize, minL,pp)

if (onlySNP) {
  simTag <- sprintf("%s,onlySNP", simTag)
}

if (!is.na(normFrac)) {
  simTagNF <- sprintf("%s,normFrac=%s", simTag, normFrac)
} else {
  simTagNF <- simTag
}

simName <- sprintf("%s,%s", dataSet, simTag)
simNameNF <- sprintf("%s,%s", dataSet, simTagNF)


## - - - - - - - - - - -
## Simulation
## - - - - - - - - - - -
simForce <- TRUE

simPath <- "simData"
simPath <- Arguments$getWritablePath(simPath)
spath <- file.path(simPath, simName)
spath <- Arguments$getWritablePath(spath)

## - - - - - - - - - - -
## Segmentation
## - - - - - - - - - - -
segForce <- TRUE

bkpPath <- "bkpData"
bkpPath <- Arguments$getWritablePath(bkpPath)

bpath <- file.path(bkpPath, simName)
bpath <- Arguments$getWritablePath(bpath)

## - - - - - - - - - - -
## Time
## - - - - - - - - - - -
timePath <- "timeData"
timePath <- Arguments$getWritablePath(timePath)

tpath <- file.path(timePath, simName)
tpath <- Arguments$getWritablePath(timePath)
