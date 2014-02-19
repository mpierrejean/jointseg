
if (!require("R.menu")) {
  source("http://aroma-project.org/hbLite.R");
  hbLite("R.menu")
  library("R.menu")
}
    
library("jointSeg")
library(DNAcopy)
dataSet <- "CRL2324,BAF"
Chip <- "HumanCNV370v1/"
pct <- c("100","79","50")
pp <- textMenu(pct, value=TRUE)
pathname <- system.file(paste("extdata/", Chip, dataSet,',', pp, ",cnRegions.xdr", sep=""), package="acnr")
dat <- loadObject(pathname)
## - - - - - - - - - - - - - - 
## Parameters of the experiment
## - - - - - - - - - - - - - - 
len <- 200000
K <- 20
B <- 50
regSize <- 0
minL <- 100
normFrac <- NA

simTag <- sprintf("ROC,n=%s,K=%s,regSize=%s,minL=%s,pct=%s", len, K, regSize, minL, pp)
print(simTag)
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

## Statistics to be compared
stats <- c("c,d|het", "c", 'd|het')
names(stats) <- stats


## candidate K
candK <- 10*K

bkpPath <- "bkpData"
bkpPath <- Arguments$getWritablePath(bkpPath)

bpath <- file.path(bkpPath, simName)
bpath <- Arguments$getWritablePath(bpath)

## - - - - - - - - - - -
## Evaluation
## - - - - - - - - - - -
evalForce <- TRUE
tols <- c(1, 2, 5, 10, 20)

evalPath <- "evalData"
evalPath <- Arguments$getWritablePath(evalPath)

epath <- file.path(evalPath, simName)
epath <- Arguments$getWritablePath(epath)


## - - - - - - - - - - -
## Evaluation
## - - - - - - - - - - -
figPath <- "fig"
figPath <- Arguments$getWritablePath(figPath)
