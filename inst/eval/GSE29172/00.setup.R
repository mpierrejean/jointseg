if (!require("R.menu")) {
  source("http://aroma-project.org/hbLite.R");
  hbLite("R.menu")
  library("R.menu")
}
library("jointseg")
dataSet <- "GSE29172,ASCRMAv2,H1395vsBL1395"
Chip <- "GenomeWideSNP_6/"
pct <- c("100","70","50","30")
pp <- textMenu(pct, value=TRUE)
dat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=as.numeric(pp)/100)


## - - - - - - - - - - - - - - 
## Parameters of the experiment
## - - - - - - - - - - - - - - 
len <- 2e5
K <- 20
B <- 50
regSize <- 0
minL <- 100
normFrac <- NA
onlySNP <- textMenu(c(TRUE,FALSE), value=TRUE)

if (onlySNP) {
  dat <- subset(dat, !is.na(b))
}
simTag <- sprintf("ROC,n=%s,K=%s,regSize=%s,minL=%s,pct=%s", len, K, regSize, minL, pp)

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

## candidate K
candK <- 10*K

## - - - - - - - - - - -
## Evaluation
## - - - - - - - - - - -
evalForce <- FALSE
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


## - - - - - - - - - - -
## Methods
## - - - - - - - - - - -
stats <- list(c("log(c)","d"), "log(c)", "d")
methTags <- c(
              sapply(stats, function(stat){sprintf("RBS+DP:%s (Kmax=%s)", paste(stat, collapse=","), candK)}),
              sapply(stats, function(stat){sprintf("GFLars+DP:%s (Kmax=%s)", paste(stat, collapse=","), candK)}),
              sprintf("DP:%s (Kmax=%s)", stats[c(2,3)], candK),
              "PSCBS",
              sprintf("CBS:%s", stats[c(2,3)])
)
