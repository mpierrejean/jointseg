if (!require("R.menu")) {
  source("http://aroma-project.org/hbLite.R");
  hbLite("R.menu")
  library("R.menu")
}
    
library("jointseg")
dataSet <- "CRL2324,BAF"
Chip <- "HumanCNV370v1/"
pct <- c("100","79","50")
pp <- textMenu(pct, value=TRUE)
dat <- loadCnRegionData(dataSet="GSE11976", tumorFraction=as.numeric(pp)/100)
## - - - - - - - - - - - - - - 
## Parameters of the experiment
## - - - - - - - - - - - - - - 
len <- 2e5
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

candK <- 5*K
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
