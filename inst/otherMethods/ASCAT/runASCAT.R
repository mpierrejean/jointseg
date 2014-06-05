library(jointseg)
## Load original data
dat <- loadCnRegionData(tumorFraction=0.7, dataSet="GSE29172")

## simplifying: remove CN probes (otherwise ASCAT gives an error because of its genotyping step)
dat <- subset(dat, !is.na(b))    

n <- 1e4 ## Length of simulated profile
K <- 10 ## Number of breakpoints

## Profile generation
sim <- getCopyNumberDataByResampling(n, 10, regData=dat)
datS <- sim$profile

res <- PSSeg(datS, method="ASCAT")
res$bestBkp

plotSeg(datS, list(sim$bkp, res$bestBkp))
