library(jointSeg)

## Load original data
dat <- loadCnRegionData(tumorFraction=1, dataSet="GSE29172")

## Length of simulated profile
n <- 100000

## Profile generation 
sim <- getCopyNumberDataByResampling(n, 50, minLength=100, regData=dat)
plotSeg(sim$profile, sim$bkp)

