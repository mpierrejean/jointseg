library(jointSeg)

## Load original data
dat <- loadCnRegionData(tumorFraction=1, dataSet="GSE29172")

## Length of simulated profile
n <- 1e5

## Profile generation
bkp <- c(5, 10, 28, 30,
         50, 52, 53, 54, 55,
         70, 73, 76, 77
         )*1e3
regions <- c("(1,1)","(0,1)","(1,1)","(1,3)","(1,1)","(1,3)","(1,2)","(1,3)","(1,2)","(1,1)", "(1,2)","(1,1)","(1,2)","(1,1)")

sim <- getCopyNumberDataByResampling(n, bkp=bkp, regData=dat, regions=regions)
plotSeg(sim$profile, sim$bkp)
