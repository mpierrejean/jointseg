library(jointSeg)
## Load original data
dat <- loadCnRegionData(tumorFraction=1, platform="Affymetrix")
## Length of simulated profile
n <- 100000
## Profile generation
bkp <- c(5000,10000,28000,30000,50000,52000,53000,54000,55000,70000,73000,76000,77000)
regions <- c("(1,1)","(0,1)","(1,1)","(1,3)","(1,1)","(1,3)","(1,2)","(1,3)","(1,2)","(1,1)", "(1,2)","(1,1)","(1,2)","(1,1)")
sim <- getCopyNumberDataByResampling(n, nBkp=14,bkp=bkp,regData=dat, connex=TRUE, regions=regions)
plotSeg(sim$profile)
