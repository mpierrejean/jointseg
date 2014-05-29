library(jointSeg)

## Load original data
dat <- loadCnRegionData(tumorFraction=1, platform="GSE29172")

## Length of simulated profile
n <- 25000

## Profile generation
regions <- c("(0,1)", "(0,2)", "(1,1)", "(1,2)")
probs <- c(0.2, 0.4, 0.2, 0.2)
regAnnot <- data.frame(region=regions, freq=probs,stringsAsFactors=FALSE)
print(regAnnot)

dat <- subset(dat, region %in% regions)
sim <- getCopyNumberDataByResampling(length=n, nBkp=20, regAnnot=regAnnot, minLength=100, regData=dat, connex=TRUE)

## Check empirical proportions
table(sim$regions)/length(sim$bkp)
plotSeg(sim$profile, sim$bkp)
