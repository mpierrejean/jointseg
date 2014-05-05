library(jointSeg)
## Load original data
dat <- loadCnRegionData(tumorFraction=1, platform="Affymetrix")
## Length of simulated profile
n <- 100000
## Profile generation
regions <- c("(1,1)", "(0,1)", "(1,2)", "(0,2)")
probs <- c(0.55, 0.10, 0.25, 0.10)
regAnnot <- data.frame(region=regions, freq=probs,stringsAsFactors=FALSE)
dat <- subset(dat, region%in%regions)
set.seed(1)
sim <- getCopyNumberDataByResampling(length=n, nBkp=10, regAnnot=regAnnot, minLength=100, regData=dat, connex=TRUE)
## Check proportions
table(sim$regions)/10
plotSeg(sim$profile, sim$bkp)
