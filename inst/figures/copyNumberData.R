library(jointSeg)
library(R.utils)

## paths
figPath <- "fig"
figPath <- Arguments$getWritablePath(figPath)
figName <- "copyNumberData"

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## Data generation
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
affyDat <- loadCnRegionData(platform="Affymetrix", tumorFraction=0.7)
bkp <- c(2334, 6121)                     ## breakpoint positions
regions <- c("(1,1)", "(1,2)", "(0,2)")  ## copy number regions
sim <- getCopyNumberDataByResampling(1e4, bkp=bkp, regions=regions, regData=affyDat)
## plotSeg(sim$profile, sim$bkp)

dat <- sim$profile
dat$d <- NA
wHet <- which(dat$genotype==1/2)
dat$d[wHet] <- 2*abs(dat$b[wHet]-1/2)

## figure setup
colG <- rep("#88888855", nrow(dat))
colG[wHet] <- "#00000088"

figHeight <- 3.5
cexLab <- 2

## c
datC <- dat
datC$b <- NULL
datC$d <- NULL

filename <- sprintf("%s,c.pdf", figName)
pathname <- file.path(figPath, filename)

pdf(pathname, height=figHeight)
par(cex.lab=cexLab)
plotSeg(datC, col=colG, sim$bkp)
dev.off()

## b
datB <- dat
datB$c <- NULL
datB$d <- NULL

filename <- sprintf("%s,b.pdf", figName)
pathname <- file.path(figPath, filename)

pdf(pathname, height=figHeight)
par(cex.lab=cexLab)
plotSeg(datB, col=colG, sim$bkp)
dev.off()

## d
datD <- dat
datD$c <- NULL
datD$b <- NULL

filename <- sprintf("%s,d.pdf", figName)
pathname <- file.path(figPath, filename)

pdf(pathname, height=figHeight)
par(cex.lab=cexLab)
plotSeg(datD, col=colG, sim$bkp)
dev.off()
