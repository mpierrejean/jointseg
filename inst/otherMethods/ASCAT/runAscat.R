library(jointSeg)
library(R.utils)
dataSet <- "GSE29172,ASCRMAv2,H1395vsBL1395"
Chip <- "GenomeWideSNP_6/"
##pct <- c("100","70","50","30")
##pp <- textMenu(pct, value=TRUE)
pathname <- system.file(paste("extdata/", Chip,dataSet,',', 100,",cnRegions.xdr", sep=""), package="acnr")
dat <- loadObject(pathname)
n <- 100000
sim <- getCopyNumberDataByResampling(n, 10, minLength=500, regData=dat, connex=TRUE)
datC <- data.frame(chrs=rep(1,n),pos=1:n,S1=log2(sim$profile$c))
row.names(datC) <- sprintf("SNP%s",1:n)
write.table(datC,"Tumor_LogR.txt", sep = "\t")
datB <- data.frame(chrs=rep(1,n),pos=1:n,S1=sim$profile$b)
row.names(datB) <- sprintf("SNP%s",1:n)
write.table(datB,"Tumor_BAF.txt", sep = "\t")

source("ascat.R")

ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt", chrs=1)
ascat.plotRawData(ascat.bc)
source("predictGG.R")
source("aspcf.R")
platform = "AffySNP6"

#ascat.gg = ascat.predictGermlineGenotypes(ascat.bc, platform)
ascat.gg2 = data.frame(S1=rep(NA,n))
ascat.gg2[!is.na(sim$profile$genotype),"S1"] <- TRUE
row.names(ascat.gg2) <- row.names(datB)

system.time(ascat.bc <- ascat.aspcf(ascat.bc,ascat.gg=ascat.gg2))
rbs <- PSSeg(sim$profile, K=50)
bkpCN <- which(diff(ascat.bc$Tumor_LogR_segmented)!=0)

getTpFp(candidates=bkpCN, relax=-1, tol = 5, trueBkp=sim$bkp)
getTpFp(candidates=rbs$bestBkp, relax=-1, tol = 5, trueBkp=sim$bkp)

plotSeg(sim$profile, breakpoints= list(sim$bkp, bkpCN))
