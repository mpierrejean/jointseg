library(jointseg)
## Load original data
dat <- loadCnRegionData(tumorFraction=1, platform="GSE29172")
dat <- subset(dat, !is.na(b))
## Length of simulated profile
n <- 10000
## Profile generation 
sim <- getCopyNumberDataByResampling(n, 10, minLength=500, regData=dat, connex=TRUE)
## Create data to match to ASCAT format
datC <- data.frame(chrs=rep(1,n),pos=1:n,S1=log2(sim$profile$c))
row.names(datC) <- sprintf("SNP%s",1:n)
colnames(datC) <- NULL
write.table(datC,"Tumor_LogR.txt", sep = "\t")
datB <- data.frame(chrs=rep(1,n),pos=1:n,S1=sim$profile$b)
row.names(datB) <- sprintf("SNP%s",1:n)
colnames(datB) <- NULL
write.table(datB,"Tumor_BAF.txt", sep = "\t")

## Source the ASCAT files
source("ascat.R")
source("predictGG.R")
source("aspcf.R")

## Load data with ASCAT function
ascat.bc = ascat.loadData("Tumor_LogR.txt","Tumor_BAF.txt", chrs=1)
## Plot data if you want
ascat.plotRawData(ascat.bc)

## Platforme definition
platform = "AffySNP6"

## Use the genotypes from the original data and format it in a ASCAT file
ascat.gg2 = data.frame(S1=rep(NA,n))
ascat.gg2[!is.na(sim$profile$genotype),"S1"] <- TRUE
row.names(ascat.gg2) <- row.names(datB)

## run ASCAT (this could take time)
ascat.bc <- ascat.aspcf(ascat.bc,ascat.gg=ascat.gg2)
## run RBS
rbs <- PSSeg(sim$profile, method="RBS", K=50, stat=c("c", "d"))

## Get bkp from ASCAT segmentation
bkpCN <- which(diff(ascat.bc$Tumor_LogR_segmented)!=0)

## Evaluation
getTpFp(candidates=bkpCN, relax=-1, tol = 5, trueBkp=sim$bkp)
getTpFp(candidates=rbs$bestBkp, relax=-1, tol = 5, trueBkp=sim$bkp)

## Plot profile with segmentation
plotSeg(sim$profile, breakpoints= list(sim$bkp, bkpCN))
