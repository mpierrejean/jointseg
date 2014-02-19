library("jointSeg")
library("R.utils")

SNRFunction <- function(regions, dat, covariance, sampleSize=1000, B=100, llind1, llind2){
  Tv <- mapply(
               function(ind1, ind2){
                 reg1 <- regions[ind1]
                 reg2 <- regions[ind2]
                 res2D <- NULL
                 resCn <- NULL
                 resBaf <- NULL
                 for(bb in 1:B){
                   sampleReg1 <- subset(dat, region==reg1)
                   sampleReg2 <- subset(dat, region==reg2)

                   ## Cn sample
                   PositionSampleReg1 <- sample(size=sampleSize, x=1:nrow(sampleReg1), replace=TRUE)
                   PositionSampleReg2 <- sample(size=sampleSize, x=1:nrow(sampleReg2), replace=TRUE)
                   cnReg1 <- sampleReg1[PositionSampleReg1, ]$c
                   cnReg2 <- sampleReg2[PositionSampleReg2, ]$c

                   ## Proportion of heterozygous SNP in dat by states
                   tt <- table(dat$region, dat$genotype, useNA="ifany")
                   prop <- tt[, 2]/rowSums(tt)

                   ## Baf heterozygous sample
                   n1star <- prop[reg1]*sampleSize
                   n2star <- prop[reg2]*sampleSize
                   PositionSampleReg1 <- sample(size=n1star, x=which(sampleReg1$genotype==0.5), replace=TRUE)
                   PositionSampleReg2 <- sample(size=n2star, x=which(sampleReg2$genotype==0.5), replace=TRUE)
                   dReg1 <- sampleReg1[PositionSampleReg1, ]$d
                   dReg2 <- sampleReg2[PositionSampleReg2, ]$d

                   y <- c(cn=mean(cnReg1, na.rm=TRUE)-mean(cnReg2, na.rm=TRUE),  d=(mean(dReg1)-mean(dReg2)))
                   S0 <- covariance[[reg1]]/matrix(c(sampleSize, n1star, n1star, n1star), ncol=2, byrow=TRUE)
                   S1 <- covariance[[reg2]]/matrix(c(sampleSize, n2star, n2star, n2star), ncol=2, byrow=TRUE)
                   
                   res2D <- c(res2D, as.numeric(t(y)%*%solve(S0+S1)%*%y))
                   resCn <- c(resCn, as.numeric(y["cn"]*(S0[1, 1]+S1[1, 1])^{-1}*y["cn"]))
                   resBaf <- c(resBaf, as.numeric(y["d"]*(S0[2, 2]+S1[2, 2])^{-1}*y["d"]))
                 }
                 name <- paste(reg1, "-", reg2)
                 return(list(cn=resCn, cb=res2D, baf=resBaf))
               }, llind1, llind2
               )
}
### SNR for Illumina data set
dataSet <- "CRL2324,BAF"
Chip <- "HumanCNV370v1/"
pct <- c(100, 79, 50)
regions <- c("(0,1)", "(0,2)", "(1,1)", "(1,2)")
TvIll<- list()
for(pp in pct){
  print(sprintf("pct tum=%s", pp))
  ## Load Illumina data set
  pathname <- system.file(paste("extdata/", Chip, dataSet, ",", pp, ",cnRegions.xdr", sep=""), package="acnr")
  dat <- loadObject(pathname)
  ## Keep only regions : loss/normal/cn-Loh/gain.
  dat <- subset(dat, region %in% regions)
  ## Transform "Baf" into "d" and "Cn" to "log(c)".
  dat$d <- 2*abs(dat$b-1/2)
  dat$c <- log2(dat$c)-1
  ## Keep only heterozygous SNP
  dat[which(dat$genotype!=0.5), ]$d <- NaN
  ## Compute variances
  varCN=as.matrix(by(dat$c, dat$region, var, na.rm=TRUE))
  vard=as.matrix(by(dat$d, dat$region, var, na.rm=TRUE))
  ## Compute covariances
  covariance <- by(cbind(dat$c, dat$d), dat$region, var, na.rm=TRUE)
  for(i in 1:4){
    covariance[[i]][1, 1] <- varCN[i, 1] 
  }
  ## To do only interesting transitions
  llind1 <- c(1, 3, 1, 2, 3)
  llind2 <- c(2, 4, 3, 4, 3)
  colnames <-  mapply(
                      function(ind1, ind2){
                        reg1 <- regions[ind1]
                        reg2 <- regions[ind2]
                        name <- sprintf("%s-%s", reg1, reg2)
                      }, llind1,  llind2
                      )
  TvIll[[which(pct==pp)]] <- SNRFunction(regions, dat, covariance,  sampleSize=500,  llind1=llind1,  llind2=llind2)
  dimnames(TvIll[[which(pct==pp)]]) <- list(c("cn", "(cn,baf)", "baf"), colnames)
}


### SNR for Affymetrix data set
dataSet <- "GSE29172,ASCRMAv2,H1395vsBL1395"
Chip <- "GenomeWideSNP_6/"
TvAff <- list()
pct=c(100, 70, 50)
for(pp in pct){
  print(sprintf("pct tum=%s", pp))
  ## Load Affymetrix data set
  pathname <- system.file(paste("extdata/", Chip, dataSet,",", pp, ",cnRegions.xdr", sep=""), package="acnr")
  dat <- loadObject(pathname)
  ## Keep only regions : loss/normal/cn-Loh/gain.
  dat <- subset(dat, region %in% regions)
  ## Transform "Baf" into "d" and "Cn" to "log(c)".
  dat$d <- 2*abs(dat$b-1/2)
  dat$c <- log2(dat$c)-1
  ## Keep only heterozygous SNP
  dat[which(dat$genotype!=0.5), ]$d <- NaN
  ## Compute variances
  varCN=as.matrix(by(dat$c, dat$region, var, na.rm=TRUE))
  vard=as.matrix(by(dat$d, dat$region, var, na.rm=TRUE))
  ## Compute covariances
  covariance <- by(cbind(dat$c, dat$d), dat$region, var, na.rm=TRUE)
  for(i in 1:4){
    covariance[[i]][1, 1] <- varCN[i, 1] 
  }
  ## To do only interesting transitions
  llind1 <- c(1, 3, 1, 2, 3)
  llind2 <- c(2, 4, 3, 4, 3)
  colnames <-  mapply(function(ind1, ind2){
                        reg1 <- regions[ind1]
                        reg2 <- regions[ind2]
                        name <- sprintf("%s-%s", reg1, reg2)
                      }, llind1, llind2)
    TvAff[[which(pct==pp)]] <- SNRFunction(regions, dat, covariance, sampleSize=500, llind1=llind1, llind2=llind2)
    dimnames(TvAff[[which(pct==pp)]]) <- list(c("cn", "(cn,baf)", "baf"), colnames)
}

### Graphics

## parameters
figPath <- "fig"
figPath <- Arguments$getWritablePath(figPath)
llind1 <- c(1, 3, 1)
llind2 <- c(2, 4, 3)
lwdLine <- 2
lwdSeg <- 1
colCn <- "blue"
colBaf <- "red"

dataSet <- "GSE29172,ASCRMAv2,H1395vsBL1395"
pct <- c(100, 70, 50)
mapply(function(ind1, ind2){
       transit <- sprintf("%s-%s", regions[ind1], regions[ind2])
       figName <- sprintf("SNRByTumorPurity,%s,%s.pdf", dataSet, transit)
       pdf(file.path(figPath, figName), width=5, height=5)
       
       mmCN <- sapply(TvAff, function(ll){mean(log(ll[,transit]$cn))})
       sdCN <- sapply(TvAff, function(ll){sd(log(ll[,transit]$cn))})
       mmBaf <- sapply(TvAff, function(ll){mean(log(ll[,transit]$baf))})
       sdBaf <- sapply(TvAff, function(ll){sd(log(ll[,transit]$baf))})
       
       plot(pct, mmCN, type="l", ylim=c(0, 12), col=colCn, ylab="log(SNR)", lwd=lwdLine, xlab="Tumour purity", lty=1)
       segments(x0=pct, y0=mmCN+sdCN, y1=mmCN-sdCN, col=colCn, lty=1, lwd=lwdSeg)
       points(pct, mmCN+sdCN, pch="-", col=colCn)
       points(pct, mmCN-sdCN, pch="-", col=colCn)
       
       lines(pct, mmBaf, type="l", col=colBaf, ylab="", lwd=lwdLine, lty=2)
       segments(x0=pct, y0=mmBaf+sdBaf, y1=mmBaf-sdBaf, col=colBaf, lty=1, lwd=lwdSeg)
       points(pct, mmBaf+sdBaf, pch="-", col=colBaf)
       points(pct, mmBaf-sdBaf, pch="-", col=colBaf)
       legend("topleft", legend=c("c", "d"), col=c(colCn, colBaf), lty=c(1, 2), bty="n")
       dev.off()
     }, llind1, llind2
       )

dataSet <- "CRL2324,BAF"
pct <- c(100, 79, 50)
mapply(function(ind1, ind2){
       transit <- sprintf("%s-%s", regions[ind1], regions[ind2])
       figName <- sprintf("SNRByTumorPurity,%s,%s.pdf", dataSet, transit)
       pdf(file.path(figPath, figName), width=5,  height=5)		
       mmCN <- sapply(TvIll, function(ll){mean(log(ll[,transit]$cn))})
       sdCN <- sapply(TvIll, function(ll){sd(log(ll[,transit]$cn))})
       mmBaf <- sapply(TvIll, function(ll){mean(log(ll[,transit]$baf))})
       sdBaf <- sapply(TvIll, function(ll){sd(log(ll[,transit]$baf))})
       
       plot(pct, mmCN, type="l", ylim=c(0, 12), col=colCn, ylab="log(SNR)", lwd=lwdLine, xlab="Tumour purity", lty=1)
       segments(x0=pct, y0=mmCN+sdCN, y1=mmCN-sdCN, col=colCn, lty=1, lwd=lwdSeg)
       points(pct, mmCN+sdCN, pch="-", col=colCn)
       points(pct, mmCN-sdCN, pch="-", col=colCn)
       
       lines(pct, mmBaf, type="l", col=colBaf, ylab="", lwd=lwdLine, lty=2)
       segments(x0=pct, y0=mmBaf+sdBaf, y1=mmBaf-sdBaf, col=colBaf, lty=1, lwd=lwdSeg)
       points(pct, mmBaf+sdBaf, pch="-", col=colBaf)
       points(pct, mmBaf-sdBaf, pch="-", col=colBaf)
       legend("topleft", legend=c("c", "d"), col=c(colCn, colBaf), lty=c(1, 2), bty="n")
       dev.off()
     }, llind1, llind2
       )

