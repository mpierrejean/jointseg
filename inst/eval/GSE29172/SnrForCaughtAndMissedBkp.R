##----------------------------------------------------------------------------------
## Compute SNR for simulated change-points For Data Set 1
##----------------------------------------------------------------------------------
library(acnr)
library(jointseg)
library("R.menu")
library(xtable)
source("myImagePlot.R")

## Function to compute SNR
SNRFunctionSample <- function(datReg1, datReg2, covariance, reg1, reg2){
  cnReg1 <- log2(datReg1$c)-1
  cnReg2 <- log2(datReg2$c)-1
  dReg1 <- 2*abs(subset(datReg1, genotype==0.5)$b-1/2)
  dReg2 <- 2*abs(subset(datReg2, genotype==0.5)$b-1/2)
  y <- c(cn=mean(cnReg1, na.rm=TRUE)-mean(cnReg2, na.rm=TRUE), d=(mean(dReg1, na.rm=TRUE)-mean(dReg2, na.rm=TRUE)))
  S0 <- covariance[[reg1]]/matrix(c(length(cnReg1), length(dReg1), length(dReg1), length(dReg1)), ncol=2, byrow=TRUE)
  S1 <- covariance[[reg2]]/matrix(c(length(cnReg2), length(dReg2), length(dReg2), length(dReg2)), ncol=2, byrow=TRUE)
  res <- as.numeric(t(y)%*%solve(S0+S1)%*%y)
  return(res)
}

## Data set parameters
dataSet <- "GSE29172,ASCRMAv2,H1395vsBL1395"
Chip <- "GenomeWideSNP_6/"
pp <- "50"
B <- 50
K <- 20
len <- 200000

## Methods which will be evaluated
candK <- 10*20
methTags <- c(sprintf("RBS+DP:log(c),d|het (Kmax=%s)", candK),
              sprintf("GFLars+DP:(log(c),d)|het (Kmax=%s)", candK),
              "PSCBS",
              sprintf("RBS+DP:log(c) (Kmax=%s)", candK),
              sprintf("GFLars+DP:log(c) (Kmax=%s)", candK),
              "CBS:log(c)",
              sprintf("cghseg:log(c) (Kmax=%s)", candK),
              sprintf("RBS+DP:d|het (Kmax=%s)", candK),
              sprintf("GFLars+DP:d|het (Kmax=%s)", candK),
              "CBS:d|het",
              sprintf("cghseg:d|het (Kmax=%s)", candK)
              )
## Studied transitions
regions <- c("(0,1)", "(0,2)", "(1,1)", "(1,2)")
SNR50Meth <- lapply(methTags, function(methTag){
  print(methTag)
  dat <- loadCnRegionData(dataSet="GSE29172", tumorFraction=as.numeric(pp)/100)
  dat$d <- 2*abs(dat$b-1/2)
  dat[which(dat$genotype!=0.5),]$d <- NaN
  dat$c = log2(dat$c)-1
  varCN = as.matrix(by(dat$c, dat$region, var))
  vard = as.matrix(by(dat$d, dat$region, var, na.rm=TRUE))
  covariance <- by(cbind(dat$c,dat$d), dat$region,var,na.rm=TRUE)
  
  for(i in 1:7){
    covariance[[i]][1,1] <- varCN[i,1] 
  }
  
  pathnameDat <-  sprintf("simData/%s,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s", dataSet, pp)
### Compute SNR for each profiles and each bkp
    SNRResults <- lapply(1:B,function(bb){
      dataSample <- loadObject(sprintf("%s/%s,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s,b=%s.xdr",pathnameDat, dataSet, pp,bb))
      start <- c(1,dataSample$bkp[-K]+1)
      inter <- dataSample$bkp
      end <- c(dataSample$bkp[-1],len)
      temp <- mapply( function(ss,ii,ee){
        ## ss:ii : index for s0 (segment 0)
        ## ii+1:ee : index for s1 (segment 1)
        datReg1 = dataSample$profile[ss:ii,]
        datReg2 = dataSample$profile[(ii+1):ee,]
        ## Region s0
        reg1 = dataSample$profile$region[ss]
        ## Region s1
        reg2 = dataSample$profile$region[ee]
        ## Compute SNR for this change-point
        res <- SNRFunctionSample(datReg1 = datReg1,
                                 datReg2 = datReg2,
                                 covariance,
                                 reg1,
                                 reg2)
        if(is.nan(res)){res <- 0}
        return(list(bkp=ii,SNR=res, reg1=reg1, reg2=reg2, len1=(ii-ss), len2=(ee-ii)))
      }, start,inter,end)
      return(temp)
    })
### Sort bkp by missed and caught for each method.
  pathnameBkp <-  sprintf("bkpData/%s,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s", dataSet, pp)
  snrbymeth=data.frame(bkp=NA,SNR=NA,status=NA, reg1 = NA, reg2=NA, len1=NA, len2=NA)
  for(bb in 1:50){
    dataBkp <- loadObject(sprintf("%s/%s,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s,b=%s,%s.xdr",pathnameBkp, dataSet, pp, bb, methTag))
    SNRs <-SNRResults[[bb]]
    trueBkp <- unlist(SNRs[1,])
    candidates <- dataBkp$bestBkp
    status= rep("missed",K)
    res=NULL
    for (cc in candidates){
      distC <- abs(trueBkp-cc)
      minC <- min(distC)
      idx <- which(distC==minC)
      if(minC<=2){
        status[idx] <- "caught"
      }
    }   
    res <- data.frame(bkp=unlist(SNRs["bkp",]),SNR=as.numeric(unlist(SNRs["SNR",])),status=status,reg1=unlist(SNRs["reg1",]),reg2=unlist(SNRs["reg2",]),len1 =unlist(SNRs["len1",]), len2 = unlist(SNRs["len2",]))
    snrbymeth <- rbind(snrbymeth,res)
  }
  snrbymeth
})
names(SNR50Meth) <- methTags

## Lists which contain for missed and caught bkp : position of bkp, region before and after bkp, SNR of bkp.
resBkpMissedBymeth <- list()
resRegMissedBymeth <- list()
resBkpCaughtBymeth <- list()
resRegCaughtBymeth <- list()
for(methTag in methTags){
  resBkpMissedBymeth[[methTag]] <- which(SNR50Meth[[methTag]]$status=="missed")
  resRegMissedBymeth[[methTag]][["reg1"]] <- SNR50Meth[[methTag]]$reg1[which(SNR50Meth[[methTag]]$status=="missed")]
  resRegMissedBymeth[[methTag]][["reg2"]] <- SNR50Meth[[methTag]]$reg2[which(SNR50Meth[[methTag]]$status=="missed")]
  resRegMissedBymeth[[methTag]][["SNR"]] <- SNR50Meth[[methTag]]$SNR[which(SNR50Meth[[methTag]]$status=="missed")]
  
  resBkpCaughtBymeth[[methTag]] <- which(SNR50Meth[[methTag]]$status=="caught")
  resRegCaughtBymeth[[methTag]][["reg1"]] <- SNR50Meth[[methTag]]$reg1[which(SNR50Meth[[methTag]]$status=="caught")]
  resRegCaughtBymeth[[methTag]][["reg2"]] <- SNR50Meth[[methTag]]$reg2[which(SNR50Meth[[methTag]]$status=="caught")]
  resRegCaughtBymeth[[methTag]][["SNR"]] <- SNR50Meth[[methTag]]$SNR[which(SNR50Meth[[methTag]]$status=="caught")]

}
## Bkp Missed
## Mean of SNR for Bkp missed form each method
ListSNR50Missed <- lapply(resRegMissedBymeth,
                          function(tt){
                            mapply(function(ind1,ind2){
                              ## Reg1
                              reg1 <- regions[ind1]
                              ## Reg2
                              reg2 <- regions[ind2]
                              ## Transitions "reg1-reg2"
                              indReg1 <- which(tt$reg1==reg1)
                              indReg2 <- which(tt$reg2==reg2)
                              indReg1Reg2 <- intersect(indReg1,indReg2)
                              ## Transitions "reg2-reg1"
                              indReg1 <- which(tt$reg1==reg2)
                              indReg2 <- which(tt$reg2==reg1)
                              indReg2Reg1 <- intersect(indReg1,indReg2)
                              ## All transitions Missed
                              indFinal <- c(indReg1Reg2,indReg2Reg1)
                              return(matrix(mean(tt$SNR[indFinal]), dimnames = list(reg1, reg2)))
                            }, 1:length(regions), rep(1:length(regions),each=length(regions)))
                          })
## Change into a matrix
matMissed <-lapply(ListSNR50Missed, function(ll){
  matrix(ll, nrow=length(regions), ncol = 4, byrow = FALSE,
         dimnames =list(regions, regions ))
})

tt1 <- c("(0,1)", "(1,1)", "(0,1)", "(0,2)")
tt2 <- c("(0,2)", "(1,2)", "(1,1)", "(1,2)")
matSNRbkpMissed <- sapply(matMissed, function(ll){
  mapply(function(reg1,reg2){
    res <- ll[reg1,reg2]
    names(res) <- sprintf("%s\n%s",reg1,reg2)
    return(res)
  },tt1, tt2, USE.NAMES=FALSE)
})

### Bkp Caught
ListSNR50Caught <- lapply(resRegCaughtBymeth,
                          function(tt){
                            mapply(
                                   function(ind1,ind2){
                                     ## Reg1
                                     reg1 <- regions[ind1]
                                     ## Reg2
                                     reg2 <- regions[ind2]
                                     ## Transitions "reg1-reg2"
                                     indReg1 <- which(tt$reg1==reg1)
                                     indReg2 <- which(tt$reg2==reg2)
                                     indReg1Reg2 <- intersect(indReg1,indReg2)
                                     ## Transitions "reg2-reg1"
                                     indReg1 <- which(tt$reg1==reg2)
                                     indReg2 <- which(tt$reg2==reg1)
                                     indReg2Reg1 <- intersect(indReg1,indReg2)
                                     ## All transitions Caught
                                     indFinal <- c(indReg1Reg2,indReg2Reg1)
                                     return(matrix(mean(tt$SNR[indFinal]), dimnames = list(reg1, reg2)))
                                   }, 1:length(regions), rep(1:length(regions),each=length(regions)))
                          })
## Change into a matrix
matSNRCaught <-lapply(ListSNR50Caught, function(ll){
  matrix(ll, nrow=length(regions), ncol = 4, byrow = FALSE,
         dimnames =list(regions, regions ))
})
matSNRbkpCaught <- sapply(matSNRCaught, function(ll){
  mapply(function(reg1,reg2){
    res <- ll[reg1,reg2]
    names(res) <- sprintf("%s\n%s",reg1,reg2)
    return(res)
  },tt1, tt2, USE.NAMES=FALSE)
})


### Graphics
## Parameters
ylabs <- c("RBS","GFLars","PSCBS","RBS","GFLars","CBS","cghseg","RBS","GFLars","CBS","cghseg")
xlabs <- sprintf("%s\n%s",tt1,tt2)
ymin <- log(min(matSNRbkpMissed, matSNRbkpCaught, na.rm=TRUE))
ymax <- log(max(matSNRbkpMissed, matSNRbkpCaught, na.rm=TRUE))
ab <- c(4.5,8.5)
mar <- c(4,8,2.5,2)

pdf(sprintf("fig/SNRForMissedBkp,%s,pct=50.pdf", dataSet),width=7, height=8.5)
zM <- t(log(matSNRbkpMissed))
myImagePlot(zM, min=ymin , max=ymax, yLabels=ylabs, xLabels=xlabs,title=c(""), ab=ab,mar=mar)
dev.off()

pdf(sprintf("fig/SNRForCaughtBkp,%s,pct=50.pdf", dataSet),width=7, height=8.5)
zC <- zM <- t(log(matSNRbkpCaught))
myImagePlot(zC, min=ymin , max=ymax, yLabels=ylabs, xLabels=xlabs,title=c(""), ab=ab,mar=mar)
dev.off()

###Table
## All transitions observed
resReg <- list(reg1=SNR50Meth[[methTag]]$reg1, reg2=SNR50Meth[[methTag]]$reg2)
TransitionsObserved <- table(resReg[["reg1"]], resReg[["reg2"]])[regions, regions]
AllTransitions <-
  mapply(function(reg1,reg2){
  res <- TransitionsObserved[reg1, reg2] + TransitionsObserved[reg2, reg1]
  names(res) <- paste(reg1, "-", reg2)
  return(res)
},tt1, tt2, USE.NAMES =FALSE)

FinalTable <- t(sapply(methTags, function(methTag){
  matReg <- table(resRegMissedBymeth[[methTag]][["reg1"]], resRegMissedBymeth[[methTag]][["reg2"]])[regions, regions]
  MissedTransitions <-
    mapply(function(reg1, reg2){
      res <- matReg[reg1, reg2] + matReg[reg2, reg1]
      names(res) <- paste(reg1, "-", reg2)
      return(res)
    },tt1, tt2, USE.NAMES =FALSE)
    MissedTransitions/AllTransitions
}))
rownames(FinalTable) <- ylabs
xtable(round(FinalTable, 2))
