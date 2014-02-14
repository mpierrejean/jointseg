##----------------------------------------------------------------------------------
## Compute SNR for simulated change-points for Data Set 2
##----------------------------------------------------------------------------------
library(acnr)
library(jointSeg)
library("R.menu")
library(xtable)
source('myImagePlot.R')

## Function to compute SNR
SNRFunctionSample <- function(datReg1, datReg2, covariance,reg1,reg2 ){
  cnReg1 <- log2(datReg1$c)-1
  cnReg2 <- log2(datReg2$c)-1
  dReg1 <- 2*abs(subset(datReg1, genotype==0.5)$b-1/2)
  dReg2 <- 2*abs(subset(datReg2, genotype==0.5)$b-1/2)
  y= c(cn=mean(cnReg1,na.rm=TRUE)-mean(cnReg2,na.rm=TRUE),  d=(mean(dReg1,na.rm=TRUE)-mean(dReg2,na.rm=TRUE)))
  S0=covariance[[reg1]]/matrix(c(length(cnReg1), length(dReg1), length(dReg1), length(dReg1)), ncol = 2, byrow=TRUE)
  S1=covariance[[reg2]]/matrix(c(length(cnReg2), length(dReg2), length(dReg2), length(dReg2)), ncol = 2, byrow=TRUE)
  res <- as.numeric(t(y)%*%solve(S0+S1)%*%y)
  return(res)
}

## Data set parameters
dataSet <- "CRL2324,BAF"
Chip <- "HumanCNV370v1/"
pp="50"
B=50
K <- 20
len <- 200000

## Methods which will be evaluated
candK <- 10*20
methTags <- c(
              sprintf("RBS+DP:log(c),d|het (Kmax=%s)", candK),
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
Transitions <- c("(0,1)","(0,2)","(1,1)","(1,2)")

SNR50Meth <- lapply(methTags, function(methTag){
  print(methTag)
  pathname <- system.file(paste("extdata/", Chip,dataSet,',', pp,",cnRegions.xdr", sep=""), package="acnr")
  dat <- loadObject(pathname)
  dat$d <- 2*abs(dat$b-1/2)
  dat[which(dat$genotype!=0.5),]$d <- NaN
  dat$c = log2(dat$c)-1
  varCN = as.matrix(by(dat$c, dat$region, var))
  vard = as.matrix(by(dat$d, dat$region, var, na.rm=TRUE))
  covariance <- by(cbind(dat$c,dat$d), dat$region,var,na.rm=TRUE)
  
  for(i in 1:7){
    covariance[[i]][1,1] <- varCN[i,1] 
  }
  
  pathnameDat <-  sprintf("../SimulatedChr1Illumina/simData/CRL2324,BAF,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s,ponderation",pp)
### Compute SNR for each profiles and each bkp
  SNRResults <- lapply(1:50,function(bb){
    dataSample <- loadObject(sprintf("%s/CRL2324,BAF,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s,ponderation,b=%s.xdr",pathnameDat,pp,bb))
    
    start <- c(1,dataSample$bkp[-20]+1)
    inter <- dataSample$bkp
    end <- c(dataSample$bkp[-1],200000)
    temp <-mapply( function(ss,ii,ee){
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
      return(list(bkp=ii,SNR=res, reg1=reg1,reg2=reg2,len1 = (ii-ss), len2 = (ee-ii)))
    }, start,inter,end)
    return(temp)
  })
### Sort bkp by missed and caught for each method.  
  pathnameBkp <-  sprintf("../SimulatedChr1Illumina/bkpData/CRL2324,BAF,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s",pp)
  snrbymeth=data.frame(bkp=NA,SNR=NA,status=NA, reg1 = NA, reg2=NA, len1=NA, len2=NA)
  for(bb in 1:50){
    dataBkp <- loadObject(sprintf("%s/CRL2324,BAF,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s,b=%s,%s.xdr",pathnameBkp,pp,bb,methTag))
    SNRs <-SNRResults[[bb]]
    trueBkp <- unlist(SNRs[1,])
    candidates <- dataBkp$bestBkp
    status= rep("missed",20)
    res=NULL
    for (cc in candidates){
      distC <- abs(trueBkp-cc)
      minC <- min(distC)
      idx <- which(distC==minC)
      if(minC<=2){
        status[idx] <- "catched"
      }
    }
    res <- data.frame(bkp=unlist(SNRs["bkp",]),SNR=as.numeric(unlist(SNRs["SNR",])),status=status,reg1=unlist(SNRs["reg1",]),reg2=unlist(SNRs["reg2",]),len1 =unlist(SNRs["len1",]), len2 = unlist(SNRs["len2",]))
    snrbymeth <- rbind(snrbymeth,res)
  }
  snrbymeth
}
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
Transitions <- c("(0,1)","(0,2)","(1,1)","(1,2)")
## Mean of SNR for Bkp missed form each method
ListSNR50Missed <- lapply(resRegMissedBymeth,
                          function(tt){
                            mapply(function(ind1,ind2){
                              ## Reg1
                              reg1 <- Transitions[ind1]
                              ## Reg2
                              reg2 <- Transitions[ind2]
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
                            }, 1:length(Transitions), rep(1:length(Transitions),each=length(Transitions)))
                          })
## Change into a matrix
matMissed <-lapply(ListSNR50Missed, function(ll){
  matrix(ll, nrow=length(Transitions), ncol = 4, byrow = FALSE,
         dimnames =list(Transitions, Transitions ))
})
matSNRbkpMissed <- sapply(matMissed, function(ll){
  mapply(function(reg1,reg2){
    res <- ll[reg1,reg2]
    names(res) <- sprintf("%s\n%s",reg1,reg2)
    return(res)
  },c("(0,1)","(1,1)","(0,1)","(0,2)"), c("(0,2)","(1,2)","(1,1)","(1,2)"), USE.NAMES=FALSE)
})

### Bkp Caught
ListSNR50Caught <- lapply(resRegCaughtBymeth,
                          function(tt){
                            mapply(
                                   function(ind1,ind2){
                                     ## Reg1
                                     reg1 <- Transitions[ind1]
                                     ## Reg2
                                     reg2 <- Transitions[ind2]
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
                                   }, 1:length(Transitions), rep(1:length(Transitions),each=length(Transitions)))
                          })
## Change into a matrix
matSNRCaught <-lapply(ListSNR50Caught, function(ll){
  matrix(ll, nrow=length(Transitions), ncol = 4, byrow = FALSE,
         dimnames =list(Transitions, Transitions ))
})
matSNRbkpCaught <- sapply(matSNRCaught, function(ll){
  mapply(function(reg1,reg2){
    res <- ll[reg1,reg2]
    names(res) <- sprintf("%s\n%s",reg1,reg2)
    return(res)
  },c("(0,1)","(1,1)","(0,1)","(0,2)"), c("(0,2)","(1,2)","(1,1)","(1,2)"), USE.NAMES=FALSE)
})

## Graphics
pdf("fig/SNRForMissedBkp,Illu,pct=50.pdf",width=7, height=8.5)
myImagePlot(t(log(matSNRbkpMissed)),min = min(log(min(matSNRbkpMissed,na.rm=TRUE)),log(min(matSNRbkpCaught,na.rm=TRUE))) ,max=max(max(log(matSNRbkpMissed),log(max(matSNRbkpCaught,na.rm=TRUE))),na.rm=TRUE)
            , yLabels= gsub("\\)","",gsub("\\(","",gsub("log","",gsub("\\+DP","",gsub("\\|het","",gsub("\\(Kmax=200\\)","",colnames(matSNRbkpMissed)))))))
            , xLabels=rownames(matSNRbkpMissed),title=c(''), ab = c(4.5,8.5),mar=c(4,8,2.5,2))
dev.off()
pdf("fig/SNRForCaughtBkp,Illu,pct=50.pdf",width=7, height=8.5)
myImagePlot(t(log(matSNRbkpCaught)),min = min(log(min(matSNRbkpMissed,na.rm=TRUE)),log(min(matSNRbkpCaught,na.rm=TRUE))) ,max=max(max(log(matSNRbkpMissed),log(max(matSNRbkpCaught,na.rm=TRUE))),na.rm=TRUE)
            , yLabels= gsub("\\)","",gsub("\\(","",gsub("log","",gsub("\\+DP","",gsub("\\|het","",gsub("\\(Kmax=200\\)","",colnames(matSNRbkpMissed)))))))
            , xLabels=rownames(matSNRbkpCaught),title=c(''), ab = c(4.5,8.5),mar=c(4,8,2.5,2))

dev.off()

###Table
## All transitions observed
resReg <- list(reg1 = SNR50Meth[[methTag]]$reg1, reg2 = SNR50Meth[[methTag]]$reg2)
TransitionsObserved <- table(resReg[["reg1"]],resReg[["reg2"]])[c("(0,1)", "(0,2)", "(1,1)" ,"(1,2)"),c("(0,1)", "(0,2)", "(1,1)" ,"(1,2)")]
AllTransitions <-
  mapply(function(reg1,reg2){
  res <- TransitionsObserved[reg1,reg2] + TransitionsObserved[reg2,reg1]
  names(res) <- paste(reg1,"-",reg2)
  return(res)
},c("(0,1)","(1,1)","(0,1)","(0,2)"), c("(0,2)","(1,2)","(1,1)","(1,2)"), USE.NAMES =FALSE)

FinalTable <- t(sapply(methTags,function(methTag){
  matReg <- table(resRegMissedBymeth[[methTag]][["reg1"]],resRegMissedBymeth[[methTag]][["reg2"]])[c("(0,1)", "(0,2)", "(1,1)" ,"(1,2)"),c("(0,1)", "(0,2)", "(1,1)" ,"(1,2)")]
  MissedTransitions <-
    mapply(function(reg1,reg2){
      res <- matReg[reg1,reg2]+ matReg[reg2,reg1]
      names(res) <- paste(reg1,"-",reg2)
      return(res)
    },c("(0,1)","(1,1)","(0,1)","(0,2)"), c("(0,2)","(1,2)","(1,1)","(1,2)"), USE.NAMES =FALSE)
    MissedTransitions/AllTransitions
}))
rownames(FinalTable) <- gsub("\\)","",gsub("\\(","",gsub("log","",gsub("\\+DP","",gsub("\\|het","",gsub("\\(Kmax=200\\)","",methTags))))))
xtable(round(FinalTable,2))
