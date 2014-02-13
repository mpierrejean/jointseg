##----------------------------------------------------------------------------------
## Compute SNR for simulated change-points
##----------------------------------------------------------------------------------

library(acnr)
library(jointSeg)
library("R.menu")
library(xtable)
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

dataSet <- "GSE29172,ASCRMAv2,H1395vsBL1395"
Chip <- "GenomeWideSNP_6/"
pct <- c("100","70","50")

SNRbypctTum <- list()
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

SNR50Meth <- list()
candK <- 10*20
for(methTag in methTags){
  print(methTag)
  pp=50
  print(pp)
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
  
  pathnameDat <-  sprintf("../SimulatedChr1/simData/GSE29172,ASCRMAv2,H1395vsBL1395,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s",pp)
    SNRResults <- lapply(1:50,function(bb){
      dataSample <- loadObject(sprintf("%s/GSE29172,ASCRMAv2,H1395vsBL1395,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s,b=%s.xdr",pathnameDat,pp,bb))
    
    start <- c(1,dataSample$bkp[-20]+1)
    inter <- dataSample$bkp
    end <- c(dataSample$bkp[-1],200000)
    testSNR <-mapply( function(ss,ii,ee){
      
      datReg1 = dataSample$profile[ss:ii,]
      datReg2 = dataSample$profile[(ii+1):ee,]
      reg1 = dataSample$profile$region[ss]
      reg2 = dataSample$profile$region[ee]
      res <- SNRFunctionSample(datReg1 = datReg1,
                               datReg2 = datReg2,
                               covariance,
                               reg1,
                               reg2)
      if(is.nan(res)){res <- 0}
      return(list(bkp=ii,SNR=res, reg1=reg1,reg2=reg2,len1 = (ii-ss), len2 = (ee-ii)))
    }, start,inter,end)
  })
  
  
  pathnameBkp <-  sprintf("../SimulatedChr1/bkpData/GSE29172,ASCRMAv2,H1395vsBL1395,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s",pp)
    snrbymeth=data.frame(bkp=NA,SNR=NA,status=NA, reg1 = NA, reg2=NA, len1=NA, len2=NA)
    for(bb in 1:50){
      dataBkp <- loadObject(sprintf("%s/GSE29172,ASCRMAv2,H1395vsBL1395,ROC,n=2e+05,K=20,regSize=0,minL=100,pct=%s,b=%s,%s.xdr",pathnameBkp,pp,bb,methTag))
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
    
    res <- data.frame(bkp=unlist(SNRs[1,]),SNR=as.numeric(unlist(SNRs[2,])),status=status,reg1=unlist(SNRs[3,]),reg2=unlist(SNRs[4,]),len1 =unlist(SNRs[5,]), len2 = unlist(SNRs[6,]))
    snrbymeth <- rbind(snrbymeth,res)
  }
  i <- which(pct==pp)
  SNR50Meth[[methTag]] <- snrbymeth
  rm(snrbymeth)
}

resBkpMissedBymeth <- list()
resLengthMissedBymeth <- list()
resRegMissedBymeth <- list()
lim <- 0
for(methTag in methTags){
  resBkpMissedBymeth[[methTag]] <- which(SNR50Meth[[methTag]]$status=="missed" & SNR50Meth[[methTag]]$SNR>=lim)
  resLengthMissedBymeth[[methTag]][["len1"]] <-
    mean(SNR50Meth[[methTag]]$len1[which(SNR50Meth[[methTag]]$status=="missed" & SNR50Meth[[methTag]]$SNR>=lim)])
  resLengthMissedBymeth[[methTag]][["len2"]] <-
    mean(SNR50Meth[[methTag]]$len2[which(SNR50Meth[[methTag]]$status=="missed" & SNR50Meth[[methTag]]$SNR>=lim)])
  resRegMissedBymeth[[methTag]][["reg1"]] <-
    SNR50Meth[[methTag]]$reg1[which(SNR50Meth[[methTag]]$status=="missed" & SNR50Meth[[methTag]]$SNR>=lim)]
  resRegMissedBymeth[[methTag]][["reg2"]] <-
    SNR50Meth[[methTag]]$reg2[which(SNR50Meth[[methTag]]$status=="missed" & SNR50Meth[[methTag]]$SNR>=lim)]
  resRegMissedBymeth[[methTag]][["SNR"]] <-
    SNR50Meth[[methTag]]$SNR[which(SNR50Meth[[methTag]]$status=="missed" & SNR50Meth[[methTag]]$SNR>=lim)]
}

resBkpCaughtBymeth <- list()
resLengthCaughtBymeth <- list()
resRegCaughtBymeth <- list()
resReg <- list()
for(methTag in methTags){
  resBkpCaughtBymeth[[methTag]] <- which(SNR50Meth[[methTag]]$status=="catched")
  resLengthCaughtBymeth[[methTag]][["len1"]] <-
    mean(SNR50Meth[[methTag]]$len1[which(SNR50Meth[[methTag]]$status=="catched")])
  resLengthCaughtBymeth[[methTag]][["len2"]] <-
    mean(SNR50Meth[[methTag]]$len2[which(SNR50Meth[[methTag]]$status=="catched")])
  resRegCaughtBymeth[[methTag]][["reg1"]] <-
    SNR50Meth[[methTag]]$reg1[which(SNR50Meth[[methTag]]$status=="catched")]
  resRegCaughtBymeth[[methTag]][["reg2"]] <-
    SNR50Meth[[methTag]]$reg2[which(SNR50Meth[[methTag]]$status=="catched")]
  resReg[["reg1"]] <- SNR50Meth[[methTag]]$reg1
  resReg[["reg2"]] <- SNR50Meth[[methTag]]$reg2
  resRegCaughtBymeth[[methTag]][["SNR"]] <-
    SNR50Meth[[methTag]]$SNR[which(SNR50Meth[[methTag]]$status=="catched" & SNR50Meth[[methTag]]$SNR>=lim)]
}
## Bkp Missed
Transitions <- c("(0,1)","(0,2)","(1,1)","(1,2)")
ListSNR50 <- lapply(resRegMissedBymeth,
       function(tt){
         mapply(
                function(ind1,ind2){
                  reg1 <- Transitions[ind1]
                  reg2 <- Transitions[ind2]
                  res <- NULL
                  indReg1 <- which(tt$reg1==reg1)
                  indReg2 <- which(tt$reg2==reg2)
                  ind <- intersect(indReg1,indReg2)
                  indReg1 <- which(tt$reg1==reg2)
                  indReg2 <- which(tt$reg2==reg1)
                  ind2 <- intersect(indReg1,indReg2)
                  indFinal <- c(ind,ind2)
                  return(matrix(mean(tt$SNR[indFinal]), dimnames = list(reg1, reg2)))
                }, 1:length(Transitions), rep(1:length(Transitions),each=length(Transitions))
                )
       })
matSNR50 <-lapply(ListSNR50, function(x){
  matrix(x, nrow=length(Transitions), ncol = 4, byrow = FALSE,
         dimnames =list(Transitions, Transitions ))
})
matSNRbkpMissed <- sapply(matSNR50, function(ll){
  mapply(function(i,j){
  res <- ll[i,j]
  reg1 <- rownames(ll)[i]
  reg2 <- rownames(ll)[j]
  names(res) <- sprintf("%s\n%s",reg1,reg2)
  return(res)
},c(1,3,1,2), c(2,4,3,4))
})

### Bkp Caught
ListSNR50Caught <- lapply(resRegCaughtBymeth,
       function(tt){
         mapply(
                function(ind1,ind2){
                  reg1 <- Transitions[ind1]
                  reg2 <- Transitions[ind2]
                  res <- NULL
                  indReg1 <- which(tt$reg1==reg1)
                  indReg2 <- which(tt$reg2==reg2)
                  ind <- intersect(indReg1,indReg2)
                  indReg1 <- which(tt$reg1==reg2)
                  indReg2 <- which(tt$reg2==reg1)
                  ind2 <- intersect(indReg1,indReg2)
                  indFinal <- c(ind,ind2)
                  return(matrix(mean(tt$SNR[indFinal]), dimnames = list(reg1, reg2)))
                }, 1:length(Transitions), rep(1:length(Transitions),each=length(Transitions))
                )
       })
matSNR50Caught <-lapply(ListSNR50Caught, function(x){
  matrix(x, nrow=length(Transitions), ncol = 4, byrow = FALSE,
         dimnames =list(Transitions, Transitions ))
})
matSNRbkpCaught <- sapply(matSNR50Caught, function(ll){
  mapply(function(i,j){
  res <- ll[i,j]
  reg1 <- rownames(ll)[i]
  reg2 <- rownames(ll)[j]
  names(res) <- sprintf("%s\n%s",reg1,reg2)
  return(res)
},c(1,3,1,2), c(2,4,3,4))
})

### Graphiques

source('myImagePlot.R')
pdf("fig/SNRForMissedBkp,Affy,pct=50.pdf",width=7, height=8.5)
myImagePlot(t(log(matSNRbkpMissed)),min = min(log(min(matSNRbkpMissed,na.rm=TRUE)),log(min(matSNRbkpCaught,na.rm=TRUE))) ,max=max(max(log(matSNRbkpMissed),log(max(matSNRbkpCaught,na.rm=TRUE))),na.rm=TRUE)
            , yLabels= gsub("\\)","",gsub("\\(","",gsub("log","",gsub("\\+DP","",gsub("\\|het","",gsub("\\(Kmax=200\\)","",colnames(matSNRbkpMissed)))))))
            , xLabels=rownames(matSNRbkpMissed),title=c(''), ab = c(4.5,8.5),mar=c(4,8,2.5,2))

dev.off()
pdf("fig/SNRForCaughtBkp,Affy,pct=50.pdf",width=7, height=8.5)
myImagePlot(t(log(matSNRbkpCaught)),min = min(log(min(matSNRbkpMissed,na.rm=TRUE)),log(min(matSNRbkpCaught,na.rm=TRUE))) ,max=max(max(log(matSNRbkpMissed),log(max(matSNRbkpCaught,na.rm=TRUE))),na.rm=TRUE)
            , yLabels= gsub("\\)","",gsub("\\(","",gsub("log","",gsub("\\+DP","",gsub("\\|het","",gsub("\\(Kmax=200\\)","",colnames(matSNRbkpMissed)))))))
            , xLabels=rownames(matSNRbkpCaught),title=c(''), ab = c(4.5,8.5),mar=c(4,8,2.5,2))
dev.off()



###Tables
matRegTot <- table( resReg[["reg1"]],resReg[["reg2"]])[c(1,2,4,5),c(1,2,4,5)]
matRegTot4Transitions <-
  mapply(function(i,j){
  res <- matRegTot[i,j]+ matRegTot[j,i]
  reg1 <- rownames(matRegTot)[i]
  reg2 <- rownames(matRegTot)[j]
  names(res) <- paste(reg1,"-",reg2)
  return(res)
},c(1,3,1,2), c(2,4,3,4))
res <- NULL
for(methTag in methTags){

  matReg <- table( resRegMissedBymeth[[methTag]][["reg1"]],resRegMissedBymeth[[methTag]][["reg2"]])[c(1,2,4,5),c(1,2,4,5)]
  matRegMissedTransitions <-
    mapply(function(i,j){
      res <- matReg[i,j]+ matReg[j,i]
      reg1 <- rownames(matReg)[i]
      reg2 <- rownames(matReg)[j]
      names(res) <- paste(reg1,"-",reg2)
      return(res)
    },c(1,3,1,2), c(2,4,3,4))
  
  res <- rbind(res,matRegMissedTransitions/matRegTot4Transitions)
}
rownames(res) <- gsub("\\)","",gsub("\\(","",gsub("log","",gsub("\\+DP","",gsub("\\|het","",gsub("\\(Kmax=200\\)","",methTags))))))
xtable(res)
