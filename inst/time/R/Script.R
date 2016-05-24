### Computing time experiences

pp <- 100
len <- 10^4
onlySNP=FALSE
source('R/00.setup.R')
source('R/01.createDataSets.R')
source("R/02e.PSCBS.R")
source("R/02a.RBS+DP.R")
source("R/02c.cghseg.R")	
source("R/02d.CBS.R")
source("R/02f.GFLars.R")
source("R/02g.CnaStruct.R")
source("R/02i.DPSeg2d.R")
## Run only with version 2.15 ofR
if(version$major<3){source("R/02b.PSCN.R")}

methTags <- c(sprintf("DP:c (Kmax=%s)", candK),
              sprintf("CnaStruct (Kmax=%s)", candK),
              sprintf("DPseg:log(c),d (Kmax=%s)", candK),
              "CBS:log(c)",
              "PSCBS",
              sprintf("RBS+DP:log(c),d (Kmax=%s)", candK),        
              sprintf("GFLars+DP:log(c) (Kmax=%s)", candK),
              sprintf("GFLars+DP:log(c),d (Kmax=%s)", candK),
              "PSCN"
              )
source('R/00.setup.R')
times <- as.matrix(sapply(methTags, function(methTag){
  print(methTag)
  tt <- sapply(1:B, function(bb){
    filename <- sprintf("%s,b=%s,%s.xdr", simNameNF, bb, methTag)
    pathname <- file.path(tpath, filename)
    if(file.exists(pathname)){
      if(methTag%in%c("DP:c (Kmax=50)","DPseg:log(c),d|het (Kmax=50)","CnaStruct (Kmax=50)")){
        time <- loadObject(pathname)
      }else{
        time <- loadObject(pathname)["segmentation"]
      }
      return(time)
    }else{
      return(NA)
    }
  })
  return(mean(tt,na.rm=TRUE))
}), ncol = 1)

library(xtable)

rownames(times) <- gsub("\\+DP","",gsub("\\(Kmax=(.*)\\)","",gsub(".segmentation", "", rownames(times))))
colnames(times) <- sprintf("Time(n = %s)", len)
xtable(times)




pp <- 100
len <- 10^5
onlySNP=FALSE
source('R/00.setup.R')
source('R/01.createDataSets.R')
source("R/02e.PSCBS.R")
source("R/02a.RBS+DP.R")
source("R/02c.cghseg.R")	
source("R/02d.CBS.R")
source("R/02f.GFLars.R")
## Run only with version 2.15 of R
if(version$major<3){source("R/02b.PSCN.R")}

methTags <- c(sprintf("DP:c (Kmax=%s)", candK),
              sprintf("CnaStruct (Kmax=%s)", candK),
              sprintf("DPseg:log(c),d (Kmax=%s)", candK),
              "CBS:log(c)",
              "PSCBS",
              sprintf("RBS+DP:log(c),d (Kmax=%s)", candK),        
              sprintf("GFLars+DP:log(c) (Kmax=%s)", candK),
              sprintf("GFLars+DP:log(c),d (Kmax=%s)", candK),
              "PSCN"
              )
B=10
times2 <- as.matrix(sapply(methTags, function(methTag){
  print(methTag)
  tt <- sapply(1:B,function(bb){
    filename <- sprintf("%s,b=%s,%s.xdr", simNameNF, bb, methTag)
    pathname <- file.path(tpath, filename)
    if(file.exists(pathname)){
      if(methTag%in%c("DP:c (Kmax=50)","DPseg:log(c),d|het (Kmax=50)","CnaStruct (Kmax=50)")){
        time <- loadObject(pathname)
      }else{
        time <- loadObject(pathname)["segmentation"]
      }
    }else{time <- NA}
    time
  })
  return(mean(tt, na.rm=TRUE))
}), ncol = 1)
library(xtable)

rownames(times2) <- gsub("\\+DP","",gsub("\\(Kmax=(.*)\\)","",gsub(".segmentation", "", rownames(times))))
colnames(times2) <- sprintf("Time(n = %s)", len)
xtable(times2)

xtable(cbind(times,times2))
