## Generation of xtable for Briefings in Bioinformatics.
dataSet <- "GSE29172,ASCRMAv2,H1395vsBL1395"
methTags <-  c("PSCBS",
               sprintf("GFLars+DP:(log(c),d)|het (Kmax=%s)", candK),
               sprintf("RBS+DP:log(c),d|het (Kmax=%s)", candK),
               "CBS:log(c)",
               sprintf("GFLars+DP:log(c) (Kmax=%s)", candK),
               sprintf("RBS+DP:log(c) (Kmax=%s)", candK),
               sprintf("cghseg:log(c) (Kmax=%s)", candK),
               "CBS:d|het",
               sprintf("GFLars+DP:d|het (Kmax=%s)", candK),
               sprintf("RBS+DP:d|het (Kmax=%s)", candK),
               sprintf("cghseg:d|het (Kmax=%s)", candK)
               )


tol <- 5;
aucMeth <- sapply(methTags, function(methTag){
  aucContamination <- NULL
  aucContaminationMom <- NULL
  for(pp in pct){
    simTag <- sprintf("ROC,n=%s,K=%s,regSize=%s,minL=%s,pct=%s,", len, K, regSize, minL, pp)
    simName <- sprintf("%s,%s", dataSet, simTag)
    AucPath <- "aucData"
    fpath <- file.path(AucPath, simName)
    fpath <- Arguments$getWritablePath(fpath)
    filename <- sprintf("%s,B=%s,%s,aucData,relax=%s.xdr", simName, B, methTag, relax)
    pathname <- file.path(fpath, filename)
    if(file.exists(pathname)){
      auc <- loadObject(pathname)
      indtol <- which(tols==tol)
      aucContamination <- c(aucContamination, auc[indtol, "meanAUC"])
      aucContaminationMom <- c(aucContaminationMom, auc[indtol, "masOmenoAUC"])
    }
  }
  return(aucContamination)
})
rownames(aucMeth) <- sprintf("%s", pct)
colnames(aucMeth) <- methTags <-  c("PSCBS", "GFLars", "RBS", "CBS", "GFLars", "RBS", "cghseg", "CBS", "GFLars", "RBS", "cghseg")

order2D <- apply(t(aucMeth)[1:3, ], 2, function(cc){order(cc, decreasing=TRUE)})[1, ]
orderCn <- apply(t(aucMeth)[4:7, ], 2, function(cc){order(cc, decreasing=TRUE)})[1, ]
orderBaf <- apply(t(aucMeth)[8:11, ], 2, function(cc){order(cc, decreasing=TRUE)})[1, ]

methBaf <- c("CBS:d|het",
             sprintf("GFLars+DP:d|het (Kmax=%s)", candK),
             sprintf("RBS+DP:d|het (Kmax=%s)", candK),
             sprintf("cghseg:d|het (Kmax=%s)", candK)
             )[orderBaf]

methCn <-  c("CBS:log(c)",
             sprintf("GFLars+DP:log(c) (Kmax=%s)", candK),
             sprintf("RBS+DP:log(c) (Kmax=%s)", candK),
             sprintf("cghseg:log(c) (Kmax=%s)", candK),
             )[orderCn]

meth2D <-  c("PSCBS",
             sprintf("GFLars+DP:(log(c),d)|het (Kmax=%s)", candK),
             sprintf("RBS+DP:log(c),d|het (Kmax=%s)", candK),
             )[order2D]

matIndOrder <- apply(rbind(meth2D, methCn, methBaf), 2, function(cc){which(methTags%in%cc)})

matCond <- sapply(1:3, function(cc){t(aucMeth)[matIndOrder[, cc], cc]})
rownames(matCond) <- c("c,d", "c", "d")
colnames(matCond) <- sprintf("%s", pct)

library(xtable)
## Product latex table
xtable(t(aucMeth), align =c("l", "r", "r", "r"), caption=sprintf("Data Set 1: pAUC by method"))
xtable(matCond, align =c("l", "r", "r", "r"), caption=sprintf("pAUC By type  : Benchmark 1 :%s", dataSet))
