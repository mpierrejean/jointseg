filenames <- sprintf("%s,b=%s.xdr", simName, 1:B)
regA <- data.frame(region = unique(dat$region), freq = c(0.90/4,0.90/4,0.10/4,0.90/4,0.90/4,0.10/4,0.10/4,0.10/4))
regA$region <- as.character(regA$region)
for (bb in 1:B) {
  filename <- filenames[bb]
  print(filename)
  pathname <- file.path(spath, filename)
  if (!file.exists(pathname) || simForce) {
    simTag100 <- sprintf("ROC,n=%s,K=%s,regSize=%s,minL=%s,pct=100", len, K, regSize, minL)
    simName100 <- sprintf("%s,%s", dataSet, simTag100)
    filename100 <- sprintf("%s,b=%s.xdr", simName100, 1:B)[bb]
    spath100 <- file.path(simPath, simName100)
    pathname100 <- file.path(spath100, filename100)
    if(file.exists(pathname100)){
      sim100 <- loadObject(pathname100)
      bkp100 <- sim100$bkp
      regions100 <- sim100$region
    }else{
      sim100 <- NULL
      bkp100 <- NULL
      regions100 <- NULL
    }
    sim <- getCopyNumberDataByResampling(len, K, bkp=bkp100, regions=regions100, minLength=minL, regionSize=regSize, regData=dat, connex=TRUE, regAnnot=regA)
    saveObject(sim, file=pathname)
  }
}
