filenames <- sprintf("%s,b=%s.xdr", simName, 1:B)
for (bb in 1:B) {
  filename <- filenames[bb]
  print(filename)
  pathname <- file.path(spath, filename)
  sim <- loadObject(pathname)
  if (!is.na(normFrac)) {
    dat <- setNormalContamination(sim$profile, normFrac)
  } else {
    dat <- sim$profile
  }
  CNA.object <- CNA(dat$c,rep(1,len),1:len)
  smoothed.CNA.obj <- smooth.CNA(CNA.object)
  dat$c <- smoothed.CNA.obj$Sample.1
  stats <- c("d", "log(c)")
  for(stat in stats){
    for (KK in candK) {
      methTag <- sprintf("cghseg:%s (Kmax=%s)", stat, KK)
      filename <- sprintf("%s,b=%s,%s.xdr", simNameNF, bb, methTag)
      pathname <- file.path(tpath, filename)
      if (!file.exists(pathname) || segForce) {
        if(stat=="log(c)"){dat$c =  log2(dat$c)-1; stat="c"}
        ## drop NA or -Inf
        dat$c[which(dat$c==-Inf)] <- NA
        indNA <- which(is.na(dat$c))
        posNotNa <-  which(!is.na(dat$c))
        datwithoutNA <- dat[posNotNa,]
        print(stat)
        res <- PSSeg(datwithoutNA, method="DP", K=KK, statistic=stat, profile=TRUE, verbose=TRUE)
        res2 <- list(bestBkp=posNotNa[res$bestBkp], 
                     initBkp=posNotNa[res$initBkp], 
                     dpBkpList=lapply(res$dpBkpList,function(bkp) posNotNa[bkp]), 
                     prof= res$prof)
        saveObject(res2$prof[, "time"], file=pathname)
        
      }
    }
  }
}
