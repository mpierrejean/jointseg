for (bb in 1:B) {
  filename <- sprintf("%s,b=%s.xdr", simName, bb)
  print(filename)
  pathname <- file.path(spath, filename)
  sim <- loadObject(pathname)
  if (!is.na(normFrac)) {
    dat <- setNormalContamination(sim$profile, normFrac)
  } else {
    dat <- sim$profile
  }
  ## drop outliers       
  CNA.object <- CNA(dat$c,rep(1,len),1:len)
  smoothed.CNA.obj <- smooth.CNA(CNA.object)
  dat$c <- smoothed.CNA.obj$Sample.1
 # stats <- c("c","c,d|het", "d|het", "log(c)", "log(c),d|het")
  stats <- c("log(c),d|het")
  for (stat in stats) {
    for (KK in candK) {
      methTag <- sprintf("RBS+DP:%s (Kmax=%s)", stat, KK)
      filename <- sprintf("%s,b=%s,%s.xdr", simNameNF, bb, methTag)
      print(filename)
      pathname <- file.path(tpath, filename)
      if (!file.exists(pathname) || segForce) {
        print(stat)
        if(length(grep("log",stat))){
          dat$c = log2(dat$c)-1;
          stat= gsub("log\\(c\\)","c", stat);
          print(stat)
        }
        dat$c[which(dat$c==-Inf)] <- NA
        indNA <- which(is.na(dat$c))
        posNotNa <-  which(!is.na(dat$c))
        datwithoutNA <- dat[posNotNa,]
        res <- PSSeg(datwithoutNA, K=KK, statistic=stat, DP = TRUE, profile=TRUE, verbose=TRUE)
        res <- list(bestBkp=posNotNa[res$bestBkp], 
                     initBkp=posNotNa[res$initBkp], 
                     dpBkpList=lapply(res$dpBkpList,function(bkp) posNotNa[bkp]), 
                     prof= res$prof)
        ##print(getTpFp(res$dpBkpList[[20]], sim$bkp, tol = 5, relax = -1))
        saveObject(res$prof[, "time"], file=pathname)
      }
    }
  }
}
