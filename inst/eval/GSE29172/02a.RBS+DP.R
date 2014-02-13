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
 stats <- c("c","c,d|het", "d|het", "log(c)", "log(c),d|het")
  for (stat in stats) {
    for (KK in candK) {
      methTag <- sprintf("RBS+DP:%s (Kmax=%s)", stat, KK)
      filename <- sprintf("%s,b=%s,%s.xdr", simNameNF, bb, methTag)
      print(filename)
      pathname <- file.path(bpath, filename)
      if (!file.exists(pathname) || segForce) {
        print(stat)
        geno <- dat
        if(length(grep("log",stat))){
          ## Log transformation
          geno$c = log2(geno$c)-1;
          stat= gsub("log\\(c\\)","c", stat);
        }
        geno$c[which(geno$c==-Inf)] <- NA
        indNA <- which(is.na(geno$c))
        posNotNa <-  which(!is.na(geno$c))
        genoWithoutNA <- geno[posNotNa,]
        res <- PSSeg(genoWithoutNA, K=KK, statistic=stat, DP = TRUE, profile=TRUE, verbose=TRUE)
        res <- list(bestBkp=posNotNa[res$bestBkp], 
                     initBkp=posNotNa[res$initBkp], 
                     dpBkpList=lapply(res$dpBkpList,function(bkp) posNotNa[bkp]), 
                     prof= res$prof)
        ##print(getTpFp(res$dpBkpList[[20]], sim$bkp, tol = 5, relax = -1))
        print(res$prof[, "time"])
        saveObject(res, file=pathname)
      }
    }
  }
}
