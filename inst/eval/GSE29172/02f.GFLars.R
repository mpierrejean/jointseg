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
  ## drop outliers
  CNA.object <- CNA(dat$c,rep(1,len),1:len)
  smoothed.CNA.obj <- smooth.CNA(CNA.object)
  dat$c <- smoothed.CNA.obj$Sample.1
  stats <- list(c("log(c)","d"), "log(c)", "d")
  lapply(stats, function(stat) {
    for (KK in candK) {
      methTag <- sprintf("GFLars+DP:%s (Kmax=%s)",  paste(stat, collapse=","), KK)
      filename <- sprintf("%s,b=%s,%s.xdr", simNameNF, bb, methTag)
      pathname <- file.path(bpath, filename)
      if (!file.exists(pathname) || segForce) {
        geno <- dat
        if(length(grep("log", stat))){
          geno$c <- log2(geno$c)-1;
          stat <- gsub("log\\(c\\)", "c", stat)
        }
        ## drop NA or -Inf
        geno$c[which(geno$c==-Inf)] <- NA
        indNA <- which(is.na(geno$c))
        posNotNa <-  which(!is.na(geno$c))
        res <- PSSeg(geno, method="GFLars", K=KK, stat=stat, profile=TRUE, verbose=FALSE)
        print(res$prof[, "time"])
        saveObject(res, file=pathname)
      }  
    }
  })
}

