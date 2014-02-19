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
  stats <- c("c", "d|het", "log(c)")
  for(stat in stats){
  methTag <- sprintf("CBS:%s", stat)
  filename <- sprintf("%s,b=%s,%s.xdr", simNameNF, bb, methTag)
  pathname <- file.path(bpath, filename)
    if (!file.exists(pathname) || segForce) {
      geno <- dat
      if(stat=="log(c)"){
        geno$c <- log2(geno$c)-1
        stat <- "c"
      }
      res <- PSSeg(geno, flavor="CBS", statistic=stat, profile=TRUE, verbose=TRUE)
      print(res$prof[, "time"])
      saveObject(res, file=pathname)
    }
  }
}
