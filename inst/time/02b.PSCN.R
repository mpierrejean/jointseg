filenames <- sprintf("%s,b=%s.xdr", simName, 1:B)
for (bb in 1:B) {
  filename <- filenames[bb]
  print(filename)
  pathname <- file.path(spath, filename)
  sim <- loadObject(pathname)

  methTag <- "PSCN"

  if (!is.na(normFrac)) {
    dat <- setNormalContamination(sim$profile, normFrac)
  } else {
    dat <- sim$profile
  }
  ## drop outliers
  CNA.object <- CNA(dat$c,rep(1,len),1:len)
  smoothed.CNA.obj <- smooth.CNA(CNA.object)
  dat$c <- smoothed.CNA.obj$Sample.1
  filename <- sprintf("%s,b=%s,%s.xdr", simNameNF, bb, methTag)
  pathname <- file.path(tpath, filename)
  if (!file.exists(pathname) || segForce) {
    tt <- try(res <- PSSeg(dat, flavor="PSCN", profile=TRUE, verbose=TRUE,platform = 'Affymetrix'))
    if (class(tt)!="try-error") {
      saveObject(res$prof[, "time"], file=pathname)
    } else {
      warning("Error caught in PSCN segmentation:", bb)
    }
  }
}
