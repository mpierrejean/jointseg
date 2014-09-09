filenames <- sprintf("%s,b=%s.xdr", simName, 1:B)
warnings("Could take a long time ...")
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
  dat$c <- log(dat$c)
  methTag <- "CnaStruct"
  for (KK in candK) {
    methTag <- sprintf("CnaStruct (Kmax=%s)", KK)
    filename <- sprintf("%s,b=%s,%s.xdr", simNameNF, bb, methTag)
    pathname <- file.path(tpath, filename)
    if (!file.exists(pathname) || segForce) {
      pathCNAstruct <- system.file("otherMethods/CnaStruct", package="jointseg")
      fileCNAstruct <- list.files(pathCNAstruct)
      pathnameCNAstruct <- file.path(pathCNAstruct, fileCNAstruct)
      source(pathnameCNAstruct)
      res <- PSSeg(dat, method="other", segFUN=doCnaStruct, K=KK, profile=TRUE, verbose=FALSE, p = 0.5,maxk=len)
      saveObject(res$prof[, "time"], file=pathname)
    }
  }
}
