filenames <- sprintf("%s,b=%s.xdr", simName, 1:B)
for (bb in 1:B) {
  filename <- filenames[bb]
  print(filename)
  pathname <- file.path(spath, filename)
  if (!file.exists(pathname) || simForce) {
    sim <- getCopyNumberDataByResampling(len, K,minLength=minL, regionSize=regSize, regData=dat, connex=TRUE)
    saveObject(sim, file=pathname)
  }
}
