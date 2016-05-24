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
    CNA.object <- CNA(dat$c, rep(1, len), 1:len)
    smoothed.CNA.obj <- smooth.CNA(CNA.object)
    dat$c <- smoothed.CNA.obj$Sample.1
    stats <- list(c("log(c)","d"), "log(c)", "d")
    lapply(stats, function(stat) {
        for (KK in candK) {
            methTag <- sprintf("RBS+DP:%s (Kmax=%s)", paste(stat, collapse=","), KK)
            filename <- sprintf("%s,b=%s,%s.xdr", simNameNF, bb, methTag)
            pathname <- file.path(bpath, filename)
            if (!file.exists(pathname) || segForce) {
                print(stat)
                geno <- dat
                if(length(grep("log", stat))){
                    ## Log transformation
                    geno$c = log2(geno$c)-1;
                    stat= gsub("log\\(c\\)", "c", stat);
                    print(stat)
                }
                res <- PSSeg(geno, method="RBS", K=KK, stat=stat, profile=TRUE, verbose=TRUE)
                print(res$prof[, "time"])
                saveObject(res, file=pathname)
            }
        }
    })
}
