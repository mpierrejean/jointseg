
################################################################
### AUC along Tolerance
################################################################
for(pp in pct){
    source('R/00.setup.R')
    pathFig <- sprintf("fig/AUC")
    Arguments$getWritablePath(pathFig)
    pathname = sprintf("%s/%s,PartialAUC,2dvs1d.pdf", pathFig, simName)
    pdf(pathname, width = 10, height=10)
    par(cex = 1.5, mar = c(4, 4, 2, 1)+0.1, mgp = c(2.5, 1, 0))
    plot(NA, NA, xlim = c(0, max(tols)), ylim = c(0, 1), xlab = "Tolerance", ylab = "Partial AUC", cex.axis= 1.2, cex.lab = 1.4)
    sapply(methTags, function(methTag){
        AucPath <- "aucData"
        fpath <- file.path(AucPath, simName)
        fpath <- Arguments$getWritablePath(fpath)
        filename <- sprintf("%s,B=%s,%s,aucData,relax=%s.xdr", simNameNF, B, methTag, relax)
        pathname <- file.path(fpath, filename)
        if(file.exists(pathname)){
            auc <- loadObject(pathname)
            mm <- which(methTags==methTag)
            xx <- tols
            lines(xx, auc[, "meanAUC"], col = cols[mm], lty=ltysDP[mm])
            mom <- auc[, "masOmenoAUC"]
            segments(xx, auc[, "meanAUC"]-mom, xx, auc[, "meanAUC"]+mom, col=cols[mm], lwd=1.8)
            points(xx, auc[, "meanAUC"]-mom, pch="-", col=cols[mm])
            points(xx, auc[, "meanAUC"]+mom, pch="-", col=cols[mm])
        }
    })
    dev.off()
}




################################################################
### AUC along Contamination
################################################################

for(tol in tols){
    print("AUC by contamination")
    ss <- sprintf("ROC,n=%s,K=%s,regSize=%s,minL=%s", len, K, regSize, minL)
    pathFig <- sprintf("fig/AUC")
    Arguments$getWritablePath(pathFig)
    pathname <- sprintf("%s/%s,tol=%s,PartialAUC,ContaminationInfluence.pdf", pathFig, figName, tol)
    pdf(pathname, width = 10, height=10)
    par(cex = 1.5, mar = c(4, 4, 2, 1)+0.1, mgp = c(2.5, 1, 0))
    plot(NA, NA, xlim = c(50, 100), ylim = c(0, 1), xlab = "Tumour purity", ylab = "Partial AUC", cex.axis= 1.2, cex.lab = 1.4)
    aucMeth <- sapply(methTags, function(methTag){
        aucContamination <- NULL
        aucContaminationMom <- NULL
        for(pp in pct){
            print(pp)
            ## source("R/00.setup.R")
            simTag <- sprintf("ROC,n=%s,K=%s,regSize=%s,minL=%s,pct=%s", len, K, regSize, minL, pp)
            simName <- sprintf("%s,%s", dataSet, simTag)
            AucPath <- "aucData"
            fpath <- file.path(AucPath, simName)
            fpath <- Arguments$getWritablePath(fpath)
            filename <- sprintf("%s,B=%s,%s,aucData,relax=%s.xdr", simName, B, methTag, relax)
            pathname <- file.path(fpath, filename)
            print(file.exists(pathname))
            if(file.exists(pathname)){
                auc <- loadObject(pathname)
                indtol <- which(tols==tol)
                aucContamination <- c(aucContamination, auc[indtol, "meanAUC"])
                print(aucContamination)
                aucContaminationMom <- c(aucContaminationMom, auc[indtol, "masOmenoAUC"])
            }
        }
        if(!is.null(aucContamination)){
            mm <- which(methTags==methTag)
            xx <- as.numeric(pct)
            print(methTag)
            print(aucContamination)
            lines(xx, aucContamination, col = cols[mm], lty=ltysDP[mm])
            mom <- aucContaminationMom
            segments(xx, aucContamination-mom, xx, aucContamination+mom , col=cols[mm], lwd=1.8)
            points(xx, aucContamination-mom, pch="-", col=cols[mm])
            points(xx, aucContamination+mom, pch="-", col=cols[mm])
            return(aucContamination)
        }
    })
    rownames(aucMeth) <- sprintf("%s", pct)
    dev.off()
}

