library("jointseg")
library("R.utils")

regData100 <- loadCnRegionData(dataSet="GSE29172", tumorFraction=1)
## regData100 <- subset(regData100, !is.na(genotype) & genotype!=0)          

len <- 1000
bkp <- c(333, 666);          ## breakpoints
tp <- c(330, 670);           ## true positives
fp <- c(475, 500, 804, 678); ## false positives

## A copy number transition:
regions <- c("(1,1)", "(1,2)", "(1,1)")
sim100 <- getCopyNumberDataByResampling(len, bkp=bkp, regData=regData100, regions=regions)

figPath <- "fig"
figPath <- Arguments$getWritablePath(figPath)
filename <- "TPFP.pdf"
pathname <- file.path(figPath, filename)

lwd <- 10

pdf(pathname, width=30, height=10)
par(mar=c(5, 5, 3, 2)+0.3, cex=2, cex.lab=2, cex.axis=1.5)
plot(NA, ylab="Copy number", xlab="Position", pch=19, ylim=c(0,4), xlim=c(0,len))
rect(c(bkp-20), c(-1,-1), c(bkp+20), c(5,5),
     col="lightgrey", border="lightgrey",
     lty=1, lwd=1)
points(sim100$profile$c, cex=0.5, ylab="copy number", xlab="position", pch=19, col="black")
abline(v=tp+0.5, col="blue", lwd=lwd, lty=1)
abline(v=fp+0.5, col="forestgreen", lwd=lwd, lty=2)

legend(0, 5,
       c("true positive", "false positive", "tolerance area", "true change-point"),
       lty=c(1, 2, 0, 0),
       col=c("blue", "forestgreen", NA, "red"),
       pch=c(NA, NA, NA, "|"),
       lwd=c(1, 1, 0, 1)*lwd,
       xpd=1, bty="n", horiz=TRUE, cex=1.5,
       fill=c(0, 0, "lightgrey", 0), border=c(0, 0, 0, 0))
axis(side=1, at=bkp, tick=TRUE, lwd.tick=5, col.ticks="red", labels=c("", ""))
mtext(side=1, at=bkp, text=c(expression(t[1]), expression(t[2])), col="red", cex=3, line=1)
dev.off()
