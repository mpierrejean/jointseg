library(jointSeg)
library(R.utils)
len= 1000
nBkp=3
seq(from = 1, to = len, length = 4)
bkp=c(333, 666)
regData100 <- loadCnRegionData(platform="Affymetrix", tumorFraction=1)
regData100 <- subset(regData100, !is.na(genotype) & genotype!=0)
## Normal LOH
region <- c('(1,1)','(1,2)','(1,1)')
sim100 <- getCopyNumberDataByResampling(len, nBkp, bkp,regData=regData100, regions = region)
figPath <- "evalTPFP"
figPath <- Arguments$getWritablePath(figPath)
pdf(file.path(figPath,'TPFP.pdf'), width = 20, height=10)
par(mar = c(5, 5, 6, 2) + 0.3,cex = 1.5,cex.lab = 1,cex.axis = 1.0)
plot(NA, cex = 0.4, ylab = 'copy number', xlab = 'position', pch = 19, ylim = c(0,4), xlim = c(0,len), cex.lab = 2.0)
rect(c(bkp-20),c(-1,-1),c(bkp+20), c(5,5),
     col = c("lightgrey"), border = "lightgrey",
     lty = 1, lwd = 1)
points(sim100$profile$c, cex = 0.5, ylab = 'copy number', xlab = 'position', pch = 19, col = 'black')
#abline(v = bkp, col = 'black', lwd = 5.0)
abline(v = c(667,332) , col = 'blue', lwd = 5.0, lty = 1)
abline(v = c(475,500,804,678) , col = 'forestgreen', lwd = 5.0, lty = 2)
legend(250,5.5,y.intersp=2, x.intersp=0.6,c("true positive","false positive","tolerance area", "true change-point") ,lty = c(1,2,0,0),col = c("blue","forestgreen",NA, "red"), pch=c(NA,NA,NA,"|"),lwd=c(5,5,0,5),xpd=1,bty = "n",ncol = 2, cex = 1.5
, fill =c(0,0,"lightgrey",0)
, border=c(0,0,0,0), merge =TRUE)
axis(side=1, at = bkp, tick = TRUE, lwd.tick = 5, col.ticks='red', labels = c("",""))
mtext(side=1, at = bkp, text = c("t1","t2"), col='red', cex = 2,line = 1)
dev.off()
