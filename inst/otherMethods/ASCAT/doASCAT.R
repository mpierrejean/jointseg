doASCAT <- function(Y, verbose=FALSE) {
  ## handle "missing values"
  Y <- subset(Y, !is.na(b))
  n <- nrow(Y)
  str(Y)
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Create data to match to ASCAT format
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rnms <- sprintf("SNP%s", 1:n)
  sampleName <- "S1"
  ch <- rep(1,n)
  po <- 1:n
  
  ## total copy numbers
  datC <- data.frame(chrs=ch, pos=po, log2(Y[, "c"]))
  colnames(datC)[3] <- sampleName
  row.names(datC) <- rnms
  pnc <- "Tumor_LogR.txt"
  write.table(datC, pnc, sep="\t")
  
  ## allelic ratios
  datB <- data.frame(chrs=ch, pos=po, Y[, "b"])
  ## datB <- data.frame(chrs=ch, pos=po, rnorm(n, sd=0.1))
  colnames(datB)[3] <- sampleName
  row.names(datB) <- rnms
  pnb <- "Tumor_BAF.txt"
  write.table(datB, pnb, sep="\t")

  ## germline genotypes
  geno <- (Y[, "genotype"] %in% c(0, 1))
  ascat.gg2 <- data.frame(geno)
  colnames(ascat.gg2) <- sampleName
  row.names(ascat.gg2) <- rnms
  str(ascat.gg2)
  
  ## Load data with ASCAT function
  ascat.bc <- ascat.loadData(pnc, pnb, chrs=1)
  file.remove(pnb, pnc)
  str(ascat.bc)
  
  ## run ASCAT (this could take time)
  res <- ascat.aspcf(ascat.bc, ascat.gg=ascat.gg2)
  str(res)
  
  ## clean up files silently created by ASCAT
  filenames <- paste(c("LogR", "BAF"), "PCFed", sampleName, sep="_")
  pathnames <- paste(filenames, ".txt", sep="")
  file.remove(pathnames)

  ## return segmentation
  bkpLogR <- which(diff(res$Tumor_LogR_segmented)!=0)
  bkpBAF <- which(diff(res$Tumor_BAF_segmented[[1]][, 1])!=0)

  list(bkp=bkpLogR)
}
