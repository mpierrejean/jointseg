url <- "http://heim.ifi.uio.no/bioinf/Projects/ASCAT/ASCAT2.1.zip"
path <- tempdir()

filename <- "ASCAT2.1.zip"
pathname <- file.path(path, filename)

download.file(url, destfile=pathname) 
unzip(pathname, exdir=path)

sourcePath <- file.path(path, "ASCAT2.1")
for (ff in list.files(sourcePath, pattern="*.R")) {
  pathname <- file.path(sourcePath, ff)
  print(paste("Sourcing", ff))
  source(pathname)
}

setwd(sourcePath) ## required by ASCAT...
