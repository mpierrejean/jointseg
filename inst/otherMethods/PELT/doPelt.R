doPelt <- function(## Run Pelt segmentation,
                                    y,
### A numeric vector or one column matrix, the signal to be segmented
                                    ...
### Parameters for cpt.mean function
                                    ){
  cpt <- changepoint::cpt.mean(y, method="PELT")
  res <- list(bkp=cpt@cpts[-length(cpt@cpts)])
  return(res) 
}

