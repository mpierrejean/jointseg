chooseK <- function(lambda, rse, n, c){
    D <- length(rse)
    pen <- lambda*(1:D)*(c+log(n/1:D))/n
    DHat <- which.min(rse/n+pen)
    return(DHat)
}
