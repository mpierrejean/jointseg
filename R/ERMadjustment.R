ERMadjustment <- function(ERM, n){
    D <- length(ERM)
    ind <- (D%/%2):D
    y <- ERM[ind]/n
    var <- ind
    model <- stats::lm(y~var)
    lambda <- -model$coefficient["var"]*n
    return(lambda)
}
