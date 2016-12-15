#' mapPositionsBack
#' 
#' Map breakpoint positions back to original space
#' 
#' 
#' @param pos A sorted list of position (or indices)
#' @author Morgane Pierre-Jean and Pierre Neuvial
mapPositionsBack <- function(pos) {
    outPos <- sapply(1:(length(pos)-1), FUN=function(ii){
        floor(stats::median(pos[c(ii, ii+1)]))
    })
    outPos <- c(outPos, pos[length(pos)])  ## make sur 'pos' and
    ## 'outPos' are of same length. (This is just formal as by
    ## definition, a break point cannot occur at the last position.)
    outPos
}
