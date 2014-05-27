mapPositionsBack <- function(
### Map breakpoint positions back to original space
    pos
    ### A sorted list of position (or indices)
    ) {
  outPos <- sapply(1:(length(pos)-1), FUN=function(ii){
    floor(median(pos[c(ii, ii+1)]))
  })
  outPos <- c(outPos, pos[length(pos)])  ## make sur 'pos' and
  ## 'outPos' are of same length. (This is just formal as by
  ## definition, a break point cannot occur at the last position.)
  outPos
}
