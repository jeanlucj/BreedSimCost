# This is not easy...
# The final two stages will be set together so that the sum is always 1
# The percentages have to include PIC + all VDP stages so length nStages+1
makeSimplexGrid <- function(bsd, justVDP=F){
  minPerc <- bsd$minPercentage
  maxPerc <- bsd$maxPercentage
  percStep <- bsd$percentageStep
  if (justVDP){
    minPerc <- minPerc[-1]
    maxPerc <- maxPerc[-1]
    percStep <- percStep[-1]
  }
  simpDim <- length(minPerc)
  if (simpDim < 2) stop("Not enough stages for simplex")
  lastTwoStep <- min(percStep[simpDim - 0:1])
  
  if (justVDP) strt <- 2 else strt <- 1
  percList <- as.list(seq(from=bsd$minPercentage[strt], 
                          to=bsd$maxPercentage[strt],
                          by=bsd$percentageStep[strt]), ncol=1)
  for (stage in setdiff(1+1:bsd$nStages, strt)){ # if justVDP don't do strt again
    newPerc <- seq(from=bsd$minPercentage[stage], to=bsd$maxPercentage[stage], 
                   by=bsd$percentageStep[stage])
    percList <- mapply(c, rep(percList, each=length(newPerc)), 
                       rep(newPerc, length(percList)), SIMPLIFY=F)
  }
  
  budgets <- lapply(percList, function(v) return(v / sum(v)))
  
  return(budgets)
}#END makeGrid
