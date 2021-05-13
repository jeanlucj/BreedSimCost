#' calcBudget function
#'
#' Once the costs are specified in bsd, this function calculates the total annual budget for the breeding scheme
#'
#' @param bsd List of breeding sheme data
#' @return The bsd list with updated intermediate and total costs
#'
#' @details Call this function once costs have been specified
#'
#' @examples
#' bsd <- calcBudget(bsd)
#'
#' @export
calcBudget <- function(bsd, nEntries=NULL, nBreedingProg=NULL){
  # Calculate the program yearly cost
  # Assumptions
  # There is no QC genotyping in the population improvement cycle
  # Variety candidates get QC genotyped at every stage
  # Variety candidates get whole genome genotyped once at the first stage
  # The number of plots in each trial is nEntries*nReps so 
  # the cost of the trial is nEntries*nReps*nLocs*plotCost

  if (is.null(nEntries)) nEntries <- bsd$nEntries
  if (is.null(nBreedingProg)) nBreedingProg <- bsd$nBreedingProg
  
  crossCosts <- nBreedingProg * bsd$nPopImpCycPerYear * bsd$crossingCost
  genoCostsPIC <- nBreedingProg * bsd$nPopImpCycPerYear *
    bsd$wholeGenomeCost
  develCosts <- nEntries[1] * bsd$candidateDevelCost
  genoCostsVDP <- sum(nEntries) * bsd$qcGenoCost + 
    nEntries[1] * bsd$wholeGenomeCost
  trialCosts <- (nEntries * bsd$nReps * bsd$nLocs) %*% bsd$plotCost
  locationCosts <- max(bsd$nLocs) * bsd$perLocationCost
  
  budget <- crossCosts + develCosts + genoCostsVDP +
    genoCostsPIC + trialCosts + locationCosts
  
  bsd$budgetVec <- c(crossCosts, genoCostsPIC, 
                     genoCostsVDP, develCosts,
                     trialCosts, locationCosts,
                     budget)
  names(bsd$budgetVec) <- c("crossCosts", "genoCostsPIC", 
                            "genoCostsVDP", "develCosts",
                            "trialCosts", "locationCosts",
                            "budget")
  return(bsd)
}

#' budgetToScheme function
#'
#' Given breeding scheme budget percentages, return numbers of crosses and plots
#' The percentages are for PIC and for allocation to each VDP stage
#' Rules: 1. For PIC, nBreedingProg has to be at least minNBreedingProg
#' 2. Each stage has to have fewer entries than the previous stage
#' If the percentages break the rules budgetToScheme returns a failure
#' The total budget allowed is given by bsd$budget minus location maintenance
#' 
#' @param percentages Numeric vector with nStages+1 cells which are 
#' alloction for PIC and for each stage.
#' @param bsd List of breeding scheme data
#'  
#' @return A revised bsd with the sizes of the stages fitting the percentages
#'
#' @details Call this function after running specifyCosts.
#'
#' @examples
#' Assume stages of CET, PYT, UYT, so percentRanges needs 4 rows
#' The first stage needs more budget because of genotyping and development costs
#' percentages <- c(0.50, 0.30, 0.10, 0.10)
#' bsd <- budgetToScheme(bsd, percentages=percentages)
#'
#' @export
budgetToScheme <- function(percentages, bsd){
  targetBudget <- bsd$budgetVec["budget"] - bsd$budgetVec["locationCosts"]
  
  # Population Improvement Cycle budget
  picBudget <- targetBudget * percentages[1]
  nBreedingProg <- round(picBudget / bsd$nPopImpCycPerYear /
                           (bsd$crossingCost + bsd$wholeGenomeCost))
  failure <- nBreedingProg < bsd$minNBreedingProg
  
  # Variety Development Pipeline budget
  costPerEntry <- bsd$plotCosts * bsd$nReps * bsd$nLocs + bsd$qcGenoCost
  costPerEntry[1] <- bsd$candidateDevelCost + bsd$wholeGenomeCost
  stageBudgets <- targetBudget * percentages[-1]
  nEntries <- round(stageBudgets / costPerEntry)
  failure <- any(diff(nEntries) > 0) # Ensure that nEntries gets smaller
  
  bsd$realizedBudget <- NA
  if (!failure){
    bsd$realizedBudget <- calcBudget(bsd, nEntries, nBreedingProg)
    bsd$nBreedingProg <- nBreedingProg
    bsd$nEntries <- nEntries
  }
  bsd$failure <- failure
  return(bsd)
}

#' gridSearch function
#'
#' Minimum and maximum budget percentages for each part and step sizes
#' These vectors are in bsd
#' Returns a list with all the outputs
#' 
#' @param bsd List of breeding scheme data
#'  
#' @return List with initial and end means and variances
#'
#' @details Call this function to set up optimization
#'
#' @examples
#' bsd <- gridSearch(bsd)
#'
#' @export
makeGrid <- function(bsd){
  
  percList <- as.list(seq(from=bsd$minPercentage[1], to=bsd$maxPercentage[1], 
                 by=bsd$percentageStep[1]), ncol=1)
  for (stage in 1+1:bsd$nStages){
    newPerc <- seq(from=bsd$minPercentage[stage], to=bsd$maxPercentage[stage], 
                   by=bsd$percentageStep[stage])
    percList <- mapply(c, rep(percList, each=length(newPerc)), 
                       rep(newPerc, length(percList)), SIMPLIFY=F)
  }
  
  budgets <- lapply(percList, function(v) return(v / sum(v)))
  
  return(budgets)
}#END gridSearch

#' runWithBudget function
#'
#' Run one replication of the scheme with specified budget percentages
#' These percentages have to be valid in terms of following rules
#' 
#' @param percentages Double vector 1 + nStages with budget percentages for 
#' the population improvement cycle and the variety development pipeline
#' @param bsd List of breeding scheme data
#'  
#' @return Vector with initial and end means and variances
#'
#' @details Call this function to test out the response of different budgets
#'
#' @examples
#' bsd <- runWithBudget(percentages, bsd)
#'
#' @export
runWithBudget <- function(percentages, bsd){
  bsd <- budgetToScheme(bsd, percentages)
  if (bsd$failure) return(c(percentages, NA))
  
  startValues <- calcCurrentStatus(bsd)
  
  for (twoPart in 1:bsd$nCyclesToRun){
    bsd$year <- bsd$year+1
    bsd <- makeVarietyCandidates(bsd)
    
    # Within the year, sequentially run the VDP...
    for (stage in 1:bsd$nStages){
      entries <- chooseTrialEntries(bsd, toTrial=bsd$stageNames[stage],
                   fromTrial=ifelse(stage == 1, NULL, bsd$stageNames[stage-1]))
      bsd <- runVDPtrial(bsd, bsd$stageNames[stage], entries)
    }

    # ... And then run the PIC. Sequencing could be made more complicated
    for (genSel in 1:bsd$nPopImpCycPerYear){
      optCont <- selectParents(bsd)
      bsd <- makeCrosses(bsd, optCont)
    }
  }
  
  endValues <- calcCurrentStatus(bsd)
  
  return(c(percentages, startValues, endValues))
}#END runWithBudget

#' runBatch function
#'
#' Run a batch of replications with specified budget percentages
#' These percentages have to be valid in terms of following rules
#' 
#' @param batchBudg List of budget percentage vectors for 
#' the population improvement cycle and the variety development pipeline
#' @param bsd List of breeding scheme data
#'  
#' @return Tibble with initial and end means and variances
#'
#' @details Call this function to test out a batch of different budgets
#'
#' @examples
#' bsd <- runBatch(batchBudg, bsd)
#'
#' @export
runBatch <- function(batchBudg, bsd){
  require(parallel)
  
  batchResults <- mclapply(batchBudg, runWithBudget,
                           bsd=bsd, mc.cores=bsd$nCores)
  # Remove results where budget was not valid
  batchResults <- batchResults[sapply(batchResults, 
                                      function(v) !any(is.na(v)))]
  nParmsInResults <- length(batchResults[[1]])
  batchResults <- batchResults %>% unlist %>% matrix(nrow=nParmsInResults) %>% 
    t %>% as_tibble(.name_repair="universal")
  # WARNING numbers here are hardcoded based on calcCurrentStatus
  colnames(batchResults) <- c(paste0("perc", c("PIC", bsd$stageNames)),
                              paste0("init", c("PopMean", "PopSD", "VarMean")),
                              paste0("end", c("PopMean", "PopSD", "VarMean")))
  batchResults <- batchResults %>% mutate(response = endVarMean - initPopMean)
  
  return(batchResults)
}#END runBatch

#' findRedoBudgets function
#'
#' Based on simulations run so far, find the budgets that gave the highest
#' gain and the budgets with high gain but also high uncertainty
#' 
#' @param bsd List of breeding scheme data
#' @param simResults List with the results of simulations so far
#'  
#' @return List of two vectors one with high gain budgets
#' and one with uncertainty budgets
#'
#' @details Call this function to choose which budgets to try again
#'
#' @examples
#' bsd <- findRedoBudgets(percentages, bsd)
#'
#' @export
findRedoBudgets <- function(simResults, bsd){
  # Non-Parametric LOESS response
  loFormula <- paste0("response ~ ", 
                      paste0(colnames(df)[1:bsd$nStages], collapse=" + "))
  loFM <- loess(loFormula, data=simResults, degree=1)
  loPred <- predict(loFM, se=T)
  loPred <- tibble(fit=loPred$fit, se=loPred$se.fit, simNum=1:nrow(simResults))
  # Budgets with the highest response
  whichBest <- order(loPred$fit, decreasing=T)[1:bsd$nHighGain]
  # Budgets with the high response and high uncertainty
  ndRows <- findNonDom(loPred, dir1Low=F, dir2Low=F, 
                     var1name="fit", var2name="se")
  # Remove ones too far from the best
  fitDist <- loPred$fit[whichBest[1]] - ndRows$fit
  keep <- which(fitDist < ndRows$se*4)
  ndRows <- ndRows[keep,]
  if (nrow(ndRows) > bsd$nUncertain){
    ndRows <- ndRows[order(ndRows$fit, decreasing=T)[1:nByPareto],]
  }
  whichUncert <- ndRows$simNum
  
  return(list(whichBest=whichBest, whichUncert=whichUncert))
}#END findRedoBudgets

#' optimizeByLOESS function
#'
#' Function to optimize a two-part strategy breeding scheme:
#' 1. Simulate a batch using given percentage ranges
#' 2. Perform LOESS fit to the gains
#' 3. Find budget with best estimated gain
#' 4. Calculate new percentage ranges: any simulation within 2*StdErr of best
#' 5. Decide on some simulations to repeat:
#'    1. Parameter space with high gain and high std err: need more info there
#'    2. Parameter space with high gain: high probability that it's best
#' Go back to 1.
#'
#' @param batchSize Integer number of simulations between LOESS fits
#' @param tolerance Numerical difference between min amd max 
#' percentage budgets for all stages
#' @param baseDir Directory if you want to have progress saved by batch. 
#' Relative to R working directory. If not empty string, include final /
#' @param maxNumBatches Integer to stop the simulations eventually if
#' the algorithm is not narrowing in on optimal parameter values
#' @return Numeric matrix with all simulation budget allocations, 
#' gen mean change, gen std dev change, total cost.
#' 
#' @details A wrapper to repeatedly simulate a scheme with different 
#' budget allocations to find optimal allocations
#' 
#' @examples
#' 
#' @export
optimizeByLOESS <- function(bsd){
  require(here)
  on.exit(expr={
      print(traceback())
      saveRDS(mget(ls()), file=here::here("data/optimizeByLOESS.rds"))
    }
  )
  
  nBatchesDone <- 0
  toleranceMet <- FALSE
  allPercRanges <- list()
  results <- tibble()
  # 1. Make grid batch
  batch <- makeGrid(bsd)
  while (batchesDone < maxNumBatches & !toleranceMet){
    # 2. Run batch
    newBatchOut <- runBatch(batch, bsd)
    # 3. Evaluate stopping rule
    results <- results %>% bind_rows(newBatchOut)
    # Non-Parametric LOESS response
    loFormula <- paste0("response ~ ", 
                   paste0(colnames(results)[1:bsd$nStages], collapse=" + "))
    loFM <- loess(loFormula, data=results, degree=1)
    loPred <- predict(loFM, se=T)
    loPred <- tibble(simNum=1:nrow(results), fit=loPred$fit, se=loPred$se.fit, results)
    # 4. Evaluate stopping rule
    bestGain <- which.max(loPred$fit)
    bestSE <- loPred$se[bestGain]
    bestClose <- which(loPred$fit[bestGain] - loPred$fit < bestSE)
    percRanges <- results %>% slice(bestClose) %>% 
      dplyr::select(contains("perc")) %>% summarise_all(range)
    toleranceMet <- all(percRanges %>% summarise_all(diff) < tolerance)
    if (!toleranceMet){
      # 5. Find which budgets to rerun
      redoSims <- findRedoBudgets(results, bsd) %>% unlist
      batch <- lapply(redoSims, function(idx) 
        results %>% slice(idx) %>% pull(1:(bsd$nStages+1)))
      combineBudg <- function(dummy, batch){
        twoBudg <- sample(length(batch), 2)
        return(0.9 * batch[[twoBudg[1]]] + 0.1 * batch[[twoBudg[2]]])
      }
      batch <- c(batch, 
                 lapply(1:(bsd$batchSize - length(batch)), combineBudg, batch))
    }
    allPercRanges <- c(allPercRanges, list(percRanges, 
      nSimClose=length(bestClose), bestGain=loPred$fit[bestGain],
      bestSE=bestSE))
    
    batchesDone <- batchesDone + 1
    
    # Save batches and results
    if (!is.null(baseDir)){
      saveRDS(results, file=here::here("data/allBatches.rds"))
      saveRDS(allPercRanges, file=here::here("data/allPercentRanges.rds"))
    }
  }#END go through batches

  on.exit()
  return(list(results=results, allPercentRanges=allPercRanges))
}
