#' calcBudget function
#'
#' Once the costs are specified in bsd, this function calculates 
#' the total annual budget for the breeding scheme
#' Assumptions
#' There is no QC genotyping in the population improvement cycle
#' Variety candidates get QC genotyped at every stage
#' Variety candidates get whole genome genotyped once at the first stage
#' The number of plots in each trial is nEntries*nReps so 
#' the cost of the trial is nEntries*nReps*nLocs*plotCost
#'
#' @param bsd List of breeding scheme data
#' @return Numeric vector with intermediate and total costs
#'
#' @details Call this function once costs have been specified
#'
#' @examples
#' bsd$budgetVec <- calcBudget(bsd)
#'
#' @export
calcBudget <- function(bsd){
  # Per task way of calculating budget
  crossCosts <- bsd$nBreedingProg * bsd$nPopImpCycPerYear * bsd$crossingCost
  genoCostsPIC <- bsd$nBreedingProg * bsd$nPopImpCycPerYear *
    bsd$wholeGenomeCost
  develCosts <- bsd$nEntries[1] * bsd$candidateDevelCost
  genoCostsVDP <- sum(bsd$nEntries) * bsd$qcGenoCost + 
    bsd$nEntries[1] * bsd$wholeGenomeCost
  trialCosts <- (bsd$nEntries * bsd$nReps * bsd$nLocs) %*% bsd$plotCost
  locationCosts <- max(bsd$nLocs) * bsd$perLocationCost

  # Per stage way of calculating budget
  # Population Improvement Cycle budget
  perPICprog <- bsd$crossingCost + bsd$wholeGenomeCost
  picBudget <- bsd$nBreedingProg * bsd$nPopImpCycPerYear * perPICprog
  
  # Variety Development Pipeline budget
  costPerEntry <- bsd$plotCosts * bsd$nReps * bsd$nLocs + bsd$qcGenoCost
  costPerEntry[1] <- costPerEntry[1] + 
    bsd$candidateDevelCost + bsd$wholeGenomeCost
  stageBudgets <- bsd$nEntries * costPerEntry
  
  # Minimal ratios to ensure that stages have fewer entries going forward
  minRatios <- c(perPICprog, costPerEntry[-bsd$nStages]) / costPerEntry

  budget <- crossCosts + develCosts + genoCostsVDP +
    genoCostsPIC + trialCosts + locationCosts
  
  # Sanity check
  if (abs(budget - sum(picBudget, stageBudgets, locationCosts)) > 1e-6){
    stop("There was a problem calculating the budget")
  }
  
  minPICbudget <- bsd$minNBreedingProg * bsd$nPopImpCycPerYear *
    sum(bsd$crossingCost + bsd$wholeGenomeCost) / sum(picBudget, stageBudgets)
  minLastStgBudget <- bsd$nToMarketingDept * 
    ((bsd$plotCosts * bsd$nReps * bsd$nLocs)[bsd$nStages] + bsd$qcGenoCost) / 
    sum(picBudget, stageBudgets)
  
  budgetVec <- c(crossCosts, genoCostsPIC, 
                 genoCostsVDP, develCosts,
                 trialCosts, locationCosts,
                 picBudget, stageBudgets,
                 c(picBudget, stageBudgets) / sum(picBudget, stageBudgets),
                 minPICbudget, minLastStgBudget, minRatios,
                 budget)
  names(budgetVec) <- c("crossCosts", "genoCostsPIC", 
                        "genoCostsVDP", "develCosts",
                        "trialCosts", "locationCosts",
                        "picBudget", paste0(bsd$stageNames, "_budget"),
                        paste0("perc_", c("PIC", bsd$stageNames)),
                        "minPICbudget", "minLastStgBudget",
                        paste0("ratio_", 1:bsd$nStages-1, 1:bsd$nStages),
                        "budget")
  return(budgetVec)
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
  if (bsd$debug){
    require(here)
    on.exit(expr={
      print(traceback())
      saveRDS(mget(ls()), file=here::here("data", "budgetToScheme.rds"))
    })
  }
  
  targetBudget <- bsd$initBudget["budget"] - bsd$initBudget["locationCosts"]
  
  numbersFromPercentages <- function(percentages){
    # Population Improvement Cycle budget
    picBudget <- targetBudget * percentages[1]
    nBreedingProg <- ceiling(picBudget / bsd$nPopImpCycPerYear /
                             (bsd$crossingCost + bsd$wholeGenomeCost))
    
    # Variety Development Pipeline budget
    costPerEntry <- bsd$plotCosts * bsd$nReps * bsd$nLocs + bsd$qcGenoCost
    costPerEntry[1] <- costPerEntry[1] + 
      bsd$candidateDevelCost + bsd$wholeGenomeCost
    stageBudgets <- targetBudget * percentages[-1]
    nEntries <- round(stageBudgets / costPerEntry)
    names(nEntries) <- bsd$stageNames
    return(list(nBreedingProg, nEntries))
  }#END numbersFromPercentages
  nBreedingProg <- numbersFromPercentages(percentages)
  nEntries <- nBreedingProg[[2]]
  nBreedingProg <- nBreedingProg[[1]]
  
  # Ensure that the population improvement cycle generates enough progeny
  failure1 <- nBreedingProg < bsd$minNBreedingProg
  # Ensure that enough variety candidates are in the final stage
  failure2 <- dplyr::last(nEntries) < bsd$nToMarketingDept
  # Ensure that budget ratios are respected
  ratios <- percentages[-length(percentages)] / percentages[-1]
  if (bsd$varietyType == "inbred") ratios[1] <- 1e9
  minRatios <- bsd$initBudget[grep("ratio", names(bsd$initBudget))]
  whichRatios <- which(minRatios - ratios > 1e-6)
  failure3 <- length(whichRatios) > 0
  if (failure1 | failure2 | failure3){
    # Decide what percentages to go toward
    if (exists("results", where=bsd)){
      percBest <- bsd$results %>% slice(which.max(fit)) %>% 
        dplyr::select(starts_with("perc")) %>% unlist
    } else{
      percBest <- bsd$initBudget[grep("perc", names(bsd$initBudget))]
    }
    # minPICbudget, minLastStgBudget
    lambda1 <- lambda2 <-lambda3 <- 0
    if (failure1){
      a <- percentages[1]
      lambda1 <- (bsd$initBudget["minPICbudget"] - a) / 
        (dplyr::first(percBest) - a)
    }
    if (failure2){
      b <- dplyr::last(percentages)
      lambda2 <- (bsd$initBudget["minLastStgBudget"] - b) / 
        (dplyr::last(percBest) - b)
    }
    if (failure3){
      for (r in whichRatios){
        p2c <- minRatios[r] * percentages[r+1]
        p2o <- minRatios[r] * percBest[r+1]
        l3 <- (p2c - percentages[r]) / 
          ((p2c - percentages[r]) - (p2o - percBest[r]))
        lambda3 <- max(l3, lambda3)
      }
    }
    lambda <- max(lambda1, lambda2, lambda3)
    percentages <- (1 - lambda) * percentages + lambda * percBest
  }
  nBreedingProg <- numbersFromPercentages(percentages)
  nEntries <- nBreedingProg[[2]]
  nBreedingProg <- nBreedingProg[[1]]
  
  # Ensure that nEntries gets smaller over stages
  failure <- any(diff(nEntries) > 0)
  
  if (!failure){
    bsd$nBreedingProg <- nBreedingProg
    bsd$nEntries <- nEntries
    bsd$realizedBudget <- calcBudget(bsd)
    bsd$percentages <- percentages
  }
  bsd$failure <- failure
  
  if (bsd$debug) on.exit()
  return(bsd)
}

#' makeGrid function
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
#' bsd <- makeGrid(bsd)
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
#' NOTE: the statements between startValues <- ... and endValues <- 
#' define the breeding scheme that is being optimized
#' 
#' @param percentages Double vector 1 + nStages with budget percentages for 
#' the population improvement cycle and the variety development pipeline
#' @param bsd List of breeding scheme data
#' @param returnBSD Logical if you want to return bsd because this is a one-off
#' rather than being called repeatedly for optimization
#'  
#' @return Vector with initial and end means and variances
#'
#' @details Call this function to test out the response of different budgets
#'
#' @examples
#' bsd <- runWithBudget(percentages, bsd)
#'
#' @export
runWithBudget <- function(percentages, bsd, returnBSD=F){
  if (bsd$debug){
    require(here)
    on.exit(expr={
      print(traceback())
      saveRDS(mget(ls()), file=here::here("data/runWithBudget.rds"))
    })
  }

  s <- Sys.time()
  percentages <- unlist(percentages)
  bsd <- budgetToScheme(percentages, bsd)
  percentages <- unlist(bsd$percentages)
  if (bsd$verbose) cat(percentages, bsd$failure, "\n")
  if (bsd$failure){
    if (bsd$verbose) print(Sys.time() - s)
    return(c(percentages, NA))
  }
  
  startValues <- calcCurrentStatus(bsd)
  
  for (twoPart in 1:bsd$nCyclesToRun){
    bsd$year <- bsd$year+1
    bsd <- makeVarietyCandidates(bsd)
    
    # Within the year, sequentially run the VDP...
    for (stage in 1:bsd$nStages){
      bsd <- chooseTrialEntries(bsd, toTrial=bsd$stageNames[stage],
                   fromTrial=ifelse(stage == 1, NULL, bsd$stageNames[stage-1]))
      bsd <- runVDPtrial(bsd, trialType=bsd$stageNames[stage])
    }

    # ... And then run the PIC. Sequencing could be made more complicated
    for (genSel in 1:bsd$nPopImpCycPerYear){
      optCont <- selectParents(bsd)
      bsd <- makeCrosses(bsd, optCont)
    }
  }
  
  endValues <- calcCurrentStatus(bsd)
  
  if (bsd$verbose) print(Sys.time() - s)
  if (bsd$debug) on.exit()
  if (returnBSD){
    return(bsd)
  } else{
    return(c(percentages, startValues, endValues, bsd$realizedBudget["budget"]))
  }
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
  if (bsd$debug){
    require(here)
    on.exit(expr={
      print(traceback())
      saveRDS(mget(ls()), file=here::here("data", "runBatch.rds"))
    })
  }

  if (bsd$debug){
    batchResults <- lapply(batchBudg, runWithBudget, bsd=bsd)
  } else{
    batchResults <- mclapply(batchBudg, runWithBudget, bsd=bsd, 
                             mc.cores=bsd$nCores)
  }
  # Remove results where budget was not valid
  batchResults <- batchResults[sapply(batchResults, 
                                      function(v) !any(is.na(v)))]
  rowNum <- length(batchResults[[1]])
  batchResults <- batchResults %>% unlist %>% matrix(nrow=rowNum) %>% t
  invisible(capture.output(
    batchResults <- as_tibble(batchResults, .name_repair="universal")
  ))
  # WARNING names here are hardcoded based on calcCurrentStatus
  colnames(batchResults) <- c(paste0("perc_", c("PIC", bsd$stageNames)),
                              paste0("init", c("PopMean", "PopSD", "VarMean")),
                              paste0("end", c("PopMean", "PopSD", "VarMean")), 
                              "Budget")
  batchResults <- batchResults %>% dplyr::mutate(response = endVarMean - initPopMean)
  
  if (bsd$debug) on.exit()
  return(batchResults)
}#END runBatch

#' findRedoBudgets function
#'
#' Based on simulations run so far, find the budgets that gave the highest
#' gain and the budgets with high gain but also high uncertainty
#' 
#' @param bsd List of breeding scheme data
#'  
#' @return List of two vectors one with high gain budgets
#' and one with uncertainty budgets
#'
#' @details Call this function to choose which budgets to try again
#'
#' @examples
#' bsd <- findRedoBudgets(bsd)
#'
#' @export
findRedoBudgets <- function(bsd){
  if (bsd$debug){  require(here)
    on.exit(expr={
      print(traceback())
      saveRDS(mget(ls()), file=here::here("data/findRedoBudgets.rds"))
    })
  }
  results <- bsd$results
  # Non-Parametric LOESS response
  loFormula <- paste0("response ~ ", 
                      paste0(colnames(results)[1:bsd$nStages], collapse=" + "))
  loFM <- loess(loFormula, data=results, degree=bsd$loessDegree)
  loPred <- predict(loFM, se=T)
  loPred <- tibble(fit=loPred$fit, se=loPred$se.fit, simNum=1:nrow(results))
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
    ndRows <- ndRows[order(ndRows$fit, decreasing=T)[1:bsd$nUncertain],]
  }
  whichUncert <- ndRows$simNum
  
  if (bsd$debug) on.exit()
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
#' @param bsd List of breeding scheme data. Has parameters for optimization too.
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
  if (bsd$debug){
    require(here)
    on.exit(expr={
      print(traceback())
      saveRDS(mget(ls()), file=here::here("data/optimizeByLOESS.rds"))
    })
  }
  
  nBatchesDone <- 0
  toleranceMet <- FALSE
  allPercRanges <- list()
  # 1. Make grid batch
  batch <- makeGrid(bsd)
  while (nBatchesDone < bsd$maxNumBatches & !toleranceMet){
    # 2. Run batch
    newBatchOut <- runBatch(batch, bsd)
    # 3. Evaluate stopping rule
    if (!exists("results", where=bsd)){
      bsd$results <- newBatchOut
    } else{
      bsd$results <- bsd$results %>% bind_rows(newBatchOut)
    }
    # Non-Parametric LOESS response
    loFormula <- paste0("response ~ ", 
                   paste0(colnames(bsd$results)[1:bsd$nStages], collapse=" + "))
    loFM <- loess(loFormula, data=bsd$results, degree=bsd$loessDegree)
    loPred <- predict(loFM, se=T)
    bsd$results <- bsd$results %>% dplyr::mutate(fit=loPred$fit, se=loPred$se.fit)
    # 4. Evaluate stopping rule
    bestGain <- which.max(bsd$results$fit)
    bestSE <- bsd$results$se[bestGain]
    bestClose <- which(bsd$results$fit[bestGain] - bsd$results$fit < 2*bestSE)
    percRanges <- bsd$results %>% slice(bestClose) %>% 
      dplyr::select(contains("perc")) %>% summarise_all(range)
    toleranceMet <- all(percRanges %>% summarise_all(diff) < bsd$tolerance)
    if (!toleranceMet){
      # 5. Find which budgets to rerun
      redoSims <- findRedoBudgets(bsd) %>% unlist
      batch <- lapply(redoSims, function(idx) 
        bsd$results %>% slice(idx) %>% dplyr::select(contains("perc")))
      # Function to combine two budgets staying closer to one
      combineBudg <- function(dummy, batch){
        b1 <- sample(length(batch), 2)
        b2 <- batch[[b1[2]]]
        b1 <- batch[[b1[1]]]
        anyNeg <- TRUE
        while(anyNeg){
          mix <- runif(1, 0.8, 1.1)
          combo <- mix * b1 + (1 - mix) * b2
          anyNeg <- any(combo <= 0)
        }
        return(combo)
      }#END combineBudg
      # batch <- c(batch, 
      #            lapply(1:(bsd$batchSize - length(batch)), combineBudg, batch))
      batch <- lapply(1:bsd$batchSize, combineBudg, batch)
    }#END !toleranceMet
    allPercRanges <- c(allPercRanges, list(list(percRanges, 
      nSimClose=length(bestClose), bestGain=bsd$results$fit[bestGain],
      bestSE=bestSE)))
    
    nBatchesDone <- nBatchesDone + 1
    
    if (bsd$saveIntermediateResults){
      # Save batches and results
      require(here)
      saveRDS(bsd$results, file=here::here("data", "allBatches.rds"))
      saveRDS(allPercRanges, file=here::here("data", "allPercentRanges.rds"))
    }
    
  }#END go through batches

  if (bsd$debug) on.exit()
  return(list(results=bsd$results, allPercentRanges=allPercRanges))
}
