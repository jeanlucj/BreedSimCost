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
#' Rules: 
#' 1. For PIC, nBreedingProg has to be at least minNBreedingProg
#' 2. Each stage has to have fewer entries than the previous stage
#' 3. At the end of the VDP, there have to be more entries than go to marketing
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
#' bsd <- budgetToScheme(percentages=percentages, bsd)
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
  }#END if failure
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
makeGrid <- function(bsd, justVDP=F){
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
  
  popParmsByCyc <- calcCurrentStatus(bsd)
  
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
    
    popParmsByCyc <- popParmsByCyc %>% bind_rows(calcCurrentStatus(bsd))
  }

  if (bsd$verbose) print(Sys.time() - s)
  if (bsd$debug) on.exit()
  if (returnBSD){
    return(bsd)
  } else{
    totalGain <- (popParmsByCyc %>% slice_tail(n=1))$varCandMean -
      (popParmsByCyc %>% slice_head(n=1))$breedPopMean
    return(list(gainBudgetPerc=c(totalGain=totalGain, 
                                 budget=bsd$realizedBudget["budget"], 
                                 percentages), 
                popParmsByCyc=popParmsByCyc))
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
#' batchResults <- runBatch(batchBudg, bsd)
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
                             mc.preschedule=F, mc.cores=bsd$nCores)
  }
  # batchResults is now a list of lists
  # Remove results where budget was not valid
  findBatchNA <- function(oneRunRes){
    !any(is.na(unlist(oneRunRes)))
  }
  batchResults <- batchResults[sapply(batchResults, 
                                      function(v) findBatchNA(v))]

  if (bsd$debug) on.exit()
  return(batchResults)
}#END runBatch
