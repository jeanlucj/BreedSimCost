
#' runVDPtrial function
#'
#' Function to run a variety development pipeline stage
#'
#' @param bsd The breeding scheme data list
#' @param trialType String name of the kind of trial to run. Must link back to
#' a trial type with parameters in bsd
#' @param entries A vector of entry IDs
#'
#' @return Updated bsd
#'
#' @details Add records to phenoRecords
#'
#' @examples
#' params <- runVDPtrial(bsd, trialType)
#'
#' @export
runVDPtrial <- function(bsd, trialType){
  entries <- bsd$entries
  # Calculate the error variance resulting from the number of locations and reps
  # Note that the number of years for a trial is always 1
  nL <- bsd$nLoc[trialType]; nR <- bsd$nRep[trialType]
  eV <- bsd$errVars[trialType]
  impliedVarE <- bsd$gxlVar/nL + bsd$gxyVar + bsd$gxyxlVar/nL + eV/nL/nR
  # Phenotypic evaluation of experimental lines
  pheno <- AlphaSimR::setPheno(bsd$varietyCandidates[entries],
                               varE=impliedVarE, simParam=bsd$SP)
  # Add the information to phenoRecords
  newRec <- tibble(year=bsd$year, trialID=bsd$nextTrialID, trialType=trialType,
                   id=pheno@id, pheno=pheno@pheno, selCrit=NA,
                   errVar=impliedVarE)
  bsd$phenoRecords <- bsd$phenoRecords %>% bind_rows(newRec)
  
  # Manage the trialID number
  bsd$nextTrialID <- bsd$nextTrialID + 1
  
  return(bsd)
}

#' chooseTrialEntries function
#'
#' Function to select which varieties to advance to the next stage assuming
#' they are IID (i.e., estimating genotypic values), by truncation selection
#'
#' @param bsd The breeding scheme data list
#' @param toTrial String the trial type the entries are going to 
#' @param fromTrial String the trial type the entries are coming from 
#' @return bsd updated with a vector of entry IDs in bsd$entries
#' @details Accesses all data in phenoRecords to pick the highest 
#' performing candidates.
#'
#' @examples
#' bsd <- chooseTrialEntries(bsd, toTrial, fromTrial)]
#'
#' @export
chooseTrialEntries <- function(bsd, toTrial, fromTrial=NULL){
  if (bsd$debug){
    require(here)
    on.exit(expr={
      print(traceback())
      saveRDS(mget(ls()), file=here::here("data", "chooseTrialEntries.rds"))
    })
  }
  
  if (toTrial == bsd$stageNames[1]){
    nInd <- nInd(bsd$varietyCandidates)
    nIndMax <- min(nInd, bsd$nEntries[1])
    entries <- bsd$varietyCandidates@id[nInd - (nIndMax - 1):0]
  } else{
    if (toTrial == "MarketingDept"){
      nToSelect <- bsd$nToMarketingDept
      fromYear <- 0
    } else{
      nToSelect <- bsd$nEntries[toTrial]
      fromYear <- -1
    }
    phenoRecords <- bsd$phenoRecords
    candidates <- phenoRecords %>% pull(id) %>% unique
    if (!is.null(fromTrial)){
      candidates <- phenoRecords %>% 
        dplyr::filter(trialType == fromTrial & year == bsd$year + fromYear) %>% 
        pull(id) %>% unique
      phenoRecords <- phenoRecords %>% dplyr::filter(id %in% candidates)
    }
    if (nrow(phenoRecords) > length(candidates)){ # There is some replication
      crit <- iidPhenoEval(phenoRecords)
    } else{
      crit <- phenoRecords %>% pull(pheno)
      names(crit) <- phenoRecords %>% pull(id)
    }
    # Add crit to phenoRecords
    for (n in names(crit)){
      bsd$phenoRecords$selCrit[bsd$phenoRecords$id == n] <- crit[n]
    }
    if (length(candidates) < nToSelect){
      # Workaround for a problem that should only happen in the transition from
      # burnin to two-part: the VDP nEntries between the two are mismatched
      print("There are too few variety candidates for trial")
      nCand <- nToSelect - length(candidates)
      stageNum <- which(bsd$stageNames == toTrial) # How many years back parents
      nBreedInd <- AlphaSimR::nInd(bsd$breedingPop)
      nIndBack <- bsd$initNBreedingProg * 
        (bsd$nPopImpCycPerYear * (stageNum - 1) + 1)
      nIndBack <- min(nIndBack, nBreedInd)
      breedPopIDs <- bsd$breedingPop@id[nBreedInd - nIndBack + 
                                          1:bsd$nBreedingProg]
      bsd <- makeVarietyCandidates(bsd, breedPopIDs, nCand) 
      nVarInd <- AlphaSimR::nInd(bsd$varietyCandidates)
      entries <- c(candidates, 
                   bsd$varietyCandidates@id[nVarInd - (nCand - 1):0])
    } else{
      crit <- crit[candidates]
      entries <- crit[(crit %>% order(decreasing=T))[1:nToSelect]] %>% names
    }
  }
  bsd$entries <- entries
  if (bsd$debug) on.exit()
  return(bsd)
}

#' makeVarietyCandidates function
#'
#' @param bsd List of breeding program data
#' @param breedPopIDs String vector IDs of breeding pop progenitors
#' of the variety candidates.  If NULL, the last nBreedingProg individuals used
#' @param nCandidates Make bespoke number of candidates if needed on the fly
#'
#' @return Updated bsd with new variety candidates in bsd$varietyCandidates
#' @details Here, creates DHs evenly distributed among the breeding progeny
#' in the last generation
#' @examples
#' bsd <- makeVarietyCandidates(bsd)]
#'
#' @export
makeVarietyCandidates <- function(bsd, breedPopIDs=NULL, nCandidates=NULL){
  if (is.null(nCandidates)) nCandidates <- bsd$nEntries[1]
  nInd <- nInd(bsd$breedingPop)
  switch(bsd$varietyType,
         clonal = {
           if (is.null(breedPopIDs)){ # Take the last breeding pop progeny
             if (nInd < nCandidates) stop("Not enough breeding pop progeny")
             if (bsd$nPopImpCycPerYear*bsd$nBreedingProg < nCandidates)
               cat("Previous-year Breeding pop progeny will be reused", bsd$nPopImpCycPerYear*bsd$nBreedingProg, nCandidates, "\n")
             newCand <- bsd$breedingPop[nInd - (nCandidates - 1):0]
           } else{
             newCand <- bsd$breedingPop[breedPopIDs]
           }
         }, 
         inbred = {
           if (is.null(breedPopIDs)){ # Take the last bsd$nBreedingProg
             nIndMax <- min(nInd, bsd$nBreedingProg)
             breedPopIDs <- bsd$breedingPop@id[nInd - (nIndMax - 1):0]
           }
           nInd <- length(breedPopIDs)
           nDH <- nCandidates %/% nInd
           if (nDH > 0){
             newCand <- makeDH(bsd$breedingPop[breedPopIDs], 
                               nDH=nDH, simParam=bsd$SP)
           }
           nExtra <- nCandidates %% nInd
           if (nExtra > 0){
             whichPar <- sample(nInd, nExtra)
             extraCand <- makeDH(bsd$breedingPop[breedPopIDs][whichPar], 
                                 nDH=1, simParam=bsd$SP)
             if (exists("newCand")){
               newCand <- c(newCand, extraCand)
             } else{
               newCand <- extraCand
             }
           }
         }#END inbred
  )#END switch
  
  # Update varietyCandidates population
  if (exists("varietyCandidates", bsd)){
    bsd$varietyCandidates <- c(bsd$varietyCandidates, newCand)
  } else{
    bsd$varietyCandidates <- newCand
  }
  
  return(bsd)
}

#' iidPhenoEval function
#'
#' Function to estimate variety candidate BLUPs assuming they are IID
#'
#' @param phenoRecords Tibble of phenotypic observations on variety candidates.
#' @return Named real vector of the BLUPs of all individuals in phenoRecords,
#' (names are the individual ids), with appropriate weights by error variance
#' @details Given all the phenotypic records calculate the best prediction of
#' the genotypic value for each individual using all its records
#'
#' @examples
#' iidBLUPs <- iidPhenoEval(phenoRecords)
#'
#' @export
iidPhenoEval <- function(phenoRecords){
  require(lme4)
  # Make error variances into weights
  phenoRecords <- phenoRecords %>% dplyr::mutate(wgt=1/errVar)
  fm <- lmer(pheno ~ (1|id), weights=wgt, data=phenoRecords)
  blup <- ranef(fm)[[1]]
  namesBlup <- rownames(blup)
  blup <- unlist(blup)
  # Ensure output has variation: needed for optimal contributions
  if (sd(blup) == 0){
    blup <- phenoRecords %>% group_by(id) %>%
      summarise(wm = weighted.mean(x=pheno, w=wgt)) %>% pull(wm)
  }
  names(blup) <- namesBlup
  return(blup)
}

#' calcVDPcovMat function
#' 
#' @param bsd List of breeding program data
#'
#' @return Covariance matrix of the true genotypic values of the variety
#' candidates with the phenotypic means from the VDP stages
#' 
#' @details This function came essentially from the selectiongain package 
#' developed by Melchinger's group: Mi et al. 2014 Crop Science. 
#' Simplified just to the "LonginII" method.
#'
#' @examples
#' vdpCovMat <- calcVDPcovMat(bsd)
#' 
#' @export
calcVDPcovMat <- function(bsd){
  dim = bsd$nStages + 1
  Vg <- bsd$genVar
  Vgy <- bsd$gxyVar
  Vgl <- bsd$gxlVar
  Vgyl <- bsd$gxyxlVar
  Ve <- bsd$errVars
  covMat <- diag(dim)*Vg
  covMat[1,] <- Vg
  covMat[,1] <- Vg
  for (i in 2:dim){
    covMat[i,i] <- Vg + Vgy +
      (Vgl+Vgyl)/bsd$nLocs[i-1] +
      Ve[i-1]/bsd$nLocs[i-1]/bsd$nReps[i-1]
  }
  for (i in 2:(dim-1)){
    for (j in (i+1):dim){
      covMat[i,j] <- Vg +
        Vgl/max(bsd$nLocs[i-1], bsd$nLocs[j-1])
      covMat[j,i] <- covMat[i,j]
    }
  }
  
  return(covMat)
}

#' multistageTrucPt function
#' 
#' @param alpha vector the selected fractions from one VDP stage to the next
#' @param corr matrix the correlation matrix from calcVDPcovMat
#' @param alg function the truncation point calculation algorithm from
#' package mvtnorm. Options are Miwa() and GenzBretz(). The latter is slower
#'
#' @return Covariance matrix of the true genotypic values of the variety
#' candidates with the phenotypic means from the VDP stages
#' 
#' @details This function came essentially from the selectiongain package 
#' developed by Melchinger's group: Mi et al. 2014 Crop Science. 
#' Simplified just to the "LonginII" method.
#'
#' @examples
#' vdpCovMat <- calcVDPcovMat(bsd)
#' 
#' @export
multistageTruncPt <- function(alpha,  corr=NA){
  if (any(is.na(corr))) stop("The corr matrix must not be NA")
  
  if (nrow(corr) == length(alpha)+1){
    corr <- corr[-1,-1]
  }
  else{
    stop ("The dimension of corr matrix must be one more than the length alpha")
  }
  
  nStg=length(alpha)
  alpha[alpha == 1] <- 0.9999
  
  lower <- rep(-Inf, nStg)
  for (stage in 1:nStg){
    truncPt <- tmvtnorm::qtmvnorm.marginal(p=1-alpha[stage], interval=c(-5,5), 
                                           n=stage, tail="lower.tail", 
                                           sigma=corr, lower=lower)
    lower[stage] <- truncPt$root
  }

  return(lower)
}

#' multistageGain function
#' 
#' @param corr matrix the correlation matrix from calcVDPcovMat
#' @param truncPts vector the truncation points that came from multistageTruncPt
#' @param Vg numeric the genotypic variance to provide a scale
#' 
#' @return numeric the expected gain from selection
#' 
#' @details This function came essentially from the selectiongain package 
#' but then I simplified it very much by just using a current R package
#'
#' @examples
#' vdpCovMat <- multistageGain(corr, truncPts)
#' 
#' @export
# mutistagecor.R author: xuefei mi, 11-03-2013, for selectiongain package v2.0.2
multistageGain <- function(corr, truncPts, Vg=1){
  require(tmvtnorm)
  # check if truncPts and corr have corresponding dimensions
  if (length(truncPts)!=nrow(corr)-1){
    stop("Dimension of truncPts must be same as nrow(corr)-1")
  }
  meanGenVal <- tmvtnorm::mtmvnorm(sigma=corr, lower=c(-Inf, truncPts))$tmean[1]
  return(meanGenVal*sqrt(Vg))
}

# Order of operations
# 0. From the full optimization, get the VDP budget
# 0.1 Calculate the implied covariance matrix in the VDP: calcVDPcovMat
# 1. With the VDP budget, calculate the min max entries: calcMinMaxEntries
# 2. With the min max entries, make the grid: makeNEntryGrid
# 3. Calculate the gains across the grid: calcVDPgain
# 4. Determine if the max is balanced more toward the end than the full two-part

#' makeNEntryGrid function
#'
#' @param nIndMin nStages vector minimum number of individuals at each stage
#' @param nIndMax nStages vector maximum number of individuals at each stage
#' @param step nStages vector: interval of number of individuals. If NULL will
#' increment for last stage and have the same number of steps for other stages
#' @param maxBudget the maximum cost of all evaluations
#' @param costProd per individual, how much does it cost to produce the seed?
#' @param costTest per individual, how much does it cost to run all the plots?
#'
#' @return List with feasible entry number schemes
#'
#' @details Call this function to set up optimization
#'
#' @examples
#' bsd <- makeNEntryGrid(bsd)
#'
#' @export
makeNEntryGrid <- function(nIndMin, nIndMax, step=NULL, 
                           minBudget=0.9*maxBudget, maxBudget, indCost){
  nStages <- length(nIndMin)
  if (is.null(step)){
    nSteps <- nIndMax[nStages] - nIndMin[nStages]
    step <- (nIndMax - nIndMin) / nSteps
  }
  nEntryList <- as.list(seq(from=nIndMin[1], to=nIndMax[1], by=step[1]) %>% round)
  for (stg in 2:nStages){
    nEntryNew <- seq(from=nIndMin[stg], to=nIndMax[stg], by=step[stg]) %>% round
    nEntryList <- mapply(c, rep(nEntryList, each=length(nEntryNew)),
                         rep(nEntryNew, length(nEntryList)), SIMPLIFY=F)
  }
  
  # Make sure that nEntries gets smaller going forward
  nIndDecreases <- function(nEntryVec){
    all(diff(nEntryVec) <= 0)
  }
  nEntryList <- nEntryList[sapply(nEntryList, nIndDecreases)]
  
  # Calculate budgets to eliminate ones that are too small or too big
  calcBudget <- function(nEntryVec, indCost){
    return(nEntryVec %*% indCost)
  }
  budgets <- sapply(nEntryList, calcBudget, indCost=indCost)
  
  return(nEntryList[budgets > minBudget & budgets < maxBudget])
}#END makeNEntryGrid

#' calcVDPGain function
#'
#' @param nEntryVec nStages vector with number of individuals at each stage
#' @param nFinal numeric how many individuals are going to the marketing dept
#' @param cov matrix covariances as came out of calcVDPcovMat
#'
#' @return numeric expected gain from the VDP
#'
#' @details Call this function on a grid of budget equal VDPs. Wrapper function 
#' to call the truncation point and gain functions with elements from the grid
#'
#' @examples
#' gains <- sapply(testGrid, calcVDPGain, 
#'                 nFinal=bsd$nToMarketingDept, cov=vdpCovMat)
#'
#' @export
calcVDPGain <- function(nEntryVec, nFinal, cov){
  nEntryVec <- c(nEntryVec, nFinal)
  genVar <- cov[1, 1]
  corr <- cov2cor(cov)
  selFrac <- nEntryVec[-1] / nEntryVec[-length(nEntryVec)]
  truncPts <- multistageTruncPt(selFrac, corr)
  return(multistageGain(corr=corr, truncPts=truncPts, Vg=genVar))
}

#' calcMinMaxEntries function
#'
#' @param bsd List of breeding program data
#' @param budgetVec vector as comes out of calcBudget. If NULL look in bsd
#' @param includePIC logical should a min and max for nProg in PIC be figured?
#'
#' @return List with vectors for minimum entries and maximum entries and cost
#' per individual entry
#'
#' @details Assume the location costs are fixed. 
#'
#' @examples
#' minMaxList <- calcMinMaxEntries(bsd)
#'
#' @export
calcMinMaxEntries <- function(bsd, budgetVec=NULL, includePIC=F){
  if (is.null(budgetVec)){
    if (exists("realizedBudget", bsd)){
      budgetVec <- bsd$realizedBudget
    }
    else{
      budgetVec <- bsd$initBudget
    }
  }
  indCost <- bsd$plotCosts * bsd$nReps * bsd$nLocs + bsd$qcGenoCost
  indCost[1] <- indCost[1] + 
    bsd$candidateDevelCost + bsd$wholeGenomeCost
  varBudget <- budgetVec["budget"] - budgetVec["locationCosts"]
  if (includePIC){
    costPerPIC <- bsd$nPopImpCycPerYear * 
      (bsd$crossingCost + bsd$wholeGenomeCost)
    indCost <- c(costPerPIC, indCost)
  }
  else{
    varBudget <- varBudget - budgetVec["picBudget"]
  }
  
  minLast <- bsd$nToMarketingDept
  minEntries <- rep(minLast, length(indCost))
  minStgCosts <- indCost * minEntries
  maxEntries <- NULL
  for (stg in 1:length(indCost)){
    costRest <- sum(minStgCosts[-(1:stg)])
    perEntToStg <- sum(indCost[1:stg])
    maxEntries <- c(maxEntries, floor((varBudget - costRest)/perEntToStg))
  }
  minEntries[1] <- maxEntries[stg]
  return(list(minEntries=minEntries, maxEntries=maxEntries, 
              indCost=indCost, varBudget=varBudget))
}
