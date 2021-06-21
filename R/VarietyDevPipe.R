
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
    nToSelect <- bsd$nEntries[toTrial]
    phenoRecords <- bsd$phenoRecords
    candidates <- phenoRecords %>% pull(id) %>% unique
    if (!is.null(fromTrial)){
      candidates <- phenoRecords %>% 
        dplyr::filter(trialType == fromTrial & year == bsd$year - 1) %>% 
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
  swtich(bsd$varietyType,
         clonal = {
           if (is.null(breedPopIDs)){ # Take the last breeding pop progeny
             if (nInd < nCandidates) stop("Not enough breeding pop progeny")
             if (bsd$nPopImpCycPerYear*bsd$nBreedingProg < nCandidates)
               print("Previous-year Breeding pop progeny will be reused")
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
  names(blup) <- namesBlup
  # Ensure output has variation: needed for optimal contributions
  if (sd(blup) == 0){
    blup <- tapply(phenoRecords$pheno, phenoRecords$id, mean)
  }
  return(blup)
}
