
#' selectParents function
#'
#' Function to select which breeding population individuals to use as parents
#' to advance the population
#'
#' @param bsd The breeding scheme data list
#' @return Data.frame with pedigree ccMat-based optimal contributions
#' @details Accesses phenoRecords to predict GEBVs among breeding population
#' then calculates optimal contributions with those GEBVs and a pedigree-based
#' coefficient of coancestry matrix.
#'
#' @examples
#' optCont <- selectParents(bsd)]
#'
#' @export
selectParents <- function(bsd){
  # Get GEBVs
  gebv <- grmPhenoEval(bsd)
  # Keep only breeding population
  nInd <- nInd(bsd$breedingPop)
  nIndMax <- min(nInd, bsd$nBreedingProg*bsd$keepNBreedingCyc)
  stopAt <- bsd$minParentAge*bsd$nPopImpCycPerYear
  selCan <- bsd$breedingPop@id[nInd - (stopAt + (nIndMax - 1):0)]
  gebv <- gebv[selCan]
  # Set up to use optiSel
  # 1. data.frame with a column named "Indiv" that has individual ids,
  # and a column named however you want with the selection criterion.
  # Here, the column is named "gebv".
  phen <- data.frame(Indiv=names(gebv), gebv=gebv)
  # 2. Need to have a coefficient of coancestry matrix
  pedigree <- data.frame(id=bsd$breedingPop@id, dam=bsd$breedingPop@mother, 
                         sire=bsd$breedingPop@father)
  pedigree <- convertNamesToRows(pedigree)
  ccMat <- pedigreeToCCmatrix(pedigree)
  rownames(ccMat) <- colnames(ccMat) <- bsd$breedingPop@id
  ccMat <- ccMat[selCan, selCan]
  # 3. the function `candes` checks that all is well.
  invisible(capture.output(cand <- optiSel::candes(phen, ccMat=ccMat, quiet=T)))
  # 4. The constraint. Here, I give it a change of inbreeding using the
  # optiContEffPop parameter in bsd
  # The constraint should be `ub.name`, where `name` is the name of the
  # cc matrix given above
  deltaF <- 1 / (2 * bsd$optiContEffPop)
  con <- list(
    ub.ccMat = deltaF + (1 - deltaF)*cand$mean$ccMat
  )
  # 5. `oc` has the optimal contributions for each individual. A number of them
  # will be close to zero and should be explicitly set to zero (see keep below)
  oc <- optiSel::opticont("max.gebv", cand, con, quiet=T, trace=F)
  oc <- oc$parent[, c("Indiv", "oc")]
  # Keep lines to that have a chance to be a parent once
  keep <- which(oc$oc > 1 / (2*bsd$nBreedingProg))
  oc <- oc[keep,]
  oc$oc <- oc$oc / sum(oc$oc)
  return(oc)
}

#' makeCrosses function
#'
#' Function to make crosses to get the next generation of breeding progeny
#' @param bsd List of breeding program data
#' @param optCont Data.frame Optimal contributions coming from selectParents
#'
#' @return Updated bsd with new S0 progeny in breedingPop
#' @details Takes the optimal contributions and generates random mating with
#' the correct number of progeny per parent. Self-fertilization is possible.
#' @examples
#' bsd <- makeCrosses(bsd, optCont)]
#'
#' @export
makeCrosses <- function(bsd, optCont){
  nProgeny <- bsd$nBreedingProg
  nParents <- nrow(optCont)
  probs <- optCont$oc
  ids <- optCont$Indiv
  if (abs(1 - sum(probs)) > 1e-9) stop("Crossing probabilites must sum to 1")
  # Two parents per progeny
  parIdx <- round(cumsum(2*nProgeny*probs))
  nEachPar <- diff(c(0, parIdx))
  parVec <- rep(ids, nEachPar)
  crossPlan <- matrix(sample(parVec), ncol=2)
  newBreedProg <- AlphaSimR::makeCross(bsd$breedingPop, crossPlan, 
                                       simParam=bsd$SP)
  bsd$breedingPop <- c(bsd$breedingPop, newBreedProg)
  return(bsd)
}

#' pedigreeToCCmatrix function
#'
#' Function to calculate a coefficient of coancestry matrix from a three-
#' column pedigree matrix:
#' The first column has to be the row number
#' Sire and dam columns refer directly to rows
#' Unknown parents need to be set to 0 (zero)
#'
#' @param threeColPed The pedigree matrix
#' @return Matrix with coefficients of coancestry
#' @details Uses the breedingPop data in bsd
#'
#' @examples
#' ccMat <- pedigreeToCCmatrix(threeColPed)]
#'
#' @export
pedigreeToCCmatrix <- function(threeColPed){
  nInd <- nrow(threeColPed)
  ccMat <- matrix(0, nInd, nInd)
  # the very first individual in the pedigree has to be a founder
  ccMat[1, 1] <- 0.5
  for (prog in 2:nInd){
    sire <- threeColPed[prog, 2]
    dam <- threeColPed[prog, 3]
    prog_1 <- prog - 1
    if (sire){
      sireRow <- ccMat[sire, 1:prog_1]
    } else{
      sireRow <- rep(0, prog_1)
    }
    if (dam){
      damRow <- ccMat[dam, 1:prog_1]
    } else{
      damRow <- rep(0, prog_1)
    }
    ccMat[prog, 1:prog_1] <- ccMat[1:prog_1, prog] <- (sireRow + damRow) / 2
    ccSelf <- 0.5
    if (sire > 0 & dam > 0) ccSelf <- ccSelf + ccMat[sire, dam] / 2
    ccMat[prog, prog] <- ccSelf
  }
  rownames(ccMat) <- colnames(ccMat) <- 1:nInd
  return(ccMat)
}

#' convertNamesToRows function
#'
#' Function to transform arbitrary name IDs to row numbers for the sake of
#' pedigreeToCCmatrix
#'
#' @param nameMat The pedigree matrix with character IDs
#' @return Matrix suitable for pedigreeToCCmatrix
#' @details Keeps everything in order
#'
#' @examples
#' threeColPed <- convertNamesToRows(nameMat)]
#'
#' @export
convertNamesToRows <- function(nameMat){
  nameToRow <- 1:nrow(nameMat)
  names(nameToRow) <- nameMat[,1]
  parVecToRow <- function(parVec){
    rowID <- integer(length(parVec))
    parVec[!(parVec %in% nameMat[,1])] <- "0"
    rowID[parVec != "0"] <- nameToRow[parVec[parVec != "0"]]
    return(rowID)
  }
  return(cbind(nameToRow, parVecToRow(nameMat[,2]), parVecToRow(nameMat[,3])))
}

#' calcGRM function
#'
#' Function to make a genomic relationship matrix to calculate GEBVs
#'
#' @param bsd List of breeding scheme data
#' @return A genomic relationship matrix
#' @details bsd contains both variety candidates that have phenotypes and
#' breeding population individuals that don't.  Both need to be in the GRM
#' for prediction
#'
#' @examples
#' grm <- calcGRM(bsd)
#'
#' @export
calcGRM <- function(bsd){
  require(sommer)
  nInd <- nInd(bsd$varietyCandidates)
  nIndMax <- min(nInd, bsd$nEntries[1]*bsd$keepNTrainingCyc)
  allVar <- bsd$varietyCandidates[nInd - (nIndMax - 1):0]
  nInd <- nInd(bsd$breedingPop)
  nIndMax <- min(nInd, bsd$nBreedingProg*bsd$keepNBreedingCyc)
  stopAt <- bsd$minParentAge*bsd$nPopImpCycPerYear
  allPop <- bsd$breedingPop[nInd - (stopAt + (nIndMax - 1):0)]
  switch(bsd$varietyType,
         clonal = {
           varNotPop <- setdiff(allVar@id, allPop@id)
           if (length(varNotPop) > 0) allPop <- c(allVar[varNotPop], allPop)
         },
         allPop <- c(allVar, allPop)
         )

  return(sommer::A.mat(AlphaSimR::pullSnpGeno(allPop, simParam=bsd$SP) - 1))
}

#' grmPhenoEval function
#'
#' Takes phenotypes from phenoRecords and combines them with the GRM across
#' variety candidates and breeding population individuals
#'
#' @param bsd List of breeding scheme data
#' @return Named real vector of the GEBVs (names are the individual ids) of all
#' individuals including breeding population and variety candidates,
#' with appropriate weights by error variance of the observations
#' @details Given all the phenotypic records calculate the GEBV
#'
#' @examples
#' grmBLUPs <- grmPhenoEval(bsd)
#'
#' @export
grmPhenoEval <- function(bsd){
  require(sommer)
  grm <- calcGRM(bsd)
  phenoDF <- bsd$phenoRecords
  phenoDF <- phenoDF %>% filter(id %in% rownames(grm)) %>% 
    dplyr::mutate(wgt=1/errVar)
  phenoDF$id <- factor(phenoDF$id, levels=rownames(grm)) # Enable prediction
  fm <- sommer::mmer(pheno ~ 1,
             random= ~ vs(id, Gu=grm),
             method="EMMA",
             rcov= ~ units,
             weights=wgt,
             data=phenoDF,
             verbose=F,
             date.warning=F)
  blup <- fm$U[[1]][[1]]
  # Ensure output has variation: needed for optimal contributions
  if (sd(blup) == 0){ # In this case, random selection
    nb <- names(blup)
    blup <- runif(length(blup))
    names(blup) <- nb
  }
  return(blup)
}
