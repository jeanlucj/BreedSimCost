#' initializeProgram function
#'
#' Function to initialize the simulation.
#' Read the founder genetic architecture and population genetic parameters.
#' Create simulation data structures.
#'
#' @param founderFile String name of the text file with parameters governing
#' the founders and the genetic architecture
#' @param schemeFile String name of the text file with parameters governing
#' the variety development pipeline and selection functions
#' @return Named list containing the `AlphaSimR` SP file, breedingPop the
#' breeding population, varietyCand an empty pop for variety candidates,
#' phenoRecords an empty tibble for phenotype records
#'
#' @details Call this function at the beginning of the simulation
#'
#' @examples
#' bsd <- initializeProgram("FounderCtrlFile.txt", "SchemeCtrlFile.txt", 
#'                          "CostsCtrlFile.txt", "OptimizationCtrlFile.txt")
#'
#' @export
initializeProgram <- function(founderFile, schemeFile, 
                              costFile, optimizationFile){
  # FounderFile
  # Read parameters to create founders
  parmNames <- c("varietyType", "nChr", "effPopSize", "quickHaplo", 
                 "segSites", "nQTL", "nSNP", "genVar", "gxeVar", 
                 "gxyVar", "gxlVar", "gxyxlVar", "meanDD", "varDD", 
                 "relAA")
  bsd <- readControlFile(founderFile, parmNames)

  # SchemeFile
  # Read parameters about the overall scheme
  parmNames <- c("nCyclesToRun", "nBurnInCycles", 
                 "nStages", "stageNames", "nEntries", "nToMarketingDept", 
                 "nReps", "nLocs", "errVars", 
                 "optiContEffPop", "nBreedingProg", "nPopImpCycPerYear", 
                 "keepNTrainingCyc", "keepNBreedingCyc", "minParentAge")
  bsdNew <- readControlFile(schemeFile, parmNames)
  bsd <- c(bsd, bsdNew)
  # Because different budget allocations will change these, store initial values
  bsd$initNEntries <- bsd$nEntries
  bsd$initNBreedingProg <- bsd$nBreedingProg
  
  # CostFile
  # Read parameters about scheme costs
  parmNames <- c("plotCosts", "perLocationCost", "crossingCost", 
                 "candidateDevelCost", "qcGenoCost", "wholeGenomeCost")
  bsdNew <- readControlFile(costFile, parmNames)
  bsd <- c(bsd, bsdNew)
  
  # OptimizationFile
  # Read parameters about optimization procedure
  parmNames <- c("nCores", "minPercentage", "maxPercentage",
                 "percentageStep", "minNBreedingProg",
                 "debug", "verbose")
  bsdNew <- readControlFile(optimizationFile, parmNames)
  bsd <- c(bsd, bsdNew)
  
  # Add miscelaneous data structures
  bsd$year <- 0
  bsd$nextTrialID <- 1
  bsd <- calcDerivedParms(bsd)
  bsd$initBudget <- calcBudget(bsd)
  
  # Create haplotypes for founder population of outbred individuals
  if (bsd$quickHaplo){
    founderHap <- AlphaSimR::quickHaplo(nInd=bsd$nBreedingProg, 
                             nChr=bsd$nChr, segSites=bsd$segSites)
  } else{
    founderHap <- AlphaSimR::runMacs2(nInd=bsd$nBreedingProg, 
                           nChr=bsd$nChr, segSites=bsd$segSites, 
                           Ne=bsd$effPopSize)
  }
  
  # Global simulation parameters from founder haplotypes
  bsd$SP <- AlphaSimR::SimParam$new(founderHap)
  bsd$SP$restrSegSites(minQtlPerChr=1, minSnpPerChr=10, overlap=FALSE)
  # Additive, dominance, and epistatic trait architecture
  bsd$SP$addTraitADE(nQtlPerChr=bsd$nQTL, var=bsd$genVar, meanDD=bsd$meanDD, varDD=bsd$varDD, relAA=bsd$relAA, useVarA=FALSE)
  # Observed SNPs per chromosome
  bsd$SP$addSnpChip(bsd$nSNP)
  
  # Create the founders
  bsd$breedingPop <- AlphaSimR::newPop(founderHap, simParam=bsd$SP)

  return(bsd)
}

#' selectParentsBurnIn function
#'
#' Function to select varieties to cross during phenotypic selection burnin
#' phase, from variety candidates, by truncation selection
#'
#' @param bsd The breeding scheme data list
#' @return Vector of entry IDs
#' @details Accesses data in phenoRecords to pick the highest among candidates
#'
#' @examples
#' entries <- selectParentsBurnIn(bsd)]
#'
#' @export
selectParentsBurnIn <- function(bsd){
  nToSelect <- min(bsd$nEntries)
  candidates <- bsd$phenoRecords$id %>% unique
  if (length(candidates) < nToSelect){
    stop("There are too few variety candidates to be parents")
  }
  if (nrow(bsd$phenoRecords) > length(candidates)){ # There is some replication
    crit <- iidPhenoEval(bsd$phenoRecords)
  } else{
    crit <- c(bsd$phenoRecords$pheno)
    names(crit) <- bsd$phenoRecords$id
  }
  parents <- crit[(crit %>% order(decreasing=T))[1:nToSelect]] %>% names
  return(parents)
}

#' makeCrossesBurnIn function
#'
#' Function to make crosses during the phenotypic selection burnin
#' @param bsd List of breeding program data
#' @param parents String vector of ids of parents
#'
#' @return Updated bsd with new S0 progeny in breedingPop
#' @details Randomly mates the parents. Self-fertilization is possible.
#' @examples
#' bsd <- makeCrossesBurnIn(bsd, parents)]
#'
#' @export
makeCrossesBurnIn <- function(bsd, parents){
  nProgeny <- max(bsd$nEntries)
  nParents <- length(parents)
  # Two parents per progeny
  parVec <- rep(parents, length.out=nProgeny*2)
  crossPlan <- matrix(sample(parVec), ncol=2)
  newBreedProg <- AlphaSimR::makeCross(bsd$varietyCandidates, crossPlan, 
                                       simParam=bsd$SP)
  bsd$breedingPop <- c(bsd$breedingPop, newBreedProg)
  return(bsd)
}

#' readControlFile function
#'
#' Function to read a text control file
#'
#' The text file should be organized as follows
#' 1. Any text after a comment symbol # will be ignored
#' 2. Control parameter names should be on their own line
#' 3. Parameter values should be on the following line. If multiple parameter values are needed they should be separated by white space but on the same line
#' @param fileName The name of the text file to be read. Must include the path to the file
#' @param parmNames A string vector with the names of the control parameters that will be searched in the text file
#' @return A named list of the parameter values read from the control file
#'
#' @details Call this function before beginning the simulation
#'
#' @examples
#' params <- readControlFile("./inputDir/ctrlFile.txt", c("nStages", "nParents", "nCrosses"))
#'
#' @export
readControlFile <- function(fileName, parmNames){
  ctrlLines <- readLines(fileName)
  ctrlLines <- sapply(ctrlLines, function(st) strsplit(st, "#", fixed=T)[[1]][1])
  getParm <- function(parmName){
    parmRow <- grep(parmName, ctrlLines)+1
    parms <- unlist(strsplit(ctrlLines[parmRow], "[[:space:]]"))
    names(parms) <- NULL
    parmsNum <- suppressWarnings(as.numeric(parms))
    if (!any(is.na(parmsNum))) parms <- parmsNum
    return(parms)
  }
  parms <- lapply(parmNames, getParm)
  names(parms) <- parmNames
  return(parms)
}

#' calcDerivedParms function
#'
#' Once you have read in parameters from a control file, or set them yourself, there are still a few derived parameters that are needed.  This function calculates them.
#'
#' @param bsd List of breeding program data
#' @return A list bsd that extends the input with a few derived parameters
#'
#' @details This function is only called internally by other functions used to specify the pipeline
#'
calcDerivedParms <- function(bsd){
  # Function to check if a parameter has no value
  nv <- function(parm){
    is.null(parm) | length(parm) == 0
  }

  # Some parms have to be logical
  makeLogical <- function(parm){
    if (nv(parm)) parm <- FALSE else parm <- as.logical(parm)
    if (is.na(parm)) parm <- FALSE
    return(parm)
  }

  # Prevent errors having to do with inconsistent parameters
  if (bsd$nSNP + bsd$nQTL >= bsd$segSites){
    print("The number of segregating sites (segSites) has to be greater than the number of SNPs (nSNP) and the number of QTL (nQTL). segSites set 10% bigger than nSNP + nQTL")
    bsd$segSites <- ceiling((bsd$nSNP + bsd$nQTL) * 1.1)
  }

  bsd$quickHaplo <- makeLogical(bsd$quickHaplo)
  bsd$debug <- makeLogical(bsd$debug)
  bsd$verbose <- makeLogical(bsd$verbose)

  # Genetic architecture defaults
  if (nv(bsd$meanDD)) bsd$meanDD <- 0
  if (nv(bsd$varDD)) bsd$varDD <- 0
  if (nv(bsd$relAA)) bsd$relAA <- 0

  # Check that these vectors are of the right length
  rightLength <- function(objName) length(get(objName, bsd)) == bsd$nStages
  vecNames <- c("stageNames", "nEntries", "nReps", "nLocs", "errVars")
  rl <- sapply(vecNames, rightLength)
  if (any(!rl)){
    stop(paste("These vectors do not have the right length:",
               paste(vecNames[!rl], collapse=" ")))
  }

  # Defaults for GxE variance
  if (any(nv(bsd$gxyVar), nv(bsd$gxlVar), nv(bsd$gxyxlVar))){
    if (!nv(bsd$gxeVar)){
      if (nv(bsd$gxyVar)) bsd$gxyVar <- bsd$gxeVar / 3
      if (nv(bsd$gxlVar)) bsd$gxlVar <- bsd$gxeVar / 3
      if (nv(bsd$gxyxlVar)) bsd$gxyxlVar <- bsd$gxeVar / 3
    } else{
      if (nv(bsd$gxyVar)) bsd$gxyVar <- 0
      if (nv(bsd$gxlVar)) bsd$gxlVar <- 0
      if (nv(bsd$gxyxlVar)) bsd$gxyxlVar <- 0
    }
  }

  # Make sure everything has names
  names(bsd$nEntries) <- names(bsd$nReps) <- names(bsd$nLocs) <- 
    names(bsd$errVars) <- bsd$stageNames

  return(bsd)
}

#' calcCurrentStatus function
#'
#' Calculate mean and variance statistics on your current breeding program.
#'
#' @param bsd List of breeding program data
#' @return A vector extracting breeding population means and variances:
#' breedPopMean, breedPopSD, varCandMean
#' 
#' @details This function is only called internally by other functions used to specify the pipeline
#'
#' @export
calcCurrentStatus <- function(bsd){
  # Mean of the last breeding population
  nInd <- nInd(bsd$breedingPop)
  currBreedPop <- bsd$breedingPop[nInd - (bsd$nBreedingProg - 1):0]
  breedPopMean <- mean(gv(currBreedPop))
  breedPopSD <- sd(gv(currBreedPop))
  breedPopAddSD <- sd(bv(currBreedPop, simParam=bsd$SP))
  
  # Mean of the variety candidates being sent to the marketing department
  bsd <- chooseTrialEntries(bsd, toTrial="MarketingDept", 
                                fromTrial=bsd$stageNames[bsd$nStages])
  varCandMean <- mean(gv(bsd$varietyCandidates[bsd$entries]))
  return(tibble(breedPopMean=breedPopMean, breedPopSD=breedPopSD, 
           breedPopAddSD=breedPopAddSD, varCandMean=varCandMean))
}

#' loessPredCount function
#'
#' Summarize data from simulations in a form to make a nice hexbin plot
#' Return both counts in bins and bin means
#'
#' @param resultMat Tibble the simulation results
#' @param nSim Integer Predictions and counts for the first nSim simulations
#' @param xlim Real max and min for plot
#' @param ylim Real max and min for plot
#' @param budg1 Integer which column of percent budget to make plot for x-axis
#' @param budg2 Integer which column of percent budget to make plot for y-axis
#' @param degree Interger the degree of the polynomial LOESS fits
#' @return A list with bin counts and means and results from LOESS predictions
#'
#' @details Use the hexbin function to tabulate 2D bins in hexagonal array
#'
#' @export
loessPredCount <- function(resultMat, nSim=nrow(resultMat), 
                      xlim=NULL, ylim=NULL, 
                      budg1=1, budg2=2, degree=1){
  require(hexbin)
  if (is.null(xlim)){
    xlim <- c(floor(min(resultMat[,budg1])*50-0.5)/50, 
              ceiling(max(resultMat[,budg1])*50+0.5)/50)
  }
  if (is.null(ylim)){
    ylim <- c(floor(min(resultMat[,budg2])*50-0.5)/50, 
              ceiling(max(resultMat[,budg2])*50+0.5)/50)
  }
  
  uptoResults <- resultMat[1:nSim,]
  # Non-Parametric LOESS response
  predictors <- resultMat %>% colnames %>% stringr::str_subset("perc")
  predictors <- predictors[-length(predictors)] # remove last predictor
  loFormula <- paste0("response ~ ", paste0(predictors, collapse=" + "))
  loFM <- stats::loess(loFormula, data=uptoResults, degree=degree)
  loPred <- predict(loFM, se=T)
  
  # Find budget with highest predicted gain
  percHiPredGain <- uptoResults %>% dplyr::slice(which.max(loPred$fit)) %>% 
    dplyr::select(contains("perc"))
  
  # Make the bins!
  bins <- hexbin::hexbin(cbind(uptoResults[,budg1], uptoResults[,budg2]), 
            xbnds=xlim, ybnds=ylim, 
            xlab=colnames(resultMat)[budg1], ylab=colnames(resultMat)[budg2])
  
  # Mean gain for each hexagon
  calcDistToCell <- function(hex){
    xcm <- bins@xcm[hex]; ycm <- bins@ycm[hex]
    sqrt((uptoResults[,budg1] - xcm)^2 + (uptoResults[,budg2] - ycm)^2)
  }
  allDist <- lapply(1:bins@ncells, calcDistToCell)
  allDist <- matrix(unlist(allDist), nSim)
  closeBin <- apply(allDist, 1, which.min)
  # If bin empty steal observation closest to that bin 
  # Not sure why this happens
  occBins <- sort(unique(closeBin))
  emptyBins <- setdiff(1:max(occBins), occBins)
  wme <- apply(allDist[,emptyBins, drop=F], 2, which.min)
  closeBin[wme] <- emptyBins
  calcBinMean <- function(bin){
    mean(loFM$fitted[closeBin == bin])
  }
  binMean <- sapply(1:bins@ncells, calcBinMean)
  
  highMean <- which.max(binMean)
  hiMeanXY <- c(bins@xcm[highMean], bins@ycm[highMean])
  highCt <- which.max(bins@count)
  hiCtXY <- c(bins@xcm[highCt], bins@ycm[highCt])
  highPred <- closeBin[which.max(loPred$fit)]
  hiPredXY <- c(bins@xcm[highPred], bins@ycm[highPred])
  
  percMean <- uptoResults %>% dplyr::filter(closeBin == highMean) %>% 
    dplyr::select(contains("perc")) %>% dplyr::summarize_all(mean)

  return(list(loPred=loPred, bins=bins, closeBin=closeBin, binMean=binMean, 
              hiMeanXY=hiMeanXY, hiCtXY=hiCtXY, hiPredXY=hiPredXY,
              percMean=percMean, percHiPredGain=percHiPredGain))
}#END loessPredCount

#' plotLoessPred function
#'
#' Make a nice hexbin plot with the summary from loessPredCount
#'
#' @param resultMat Tibble the simulation results
#' @param nSim Integer Predictions and counts for the first nSim simulations
#' @param xlim Real max and min for plot
#' @param ylim Real max and min for plot
#' @param budg1 Integer which column of percent budget to make plot for x-axis
#' @param budg2 Integer which column of percent budget to make plot for y-axis
#' @param binMeanContrast Numeric a higher value makes the peak stand out more
#' @param plotMn Logical if TRUE, plot mean of gain for bin not bin counts
#' @param plotHiMn Logical whether to put a red circle at the hexbin with the
#' highest mean
#' @param plotHiCt Logical whether to put a green circle where
#' the most simulations were run
#' @param plotHiPredGain Logical whether to put a blue circle where on the point
#' with the highest LOESS predicted gain
#' @param giveRange Logical whether to put the range of gains in the main title
#' @return A list with bin counts and means and results from LOESS predictions
#'
#' @details Makes a plot
#'
#' @export
plotLoessPred <- function(resultMat, nSim=nrow(resultMat), 
                     xlim=NULL, ylim=NULL,
                     budg1=1, budg2=2, binMeanContrast=3, 
                     plotMn=T, plotHiMn=T, plotHiCt=F, plotHiPredGain=F,
                     giveRange=T){
  require(hexbin)
  require(grid)
  lpc <- loessPredCount(resultMat=resultMat, nSim=nSim, xlim=xlim, 
                        ylim=ylim, budg1=budg1, budg2=budg2)
  
  rangeTitle <- " Simulations"
  if (giveRange) rangeTitle <- paste0(" Sims. Range: ", paste(round(range(lpc$binMean), 1), collapse=" to "))
  prefTitle <- paste0(rep(" ", 4 - ceiling(log10(nSim))), collapse="")
  main <- paste0(prefTitle, nSim, rangeTitle)
  if (plotMn){
    bmc <- binMeanContrast
    binRange <- diff(range(lpc$binMean))^bmc
    meanAsCount <- round(99*(lpc$binMean - min(lpc$binMean))^bmc / binRange) + 1
    lpc$bins@count <- meanAsCount
  }
  p <- hexbin::gplot.hexbin(lpc$bins, main=main, legend=ifelse(plotMn, 0, 1))
  pushHexport(p$plot.vp)
  if (plotHiMn){
    grid::grid.points(lpc$hiMeanXY[1], lpc$hiMeanXY[2], gp=gpar(col="red"), 
                      pch=16)
  }
  if (plotHiCt){
    grid::grid.points(lpc$hiCtXY[1], lpc$hiCtXY[2], gp=gpar(col="green"), 
                      pch=16)
  }
  if (plotHiPredGain){
    grid::grid.points(lpc$hiPredXY[1], lpc$hiPredXY[2], gp=gpar(col="blue"), 
                      pch=16)
  }
  upViewport()
  return(lpc)
}#END plotLoessPred
