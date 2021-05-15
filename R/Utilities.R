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
#' phenoRecords an empty tibble for phenotype records, and inventory an empty
#' tibble for variety candidates
#'
#' @details Call this function at the beginning of the simulation
#'
#' @examples
#' bsd <- initializeProgram("FounderCtrlFile.txt")
#'
#' @export
initializeProgram <- function(founderFile, schemeFile, 
                              costFile, optimizationFile){
  # Read parameters to create founders
  parmNames <- c("nChr", "effPopSize", "quickHaplo", 
                 "segSites", "nQTL", "nSNP", "genVar", "gxeVar", 
                 "gxyVar", "gxlVar", "gxyxlVar", "meanDD", "varDD", 
                 "relAA")
  bsd <- readControlFile(founderFile, parmNames)

  # Read parameters about the overall scheme
  parmNames <- c("nCyclesToRun", "nBurnInCycles", "nStages", "stageNames", 
                 "nEntries", "nReps", "nLocs", "errVars",
                 "seedNeeded", "seedProduced", "optiContEffPop",
                 "nBreedingProg", "nPopImpCycPerYear", "keepNTrainingCyc", 
                 "keepNBreedingCyc", "varietiesCanBeParents")
  bsdNew <- readControlFile(schemeFile, parmNames)
  bsd <- c(bsd, bsdNew)
  
  # Read parameters about scheme costs
  parmNames <- c("plotCosts", "perLocationCost", "crossingCost", 
                 "candidateDevelCost", "qcGenoCost", "wholeGenomeCost")
  bsdNew <- readControlFile(costFile, parmNames)
  bsd <- c(bsd, bsdNew)
  
  # Read parameters about optimization procedure
  parmNames <- c("nCores", "minPercentage", "maxPercentage",
                 "percentageStep", "minNBreedingProg",
                 "tolerance", "batchSize", "maxNumBatches",
                 "nHighGain", "nUncertain", "debug", 
                 "verbose", "saveIntermediateResults")
  bsdNew <- readControlFile(optimizationFile, parmNames)
  bsd <- c(bsd, bsdNew)
  
  # Add miscelaneous data structures
  bsd$year <- 0
  bsd$nextTrialID <- 1
  bsd <- calcDerivedParms(bsd)
  bsd <- calcBudget(bsd)
  
  # Create haplotypes for founder population of outbred individuals
  if (bsd$quickHaplo){
    founderHap <- quickHaplo(nInd=bsd$nBreedingProg, 
                             nChr=bsd$nChr, segSites=bsd$segSites)
  } else{
    founderHap <- runMacs2(nInd=bsd$nBreedingProg, 
                           nChr=bsd$nChr, segSites=bsd$segSites, 
                           Ne=bsd$effPopSize)
  }
  
  # Global simulation parameters from founder haplotypes
  bsd$SP <- SimParam$new(founderHap)
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
  bsd$saveIntermediateResults <- makeLogical(bsd$saveIntermediateResults)
  bsd$varietiesCanBeParents <- makeLogical(bsd$varietiesCanBeParents)
  
  # Genetic architecture defaults
  if (nv(bsd$meanDD)) bsd$meanDD <- 0
  if (nv(bsd$varDD)) bsd$varDD <- 0
  if (nv(bsd$relAA)) bsd$relAA <- 0

  # Check that these vectors are of the right length
  rightLength <- function(objName) length(get(objName, bsd)) == bsd$nStages
  vecNames <- c("stageNames", "nEntries", "nReps", "nLocs", "errVars",
                "seedNeeded", "seedProduced")
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
  names(bsd$nEntries) <- names(bsd$nReps) <- names(bsd$nLocs) <- names(bsd$errVars) <-
    names(bsd$seedNeeded) <- names(bsd$seedProduced)<- bsd$stageNames

  return(bsd)
}

#' calcCurrentStatus function
#'
#' Calculate mean and variance statistics on your current breeding program.
#'
#' @param bsd List of breeding program data
#' @return A vector extracting breeding population means and variances
#'
#' @details This function is only called internally by other functions used to specify the pipeline
#'
#' @export
calcCurrentStatus <- function(bsd){
  nInd <- nInd(bsd$breedingPop)
  currBreedPop <- bsd$breedingPop[nInd - (bsd$nBreedingProg - 1):0]
  breedPopMean <- mean(gv(currBreedPop))
  breedPopSD <- sd(gv(currBreedPop))
  nInd <- nInd(bsd$varietyCandidates)
  bestVarCand <- bsd$varietyCandidates[nInd - (bsd$nEntries[bsd$nStages] - 1):0]
  varCandMean <- mean(gv(bestVarCand))
  return(c(breedPopMean=breedPopMean, breedPopSD=breedPopSD, 
           varCandMean=varCandMean))
}

#' loessForNsim function NOT FINISHED
#'
#' Make a nice hexbin plot with results from optimization.
#'
#' @param bsd List of breeding program data
#' @return A vector extracting breeding population means and variances
#'
#' @details This function is only called internally by other functions used to specify the pipeline
#'
#' @export
loessForNsim <- function(nSim, allBatches, nToRepeat, rand=FALSE, 
                         xlim=NULL, ylim=NULL, budg1=1, budg2=2){
  # Frame for hexbin plots
  if (is.null(xlim)){
    xlim <- c(floor(min(allBatches[,budg1])*50-0.5)/50, 
              ceiling(max(allBatches[,budg1])*50+0.5)/50)
  }
  if (is.null(ylim)){
    ylim <- c(floor(min(allBatches[,budg2])*50-0.5)/50, 
              ceiling(max(allBatches[,budg2])*50+0.5)/50)
  }

  sims <- 1:nSim
  if (rand) sims <- sort(sample(nrow(allBatches), nSim))
  uptoBatches <- allBatches[sims,]
  
  stageNames <- colnames(allBatches)[grep("budget", colnames(allBatches))]
  stageNames <- stageNames %>% substring(8)
  
  # Non-Parametric LOESS response
  loFormula <- paste0("response ~ ", paste0(paste0("budget.", stageNames[-1]), collapse=" + "))
  loFM <- loess(loFormula, data=uptoBatches)
  loPred <- predict(loFM, se=T)
  
  # choose which have high response and high std. err. of response
  # Use Pareto front, looking for high fit and high std. err.
  paretoRepeat <- NULL
  fitStdErr <- tibble(batchID=1:nrow(uptoBatches), fit=loPred$fit, se=loPred$se.fit)
  while (length(paretoRepeat) < nToRepeat){
    nonDomSim <- findNonDom(fitStdErr, dir1Low=F, dir2Low=F, var1name="fit", var2name="se")
    nds <- nonDomSim$batchID
    if (length(nds) > nToRepeat - length(paretoRepeat)){
      nds <- sample(nds, nToRepeat - length(paretoRepeat))
    }
    paretoRepeat <- c(paretoRepeat, nds)
    fitStdErr <- fitStdErr %>% dplyr::filter(!(batchID %in% nds))
  }
  fitStdErr <- tibble(batchID=1:nrow(uptoBatches), fit=loPred$fit, se=loPred$se.fit)
  plot(fitStdErr$fit, fitStdErr$se)
  points(fitStdErr$fit[paretoRepeat], fitStdErr$se[paretoRepeat], pch=16, col=2, cex=0.8)
  # Use ones that have the best chance of beating the current best
  bestGain <- max(fitStdErr$fit)
  probBetter <- pnorm(bestGain, fitStdErr$fit, fitStdErr$se, lower.tail=FALSE)
  hiProb <- order(probBetter, decreasing=T)[1:nToRepeat]
  hist(probBetter[hiProb])
  hist(hiProb)
  plot(fitStdErr$fit, fitStdErr$se)
  points(fitStdErr$fit[hiProb], fitStdErr$se[hiProb], pch=16, col=2, cex=0.8)
  
  bestFit <- which.max(loPred$fit)
  bestSE <- loPred$se.fit[bestFit]
  
  # Plot all and circle best
  bins <- hexbin(cbind(uptoBatches[,budg1], uptoBatches[,budg2]), xbnds=xlim, ybnds=ylim)
  p <- plot(bins, xlab=names(budg1), ylab=names(budg2), main=nSim)
  pushHexport(p$plot.vp)
  grid.points(uptoBatches[bestFit, budg1], uptoBatches[bestFit, budg2], gp=gpar(col="red"))
  upViewport()
  
  # Plot Pareto repeats
  bins <- hexbin(cbind(uptoBatches[paretoRepeat, budg1], uptoBatches[paretoRepeat, budg2]), xbnds=xlim, ybnds=ylim)
  pdf(paste0("_pdfs/RepPareto", nSim, ".pdf"))
  p <- plot(bins, xlab=names(budg1), ylab=names(budg2), main=paste("Rep. Pareto", nSim), legend=FALSE)
  pushHexport(p$plot.vp)
  grid.points(uptoBatches[bestFit, budg1], uptoBatches[bestFit, budg2], gp=gpar(col="red"))
  upViewport()
  dev.off()
  
  # Plot high probability repeats
  bins <- hexbin(cbind(uptoBatches[hiProb, budg1], uptoBatches[hiProb, budg2]), xbnds=xlim, ybnds=ylim)
  pdf(paste0("_pdfs/RepBest", nSim, ".pdf"))
  p <- plot(bins, xlab=names(budg1), ylab=names(budg2), main=paste("Rep. Best", nSim), legend=FALSE)
  pushHexport(p$plot.vp)
  grid.points(uptoBatches[bestFit, budg1], uptoBatches[bestFit, budg2], gp=gpar(col="red"))
  upViewport()
  dev.off()
  
  whichClose <- which(max(loPred$fit) - loPred$fit < bestSE)
  print(paste("Number of Sims <1 StdErr to best", length(whichClose)))
  nClose <- sum(max(loPred$fit) - loPred$fit < 2*bestSE)
  return(c(bestFit=loPred$fit[bestFit], bestSE=bestSE, nClose=nClose, loPred=loPred, probBetter=probBetter))
}
