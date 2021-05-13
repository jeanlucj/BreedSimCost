# GOING FOR LOW on both criteria

#' idx1 function
#'
#' function to return the number of points lower than the query point
#' in the first column of a matrix
#'
#' @param pt Real query point
#' @param nonDom Real matrix 
#' @param var1name Character: name of the column of column 1
#' @return Integer number of values in nonDom <= to pt
#' 
#' @details Utility for Pareto front
idx1 <- function(pt, nonDom, var1name){
  sum(nonDom[, var1name] <= c(pt[,var1name]))
}

#' insertRow function
#'
#' Insert a vector (vec) into a matrix (mat) at row idx
#'
#' @param vec Real row vector
#' @param idx Integer insertion point 
#' @param mat Real matrix 
#' @return Augmented matrix
#' 
#' @details Utility for Pareto front
insertRow <- function(vec, idx, mat){
  top <- NULL; if (idx) top <- mat[1:idx,,drop=F]
  return(rbind(top, vec, mat[-(1:idx),,drop=F]))
}

#' testDom function
#'
#' Test if any of the points in nonDom are dominated by pt
#'
#' @param pt Real double
#' @param nonDom Real two-column matrix 
#' @param var1name Character: name of the column of column 1
#' @param var2name Character: name of the column of column 2
#' @return Logical vector size number of matrix rows
#' 
#' @details Utility for Pareto front
testDom <- function(pt, nonDom, var1name, var2name){
  idx <- sum(nonDom[, var1name] < c(pt[,var1name]))
  c(rep(F, idx), nonDom[-(1:idx), var2name] > c(pt[,var2name]))
}

#' checkAddRemove function
#'
#' Test if any of the points in nonDom are dominated by pt
#'
#' @param cand Real double that you want to check if it is dominated or not
#' @param nonDom Real two-column matrix, with vectors in rows that are 
#' non-dominated and sorted from low to high for their first column
#' @param var1name Character: name of the column of column 1
#' @param var2name Character: name of the column of column 2
#' @return A new nonDom matrix, adding cand if it is non-dominated and
#' removing any rows that are dominated by it
#' 
#' @details Utility for Pareto front
checkAddRemove <- function(cand, nonDom, var1name, var2name){
  idx <- idx1(cand, nonDom, var1name)
  # if idx == 0 cand is lower than all nonDom in the first column
  # so it cannot be dominated
  # else, there are some points lower than it, and if any are also lower
  # in the second coordinate, the point is dominated
  dominated <- idx > 0 & any(nonDom[1:idx,var2name] < c(cand[,var2name]))
  if (!dominated){
    nonDom <- insertRow(cand, idx, nonDom)
    newDom <- testDom(cand, nonDom, var1name, var2name)
    if (any(newDom)){
      nonDom <- nonDom[!newDom,]
    }
  }
  return(nonDom)
}

#' findNonDom function
#'
#' Find the non-dominated points in a matrix
#'
#' @param mat Real matrix with named columns
#' @param dir1Low Logical: desired direction in down for column 1
#' @param dir2Low Logical: desired direction in down for column 2
#' @param var1name Character: name of the column of column 1
#' @param var2name Character: name of the column of column 2
#' @return nonDom Real matrix, with rows that are non-dominated with 
#' respect to var1 and var2 and sorted for var1
#' 
#' @details Utility for Pareto front
#' @export
findNonDom <- function(mat, dir1Low=T, dir2Low=T, var1name=NULL, var2name=NULL){
  if (is.null(var1name)) var1name <- 1
  if (is.null(var2name)) var2name <- 2
  # function tries to go low for both columns so reverse columns if F
  if (!dir1Low) mat[, var1name] <- -mat[, var1name]
  if (!dir2Low) mat[, var2name] <- -mat[, var2name]
  # find the minima of columns 1 and 2
  min1 <- which(mat[, var1name] == min(mat[, var1name]))
  if (length(min1) > 1) min1 <- min1[which.min(unlist(mat[min1, var2name]))]
  nonDom <- mat[min1,,drop=F]
  mat <- mat[-min1,,drop=F]
  min2 <- which(mat[, var2name] == min(mat[, var2name]))
  if (length(min2) > 1) min2 <- min2[which.min(unlist(mat[min2, var1name]))]
  nonDom <- rbind(nonDom, mat[min2,,drop=F])
  mat <- mat[-min2,,drop=F]
  for (i in 1:nrow(mat)){
    nonDom <- checkAddRemove(mat[i,,drop=F], nonDom, var1name, var2name)
  }
  if (!dir1Low) nonDom[, var1name] <- -nonDom[, var1name]
  if (!dir2Low) nonDom[, var2name] <- -nonDom[, var2name]
  return(nonDom)
}
