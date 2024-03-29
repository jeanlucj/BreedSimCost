% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CostOptimization.R
\name{optimizeByLOESS}
\alias{optimizeByLOESS}
\title{optimizeByLOESS function}
\usage{
optimizeByLOESS(bsd)
}
\arguments{
\item{bsd}{List of breeding scheme data. Has parameters for optimization too.}
}
\value{
Numeric matrix with all simulation budget allocations, 
gen mean change, gen std dev change, total cost.
}
\description{
Function to optimize a two-part strategy breeding scheme:
1. Simulate a batch using given percentage ranges
2. Perform LOESS fit to the gains
3. Find budget with best estimated gain
4. Calculate new percentage ranges: any simulation within 2*StdErr of best
5. Decide on some simulations to repeat:
   1. Parameter space with high gain and high std err: need more info there
   2. Parameter space with high gain: high probability that it's best
Go back to 1.
}
\details{
A wrapper to repeatedly simulate a scheme with different 
budget allocations to find optimal allocations
}
