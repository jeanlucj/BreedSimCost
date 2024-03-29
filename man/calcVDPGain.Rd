% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VarietyDevPipe.R
\name{calcVDPGain}
\alias{calcVDPGain}
\title{calcVDPGain function}
\usage{
calcVDPGain(nEntryVec, nFinal, cov)
}
\arguments{
\item{nEntryVec}{nStages vector with number of individuals at each stage}

\item{nFinal}{numeric how many individuals are going to the marketing dept}

\item{cov}{matrix covariances as came out of calcVDPcovMat}
}
\value{
numeric expected gain from the VDP
}
\description{
calcVDPGain function
}
\details{
Call this function on a grid of budget equal VDPs. Wrapper function 
to call the truncation point and gain functions with elements from the grid
}
\examples{
gains <- sapply(testGrid, calcVDPGain, 
                nFinal=bsd$nToMarketingDept, cov=vdpCovMat)

}
