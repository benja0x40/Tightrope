% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BRD_V3.R
\name{BackgroundCandidates}
\alias{BackgroundCandidates}
\title{Index of background candidates}
\usage{
BackgroundCandidates(brd)
}
\arguments{
\item{brd}{result of a prior call to the \link{BRD} function.}
}
\value{
\code{BackgroundCandidates} returns an integer vector.
}
\description{
Provides the index of background candidates within rows of the read count
matrix used to identify these candidates.
}
\seealso{
\link{BRD},
  \link{PlotBRD},
  \link{PlotCountDistributions}
}
