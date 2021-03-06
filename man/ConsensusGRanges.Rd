% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenomicRanges.R
\name{ConsensusGRanges}
\alias{ConsensusGRanges}
\title{Consensus genomic ranges}
\usage{
ConsensusGRanges(grl, minoverlap, ignore.strand = TRUE, seqinfo = NULL)
}
\arguments{
\item{grl}{a list of \link{GRanges} objects refering to the same genome build.}

\item{minoverlap}{
    See \code{?\link[IRanges]{findOverlaps}} in the \pkg{IRanges} package
    for a description of these arguments.
  }

\item{ignore.strand}{
    When set to \code{TRUE}, the strand information is ignored in the
    overlap calculations.
  }

\item{seqinfo}{a \link{Seqinfo} object.}
}
\value{
\code{ConsensusGRanges} returns a \link{GRanges} object.
}
\description{
Consensus genomic ranges
}
\keyword{internal}
