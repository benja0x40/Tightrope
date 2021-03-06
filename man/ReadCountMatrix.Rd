% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadCountMatrix.R
\name{ReadCountMatrix}
\alias{ReadCountMatrix}
\title{Compute read counts over genomic intervals for multiple bam files}
\usage{
ReadCountMatrix(bam.files, grg, paired, names = NULL, maxgap = -1L,
  minoverlap = 0L, type = "any", ignore.strand = TRUE, ...)
}
\arguments{
\item{bam.files}{character vector containing paths to the bam files (mapped reads).}

\item{grg}{a \link{GRanges} object defining the considered genomic intervals.}

\item{paired}{logical vector, must be set to TRUE for files with paired-end reads
and to FALSE for files with single-end reads.
For conciseness, a single logical value can be provided when all bam files
are consistent with this value.}

\item{names}{user provided column names for the resulting count matrix.
When omitted, the bam file names are used.}

\item{maxgap}{
    See \code{?\link[IRanges]{findOverlaps}} in the \pkg{IRanges} package
    for a description of these arguments.
  }

\item{minoverlap}{
    See \code{?\link[IRanges]{findOverlaps}} in the \pkg{IRanges} package
    for a description of these arguments.
  }

\item{type}{
    See \code{?\link[IRanges]{findOverlaps}} in the \pkg{IRanges} package
    for a description of these arguments.
  }

\item{ignore.strand}{
    When set to \code{TRUE}, the strand information is ignored in the
    overlap calculations.
  }

\item{...}{optional arguments forwarded to
\link{readGAlignments} or \link{readGAlignmentPairs}.}
}
\value{
\code{ReadCountMatrix} returns a read counts matrix where
rows correspond to the provided genomic intervals
and columns to the provided bam files.
}
\description{
Compute read counts over genomic intervals for multiple bam files
}
\seealso{
\link{DitherCounts},
  \link{NormalizeCountMatrix},
  \link{PlotCountDistributions}
}
