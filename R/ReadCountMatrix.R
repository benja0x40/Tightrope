# =============================================================================.
#' Compute read counts over genomic intervals for multiple bam files
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{DitherCounts},
#'   \link{NormalizeCountMatrix},
#'   \link{PlotCountDistributions}
# -----------------------------------------------------------------------------.
#' @inheritParams GenomicRanges::countOverlaps
#'
#' @param bam.files
#' character vector containing paths to the bam files (mapped reads).
#'
#' @param grg
#' a \link{GRanges} object defining the considered genomic intervals.
#'
#' @param paired
#' logical vector, must be set to TRUE for files with paired-end reads
#' and to FALSE for files with single-end reads.
#' For conciseness, a single logical value can be provided when all bam files
#' are consistent with this value.
#'
#' @param names
#' user provided column names for the resulting count matrix.
#' When omitted, the bam file names are used.
#'
#' @param ...
#' optional arguments forwarded to
#' \link{readGAlignments} or \link{readGAlignmentPairs}.
#'
#' @return
#' \code{ReadCountMatrix} returns a read counts matrix where
#' rows correspond to the provided genomic intervals
#' and columns to the provided bam files.
# -----------------------------------------------------------------------------.
#' @export
ReadCountMatrix <- function(
  bam.files, grg, paired, names = NULL,
  maxgap = -1L, minoverlap = 0L, type = "any", ignore.strand = TRUE, ...
) {

  rex <- "\\.bam$"
  if(! all(grepl(rex, bam.files, ignore.case = TRUE, perl = TRUE))) {
    stop('some file names lack the ".bam" extension')
  }
  if(! all(sapply(bam.files, file.exists))) {
    stop("some bam files are missing")
  }
  if(is.null(names)) {
    names <- gsub(rex, "", basename(bam.files), ignore.case = TRUE, perl = TRUE)
  }

  nbr <- length(bam.files)
  paired <- rep(paired, nbr, length.out = nbr)
  cnt <- matrix(0, length(grg), nbr)
  for(i in 1:nbr) {
    message("ReadCountMatrix | ", names[i], " => ", bam.files[i])
    if(paired[i]) {
      aln <- readGAlignmentPairs(bam.files[i], ...)
    } else {
      aln <- readGAlignments(bam.files[i], ...)
    }
    cnt[, i] <- countOverlaps(
      grg, as(aln, "GRanges"), maxgap = maxgap, minoverlap = minoverlap,
      type = type, ignore.strand = ignore.strand
    )
  }
  colnames(cnt) <- names

  cnt
}
