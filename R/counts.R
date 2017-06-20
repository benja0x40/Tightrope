# FUNCTIONS | READ COUNTS ######################################################

# =============================================================================.
#' Detect computable values (i.e. not NA nor Inf)
# -----------------------------------------------------------------------------.
#' @param x
#' numeric vector or matrix
#'
#' @return
#' FiniteValues returns a \code{logical} vector
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
FiniteValues <- function(x) {
  if(is.null(dim(x))) {
    x <- sapply(x, FUN = is.finite)
  } else {
    n <- ncol(x)
    x <- t(apply(x, MARGIN = 1, FUN = is.finite))
    x <- rowSums(x) == n
  }
  x
}

# =============================================================================.
#' Read count matrix summary
# -----------------------------------------------------------------------------.
#' @param cnt
#' matrix of read counts (rows = observations, columns = samples or conditions).
#'
#' @param detailed
#' logical, if TRUE return min and max read counts in each interval
#'
#' @return
#' DetectCounts returns a \code{data.frame} with the following columns:
#' \item{all}{logical vector indicating rows of the \code{cnt} matrix with at least one count in all columns}
#' \item{none}{logical vector indicating rows of the \code{cnt} matrix with no counts in all columns}
#' \item{nbr}{integer vector representing the number of columns with at least one counts for each row of the \code{cnt} matrix}
#' If \code{detailed} is equal to \code{TRUE}, the returned \code{data.frame} also includes the following columns:
#' \item{min}{integer vector representing the maximum count value for each row of the \code{cnt} matrix}
#' \item{max}{integer vector representing the minimum count value for each row of the \code{cnt} matrix}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
DetectCounts <- function(cnt, detailed = F) {
  cnt <- as.matrix(cnt)
  n <- ncol(cnt)
  k <- apply(cnt, MARGIN = 1, FUN = function(x) { sum(x > 0) })
  res <- data.frame(
    all  = k == n,
    none = k == 0,
    nbr  = k
  )
  if(detailed) {
    idx <- which(k > 0)
    res$min <- 0
    res$min[idx] <- apply(cnt[idx, ], 1, FUN=function(x) { min(x[x > 0]) })
    res$max <- apply(cnt, MARGIN = 1, FUN = max)
  }
  res
}

# =============================================================================.
#' Dither counts
# -----------------------------------------------------------------------------.
#' @param x
#' matrix of read counts (rows = observations, columns = samples or conditions).
#'
#' @return DitherCounts returns a \code{matrix} of dithered counts
# -----------------------------------------------------------------------------.
#' @export
DitherCounts <- function(x) {

  zero <- x == 0
  xmin <- min(x[! zero], na.rm = T)
  xmax <- max(x[! zero], na.rm = T)

  # dithering
  x <- x + rtriangle(length(x), a = -1, b = 1)
  # lower limit
  k <- x < xmin & ! zero
  x[k] <- runif(sum(k), xmin - 0.5, xmin)
  # upper limit
  k <- x > xmax
  x[k] <- runif(sum(k), xmax, xmax + 0.5)
  # no count
  x[zero] <- 0

  x
}

# =============================================================================.
#' Merge count data
# -----------------------------------------------------------------------------.
#' @param x
#' list of read count matrices
#' (rows = observations, columns = samples or conditions).
#'
#' @param y
#' list of read count matrices.
#'
#' @return JoinColumns returns a \code{list}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
JoinColumns <- function(x, y) {
  for(lbl in names(x)) {
    x[[lbl]] <- cbind(x[[lbl]], y[[lbl]])
  }
  x
}
# =============================================================================.
#' Extract count data
# -----------------------------------------------------------------------------.
#' @param x
#' list of read count matrices
#' (rows = observations, columns = samples or conditions).
#'
#' @param lst
#' column names or indices to be extracted
#'
#' @return ExtractColumns returns a \code{list}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ExtractColumns <- function(x, lst) {
  for(lbl in names(x)) {
    x[[lbl]] <- x[[lbl]][, lst]
  }
  x
}

# =============================================================================.
#' Count mapped reads overlapping with genomic intervals
# -----------------------------------------------------------------------------.
# TODO: implement skipping existing count columns with provided counts
# TODO: check if this can be faster and as/more flexible using summerizeOverlaps
# -----------------------------------------------------------------------------.
#' @param aln
#' list of \link{GAlignments} objects with mapped reads from different samples
#'
#' @param grg
#' \link{GRanges} object defining the considered genomic intervals
#'
#' @return
#' makeReadCountMatrix matrix of read counts
#' (rows = observations, columns = samples or conditions).
# -----------------------------------------------------------------------------.
#' @export
makeReadCountMatrix <- function(aln, grg, ...) {
  cnt <- matrix(0, length(grg), length(aln))
  for(i in 1:length(aln)) {
    message("Read counts for sample ", i)
    cnt[,i] <- countOverlaps(
      query=grg, subject=aln[[i]], ...
    )
  }
  cnt
}

# =============================================================================.
#' makeReadCountsInROI
# -----------------------------------------------------------------------------.
# TODO: implement read.extension
# TODO: implement different inputs (bam files, GenomicAlignment)
# -----------------------------------------------------------------------------.
#'
#' @param ALN
#' mapped reads
#'
#' @param ROI
#' genomic intervals
#'
#' @param CNT
#' previsouly computed read counts
#'
#' @param labels
#' sample names
#'
#' @param ignore.strand
#' default = T
#'
#' @param read.extension
#' not yet implemented
#'
#' @param ...
#'
#' @return makeReadCountsInROI returns a \code{list} of read count matrixes
# -----------------------------------------------------------------------------.
#' @export
makeReadCountsInROI <- function(
  ALN, ROI, CNT = NULL, labels = NULL, ignore.strand = T, read.extension = 0,
  ...
) {

  ROI.CNT        <- vector("list", length(ROI))
  names(ROI.CNT) <- names(ROI)

  if(is.null(labels)) labels <- names(ALN)

  for(roi in names(ROI)) {

    if(! is.null(CNT[[roi]])) { # Skip if counts already provided in CNT
      ROI.CNT[[roi]] <- CNT[[roi]]
    } else {

      txt_out(x = "=")
      message(roi)
      txt_out(x = "-")
      cnt <- makeReadCountMatrix(
        ALN, ROI[[roi]], ignore.strand = ignore.strand, ...
      )
      colnames(cnt) <- labels
      ROI.CNT[[roi]] <- cnt
    }
  }

  ROI.CNT
}
