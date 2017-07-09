# =============================================================================.
#
# -----------------------------------------------------------------------------.
VerifyInputs.RCM <- function(cnt, col, aln, grg) {
  msg <- NULL
  if(! is.character(col)) msg <- c(msg, "col is not of a character string")
  if(! is.matrix(cnt)) msg <- c(msg, "cnt is not a matrix")
  if(! is(grg, "GRanges")) msg <- c(msg, "grg is not a GRanges object")
  if(is.character(aln)) {
    if(! grepl("\\.(bam|rdata)$", aln, perl = T, ignore.case = T)) {
      msg <- c(msg, "aln file path is not ending by .bam or .rdata")
    } else {
      if(! file.exists(aln)) msg <- c(msg, "aln file not found")
    }
  } else {
    chk <- is(aln)
    if(! any(c("GAlignments", "GRanges") %in% chk)) {
      msg <- c(msg, "aln is not a GAlignments or GRanges object")
    }
  }
  if(is.null(msg)) {
    if(is.null(colnames(cnt))) msg <- c(msg, "missing column names in cnt")
    if(length(grg) != nrow(cnt)) {
      msg <- c(msg, "the length of grg doesn't match the number of rows in cnt")
    }
    if(! col %in% colnames(cnt)) {
      msg <- c(msg, "col is not a column name of cnt")
    }
  }
  if(! is.null(msg)) msg <- paste(msg, collapse = "\n", sep = "")
  msg
}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
AppendCounts <- function(cnt, col, aln, grg, ...) {

  lbl <- colnames(cnt)
  if(col %in% lbl) stop(col, " column already exist in cnt matrix")

  i <- ncol(cnt) + 1
  cnt <- cbind(cnt, countOverlaps(query = grg, subject = aln, ...))
  colnames(cnt)[i] <- col

  cnt
}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
UpdateCounts <- function(cnt, col, aln, grg, ...) {

  lbl <- colnames(cnt)
  if(! col %in% lbl) stop(col, " column doesn't exist in cnt matrix")

  i <- match(col, colnames(cnt))
  cnt[, i] <- countOverlaps(query = grg, subject = aln, ...)

  cnt
}

# =============================================================================.
#' count reads overlapping with genomic intervals
# -----------------------------------------------------------------------------.
# TODO: implement skipping existing count columns with provided counts
# TODO: check if this can be faster and as/more flexible using summerizeOverlaps
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{NonZeroCounts},
#' \link{DitherCounts},
#' \link{JoinColumns},
#' \link{ExtractColumns},
#' \link{makeReadCountsInROI}
# -----------------------------------------------------------------------------.
#' @param lst
#' list of \link{GAlignments} objects with mapped reads from different samples
#'
#' @param grg
#' \link{GRanges} object defining the considered genomic intervals
#'
#' @return
#' \code{ReadCountMatrix} returns a matrix of read counts with
#' rows = observations and columns = measurement conditions.
# -----------------------------------------------------------------------------.
#' @export




# Inputs: list of GAlignments, GRanges, rdata or bam files paths
ReadCountMatrix <- function(lst, grg, cnt = NULL, ...) {


  lst <- dir(
    "/Volumes/USB16GB/LT_WORKS/NextSeq_ESC_BRD_R1_TEST/_MAPPEDREADS_/bowtie2/genomic_dm6/",
    full.names = T
  )
  names(lst) <- gsub("\\.[^\\.]+$", "", basename(lst), perl = T)


  if(! is.null(cnt)) {
  } else {
    cnt <-  matrix(0, length(grg), 0)
  }
  # File paths
  chk <- all(sapply(lapply(lst, is), "[[", 1) == "character")
  if(chk) chk <- all(file.exists(lst))
  if(chk) chk <- all(grepl("\\.(bam|rdata)$", lst, perl = T, ignore.case = T))
  if(chk) {
  }

  cnt <- matrix(0, length(grg), length(lst))
  for(i in 1:length(lst)) {
    message("Read counts for sample ", i)
    cnt[,i] <- countOverlaps(
      query=grg, subject=lst[[i]], ...
    )
  }
  cnt
}

# =============================================================================.
#' Merge count data
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{ExtractColumns},
#' \link{RemoveColumns},
#' \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @param x
#' list of read count matrices
#' (rows = observations, columns = measurement conditions).
#'
#' @param y
#' list of read count matrices.
#'
#' @return
#' \code{JoinColumns} returns a list.
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
#' @seealso
#' \link{JoinColumns},
#' \link{RemoveColumns},
#' \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @param x
#' list of read count matrices
#' (rows = observations, columns = measurement conditions).
#'
#' @param lst
#' column names or indices to be extracted.
#'
#' @return
#' \code{ExtractColumns} returns a list.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ExtractColumns <- function(x, lst) {
  for(lbl in names(x)) {
    x[[lbl]] <- x[[lbl]][, idx, drop = F]
  }
  x
}

# =============================================================================.
#' Remove count data
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{JoinColumns},
#' \link{ExtractColumns},
#' \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @param x
#' list of read count matrices
#' (rows = observations, columns = measurement conditions).
#'
#' @param lst
#' names or indices of columns to be removed.
#'
#' @return
#' \code{RemoveColumns} returns a list.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
RemoveColumns <- function(x, lst) {
  idx <- lst
  for(lbl in names(x)) {
    if(is.character(lst)) idx <- match(lst, colnames(x[[lbl]]))
    x[[lbl]] <- x[[lbl]][, - idx, drop = F]
  }
  x
}
