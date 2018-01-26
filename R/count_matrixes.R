# =============================================================================.
#' StopOnError
# -----------------------------------------------------------------------------.
#' @param msg
#' error message.
#'
#' @param chk
#' logical vector.
#'
#' @param x
#' character vector.
#'
#' @return
#' \code{StopOnError} stops execution in case of error or returns FALSE.
# -----------------------------------------------------------------------------.
#' @keywords internal
StopOnError <- function(msg, chk, x) {
  if(! all(chk)) {
    chk <- paste(paste0("[", which(! chk), "] ", x[! chk]), collapse = "\n")
    msg <- paste(msg, chk, sep = "\n")
    message("Context: ", format(sys.call(-1)))
    stop(msg, call. = F, domain = NA)
  }
  F
}

# =============================================================================.
#' VerifyInputs
# -----------------------------------------------------------------------------.
#' @param cnt
#' a numeric matrix is expected.
#'
#' @param grg
#' a GRanges object is expected.
#'
#' @return
#' \code{VerifyInputs} stops execution in case of invalid input or returns TRUE.
# -----------------------------------------------------------------------------.
#' @keywords internal
VerifyInputs <- function(cnt, grg) {
  chk <- rep(T,  2)
  msg <- rep("", 2)
  if(! (is.matrix(cnt) | is.null(cnt))) {
    chk[1] <- F
    msg[1] <- "matrix required"
  }
  if(! is(grg, "GRanges")) {
    chk[2] <- F
    msg[2] <- "GRanges object required"
  }
  if(all(chk) & ! is.null(cnt)) {
    if(is.null(colnames(cnt))) {
      chk[1] <- F
      msg[1] <- "missing column names in count matrix"
    }
    if(length(grg) != nrow(cnt)) {
      chk[2] <- F
      msg[2] <- "the number of observations in the count matrix doesn't match the number of genomic ranges"
    }
  }
  if(! all(chk)) {
    chk <- paste(paste0("[", which(! chk), "] ", msg[! chk]), collapse = "\n")
    msg <- paste("invalid arguments", chk, sep = "\n")
    message("Context: ", format(sys.call(-1)))
    stop(msg, call. = F, domain = NA)
  }
  T
}

# =============================================================================.
# TODO: paired reads
# =============================================================================.
#' load/verify a list of mapped read files
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{MakeReadCounts},
#' \link{ReadCountMatrix},
#' \link{AppendReadCounts},
#' \link{UpdateReadCounts}
# -----------------------------------------------------------------------------.
#' @param aln
#' a vector of bam file paths or a list with named elements and containing
#' exclusively \link{GAlignments} objects.
#'
#' @param NameByFile
#' use file names to generate the name of each measurement condition, meaning
#' each sequenced biological sample (logical, default = F).
#'
#' @param ...
#' additional arguments passed to the \link{readGAlignments} function.
#'
#' @return
#' \code{LoadMappedReads} returns a list of \link{GAlignments} objects.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
LoadMappedReads <- function(aln, NameByFile = F, ...) {

  rex <- "^(.*)\\.(bam|rdata)$"

  if(is.null(names(aln)) & ! (is.character(aln) & NameByFile)) {
    stop("names are required to designate each measurement condition")
  }
  if(is.character(aln)) {

    chk <- grepl(rex, aln, perl = T, ignore.case = T)
    if(! StopOnError("file path not ending by .bam or .rdata", chk, aln)) {

      chk <- file.exists(aln)
      if(! StopOnError("mapped reads file not found", chk, aln)) {

        if(NameByFile) {
          names(aln) <- gsub(
            rex, "\\1", basename(aln), perl = T, ignore.case = T
          )
        }
        ext <- gsub(rex, "\\2", aln, perl = T, ignore.case = T)
        message("loading mapped reads:\n", paste(names(aln), collapse = "\n"))

        if(all(ext == "bam"))   aln <- lapply(aln, readGAlignments, ...)
        if(all(ext == "rdata")) aln <- lapply(aln, readRDS)

      }
    }
  } else {
    if(! is.list(aln)) {
      msg <- paste(
        "mapped reads must be provided as a vector of bam file paths",
        "or as a list containing GAlignments objects"
      )
      stop(msg)
    }
    obj <- sapply(lapply(aln, is), "[[", 1)
    chk <- obj %in% c("GAlignments", "GRanges")
    StopOnError("GAlignments or GRanges object(s) required", chk, obj)
  }

  aln
}

# =============================================================================.
#' append read counts for additional conditions
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{LoadMappedReads},
#' \link{UpdateReadCounts},
#' \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @inheritParams GenomicRanges::countOverlaps
#'
#' @param cnt
#' a numeric matrix representing read counts,
#' where rows = observations and columns = measurement conditions.
#'
#' @param aln
#' a vector of bam file paths or a list with named elements and containing
#' exclusively \link{GAlignments} objects.
#'
#' @param grg
#' a \link{GRanges} object defining the genomic intervals of interest.
#'
#' @param NameByFile
#' use file names to generate the name of each measurement condition, meaning
#' each sequenced biological sample (logical, default = F).
#'
#' @param ...
#' additional arguments passed to the \link{readGAlignments} function.
#'
#' @return
#' functions \code{ReadCountMatrix}, \code{AppendReadCounts}, and
#' \code{UpdateReadCounts} return a read counts matrix where
#' rows = observations and columns = measurement conditions.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
AppendReadCounts <- function(
  cnt, aln, grg, NameByFile = F,
  maxgap = 0L, minoverlap = 1L, type = "any", ignore.strand = T, ...
) {

  VerifyInputs(cnt, grg)
  cnt.col <- colnames(cnt)

  aln <- LoadMappedReads(aln, NameByFile = NameByFile, ...)
  aln.col <- names(aln)

  chk <- ! aln.col %in% cnt.col
  StopOnError("column already exist in read count matrix", chk, aln.col)

  cnt <- cbind(
    cnt, matrix(0, length(grg), length(aln), dimnames = list(NULL, aln.col))
  )
  for(i in aln.col) {
    # message("creating read counts for ", aln.col)
    cnt[, i] <- countOverlaps(
      query = grg, subject = aln[[i]],
      maxgap = maxgap, minoverlap = minoverlap, type = type,
      ignore.strand = ignore.strand
    )
  }

  cnt
}

# =============================================================================.
#' update read counts for specific conditions
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{LoadMappedReads},
#' \link{AppendReadCounts},
#' \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @inherit AppendReadCounts
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
UpdateReadCounts <- function(
  cnt, aln, grg, NameByFile = F,
  maxgap = 0L, minoverlap = 1L, type = "any", ignore.strand = T, ...
) {

  VerifyInputs(cnt, grg)
  cnt.col <- colnames(cnt)

  aln <- LoadMappedReads(aln, NameByFile = NameByFile, ...)
  aln.col <- names(aln)

  chk <- aln.col %in% cnt.col
  StopOnError("column doesn't exist in read count matrix", chk, aln.col)

  for(i in aln.col) {
    # message("updating read counts for ", aln.col)
    cnt[, i] <- countOverlaps(
      query = grg, subject = aln[[i]],
      maxgap = maxgap, minoverlap = minoverlap, type = type,
      ignore.strand = ignore.strand
    )
  }

  cnt
}

# =============================================================================.
#' count reads in genomic intervals for different conditions
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{LoadMappedReads},
#' \link{AppendReadCounts},
#' \link{UpdateReadCounts},
#' \link{NonZeroCounts},
#' \link{DitherCounts}
# -----------------------------------------------------------------------------.
#' @description
#' \code{ReadCountMatrix} computes the number of reads overlapping with provided
#' genomic intervals for different measurement conditions.
#'
#' @inherit AppendReadCounts
#'
#' @param update
#' instead of producing new count matrixes from scratch, this argument allows
#' to update targeted count columns when existing count matrixes are provided
#' (logical, default = F).
# -----------------------------------------------------------------------------.
#' @export
ReadCountMatrix <- function(
  aln, grg, cnt = NULL, update = F, NameByFile = F,
  maxgap = 0L, minoverlap = 1L, type = "any", ignore.strand = T, ...
) {

  if(update & is.null(cnt)) stop("missing the read count matrix to be updated")

  aln <- LoadMappedReads(aln, NameByFile = NameByFile, ...)

  if(is.null(cnt) | ! update) {
    cnt <- matrix(
      0, length(grg), length(aln), dimnames = list(NULL, names(aln))
    )
  }

  cnt <- UpdateReadCounts(
    cnt, aln, grg,
    maxgap = maxgap, minoverlap = minoverlap, type = type,
    ignore.strand = ignore.strand
  )

  cnt
}

# =============================================================================.
# TEST: with a subset of aln and grg, only those should be updated...
# TODO: implement read.extension
# =============================================================================.
#' make read counts for different locations and conditions
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{LoadMappedReads},
#' \link{JoinColumns},
#' \link{ExtractColumns},
#' \link{RemoveColumns}
# -----------------------------------------------------------------------------.
#' @description
#' \code{MakeReadCounts} computes read counts for different sets of genomic
#' locations and different measurement conditions.
#'
#' @inheritParams ReadCountMatrix
#'
#' @param locations
#' a list with named elements, each one being a \link{GRanges} object defining
#' a set of genomic intervals.
#'
#' @param ...
#' additional arguments passed to the \link{readGAlignments} functions.
#'
#' @return
#' \code{MakeReadCounts} returns a list fo read count matrices.
# -----------------------------------------------------------------------------.
#' @export
MakeReadCounts <- function(
  aln, locations, cnt = NULL, update = F, NameByFile = F,
  maxgap = 0L, minoverlap = 1L, type = "any", ignore.strand = T, ...
) {

  if(is.null(names(locations))) {
    stop("names are required to designate each set of genomic locations")
  }
  if(update & is.null(cnt)) stop("missing read count matrix to be updated")

  if(is.null(cnt)) cnt <- list()

  for(spl in names(aln)) {

    x <- LoadMappedReads(aln[spl], NameByFile = NameByFile, ...)

    for(roi in names(locations)) {
      message("counting reads ", spl, " | ", roi)

      if(update) f <- UpdateReadCounts else  f <- AppendReadCounts

      cnt[[roi]] <- f(
        cnt[[roi]], aln = x, grg = locations[[roi]],
        maxgap = maxgap, minoverlap = minoverlap, type = type,
        ignore.strand = ignore.strand
      )
    }
  }

  cnt
}

# =============================================================================.
#' Merge count data
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{ExtractColumns},
#' \link{RemoveColumns},
#' \link{MakeReadCounts}
# -----------------------------------------------------------------------------.
#' @param x
#' list of read count matrices
#' (rows = observations, columns = measurement conditions).
#'
#' @param y
#' list of read count matrices.
#'
#' @return
#' \code{JoinColumns} returns a list of read count matrices.
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
#' \link{MakeReadCounts}
# -----------------------------------------------------------------------------.
#' @param x
#' list of read count matrices
#' (rows = observations, columns = measurement conditions).
#'
#' @param lst
#' names or indices of columns to be extracted.
#'
#' @return
#' \code{ExtractColumns} returns a list of read count matrices.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ExtractColumns <- function(x, lst) {
  for(lbl in names(x)) {
    x[[lbl]] <- x[[lbl]][, lst, drop = F]
  }
  x
}

# =============================================================================.
#' Remove count data
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{JoinColumns},
#' \link{ExtractColumns},
#' \link{MakeReadCounts}
# -----------------------------------------------------------------------------.
#' @param x
#' list of read count matrices
#' (rows = observations, columns = measurement conditions).
#'
#' @param lst
#' names or indices of columns to be removed.
#'
#' @return
#' \code{RemoveColumns} returns a list of read count matrices.
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
