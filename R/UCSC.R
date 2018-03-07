# =============================================================================.
#' Retrieve chromosome information from UCSC and build Seqinfo object
# -----------------------------------------------------------------------------.
#' @param x
#' character.
#'
#' @return
#' \code{BuildSeqInfo} returns a \link{Seqinfo} object.
# -----------------------------------------------------------------------------.
#' @export
BuildSeqInfo <- function(x) {
  s <- fetchExtendedChromInfoFromUCSC(genome = x)
  s <- with(
    s, Seqinfo(
      seqnames = UCSC_seqlevel, seqlengths = UCSC_seqlength,
      isCircular = circular, genome = x
    )
  )
  s
}

# HIDDEN #######################################################################

# =============================================================================.
#' Import UCSC cpgIslandExt.txt file format
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{ImportGenomicRanges}
# -----------------------------------------------------------------------------.
#' @param fpath
#' file path or url.
#'
#' @param ...
#' optional parameters forwarded to the \link{ImportGenomicRanges} function.
#'
#' @return
#' \code{ImportCpGIslandExt} returns a \link{GRanges} object.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ImportCpGIslandExt <- function (fpath, ...) {
  grg <- ImportGenomicRanges(
    fpath, chr=2, start=3, end=4,
    xidx = c(1,5:11),
    xlbl = c(
      "bin",
      "name",
      "length",
      "cpgNum",
      "gcNum",
      "perCpg",
      "perGc",
      "obsExp"
    ),
    ...
  )
  grg
}

# =============================================================================.
#' Download and import UCSC's liftOver chain files
# -----------------------------------------------------------------------------.
#' @param url
#' character.
#'
#' @param path
#' character.
#'
#' @return
#' \code{ImportLiftOverChain} returns a \link{`Chain-class`} object.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ImportLiftOverChain <- function(url, path = "") {
  lfp <- MakePath(path, basename(url))
  download.file(url, lfp)
  system(paste("gunzip", lfp))
  lfp <- gsub("\\.gz$", "", lfp)
  chain <- import.chain(lfp)
  file.remove(lfp)
  chain
}
