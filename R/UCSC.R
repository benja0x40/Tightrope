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
