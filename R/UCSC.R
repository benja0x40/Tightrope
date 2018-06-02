# =============================================================================.
#' Retrieve chromosome information from UCSC and build Seqinfo object
# -----------------------------------------------------------------------------.
#' @param genome
#' name of an UCSC genome build.
#'
#' @return
#' \code{BuildSeqInfo} returns a \link{Seqinfo} object.
# -----------------------------------------------------------------------------.
#' @export
BuildSeqInfo <- function(genome) {
  s <- fetchExtendedChromInfoFromUCSC(genome = genome)
  s <- with(
    s, Seqinfo(
      seqnames = UCSC_seqlevel, seqlengths = UCSC_seqlength,
      isCircular = circular, genome = genome
    )
  )
  s
}

# =============================================================================.
#' Import CpG-islands from UCSC genomes
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{ImportGenomicRanges},
#'   \link{CleanupGRanges},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @param genome
#' name of an UCSC genome build. When \code{genome} is unspecified, the
#' \code{fpath} argument must be specified.
#'
#' @param fpath
#' path to a cpgIslandExt.txt.gz file. When \code{fpath} is unspecified, the
#' \code{genome} argument must be specified.
#'
#' @param goldenPath
#' base URL of UCSC genome repositories
#' (default = "http://hgdownload.soe.ucsc.edu/goldenPath").
#'
#' @param cpgIslandExt
#' path to cpgIslandExt.txt.gz files in the UCSC genome repositories
#' (default = "database/cpgIslandExt.txt.gz").
#'
#' @param keep.file
#' logical value controlling whether the downloaded file should be kept or
#' deleted (default = FALSE, delete).
#'
#' @param ...
#' optional parameters forwarded to the \link{ImportGenomicRanges} function.
#'
#' @return
#' \code{ImportCpGIslandExt} returns a \link{GRanges} object.
# -----------------------------------------------------------------------------.
#' @export
ImportCpGIslandExt <- function (
  genome = NULL, fpath = NULL,
  goldenPath = "http://hgdownload.soe.ucsc.edu/goldenPath",
  cpgIslandExt = "database/cpgIslandExt.txt.gz",
  keep.file = FALSE, quiet = FALSE
) {

  if(is.null(fpath) & is.null(genome)) stop("missing genome or file path")
  if(is.null(genome) & ! is.null(fpath)) {
    if(! file.exists(fpath)) stop("file not found")
  }

  sinfo <- NULL
  if(is.null(fpath) & ! is.null(genome)) {
    fpath <- paste0("UCSC_cpgIslandExt_", genome, ".txt.gz")
    utils::download.file(
      paste0(goldenPath, "/", genome, "/", cpgIslandExt), destfile = fpath,
      quiet = quiet
    )
    sinfo <- BuildSeqInfo(genome)
  }

  args <- list(
    fpath = fpath, chr = 2, start = 3, end = 4,
    xidx = c(1, 5:11),
    xlbl = c(
      "bin", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp"
    )
  )
  if(! ("seqinfo" %in% names(args) | is.null(sinfo))) {
    args <- c(args, seqinfo = sinfo)
  }
  grg <- do.call(ImportGenomicRanges, args = args)

  if(! keep.file) file.remove(fpath)
  grg
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
  utils::download.file(url, lfp)
  system(paste("gunzip", lfp))
  lfp <- gsub("\\.gz$", "", lfp)
  chain <- import.chain(lfp)
  file.remove(lfp)
  chain
}
