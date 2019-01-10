# =============================================================================.
#' ** RESERVED FOR INTERNAL USE **
# -----------------------------------------------------------------------------.
#' @description
#' List of global options for the Tightrope package
# -----------------------------------------------------------------------------.
# from Barbouille
# FiniteValues RowSampler RankScore ASH2D
# ColorMapper ColorLegend ScatterMaps SideBySide
# -----------------------------------------------------------------------------.
#' @import methods
#' @import grDevices
#' @import graphics
#' @importFrom shades gradient
#' @importFrom stats median prcomp approx rnbinom runif
#' @importFrom utils download.file read.delim setTxtProgressBar txtProgressBar
#' @importFrom matrixStats rowMeans2 colMeans2 colSds
#' @importFrom triangle rtriangle
#' @importFrom ica icafast
#' @importFrom mixtools dmvnorm rmvnorm
#' @importFrom biomaRt useMart getBM
#' @import Rsamtools
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import GenomicFeatures
#' @importFrom QuickShift QuickShift
#' @importFrom Barbouille xylim FiniteValues RankScore ColorMapper ColorLegend ScatterMaps SideBySide
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
DefaultOptions <- function() {
  list(
    # BRD
    sampling  = NA,
    smobs     = TRUE,
    dither    = 5,
    zscore    = TRUE,
    bins      = 705,
    smoothing = 25,
    bdt       = c(0.15, 0.05),
    ncl       = 1,
    mincs     = 50,

    # PlotBRD
    palette    = "magma",
    gradient   = "hcl.duo.light",
    saturation = 1.0,

    messages  = TRUE
  )
}

# =============================================================================.
#' Global options for Tightrope functions
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{ScatterPlot},
#'   \link{SideBySide}
# -----------------------------------------------------------------------------.
#' @description
#' This function sets the default value of arguments used by the main functions
#' of the Tightrope package.
#'
#' @param ...
#' Any of the following arguments:
# -----------------------------------------------------------------------------.
#' @export
Tightrope <- function(...) {

  opt <- names(Tightrope::DefaultOptions())

  cfg <- list(...)
  cfg <- cfg[names(cfg) %in% opt]
  cfg <- cfg[! sapply(cfg, is.null)]

  if(length(cfg) > 0) {
    names(cfg) <- paste0("Tightrope.", names(cfg))
    options(cfg)
  } else {
    cfg <- options()[paste0("Tightrope.", opt)]
    if(any(is.na(names(cfg)))) {
      stop("missing global options")
    } else {
      names(cfg) <- gsub("^Tightrope\\.", "", names(cfg))
    }
    cfg
  }
}

# =============================================================================.
#' ** RESERVED FOR INTERNAL USE **
# -----------------------------------------------------------------------------.
#' @description
#' Remove global options of the Tightrope package from the R environment
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
RemoveOptions <- function() {
  cfg <- options()
  cfg <- cfg[grepl( "^Tightrope\\.", names(cfg))]
  cfg[] <- vector("list", length(cfg))
  options(cfg)
}

# =============================================================================.
#' ** RESERVED FOR INTERNAL USE **
# -----------------------------------------------------------------------------.
#' @description
#' Reinitialize global options of the Tightrope package
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ResetOptions <- function() {
  do.call(Tightrope, Tightrope::DefaultOptions())
}
