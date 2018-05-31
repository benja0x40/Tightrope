# =============================================================================.
#' Density estimation in reduced dimensions
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{BRD},
#'   \link{PlotBRD},
#'   \link{DitherCounts},
#'   \link{Barbouille::ASH2D}
# -----------------------------------------------------------------------------.
#' @description
#' Count Density after Dithering and Dimensionality Reduction (CDaDaDR).
#'
#' @param cnt
#' matrix of read counts
#' (rows = observations, columns = measurement conditions).
#'
#' @param movs
#' precomputed mean count values to be subtracted to each observation
#' (default = NULL, zero).
#'
#' @param bins
#' number of bins per principal component for density estimations
#' (default = 500).
#'
#' @param smoothing
#' number of consecutive bins for local average smoothing of estimated densities
#' (default = 10).
#'
#' @param dither
#' number of replicates of the count dithering performed by \link{DitherCounts}:
#' 1 = single, 2 = duplicate, 3 = triplicate, etc. (default = 5).
#'
#' @param zscore
#' transform read count projections into z-scores (default = TRUE, recommended).
#'
#' @param method
#' projection method, either "pca" (default) performed by \link{prcomp},
#' or "ica" performed by \link{icafast}.
#'
#' @return
#' \code{CDaDaDR.2D} returns a list with the following elements:
#' \item{parameters}{list with the value of each CDaDaDR.2D argument}
#' \item{status}{list with informations on the dimensionality reduction}
#' \item{density}{estimated density at each observation}
#' \item{projection}{transformed count matrix}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
CDaDaDR.2D <- function(
  cnt, movs = NULL, bins = 500, smoothing = 10, dither = 5, zscore = TRUE,
  method = c("pca", "ica"), progress = F
) {

  method <- match.arg(method)
  bins      <- rep(bins, length.out = 2)
  smoothing <- rep(smoothing, length.out = 2)

  parameters <- list(
    bins = bins, smoothing = smoothing, dither = dither, zscore = zscore
  )

  status <- list(
    dimensions = paste(method, "projection"),
    original   = ncol(cnt),
    reduced    = 2
  )

  if(is.null(movs)) movs <- 0

  p <- 0

  if(progress) pb <- utils::txtProgressBar(min = 0, max = dither, char = "|", style = 3)
  for(i in 1:dither) {

    # Count transformations
    x <- log2(DitherCounts(cnt)) # Dithering
    x <- x - movs                # Subtract mean of each observation  (if != 0)

    # Principal components (dimensionality reduction)
    if(method == "pca") x <- prcomp(x, retx = TRUE, center = TRUE, scale. = TRUE)$x[, 1:2]
    # Independent components (dimensionality reduction)
    if(method == "ica") x <- icafast(x, nc = 2, center = TRUE)$Y[, 1:2]

    colnames(x) <- paste0("C", 1:2)

    # z-score transformation
    if(zscore) x <- t((t(x) - matrixStats::colMeans2(x)) / matrixStats::colSds(x))

    # Density estimation
    p <- p + Barbouille::ASH2D(x, n = bins, k = smoothing) / dither

    if(progress) utils::setTxtProgressBar(pb, i)
  }
  if(progress) close(pb)

  list(
    parameters = parameters,
    status     = status,
    density    = p,
    projection = x
  )
}
