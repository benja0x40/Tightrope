# =============================================================================.
#' Count Density After Dithering And Dimensionality Reduction
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{BRD},
#'   \link{DitherCounts},
#'   \link{knn_density}
# -----------------------------------------------------------------------------.
#' @param cnt
#' matrix of read counts
#' (rows = observations, columns = measurement conditions).
#'
#' @param knn
#' number of nearest neighbors, which corresponds to the smoothing parameter
#' of estimated densities (larger values = smoother).
#'
#' @param smobs
#' subtract the mean count value of each observation (logical, default = F).
#'
#' @param cvt
#' percentage of variance to be ignored.
#' The value of \code{cvt} determines the stringency of the dimensionality
#' reduction which is performed by principal component analysis using
#' \link{prcomp}.
#' Dimensions retained after projection of the read count dataset on its
#' principal components will preserve \code{100 * (1 - cvt)} percents
#' of the initial variance.
#' For instance, after dimensionality reduction with the default setting
#' (\code{cvt} = 0.5), approximately 50 percents of the initial variance
#' will remain in the reduced data.
#' Setting \code{cvt} = 0 bypasses the dimensionality reduction and performs
#' density estimations in the original multivariate measurement space
#' of the dataset.
#'
#' @param npc
#' number of dimensions retained after projection on principal components.
#' When \code{npc} is specified, the \code{cvt} parameter is ignored.
#'
#' @param zscore
#' transform counts into z-scores (logical, default = T).
#'
#' @param dither
#' number of replicates for count dithering:
#' 1 = single (default), 2 = duplicate, 3 = triplicate, etc.
#'
#' @param progress
#' show progress bar (logical, default = F).
#'
#' @return
#' cdadadr returns a \code{list} with the following elements:
#' \item{parameters}{list with the value of each cdadadr argument}
#' \item{status}{list with informations on the dimensionality reduction}
#' \item{density}{knn density}
#' \item{projection}{transformed count matrix}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
cdadadr <- function(
  cnt, knn, smobs = F, cvt = 0.5, npc = NA, zscore = T, dither = 1, progress = F
) {

  parameters <- list(
    knn = knn, smobs = smobs, cvt = cvt, npc = npc, zscore = zscore, dither = dither
  )
  status <- list(
    dimensions = "initial measurements",
    original   = ncol(cnt),
    reduced    = NA,
    dimcut     = "",
    pcsd       = NA, # standard deviation of principal components
    cvpc       = NA  # cumulated percentage of variance
  )

  if(progress) pb <- txtProgressBar(min = 0, max = dither, char = "|", style = 3)

  p <- 0

  for(i in 1:dither) {

    # Count dithering
    x <- log2(DitherCounts(cnt))

    # Subtract mean of each observation
    if(smobs) x <- x - rowMeans(x)

    # Dimensionality reduction
    if(cvt > 0 | ! is.na(npc)) {
      status$dimensions <- "principal components"

      if(! is.na(npc)) {
        # Retain PC as directly specified
        status$dimcut <- paste("fixed number of PC =", npc)
        if(npc < 2) {
          npc <- 2
          status$dimcut <- "inapplicable PC number threshold"
        }
        idx <- 1:npc
      }

      # Principal components
      x <- prcomp(x, retx = T, center = T, scale. = T)

      status$pcsd <- x$sdev
      status$cvpc <- 1 - cumsum(x$sdev / sum(x$sdev))
      names(status$pcsd) <- paste0("PC", 1:ncol(cnt))
      names(status$cvpc) <- names(status$pcsd)

      if(is.na(npc)) {
        # Retain PC conserving 100 * (1 - cvt) % of the variance
        status$dimcut <- paste("variance threshold = ", cvt)
        idx <- which(status$cvpc > cvt)
        if(length(idx) < 2) {
          idx <- 1:2
          status$dimcut <- "inapplicable component variance threshold"
        }
      }
      status$reduced <- length(idx)
      x <- x$x[, idx]
    }

    # z-score transformation
    if(zscore) x <- t((t(x) - colMeans(x)) / apply(x, MARGIN = 2, sd))

    # Density estimation
    p <- p + knn_density(x, k = knn) / dither

    if(progress) setTxtProgressBar(pb, i)
  }
  if(progress) close(pb)

  list(
    parameters = parameters,
    status     = status,
    density    = p,
    projection = x
  )
}
