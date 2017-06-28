# =============================================================================.
#' Count Density after Dithering and Dimensionality Reduction
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
#' @param movs
#' precomputed mean count values to be subtracted to each observation
#' (default = NULL, zero).
#'
#' @param dither
#' number of replicates for count dithering:
#' 1 = single (default), 2 = duplicate, 3 = triplicate, etc.
#'
#' @param npc
#' number of dimensions retained after projection on principal components.
#' When \code{npc} is specified, the \code{cvt} parameter is ignored.
#'
#' @param zscore
#' transform counts into z-scores (logical, default = T).
#'
#' @param rare
#' defines exceptionally low density values (default = 0.01).
#'
#' @param method
#' projection method, either "pca" (default) or "ica".
#'
#' @param progress
#' show progress bar (logical, default = F).
#'
#' @return
#' CDaDaDR returns a \code{list} with the following elements:
#' \item{parameters}{list with the value of each CDaDaDR argument}
#' \item{status}{list with informations on the dimensionality reduction}
#' \item{density}{knn density}
#' \item{projection}{transformed count matrix}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
CDaDaDR <- function(
  cnt, knn, smobs = F, movs = NULL,
  dither = 1, npc = NA, zscore = T, rare = 0.01,
  method = c("pca", "ica"),
  progress = F
) {

  parameters <- list(
    knn = knn, smobs = smobs,
    dither = dither, npc = npc, zscore = zscore, rare = rare
  )
  status <- list(
    dimensions = "initial measurements",
    original   = ncol(cnt),
    reduced    = NA,
    sdpc       = NA, # standard deviation per component
    cvpc       = NA  # cumulated percentage of variance per component
  )

  if(progress) pb <- txtProgressBar(min = 0, max = dither, char = "|", style = 3)

  method <- method[1]

  if(smobs & is.null(movs)) {
    movs <- knn_mean(log2(cnt), k = 5) # ???
    movs <- rowMeans(movs)
  }

  p <- 0

  for(i in 1:dither) {

    # Count transformations
    x <- log2(DitherCounts(cnt)) # Dithering
    if(smobs) x <- x - movs      # Subtract mean of each observation  (if != 0)

    # Dimensionality reduction
    if(! is.na(npc)) {

      status$dimensions <- paste(method, "projection")
      if(npc < 2) {
        npc <- 2
        warning("minimum npc value is 2")
      }
      idx <- 1:npc

      # Principal components
      if(method == "pca") {
        x <- prcomp(x, retx = T, center = T, scale. = T)
        status$sdpc <- x$sdev
        names(status$sdpc) <- paste0("C", 1:ncol(cnt))
        x <- x$x[, idx]
      }
      # Independent components
      if(method == "ica") {
        x <- icafast(x, nc = npc, center = T)
        status$sdpc <- x$vafs
        names(status$sdpc) <- paste0("C", 1:npc)
        x <- x$Y[, idx]
      }
      colnames(x) <- names(status$sdpc)[1:npc]

      status$cvpc <- with(status, 1 - cumsum(sdpc / sum(sdpc)))
      names(status$cvpc) <- names(status$sdpc)

      status$reduced <- length(idx)
    }

    # z-score transformation
    if(zscore) x <- t((t(x) - colMeans(x)) / apply(x, MARGIN = 2, sd))

    # Density estimation
    p <- p + knn_density(x, k = knn) / dither

    if(progress) setTxtProgressBar(pb, i)
  }
  if(progress) close(pb)

  rare <- p < mean(p[rankstat(p) < rare])

  list(
    parameters = parameters,
    status     = status,
    density    = p,
    rare       = rare,
    projection = x
  )
}
