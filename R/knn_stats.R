# FUNCTIONS | KNN STATS ########################################################

# # =============================================================================.
# #' k-nearest neighbor statistics - DEPRECATED
# # -----------------------------------------------------------------------------.
# # TODO: qry and sbj are not handled properly
# # TODO: implement the case were qry and sbj are one-dimension (single samples)
# # -----------------------------------------------------------------------------.
# #' @param qry
# #' numeric matrix representing multivariate observations
# #'
# #' @param k
# #' number of nearest neighbors
# #'
# #' @param sbj
# #'
# #' @param FUN
# #' statistical function (default = \link{knn_density})
# #'
# #' @return knn_stat returns a numeric vector
# # -----------------------------------------------------------------------------.
# #' @keywords internal
# #' @export
# knn_stat <- function(qry, k, sbj = NULL, FUN = knn_density) {
#
#   N  <- nrow(qry)                 # total number of example values in the sample
#   idx <- which(FiniteValues(qry)) # detect NA and Inf
#
#   # find KNN
#   if(is.null(sbj)) {
#     # default: qry = sbj
#     v <- get.knn(data = qry[idx, ], k = k, algorithm = "kd_tree")
#   } else {
#     # assume that sbj only contains usable values
#     v <- get.knnx(data = sbj, query = qry[idx, ], k = k, algorithm = "kd_tree")
#   }
#
#   # compute stat
#   d <- FUN(qry[idx, ], v$nn.index, v$nn.dist, k)
#   d_all <- rep(NA, N)
#   d_all[idx] <- d
#
#   d
# }

# # =============================================================================.
# #' knn statistics - DEPRECATED
# # -----------------------------------------------------------------------------.
# #' @seealso
# #'   \link{knn_density},
# #'   \link{knn_musigma2}
# # -----------------------------------------------------------------------------.
# #' @param data
# #' numeric matrix representing multivariate observations
# #'
# #' @param k
# #' number of nearest neighbors which corresponds to the smoothing parameter
# #' of estimated statistics (larger values = smoother estimations).
# #'
# #' @return
# #' knn_stats returns a \code{list} with the following elements:
# #' \item{density}{knn density estimator}
# #' \item{entropy}{knn entropies}
# #' and additional elements resulting from \link{knn_musigma2}
# # -----------------------------------------------------------------------------.
# #' @keywords internal
# #' @export
# knn_stats <- function(data, k) {
#
#   entropy <- function(p) { - sum(p * log(p)) }
#
#   # get nearest neighbors
#   r <- get.knn(data = data, k = k)
#
#   # compute local density
#   p <- knn_density(data, k, d = r$nn.dist)
#
#   # compute local centroids and corresponding variance
#   musigma2 <- knn_musigma2(data, k, i = r$nn.index, d = r$nn.dist)
#
#   # compute local entropy
#   e <- knn_values(p, r$nn.index)
#   e <- apply(e, 1, entropy)
#
#   c(list(density = p, entropy = e), musigma2)
# }

# =============================================================================.
#' knn_values
# -----------------------------------------------------------------------------.
#' @param v
#' numeric vector.
#'
#' @param i
#' precomputed matrix of nearest neighbor indices.
#'
#' @return
#' knn_values returns a matrix.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
knn_values <- function(v, i) {
  matrix(v[as.vector(i)], nrow(i), ncol(i))
}

# =============================================================================.
#' knn_smoothing
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{knn_density},
#'   \link{knn_musigma2}
# -----------------------------------------------------------------------------.
#' @param v
#' numeric vector.
#'
#' @param i
#' precomputed matrix of nearest neighbor indices.
#'
#' @param f
#' smoothing function (default = mean).
#'
#' @return
#' knn_smoothing returns a numeric vector.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
knn_smoothing <- function(v, i, f = mean) {
  v <- knn_values(v, i)
  v <- apply(v, 1, f)
  v
}

# =============================================================================.
#' knn_density
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{knn_musigma2}
# -----------------------------------------------------------------------------.
#' @description
#' k-nearest neighbor (knn) estimator of the density \deqn{P(Xi) ~ k / N Vi}
#' where \eqn{N} is the number of observations and \eqn{Vi} the volume of a
#' sphere of radius equal to the distance between \eqn{Xi} and its k-nearest
#' neighbor.
#'
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param k
#' number of nearest neighbors, which corresponds to the smoothing parameter
#' of estimated densities (larger values = smoother).
#'
#' @param i
#' precomputed matrix of nearest neighbor indices (optional).
#'
#' @param d
#' precomputed matrix of nearest neighbor distances (optional).
#'
#' @param smoothing
#' perfom a knn average smoothing of the density (default = T, recommended).
#'
#' @return
#' knn_density returns a numeric vector.
# -----------------------------------------------------------------------------.
#' @export
knn_density <- function(x, k, i = NULL, d = NULL, smoothing = T) {

  if(is.null(i) | is.null(d)) {
    r <- get.knn(data = x, k = k)
    i <- r$nn.index
    d <- r$nn.dist
  }

  # compute density in D = ncol(x) dimensions
  N  <- nrow(x)                 # total number of observations in the sample
  D  <- ncol(x)                 # number of dimensions of each observation
  V1 <- pi ^ (D / 2) / factorial(D / 2) # volume of the unit sphere
  p <- k / (N * V1 * d[, k] ^ D)        # standard knn density estimator

  if(smoothing) {
    p <- knn_smoothing(p, i)
  }

  p <- p / sum(p)
}

# =============================================================================.
#' knn_musigma2
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{knn_density}
# -----------------------------------------------------------------------------.
#' @description
#' k-nearest neighbor (knn) centroid and corresponding standard deviation.
#'
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param k
#' number of nearest neighbors (i.e. smoothing factor).
#'
#' @param i
#' precomputed matrix of nearest neighbor indices (optional).
#'
#' @param d
#' precomputed matrix of nearest neighbor distances (optional).
#'
#' @param smoothing
#' perfom a knn average smoothing of the standard deviations (default = F).
#'
#' @return
#' knn_musigma2 returns a \code{list} with the following elements:
#' \item{mu}{knn centroids}
#' \item{sigma2}{knn variances (squared standard deviations)}
# -----------------------------------------------------------------------------.
#' @export
knn_musigma2 <- function(x, k, i = NULL, d = NULL, smoothing = F) {

  if(is.null(i) | is.null(d)) {
    r <- get.knn(data = x, k = k)
    i <- r$nn.index
    d <- r$nn.dist
  }

  # compute centroids
  m <- matrix(0, nrow(x), ncol(x))
  for(j in 1:ncol(x)) {
    m[, j] <- rowMeans(knn_values(x[, j], i))
  }

  # compute average distances to centroids
  v <- rowMeans(knnx.dist(data = x, query = m, k = k)^2)

  if(smoothing) {
    v <- knn_smoothing(v, i)
  }

  list(mu = m, sigma2 = v)
}

# =============================================================================.
#' Count Density After Dithering And Dimensionality Reduction
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{DitherCounts},
#'   \link{prcomp},
#'   \link{knn_density}
# -----------------------------------------------------------------------------.
#' @param cnt
#' matrix of read counts (rows = observations, columns = samples or conditions).
#'
#' @param k
#' number of nearest neighbors, which corresponds to the smoothing parameter
#' of estimated densities (larger values = smoother).
#'
#' @param smobs
#' subtract the mean count value of each observation (logical, default = F).
#'
#' @param dred
#' stringency of the dimensionality reduction which is performed by
#' Principal Component Analysis. Dimensions retained after projection of the
#' read count dataset on its principal components will preserve
#' 100 * (1 - dred) % of the total variance.
#' For instance with the default value (dred = 0.05), 95% of the total variance
#' will be retained after dimensionality reduction.
#' Setting dred = 0 bypasses the dimensionality reduction and performs
#' density estimations in the original multivariate space of the dataset.
#'
#' @param npc
#' number of dimensions retained after projection on principal components.
#' When \code{npc} is specified, the \code{dred} parameter is ignored.
#'
#' @param zscore
#' transform counts into z-scores (logical, default = T).
#'
#' @param dither
#' number of replicates of the dithered counts: 0 = no dithering,
#' 1 = single (default), 2 = duplicate, 3 = triplicate, etc.
#'
#' @param progress
#' show progress bar (logical, default = F).
#'
#' @return
#' cdadadr returns a \code{list} with the following elements:
#' \item{density}{knn density}
#' \item{parameters}{list with the value of each cdadadr argument}
#' \item{projection}{list}
# -----------------------------------------------------------------------------.
#' @export
cdadadr <- function(
  cnt, k, smobs = F, dred = 0.05, npc = NA, zscore = T, dither = 1, progress = F
) {

  p <- 0

  msg <- list(
    status   = "original dimensions",
    original = ncol(cnt),
    reduced  = NA,
    dred     = ""
  )

  if(progress) pb <- txtProgressBar(min = 0, max = dither, char = "|", style = 3)

  for(i in 1:dither) {

    # Count dithering
    x <- log2(DitherCounts(cnt))

    # Subtract mean of each observation
    if(smobs) x <- x - rowMeans(x)

    # Dimensionality reduction
    if(dred > 0 | ! is.na(npc)) {
      msg$status <- "principal components"

      if(! is.na(npc)) {
        # Retain PC as directly specified
        msg$dred <- "fixed number of PC"
        if(npc < 2) {
          npc <- 2
          msg$dred <- "irrelevant number of PC"
        }
        idx <- 1:npc
      }

      # Compute Principal Components (PC)
      x <- prcomp(x, retx = T, center = T, scale. = T)

      if(is.na(npc)) {
        # Retain PC conserving 100 * (1 - dred) % of the variance
        msg$dred <- "variance threshold"
        a <- 1 - cumsum(x$sdev / sum(x$sdev))
        idx <- which(a > dred)
        if(length(idx) < 2) {
          idx <- 1:2
          msg$dred <- "irrelevant variance threshold"
        }
      }
      msg$reduced <- length(idx)
      x <- x$x[, idx]
    }

    # z-score transformation
    if(zscore) x <- t((t(x) - colMeans(x)) / apply(x, MARGIN = 2, sd))

    # Density estimation
    p <- p + knn_density(x, k = k) / dither

    if(progress) setTxtProgressBar(pb, i)
  }
  if(progress) close(pb)

  list(
    density = p,
    parameters = list(
      k = k, dither = dither, dred = dred, smobs = smobs, zscore = zscore
    ),
    projection = list(msg = msg, x = x)
  )
}
