# FUNCTIONS | KNN STATS ########################################################

# =============================================================================.
#
# -----------------------------------------------------------------------------.
knn_values <- function(v, i) {
  matrix(v[as.vector(i)], nrow(i), ncol(i))
}

# =============================================================================.
#' k-nearest neighbors standard deviation
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate observations
#'
#' @param nn.idx
#' matrix of nearest neighbor indexes
#'
#' @param nn.dist
#' matrix of nearest neighbor distances
#'
#' @param k
#' number of nearest neighbors
#'
#' @return knn_sd returns a vector of standard deviations
# -----------------------------------------------------------------------------.
#' @export
knn_sd <- function(x, nn.idx, nn.dist, k) {
  apply(nn.dist, MARGIN = 1, FUN = sd, na.rm = T)
}

# =============================================================================.
#' k-nearest neighbors estimator of density P(x) ~ k / N V
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate observations
#'
#' @param nn.idx
#' matrix of nearest neighbor indexes
#'
#' @param nn.dist
#' matrix of nearest neighbor distances
#'
#' @param k
#' number of nearest neighbors
#'
#' @return knn_density returns a vector of densities
# -----------------------------------------------------------------------------.
#' @export
knn_density <- function(x, nn.idx, nn.dist, k) {
  N  <- nrow(x)                 # total number of example values in the sample
  D  <- ncol(x)                 # number of dimensions
  V1 <- pi ^ (D / 2) / factorial(D / 2) # volume of the unit sphere in D dim.
  k / (N * V1 * nn.dist[, k] ^ D)       # standard density estimator
}

# =============================================================================.
#' k-nearest neighbor statistics
# -----------------------------------------------------------------------------.
# TODO: qry and sbj are not handled properly
# TODO: implement the case were qry and sbj are one-dimension (single samples)
# -----------------------------------------------------------------------------.
#' @param qry
#' numeric matrix representing multivariate observations
#'
#' @param k
#' number of nearest neighbors
#'
#' @param sbj
#'
#' @param FUN
#' statistical function (default = \link{knn_density})
#'
#' @return knn_stat returns a numeric vector
# -----------------------------------------------------------------------------.
#' @export
knn_stat <- function(qry, k, sbj = NULL, FUN = knn_density) {

  N  <- nrow(qry)                 # total number of example values in the sample
  idx <- which(finiteValues(qry)) # detect NA and Inf

  # find KNN
  if(is.null(sbj)) {
    # default: qry = sbj
    v <- get.knn(data = qry[idx, ], k = k, algorithm = "kd_tree")
  } else {
    # assume that sbj only contains usable values
    v <- get.knnx(data = sbj, query = qry[idx, ], k = k, algorithm = "kd_tree")
  }

  # compute stat
  d <- FUN(qry[idx, ], v$nn.index, v$nn.dist, k)
  d_all <- rep(NA, N)
  d_all[idx] <- d

  d
}
# =============================================================================.
#' knn statistics
# -----------------------------------------------------------------------------.
#'
#' @param data
#' numeric matrix representing multivariate observations
#'
#' @param k
#' number of nearest neighbors which corresponds to the smoothing parameter
#' of estimated statistics (larger values = smoother estimations).
#'
#' @param smoothing
#' default = T, recommended
#'
#' @return
#' cdadadr returns a \code{list} with the following elements:
#' \item{density}{knn density}
#' \item{variance}{knn variance}
# -----------------------------------------------------------------------------.
#' @export
knn_stats <- function(data, k, smoothing = T) {

  entropy <- function(p) { - sum(p * log(p)) }

  # get nearest neighbors
  r <- get.knn(data = data, k = k)

  # compute local density
  p <- knn_density(data, r$nn.index, r$nn.dist, k)
  if(smoothing) {
    p <- knn_values(p, r$nn.index)
    p <- apply(p, 1, mean)
  }
  p <- p / sum(p)

  # compute local variance
  m <- matrix(0, nrow(data), ncol(data))
  for(i in 1:ncol(data)) m[, i] <- rowMeans(knn_values(data[, i], r$nn.index))
  v <- rowMeans(knnx.dist(data = data, query = m, k = k)^2)

  # compute local entropy
  e <- knn_values(p, r$nn.index)
  e <- apply(e, 1, entropy)

  list(density = p, variance = v, entropy = e)
}
# =============================================================================.
#' Count Density After Dithering And Dimensionality Reduction
# -----------------------------------------------------------------------------.
#' @param cnt
#' matrix of read counts (row = observations, columns = samples or conditions).
#'
#' @param k
#' number of nearest neighbors which corresponds to the smoothing parameter
#' of estimated densities (larger values = smoother).
#'
#' @param dither
#' number of replicates of the dithered counts: 0 = no dithering,
#' 1 = single (default), 2 = duplicate, 3 = triplicate, etc.
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
#' @param smobs
#' subtract the mean value of each observation
#'
#' @param zscore
#'
#' @return
#' cdadadr returns a \code{list} with the following elements:
#' \item{density}{knn density}
#' \item{variance}{knn variance}
#' \item{entropy}{knn entropy}
#' \item{msg}{status of the dimensionality reduction}
#' \item{x}{transformed data}
# -----------------------------------------------------------------------------.
#' @export
cdadadr <- function(cnt, k, dither = 1, dred = 0.05, smobs = F, zscore = F) {

  r <- NULL

  pb <- txtProgressBar(min = 0, max = dither, char = "|", style = 3)
  for(i in 1:dither) {

    # Count dithering
    x <- log2(ditherCounts(cnt))

    # Centering
    if(centered) x <- x - rowMeans(x)

    # Dimensionality reduction
    msg <- "original dimensions"
    if(dred > 0) {
      # Compute Principal Components (PC)
      x <- prcomp(x, retx = T, center = T, scale. = T)

      # Retain PC conserving 100 * (1 - dred) % of the variance
      idx <- which(cumsum(x$sdev) / sum(x$sdev) > 1 - dred)
      if(length(idx) < 2) {
        idx <- 1:2
        warning("dimensionality reduction stringency (dred value) is too high")
      }
      msg <- paste("components =", length(idx))
      x <- x$x[, idx]
    }

    # z-score transformation
    if(zscore) x <- t((t(x) - colMeans(x)) / apply(x, MARGIN = 2, sd))

    # Density estimation
    l <- knn_stats(x, k = k)
    if(is.null(r)) {
      r <- l
    } else {
      for(id in names(r)) r[[id]] <- r[[id]] + l[[id]] / dither
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  c(r, list(msg = msg, x = x))
}
