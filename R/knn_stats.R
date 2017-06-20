# FUNCTIONS | KNN STATS ########################################################

# # =============================================================================.
# #' knn in 1D (not so fast, NOT FULLY WORKING)
# # -----------------------------------------------------------------------------.
# #' @param data
# #' numeric matrix representing multivariate data where rows = observations
# #' and columns = measurement conditions.
# #'
# #' @param k
# #' number of nearest neighbors.
# #'
# #' @return
# # -----------------------------------------------------------------------------.
# #' @keywords internal
# #' @export
# get.knn_1D <- function(data, k) {
#
#   N  <- length(data) # number of observations
#
#   o <- order(data)
#   h <- data[o]
#
#   d  <- matrix(0, N, k)
#   da.1 <- c(h[2:N] - h[1:(N - 1)], Inf)
#   db.1 <- c(Inf, h[2:N] - h[1:(N - 1)])
#
#   da <- da.1
#   db <- db.1
#
#   i  <- matrix(0, N, k)
#   ia <- 1:N
#   ib <- 1:N
#
#   for(j in 1:k) {
#     # find nn either below or above
#     ma <- which(da <  db)
#     mb <- which(db <= da)
#     # update nn index
#     ia[ma] <- ia[ma] + 1
#     ib[mb] <- ib[mb] - 1
#     # save nn index
#     i[ma, j] <- ia[ma]
#     i[mb, j] <- ib[mb]
#     # save nn distance
#     d[ma, j] <- da[ma]
#     d[mb, j] <- db[mb]
#     # update nn distances
#     da[ma] <- da[ma] + da.1[ia[ma]]
#     db[mb] <- db[mb] + db.1[ib[mb]]
#   }
#   i[o, ] <- i
#   d[o, ] <- d
#
#   list(nn.index = i, nn.dist = d)
# }

# =============================================================================.
#' Extract local values from knn index matrix
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
#' Apply smoothing function to local values
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
#' knn density estimator
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
#' and columns = measurement conditions.
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

  x <- as.matrix(x)

  N  <- nrow(x) # number of observations
  D  <- ncol(x) # number of dimensions of each observation

  if(is.null(i) | is.null(d)) {
    if(D == 1) r <- get.knn(data = x, k = k)
    else       r <- get.knn(data = x, k = k)

    i <- r$nn.index
    d <- r$nn.dist
  }

  # Compute density in any dimensions D
  V1 <- pi ^ (D / 2) / factorial(D / 2) # volume of the unit sphere
  p <- k / (N * V1 * d[, k] ^ D)        # standard knn density estimator

  if(smoothing) {
    p <- knn_smoothing(p, i)
  }

  p <- p / sum(p)
}

# =============================================================================.
#' Local mean and variance
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{knn_density}
# -----------------------------------------------------------------------------.
#' @description
#' k-nearest neighbor (knn) centroid and corresponding standard deviation.
#'
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = measurement conditions.
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

  x <- as.matrix(x)

  N  <- nrow(x) # number of observations
  D  <- ncol(x) # number of dimensions of each observation

  if(is.null(i) | is.null(d)) {
    if(D == 1) r <- get.knn(data = x, k = k)
    else       r <- get.knn(data = x, k = k)

    i <- r$nn.index
    d <- r$nn.dist
  }

  # compute centroids
  m <- matrix(0, N, D)
  for(j in 1:D) {
    m[, j] <- rowMeans(knn_values(x[, j], i))
  }

  # compute average distances to centroids
  v <- rowMeans(knnx.dist(data = x, query = m, k = k)^2)

  if(smoothing) {
    v <- knn_smoothing(v, i)
  }

  list(mu = m, sigma2 = v)
}
