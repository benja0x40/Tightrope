# FUNCTIONS | KNN STATS ########################################################

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
