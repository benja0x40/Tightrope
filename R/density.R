# =============================================================================.
#' Naive mode estimation
# -----------------------------------------------------------------------------.
#' @param \code{x} \code{vector}
#' @param \code{na.rm} \code{logical}
#'
#' @return scalar
# -----------------------------------------------------------------------------.
findMode <- function(x, na.rm) {
  d <- density(x, na.rm=na.rm)
  m <- d$x[which.max(d$y)]
  m
}
# =============================================================================.
# TODO: densityByKernel
# -----------------------------------------------------------------------------.
densityByKernel <- function(qry, k, sbj=NULL, FUN=NULL) {
}
# =============================================================================.
#' Use k-nearest neighbors as an estimator of the density P(x) ~ k / N V
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
# -----------------------------------------------------------------------------.
#' @param qry
#' numeric matrix for which the empirical density will be estimated.
#'
#' @param k
#' index of the nearest neighbor used for local density estimation.
#'
#' @param sbj
#' numeric matrix providing the multivariate samples from which an empirical
#' density will be derived (default = qry).
#'
#' @param FUN
#'
#' @return \code{list}
# -----------------------------------------------------------------------------.
# TODO: implement the case were qry and sbj are one-dimension (single samples)
# TEST:
# x <- rnorm(100000, mean=c(0,5))
# y <- rnorm(100000, mean=c(0,5))
# d <- densityByKNN(cbind(x,y), k=100)
# plot(x, y, col=rgb(d$prc/100,0,0), pch=20, cex=0.5)
# NOTE: This illustrates nicely the curse of dimensionality and practically
# this KNN based density should be fine up to 16 dimensions. Beyond that, a
# simple way to circumvent the exponential increase of sparsity is to reduce
# dimensions (PCA or MDS) before estimating density.
# D <- 1:100
# V1 <- pi^(D/2)/factorial(D/2)
# layout(matrix(1:2,1,2))
# plot(D, V1)
# plot(D, V1, log='y')
# -----------------------------------------------------------------------------.
densityByKNN <- function(qry, k, sbj=NULL, FUN=NULL) {

  idx <- which(finiteValues(qry)) # detect NA and Inf

  N  <- nrow(qry)               # total number of example values in the sample
  D  <- ncol(qry)               # number of dimensions
  V1 <- pi^(D/2)/factorial(D/2) # volume of the unit sphere in D dimensions

  # find KNN
  if(is.null(sbj)) {
    # default: qry = sbj
    v <- get.knnx(data=qry[idx,], query=qry[idx,], k=k, algorithm="kd_tree")
  } else {
    # assume that sbj only contains usable values
    v <- get.knnx(data=sbj, query=qry[idx,], k=k, algorithm="kd_tree")
  }

  # Compute density
  if(! is.null(FUN)) {
    d <- apply(v$nn.dist, MARGIN=1, FUN=FUN) # user provided estimator function
  } else {
    d <- k/(N * V1 * v$nn.dist[,k]^D)        # standard estimator
  }
  d_all <- rep(NA, N)
  d_all[idx] <- d

  # Find example values with minimum and maximum density
  idx.min <- idx[which.min(d)]
  idx.max <- idx[which.max(d)]

  # Compute percentiles of the density distribution
  q <- quantile(d, probs=0:100/100)

  # Compute the percentile of density corresponding to each example values
  prc_all <- rep(NA, N)
  prc_all[idx] <- 100*rank(d)/length(d)
  # Result
  list(d=d_all, q=q, prc=prc_all, idx.min=idx.min, idx.max=idx.max)
}
