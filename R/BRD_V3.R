# =============================================================================.
#' Convert (x, y) to (a = mean, m = difference) coordinates
# -----------------------------------------------------------------------------.
#' @description
#' compute \eqn{a = (y + x)/2}, and \eqn{m = (y - x)}
#'
#' @param x
#' numeric vector, or a matrix with 2 columns or a list containing two vectors.
#'
#' @param y
#' numeric vector.
#'
#' @return
#' \code{xy2md} returns a list or matrix with (a = mean, m = difference) values.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
xy2md <- function(x, y = NULL) {

  chk  <- 0

  if(is.null(y) & ! is.null(dim(x))) {
    y <- x[,2]
    x <- x[,1]
    chk <- 1
  }
  if(is.null(y) & is.list(x)) {
    y <- x[[2]]
    x <- x[[1]]
  }

  a <- (y + x)/2
  m <- (y - x)

  if(chk == 0) r <- list(a = a, m = m)
  if(chk == 1) r <- cbind(a = a, m = m)

  r
}

# =============================================================================.
#' Average ranking for control measurements
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{BRD}
# -----------------------------------------------------------------------------.
#' @description
#' Average ranking for control measurements
#' (e.g. input DNA, ChIP-seq performed with IgG or in KO conditions).
#'
#' @param x
#' matrix of read counts
#' (rows = observations, columns = measurement conditions).
#'
#' @param by
#' grouping index.
#'
#' @param controls
#' columns corresponding to control measurements (Input, IgG, etc.)
#'
#' @return
#' \code{bg_ranking} returns group identifiers ordered by decreasing background
#' intensity, meaning the average difference between control enrichment
#' conditions.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
bg_ranking <- function(x, by, controls) {

  bg <- split(as.data.frame(as.matrix(x)), f = by)
  ng <- length(bg)
  z <- t(sapply(bg, colMeans))
  a <- matrixStats::rowMeans2(as.matrix(z[, - controls]))
  b <- matrixStats::rowMeans2(as.matrix(z[, controls]))
  bg <- (1:ng)[order(b - a, decreasing = TRUE)]

  bg
}

# =============================================================================.
#' Background Reads Density
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{PlotBRD},
#'   \link{BackgroundCandidates},
#'   \link{ScalingFactors},
#'   \link{NormalizeCountMatrix},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @description
#' The \code{BRD} function searches for background candidates based on the
#' provided read count matrix and uses these background candidates to estimate
#' normalization factors.
#'
#' @inheritParams CDaDaDR.2D
#'
#' @param controls
#' count columns corresponding to control measurements (Input, IgG, etc.).
#'
#' @param smobs
#' subtract a mean count value for each observation (default = TRUE, recommended).
#'
#' @param bdt
#' numeric vector of length 2 defining background density thresholds both
#' expressed as proportions between 0 and 1.
#' \code{bdt[1]} specifies the global threshold used to discard observations
#' with low density prior to clustering.
#' \code{bdt[2]} determines the maximum density loss allowed
#' when selecting core observations relatively to the local maximum density
#' in each cluster, and thus defining the background candidates.
#' By default the value of \code{bdt} is \code{c(0.2, 0.05)} meaning that,
#' in terms of density percentiles, the bottom 20 percents will be filtered out
#' before clustering and only the top 5 percents can be selected as background
#' candidates among each cluster.
#'
#' @param ncl
#' number of clusters (density modes) to be distinguished.
#'
#' @param mincs
#' minimum size of cluster cores, as number of observations.
#'
#' @return
#' \code{BRD} returns a list with the following elements:
#' \item{parameters}{call parameters of the function.}
#' \item{status}{execution status.}
#' \item{nonzero}{indices of initial observations with count > 0.}
#' \item{dred}{
#'   dimensionality-reduced non-zero observations.
#' }
#' \item{subsets}{
#'   partition of non-zero observations into background candidate subsets.
#' }
#' \item{populations}{
#'   summary of the core population in each subset.
#' }
#' \item{theta}{
#'   fitted distribution parameters for each core population.
#' }
#' \item{log2counts}{
#'   dithered and log2 transformed counts.
#' }
#' \item{normfactors}{
#'   BRD normalization factors.
#' }
# -----------------------------------------------------------------------------.
#' @export
BRD <- function(
  cnt, controls, sampling = NULL, smobs = NULL, dither = NULL, zscore = NULL,
  bins = NULL, smoothing = NULL, bdt = NULL, ncl = NULL, mincs = NULL
) {

  if(is.null(colnames(cnt))) stop("missing column names in the count matrix")
  if(length(controls) < 2) stop("BRD requires at least 2 controls")

  cfg <- Tightrope() # Global options
  DefaultArgs(cfg, ignore = c("cnt", "controls"), from = BRD)


  xps <- colnames(cnt)
  inp <- match(controls, xps)

  parameters <- list(
    experiments = xps, controls = xps[inp], smobs = smobs,
    bins = bins, smoothing = smoothing, dither = dither, zscore = zscore,
    bdt = bdt, ncl = ncl, mincs = mincs
  )

  # Sampling (optional) and cleanup of missing values
  cleanup <- TRUE
  if(! is.null(sampling)) {
    cleanup <- rep(F, nrow(cnt))
    cleanup[RowSampler(cnt, max = sampling)] <- TRUE
  }
  cleanup <- cleanup & FiniteValues(log2(cnt))
  ldc <- log2(DitherCounts(cnt))
  ldc <- ldc[cleanup, ]     # Dithered and log2 transformed counts
  cnt <- cnt[cleanup, ]     # Raw counts

  ldc <- RemoveInputBias(ldc, controls = inp, as.log2 = TRUE, safe = TRUE)

  # Precompute the mean of each observation in control and enrichment conditions
  movs <- matrixStats::rowMeans2(as.matrix(log2(cnt[,   inp, drop = FALSE])))
  menr <- matrixStats::rowMeans2(as.matrix(log2(cnt[, - inp, drop = FALSE])))
  intensity <- xy2md(menr, movs)$m

  # Compute count density after dithering and dimensionality reduction
  dred <- CDaDaDR.2D(
    cnt, movs = movs,
    bins = bins, smoothing = smoothing, dither = dither, zscore = zscore,
    method = "pca"
  )
  dred <- c(dred, list(ctl = movs, enr = menr, intensity = intensity))

  # Filter out observations with lower density than specified
  rd <- RankScore(dred$density)
  sbs <- which(rd > bdt[1])
  sbs <- with(
    dred, list(
      i = sbs,
      b = intensity[sbs], d = density[sbs], r = rd[sbs], x = projection[sbs, ]
    )
  )

  # Find clusters by density gradient ascent in the dred space
  qsc <- with(sbs, QuickShiftClustering(x, d, n = ncl))
  sbs$cluster <- qsc$membership

  # Select top density (core) observations in each cluster and count members
  pop <- data.frame(cluster = rep(0, ncl), core = rep(0, ncl))
  pop$max_density <- as.vector(by(sbs$r, sbs$cluster, max))
  sbs$core  <- FALSE
  for(grp in 1:ncl) {
    k <- with(sbs, cluster == grp)
    sbs$core[k] <- sbs$r[k] > max(sbs$r[k]) - bdt[2]
    pop$cluster[grp] <- sum(k)
    pop$core[grp]    <- sum(sbs$core[k])
  }

  # Rank cluster cores by background levels (enrichment in control samples)
  if(ncl > 1) {
    bg_rnk <- with(
      sbs, bg_ranking(ldc[i[core], ], by = cluster[core], controls = inp)
    )
  } else {
    bg_rnk <- 1
  }
  bg_clu <- bg_rnk[1]
  bg_idx <- with(sbs, i[core & cluster == bg_clu])
  pop$background <- 1:ncl == bg_clu

  # Fit a multivariate gaussian model to each cluster core of minimum size
  theta <- list()
  for(grp in 1:ncl) {
    idx <- with(sbs, i[core & cluster == grp]) # select cluster core
    dns <- with(sbs, d[core & cluster == grp]) # core density
    if(length(idx) > mincs) {
      # theta[[grp]] <- mv_mle(ldc[idx, ], w = dns)
      theta[[grp]] <- list(mu = 0, sigma = 0)
      for(i in 1:dither) {
        tst <- log2(DitherCounts(cnt[idx, ]))
        tst <- RemoveInputBias(tst, controls = inp, as.log2 = TRUE, safe = TRUE)
        tst <- mv_mle(tst, w = dns)
        theta[[grp]]$mu    <- theta[[grp]]$mu + tst$mu / dither
        theta[[grp]]$sigma <- theta[[grp]]$sigma + tst$sigma / dither
      }
    } else {
      theta[[grp]] <- NA
    }
  }

  res <- list(
    parameters  = parameters,
    status      = "clustering",
    nonzero     = which(cleanup), # observation with count > 0
    dred        = dred,           # CDaDaDR results
    subsets     = sbs,            # clustering of non-zero observations
    populations = pop,            # summary of cluster/core populations
    theta       = theta,          # distribution parameters for each cluster
    log2counts  = ldc             # dithered and log2 transformed counts
  )

  # Background cluster
  res$status <- "successful selection fo background candidates"
  bg_theta <- theta[[bg_clu]]
  if(length(bg_theta) > 1) {
    res <- c(
      res, list(
        normfactors = mean(bg_theta$mu) - bg_theta$mu,
        bg_clurank  = bg_rnk,
        bg_cluster  = bg_clu,
        bg_theta    = bg_theta,
        bg_members  = bg_idx
      )
    )
  } else {
    res$status <- "selection of background candidates has failed"
    warning(res$status)
  }

  res
}

# =============================================================================.
#' Index of background candidates
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{BRD},
#'   \link{PlotBRD},
#'   \link{PlotCountDistributions}
# -----------------------------------------------------------------------------.
#' @description
#' Provides the index of background candidates within rows of the read count
#' matrix used to identify these candidates.
#'
#' @param brd
#' result of a prior call to the \link{BRD} function.
#'
#' @return
#' \code{BackgroundCandidates} returns an integer vector.
# -----------------------------------------------------------------------------.
#' @export
BackgroundCandidates <- function(brd) {
  brd$bg_members
}

# =============================================================================.
#' Get BRD scaling factors
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{BRD},
#'   \link{NormalizeCountMatrix}
# -----------------------------------------------------------------------------.
#' @description
#' Extract the vector of scaling factors from the result of a prior call to the
#' \link{BRD} function.
#'
#' @param brd
#' result of a prior call to the \link{BRD} function.
#'
#' @param as.log2
#' logical (default = FALSE, no).
#'
#' @return
#' \code{ScalingFactors} returns a numeric vector.
# -----------------------------------------------------------------------------.
#' @export
ScalingFactors <- function(brd, as.log2 = FALSE) {

  # Retrieve normalization factors
  f <- brd$normfactors
  if(is.null(f)) warning("invalid argument, this function requires the object resulting from the BRD function as argument")

  # Revert log2 transformation if needed
  if(! as.log2) f <- 2^f

  f
}

# =============================================================================.
#' Apply BRD scaling factors
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{BRD},
#'   \link{ScalingFactors},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @description
#' Apply precomputed BRD normalization factors to a read count matrix.
#' Use \code{as.log2 = TRUE} if you want to provide a matrix of log2
#' transformed read counts instead of raw read counts.
#'
#' @inheritParams ScalingFactors
#'
#' @param x
#' matrix of read counts (rows = observations, columns = samples or conditions).
#'
#' @return
#' \code{NormalizeCountMatrix} returns a \code{matrix}.
# -----------------------------------------------------------------------------.
#' @export
NormalizeCountMatrix <- function(x, brd, as.log2 = FALSE) {

  if(any(x < 0) & ! as.log2) {
    message("use as.log2 = TRUE to normalize log2 transformed read counts")
    stop("the read count matrix has negative values")
  }

  # Retrieve normalization factors
  f <- brd$normfactors
  if(is.null(f)) {
    message("provide a result from calling BRD function as argument")
    stop("invalid brd object")
  }
  f <- f[colnames(x)]

  if(as.log2) { # Normalization of log2 transformed counts = t(t(x) + f)
    x <- AddByRow(x, f)
  } else { # Normalization of raw counts = t(t(x) * 2^f)
    x <- MulByRow(x, 2^f)
  }

  x
}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
#' @export
RemoveInputBias <- function(X, controls, as.log2 = FALSE, safe = FALSE) {
  if(! as.log2) X <- log2(X)
  if(! safe) X[X == -Inf] <- NA
  bias <- X[, controls]
  # TODO: could be improved by adding GC content and mappability
  bias <- matrixStats::rowMeans2(bias, na.rm = T) - median(bias, na.rm = T)
  X <- X - bias
  if(! as.log2) X <- 2^X
  X
}
