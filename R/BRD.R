# =============================================================================.
#' Average ranking for control measurements
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{BRD}
# -----------------------------------------------------------------------------.
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
#' group identifiers ordered by decreasing average difference
#' between controls and other measurement conditions.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
bg_ranking <- function(x, by, controls) {

  bg <- split(as.data.frame(as.matrix(x)), f = by)
  ng <- length(bg)
  z <- t(sapply(bg, colMeans))
  a <- rowMeans(as.matrix(z[, - controls]))
  b <- rowMeans(as.matrix(z[, controls]))
  bg <- (1:ng)[order(b - a, decreasing = T)]

  bg
}

# =============================================================================.
#' Background Read Density
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{CPSES},
#'   \link{QuickShiftClustering},
#'   \link{cdadadr}
# -----------------------------------------------------------------------------.
#' @inheritParams cdadadr
#'
#' @param controls
#' columns corresponding to control measurements (Input, IgG, etc.).
#'
#' @param cst
#' core subset density threshold.
#'
#' @param ncl
#' number of clusters used to partition the core subset.
#'
#' @param mincs
#' minimum size of each cluster, as number of observations.
#'
#' @return
#' \code{BRD} returns a \code{list} with the following elements:
#' \item{parameters}{call parameters of the function}
#' \item{status}{execution status}
#' \item{nonzero}{indices of initial observations with count > 0}
#' \item{dred}{result of \link{cdadadr} applied to non-zero observations}
#' \item{coreset}{indices of non-zero observations selected as core subset}
#' \item{clusters}{result of \link{QuickShiftClustering} applied to the core subset}
#' \item{theta}{distribution parameters for each cluster in the core subset}
#' \item{log2counts}{dithered and log2 transformed counts}
# -----------------------------------------------------------------------------.
#' @export
BRD <- function(
  cnt, controls = NULL,
  dither = 3, smobs = T, cvt = 0.5, npc = NA, zscore = T,
  knn = 200, cst = 0.5, ncl = 3, mincs = 300, progress = F
) {

  if(is.null(colnames(cnt))) stop("missing column names in the count matrix")

  xps <- colnames(cnt)
  inp <- match(controls, xps)

  parameters <- list(
    experiments = xps, controls = xps[inp],
    dither = dither, smobs = smobs, cvt = cvt, npc = npc, zscore = zscore,
    knn = knn, cst = cst, ncl = ncl, mincs = mincs
  )

  # Cleanup of missing values
  cleanup <- FiniteValues(log(cnt))
  ldc <- log2(DitherCounts(cnt))
  ldc <- ldc[cleanup, ] # Dithered and log2 transformed counts
  cnt <- cnt[cleanup, ] # Raw counts

  # Count density after dithering and dimensionality reduction
  dred <- cdadadr(
    cnt, dither = dither, smobs = smobs, cvt = cvt, npc = npc,
    zscore = zscore, knn = knn
  )

  # Select the core subset using density threshold
  core <- which(rankstat(dred$density) > cst)
  core <- with(dred, list(i = core, d = density[core], x = projection[core, ]))

  # Clustering by density gradient ascent in the dred space
  qsc <- with(core, QuickShiftClustering(x, d, n = ncl, progress = progress))

  # Background cluster
  if(! is.null(controls)) {
    if(ncl > 1) {
      bg_rnk <- bg_ranking(ldc[core$i, ], by = qsc$membership, controls = inp)
    } else {
      bg_rnk <- 1
    }
    bg_clu <- bg_rnk[1]
    bg_idx <- core$i[qsc$membership == bg_clu]
  }

  # Fit multivariate gaussian model to each cluster of minimum size
  theta <- list()
  for(grp in 1:ncl) {
    idx <- core$i[qsc$membership == grp] # select cluster
    if(length(idx) > mincs) {
      mm <- mv_gmm(ldc[idx, ], ns = 1)   # EM fitting
      theta[[grp]] <- mm$theta[[1]]
    } else {
      theta[[grp]] <- NA
    }
  }

  res <- list(
    parameters = parameters,
    status     = "clustering",
    nonzero    = which(cleanup), # observation with count > 0
    dred       = dred,           # cdadadr results
    coreset    = core,           # core subset of observations
    clusters   = qsc,            # QuickShiftClustering results
    theta      = theta,          # distribution parameters for each cluster
    log2counts = ldc             # dithered and log2 transformed counts
  )

  # Background cluster
  if(! is.null(controls)) {
    res$status <- "BRD candidate cluster selected"
    bg_theta <- theta[[bg_clu]]
    if(! is.na(bg_theta)) {
      res <- c(
        res, list(
          normfactors = mean(bg_theta$mu) - bg_theta$mu,
          bg_cluster  = bg_clu,
          bg_theta    = bg_theta,
          bg_members  = bg_idx
        )
      )
    } else {
      res$status <- "BRD candidate cluster has few observations"
      warning(res$status)
    }
  }

  res
}
