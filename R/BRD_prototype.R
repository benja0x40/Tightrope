# =============================================================================.
#' Background Read Density
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{CPSES},
#'   \link{QuickShiftClustering},
#'   \link{CDaDaDR.prototype}
# -----------------------------------------------------------------------------.
#' @inheritParams CDaDaDR.prototype
#'
#' @param controls
#' columns corresponding to control measurements (Input, IgG, etc.).
#'
#' @param ignore
#' logical vector indicating observations (rows) to be ignored
#' (default = F, none).
#'
#' @param bdt
#' numeric vector of length 2 defining background density thresholds both
#' expressed as proportions between 0 and 1.
#' (\code{bdt[1]}) specifies the global threshold used to discard observations
#' with low density prior to clustering.
#' (\code{bdt[2]}) determines the maximum density loss allowed
#' when selecting core observations relatively to the local maximum density
#' in each cluster, and thus defining the candidate background populations.
#' By default the value of \code{bdt} is \code{c(0.3, 0.1)} meaning that,
#' in terms of density percentiles, the bottom 30% will be filtered out before
#' clustering and only the top 10% can be selected as background candidates
#' among each cluster.
#'
#' @param ncl
#' number of clusters for partitioning observations.
#'
#' @param mincs
#' minimum size of each cluster, as number of observations.
#'
#' @return
#' \code{BRD} returns a \code{list} with the following elements:
#' \item{parameters}{call parameters of the function.}
#' \item{status}{execution status.}
#' \item{nonzero}{indices of initial observations with count > 0.}
#' \item{dred}{
#'   dimensionality-reduced non-zero observations
#'   (result from the \link{CDaDaDR.prototype} function).
#' }
#' \item{subsets}{
#'   partition of non-zero observations into background candidate subsets.
#' }
#' \item{populations}{
#'   summary of core populations.
#' }
#' \item{theta}{
#'   fitted distribution parameters for each core candidate populations.
#' }
#' \item{log2counts}{
#'   dithered and log2 transformed counts.
#' }
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
BRD.prototype <- function(
  cnt, controls = NULL, widths = 1, ignore = F,
  dither = 3, smobs = F, cvt = 0, npc = 2, zscore = T, knn = 200, rare = 0.01,
  bdt = c(0.4, 0.1), ncl = 2, mincs = 100, progress = F
) {

  if(is.null(colnames(cnt))) stop("missing column names in the count matrix")

  if(length(ignore) == 1) ignore <- rep(ignore, nrow(cnt))
  if(length(widths) == 1) widths <- rep(widths, nrow(cnt))

  if(length(ignore) != nrow(cnt)) {
    stop("inconsistent count matrix and ignore vector")
  }
  if(length(widths) != nrow(cnt)) {
    stop("inconsistent count matrix and widths vector")
  }

  xps <- colnames(cnt)
  inp <- match(controls, xps)

  parameters <- list(
    experiments = xps, controls = xps[inp],
    dither = dither, smobs = smobs, cvt = cvt, npc = npc, zscore = zscore,
    knn = knn, bdt = bdt, ncl = ncl, mincs = mincs
  )

  # Cleanup of missing values and ignored if provided
  cleanup <- FiniteValues(log2(cnt)) & ! ignore
  ldc <- log2(DitherCounts(cnt))
  ldc <- ldc[cleanup, ]     # Dithered and log2 transformed counts
  cnt <- cnt[cleanup, ]     # Raw counts
  widths <- widths[cleanup] # Region sizes

  # Precompute the mean of each observation for controls and experiments
  movs <- rowMeans(as.matrix(log2(cnt[,   inp])))

  # Count density after dithering and dimensionality reduction
  dred <- CDaDaDR.prototype(
    cnt, knn = knn,
    widths = widths, smobs = smobs, movs = movs,
    dither = dither, cvt = cvt, npc = npc, zscore = zscore, rare = rare
  )

  # Filter out low density observations
  rd <- rankstat(dred$density)
  sbs <- which(rd > bdt[1])
  sbs <- with(
    dred, list(i = sbs, d = density[sbs], r = rd[sbs], x = projection[sbs, ])
  )

  # Find clusters by density gradient ascent in the dred space
  qsc <- with(sbs, QuickShiftClustering(x, d, n = ncl, progress = progress))
  sbs$cluster <- qsc$membership

  # Select top density (core) observations in each cluster and count members
  pop <- data.frame(cluster = rep(0, ncl), core = rep(0, ncl))
  pop$max_density <- as.vector(by(sbs$r, sbs$cluster, max))
  sbs$core  <- F
  for(grp in 1:ncl) {
    k <- with(sbs, cluster == grp)
    sbs$core[k] <- sbs$r[k] > max(sbs$r[k]) - bdt[2]
    pop$cluster[grp] <- sum(k)
    pop$core[grp]    <- sum(sbs$core[k])
  }

  # Rank cluster cores by background levels (enrichment in control samples)
  if(! is.null(controls)) {
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
  }

  # Fit a multivariate gaussian model to each cluster core of minimum size
  theta <- list()
  for(grp in 1:ncl) {
    idx <- with(sbs, i[core & cluster == grp]) # select cluster core
    if(length(idx) > mincs) {
      mm <- mv_gmm(ldc[idx, ], ns = 1)          # EM fitting
      theta[[grp]] <- mm$theta[[1]]
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
  if(! is.null(controls)) {
    res$status <- "BRD candidate cluster selected"
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
      res$status <- "Background candidate population is below the expected minimum size"
      warning(res$status)
    }
  }

  res
}
