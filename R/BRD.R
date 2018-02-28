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
  a <- rowMeans(as.matrix(z[, - controls]))
  b <- rowMeans(as.matrix(z[, controls]))
  bg <- (1:ng)[order(b - a, decreasing = T)]

  bg
}

# =============================================================================.
#' DensityCorrectedByIntensity
# -----------------------------------------------------------------------------.
#' @param d
#' density vector
#'
#' @param i
#' intensity vector
#'
#' @return
#' \code{codeDensityCorrectedByIntensity} returns the vector of corrected
#' densities.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
DensityCorrectedByIntensity <- function(d, i, k) {
  ab <- function(rx, ry) {
    f  <- diff(ry) / diff(rx)
    f <- c(a = ry[1] - f * rx[1], b = f)
    f
  }
  yabx <- function(x, f) {
    f[1] + x * f[2]
  }
  n <- length(d)
  o <- order(i)
  k <- round(0.1 * n)
  s <- caTools::runmean(d[o], k = k)
  ri <- range(i)
  rd <- s[c(1, n)]
  if(rd[1] > rd[2]) {
    rd <- c(rd[2] / rd[1], 1)
    d <- d * yabx(i, ab(ri, rd))
  } else {
    rd <- c(1, rd[1] / rd[2])
  }
  d
}

# =============================================================================.
#' Background Read Density
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{PlotBRD},
#'   \link{CDaDaDR},
#'   \link{QuickShiftClustering},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @inheritParams CDaDaDR
#'
#' @param controls
#' count columns corresponding to control measurements (Input, IgG, etc.).
#'
#' @param bdt
#' numeric vector of length 2 defining background density thresholds both
#' expressed as proportions between 0 and 1.
#' \code{bdt[1]} specifies the global threshold used to discard observations
#' with low density prior to clustering.
#' \code{bdt[2]} determines the maximum density loss allowed
#' when selecting core observations relatively to the local maximum density
#' in each cluster, and thus defining the candidate background populations.
#' By default the value of \code{bdt} is \code{c(0.2, 0.05)} meaning that,
#' in terms of density percentiles, the bottom 20 percents will be filtered out
#' before clustering and only the top 5 percents can be selected as background
#' candidates among each cluster.
#'
#' @param ncl
#' number of clusters partitioning the population of observations with density
#' above \code{bdt[1]} percents.
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
#'   dimensionality-reduced non-zero observations
#'   (result from the \link{CDaDaDR} function).
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
# -----------------------------------------------------------------------------.
#' @export
BRD <- function(
  cnt, controls, smobs = T, dither = 5, zscore = T, knn = 300,
  bdt = c(0.2, 0.05), ncl = 1, mincs = 100, progress = F
) {

  if(is.null(colnames(cnt))) stop("missing column names in the count matrix")
  if(length(controls) < 2) stop("BRD requires at least 2 controls")

  xps <- colnames(cnt)
  inp <- match(controls, xps)

  parameters <- list(
    experiments = xps, controls = xps[inp],
    smobs = smobs, dither = dither, zscore = zscore, knn = knn,
    bdt = bdt, ncl = ncl, mincs = mincs
  )

  # Cleanup of missing values
  cleanup <- FiniteValues(log2(cnt))
  ldc <- log2(DitherCounts(cnt))
  ldc <- ldc[cleanup, ]     # Dithered and log2 transformed counts
  cnt <- cnt[cleanup, ]     # Raw counts

  # Precompute the mean of each observation in control and enrichment conditions
  movs <- rowMeans(as.matrix(log2(cnt[,   inp]))) # , drop = F
  menr <- rowMeans(as.matrix(log2(cnt[, - inp]))) # , drop = F
  intensity <- xy2md(menr, movs)$m

  # Compute count density after dithering and dimensionality reduction
  dred <- CDaDaDR(
    cnt, knn = knn, smobs = smobs, movs = movs,
    dither = dither, npc = 2, zscore = zscore, method = "pca"
  )

  # Correct density
  ds <- DensityCorrectedByIntensity(dred$density, intensity, k = knn)
  dred <- c(dred, list(ctl = movs, enr = menr, intensity = intensity))
  dred$knn_density <- dred$density
  dred$density <- ds

  # Filter out observations with lower density than specified
  rd <- rankstat(dred$density)
  sbs <- which(rd > bdt[1])
  sbs <- with(
    dred, list(
      i = sbs,
      b = intensity[sbs], d = density[sbs], r = rd[sbs], x = projection[sbs, ]
    )
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
    if(length(idx) > mincs) {
      mm <- mv_gmm(ldc[idx, ], ns = 1)         # EM fitting
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

  res
}
