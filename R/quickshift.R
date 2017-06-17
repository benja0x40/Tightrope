# FUNCTIONS | QUICKSHIFT #######################################################

# =============================================================================.
#' QuickShift algorithm (Vedaldi & Soatto, 2008)
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param d
#' numeric vector representing the density estimation at each observation.
#'
#' @param progress
#' progress bar (logical, default = F)
#'
#' @return QuickShift returns a graph object (igraph package)
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
QuickShift <- function (x, d, progress = F) {

  if(progress) pb <- txtProgressBar(min = 1, max = nrow(x), char = "|", style= 3)

  g <- graph.empty(n = nrow(x))
  i.a <- which(FiniteValues(x) & ! is.na(d))

  while(length(i.a) > 1) {

    knn <- get.knnx(
      data = x[i.a,], query = x[i.a,], k = 2, algorithm = "kd_tree"
    )
    i.b <- i.a[knn$nn.index[, 2]]

    chk <- d[i.b] >= d[i.a]
    g <- g + igraph::edges(
      rbind(i.a, i.b)[, chk], distance = knn$nn.dist[chk, 2]
    )
    i.a <- i.a[! chk]

    if(progress) setTxtProgressBar(pb, nrow(x) - length(i.a) + 1)
  }
  if(progress) close(pb)

  g
}

# =============================================================================.
#' QuickShiftCutClusters
# -----------------------------------------------------------------------------.
# TODO: fix bug with ordering cluster sizes and memberships
# -----------------------------------------------------------------------------.
#' @param g
#' QuickShift graph resulting from the \link{QuickShift} function
#'
#' @param n
#' expected number of clusters
#'
#' @param ecut
#' maximum branch length
#'
#' @return
#' QuickShiftCutClusters returnd a \code{list} with the following elements
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
QuickShiftCutClusters <- function(g, n = NULL, ecut = NULL) {
  if(is.null(ecut) & ! is.null(n)) {
    ecut <- mean(sort(E(g)$distance, decreasing = T)[c(n, n+1) - 1])
  }
  if(is.null(ecut)) stop("Missing value for n or ecut")
  V(g)$id <- 1:length(V(g))
  g <- g - E(g)[distance > ecut]
  r <- which(igraph::degree(g, mode = "out") == 0)
  grp <- list(
    membership = rep(NA, length(V(g))),
    csize      = rep(0, length(r)),
    center     = r,
    nbr        = length(r),
    ecut       = ecut
  )
  for(k in 1:length(r)) {
    sg <- subcomponent(g, v = r[k], mode = "in")
    grp$membership[sg$id] <- k
  }
  grp$csize <- as.vector(table(grp$membership))
  o <- order(grp$csize, decreasing = T, na.last = T)
  grp$csize <- grp$csize[o]
  grp$membership <- o[grp$membership]
  grp
}

# =============================================================================.
#' QuickShiftClustering
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param d
#' numeric vector representing the density estimation at each observation.
#'
#' @param n
#' number of clusters
#'
#' @param ...
#'
#' @return
#' QuickShiftClustering returns a \code{list} with the following elements
# -----------------------------------------------------------------------------.
#' @export
QuickShiftClustering <- function (x, d, n, ...) {
  qs <- QuickShift(x, d, ...)
  qs <- QuickShiftCutClusters(qs, n = n)
  qs
}
