# FUNCTIONS | QUICKSHIFT #######################################################

# =============================================================================.
#' QuickShift algorithm (Vedaldi & Soatto, 2008)
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param d
#' statistical score associated to observations (e.g. estimated density)
#'
#' @param plot
#' logical (default = F)
#'
#' @return QuickShift returns a graph object (igraph package)
# -----------------------------------------------------------------------------.
#' @export
QuickShift <- function (x, d, plot = F) {
  g <- graph.empty(n = nrow(x))
  V(g)$d <- d
  i.a <- which(FiniteValues(x) & ! is.na(d))
  if(plot) {
    plot(x[, 1], x[, 2], pch = 20, col=rgb(0, 0, 0, 0.1))
  }
  pb <- txtProgressBar(min = 1, max = nrow(x), char = "|", style= 3)
  while(length(i.a) > 1) {
    knn <- get.knnx(
      data = x[i.a, ], query = x[i.a, ], k = 2, algorithm = "kd_tree"
    )
    i.b <- i.a[knn$nn.index[, 2]]
    chk <- d[i.b] >= d[i.a]
    g <- g + igraph::edges(rbind(i.a, i.b)[, chk], distance = knn$nn.dist[chk, 2])
    if(plot) {
      arrows(
        x[i.a[chk], 1], x[i.a[chk], 2], x[i.b[chk], 1], x[i.b[chk], 2],
        length = 0.1, col = rgb(0, 0, 0, 0.5)
      )
    }
    i.a <- i.a[! chk]
    setTxtProgressBar(pb, nrow(x) - length(i.a) + 1)
  }
  close(pb)
  V(g)$weight <- ego_size(g, order = nrow(x), mode = "in")
  E(g)$weight <- V(g)[ends(g, E(g))[, 1]]$weight
  g
}

# =============================================================================.
#' QuickShift root and top level branches
# -----------------------------------------------------------------------------.
#' @param i
#' node index
#'
#' @param d
#' depth level
#'
#' @param dmax
#' maximum depth
#'
#' @return QuickShiftSearchClusters returns a \code{data.frame}
# -----------------------------------------------------------------------------.
#' @export
QuickShiftSearchClusters <- function(i = 0, d = 0, dmax) {
  res <- NULL
  if(d < dmax) {
    if(i==0) {
      i <- which(igraph::degree(g, mode = "out") == 0)
    }
    e <- incident_edges(g, i, mode = "in")[[1]]
    if(length(e)>0) {
      j <- ends(g, e)[,1]
      # arrows(
      #   x[j], y[j], x[i], y[i], col=rgb(1, 0, 0, 0.75), length = 0.1, lwd = 1
      # )
      res <- cbind(
        edge = as.numeric(e), distance = e$distance, weight = V(g)[j]$weight
      )
      for(k in j) {
        res <- rbind(
          res, QuickShiftSearchClusters(k, d + 1, dmax)
        )
      }
    }
  }
  res
}

# =============================================================================.
#' QuickShift graph partition
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
#' QuickShiftCutClusters return a \code{list} with the following elements
# -----------------------------------------------------------------------------.
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
