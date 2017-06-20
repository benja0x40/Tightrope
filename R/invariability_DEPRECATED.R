# FUNCTIONS | INVARIABILITY ####################################################

# =============================================================================.
#' partition of multivariate observations
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param tx
#' numeric matrix representing transformed observations (default = x).
#'
#' @param p
#' statistical score associated to observations (e.g. estimated density).
#'
#' @param ns
#' number of clusters.
#'
#' @param max_iter
#' maximum number of iterations for expectation-maximization
#'
#' @return
#' populations returns a \code{list} with the following elements:
#' \item{grp}{list of clusters resulting from function \link{QuickShiftCutClusters}}
#' \item{qstree}{QuickShift tree}
#' \item{theta}{list of multivariate distribution parameters}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
populations <- function(x, tx = NULL, p, ns, max_iter = 100) {

  if(is.null(tx)) g <- QuickShift(x, d = p)
  else g <- QuickShift(tx, d = p)

  grp <- QuickShiftCutClusters(g, n = ns)

  bg <- list()
  for(i in 1:ns) {
    y <- x[grp$membership == i, ]
    theta <- mv_init_param(y, ns = 1)
    EM <- mv_emloop(y, theta, max_iter = max_iter)
    bg[[i]] <- EM$theta[[1]]
  }

  c(grp, list(qstree = g, theta = bg))
}

# =============================================================================.
#' lowVariability
# -----------------------------------------------------------------------------.
#' @param cnt
#' matrix of read counts (rows = observations, columns = samples or conditions).
#'
#' @param bg_cols
#' columns with background/control profiles (Input, IgG, etc.)
#'
#' @param fit_idx
#' columns considered for low variability analysis
#'
#' @param viz_idx
#' columns for 2D plots
#'
#' @param mincv
#' minimal component variance
#'
#' @param n
#' number of clusters
#'
#' @param k
#' number of nearest neighbors
#'
#' @param dithering
#' number of dithering passes
#'
#' @param max_iter
#' maximum number of iterations
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
lowVariability <- function(
  cnt, bg_cols, fit_idx, viz_idx, mincv = 0.05,
  n = 3, k = 64, dithering = 3, max_iter = 100
) {

  chk <- FiniteValues(log(cnt))
  i <- 0
  d <- 0
  s <- 0
  y <- 0
  pb <- txtProgressBar(0, dithering, char = "|", style = 3)
  while(i < dithering) {

    x <- cnt
    x <- DitherCounts(x)
    x <- log2(x)
    x <- x[chk, ]

    y <- x[, fit_idx]

    # optional dimensionality reduction (pca projection)
    pca_msg <- ""
    if(mincv > 0) {
      y <- prcomp(x, retx = T, center = T, scale. = T)
      idx <- which(y$sdev / sum(y$sdev) > mincv)
      if(length(idx) < 2) stop("mincv threshold is too stringent")
      pca_msg <- paste("components =", length(idx))
      y <- y$x[, idx]
    }
    # center and reduce count variables
    # y <- t((t(y) - colMeans(y)) / apply(y, MARGIN = 2, sd))
    # estimate local density
    d <- d + knn_stat(y, k = k, FUN = knn_density) / dithering
    # estimate local deviation
    s <- s + knn_stat(y, k = k, FUN = knn_sd) / dithering

    i <- i + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)

  # combine density and deviation
  s <- 1 / s
  p.d <- rankstat(d)
  p.s <- rankstat(s)
  v <- sqrt(p.d * p.s)

  grp_d <- populations(x, tx = y, p = d, ns = n, max_iter = max_iter)
  grp_s <- populations(x, tx = y, p = s, ns = n, max_iter = max_iter)
  grp_v <- populations(x, tx = y, p = v, ns = n, max_iter = max_iter)

  grp_d$bg_grp <- bg_grp(x, grp_d, bg_cols)
  grp_s$bg_grp <- bg_grp(x, grp_s, bg_cols)
  grp_v$bg_grp <- bg_grp(x, grp_v, bg_cols)


  control_graphs <- function(x, viz_idx, y, d, grp) {
    clr <- rainbow(grp$nbr, alpha = 0.2)
    plot(y[, 1:2], pch = 19, cex = 0.3, col = colorize(d), main = pca_msg)
    plot(x[, viz_idx], pch = 19, cex = 0.3, col = colorize(d))
    clr <- plot_groups_2D(
      x[, viz_idx], clr = clr[grp$membership], pch = 19, cex = 0.3
    )
    plot_sources_2D(
      grp$theta, components = viz_idx, clr = grey(0), lwd = 2, center = F
    )
    a <- diff(grp$theta[[grp$bg_grp]]$mu[viz_idx])
    abline(a = a, b = 1)
    for(i in 1:grp$nbr) {
      text(matrix(grp$theta[[i]]$mu[viz_idx], nrow = 1), labels = i)
    }
  }

  layout(matrix(1:9, 3, 3, byrow = T))
  control_graphs(x, viz_idx, y, d = d, grp_d)
  control_graphs(x, viz_idx, y, d = s, grp_s)
  control_graphs(x, viz_idx, y, d = v, grp_v)

  list(grp_d = grp_d, grp_s = grp_s, grp_v = grp_v)
}
