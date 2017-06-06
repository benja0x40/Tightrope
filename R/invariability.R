# FUNCTIONS | INVARIABILITY ####################################################

# =============================================================================.
# ranking
# -----------------------------------------------------------------------------.
ranked_p <- function(x) { (rank(x) - 0.5) / length(x) }

# =============================================================================.
#' partition of multivariate observations
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate observations
#'
#' @param p
#' statistical score associated to observations (e.g. estimated density)
#'
#' @param ns
#' number of clusters
#'
#' @param columns
#' components to be considered
#'
#' @param max_iter
#' maximum number of iterations for expectation-maximization
#'
#' @return
#' populations returns a \code{list} with the following elements:
#' \item{grp}{list of clusters resulting from function \link{quickShiftCutClusters}}
#' \item{qstree}{quickshift tree}
#' \item{theta}{list of multivariate distribution parameters}
# -----------------------------------------------------------------------------.
#' @export
populations <- function(x, p, ns, columns = NULL, max_iter = 100) {

  if(is.null(columns)) columns <- 1:ncol(x)

  g <- quickShift(x[, columns], d = p)
  grp <- quickShiftCutClusters(g, n = ns)
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
#' background population
# -----------------------------------------------------------------------------.
#' @param x
#'
#' @param grp
#'
#' @param bg_cols
#'
#' @return cluster id
# -----------------------------------------------------------------------------.
#' @export
bg_grp <- function(x, grp, bg_cols) {
  bg <- x[grp$center, ]
  bg <- cbind(rowMeans(bg[, - bg_cols]), rowMeans(bg[, bg_cols]))
  print(bg)
  bg <- which(bg[, 2] > bg[, 1])
  bg
}

# =============================================================================.
#' lowVariability
# -----------------------------------------------------------------------------.
#' @param cnt
#' matrix of read counts (columns = samples, rows = intervals)
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
#' @export
lowVariability <- function(
  cnt, bg_cols, fit_idx, viz_idx,
  n = 3, k = 64, dithering = 3, max_iter = 100
) {

  chk <- finiteValues(log(cnt))
  i <- 0
  d <- 0
  s <- 0
  y <- 0
  pb <- txtProgressBar(1, dithering, char = "|", style = 3)
  while(i < dithering) {
    x <- cnt
    x <- ditherCounts(x)
    x <- log2(x)
    x <- x[chk, ]
    y <- y + x / dithering
    # center and reduce count variables
    x <- t((t(x) - colMeans(x)) / apply(x, MARGIN = 2, sd))
    # estimate local density
    d <- d + knn_stat(x, k = k, FUN = knn_density) / dithering
    # estimate local deviation
    s <- s + knn_stat(x, k = k, FUN = knn_sd) / dithering

    i <- i + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)

  x <- y

  # combine density and deviation
  p.d <- ranked_p(d)
  p.s <- ranked_p(1 / s)
  v <- sqrt(p.d * p.s)

  layout(matrix(1:6, 3, 2, byrow = T))

  grp_d <- populations(x, p = d, ns = n, columns = fit_idx, max_iter = max_iter)
  grp_s <- populations(x, p = 1 / s, ns = n, columns = fit_idx, max_iter = max_iter)
  grp_v <- populations(x, p = v, ns = n, columns = fit_idx, max_iter = max_iter)

  plot(x[, viz_idx], pch = 19, cex = 0.3, col = colorize(d))
  clr <- rainbow(grp_d$nbr, alpha = 0.2)
  clr <- plot_groups_2D(
    x[, viz_idx], clr = clr[grp_d$membership], pch = 19, cex = 0.3
  )
  plot_sources_2D(
    grp_d$theta, components = viz_idx, clr = grey(0), lwd = 2, center = F
  )
  abline(a = diff(grp_d$theta[[bg_grp(x, grp_d, bg_cols)]]$mu[viz_idx]), b = 1)

  plot(x[, viz_idx], pch = 19, cex = 0.3, col = colorize(1 / s))
  clr <- rainbow(grp_s$nbr, alpha = 0.2)
  clr <- plot_groups_2D(
    x[, viz_idx], clr = clr[grp_s$membership], pch = 19, cex = 0.3
  )
  plot_sources_2D(
    grp_s$theta, components = viz_idx, clr = grey(0), lwd = 2, center = F
  )
  abline(a = diff(grp_s$theta[[bg_grp(x, grp_s, bg_cols)]]$mu[viz_idx]), b = 1)

  plot(x[, viz_idx], pch = 19, cex = 0.3, col = colorize(v))
  clr <- rainbow(grp_v$nbr, alpha = 0.2)
  clr <- plot_groups_2D(
    x[, viz_idx], clr = clr[grp_v$membership], pch = 19, cex = 0.3
  )
  plot_sources_2D(
    grp_v$theta, components = viz_idx, clr = grey(0), lwd = 2, center = F
  )
  abline(a = diff(grp_v$theta[[bg_grp(x, grp_v, bg_cols)]]$mu[viz_idx]), b = 1)
}
