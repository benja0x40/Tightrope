# FUNCTIONS | 2D PLOTS #########################################################

# =============================================================================.
#' empty.plot
# -----------------------------------------------------------------------------.
#' @param axes logical
#' @param xlab character
#' @param ylab character
#' @param ...
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
empty.plot <- function(axes = T, xlab = '', ylab = '', ...) {
  plot(0, type = 'n', axes = axes, xlab = xlab, ylab = ylab, ...)
}

# =============================================================================.
#' PlotQuickShift
# -----------------------------------------------------------------------------.
#' @param x matrix
#' @param g graph
#' @param ...
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
PlotQuickShift <- function(x, g, new = T, ...) {

  if(new) plot(x[, 1], x[, 2], type='n')

  el <- as_edgelist(g)
  suppressWarnings(
    arrows(
      x[el[, 1], 1], x[el[, 1], 2], x[el[, 2], 1], x[el[, 2], 2],
      length = 0.05, col = rgb(0, 0, 0, 0.2)
    )
  )
}

# =============================================================================.
#' plotHistograms
# -----------------------------------------------------------------------------.
#' @param x numeric vector
#' @param y numeric vector
#' @param bins integer
#' @param xlim range
#' @param log logical
#' @param rel logical
#' @param ...
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
plotHistograms <- function(x, y, bins = 100, xlim = NULL, log = F, rel = F, ...) {

  n <- length(x)

  if(log) {
    x <- log2(x)
    y <- log2(y)
  }

  bins <- bins - 1
  r <- diff(range(x, y, na.rm = T))
  brk <- 0:bins/bins * (1.2 * r) - 0.1 * r + min(x, y, na.rm = T)
  h.x <- hist(x, breaks = brk, plot = F)
  h.y <- hist(y, breaks = brk, plot = F)

  h.x$counts <- h.x$counts / n
  h.y$counts <- h.y$counts / n

  if(rel) {
    chk <- h.x$counts > 0
    h.x$counts[chk] <- with(h.x, counts[chk] / sum(diff(breaks[c(chk, T)]) * counts[chk]))
    chk <- h.y$counts > 0
    h.y$counts[chk] <- with(h.y, counts[chk] / sum(diff(breaks[c(chk, T)]) * counts[chk]))
  }

  if(is.null(xlim)) xlim = range(x, y)
  ylim <- range(h.x$counts, h.y$counts)
  ylim <- min(ylim) + c(0, 1.1 * diff(ylim))
  empty.plot(xlim = xlim, ylim = ylim, yaxs = 'i', ...)
  chk <- h.x$counts > 0
  k <- length(h.y$mids)
  x <- c(h.y$mids[1], h.y$mids, h.y$mids[k])
  y <- c(0, h.y$counts, 0)
  polygon(x, y, border = NA, col = rgb(1, 0.6, 0, 0.5))
  chk <- h.y$counts > 0
  points(h.y$mids[chk], h.y$counts[chk], type = 'l', col = rgb(1, 0, 0, 0.75), lwd = 2)
  points(h.x$mids[chk], h.x$counts[chk], type = 'h', lwd = 2, col = grey(0.3))
}

# =============================================================================.
#' plot_samples
# -----------------------------------------------------------------------------.
#' @param x matrix
#' @param idx vector
#' @param symetric logical
#' @param cex numeric
#' @param ...
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
plot_samples <- function(
  x, idx = 1:2, symetric = T, cex = 0.5, xlim = NULL, ylim = NULL, ...
) {

  rng = matrix(range(x), 2, 2)
  x <- x[, idx]
  xlab <- colnames(x)[1]
  ylab <- colnames(x)[2]
  if(! symetric) rng <- apply(x, 2, range)

  if(! is.null(xlim)) rng[, 1] <- xlim
  if(! is.null(ylim)) rng[, 2] <- ylim

  ScatterPlot(
    x, xlim = rng[, 1], ylim = rng[, 2], cex = cex,
    xlab = xlab, ylab = ylab, ...
  )

  # for(g in grp_order) points(x[grp == g, ], col = grp_clr[g], cex = cex)
}

# =============================================================================.
#' 2D scatter plot with groups
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix, the two first columns being used as x and y coordinates
#'
#' @param p
#' matrix of group membership probabilities
#'
#' @param posterior
#' membership probability threshold
#'
#' @param contrast
#' membership probability contrast factor
#'
#' @param clr
#' vector of colors
#'
#' @param alpha
#' color transparency
#'
#' @param ...
#'
#' @return plot_groups_2D returns a vector of colors
# -----------------------------------------------------------------------------.
#' @export
plot_groups_2D <- function(x, p = NULL, posterior = 0, contrast = 1, clr = NULL, alpha = 0.1, ...) {

  nx <- nrow(x) # number of observations

  if(is.null(p)) p <- matrix(1, nx)
  else p <- as.matrix(p)
  ns <- ncol(p) # number of sources

  if(is.null(clr)) clr <- rainbow(ns, alpha = alpha)
  if(length(clr) == 1) clr = rep(clr, nx)

  if(ns > 1) {
    low <- t(apply(p, MARGIN = 1, FUN = sort, decreasing = T))
    low <- (low[, 1] < posterior) | (low[, 1] / low[, 2] < contrast)
    grp <- apply(p, MARGIN = 1, FUN = which.max)
    clr <- clr[grp]
    clr[low] <- grey(0.5, 0.2)
  }

  plot(x, col = clr, ...)

  clr
}

# =============================================================================.
#' plot 2D normal distributions
# -----------------------------------------------------------------------------.
#' @param theta
#' list of multivariate distribution parameters
#'
#' @param components
#' components to be represented (default = 1:2)
#'
#' @param p
#' probability boundary of represented distributions
#'
#' @param clr
#' vector of colors
#'
#' @param center
#' logical (default = F)
#'
#' @param ...
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @export
plot_sources_2D <- function(theta, components = NULL, p = 0.5, clr = NULL, center = F, ...) {

  ns <- length(theta) # number of sources
  if(is.null(clr)) clr <- grey(0, alpha = 0.5)
  if(length(clr) == 1) clr = rep(clr, ns)
  if(is.null(components)) components <- 1:2

  for(i in 1:ns) {
    mu <- theta[[i]][[2]][components]
    ellipse(
      mu = mu, sigma = theta[[i]][[3]][components, components],
      alpha = p, col = clr[i], ...
    )
    if(center) points(mu[1], mu[2], col = clr[i], pch = 19)
  }
}
