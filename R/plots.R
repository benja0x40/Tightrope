# =============================================================================.
#' Plot empirical distributions of read counts
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{DitherCounts},
#'   \link{SideBySideDensity},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @param x
#' a read count matrix, preferably dithered and log2 transformed.
#'
#' @param ...
#' optional arguments forwarded to the \link{SideBySideDensity} function.
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @examples
#' # Simulation of a read count matrix with 3 populations of 10000 observations
#' p <- c(10000, 10000, 10000)
#' m <- DefineSimulation(
#'   chip = 5, patterns = c("^", "v"), enrichment = c(1.0, 3), replicate = 2
#' )
#' r <-  MakeSimulation(p = p, m = m)
#'
#' grp <- r$group # Population memberships
#' cnt <- r$data  # Simulated counts
#'
#' # Prepare figure layout and graphic options
#' layout(matrix(1:4, 2, 2, byrow = T))
#' par(pch = 20)
#'
#' # Show the empirical distribution of simulated populations
#' l2c <- log2(DitherCounts(cnt)) # Dithering and log2 transformation
#' xyl <- range(l2c[FiniteValues(l2c), ])
#'
#' r <- PlotCountDistributions(l2c, ylim = xyl, main = "Total")
#' for(i in sort(unique(grp))) {
#'   main <- ifelse(i == 1, "Invariable subset", paste("Variable subset", i - 1))
#'   r <- PlotCountDistributions(l2c[grp == i, ], ylim = xyl, main = main)
#' }
# -----------------------------------------------------------------------------.
#' @export
PlotCountDistributions <- function(x, ...) {
  r <- SideBySideDensity(
    x, method = "ash", parameters = list(color = "Wry"),
    ylab = "log2(counts)", las = 2, ...
  )
}

# HIDDEN #######################################################################

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
#' @keywords internal
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
#' @keywords internal
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
PlotQuickShift <- function(
  x, g, new = T, length = 0.05, col = rgb(0, 0, 0, 0.2), ...
) {

  if(new) plot(x[, 1], x[, 2], type='n')

  el <- as_edgelist(g)
  suppressWarnings(
    arrows(
      x[el[, 1], 1], x[el[, 1], 2], x[el[, 2], 1], x[el[, 2], 2],
      length = length, col = col, ...
    )
  )
}
