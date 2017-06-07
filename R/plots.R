# FUNCTIONS | 2D PLOTS #########################################################

# =============================================================================.
#' quantile color mapping
# -----------------------------------------------------------------------------.
#' @param x
#' numeric vector
#'
#' @return colorize returns a vector of colors
# -----------------------------------------------------------------------------.
#' @export
colorize <- function(x) {

  q <- c(0, 0.33, 0.66, 1.0)
  k <- c("grey", "black", "red", "yellow")
  x <- rankstat(x)

  clr_prm <- defineColors(thresholds = q, colors = k)
  clr <- makeColors(x, parameters = clr_prm)
  clr
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
