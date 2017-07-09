# =============================================================================.
# NOT WORKING
# =============================================================================.
#' Local curvature (2D)
# -----------------------------------------------------------------------------.
#' @param p density vector
#'
#' @inheritParams knn_density
#'
#' @return numeric vector
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
knn_curvature <- function(p, k, i, d, smoothing = F, adaptative = T) {

  # Second Order Derivative Gaussian Kernel (1 dimension)
  # sodgk <- function(x, s = max(x)) {
  #   s <- s / 3
  #   x2 <- x * x
  #   s2 <- s * s
  #   g <- - ((s2 - x2) / s^5 * sqrt(2 * pi)) * exp(- x2 / (2 * s2))
  #   g * 10^(3*log10(s)) # normalized coefficients
  #   g
  # }
  # kernel coefficients at k-nearest neighbors
  # if(adaptative) {
  #   g <- t(apply(cbind(0, d), 1, sodgk))
  # } else {
  #   g <- t(apply(cbind(0, d), 1, sodgk, s = min(d[, k])))
  # }
  # p <- cbind(p, knn_values(p, i))
  # g <- rowSums(g * p) # BUG/TODO: needs to be normalized by disc areas

  # p <- S01(p)
  q <- knn_values(p, i)
  h0 <- k %/% 2
  h1 <- h0 + 1
  a <- d[, h0]^2 / (d[, k]^2 - d[, h1]^2)
  g <- rowSums(a * q[, h1:k]) - rowSums(cbind(p, q[, 1:h0]))

  if(smoothing) {
    g <- knn_smoothing(g, i)
  }

  g
}

