# =============================================================================.
#' rankstat
# -----------------------------------------------------------------------------.
#' @description
#' rank statistics
#'
#' @param x
#' numeric vector
#'
#' @return
#' \eqn{(rank(x) - 0.5) / N} where N = length(x)
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
rankstat <- function(x) { (rank(x) - 0.5) / length(x) }

# =============================================================================.
#' S01
# -----------------------------------------------------------------------------.
#' @description
#' rescale x into [0, 1]
#'
#' @param x
#' numeric vector
#'
#' @return
#' S01 returns x rescaled such that range(x) = [0, 1]
# -----------------------------------------------------------------------------.
#' @keywords internal
S01 <- function(x) { (x - min(x)) / diff(range(x)) }

# =============================================================================.
#' SX2Y
# -----------------------------------------------------------------------------.
#' @description
#' rescale x into range(y)
#'
#' @param x
#' numeric vector
#'
#' @param y
#' numeric vector
#'
#' @return
#' SX2Y returns x rescaled such that range(x) = range(y)
# -----------------------------------------------------------------------------.
#' @keywords internal
SX2Y <- function(x, y) { S01(x) * diff(range(y)) + min(y) }

# =============================================================================.
#' SkewedProbabilities
# -----------------------------------------------------------------------------.
#' @param k
#' resolution
#'
#' @param beta
#' exponent
#'
#' @param gamma
#' mix
#'
#' @param sx
#' +1 or -1
#'
#' @param sy
#' +1 or -1
#'
#' @param plot
#' default = F
#'
#' @return
#' SkewedProbabilities returns a vector of probabilities
# -----------------------------------------------------------------------------.
#' @keywords internal
SkewedProbabilities <- function(
  k, beta = - 1/8, gamma = 0.5, sx = -1, sy = 1, plot = F
) {

  f <- function(a, b) { (a:b - 0.5) / max(a, b) }
  r <- function(a, b) { S01(f(a, b)^beta) }

  if(sx > 0) x <- f(1, k)
  if(sx < 0) x <- f(k, 1)
  if(sy > 0) y <- min(x) + r(1, k)
  if(sy < 0) y <- min(x) + r(k, 1)

  p <- gamma * x + (1 - gamma) * y
  p <- SX2Y(p, x)

  if(plot) {
    plot(
      0, type = 'n', xlim = c(1, k), ylim = range(x, y, p),
      xlab = "", ylab = ""
    )
    points(1:k, x, type = "l", col = grey(0.5))
    points(1:k, y, type = "l", col = grey(0.5))
    points(1:k, p, type = "l", col = rgb(1, 0, 0))
  }

  p
}

# =============================================================================.
#' MatchLogSpace
# -----------------------------------------------------------------------------.
#' @param p
#' numeric vector
#'
#' @param x
#' sorted vector of unique and strictly positive values
#'
#' @param b
#' base of the logarithm (default = 2)
#'
#' @return
#' MatchLogSpace returns a numeric vector
# -----------------------------------------------------------------------------.
#' @keywords internal
MatchLogSpace <- function(p, x, b = 2) {
  n <- length(p)
  q <- approx(x = b^(SX2Y(1:n, log(x, base = b))), y = p, xout = x, rule = 2)$y
  q
}

