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
#' @export
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

