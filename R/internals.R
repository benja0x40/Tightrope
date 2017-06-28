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
#' @seealso
#'   \link{SX2Y}
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
#' @seealso
#'   \link{S01}
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
#' @export
SX2Y <- function(x, y) { S01(x) * diff(range(y)) + min(y) }

# =============================================================================.
#' Converts log2(counts) to AM values
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{ma2lc}
# -----------------------------------------------------------------------------.
#' @description
#' a = (y + x)/2, m = (y - x)
#' @param x numeric vector, or matrix or list
#' @param y numeric vector
#'
#' @return \code{list} with a and m values
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
lc2ma <- function(x, y = NULL) {

  chk  <- 0

  if(is.null(y) & ! is.null(dim(x))) {
    y <- x[,2]
    x <- x[,1]
    chk <- 1
  }
  if(is.null(y) & is.list(x)) {
    y <- x[[2]]
    x <- x[[1]]
  }

  a <- (y + x)/2
  m <- (y - x)

  if(chk == 0) r <- list(a = a, m = m)
  if(chk == 1) r <- cbind(a = a, m = m)

  r
}

# =============================================================================.
#' Converts AM values to log2(counts)
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{lc2ma}
# -----------------------------------------------------------------------------.
#' @description
#' x = a + m/2, y = a - m/2
#' @param a numeric vector, or matrix or list
#' @param m numeric vector
#'
#' @return \code{list} with x and y values
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ma2lc <- function(a, m = NULL) {

  chk  <- 0

  if(is.null(m) & ! is.null(dim(a))) {
    m <- a[,2]
    a <- a[,1]
    chk <- 1
  }
  if(is.null(m) & is.list(a)) {
    m <- a[[2]]
    a <- a[[1]]
  }

  x <- a + m/2
  y <- a - m/2

  if(chk == 0) r <- list(x = x, y = y)
  if(chk == 1) r <- cbind(x = x, y = y)

  r
}
