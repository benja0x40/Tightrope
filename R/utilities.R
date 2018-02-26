# HIDDEN #######################################################################

# =============================================================================.
#' HostName
# -----------------------------------------------------------------------------.
#' @return character
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
HostName <- function() {
  as.character(Sys.info()["nodename"])
}

# =============================================================================.
#' Current time
# -----------------------------------------------------------------------------.
#' @return character
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
clock <- function() {
  format(Sys.time(), '%H:%M:%S')
}

# =============================================================================.
#' Alias for mapply with USE.NAMES = F
# -----------------------------------------------------------------------------.
#' @param ...
#' see \link{mapply}.
#'
#' @return
#' see \link{mapply}.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
fmap <- function(...) {
  mapply(..., USE.NAMES = F)
}

# =============================================================================.
#' Make a destination path creating folders when necessary
# -----------------------------------------------------------------------------.
#' @param x
#' character value indicating a file system path.
#'
#' @return
#' see \link{system}.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
mkdir <- function(x) {
  system(paste0('mkdir -p "', x, '"'))
}

# =============================================================================.
#' Rescale x non-linearly into ]0, 1[
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{S01},
#'   \link{SX2Y}
# -----------------------------------------------------------------------------.
#' @description
#' rescale values non-linearly to the unit interval using rank scores calculated
#' by \eqn{q = (rank(x) - 0.5) / N} where N = length(x).
#'
#' @param x
#' numeric vector
#'
#' @return
#' \code{rankstat} returns rank scores
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
rankstat <- function(x) { (rank(x) - 0.5) / length(x) }

# =============================================================================.
#' Rescale x linearly into [0, 1]
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{rankstat},
#'   \link{SX2Y}
# -----------------------------------------------------------------------------.
#' @description
#' rescale values linearly to the unit interval.
#'
#' @param x
#' numeric vector.
#'
#' @return
#' \code{S01} returns x linearly rescaled such that range(x) = [0, 1]
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
S01 <- function(x) { (x - min(x)) / diff(range(x)) }

# =============================================================================.
#' SX2Y
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{rankstat},
#'   \link{S01}
# -----------------------------------------------------------------------------.
#' @description
#' rescale x values linearly to match the range of y.
#'
#' @param x
#' numeric vector
#'
#' @param y
#' numeric vector
#'
#' @return
#' \code{SX2Y} returns x rescaled such that range(x) = range(y).
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
SX2Y <- function(x, y) { S01(x) * diff(range(y)) + min(y) }

# =============================================================================.
#' Convert (x, y) to (a = mean, m = difference) coordinates
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{md2xy}
# -----------------------------------------------------------------------------.
#' @description
#' compute \eqn{a = (y + x)/2}, and \eqn{m = (y - x)}
#' The inverse transformation of \code{xy2md} is \link{md2xy}.
#'
#' @param x
#' numeric vector, or a matrix with 2 columns or a list containing two vectors.
#'
#' @param y
#' numeric vector.
#'
#' @return
#' \code{xy2md} returns a list or matrix with (a = mean, m = difference) values.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
xy2md <- function(x, y = NULL) {

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
#' Convert (a = mean, m = difference) to (x, y) coordinates
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{xy2md}
# -----------------------------------------------------------------------------.
#' @description
#' compute \eqn{x = a + m/2}, and \eqn{y = a - m/2}.
#' The inverse transformation of \code{md2xy} is \link{xy2md}.
#'
#' @param a
#' numeric vector, or a matrix with 2 columns or a list containing two vectors.
#'
#' @param m
#' numeric vector
#'
#' @return
#' \code{md2xy} returns a list or matrix with (x, y) values.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
md2xy <- function(a, m = NULL) {

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
