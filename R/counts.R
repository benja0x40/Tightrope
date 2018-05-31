# =============================================================================.
#' Apply dithering to read counts
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{NonZeroCounts},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @description
#' The \code{DitherCounts} function applies a triangular dithering filter
#' which eliminates the presence of identical values in a read count matrix.
#'
#' @param x
#' matrix of read counts (rows = observations, columns = samples or conditions).
#'
#' @return
#' \code{DitherCounts} returns a matrix of dithered counts.
# -----------------------------------------------------------------------------.
#' @export
DitherCounts <- function(x) {

  zero <- x == 0
  xmin <- min(x[! zero], na.rm = TRUE)
  xmax <- max(x[! zero], na.rm = TRUE)

  # dithering
  x <- x + rtriangle(length(x), a = -1, b = 1)
  # lower limit
  k <- x < xmin & ! zero
  x[k] <- runif(sum(k), xmin - 0.5, xmin)
  # upper limit
  k <- x > xmax
  x[k] <- runif(sum(k), xmax, xmax + 0.5)
  # no count
  x[zero] <- 0

  x
}

# HIDDEN #######################################################################

# =============================================================================.
#' Remove observations with missing values
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{DetectCounts},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @param cnt
#' matrix of read counts
#' (rows = observations, columns = measurement conditions).
#'
#' @return
#' \code{NonZeroCounts} returns a matrix.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
NonZeroCounts <- function(cnt) {
  chk <- FiniteValues(log2(cnt))
  cnt <- cnt[chk, ]
  cnt
}

# =============================================================================.
#' DetectCounts
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{NonZeroCounts},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @param cnt
#' matrix of read counts
#' (rows = observations, columns = measurement conditions).
#'
#' @param detailed
#' logical, if TRUE return min and max read counts in each interval
#'
#' @return
#' \code{DetectCounts} returns a data.frame with the following columns:
#' \item{all}{logical vector indicating rows of the \code{cnt} matrix with at least one count in all columns}
#' \item{none}{logical vector indicating rows of the \code{cnt} matrix with no counts in all columns}
#' \item{nbr}{integer vector representing the number of columns with at least one counts for each row of the \code{cnt} matrix}
#' If \code{detailed} is equal to \code{TRUE}, the returned \code{data.frame} also includes the following columns:
#' \item{min}{integer vector representing the maximum count value for each row of the \code{cnt} matrix}
#' \item{max}{integer vector representing the minimum count value for each row of the \code{cnt} matrix}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
DetectCounts <- function(cnt, detailed = FALSE) {
  cnt <- as.matrix(cnt)
  n <- ncol(cnt)
  k <- apply(cnt, MARGIN = 1, FUN = function(x) { sum(x > 0) })
  res <- data.frame(
    all  = k == n,
    none = k == 0,
    nbr  = k
  )
  if(detailed) {
    idx <- which(k > 0)
    res$min <- 0
    res$min[idx] <- apply(cnt[idx, ], 1, FUN=function(x) { min(x[x > 0]) })
    res$max <- apply(cnt, MARGIN = 1, FUN = max)
  }
  res
}
