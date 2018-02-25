# =============================================================================.
#' Random counts
# -----------------------------------------------------------------------------.
#' @description
#' Generate pseudo random count values
#'
#' @param n
#' number of observations (rows).
#'
#' @param d
#' number of dimensions (columns).
#'
#' @param k
#' by default count values are sampled in [1, \code{k}].
#' This argument is ignored when arguments \code{p} or \code{v} are provided.
#'
#' @param p
#' vector of count probabilities (default = uniform).
#'
#' @param v
#' vector of count values (default = \code{1:k})
#'
#' @param extended
#' logical (default = F)
#'
#' @return
#' \code{RandomCounts} returns a list.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
RandomCounts <- function(n, d, k = NULL, p = NULL, v = NULL, extended = F) {

  chk <- sum(2^(0:2) * (! sapply(list(k, p, v), is.null)))

  if(chk == 0) {
    stop("provide at least one of the 'k', 'p' or 'v' arguments")
  }

  # p != NULL, v == NULL
  if(chk %in% (0:1 + 2)) k <- length(p)

  # v != NULL, p == NULL
  if(chk %in% (0:1 + 4)) k <- length(v)

  # k == NULL
  if(chk %in% (0:1 + 6)) k <- length(p) # assuming length(p) == length(v)

  # p == NULL
  if(! bitwAnd(chk, 2)) p <- rep(1, k)

  # v == NULL
  if(! bitwAnd(chk, 4)) v <- 1:k

  if(k != length(p) | k != length(v)) {
    stop("arguments 'p' and 'v' have different lengths")
  }
  if(k < 2) {
    stop("generated counts must have at least two distinct values")
  }

  p <- p / sum(p)

  x <- matrix(
    sample(1:k, size = n * d, replace = T, prob = p),
    nrow = n, ncol = d
  )

  # Value substitution and result
  x <- matrix(v[x], nrow = n, ncol = d)
  r <- list(x = x, k = k, p = p, v = v)

  if(extended) {
    m <- cbind(
      rep(1:k, k),
      rep(1:k, each = k)
    )
    z <- rowMeans(apply(m, 2, function(i) p[i]))
    z <- z / sum(z)
    # Value substitution and result
    m <- matrix(v[m], k^2, 2)
    r <- c(r, list(m = m, d = z))
  }

  r
}
