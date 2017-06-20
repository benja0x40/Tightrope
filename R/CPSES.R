# =============================================================================.
#' Signal Extraction Scaling (Diaz et al., 2012)
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{CPSES},
#'   \link{BRD}
# -----------------------------------------------------------------------------.
#' @param x
#' vector of read counts in the control condition (Input, IgG, etc.)
#'
#' @param y
#' vector of read counts in the specific condition
#'
#' @return
#' scaling factor
# -----------------------------------------------------------------------------.
#' @export
Diaz.SES <- function(x, y) {

  i <- order(y)
  Y <- cumsum(y[i])
  X <- cumsum(x[i])

  n <- length(i)
  p <- Y/Y[n]
  q <- X/X[n]

  k <- which.max(abs(q - p))

  Y[k]/X[k]
}

# =============================================================================.
#' Combined Pairwise SES
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{Diaz.SES},
#'   \link{BRD}
# -----------------------------------------------------------------------------.
#' @param x
#' matrix of read counts
#' (rows = observations, columns = measurement conditions).
#'
#' @param dither
#' number of replicates for count dithering:
#' 1 = single (default), 2 = duplicate, 3 = triplicate, etc.
#'
#' @return
#' scaling factors
# -----------------------------------------------------------------------------.
#' @export
CPSES <- function(cnt, dither = 1) {

  n <- ncol(cnt)

  alpha <- matrix(diag(n), n, n, dimnames = list(colnames(x), colnames(x)))

  for(d in 1:dither) {
    x <- DitherCounts(cnt)

    # SES for all pairs
    for(i in 1:n) {
      for(j in 1:n) {
        if(i != j) {
          alpha[i, j] <- alpha[i, j] + Diaz.SES(x[, i], x[, j]) / dither
        }
      }
    }
  }

  # Combined SES
  alpha <- log(alpha)
  alpha <- colMeans(alpha)
  alpha <- exp(- alpha)

  alpha
}
# -----------------------------------------------------------------------------.
# Validations
# layout(matrix(1:6, 6, 1))
# par(mar = c(4, 4, 2, 1))
# n <- 300
#
# # TEST 1 ////
# x <- round(rnorm(n) + 15)
# y <- round(c(rnorm(n/2), 10 + rnorm(n/2)) + 10)
#
# plot(x, type = "l", ylim = range(x, y))
# points(y, type = "l", col = rgb(1, 0.5, 0))
#
# x <- x * Diaz.SES(x, y)
# plot(x, type = "l", ylim = range(x, y))
# points(y, type = "l", col = rgb(1, 0.5, 0))
#
# # TEST 2 ////
# x <- round(rnorm(n) + 15)
# z <- round(c(rnorm(n/3), 20 + rnorm(n/3), rnorm(n/3)) + 5)
#
# plot(x, type = "l", ylim = range(x, z))
# points(z, type = "l", col = rgb(1, 0.0, 0))
#
# x <- x * Diaz.SES(x, z)
# plot(x, type = "l", ylim = range(x, z))
# points(z, type = "l", col = rgb(1, 0.0, 0))
#
# # TEST 3 ////
# x <- round(rnorm(n) + 15)
# y <- round(c(rnorm(n/2), 15 + rnorm(n/2)) + 10)
# z <- round(c(rnorm(n/3), 15 + rnorm(n/3), rnorm(n/3)) + 5)
#
# plot(x, type = "l", ylim = range(x, y, z))
# points(y, type = "l", col = rgb(1, 0.5, 0))
# points(z, type = "l", col = rgb(1, 0.0, 0))
#
# alpha <- CPSES(cbind(x, y, z))
# x <- x * alpha[1]
# y <- y * alpha[2]
# z <- z * alpha[3]
# plot(x, type = "l", ylim = range(x, y, z))
# points(y, type = "l", col = rgb(1, 0.5, 0))
# points(z, type = "l", col = rgb(1, 0.0, 0))
