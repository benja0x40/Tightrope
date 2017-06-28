# =============================================================================.
#' knn in 1D (SLOW + ALGORITHM NOT FULLY WORKING)
# -----------------------------------------------------------------------------.
#' @param data
#' numeric matrix representing multivariate data where rows = observations
#' and columns = measurement conditions.
#'
#' @param k
#' number of nearest neighbors.
#'
#' @return
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
get.knn_1D <- function(data, k) {

  N  <- length(data) # number of observations

  o <- order(data)
  h <- data[o]

  d  <- matrix(0, N, k)
  da.1 <- c(h[2:N] - h[1:(N - 1)], Inf)
  db.1 <- c(Inf, h[2:N] - h[1:(N - 1)])

  da <- da.1
  db <- db.1

  i  <- matrix(0, N, k)
  ia <- 1:N
  ib <- 1:N

  for(j in 1:k) {
    # find nn either below or above
    ma <- which(da <  db)
    mb <- which(db <= da)
    # update nn index
    ia[ma] <- ia[ma] + 1
    ib[mb] <- ib[mb] - 1
    # save nn index
    i[ma, j] <- ia[ma]
    i[mb, j] <- ib[mb]
    # save nn distance
    d[ma, j] <- da[ma]
    d[mb, j] <- db[mb]
    # update nn distances
    da[ma] <- da[ma] + da.1[ia[ma]]
    db[mb] <- db[mb] + db.1[ib[mb]]
  }
  i[o, ] <- i
  d[o, ] <- d

  list(nn.index = i, nn.dist = d)
}
