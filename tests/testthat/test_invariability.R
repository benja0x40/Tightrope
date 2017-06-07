context("invariability")

# FUNCTIONS ####################################################################

# =============================================================================.
#
# -----------------------------------------------------------------------------.
translation_2D <- function(v, dv) {
  v + matrix(dv, nrow(v), length(dv), byrow = T)
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
make_square_grid_2D <- function(n, normalized = T) {
  x <- rep(1:n, n) - (n + 1) / 2
  y <- rep(1:n, each = n) - (n + 1) / 2
  xy <- cbind(x, y)
  if(normalized) xy <- 2 / n * xy
  xy
}
make_disk_grid_2D <- function(n, normalized = T) {
  x <- rep(1:n, n) - (n + 1) / 2
  y <- rep(1:n, each = n) - (n + 1) / 2
  d <- sqrt(x^2 + y^2)
  xy <- cbind(x, y)[d < n / 2, ]
  if(normalized) xy <- 2 / n * xy
  xy
}
make_ring_grid_2D <- function(n, normalized = T) {
  x <- rep(1:n, n) - (n + 1) / 2
  y <- rep(1:n, each = n) - (n + 1) / 2
  d <- sqrt(x^2 + y^2)
  xy <- cbind(x, y)[d < n / 2 & d > n / 4, ]
  if(normalized) xy <- 2 / n * xy
  xy
}
# =============================================================================.
#
# -----------------------------------------------------------------------------.
make_square_unif_2D <- function(n, normalized = T) {
  x <- n * (runif(n * n) - 1 / 2)
  y <- n * (runif(n * n) - 1 / 2)
  xy <- cbind(x, y)
  if(normalized) xy <- 2 / n * xy
  xy
}
make_disk_unif_2D <- function(n, normalized = T) {
  x <- n * (runif(n * n) - 1 / 2)
  y <- n * (runif(n * n) - 1 / 2)
  d <- sqrt(x^2 + y^2)
  xy <- cbind(x, y)[d < n / 2, ]
  if(normalized) xy <- 2 / n * xy
  xy
}
make_ring_unif_2D <- function(n, normalized = T) {
  x <- n * (runif(n * n) - 1 / 2)
  y <- n * (runif(n * n) - 1 / 2)
  d <- sqrt(x^2 + y^2)
  xy <- cbind(x, y)[d < n / 2 & d > n / 4, ]
  if(normalized) xy <- 2 / n * xy
  xy
}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
# n <- 100
# layout(matrix(1:9, 3, 3, byrow = T))
# plot(make_square_grid_2D(n), pch = 20, col = grey(0, alpha = 0.5))
# plot(make_disk_grid_2D(n), pch = 20, col = grey(0, alpha = 0.5))
# plot(make_ring_grid_2D(n), pch = 20, col = grey(0, alpha = 0.5))
# plot(make_square_unif_2D(n), pch = 20, col = grey(0, alpha = 0.5))
# plot(make_disk_unif_2D(n), pch = 20, col = grey(0, alpha = 0.5))
# plot(make_ring_unif_2D(n), pch = 20, col = grey(0, alpha = 0.5))

# TESTS ########################################################################

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("populations", {

  layout(matrix(1:9, 3, 3, byrow = T))

  # Test ////
  n <- 30
  lst <- list(
    translation_2D(make_disk_unif_2D(n / 3), c( 0,  0)),
    translation_2D(make_disk_unif_2D(n / 2), c(-3, -3)),
    translation_2D(make_disk_unif_2D(n / 1), c( 4,  4))
  )
  x <- do.call(rbind, lst)
  nbr <- rev(sapply(lst, nrow))
  mbr <- c(rep(3, nbr[3]), rep(2, nbr[2]), rep(1, nbr[1]))
  ctr <- sapply(rev(lst), colMeans)

  d <- knn_stat(x, k = 100)
  grp <- populations(x, p = d, ns = 3)
  mu <- sapply(grp$theta, "[[", "mu")

  plot(x, pch = 20, col = grey(0, alpha = 0.5))
  plot(x, pch = 20, col = colorize(d))
  clr <- rainbow(grp$nbr, alpha = 0.5)
  clr <- plot_groups_2D(x, clr = clr[grp$membership], pch = 20)
  plot_sources_2D(
    grp$theta, components = 1:2, clr = grey(0), lwd = 2, center = F
  )

  expect_equal(grp$csize, nbr)
  expect_true(all(grp$membership == mbr))
  expect_true(mean(abs(ctr - mu)) < 1E-6)

  # Test ////
  n <- 20
  lst <- list(
    translation_2D(make_disk_unif_2D(n) / 1, c( 0,  0)),
    translation_2D(make_disk_unif_2D(n) / 2, c(-2, -2)),
    translation_2D(make_disk_unif_2D(n) / 3, c( 3,  3))
  )
  nbr <- min(sapply(lst, nrow))
  lst <- lapply(lst, function(x) x[1:nbr, ])
  x <- do.call(rbind, lst)
  nbr <- sapply(lst, nrow)
  mbr <- c(rep(1, nbr[1]), rep(2, nbr[2]), rep(3, nbr[3]))
  ctr <- sapply(lst, colMeans)

  d <- knn_stat(x, k = 100)
  grp <- populations(x, p = d, ns = 3)
  mu <- sapply(grp$theta, "[[", "mu")

  plot(x, pch = 20, col = grey(0, alpha = 0.5))
  plot(x, pch = 20, col = colorize(d))
  clr <- rainbow(grp$nbr, alpha = 0.5)
  clr <- plot_groups_2D(x, clr = clr[grp$membership], pch = 20)
  plot_sources_2D(
    grp$theta, components = 1:2, clr = grey(0), lwd = 2, center = F
  )

  expect_equal(grp$csize, nbr)
  expect_true(all(grp$membership == mbr))
  expect_true(mean(abs(ctr - mu)) < 1E-6)
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("bg_grp", {
  layout(matrix(1:9, 3, 3, byrow = T))

  # Test ////
  n <- 5000
  # Generate random samples from a mixture of normal distributions
  prm <- list(
    list(alpha = 0.3, mu = c( 1.5, -1.5), sigma = matrix(c(1.5, 0.5, 0.5, 0.25), 2)),
    list(alpha = 0.5, mu = c( 0.0,  0.0), sigma = matrix(c(3.0, 1.0, 1.0, 1.0), 2)),
    list(alpha = 0.2, mu = c(-2.0,  2.0), sigma = matrix(c(2.0, 0.5, 0.5, 1.0), 2))
  )
  spl <- mv_rsg(n = n, theta = prm)
  #
  lst <- list(
    spl[spl$group == 1, -1],
    spl[spl$group == 2, -1],
    spl[spl$group == 3, -1]
  )
  x <- do.call(rbind, lst)
  nbr <- sapply(lst, nrow)[c(2, 1, 3)]
  mbr <- c(rep(2, nbr[2]), rep(1, nbr[1]), rep(3, nbr[3]))
  ctr <- sapply(lst, colMeans)[, c(2, 1, 3)]

  d <- knn_stat(x, k = 100)
  grp <- populations(x, p = d, ns = 3)
  mu <- sapply(grp$theta, "[[", "mu")
  bg <- bg_grp(x, grp, bg_cols = 1)

  plot(x, pch = 20, col = grey(0, alpha = 0.5))
  plot(x, pch = 20, col = colorize(d))
  clr <- rainbow(grp$nbr, alpha = 0.5)
  clr <- plot_groups_2D(x, clr = clr[grp$membership], pch = 20)
  plot_sources_2D(
    grp$theta, components = 1:2, clr = grey(0), lwd = 2, center = F
  )
  points(matrix(grp$theta[[bg]]$mu, nrow = 1), pch = 19)

  expect_true(mean(abs(ctr - mu)) < 0.5)

})
