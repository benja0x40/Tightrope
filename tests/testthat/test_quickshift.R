context("quickshift")

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
test_that("quickShift", {

  layout(matrix(1:9, 3, 3, byrow = T))

  # Test ////
  n <- 30
  x <- rbind(
    translation_2D(make_square_unif_2D(n / 3), c( 0,  0)),
    translation_2D(make_square_unif_2D(n / 2), c(-3, -3)),
    translation_2D(make_square_unif_2D(n / 1), c( 4,  4))
  )
  mbr <- c(rep(3, (n / 3)^2), rep(2, (n / 2)^2), rep(1, (n / 1)^2))
  nbr <- (n / (1:3))^2

  d <- knn_density(x, k = 75)
  plot(x, pch = 20, col = grey(0, alpha = 0.5))
  plot(x, pch = 20, col = colorize(d))
  g <- suppressWarnings(quickShift(x, d, plot = F))
  grp <- quickShiftCutClusters(g, n = 3)
  clr <- rainbow(grp$nbr, alpha = 0.5)
  clr <- plot_groups_2D(x, clr = clr[grp$membership], pch = 20)

  expect_equal(grp$csize, nbr)
  expect_true(all(grp$membership == mbr))

  # Test ////
  n <- 20
  x <- rbind(
    translation_2D(make_square_unif_2D(n) / 1, c( 0,  0)),
    translation_2D(make_square_unif_2D(n) / 2, c(-2, -2)),
    translation_2D(make_square_unif_2D(n) / 3, c( 3,  3))
  )
  mbr <- rep(1:3, each = n^2)
  nbr <- rep(n^2, 3)

  d <- knn_density(x, k = 75)
  plot(x, pch = 20, col = grey(0, alpha = 0.5))
  plot(x, pch = 20, col = colorize(d))
  g <- suppressWarnings(quickShift(x, d, plot = F))
  grp <- quickShiftCutClusters(g, n = 3)
  clr <- rainbow(grp$nbr, alpha = 0.5)
  clr <- plot_groups_2D(x, clr = clr[grp$membership], pch = 20)

  expect_equal(grp$csize, nbr)
  expect_true(all(grp$membership == mbr))
})
