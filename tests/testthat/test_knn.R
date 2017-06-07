context("knn")

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
test_that("knn_sd", {
  # Test
  x       <- NULL
  nn.idx  <- NULL
  nn.dist <- cbind(-(1:10), 0, 1:10) # negative values are ok for testing
  k       <- NULL
  r <- knn_sd(x, nn.idx, nn.dist, k)
  expect_true(all(r == 1:10))
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("knn_density", {
  # Test
  x       <- NULL
  nn.idx  <- NULL
  nn.dist <- NULL
  k       <- NULL
})
