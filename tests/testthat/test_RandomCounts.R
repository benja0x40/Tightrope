# =============================================================================.
context("RandomCounts")
# -----------------------------------------------------------------------------.
test_that("RandomCounts", {

  # Test: errors ////
  msg = "provide at least one of the 'k', 'p' or 'v' arguments"
  expect_error(RandomCounts(n = 10, d = 1, k = NULL, p = NULL, v = NULL), msg)

  msg = "generated counts must have at least two distinct values"
  expect_error(RandomCounts(n = 10, d = 1, k = 1,    p = NULL, v = NULL), msg)
  expect_error(RandomCounts(n = 10, d = 1, k = NULL, p = 1,    v = NULL), msg)
  expect_error(RandomCounts(n = 10, d = 1, k = NULL, p = NULL, v = 1), msg)

  msg = "arguments 'p' and 'v' have different lengths"
  expect_error(RandomCounts(n = 10, d = 1, k = NULL, p = 1:3,  v = 1:9), msg)
  expect_error(RandomCounts(n = 10, d = 1, k = NULL, p = 1:9,  v = 1:3), msg)
  expect_error(RandomCounts(n = 10, d = 1, k = 100,  p = 1:3,  v = 1:9), msg)
  expect_error(RandomCounts(n = 10, d = 1, k = 100,  p = 1:9,  v = 1:3), msg)

  expect_equal(RandomCounts(n = 10, d = 1, k = 100,  p = 1:9,  v = 1:9)$k, 9)

  # Test: p == NULL, v == NULL ////
  n <- 5
  k <- n
  p <- NULL
  v <- NULL
  for(d in 1:3) {
    r <- RandomCounts(n = n, d = d, k = k, p = p, v = v)
    expect_identical(r$k, n)
    expect_identical(r$p, rep(1 / n, n))
    expect_identical(r$v, 1:n)
  }

  # Test: p != NULL, v == NULL ////
  n <- 5
  k <- NULL
  p <- rep(c(0.1, 0.9), length.out = n)
  v <- NULL
  for(d in 1:3) {
    r <- RandomCounts(n = n, d = d, k = k, p = p, v = v)
    expect_equal(nrow(r$x), n)
    expect_equal(ncol(r$x), d)
    expect_equal(r$k, n)
    expect_identical(r$p, p / sum(p))
    expect_identical(r$v, 1:n)
  }

  # Test: p == NULL, v != NULL ////
  n <- 5
  k <- NULL
  p <- NULL
  v <- 1:n + 5
  for(d in 1:3) {
    r <- RandomCounts(n = n, d = d, k = k, p = p, v = v)
    expect_equal(nrow(r$x), n)
    expect_equal(ncol(r$x), d)
    expect_equal(r$k, n)
    expect_identical(r$p, rep(1 / n, n))
    expect_identical(r$v, v)
  }

  # Test: p != NULL, v != NULL ////
  n <- 5
  k <- NULL
  p <- rep(c(0.1, 0.9), length.out = n)
  v <- 1:n + 5
  for(d in 1:3) {
    r <- RandomCounts(n = n, d = d, k = k, p = p, v = v)
    expect_equal(nrow(r$x), n)
    expect_equal(ncol(r$x), d)
    expect_equal(r$k, n)
    expect_identical(r$p, p / sum(p))
    expect_identical(r$v, v)
  }

})
