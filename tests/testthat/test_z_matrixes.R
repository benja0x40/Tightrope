# > matrixes ===================================================================
context("matrixes")

# + m2v ------------------------------------------------------------------------
test_that("m2v", {
  a <- m2v(1:3, 1:3, nrow = 3)
  b <- c(1, 5, 9)
  expect_equal(a, b)
  a <- m2v(3:1, 1:3, nrow = 3)
  b <- c(3, 5, 7)
  expect_equal(a, b)
})

# + v2m ------------------------------------------------------------------------
test_that("v2m", {
  a <- v2m(c(1, 5, 9), nrow = 3)
  b <- cbind(1:3, 1:3)
  expect_equal(a, b)
  a <- v2m(c(3, 5, 7), nrow = 3)
  b <- cbind(3:1, 1:3)
  expect_equal(a, b)
})

# + AddByRow -------------------------------------------------------------------
test_that("AddByRow", {
  v <- 5:1
  m <- matrix(1:50, 10, 5)
  r <- t(t(m) + v)
  expect_equal(AddByRow(m, v), r)
})

# + MulByRow -------------------------------------------------------------------
test_that("MulByRow", {
  v <- 5:1
  m <- matrix(1:50, 10, 5)
  r <- t(t(m) * v)
  expect_equal(MulByRow(m, v), r)
})
