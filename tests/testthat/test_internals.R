# > internals ==================================================================
context("internals.R")

# + S01 ------------------------------------------------------------------------
test_that("S01", {
  expect_equal(range(S01(0:10)), c(0, 1))
  expect_equal(range(S01(-5:5)), c(0, 1))
})

# + SX2Y -----------------------------------------------------------------------
test_that("SX2Y", {
  a <- 0:10
  b <- -5:5
  expect_equal(range(SX2Y(a, b)), range(b))
  expect_equal(range(SX2Y(b, a)), range(a))
})

# + SkewedProbabilities --------------------------------------------------------
test_that("SkewedProbabilities", {

  # Work as expected
  # layout(matrix(1:4, 2, 2, byrow = T))
  # par(mar = c(5, 4, 3, 1))
  # p <- SkewedProbabilities(1000, beta =  1/8, gamma = 0.5, sx = +1, sy = +1, plot = T)
  # p <- SkewedProbabilities(1000, beta =  1/8, gamma = 0.5, sx = -1, sy = -1, plot = T)
  # p <- SkewedProbabilities(1000, beta = -1/8, gamma = 0.5, sx = +1, sy = -1, plot = T)
  # p <- SkewedProbabilities(1000, beta = -1/8, gamma = 0.5, sx = -1, sy = +1, plot = T)
  # Does not work as expected when gamma != 0.5
  # layout(matrix(1:4, 2, 2, byrow = T))
  # par(mar = c(5, 4, 3, 1))
  # p <- SkewedProbabilities(1000, beta =  1/2, gamma = 0.5, sx = +1, sy = -1, plot = T)
  # p <- SkewedProbabilities(1000, beta =  1/2, gamma = 0.5, sx = -1, sy = +1, plot = T)
  # p <- SkewedProbabilities(1000, beta = -1/8, gamma = 0.5, sx = +1, sy = +1, plot = T)
  # p <- SkewedProbabilities(1000, beta = -1/8, gamma = 0.5, sx = -1, sy = -1, plot = T)

})

