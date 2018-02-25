# > utilities ==================================================================
context("utilities")

# + rankstat -------------------------------------------------------------------
test_that("rankstat", {
  expect_equal(rankstat(1), 0.5)
  expect_equal(rankstat(1:2), c(0.25, 0.75))
  expect_equal(rankstat(2:1), c(0.75, 0.25))
  expect_equal(range(rankstat(1:10)), c(0.05, 0.95))
})

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

# + xy2md ----------------------------------------------------------------------
test_that("xy2md", {

})

# + md2xy ----------------------------------------------------------------------
test_that("md2xy", {

})
