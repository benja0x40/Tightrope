context("counts")

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("DitherCounts", {
  # Test ////
  x <- rep(1, 1000)
  r <- DitherCounts(x)
  d <- abs(x - r)
  expect_false(any(duplicated(r)))
  expect_true(max(d) < 0.5)
  expect_equal(round(min(r), 1), min(x) - 0.5)
  expect_equal(round(max(r), 1), max(x) + 0.5)
  # Test ////
  x <- rep(5, 1000)
  r <- DitherCounts(x)
  d <- abs(x - r)
  expect_false(any(duplicated(r)))
  expect_true(max(d) < 0.5)
  expect_equal(round(min(r), 1), min(x) - 0.5)
  expect_equal(round(max(r), 1), max(x) + 0.5)
  # Test ////
  x <- rep(1:5, 1000)
  r <- DitherCounts(x)
  d <- abs(x - r)
  expect_false(any(duplicated(r)))
  expect_true(max(d) < 1.0)
  expect_equal(round(min(r), 1), min(x) - 0.5)
  expect_equal(round(max(r), 1), max(x) + 0.5)
  # Test ////
  x <- rep(5:9, 1000)
  r <- DitherCounts(x)
  d <- abs(x - r)
  expect_false(any(duplicated(r)))
  expect_true(max(d) < 1.0)
  expect_equal(round(min(r), 1), min(x) - 0.5)
  expect_equal(round(max(r), 1), max(x) + 0.5)
  # Test ////
  x <- rep(c(0, 1:5), 1000)
  k <- x == 0
  r <- DitherCounts(x)
  d <- abs(x - r)
  expect_false(any(duplicated(r[!k])))
  expect_true(max(d) < 1.0)
  expect_true(all(r[k] == 0))
  expect_equal(round(min(r[!k]), 1), min(x[!k]) - 0.5)
  expect_equal(round(max(r[!k]), 1), max(x[!k]) + 0.5)
  # Test ////
  x <- rep(c(0, 5:9), 1000)
  k <- x == 0
  r <- DitherCounts(x)
  d <- abs(x - r)
  expect_false(any(duplicated(r[!k])))
  expect_true(max(d) < 1.0)
  expect_true(all(r[k] == 0))
  expect_equal(round(min(r[!k]), 1), min(x[!k]) - 0.5)
  expect_equal(round(max(r[!k]), 1), max(x[!k]) + 0.5)
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("FiniteValues", {
  tst <- log(c(0,NA,1,Inf))
  res <- c(F,F,T,F)
  expect_equal(FiniteValues(tst), res)
  expect_equal(FiniteValues(cbind(tst, tst)), res)
  expect_equal(FiniteValues(cbind(1, tst)), res)
  expect_equal(FiniteValues(cbind(tst, 1)), res)
  expect_equal(FiniteValues(cbind(1, tst, 1)), res)
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("detectCounts", {
  tst <- rbind(
    c(0,0,0),
    c(0,0,1),
    c(0,1,0),
    c(1,0,0),
    c(0,1,1),
    c(1,1,0),
    c(1,0,1),
    c(1,1,1)
  )
  res <- cbind(
    all  = c(F,F,F,F,F,F,F,T),
    none = c(T,F,F,F,F,F,F,F),
    nbr  = c(0,1,1,1,2,2,2,3)
  )
  expect_equal(detectCounts(tst)$all, as.logical(res[,"all"]))
  expect_equal(detectCounts(tst)$none, as.logical(res[,"none"]))
  expect_equal(detectCounts(tst)$nbr, res[,"nbr"])
})
