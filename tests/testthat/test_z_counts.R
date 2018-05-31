# > counts =====================================================================
context("counts")

# + DitherCounts ---------------------------------------------------------------
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

# + NonZeroCounts ---------------------------------------------------------------
test_that("NonZeroCounts", {
  cnt <- rbind(
    c(3, 3, 3),
    c(0, 0, 0),
    c(0, 0, 1),
    c(0, 1, 0),
    c(1, 0, 0),
    c(2, 2, 2),
    c(0, 1, 1),
    c(1, 1, 0),
    c(1, 0, 1),
    c(1, 1, 1)
  )
  tst <- rbind(
    c(3, 3, 3),
    c(2, 2, 2),
    c(1, 1, 1)
  )
  expect_equal(NonZeroCounts(cnt), tst)
})

# + DetectCounts ---------------------------------------------------------------
test_that("DetectCounts", {
  cnt <- rbind(
    c(0, 0, 0),
    c(0, 0, 1),
    c(0, 1, 0),
    c(1, 0, 0),
    c(0, 2, 3),
    c(4, 5, 0),
    c(6, 0, 7),
    c(1, 1, 1)
  )
  tst <- cbind(
    all  = c(F, F, F, F, F, F, F, T),
    none = c(T, F, F, F, F, F, F, F),
    nbr  = c(0, 1, 1, 1, 2, 2, 2, 3),
    min  = c(0, 1, 1, 1, 2, 4, 6, 1),
    max  = c(0, 1, 1, 1, 3, 5, 7, 1)
  )

  r <- DetectCounts(cnt)
  expect_equal(r$all, as.logical(tst[,"all"]))
  expect_equal(r$none, as.logical(tst[,"none"]))
  expect_equal(r$nbr, tst[,"nbr"])

  r <- DetectCounts(cnt, detailed = TRUE)
  expect_equal(r$min, tst[,"min"])
  expect_equal(r$max, tst[,"max"])
})
