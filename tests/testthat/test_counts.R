context("Basic functions on read counts")

# =============================================================================.
tst <- log(c(0,NA,1,Inf))
res <- c(F,F,T,F)
# -----------------------------------------------------------------------------.
test_that("finiteValues properly detects finite values", {
  expect_equal(finiteValues(tst), res)
  expect_equal(finiteValues(cbind(tst, tst)), res)
  expect_equal(finiteValues(cbind(1, tst)), res)
  expect_equal(finiteValues(cbind(tst, 1)), res)
  expect_equal(finiteValues(cbind(1, tst, 1)), res)
})

# =============================================================================.
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
# -----------------------------------------------------------------------------.
test_that("detectCounts properly detects counts", {
  expect_equal(detectCounts(tst)$all, as.logical(res[,"all"]))
  expect_equal(detectCounts(tst)$none, as.logical(res[,"none"]))
  expect_equal(detectCounts(tst)$nbr, res[,"nbr"])
})
