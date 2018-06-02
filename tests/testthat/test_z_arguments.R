# > Arguments ==================================================================
context("Arguments")

# + DefaultArgs -----------------------------------------------------------------
test_that("DefaultArgs", {

  cfg <- list(x = 1, y = 2, z = 3)

  src <- new.env()
  dst <- new.env()

  DefaultArgs(cfg, from = src)
  expect_identical(names(src), character(0))

  DefaultArgs(cfg, to = dst)
  expect_identical(as.list(dst), cfg)

  dst <- new.env()
  dst$y <- 0
  DefaultArgs(cfg, to = dst)
  expect_identical(as.list(dst), list(x = 1, y = 0, z = 3))

  alt <- list(x = 3, y = 2, z = 1)
  dst <- as.environment(alt)

  DefaultArgs(cfg, to = dst)
  expect_identical(as.list(dst, sorted = TRUE), alt)

  src <- new.env()
  src$y <- 0

  DefaultArgs(cfg, from = src, to = dst)
  expect_identical(as.list(dst, sorted = TRUE), alt)

  dst <- new.env()
  DefaultArgs(cfg, from = src, to = dst)
  expect_identical(as.list(dst), list(x = 1, y = 0, z = 3))

  f <- function(x = NULL, y = NULL, z = NULL, ...) {
    DefaultArgs(cfg)
    list(x = x, y = y, z = z)
  }

  expect_identical(f(), cfg)
  expect_identical(f(x = 0)$x, 0)
  expect_identical(f(y = 0)$y, 0)
  expect_identical(f(z = 0)$z, 0)
  expect_identical(f(i = 1:10), cfg)

  f <- function(x = NULL, y = NULL, z = NULL, ...) {
    DefaultArgs(cfg, ignore = "...", from = f)
    list(x = x, y = y, z = z)
  }

  expect_identical(f(), cfg)
  expect_identical(f(x = 0)$x, 0)
  expect_identical(f(y = 0)$y, 0)
  expect_identical(f(z = 0)$z, 0)
  expect_identical(f(i = 1:10), cfg)

})

# + VectorArgs ------------------------------------------------------------------
test_that("VectorArgs", {

  x <- 0
  y <- 1:10

  VectorArgs(c("x", "y"))
  expect_identical(x, rep(0, 10))
  expect_identical(y, 1:10)

  VectorArgs(c("x", "y"), size = 15)
  expect_identical(x, rep(0, 15))
  expect_identical(y[11:15], 1:5)

  a <- list(x = 0, y = 1:10)

  r <- VectorArgs(c("x", "y"), from = a)
  expect_identical(r$x, rep(0, 10))
  expect_identical(r$y, 1:10)

  r <- VectorArgs(c("x", "y"), from = a, size = 15)
  expect_identical(r$x, rep(0, 15))
  expect_identical(r$y[11:15], 1:5)

})


# + ClonalArg ------------------------------------------------------------------
test_that("ClonalArg", {

  a <- c("x", "y")

  d <- "i"
  v <- "j"
  r <- ClonalArg(u = v, a, d)
  expect_identical(r, list(x = v, y = v))

  d <- LETTERS[1:5]
  v <- LETTERS[5:1]
  r <- ClonalArg(u = v, a, d)
  expect_identical(r, list(x = v, y = v))

  d <- matrix("I", 2, 2)
  v <- matrix("J", 2, 2)
  r <- ClonalArg(u = v, a, d)
  expect_identical(r, list(x = v, y = v))

  d <- 1
  v <- 0
  r <- ClonalArg(u = 0, a, d)
  expect_identical(r, list(x = v, y = v))

  r <- ClonalArg(u = c(x = 0), a, d)
  expect_identical(r, list(x = v, y = d))

  r <- ClonalArg(u = c(y = 0), a, d)
  expect_identical(r, list(x = d, y = v))

  r <- ClonalArg(u = c(x = v, y = v), a, d)
  expect_identical(r, list(x = v, y = v))

  r <- ClonalArg(u = list(x = v), a, d)
  expect_identical(r, list(x = v, y = d))

  r <- ClonalArg(u = list(y = v), a, d)
  expect_identical(r, list(x = d, y = v))

  r <- ClonalArg(u = list(x = v, y = v), a, d)
  expect_identical(r, list(x = v, y = v))

  d <- 1:2
  v <- c(0, 0)

  r <- ClonalArg(u = 0, a, d)
  expect_identical(r, list(x = v, y = v))

  r <- ClonalArg(u = c(x = 0), a, d)
  expect_identical(r, list(x = v, y = d))

  r <- ClonalArg(u = c(y = 0), a, d)
  expect_identical(r, list(x = d, y = v))

  r <- ClonalArg(u = list(x = v), a, d)
  expect_identical(r, list(x = v, y = d))

  r <- ClonalArg(u = list(y = v), a, d)
  expect_identical(r, list(x = d, y = v))

  d <- matrix(1, 2, 2)
  v <- matrix(0, 2, 2)

  r <- ClonalArg(u = 0, a, d)
  expect_identical(r, list(x = v, y = v))

  r <- ClonalArg(u = c(x = 0), a, d)
  expect_identical(r, list(x = v, y = d))

  r <- ClonalArg(u = c(y = 0), a, d)
  expect_identical(r, list(x = d, y = v))

  r <- ClonalArg(u = list(x = v), a, d)
  expect_identical(r, list(x = v, y = d))

  r <- ClonalArg(u = list(y = v), a, d)
  expect_identical(r, list(x = d, y = v))

})

# Tightrope ####################################################################
