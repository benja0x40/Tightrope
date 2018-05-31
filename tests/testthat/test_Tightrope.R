# > Tightrope ================================================================
context("Tightrope")

# + DefaultOptions -------------------------------------------------------------
test_that("DefaultOptions", {

  expect_identical(Tightrope::DefaultOptions(), Tightrope())

})

# + Arbitrary options ----------------------------------------------------------
test_that("Arbitrary options", {

  cfg <- Tightrope()
  opt <- names(cfg)

  x <- "SomeValue"
  a <- list()
  for(k in opt) {
    a[[k]] <- x
    do.call(Tightrope, a)
    expect_identical(Tightrope()[[k]], x, info = k)
    a[[k]] <- NULL
  }

  do.call(Tightrope, cfg) # Restore default values
  expect_identical(cfg, Tightrope())

})

# + Remove & Reset -------------------------------------------------------------
test_that("Remove & Reset", {

  Tightrope::RemoveOptions()
  expect_error(Tightrope())
  Tightrope::ResetOptions()
  expect_identical(Tightrope::DefaultOptions(), Tightrope())

})
