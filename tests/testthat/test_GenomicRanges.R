# > GenomicRanges ==============================================================
context("GenomicRanges")

# + GenomicTiling --------------------------------------------------------------
test_that("GenomicTiling", {

  asm <- "AlphaBetaGamma"

  gnm <-Seqinfo(
    seqnames = LETTERS[1:3],
    seqlengths = c(1000, 2000, 3000),
    genome = asm
  )

  grg <- GenomicTiling(gnm, 100)

  expect_is(grg, "GRanges")
  expect_identical(unique(genome(grg)), asm)
  expect_equal(length(grg), 60)
  expect_true(all(width(grg) == 100))
  expect_identical(
    as.vector(seqnames(grg)), c(rep("A", 10), rep("B", 20), rep("C", 30))
  )

  grg <- GenomicTiling(gnm, 100, w = 50)

  expect_equal(length(grg), 60)
  expect_true(all(width(grg) == 50))

})

# + CleanupGRanges --------------------------------------------------------------
test_that("CleanupGRanges", {

  asm <- "hg38"
  cgi <- ImportCpGIslandExt(genome = asm, quiet = TRUE)
  grg <- CleanupGRanges(cgi, seqinfo(cgi), "Homo sapiens")

  expect_true(any(grepl("_", seqlevels(cgi))))
  expect_false(any(grepl("_", seqlevels(grg))))
  expect_equal(length(seqlevels(grg)), 24)

})

