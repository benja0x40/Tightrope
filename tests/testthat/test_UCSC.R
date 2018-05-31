# > UCSC =======================================================================
context("UCSC")

# + BuildSeqInfo ---------------------------------------------------------------
test_that("BuildSeqInfo", {

  expect_error(BuildSeqInfo(genome = "A_Non_Existent_Genome"))

  for(asm in c("mm10", "hg38")) {
    sinf <- BuildSeqInfo(genome = asm)
    expect_is(sinf, "Seqinfo")
    expect_identical(unique(genome(sinf)), asm)
  }
})

# + ImportCpGIslandExt ---------------------------------------------------------
test_that("ImportCpGIslandExt", {

  expect_error(
    ImportCpGIslandExt(), regexp = "missing"
  )
  expect_error(
    ImportCpGIslandExt(fpath = "A.Non.Existent.File"),
    regexp = "file not found"
  )

  asm <- "mm10"
  flp <- paste0("UCSC_cpgIslandExt_", asm, ".txt.gz")

  cgi <- ImportCpGIslandExt(genome = asm, keep.file = TRUE, quiet = TRUE)
  expect_true(file.exists(flp))
  expect_true(file.remove(flp))

  cgi <- ImportCpGIslandExt(genome = asm, quiet = TRUE)
  expect_false(file.exists(flp))

  lst <- c(
    "bin", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp"
  )

  for(asm in c("mm10", "hg38")) {
    cgi <- ImportCpGIslandExt(genome = asm, quiet = TRUE)
    expect_is(cgi, "GRanges")
    expect_identical(unique(genome(cgi)), asm)
    expect_true(length(cgi) > 16000)
    expect_identical(colnames(mcols(cgi)), lst)
  }
})

