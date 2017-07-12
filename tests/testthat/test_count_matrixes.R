context("count_matrixes")

# =============================================================================.
#
# -----------------------------------------------------------------------------.
mapped_reads      <- "data/mapped_reads"
genomic_intervals <- "data/genomic_intervals"

mapped_reads      <- "../../data/mapped_reads"
genomic_intervals <- "../../data/genomic_intervals"

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("VerifyInputs", {

  # Preparation ////
  lst <- dir(genomic_intervals, pattern = ".rdata$", full.names = T)
  names(lst) <- gsub("\\.rdata", "", basename(lst), perl = T, ignore.case = T)
  grg <- lapply(lst, readRDS)
  cnt <- matrix(0, length(grg$GN), 5, dimnames = list(NULL, LETTERS[1:5]))
  # Test ////
  expect_true(VerifyInputs(NULL, grg$GN))
  expect_true(VerifyInputs(cnt,  grg$GN))
  expect_error(VerifyInputs(cnt, NULL))
  expect_error(VerifyInputs(cnt, grg))
  expect_error(VerifyInputs(cnt, grg$GN_ENS))

  # Preparation ////
  colnames(cnt) <- NULL
  # Test ////
  expect_error(VerifyInputs(cnt,  grg$GN))

})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("LoadMappedReads", {

  # Preparation ////
  lst <- dir(mapped_reads, pattern = ".bam$", full.names = T)
  aln <- lapply(lst, readGAlignments)
  lbl <- gsub("\\.bam", "", basename(lst), perl = T, ignore.case = T)
  # Test ////
  expect_error(LoadMappedReads(aln = NULL))
  expect_error(LoadMappedReads(aln = lbl))
  expect_error(LoadMappedReads(aln = lst))
  expect_error(LoadMappedReads(aln = aln))

  # Preparation ////
  names(lst) <- lbl
  names(aln) <- lbl
  # Test ////
  expect_identical(LoadMappedReads(lst), aln)
  expect_identical(LoadMappedReads(aln), aln)

  # Preparation ////
  names(lst) <- NULL
  # Test ////
  expect_identical(LoadMappedReads(lst, NameByFile = T), aln)

  # Preparation ////
  lst <- dir(mapped_reads, pattern = ".rdata$", full.names = T)
  names(lst) <- lbl
  names(aln) <- lbl
  # Test ////
  expect_identical(LoadMappedReads(lst), aln)
  expect_identical(LoadMappedReads(aln), aln)

  # Preparation ////
  names(lst) <- NULL
  # Test ////
  expect_identical(LoadMappedReads(lst, NameByFile = T), aln)

  # Preparation ////
  aln <- lapply(aln, as, Class = "GRanges")
  # Test ////
  expect_identical(LoadMappedReads(aln), aln)
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("AppendReadCounts", {

  # Preparation ////
  lst <- dir(genomic_intervals, pattern = ".rdata$", full.names = T)
  grg <- LoadMappedReads(lst, NameByFile = T) # ;-)
  lst <- dir(mapped_reads, pattern = ".bam$", full.names = T)
  aln <- LoadMappedReads(lst, NameByFile = T)
  nbr <- length(grg$GN)

  # Test ////
  x <- NULL
  x <- AppendReadCounts(x, aln[1:2], grg$GN)
  expect_identical(colnames(x), names(aln[1:2]))
  expect_equal(dim(x), c(nbr, 2))
  x <- AppendReadCounts(x, aln[3], grg$GN)
  expect_identical(colnames(x), names(aln[1:3]))
  expect_equal(dim(x), c(nbr, 3))

  # Test ////
  x <- NULL
  x <- AppendReadCounts(x, lst[1:2], grg$GN, NameByFile = T)
  expect_identical(colnames(x), names(aln[1:2]))
  expect_equal(dim(x), c(nbr, 2))
  x <- AppendReadCounts(x, lst[3], grg$GN, NameByFile = T)
  expect_identical(colnames(x), names(aln[1:3]))
  expect_equal(dim(x), c(nbr, 3))

  # Test ////
  x <- NULL
  x <- AppendReadCounts(x, aln[1:2], grg$GN)
  expect_error(AppendReadCounts(x, aln[2], grg$GN))
  expect_error(AppendReadCounts(x, aln[2], grg$GN_ENS))
  expect_error(AppendReadCounts(x, lst[2], grg$GN, NameByFile = T))
  expect_error(AppendReadCounts(x, lst[2], grg$GN_ENS, NameByFile = T))
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("UpdateReadCounts", {

  # Preparation ////
  lst <- dir(genomic_intervals, pattern = ".rdata$", full.names = T)
  grg <- LoadMappedReads(lst, NameByFile = T) # ;-)
  lst <- dir(mapped_reads, pattern = ".bam$", full.names = T)
  aln <- LoadMappedReads(lst, NameByFile = T)
  nbr <- length(grg$GN)
  cnt <- AppendReadCounts(NULL, aln, grg$GN)

  # Test ////
  x <- matrix(0, nbr, length(aln), dimnames = list(NULL, names(aln)))
  x <- UpdateReadCounts(x, aln[c(1, 3)], grg$GN)
  expect_identical(x[, c(1, 3)], cnt[, c(1, 3)])
  expect_true(all(x[, 2] == 0))

  # Test ////
  x <- matrix(0, nbr, length(aln), dimnames = list(NULL, names(aln)))
  x <- UpdateReadCounts(x, lst[c(1, 3)], grg$GN, NameByFile = T)
  expect_identical(x[, c(1, 3)], cnt[, c(1, 3)])
  expect_true(all(x[, 2] == 0))
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("ReadCountMatrix", {

  # Preparation ////
  lst <- dir(genomic_intervals, pattern = ".rdata$", full.names = T)
  grg <- LoadMappedReads(lst, NameByFile = T) # ;-)
  lst <- dir(mapped_reads, pattern = ".bam$", full.names = T)
  aln <- LoadMappedReads(lst, NameByFile = T)
  nbr <- length(grg$GN)
  cnt <- ReadCountMatrix(aln, grg$GN)

  # Test ////
  x <- matrix(0, nbr, length(aln), dimnames = list(NULL, names(aln)))
  x <- ReadCountMatrix(aln[c(1, 3)], grg$GN, cnt = x, update = T)
  expect_identical(x[, c(1, 3)], cnt[, c(1, 3)])
  expect_true(all(x[, 2] == 0))

  # Test ////
  x <- ReadCountMatrix(aln[2], grg$GN, cnt = x, update = T)
  expect_identical(x, cnt)

  # Test ////
  x <- matrix(0, nbr, length(aln), dimnames = list(NULL, names(aln)))
  x <- ReadCountMatrix(
    lst[c(1, 3)], grg$GN, cnt = x, update = T, NameByFile = T
  )
  expect_identical(x[, c(1, 3)], cnt[, c(1, 3)])
  expect_true(all(x[, 2] == 0))

  # Test ////
  x <- ReadCountMatrix(
    lst[2], grg$GN, cnt = x, update = T, NameByFile = T
  )
  expect_identical(x, cnt)
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("MakeReadCounts", {

  # Preparation ////
  lst <- dir(genomic_intervals, pattern = ".rdata$", full.names = T)
  grg <- LoadMappedReads(lst, NameByFile = T) # ;-)
  lst <- dir(mapped_reads, pattern = ".bam$", full.names = T)
  aln <- LoadMappedReads(lst, NameByFile = T)

  # Preparation ////
  x <- MakeReadCounts(aln[1], grg[1])
  x <- MakeReadCounts(aln[2], grg[2], cnt = x)
  x <- MakeReadCounts(aln[3], grg[3], cnt = x)

  # Test ////
  expect_identical(names(x), names(grg[1:3]))
  expect_identical(colnames(x[[1]]), names(aln[1]))
  expect_identical(colnames(x[[2]]), names(aln[2]))
  expect_identical(colnames(x[[3]]), names(aln[3]))

  # Preparation ////
  x <- MakeReadCounts(aln[1], grg[1:3])
  x <- MakeReadCounts(aln[2], grg[1:3], cnt = x)
  x <- MakeReadCounts(aln[3], grg[1:3], cnt = x)

  # Test ////
  expect_identical(names(x), names(grg[1:3]))
  expect_identical(colnames(x[[1]]), names(aln[1:3]))
  expect_identical(colnames(x[[2]]), names(aln[1:3]))
  expect_identical(colnames(x[[3]]), names(aln[1:3]))

  # Preparation ////
  y <- x
  y[[1]][, 2] <- 0
  y[[2]][, 2] <- 0
  y[[3]][, 2] <- 0
  y <- MakeReadCounts(aln[2], grg[1:3], cnt = y, update = T)

  # Test ////
  expect_identical(x[[1]][, 2], y[[1]][, 2])
  expect_identical(x[[2]][, 2], y[[2]][, 2])
  expect_identical(x[[3]][, 2], y[[3]][, 2])
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("JoinColumns", {

  # Preparation ////
  lst <- dir(genomic_intervals, pattern = ".rdata$", full.names = T)
  grg <- LoadMappedReads(lst, NameByFile = T) # ;-)
  lst <- dir(mapped_reads, pattern = ".bam$", full.names = T)
  aln <- LoadMappedReads(lst, NameByFile = T)

  # Preparation ////
  x <- MakeReadCounts(aln[1], grg[1:3])
  y <- MakeReadCounts(aln[3], grg[1:3])
  a <- JoinColumns(x, y)
  b <- MakeReadCounts(aln[c(1, 3)], grg[1:3])

  # Test ////
  expect_identical(a, b)
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("ExtractColumns", {

  # Preparation ////
  lst <- dir(genomic_intervals, pattern = ".rdata$", full.names = T)
  grg <- LoadMappedReads(lst, NameByFile = T) # ;-)
  lst <- dir(mapped_reads, pattern = ".bam$", full.names = T)
  aln <- LoadMappedReads(lst, NameByFile = T)

  # Preparation ////
  x <- MakeReadCounts(aln[1:3], grg[1:3])
  a <- ExtractColumns(x, names(aln[2]))
  b <- MakeReadCounts(aln[2], grg[1:3])

  # Test ////
  expect_identical(a, b)
})

# =============================================================================.
#
# -----------------------------------------------------------------------------.
test_that("RemoveColumns", {

  # Preparation ////
  lst <- dir(genomic_intervals, pattern = ".rdata$", full.names = T)
  grg <- LoadMappedReads(lst, NameByFile = T) # ;-)
  lst <- dir(mapped_reads, pattern = ".bam$", full.names = T)
  aln <- LoadMappedReads(lst, NameByFile = T)

  # Preparation ////
  x <- MakeReadCounts(aln[1:3], grg[1:3])
  a <- RemoveColumns(x, names(aln[c(1, 3)]))
  b <- MakeReadCounts(aln[2], grg[1:3])

  # Test ////
  expect_identical(a, b)
})
