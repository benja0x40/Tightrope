# =============================================================================.
# TODO: check if this can be faster and as/more flexible using summerizeOverlaps
# -----------------------------------------------------------------------------.
#' Count mapped reads overlapping with genomic intervals
# -----------------------------------------------------------------------------.
#' @param \code{aln} list of \link{GAlignments} objects with mapped reads from different samples
#' @param \code{grg} \link{GRanges} object defining the considered genomic intervals
# -----------------------------------------------------------------------------.
#' @return \code{matrix} of read counts (columns = samples, rows = genomic intervals)
# -----------------------------------------------------------------------------.
makeReadCountMatrix <- function(aln, grg, ...) {
  cnt <- matrix(0, length(grg), length(aln))
  for(i in 1:length(aln)) {
    message("Read counts for sample ", i)
    cnt[,i] <- countOverlaps(
      query=grg, subject=aln[[i]], ...
    )
  }
  cnt
}
# =============================================================================.
#' Read count matrix summary
# -----------------------------------------------------------------------------.
#' @param \code{cnt} matrix of read counts (columns = samples, rows = genomic intervals)
#' @param \code{detailed} logical, if TRUE return min and max read counts in each interval
# -----------------------------------------------------------------------------.
#' @return
#' detectCounts returns a \code{data.frame} with the following columns:
#' \item{all}{logical vector indicating rows of the \code{cnt} matrix with at least one count in all columns}
#' \item{none}{logical vector indicating rows of the \code{cnt} matrix with no counts in all columns}
#' \item{nbr}{integer vector representing the number of columns with at least one counts for each row of the \code{cnt} matrix}
#' If \code{detailed} is equal to \code{TRUE}, the returned \code{data.frame} also includes the following columns:
#' \item{min}{integer vector representing the maximum count value for each row of the \code{cnt} matrix}
#' \item{max}{integer vector representing the minimum count value for each row of the \code{cnt} matrix}
# -----------------------------------------------------------------------------.
detectCounts <- function(cnt, detailed=F) {
  cnt <- as.matrix(cnt)
  n <- ncol(cnt)
  k <- apply(cnt, MARGIN=1, FUN=function(x) { sum(x>0) })
  res <- data.frame(
    all  = k == n,
    none = k == 0,
    nbr  = k
  )
  if(detailed) {
    idx <- which(k > 0)
    res$min <- 0
    res$min[idx] <- apply(cnt[idx,], 1, FUN=function(x) { min(x[x>0]) })
    res$max <- apply(cnt, MARGIN=1, FUN=max)
  }
  res
}
# =============================================================================.
#' Find finite values consistently in all samples (log transformed counts)
# -----------------------------------------------------------------------------.
#' @param \code{x}
#'
#' @return \code{logical} vector
# -----------------------------------------------------------------------------.
finiteValues <- function(x) {
  if(is.null(dim(x))) {
    x <- sapply(x, FUN=is.finite)
  } else {
    n <- ncol(x)
    x <- t(apply(x, MARGIN=1, FUN=is.finite))
    x <- rowSums(x) == n
  }
  x
}
# =============================================================================.
#' Dither read counts
# -----------------------------------------------------------------------------.
#' @param \code{cnt} matrix of read counts (columns = samples, rows = intervals)
#' @param \code{pdf} probability density of dithering values, either "triangle" (default) or "normal"
#' @param \code{fwhm} full width at half maximum for normal pdf, default = 1
#' @param \code{min.cnt} minimum count value after dithering, default = 0.5
#' @param \code{missing.cnt} logical, if TRUE generate random values for a fraction of null counts
#'
#' @return \code{matrix}
# -----------------------------------------------------------------------------.
ditherCounts <- function(
  cnt, pdf=c("triangle", "normal"), fwhm=1, min.cnt=0.5, missing.cnt=F
) {
  pdf <- pdf[1]
  x <- cnt
  n <- ncol(x)
  m <- nrow(x)
  nd <- matrix(detectCounts(x)$none, m, n)
  chk0 <- x == 0
  chk1 <- x == 1

  if(pdf == "triangle") { # Linearly distributed noise
    x <- x + rtriangle(length(x), a=-1, b=1)
    # No data
    x[chk0] <- 0
    # Random missing values
    if(missing.cnt) {
      idx <- which(chk0 & ! nd)
      idx <- sample(idx, size = sum(chk1)%/%2)
      x[idx] <- rtriangle(length(idx), a=0, b=1, c=0)
    }
  }
  if(pdf == "normal") { # Normally distributed noise
    sigma <- fwhm / (2 * sqrt(2 * log(2)))
    x <- x + rnorm(n * m, 0, sigma)
    # No data
    x[chk0] <- 0
    # Random missing values
    if(missing.cnt) {
      idx <- which(chk0 & ! nd)
      idx <- sample(idx, size = sum(chk1)%/%2)
      x[idx] <- abs(rnorm(length(idx), 0, sigma))
    }
  }
  x[x < min.cnt] <- 0
  cnt <- x
  cnt
}
# =============================================================================.
#' Compute log2 transformed counts with different methods
# -----------------------------------------------------------------------------.
#' @param \code{cnt} matrix of read counts (columns = samples, rows = intervals)
#' @param \code{lib.size}
#' @param \code{chip.con} ChIP DNA concentration, default = 1
#' @param \code{feat.len}
#' @param \code{method} either "raw" (default), "cpm", "rpkm", "voom"
#' @param \code{pc} pseudo count, default = 0
#'
#' @return \code{matrix} of log-transformed counts
# -----------------------------------------------------------------------------.
logCounts <- function(cnt, lib.size=1E6, chip.con=1, feat.len=1E3, method="raw", pc=0) {
  cnt <- as.matrix(cnt)
  lbl <- colnames(cnt)
  k <- as.vector(as.matrix(chip.con * 1e+06/lib.size))
  w <- as.vector(as.matrix(1e+03/feat.len))
  if(method=="raw") {
    cnt <- cnt + pc
  }
  if(method=="cpm") {
    cnt <- t(k * t(cnt + pc))
  }
  if(method=="rpkm") {
    cnt <- w * t(k * t(cnt + pc))
  }
  if(method=="voom") {
    s <- as.vector(as.matrix(1e+06/(lib.size + 1)))
    cnt <- t(s * t(cnt + 0.5))
  }
  cnt <- log2(cnt)
  colnames(cnt) <- lbl
  cnt
}
# =============================================================================.
#' Converts log2(counts) to AM values
#' @description
#' a = (y + x)/2, m = (y - x)
# -----------------------------------------------------------------------------.
#' @param \code{x}
#' @param \code{y}
#'
#' @return \code{list} with a and m values
# -----------------------------------------------------------------------------.
lc2ma <- function(x, y = NULL) {

  chk  <- 0

  if(is.null(y) & ! is.null(dim(x))) {
    y <- x[,2]
    x <- x[,1]
    chk <- 1
  }
  if(is.null(y) & is.list(x)) {
    y <- x[[2]]
    x <- x[[1]]
  }

  a <- (y + x)/2
  m <- (y - x)

  if(chk == 0) r <- list(a = a, m = m)
  if(chk == 1) r <- cbind(a = a, m = m)

  r
}
# =============================================================================.
#' Converts AM values to log2(counts)
#' @description
#' x = a + m/2, y = a - m/2
# -----------------------------------------------------------------------------.
#' @param \code{a}
#' @param \code{m}
#'
#' @return \code{list} with x and y values
# -----------------------------------------------------------------------------.
ma2lc <- function(a, m = NULL) {

  chk  <- 0

  if(is.null(m) & ! is.null(dim(a))) {
    m <- a[,2]
    a <- a[,1]
    chk <- 1
  }
  if(is.null(m) & is.list(a)) {
    m <- a[[2]]
    a <- a[[1]]
  }

  x <- a + m/2
  y <- a - m/2

  if(chk == 0) r <- list(x = x, y = y)
  if(chk == 1) r <- cbind(x = x, y = y)

  r
}
# =============================================================================.
#' Histogram
# -----------------------------------------------------------------------------.
#' @param x
#' @param y
#' @param breaks
#' @param x.name
#' @param y.name
#' @param x.col
#' @param y.col
#' @param ...
# -----------------------------------------------------------------------------.
#' @return produce an histogram in the active graphic device
# -----------------------------------------------------------------------------.
dualHist <- function(x, y, breaks=100, x.name=NULL, y.name=NULL, x.col = rgb(0,0,0,0.5), y.col = rgb(1,0,0,0.5), ...) {
  if(is.null(x.name)) x.name <- deparse(substitute(x))
  if(is.null(y.name)) y.name <- deparse(substitute(y))
  hx <- hist(x, breaks=breaks, plot=F)
  hy <- hist(y, breaks=breaks, plot=F)
  xlim <- range(c(hx$breaks, hy$breaks))
  ylim <- c(0, max(hy$counts, hx$counts, na.rm = T))
  hist(x, breaks=breaks, xlim=xlim, ylim=ylim, col=x.col, border=rgb(0,0,0,0), ...)
  hist(y, breaks=breaks, xlim=xlim, ylim=ylim, col=y.col, border=rgb(0,0,0,0), ..., add=T)
  legend("topright", c(x.name, y.name), fill=c(x.col, y.col), bty='n')
}
# =============================================================================.
#' Scatterplot with points colored according to kNN density estimation
# -----------------------------------------------------------------------------.
#' @param \code{X} matrix with two columns of numeric values
#' @param \code{col} index of color gradient (0 = grey, 1 = red, 2 = green, 3 = blue)
#' @param \code{grp} logical vector defining an highlighted subset of points
#' @param \code{grp.col} index of color gradient used for the grp subset
#' @param \code{zline} logical, activates the drawing of line y = 0 when TRUE
#' @param \code{xyline} logical, activates the drawing of line y = x when TRUE
#' @param \code{k} number of nearest neighbors for kNN density estimation
#' @param \code{...} optional parameters passed to the plot function
#'
#' @return produce a scatter plot in the active graphic device
# -----------------------------------------------------------------------------.
scPlotCounts <- function(X, col=0, grp=NULL, grp.col=0, zline=F, xyline=F, k=256, ...) {

  o <- 0:49/49
  s <- 0.5
  point.colors <- list(
    grey(99:0/132),
    c(
      rgb(o, 0, 0, s),
      rgb(1, o, 0, s)
    ),
    c(
      rgb(0, o, 0, s),
      rgb(o, 1, 0, s)
    ),
    c(
      rgb(0, 0, o, s),
      rgb(0, o, 1, s)
    ),
    rep(rgb(0,0.5,1,0.5), 100)
  )

  Y <- X[finiteValues(X),]
  dens <- densityByKNN(X, k=k)
  clrs <- point.colors[[col+1]][ceiling(dens$prc)]
  plot(X, col=clrs, ...)
  if(! is.null(grp)) {
    dens <- densityByKNN(X[grp,], k=k/2, sbj = Y)
    clrs <- point.colors[[grp.col+1]][ceiling(dens$prc)]
    points(X[grp,], col=clrs, pch=20, cex=0.3)
  }
  if(zline) {
    abline(h=0, col=rgb(0,0,0,0.5), lwd=1)
  }
  if(xyline) {
    abline(a=0, b=1, col=rgb(0,0,0,0.5), lwd=1)
  }
}
