# =============================================================================.
#' Signal Extraction Scaling
# -----------------------------------------------------------------------------.
#' @param x
#' Input or IgG raw read counts
#'
#' @param y
#' ChIP raw read counts
# -----------------------------------------------------------------------------.
#' @return list
# -----------------------------------------------------------------------------.
diaz.SES <- function(x, y) {

  i <- order(y)
  Y <- cumsum(y[i])
  X <- cumsum(x[i])

  n <- length(i)
  p <- Y/Y[n]
  q <- X/X[n]

  k <- which.max(abs(q - p))

  list(alpha = Y[k]/X[k], idx = (1:n)[i[1:k]], x.p = q, y.p = p)
}
# =============================================================================.
#' Scaling factors between read counts in different samples
# -----------------------------------------------------------------------------.
#' @param cnt
#' matrix of log transformed read counts
#' (columns = samples, rows = genomic intervals)
#'
#' @param rpkb
#' matrix of log transformed reads per kb associated to the \code{cnt} matrix
#' (columns = samples, rows = genomic intervals)
#'
#' @param background
#' logical vector defining the background subset of counts (default = all)
#'
#' @param enriched
#' logical vector defining the enriched subset of counts (default = none)
#'
#' @param replicates
#' logical matrix indicating replicate samples (default = none)
#'
#' @param center
#' function for center of difference calculations.
#' The default center function is \link{findMode} which identifies the mode in
#' a kernel density estimation of the distribution of read count differences
#' using the generic R function \code{density}.
#' This default function can be replaced by \code{median}, \code{mean}, or
#' similar functions.
#'
#' @param pairwise.density
#' logical (default = T, recommended).
#' If \code{TRUE} count density estimations and subsequent selections of
#' counts assumed to exhibit minimal variations are performed in bi-variate
#' spaces corresponding to all pairwise sample comparisons.
#' If \code{FALSE} count density estimation is performed only once, in the
#' multi-variate space of all samples.
#'
#' @param d
#' density threshold (in percentage) used to select a subset of counts assumed
#' to reflect minimal variations
#'
#' @param k
#' number of nearest neighbors for kNN density estimation
#'
#' @param plot
#' @param as.png
#' @param xylim
#' @param xylim.rpkb
#' @param ...
#' optional parameters passed to the \link{scalingFactors} function.
#' These include:
# -----------------------------------------------------------------------------.
#' @return
#' scalingFactors returns a \code{list} with the following elements:
#' \item{pairwise.factors}{
#'   matrix of pairwise normalization factors between samples
#'   (\code{cnt} and \code{rpkb} columns)
#' }
#' \item{low.variation}{
#'   logical vector indicating the subset of rows assumed invariable by the
#'   normalization
#' }
# -----------------------------------------------------------------------------.
scalingFactors <- function(cnt, rpkb, background=NULL, enriched=NULL, replicates=0, center=NULL, pairwise.density=T, d=80, k=128, plot=F, as.png=T, xylim=NULL, xylim.rpkb=NULL, ...) {

  control.plots <- function(i, j) {
    xlab <- i
    ylab <- j
    if(! is.null(colnames(cnt))) {
      xlab <- colnames(cnt)[i]
      ylab <- colnames(cnt)[j]
    }
    if(as.png) {
      png(paste("scalingFactors.", xlab, "_x_", ylab, ".png", sep=""), width=1500, height=550)
      layout(matrix(1:3,1,3,byrow=T))
    }
    # Density plot
    plot(density(rpkb[chk,j] - rpkb[chk,i]), main="", xlab="delta",lwd=2)
    dns <- density(rpkb[idx,j] - rpkb[idx,i])
    dns$y <- dns$y * sum(chk & invariable)/sum(chk)
    points(dns, type='l', col=rgb(1,0,1,0.7),lwd=2)
    abline(v=pwf[i,j], col=rgb(0.0,0.5,1),lwd=2)
    # RPKB scatterplot
    scPlotCounts(
      rpkb[,c(i,j)], k=128, xyline=T, xlim=xylim.rpkb, ylim=xylim.rpkb,
      pch=20, cex=0.3, xlab=xlab, ylab=ylab, ...
    )
    abline(h=0, col=rgb(0,0,0,0.5))
    abline(v=0, col=rgb(0,0,0,0.5))
    if(plot.en) {
      points(rpkb[chk & en, c(i,j)], pch=20, cex=0.3, col=rgb(1,0,0,0.2))
    }
    if(plot.bg) {
      points(rpkb[chk & bg, c(i,j)], pch=20, cex=0.3, col=rgb(0,0,1,0.2))
    }
    legend(
      "bottomright", paste(c("enriched", "background"), c(sum(en), sum(bg))),
      fill=c(rgb(1,0,0,0.7), rgb(0,0,1,0.7)), bty='n'
    )
    if(replicates[i,j]) {
      legend("topleft", "replicates", bty='n')
    } else {
      points(rpkb[idx,c(i,j)], pch=20, cex=0.3, col=rgb(0,1,0,0.1))
    }
    abline(a=pwf[i,j], b=1, col=rgb(0.0,0.5,1),lwd=2)
    # Read count scatterplot
    scPlotCounts(
      cnt[,c(i,j)], k=128, xyline=T, xlim=xylim, ylim=xylim,
      pch=20, cex=0.3, xlab=xlab, ylab=ylab, ...
    )
    if(plot.en) {
      points(cnt[chk & en, c(i,j)], pch=20, cex=0.3, col=rgb(1,0,0,0.2))
    }
    if(plot.bg) {
      points(cnt[chk & bg, c(i,j)], pch=20, cex=0.3, col=rgb(0,0,1,0.2))
    }
    legend(
      "bottomright", paste(c("enriched", "background"), c(sum(en), sum(bg))),
      fill=c(rgb(1,0,0,0.7), rgb(0,0,1,0.7)), bty='n'
    )
    if(replicates[i,j]) {
      legend("topleft", "replicates", bty='n')
    } else {
      points(cnt[idx,c(i,j)], pch=20, cex=0.3, col=rgb(0,1,0,0.1))
    }
    abline(a=pwf[i,j], b=1, col=rgb(0.0,0.5,1),lwd=2)
    if(as.png) {
      dev.off()
    }
  }

  n <- ncol(cnt)
  m <- nrow(cnt)

  if(is.null(center))     center <- findMode
  if(is.null(xylim))      xylim <- c(0, round(max(cnt, na.rm = T), 1))
  if(sum(replicates)==0)  replicates <- matrix(F, n, n)

  # Define counts used for density estimations (eligible as invariable)
  if(is.null(enriched))   enriched   <- matrix(F, m, n)
  if(is.null(background)) background <- matrix(T, m, n)
  if(is.null(ncol(enriched)))   enriched   <- matrix(enriched, m, n)
  if(is.null(ncol(background))) background <- matrix(background, m, n)

  plot.en <- any(enriched)
  plot.bg <- ! all(background)

  low.variation <- rep(T, m)

  chk <- finiteValues(cnt)
  # Estimate invariable counts based on global density
  if(! pairwise.density) {
    en <- apply(enriched, 1, any)
    bg <- apply(background, 1, all)
    invariable <- rep(F, m)
    idx <- which(chk & (bg & ! en))
    dns <- densityByKNN(rpkb[idx,], k = k)
    invariable[idx] <- dns$prc > d
    low.variation <- low.variation & invariable
  }

  if(is.null(xylim.rpkb)) xylim.rpkb <- round(range(rpkb[chk,], na.rm = T), 1)

  pwf <- matrix(0, n, n)

  pb <- txtProgressBar(min = 1, max = n, char = "|", style = 3)
  for(i in 1:n) {
    for(j in 1:n) {
      if(i!=j) {
        chk <- finiteValues(cnt[,c(i,j)])
        # Estimate invariable counts based on pairwise density
        if(pairwise.density) {
          en <- apply(enriched[,c(i,j)], 1, any)
          bg <- apply(background[,c(i,j)], 1, all)
          invariable <- rep(F, m)
          idx <- which(chk & (bg & ! en))
          dns <- densityByKNN(rpkb[idx,c(i,j)], k = k)
          invariable[idx] <- dns$prc > d
          if(! replicates[i,j]) low.variation <- low.variation & invariable
        }
        # Center of differences based on subset of rows assumed to be invariable
        idx <- which(chk & invariable)
        # In case of replicates, all rows are assumed to be invariable
        if(replicates[i,j]) idx <- which(chk)
        # Compute center of differences (rpkb is better than cnt)
        pwf[i,j] <- center(rpkb[idx,j] - rpkb[idx,i], na.rm=T)
        # Make optional control scatterplots (avoid redundant plots)
        if(plot & j > i) {
          control.plots(i, j)
        }
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  list(pairwise.factors = pwf, low.variation = low.variation)
}
# =============================================================================.
#' Scaling factors between read counts in different samples
# -----------------------------------------------------------------------------.
#' @description
#' Compute normalization factors based on control regions (for instance genes),
#' possibly using pre-defined background or enriched regions to further control
#' the set of genomic intervals that are assumed to exhibit minimal variations.
# -----------------------------------------------------------------------------.
#' @param aln
#' list of \link{GAlignments} objects
#'
#' @param ctr
#' \link{GRanges} defining control regions
#'
#' @param bgr
#' \link{GRanges} defining background regions (default = all control regions)
#'
#' @param enr
#' \link{GRanges} defining enriched regions to be avoided (default = none)
#'
#' @param use.bgr
#' logical vector defining the samples from which background regions are taken into account (default = all)
#'
#' @param use.enr
#' logical vector defining the samples from which enriched regions are taken into account (default = all)
#'
#' @param replicates
#' logical matrix indicating replicate samples (default = none)
#'
#' @param lbls
#' character vector of samples names (default = none)
#'
#' @param dither
#' logical, activates count dithering (default = T, recommended)
#'
#' @param plot
#' logical, activates the drawing of control plots (default = F)
#'
#' @param ...
#' optional parameters passed to the \link{scalingFactors} function.
#' These include:
#' \itemize{
#'
#'   \item center
#'   function (see \link{scalingFactors})
#'
#'   \item pairwise.density
#'   logical controling the use of pairwise or global count density estimations
#'   and invariable count selections
#'
#'   \item d
#'   density threshold
#'
#'   \item k
#'   number of nearest neighbors in kNN density estimation
#' }
#' as well as other parameters controlling the optional scatterplots.
# -----------------------------------------------------------------------------.
#' @return
#' makeScalingFactors returns a \code{list} with the following elements:
#' \item{counts}{
#'   matrix of read counts in the control regions
#'   (columns = samples, rows = genomic intervals)
#' }
#' \item{scaling}{
#'   vector of normalization factor for read counts in each sample
#'   (\code{counts} columns)
#' }
#' \item{pairwise}{
#'   matrix of pairwise normalization factors between samples
#'   (\code{counts} columns)
#' }
#' \item{low.variation}{
#'   logical vector indicating the subset of rows assumed invariable by the
#'   normalization
#' }
# -----------------------------------------------------------------------------.
makeScalingFactors <- function(aln, ctr, bgr=NULL, enr=NULL, use.bgr=NULL, use.enr=NULL, replicates=0, lbls=NULL, dither=T, plot=F, ...) {

  # Detect sample names
  if(is.null(lbls)) {
    if(! is.null(colnames(replicates))) lbls <- colnames(replicates)
    if(! is.null(names(aln)))           lbls <- names(aln)
  }

  bg <- NULL
  en <- NULL

  # Find control regions contained in pre-defined background regions
  if( is.logical(bgr)) {
    bg <- bgr
  }
  if( is(bgr, "GRanges")) {
    bg <- countOverlaps(ctr, bgr, type = "within", ignore.strand=T)
    bg <- bg > 0
  }
  if(is.list(bgr)) {
    bg <- matrix(T, length(ctr), length(bgr))
    for(i in 1:length(bgr)) {
      ovl <- countOverlaps(ctr, bgr[[i]], type = "within", ignore.strand=T)
      bg[,i] <- ovl > 0
    }
    if(! is.null(use.bgr)) bg[,! use.bgr] <- T
  }
  # Find control regions that overlap with pre-defined enriched regions
  if( is.logical(enr)) {
    en <- enr
  }
  if( is(enr, "GRanges")) {
    en <- countOverlaps(ctr, enr, type = "any", ignore.strand=T)
    en <- en > 0
  }
  if(is.list(enr)) {
    en <- matrix(F, length(ctr), length(enr))
    for(i in 1:length(enr)) {
      ovl <- countOverlaps(ctr, enr[[i]], type = "any", ignore.strand=T)
      en[,i] <- ovl > 0
    }
    if(! is.null(use.enr)) en[,! use.enr] <- F
  }

  # Make read counts in control regions
  message("Counting reads in control regions")
  cnt <- makeReadCountMatrix(
    aln, ctr, type = "any", ignore.strand=T
  )
  if(! is.null(lbls)) colnames(cnt) <- lbls
  rpkb <- 1E3 * cnt / width(ctr)

  # Apply dithering to read counts to make more robust density estimations
  x <- cnt
  if(dither) {
    x <- ditherCounts(x)
  }
  y <- 1E3 * x / width(ctr)

  # Log transform read counts
  x <- log10(x)
  y <- log10(y)

  # Optionally: make control graph of the dithering procedure
  if(dither & plot) {
    dualHist(
      log10(cnt), x, breaks=100, x.name="raw", y.name="dithered",
      xlab="log10(counts)", ylab="frequency", main="read counts"
    )
    dualHist(
      log10(rpkb), y, breaks=100, x.name="raw", y.name="dithered",
      xlab="log10(rpkb)", ylab="frequency", main="reads per kb"
    )
    plot(0, type='n', axes=F, xlab="", ylab="") # Empty plot
  }

  # Compute normalization factors
  res <- scalingFactors(
    cnt = x, rpkb = y,
    background = bg, enriched = en, replicates = replicates, plot = plot, ...
  )

  pwnf         <- res$pairwise.factors
  norm.factors <- apply(pwnf, MARGIN = 2, FUN = mean, na.rm=T)

  # Revert log transformation
  norm.factors <- 10^norm.factors
  pwnf         <- 10^pwnf

  list(
    counts        = cnt,
    background    = bg,
    enriched      = en,
    scaling       = norm.factors,
    pairwise      = pwnf,
    low.variation = res$low.variation
  )
}
# =============================================================================.
#' Center of differences between log transformed counts in different samples
# -----------------------------------------------------------------------------.
#' @param x
#' matrix of log transformed read counts (columns = samples, rows = genomic intervals)
#'
#' @param invariable
#' logical vector defining the set of rows in the read count matrix x that are
#' assumed invariable (default = all).
#'
#' @param center
#' function for center of difference calculations.
#' The default center function is \link{findMode} which identifies the mode in
#' a kernel density estimation of the distribution of read count differences
#' using the generic R function \code{density}.
#' This default function can be replaced by \code{median}, \code{mean}, or
#' similar functions.
#'
#' @param replicates
#' logical matrix indicating replicate columns (default = none)
#'
#' @param plot
#' logical, activates the drawing of control scatterplots (default = F)
#'
#' @param xylim
#' defining ranges of x and y axes for scatterplots (default = auto)
#'
#' @param ...
#' graphical parameters passed to the \link{scPlotCounts} function.
# -----------------------------------------------------------------------------.
#' @return
#' matrix of pairwise center of differences between all samples
# -----------------------------------------------------------------------------.
centerOfDifferences <- function(x, invariable = NULL, center = NULL, replicates=0, plot=F, xylim=NULL, ...) {

  n <- ncol(x)
  m <- nrow(x)

  if(is.null(invariable)) invariable <- rep(T, m)
  if(is.null(center))     center <- findMode
  if(is.null(xylim))      xylim <- c(0, round(max(x, na.rm = T), 1))
  if(sum(replicates)==0)  replicates <- matrix(F, n, n)

  f <- matrix(0, n, n)

  pb <- txtProgressBar(min = 1, max = n, char = "|", style = 3)
  for(i in 1:n) {
    for(j in 1:n) {
      chk <- is.finite(x[,j]) & is.finite(x[,i])
      # Center of differences based on subset of rows assumed to be invariable
      idx <- which(chk & invariable)
      # In case of replicate samples, all rows are assumed to be invariable
      if(replicates[i,j]) idx <- which(chk)
      if(i!=j) {
        # Compute center of differences
        f[i,j] <- center(x[idx,j] - x[idx,i], na.rm=T)
        # Make optional control scatterplots
        if(plot) {
          xlab <- i
          ylab <- j
          if(! is.null(colnames(x))) {
            xlab <- colnames(x)[i]
            ylab <- colnames(x)[j]
          }
          scPlotCounts(
            x[,c(i,j)], k=128, xyline=T, xlim=xylim, ylim=xylim,
            pch=20, cex=0.3, xlab=xlab, ylab=ylab, ...
          )
          if(replicates[i,j]) {
            legend("topleft", "replicates", bty='n')
          } else {
            points(x[idx,c(i,j)], pch=20, cex=0.3, col=rgb(1,0,1,0.3))
          }
          abline(a=f[i,j], b=1, col=rgb(0.0,0.5,1),lwd=2)
        }
      }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  f
}

# =============================================================================.
#' Normalization based on center of differences
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{normalizeByPrincipalComponent}
# -----------------------------------------------------------------------------.
#' @param M
#' matrix of log transformed read counts (columns = samples, rows = genomic intervals)
#'
#' @param invariable
#' logical vector defining the set of rows in the read count matrix M that are
#' assumed invariable (default = all).
#'
#' @param method
#' either "mode" (default), "median" or "mean"
#'
#' @note
#' This is a deprecated function that should be removed in future versions
# -----------------------------------------------------------------------------.
#' @return
#' normalizeByCenterOfDifferences returns a \code{list} with the following elements:
#' \item{X}{matrix of normalized read counts}
#' \item{f}{matrix of pairwise normalization factors}
# -----------------------------------------------------------------------------.
normalizeByCenterOfDifferences <- function(M, invariable = NULL, method=c("mode", "median", "mean"), cod.matrix=NULL) {
  if(is.null(cod.matrix)) {
    f <- centerOfDifferences(M=M, invariable=invariable, method=method)
  } else {
    f <- cod.matrix
  }
  f <- apply(f, MARGIN = 2, FUN = mean)
  list(X=t(t(M) - f), f=f)
}
# =============================================================================.
#' Estimate normalization parameters based on principal component analysis
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{normalizeByPrincipalComponent}
# -----------------------------------------------------------------------------.
#' @param M
#' matrix of log transformed read counts (columns = samples, rows = genomic intervals)
#'
#' @param invariable
#' logical vector defining the set of rows in the read count matrix M that are
#' assumed invariable (default = all).
#'
#' @param pc
#' character defining which principal component should be extracted
#' (default = "PC1").
# -----------------------------------------------------------------------------.
#' @return
#' getPrincipalComponent returns a \code{list} with the following elements:
#' \item{SLOPES}{matrix of b values according to the linear model y = b x + a}
#' \item{INTRCP}{matrix of a values according to the linear model y = b x + a}
# -----------------------------------------------------------------------------.
getPrincipalComponent <- function(M, invariable=NULL, pc="PC1") {
  n <- ncol(M)
  m <- nrow(M)
  if(is.null(invariable)) {
    invariable <- rep(T, m)
  }
  # Use PCA to extract the first principal component (axis of maximal variance)
  r <- prcomp(M[which(finiteValues(M) & invariable),], retx = F)
  # Derive slope and intercept for each variable as a function of each other
  SLOPES <- matrix(NA, n, n)
  INTRCP <- matrix(NA, n, n)
  for(i in 1:n) {
    for(j in 1:n) {
      SLOPES[i,j] <- r$rotation[j,pc]/r$rotation[i,pc]
      INTRCP[i,j] <- r$center[j] - SLOPES[i,j] * r$center[i]
    }
  }
  list(SLOPES=SLOPES, INTRCP=INTRCP, PCA=r)
}
# =============================================================================.
#' PCA-based normalization
# -----------------------------------------------------------------------------.
#' @author Benjamin Leblanc
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{getPrincipalComponent}
# -----------------------------------------------------------------------------.
#' @param M
#' matrix of log transformed read counts (columns = samples, rows = genomic intervals)
#'
#' @param pc
#' normalization parameters resulting from function \link{getPrincipalComponent}
#'
# -----------------------------------------------------------------------------.
#' @return
#' matrix of normalized read counts
# -----------------------------------------------------------------------------.
normalizeByPrincipalComponent <- function(M, pc) {
  n <- ncol(M)
  m <- nrow(M)
  M.nrm <- matrix(0, m, n)
  pb <- txtProgressBar(min = 1, max = n, char = "|", style = 3)
  for(i in 1:n) {
    x <- matrix(M[,i], n, m, byrow=T)
    x <- apply((x - pc$INTRCP[,i])/pc$SLOPES[,i], MARGIN=2, FUN=mean)
    M.nrm[,i] <- x
    setTxtProgressBar(pb, i)
  }
  close(pb)
  colnames(M.nrm) <- colnames(M)
  M.nrm
}
