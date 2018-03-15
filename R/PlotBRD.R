# =============================================================================.
#' Control plots from BRD results
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{BRD},
#'   \link{CDaDaDR}
# -----------------------------------------------------------------------------.
#' @param brd
#' result from the \link{BRD} function.
#'
#' @param with.axes
#' logical, include plot axes (default = T).
#'
#' @param with.legend
#' logical, include legends (default = T).
#'
#' @param res
#' number of bins used to make the intensity-density scatterplot
#' (default = 400).
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @export
PlotBRD <- function(brd, with.axes = T, with.legend = T, res = 400) {

  clr.prm <- AutoColorParameters("ry")
  clr.prm$thresholds <- round(clr.prm$thresholds, 2)

  ncl <- brd$parameters$ncl
  grp.clr <- rgb(1, 0.8, 0)
  if(ncl > 1) {
    clr <- TransformColors(SuperRainbow(3), delta.H = 30)[2:ncl]
    grp.clr <- c(grp.clr, clr)
  }

  bg <- 0
  if(! is.null(brd$bg_clurank)) {
    grp.clr <- grp.clr[brd$bg_clurank]
    bg <- brd$bg_cluster
  }

  idx   <- which(brd$subsets$core)
  grp   <- brd$subsets$cluster[idx]
  alpha <- rep(0.1, length(idx))

  spc <- 0
  if(with.legend) spc <- c(0, -0.15)
  lim <- with(
    brd$dred, xylim(projection[! rare, 1:2], spacing = spc, margin = 0.1)
  )

  # Plot CDaDaDR ///////////////////////////////////////////////////////////////
  # o <- order(brd$dred$knn_density)
  # with(
  #   brd$dred, plot_samples(
  #     projection[o, ], xlim = lim$x, ylim = lim$y,
  #     col = colorize(knn_density[o], mode = "rank", colors = "ry"),
  #     alpha = ! rare[o], main = "CDaDaDR"
  #   )
  # )
  # Plot corrected density /////////////////////////////////////////////////////
  o <- order(brd$dred$density)
  with(
    brd$dred, plot_samples(
      projection[o, ], xlim = lim$x, ylim = lim$y, axes = with.axes,
      col = colorize(density[o], mode = "rank", colors = "ry"),
      alpha = ! rare[o], main = "density"
    )
  )
  if(with.legend) {
    ColorLegend(
      "b", horiz = T, size = c(70, 3), parameters = clr.prm,
      ticks = 0:4/4, tick.pos = -1, cex = 0.75
    )
  }
  # Plot subsets ///////////////////////////////////////////////////////////////
  o <- order(brd$dred$density)
  with(
    brd$dred, plot_samples(
      projection[o, ], xlim = lim$x, ylim = lim$y, axes = with.axes,
      col = colorize(density[o], mode = "rank"), alpha = ! rare[o],
      main = "clusters"
    )
  )
  plot_samples(brd$subsets$x[idx, ], col = grp.clr[grp], alpha = alpha, add = T)
  if(is.null(brd$bg_cluster)) {
    legend("bottomleft", legend = brd$status, cex = 0.75, bty = 'n')
  } else {
    if(with.legend) {
      legend(
        "bottomleft", legend = "background candidates",
        fill = grp.clr[bg], cex = 0.75, bty = 'n'
      )
    }
  }
  # Plot intensity versus density //////////////////////////////////////////////
  cmf <- function(k) colorize(k, mode = "rank")
  h <- with(
    brd$dred, BivariateDensity(
      intensity, density, nx = res, method = "ash", ash = list(m = c(3, 3)),
      plot = T, axes = with.axes, # mapper = cmf,
      xlab = "background intensity", ylab = "density", main = "thresholds"
    )
  )
  # Overlay subsets
  ColorChannel(grp.clr, "a") <- 0.15
  with(
    brd$subsets, points(b[idx], d[idx], pch = 20, cex = 0.5, col = grp.clr[grp])
  )
  abline(h = min(brd$subsets$d), col = "red", lwd = 1)
  # Plot count distribution ////////////////////////////////////////////////////
  y <- brd$log2counts
  i <- brd$subsets$i
  h <- SideBySideDensity(
    y, method = "ash",
    xlab = "columns of the read count matrix", ylab = "log2(counts)",
    main = "distributions", las = 2, x.labels = F
  )
  y <- y[i, ]
  if(! is.null(brd$bg_cluster)) {
    k <- which(brd$subsets$cluster == brd$bg_cluster)
    for(j in 1:ncol(y)) {
      x <- rep(j, length(idx)) + runif(length(idx), -0.45, 0.45)
      points(x, y[idx, j], pch = 20, cex = 0.5, col = grp.clr[grp])
    }
  }
}
