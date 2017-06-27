# =============================================================================.
#' PlotBRD
# -----------------------------------------------------------------------------.
#' @param brd result from \link{BRD}
#' @param core.only logical (default = F)
#' @param with.legend logical (default = T)
#' @return NULL
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
PlotBRD <- function(brd, core.only = F, with.legend = T) {

  clr.prm <- AutoColorParameters("ry")
  clr.prm$thresholds <- round(clr.prm$thresholds, 2)

  ncl <- brd$parameters$ncl
  grp.clr <- rgb(1, 0.9, 0)
  if(ncl > 1) {
    grp.clr <- c(
      grp.clr, transformColors(SuperRainbow(3), delta.H = 30)[2:ncl]
    )
  }
  bg <- 0
  if(! is.null(brd$bg_clurank)) {
    grp.clr <- grp.clr[brd$bg_clurank]
    bg <- brd$bg_cluster
  }
  if(core.only) {
    idx   <- which(brd$subsets$core)
    grp   <- brd$subsets$cluster[idx]
    alpha <- rep(0.1, length(idx))
  } else {
    grp.clr <- rbind(
      transformColors(grp.clr, S.range = 0.2),
      transformColors(grp.clr, S.range = 1.0)
    )
    grp.clr <- as.vector(grp.clr)
    idx   <- 1:nrow(brd$subsets$x)
    grp   <- with(brd$subsets, 2 * cluster + core - 1)
    alpha <- ifelse(brd$subsets$core, 0.2, 0.05)
    bg    <- 2 * bg
  }

  spc <- 0
  if(with.legend) spc <- c(0, -0.15)
  lim <- with(
    brd$dred, xylim(projection[! rare, 1:2], spacing = spc, margin = 0.1)
  )

  o <- order(brd$dred$density)

  # Plot Count Density after Dithering and Dimensionality Reduction
  with(
    brd$dred, plot_samples(
      projection[o, ], xlim = lim$x, ylim = lim$y,
      col = colorize(density[o], mode = "rank", colors = "ry"),
      alpha = ! rare[o], main = "CDaDaDR"
    )
  )
  if(with.legend) {
    colorLegend(
      "b", horiz = T, size = c(70, 3), parameters = clr.prm,
      ticks = 0:4/4, tick.pos = -1, cex = 0.75
    )
  }
  # Plot CDaDaDR and overlay subsets
  with(
    brd$dred, plot_samples(
      projection[o, ], xlim = lim$x, ylim = lim$y,
      col = colorize(density[o], mode = "rank"), alpha = ! rare[o],
      main = "segmentation and clustering"
    )
  )
  plot_samples(brd$subsets$x[idx, ], col = grp.clr[grp], alpha = alpha, add = T)
  if(with.legend ) {
    if(is.null(brd$bg_cluster)) {
      legend("bottomleft", legend = brd$status, cex = 0.75, bty = 'n')
    } else {
      legend(
        "bottomleft", legend = "background candidates",
        fill = grp.clr[bg], cex = 0.75, bty = 'n'
      )
    }
  }
  # Controls versus experiments

  k <- with(brd$parameters, match(controls, experiments))
  x <- with(brd, rowMeans(as.matrix(log2counts[, - k])))
  y <- with(brd, rowMeans(as.matrix(log2counts[,   k])))
  am <- cbind(
    `average log2(counts)` = (x + y) / 2,
    `average log2ratio`    = (y - x)
  )

  with(
    brd$dred, plot_samples(
      am[o, ], symetric = F,
      col = colorize(density[o], mode = "rank", colors = "ry"),
      alpha = ! rare[o], main = "controls versus experiments"
    )
  )

  plot(density(am[, 2]))
  # plot_samples(
  #   am[brd$subsets$i, ][idx, ], col = grp.clr[grp], alpha = alpha, add = T
  # )
}
