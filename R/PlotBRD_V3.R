# =============================================================================.
#' Plot empirical distributions of read counts
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{DitherCounts},
#'   \link{ReadCountMatrix}
# -----------------------------------------------------------------------------.
#' @param x
#' a read count matrix, preferably dithered and log2 transformed.
#'
#' @param ...
#' optional arguments forwarded to the \link{SideBySide} function.
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @examples
#' # Simulation of a read count matrix with 3 populations of 10000 observations
#' p <- c(10000, 10000, 10000)
#' m <- DefineSimulation(
#'   chip = 5, patterns = c("^", "v"), enrichment = c(1.0, 3), replicate = 2
#' )
#' r <-  MakeSimulation(p = p, m = m)
#'
#' grp <- r$group # Population memberships
#' cnt <- r$data  # Simulated counts
#'
#' # Prepare figure layout and graphic options
#' layout(matrix(1:4, 2, 2, byrow = TRUE))
#' par(pch = 20)
#'
#' # Show the empirical distribution of simulated populations
#' l2c <- log2(DitherCounts(cnt)) # Dithering and log2 transformation
#' xyl <- range(l2c[FiniteValues(l2c), ])
#'
#' r <- PlotCountDistributions(l2c, ylim = xyl, main = "Total")
#' for(i in sort(unique(grp))) {
#'   main <- ifelse(i == 1, "Invariable subset", paste("Variable subset", i - 1))
#'   r <- PlotCountDistributions(l2c[grp == i, ], ylim = xyl, main = main)
#' }
# -----------------------------------------------------------------------------.
#' @export
PlotCountDistributions <- function(
  X, pops = NULL, sampling = NULL, vb = 0, names = T,
  palette = "magma", gradient = "hcl.mono.grey", # "hcl.duo.light"
  legend = T, inset = c(-0.25, 0), saturation = 0.9, ...
) {
  if(is.null(pops)) {
    par(mar = c( 5 + 11 * names, 5, 4, 2))
    SideBySide(
      X, vb = vb, sampling = sampling,
      colors = "Wry", gradient = "colorize", saturation = saturation,
      las = 2, label = "log2(counts)", names = names, ...
    )
  } else{
    par(mar = c( 5 + 11 * names, 5, 4, 2 + 6 * legend))
    clr <- list(
      d = "grey",
      p = c("grey", shades::gradient(palette, steps = nlevels(pops)))
    )
    SideBySide(
      X, pops = pops, db = nlevels(pops) * 10, vb = vb,
      sampling = sampling, smoothing = c(3, 5),
      proportions = "static", ordering = "static", scales = "maximized",
      colors = clr["p"], saturation = saturation, gradient = gradient,
      las = 2, label = "log2(counts)", names = names, ...
    )
    if(legend) {
      mpr <- lapply(
        clr$p, ColorMapper, gradient = gradient, saturation = saturation
      )
      mpr.clr <- sapply(mpr, function(x) x$cmf(0.8))
      k <- which(summary(pops) > 0)
      legend(
        "topright", inset = inset,
        legend = levels(pops)[k], fill = mpr.clr[k], bty = "n", xpd = T
      )
    }
  }
}

# =============================================================================.
#' Control plots from BRD results
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{BRD},
#'   \link{CDaDaDR.2D}
# -----------------------------------------------------------------------------.
#' @description
#' The \code{PlotBRD} function produces control graphs representing intermediate
#' results of the \link{BRD} function. These control graphs are useful to adapt
#' parameters of the BRD method to the considered read count matrix.
#'
#' @param brd
#' result from calling the \link{BRD} function.
#'
#' @param title
#' common title for the generated plots (default = none).
#'
#' @param plots
#' character vector specifying which plots are being generated. By default
#' \code{PlotBRD} generates 4 plots named "density", "clusters", "thresholds",
#' and "distributions". It is possible to generate only a subset of these plots
#' by using the corresponding vector of plot names, like for instance
#' \code{plots = c("density", "thresholds")}.
#'
#' @param with.axes
#' logical, include plot axes or not (default = TRUE, yes).
#'
#' @param with.legend
#' logical, include legends or not (default = TRUE, yes).
#'
#' @param res
#' resolution in number of bins for the "thresholds" plot (default = 400).
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @export
PlotBRD <- function(
  brd, title = "",
  plots = c("density", "clusters", "thresholds", "distributions"),
  palette = "magma", gradient = "hcl.duo.light", saturation = 1.0,
  with.axes = TRUE, with.legend = TRUE, res = 400
) {

  plots <- tolower(plots)

  ncl <- brd$parameters$ncl
  idx   <- which(brd$subsets$core)
  grp   <- brd$subsets$cluster[idx]

  sbs <- rep(0, nrow(brd$dred$projection))
  sbs[brd$subsets$i[idx]] <- grp
  sbs <- sbs + 1

  bgc <- brd$bg_cluster
  clr <- hsv(seq(180, 300, length.out = ncl) / 360)
  clr[bgc] <- rgb(1.0, 1.0, 0.0)
  clr <- list(d = "grey", p = c("grey", clr))
  mpr <- lapply(
    clr$p, ColorMapper, gradient = gradient, saturation = saturation
  )
  mpr.clr <- sapply(mpr, function(x) x$cmf(0.95))

  # Observations with lowest density
  low <- log10(nrow(brd$dred$projection)) - 1
  low <- max(1, min(4, low))
  low <- round(approx(1:4, y = c(10, 50, 250, 2500), xout = low)$y)
  low <- low / nrow(brd$dred$projection)
  low <- RankScore(brd$dred$density) < low

  # ----------------------------------------------------------------------------
  mklim <- function(spc) {
    xylim(
      brd$dred$projection[! low, 1:2], spacing = with.legend * spc, margin = 0.1
    )
  }
  # ----------------------------------------------------------------------------
  lgsts <- function(pos) {
    if(with.legend) {
      if(! is.null(bgc)) {
        legend(
          pos, legend = "background candidates",
          fill = mpr.clr[bgc + 1], cex = 0.75, bty = 'n'
        )
      } else {
        legend(pos, legend = brd$status, cex = 0.75, bty = 'n')
      }
    }
  }
  # ----------------------------------------------------------------------------
  lgclr <- function(pos, clr, gradient = gradient, saturation = saturation, ...) {
    if(with.legend) {
      mpr <- ColorMapper(clr, gradient = gradient, saturation = saturation)
      ColorLegend(pos, parameters = mpr$cmp, ticks = 0:4/4, cex = 0.75, ...)
    }
  }
  # ----------------------------------------------------------------------------
  lim <- mklim(c(-0.2, -0.2))
  # Plot density ///////////////////////////////////////////////////////////////
  if("density" %in% plots) {
    with(
      brd,
      ScatterMaps(
        dred$projection, rng = cbind(lim$x, lim$y), x = "C1", y = "C2",
        bins = parameters$bins, smoothing = parameters$smoothing,
        colors = "Wry", gradient = "colorize", # "red", "hsv.mono.grey"
        scales = "absolute", render = "maximum",
        main = paste("density", title, sep = "\n")
      )
    )
    # points(brd$dred$projection[low, ], col = grey(0, 0.1), pch = 20, cex = 0.5)
    lgclr(
      clr = "Wry", gradient = "colorize", saturation = saturation,
      pos = "b", horiz = TRUE, size = c(70, 3), tick.pos = -1
    )
  }
  # Plot subsets ///////////////////////////////////////////////////////////////
  if("clusters" %in% plots) {
    with(
      brd,
      ScatterMaps(
        dred$projection, rng = cbind(lim$x, lim$y), x = "C1", y = "C2",
        bins = parameters$bins, smoothing = parameters$smoothing,
        pops = sbs, colors = clr["p"],
        gradient = gradient, scales = "absolute", render = "maximum",
        main = paste("density clusters", title, sep = "\n")
      )
    )
    # points(brd$dred$projection[low, ], col = grey(0, 0.1), pch = 20, cex = 0.5)
    lgclr(
      clr = "grey", gradient = gradient, saturation = saturation,
      pos = "b", horiz = TRUE, size = c(70, 3), tick.pos = -1
    )
    lgsts("topleft")
  }
  # Plot intensity versus density //////////////////////////////////////////////
  if("thresholds" %in% plots) {
    with(
      brd$dred,
      ScatterMaps(
        cbind(background = intensity, density = density),
        rng = cbind(range(intensity), range(density)), extend = c(1, 1.15),
        bins = brd$parameters$bins, smoothing = brd$parameters$smoothing,
        pops = sbs, colors = clr["p"],
        gradient = gradient, scales = "balanced", render = "maximum",
        main = paste("density thresholds", title, sep = "\n")
      )
    )
    # Overlay subsets
    # with(
    #   brd$subsets, points(
    #     b[idx], d[idx], pch = 20, cex = 0.5,
    #     col = shades::opacity(mpr.clr, 0.1)[grp + 1]
    #   )
    # )
    abline(h = min(brd$subsets$d), col = rgb(1, 0.0, 0), lwd = 1)
    lgclr(
      clr = "grey", gradient = gradient, saturation = saturation,
      pos = "l", horiz = FALSE, size = c(50, 3)
    )
    lgsts("topleft")
  }
  # Plot count distribution ////////////////////////////////////////////////////
  if("distributions" %in% plots) {
    ldc <- with(
      brd, RemoveInputBias(
        log2counts, controls = parameters$controls, as.log2 = TRUE, safe = TRUE
      )
    )
    ext <- c(1.0, 1 + 0.15 * with.legend)
    SideBySide(
      ldc, extend = ext, db = 10, vb = 0, sampling = 5E5,
      pops = sbs, scales = "balanced", render = "maximum",
      colors = clr["p"], gradient = gradient,
      label = "log2(counts)", names = FALSE, grid = "white",
      main = paste("count distributions", title, sep = "\n")
    )
    lgsts("topleft")
  }
}
