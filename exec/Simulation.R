library(Tightrope)

# FUNCTIONS ####################################################################

# =============================================================================.
#
# -----------------------------------------------------------------------------.
Preview <- function(brd.args, sim.prm) {
  r <- with(sim.prm, MakeSimulation(p = p, m = m))

  grp <- r$group
  cnt <- r$data

  chip <- grep("ChIP", colnames(cnt), value = T)
  ctrl <- grep("Input", colnames(cnt), value = T)

  l2c <- log2(DitherCounts(cnt))
  xyl <- range(l2c[FiniteValues(l2c), ])

  brd <- do.call(BRD, args = c(list(cnt = cnt, controls = ctrl), brd.args))

  layout(matrix(1:9, 3, 3, byrow = T))
  r <- PlotCountDistributions(l2c, ylim = xyl, main = "Total")
  for(i in sort(unique(grp))) {
    main <- ifelse(i == 1, "Invariable subset", paste("Variable subset", i - 1))
    r <- PlotCountDistributions(l2c[grp == i, ], ylim = xyl, main = main)
  }
  input <- rowMeans(l2c[, ctrl])
  for(lbl in chip) {
    r <- BivariateDensity(
      input, l2c[, lbl], xlim = xyl, ylim = xyl, parameters = list(color = "Wry"),
      method = "ash", ash = list(m = c(3, 3)), xlab = "average of Inputs",
      ylab = lbl
    )
    # points(cbind(input, l2c[, lbl])[brd$bg_members, ], pch = 20, cex = 0.5, col = rgb(0, 0.5, 1))
    for(i in -10:10) abline(a = i, b = 1, col = grey(0.5, alpha = 0.15))
    abline(a = 0, b = 1, col = grey(0.5))
  }

  layout(matrix(1:4, 2, 2, byrow = T))
  par(pch = 20)
  PlotBRD(brd, with.legend = F)

  reassign(chip, pos = globalenv(), src = environment())
  reassign(ctrl, pos = globalenv(), src = environment())
  reassign(cnt, pos = globalenv(), src = environment())
  reassign(grp, pos = globalenv(), src = environment())
  reassign(l2c, pos = globalenv(), src = environment())
  reassign(xyl, pos = globalenv(), src = environment())
  reassign(brd, pos = globalenv(), src = environment())
}


# =============================================================================.
#
# -----------------------------------------------------------------------------.
RunBenchmark <- function(brd.args, sim.prm, reproduce) {

  res <- list(
    passed = c(),
    failed = rep(F, reproduce)
  )
  for(i in 1:reproduce) {
    r <- with(sim.prm, MakeSimulation(p = p, m = m))

    grp <- r$group
    cnt <- r$data
    ctrl <- grep("Input", colnames(cnt), value = T)

    brd <- do.call(BRD, args = c(list(cnt = cnt, controls = ctrl), brd.args))

    if(grepl("successful", brd$status)) {
      res$passed <- rbind(res$passed, ScalingFactors(brd))
    } else {
      res$failed[i] <- T
    }
  }

  res
}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
RunTest <- function(test, values, brd.args, sim.prm, reproduce) {

  chk <- is.null(dim(values))
  if(chk) {
    n <- length(values)
  } else {
    n <- nrow(values)
  }

  r <- vector("list", n)
  for(i in 1:n) {
    message("[ benchmark ] ", test, " | ", i, " / ", n)
    if(chk) {
      brd.args[[test]] <- values[i]
    } else {
      brd.args[[test]] <- values[i, ]
    }
    r[[i]] <- RunBenchmark(brd.args, sim.prm, reproduce)
  }

  r
}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
PlotBenchmark <- function(r, test, values) {
  f <- sapply(lapply(r, "[[", "failed"), sum)
  names(f) <- values
  barplot(100 * (n - f) / n, xlab = test, ylab = "BRD success (%)")

  p <- lapply(lapply(r, "[[", "passed"), as.vector)
  BoxPlot(p, xlab = test, ylab = "scaling factors")
  y <- rep(quantile(unlist(p), probs = 0.005), n)
  text(x = 1:n, y = y, labels = values, pos = 1, srt = 45, xpd = T)
}

# CONFIGURATION ################################################################

brd.args <- list(ncl = 3, dither = 3, knn = 300, bdt = c(0.2, 0.05))
sim.prm <- list(
  p = c(2000, 10000, 10000),
  m = DefineSimulation(
    chip = 5, patterns = c("^", "v"), enrichment = c(1.0, 4.0), replicate = 2
  )
)

brd.args <- list(ncl = 3, dither = 3, knn = 300, bdt = c(0.2, 0.05))
sim.prm <- list(
  p = c(2000, 10000, 10000),
  m = DefineSimulation(
    chip = 7, patterns = c("^", "v"), enrichment = c(1.0, 3.0), replicate = 1
  )
)

brd.args <- list(ncl = 3, dither = 3, knn = 300, bdt = c(0.2, 0.05))
sim.prm <- list(
  p = c(2000, 10000, 10000),
  m = DefineSimulation(
    chip = 5, patterns = c("+", "-"), enrichment = c(1.0, 4.0), replicate = 1
  )
)

# TESTS ########################################################################

if(F) {
  # Simulation of a read count matrix with 3 populations of observations
  # The first population is globally invariable (N = 1000) whereas the two
  # others present negatively correlated variations.
  # By design, the simulation of the invariable population is such that the
  # true value of normalization factors is always equal to 1.0
  p <- c(1000, 4000, 10000)
  m <- DefineSimulation(
    chip = 5, patterns = c("^", "v"), enrichment = c(1.0, 2.5), replicate = 2
  )
  r <-  MakeSimulation(p = p, m = m)

  grp <- r$group # Population memberships
  cnt <- r$data  # Simulated counts

  # Names of the simulated ChIP and Input samples
  chip <- grep("ChIP", colnames(cnt), value = T)
  ctrl <- grep("Input", colnames(cnt), value = T)

  # Prepare figure layout and graphic options
  layout(matrix(1:4, 2, 2, byrow = T))
  par(pch = 20)

  # Show the empirical distribution of simulated populations
  l2c <- log2(DitherCounts(cnt))
  xyl <- range(l2c[FiniteValues(l2c), ])

  r <- PlotCountDistributions(l2c, ylim = xyl, main = "Total")
  for(i in sort(unique(grp))) {
    main <- ifelse(i == 1, "Invariable subset", paste("Variable subset", i - 1))
    r <- PlotCountDistributions(l2c[grp == i, ], ylim = xyl, main = main)
  }

  # Search background candidates knowing only that there should be 3 populations
  brd <- BRD(cnt = cnt, controls = ctrl, ncl = 3, bdt = c(0.1, 0.05))

  PlotBRD(brd) # Show control graphs

  barplot(ScalingFactors(brd), las = 2) # Show normalization factors
}

# BENCHMARKS ###################################################################

benchmarks <- list()
values     <- list()

# =============================================================================.
#
# -----------------------------------------------------------------------------.
sim.prm <- list(
  p = c(2000, 10000),
  m = DefineSimulation(
    chip = 3, patterns = "+", enrichment = c(1.0, 4.0), replicate = 1
  )
)
brd.args <- list(ncl = 2, dither = 5, knn = 150, bdt = c(0.2, 0.1))
# -----------------------------------------------------------------------------.
if(F) Preview(brd.args, sim.prm)
# -----------------------------------------------------------------------------.
if(F) {
  test <- "dither"
  values[[test]] <- 1:10
  n <- length(values[[test]])
  benchmarks[[test]] <- RunTest(
    test, values[[test]], brd.args, sim.prm, reproduce = 15
  )
  layout(matrix(1:4, 2, 2, byrow = T))
  PlotBenchmark(benchmarks[[test]], test, values[[test]])
}
# -----------------------------------------------------------------------------.
if(F) {
  test <- "knn"
  values[[test]] <- c(10, 25, 50, 100, 150, 200, 250, 300, 400, 500)
  n <- length(values[[test]])
  r <- vector("list", n)
  for(i in 1:n) {
    message("Test ", i, " / ", n)
    brd.args[[test]] <- values[[test]][i]
    r[[i]] <- RunBenchmark(brd.args, sim.prm, reproduce = 10)
  }
  benchmarks[[test]] <- r

  layout(matrix(1:4, 2, 2, byrow = T))
  PlotBenchmark(r, test, values = values[[test]])
}

