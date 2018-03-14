library(Tightrope)

# FUNCTIONS ####################################################################

# =============================================================================.
#
# -----------------------------------------------------------------------------.
PlotDistribution <- function(x, ...) {
  SideBySideDensity(
    x, las = 2, parameters = list(color = "Wry"),
    method = "ash", ash = list(m = c(3, 3)), ...
  )
}

# TESTS ########################################################################

p <- c(2000, 10000, 10000)
m <- DefineSimulation(
  chip = 5, patterns = c("^", "v"), enrichment = c(1.0, 4.0), replicate = 2
)
ncl <- 3

p <- c(2000, 10000, 10000)
m <- DefineSimulation(
  chip = 7, patterns = c("^", "v"), enrichment = c(1.0, 3.0), replicate = 1
)
ncl <- 3

p <- c(2000, 10000, 10000)
m <- DefineSimulation(
  chip = 5, patterns = c("+", "-"), enrichment = c(1.0, 3.0), replicate = 1
)
ncl <- 2


r <- MakeSimulation(p = p, m = m)

grp <- r$group
cnt <- r$data

chip <- grep("ChIP", colnames(cnt), value = T)
ctrl <- grep("Input", colnames(cnt), value = T)

l2c <- log2(DitherCounts(cnt))
l   <- xylim(l2c, symetric = T)

brd <- BRD(cnt, controls = ctrl, ncl = ncl, dither = 10, knn = 300)

layout(matrix(1:4, 2, 2, byrow = T))
par(pch = 20)
PlotBRD(brd, with.legend = F)

layout(matrix(1:4, 2, 2, byrow = T))
r <- PlotDistribution(l2c, ylim = l$y, main = "Total")
for(i in sort(unique(grp))) {
  r <- PlotDistribution(l2c[grp == i, ], ylim = l$y, main = paste("Group", i))
}

# layout(matrix(1:9, 3, 3, byrow = T))
# input <- rowMeans(l2c[, ctrl])
# for(lbl in chip) {
#   r <- BivariateDensity(
#     input, l2c[, lbl], xlim = l$x, ylim = l$y, parameters = list(color = "Wry"),
#     method = "ash", ash = list(m = c(3, 3)), xlab = "average Input", ylab = lbl
#   )
#   abline(a = 0, b = 1, col = grey(0.5))
# }

# 2^brd$normfactors

