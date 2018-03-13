library(Tightrope)

# x <- c(
#   H3K27me3_WT   = -0.08283696,
#   H3K27me3_K27M = 1.08049920,
#   Input_WT      = -0.63129889,
#   Input_K27M    = -0.36636335
# )
#

nx <- 10000
ny <- 20000

x <- 0.0 + rnbinom(nx, size = 1, mu = 15)
y <- 1.0 + rnbinom(ny, size = 5, mu = 10)
k <- cbind(
  ChIP_A  = c(x, 3.0 * y) + rnbinom(nx + ny, size = 10, mu = 5),
  ChIP_B  = c(x, 2.0 * y) + rnbinom(nx + ny, size = 10, mu = 5),
  Input_A = c(x, 1.0 * y) + rnbinom(nx + ny, size = 10, mu = 5),
  Input_B = c(x, 1.0 * y) + rnbinom(nx + ny, size = 10, mu = 5)
)
k <- round(k)

k <- t(t(k) * 2:5/2)


layout(matrix(1:4, 2, 2, byrow = T))
brd <- BRD(k, controls = colnames(k)[2:3], ncl = 2)
par(pch = 20)
PlotBRD(brd, with.legend = F)


k <- log2(DitherCounts(k))
l <- xylim(k, symetric = T)

r <- SideBySideDensity(k, parameters = list(color = "Wry"))
r <- BivariateDensity(k[1:n, ], xlim = l$x, ylim = l$y, parameters = list(color = "Wry"))
abline(a = 0, b = 1, col = grey(0.5))
r <- BivariateDensity(k[1:n+n, ], xlim = l$x, ylim = l$y, parameters = list(color = "Wry"))
abline(a = 0, b = 1, col = grey(0.5))
r <- BivariateDensity(k, xlim = l$x, ylim = l$y, parameters = list(color = "Wry"))
abline(a = 0, b = 1, col = grey(0.5))




y <- x + rnbinom(n, size = 10, mu = 5)
y <- log2(DitherCounts(y))
k <- sample(which(y > 3 & y < 6), size = h)

a <- x
b <- x
a[k] <- a[k] + abs(rnorm(h, 20, 10))
b[k] <- b[k] + abs(rnorm(h, 10, 10))

cnt <- cbind(
  ChIP_A  = 0.9 * (a + rnbinom(n, size = 10, mu = 5)),
  ChIP_B  = 1.8 * (b + rnbinom(n, size = 10, mu = 5)),
  Input_A = 0.6 * (x + rnbinom(n, size = 10, mu = 5)),
  Input_B = 0.8 * x + rnbinom(n, size = 10, mu = 5)
)
cnt <- round(cnt)

brd <- BRD(cnt, controls = colnames(cnt)[2:3], ncl = 2)
PlotBRD(brd)

x <- log2(DitherCounts(cnt))

layout(matrix(1:4, 2, 2, byrow = T))
r <- SideBySideDensity(x, parameters = list(color = "Wry"))
r <- BivariateDensity(x[, c(3, 4)], parameters = list(color = "Wry"))
r <- BivariateDensity(x[, c(1, 2)], parameters = list(color = "Wry"))
r <- BivariateDensity(x[, c(1, 3)], parameters = list(color = "Wry"))





