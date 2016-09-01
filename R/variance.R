# =============================================================================.
#' voom-like variance weights on log-transformed counts
# -----------------------------------------------------------------------------.
#' @param \code{cnt} matrix of read counts (columns = samples, rows = intervals)
#' @param \code{span}
#' @param \code{plot}
#'
#' @return \code{list}
# -----------------------------------------------------------------------------.
estimateVarianceWeights <- function(cnt, span=0.5, plot=F) {

  # Compute mean and variance per observation
  chk <- finiteValues(cnt)
  mu  <- rowMeans(cnt[chk,], na.rm = T)
  sigma <- apply(cnt[chk,], MARGIN = 1, FUN = sd, na.rm=T)

  # Fit mean-variance relation by LOWESS regression
  nrm <- lowess(mu, sqrt(sigma), f=span)
  nrm$y <- nrm$y^2

  # Normalize variance weights
  refMean    <- nrm$x / mean(mu)
  refWeights <- nrm$y / min(nrm$y, na.rm=T)

  if(plot) {
    plot(mu, sqrt(sigma), pch=20, cex=0.3, col=rgb(0,0,0,0.1))
    lines(nrm$x, sqrt(nrm$y), col=rgb(0.5,0,1), lwd=2)
    plot(refMean, refWeights, type='l', col=rgb(0.5,0,1), lwd=2)
  }
  list(refMean = refMean, refWeights=refWeights)
}

# =============================================================================.
#' variance stabilization on log-transformed counts using voom-like weights
# -----------------------------------------------------------------------------.
#' @param \code{cnt} matrix of read counts (columns = samples, rows = intervals)
#' @param \code{refMean}
#' @param \code{refWeights}
#'
#' @return \code{matrix}
# -----------------------------------------------------------------------------.
correctByVarianceWeights <- function(cnt, refMean, refWeights) {

  # Mask infinite values
  isi <- is.infinite(cnt)
  inf <- cnt[isi]
  cnt[isi] <- NA

  # Compute mean per observation
  mu <- rowMeans(cnt, na.rm=T)

  # Mask NA values
  usa <- ! is.na(mu)
  mu <- mu[usa]

   # Interpolate variance weights scaled to observations
  w <- approx(x = refMean * mean(mu), y = refWeights, xout = mu)$y

  # Apply variance stabilization
  cnt[usa,] <- (cnt[usa,]  - mu) / w + mu

  # Restore infinite values
  cnt[isi]  <- inf
  cnt
}

# =============================================================================.
# DRAFT: generalization to n-dimensions should work but not usable yet.
# -----------------------------------------------------------------------------.
# Observed issues:
# + Projected point P can be slighlty off diagonal
# + baseMean and baseStdv seem to be affected by numerical errors
#   for values close to zero (=> catastrophic cancelation quite likely)
# -----------------------------------------------------------------------------.
# TESTS
# x <- rnorm(100)
# m <- matrix(x, 100, 4) + rnorm(400)/10
# h <- projectionStatistics(m)
# si <- sqrt(apply((m - h$P)^2, MARGIN = 1, FUN = mean))
# s <- (m - h$P) * sqrt(si) + h$P
#
# layout(matrix(1:9,3,3,byrow=T))
# tst <- estimateMeanVarianceRelation(m, plot = T)
#
# plot(m, xlim=range(m), ylim=range(m))
# abline(0,1)
# segments(m[,1], m[,2], h$P[,1], h$P[,2])
# segments(h$P[,1], h$P[,2], s[,1], s[,2], col='red')
# arrows(h$P[,1], h$P[,2], s[,1], s[,2], col='red')
# -----------------------------------------------------------------------------.
# projectionStatistics <- function(M, v=NULL) {
#   n <- ncol(M)
#   m <- nrow(M)
#   # Default projection = n-dimension diagonal (e.g. y = x in 2D)
#   if(is.null(v)) {
#     v <- rep(1/sqrt(n), n) # sqrt(n*(x^2))) = 1 => x = 1/sqrt(n)
#   }
#   M.avg <- apply(M, MARGIN = 2, FUN = mean)
#   M_ctr <- t(t(M) - M.avg) # Centered matrix
#   P <- (M_ctr %*% v) %*% v # Projection
#   P <- scale(P, center = - M.avg, scale=F)
#   s <- sign(apply(P, MARGIN = 1, FUN = min))
#   baseMean <- s * sqrt(apply(P^2, MARGIN = 1, FUN = sum))
#   baseStdv <- sqrt(apply((M - P)^2, MARGIN = 1, FUN = mean))
#   list(P=P, baseMean=baseMean, baseStdv=baseStdv)
# }

# =============================================================================.
# DRAFT: generalization to n-dimensions should work but not usable yet.
# -----------------------------------------------------------------------------.
# stabilizeVariance <- function(cnt, refMean, refWeights, sigma=1) {
#
#   mu <- apply(cnt, MARGIN = 1, FUN = mean, na.rm=T)
#   cnt_proj <- projectionStatistics(cnt)
#   baseMean <- cnt_proj$baseMean
#
#   # Compute the variance correction factor for all M values
#   refWeights <- pmax(refWeights, sigma)
#   w <- approx(x=refMean, y=refWeights, xout=mu, rule=2, ties='ordered')$y
#   # w <- refWeights
#   # No variance correction applises below the mode of stddev ?
#   # Apply variance correction
#   # cnt <- sigma/w * (cnt - cnt_proj$P) + cnt_proj$P
#   cnt <- sigma/w * (cnt - mu) + mu
#   cnt
# }

# =============================================================================.
# DRAFT: generalization to n-dimensions should work but not usable yet.
# -----------------------------------------------------------------------------.
# estimateMeanVarianceRelation <- function(cnt, smooth=0.4, sqroot=T, plot=F) {
#   cnt_proj <- projectionStatistics(cnt)
#   baseMean <- cnt_proj$baseMean
#   baseStdv <- cnt_proj$baseStdv
#
#   mu <- apply(cnt, MARGIN = 1, FUN = mean, na.rm=T)
#   si <- apply(cnt, MARGIN = 1, FUN = sd, na.rm=T)
#
#   # Detect mode of the standard deviation
#   d <- density(baseStdv)
#   smod <- d$x[which.max(d$y)]
#   alt.d <- density(si)
#   alt.smod <- alt.d$x[which.max(alt.d$y)]
#   # LOWESS regression
#   if(sqroot) {
#     alt <- lowess(mu, sqrt(si), f=smooth)
#     alt$y <- alt$y^2
#     nrm <- lowess(baseMean, sqrt(baseStdv), f=smooth)
#     nrm$y <- nrm$y^2
#   } else {
#     alt <- lowess(mu, si, f=smooth)
#     nrm <- lowess(baseMean, baseStdv, f=smooth)
#   }
#   if(plot) {
#     plot(mu, si, pch=20, cex=0.3, col=rgb(0,0,0,0.5))
#     abline(h=smod, col=rgb(0,0,0,0.5), lwd=2)
#     lines(alt$x, alt$y, col=rgb(1,0,0.5,0.8), lwd=3)
#     legend("topright", c("lowess", "mode"), fill=c(rgb(1,0,0.5), rgb(0,0,0,0.5)), bty='n')
#
#     plot(baseMean, baseStdv, pch=20, cex=0.3, col=rgb(0,0,0,0.5))
#     abline(h=smod, col=rgb(0,0,0,0.5), lwd=2)
#     lines(nrm$x, nrm$y, col=rgb(1,0,0.5,0.8), lwd=3)
#     legend("topright", c("lowess", "mode"), fill=c(rgb(1,0,0.5), rgb(0,0,0,0.5)), bty='n')
#
#     plot(mu, baseMean, pch=20, cex=0.3, col=rgb(0,0,0,0.5))
#     abline(0,1,col=rgb(1,0,0,0.5))
#
#     plot(alt.d, xlab="Standard deviation", main="si")
#     abline(v=alt.smod, col=rgb(1,0,0,0.5), lwd=2)
#     legend("topright", paste("mode =", smod), fill="red", bty='n')
#
#     plot(d, xlab="Standard deviation", main="baseStdv")
#     abline(v=smod, col=rgb(1,0,0,0.5), lwd=2)
#     legend("topright", paste("mode =", smod), fill="red", bty='n')
#
#     plot(si, baseStdv, pch=20, cex=0.3, col=rgb(0,0,0,0.5))
#     abline(0,1,col=rgb(1,0,0,0.5))
#
#
#   }
#   list(baseMean=mu, baseStdv=si, sigma=alt.smod, x=alt$x, y=alt$y)
# }
