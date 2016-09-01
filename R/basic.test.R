# =============================================================================.
#' Fit normal distribution to the data and compute associated p-values
# -----------------------------------------------------------------------------.
#' @param m
#' @param invariant
#' @param txt
#' @param plot
#'
#' @return \code{list}
# -----------------------------------------------------------------------------.
fitNormal <- function(m, invariant=NULL, txt="", plot=F, brks.len=70) {

  ok    <- finiteValues(m)

  # fit normal model using mixture of 2 normal distributions
  mfit <- normalmixEM(m[ok])
  idx <- which.min(abs(mfit$mu))
  mu <- mfit$mu[idx]
  sigma <- mfit$sigma[idx]

  # Compute p-value and FDR based on the fitted normal model
  pval <- rep(NA, length(m))
  padj <- rep(NA, length(m))
  pval[ok] <- pnorm(abs(m[ok]), mean=mu, sd=sigma, log.p = F, lower.tail = F)
  padj[ok] <- p.adjust(pval[ok], method = "fdr")

  if(plot) {
    mlim  <- range(m[ok])
    blen <- brks.len
    brks <- seq(mlim[1], mlim[2], length.out = blen)
    h <- hist(m[ok], breaks=brks, col=rgb(0,0,0,0.5), border=rgb(0,0,0,0), xlab="log2 ratio", main=txt)
    x <- seq(-10, 10, length.out = 2000)
    y <- dnorm(x, mean=mu, sd=sigma)
    y <- y/max(y) * quantile(h$counts, probs=(blen-2)/blen)
    lines(x, y, col=rgb(1,0,0,0.7), lwd=3)
    legend("topright", c("observed", "normal model"), fill=c("grey","red"), bty='n')

    # Make qqplot of model versus data
    ymir  <- c(m[ok & invariant], - m[ok & invariant])
    ymir <- ymir[ymir!=0]
    qmod <- qnorm((1:100-0.5)/100, mean = mu, sd=sigma)
    qdat <- quantile(ymir, probs=(1:100-0.5)/100)
    plot(qmod, qdat, xlab="model", ylab="data", main="Q-Q plot")
    abline(0,1, col=rgb(0,0,0,0.5), lwd=2)
  }

  list(mu=mu, sigma=sigma, pval=pval, padj=padj)
}
# =============================================================================.
#' Categorize differences according to fold-changes
# -----------------------------------------------------------------------------.
#' @param a
#' @param m
#' @param pval
#' @param padj
#' @param bg.tail
#' @param cmp
#' @param bg
#' @param h
#' @param amin
#'
#' @return \code{data.frame}
# -----------------------------------------------------------------------------.
statusByFoldChange <- function(a, m, pval, padj, bg.tail, cmp, bg=log2(1.5), h=log2(2), amin=4) {
  status <- rep("undefined", length(m))
  if(bg.tail=="none") {
    status[abs(m) < bg] <- "background"
  }
  if(bg.tail=="positive") {
    status[m > - bg] <- "background"
  }
  if(bg.tail=="negative") {
    status[m <   bg] <- "background"
  }
  status[a > amin & m >   h] <- "higher"
  status[a > amin & m < - h] <- "lower"
  res <- data.frame(a, m, pval, padj, status)
  colnames(res) <- paste(cmp, c("avglog2", "log2ratio", "pval", "padj", "status"), sep="_")
  res
}
# =============================================================================.
#' Categorize differences according to p-values
# -----------------------------------------------------------------------------.
#' @param a
#' @param m
#' @param pval
#' @param padj
#' @param bg.tail
#' @param cmp
#' @param bg
#' @param h
#' @param amin
#'
#' @return \code{data.frame}
# -----------------------------------------------------------------------------.
statusByPValue <- function(a, m, pval, padj, bg.tail, cmp, bg=0.2, h=1E-02, amin=4) {
  status <- rep("undefined", length(m))
  if(bg.tail=="none") {
    status[padj > bg] <- "background"
  }
  if(bg.tail=="positive") {
    k <-  - max(abs(m[padj > bg]), na.rm=T)
    status[m > k] <- "background"
  }
  if(bg.tail=="negative") {
    k <-  max(abs(m[padj > bg]), na.rm=T)
    status[m < k] <- "background"
  }

  status[a > amin & m > 0 & padj < h] <- "higher"
  status[a > amin & m < 0 & padj < h] <- "lower"
  res <- data.frame(a, m, pval, padj, status)
  colnames(res) <- paste(cmp, c("avglog2", "log2ratio", "pval", "padj", "status"), sep="_")
  res
}
