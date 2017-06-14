# =============================================================================.
#' Random counts
# -----------------------------------------------------------------------------.
#' @description
#' Generate pseudo random count values
#'
#' @param n
#' number of observations (rows).
#'
#' @param d
#' number of dimensions (columns).
#'
#' @param k
#' by default count values are sampled in [1, \code{k}].
#' This argument is ignored when arguments \code{p} or \code{v} are provided.
#'
#' @param p
#' vector of count probabilities (default = uniform).
#'
#' @param v
#' vector of count values (default = \code{1:k})
#'
#' @param extended
#' logical (default = F)
#'
#' @return
#' RandomCounts returns a \code{list}
# -----------------------------------------------------------------------------.
#' @export
RandomCounts <- function(n, d, k = NULL, p = NULL, v = NULL, extended = F) {

  chk <- sum(2^(0:2) * (! sapply(list(k, p, v), is.null)))

  if(chk == 0) {
    stop("provide at least one of the 'k', 'p' or 'v' arguments")
  }

  # p != NULL, v == NULL
  if(chk %in% (0:1 + 2)) k <- length(p)

  # v != NULL, p == NULL
  if(chk %in% (0:1 + 4)) k <- length(v)

  # k == NULL
  if(chk %in% (0:1 + 6)) k <- length(p) # assuming length(p) == length(v)

  # p == NULL
  if(! bitwAnd(chk, 2)) p <- rep(1, k)

  # v == NULL
  if(! bitwAnd(chk, 4)) v <- 1:k

  if(k != length(p) | k != length(v)) {
    stop("arguments 'p' and 'v' have different lengths")
  }
  if(k < 2) {
    stop("generated counts must have at least two distinct values")
  }

  p <- p / sum(p)

  x <- matrix(
    sample(1:k, size = n * d, replace = T, prob = p),
    nrow = n, ncol = d
  )

  # Value substitution and result
  x <- matrix(v[x], nrow = n, ncol = d)
  r <- list(x = x, k = k, p = p, v = v)

  if(extended) {
    m <- cbind(
      rep(1:k, k),
      rep(1:k, each = k)
    )
    z <- rowMeans(apply(m, 2, function(i) p[i]))
    z <- z / sum(z)
    # Value substitution and result
    m <- matrix(v[m], k^2, 2)
    r <- c(r, list(m = m, d = z))
  }

  r
}

# =============================================================================.
#' Random counts representing 3 populations
# -----------------------------------------------------------------------------.
# 1 = poorcounts : low mappability
# 2 = invariable : background or low variation
# 3 = variable   : specific signal with significant variations
# -----------------------------------------------------------------------------.
#' @param n
#' number of observations
#'
#' @param d
#' number of conditions (dimensions)
#'
#' @param mu
#' average read counts (provided as log transformed value)
#'
#' @param sigma
#'
#' @param poorcounts
#' parameters of the barely detectable subpopulation
#'
#' @param invariable
#' parameters of the subpopulation generating background or low variability
#' counts
#'
#' @param variable
#' parameters of the subpopulation showing significant variations
#'
#' @param as.counts
#' return counts (default = T, yes)
#'
#' @return
#' SimulatedCounts returns a list
# -----------------------------------------------------------------------------.
#' @export
SimulatedCounts <- function(
  n, d = 3, mu = 5, sigma = 3,
  poorcounts = list(), invariable = list(), variable = list(), as.counts = T
) {

  prm <- data.frame(
    row.names = c("poorcounts", "invariable", "variable"),
    id     = 1:3,
    n      = 0,
    alpha  = c(0.05, 0.35, 0.60),
    # beta   = c(0.00, -2/3, -1/10),
    beta   = c(-1/4, -1/4,  1/2),
    gamma  = c(0.50, 0.50, 0.50),
    sx     = c(  -1,   -1,   -1),
    sy     = c(   1,    1,    1),
    mu     = c(0.50, 0.00, 2.00),
    sigma  = c(0.50, 1.00, 0.50)
  )

  for(rn in rownames(prm)) {
    for(cn in colnames(prm)) {
      xp <- eval(parse(text = paste0(rn, "$", cn)))
      if(! is.null(xp)) prm[rn, cn] <- xp
    }
  }
  prm$alpha <- prm$alpha / sum(prm$alpha)
  grp <- sample(1:3, size = n, replace = T, prob = prm$alpha)
  prm$n <- table(grp)

  k <- sapply(1:3, "==", grp)

  p <- cbind(
    with(
      prm["poorcounts", ],
      SkewedProbabilities(
        k = 1000, beta = beta, gamma = gamma, sx = sx, sy = sy
      )
    ),
    with(
      prm["invariable", ],
      SkewedProbabilities(
        k = 1000, beta = beta, gamma = gamma, sx = sx, sy = sy
      )
    ),
    with(
      prm["variable", ],
      SkewedProbabilities(
        k = 1000, beta = beta, gamma = gamma, sx = sx, sy = sy
      )
    )
  )

  a <- sort(1 + abs(rnorm(n, mean = mu, sd = sigma)))
  s <- seq.int(1, ceiling(2^max(a)))
  p <- cbind(
    MatchLogSpace(p[, 1], s), MatchLogSpace(p[, 2], s), MatchLogSpace(p[, 3], s)
  )
  p <- p[ceiling(2^a), ]

  u <- (1.1 * max(a) - a) / (1.1 * max(a) - min(a))
  x <- matrix(0, n, d)
  for(i in 1:d) {
    x[k[, 1], i] <- a[k[, 1]] / mu + with(
      prm["poorcounts", ], p[k[, 1], 1] * rnorm(n, mu, sigma)
    )
    x[k[, 2], i] <- a[k[, 2]] + with(
      prm["invariable", ], p[k[, 2], 2] * rnorm(n, mu * 1, sigma * 1)
    )
    x[k[, 3], i] <- a[k[, 3]] + with(
      prm["variable", ], u[k[, 3]] * p[k[, 3], 3] * rnorm(n, mu * i, sigma * i)
    )
    # varying sequencing yields
    x[, i] <- x[, i] + i * (- 1)^i
  }
  colnames(x) <- c("Input", paste0("ChIP_", 1:(d-1)))

  if(as.counts) x <- ceiling(2^x)
  list(prm = prm, grp = grp, x = x)
}
