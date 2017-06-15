# GEOMETRY #####################################################################

# =============================================================================.
#' translation_2D
# -----------------------------------------------------------------------------.
#' @param v matrix
#' @param dv vector
#'
#' @return matrix
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
translation_2D <- function(v, dv) {
  v + matrix(dv, nrow(v), length(dv), byrow = T)
}

# =============================================================================.
#' RotationMatrix2D
# -----------------------------------------------------------------------------.
#' @param alpha angle
#'
#' @return matrix
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
RotationMatrix2D <- function(alpha) {
  r <- diag(rep(0, 2))
  r[1, ] <- c(cos(alpha), -sin(alpha))
  r[2, ] <- c(sin(alpha),  cos(alpha))
  r
}

# PARAMETERS ###################################################################

# =============================================================================.
#' CovMat2D
# -----------------------------------------------------------------------------.
#' @param alpha angle in radians
#' @param sx sigma2(x)
#' @param sy sigma2(y)
#'
#' @return covariance matrix
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
CovMat2D <- function(alpha, sx, sy) {
  r <- RotationMatrix2D(alpha)
  s <- diag(c(sx, sy))
  sigma <- r %*% s %*% t(r)
  sigma
}

# 2D PATTERNS ##################################################################

# =============================================================================.
#' make_square_grid_2D
# -----------------------------------------------------------------------------.
#' @param n number
#' @param normalized logical (default = T)
#'
#' @return matrix
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
make_square_grid_2D <- function(n, normalized = T) {
  x <- rep(1:n, n) - (n + 1) / 2
  y <- rep(1:n, each = n) - (n + 1) / 2
  xy <- cbind(x, y)
  if(normalized) xy <- 2 / n * xy
  xy
}

# =============================================================================.
#' make_disk_grid_2D
# -----------------------------------------------------------------------------.
#' @param n number
#' @param normalized logical (default = T)
#'
#' @return matrix
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
make_disk_grid_2D <- function(n, normalized = T) {
  x <- rep(1:n, n) - (n + 1) / 2
  y <- rep(1:n, each = n) - (n + 1) / 2
  d <- sqrt(x^2 + y^2)
  xy <- cbind(x, y)[d < n / 2, ]
  if(normalized) xy <- 2 / n * xy
  xy
}

# =============================================================================.
#' make_ring_grid_2D
# -----------------------------------------------------------------------------.
#' @param n number
#' @param normalized logical (default = T)
#'
#' @return matrix
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
make_ring_grid_2D <- function(n, normalized = T) {
  x <- rep(1:n, n) - (n + 1) / 2
  y <- rep(1:n, each = n) - (n + 1) / 2
  d <- sqrt(x^2 + y^2)
  xy <- cbind(x, y)[d < n / 2 & d > n / 4, ]
  if(normalized) xy <- 2 / n * xy
  xy
}

# =============================================================================.
#' make_square_unif_2D
# -----------------------------------------------------------------------------.
#' @param n number
#' @param normalized logical (default = T)
#'
#' @return matrix
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
make_square_unif_2D <- function(n, normalized = T) {
  x <- n * (runif(n * n) - 1 / 2)
  y <- n * (runif(n * n) - 1 / 2)
  xy <- cbind(x, y)
  if(normalized) xy <- 2 / n * xy
  xy
}

# =============================================================================.
#' make_disk_unif_2D
# -----------------------------------------------------------------------------.
#' @param n number
#' @param normalized logical (default = T)
#'
#' @return matrix
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
make_disk_unif_2D <- function(n, normalized = T) {
  x <- n * (runif(n * n) - 1 / 2)
  y <- n * (runif(n * n) - 1 / 2)
  d <- sqrt(x^2 + y^2)
  xy <- cbind(x, y)[d < n / 2, ]
  if(normalized) xy <- 2 / n * xy
  xy
}

# =============================================================================.
#' make_ring_unif_2D
# -----------------------------------------------------------------------------.
#' @param n number
#' @param normalized logical (default = T)
#'
#' @return matrix
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
make_ring_unif_2D <- function(n, normalized = T) {
  x <- n * (runif(n * n) - 1 / 2)
  y <- n * (runif(n * n) - 1 / 2)
  d <- sqrt(x^2 + y^2)
  xy <- cbind(x, y)[d < n / 2 & d > n / 4, ]
  if(normalized) xy <- 2 / n * xy
  xy
}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
# n <- 100
# layout(matrix(1:9, 3, 3, byrow = T))
# plot(make_square_grid_2D(n), pch = 20, col = grey(0, alpha = 0.5))
# plot(make_disk_grid_2D(n), pch = 20, col = grey(0, alpha = 0.5))
# plot(make_ring_grid_2D(n), pch = 20, col = grey(0, alpha = 0.5))
# plot(make_square_unif_2D(n), pch = 20, col = grey(0, alpha = 0.5))
# plot(make_disk_unif_2D(n), pch = 20, col = grey(0, alpha = 0.5))
# plot(make_ring_unif_2D(n), pch = 20, col = grey(0, alpha = 0.5))

# ND GENERATORS ################################################################

# =============================================================================.
#' SquareGrid
# -----------------------------------------------------------------------------.
#' @param k integer
#' @param d dimension
#' @param radius distance of exclusion
#' @param subset vector of indices
#'
#' @return list
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
SquareGrid <- function(k, d, radius = -1, subset = NULL) {
  x <- matrix(0, k^d, d)
  if(k > 1) {
    for(i in 1:d) {
      x[, i] <- as.numeric(gl(k, k^(i-1), length = k^d))
    }
    x <- 2 * ((x - 1) / (k - 1) - 1 / 2)
  }
  if(! is.null(subset)) {
    x <- x[subset, , drop = F]
  }
  r <- apply(abs(x), 1, max)
  chk <- r > radius
  x <- x[chk, , drop = F]
  list(n = nrow(x), d = d, k = k, x = x)
}
# =============================================================================.
#' RoundGrid
# -----------------------------------------------------------------------------.
#' @param k integer
#' @param d dimension
#' @param radius distance of exclusion
#' @param subset vector of indices
#'
#' @return list
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
RoundGrid <- function(k, d, radius = -1, subset = NULL) {
  m <- SquareGrid(k, d)
  r <- sqrt(rowSums(m$x^2))
  chk <- r <= 1 & r > radius
  m$x <- m$x[chk, , drop = F]
  if(! is.null(subset)) {
    m$x <- m$x[subset, , drop = F]
  }
  m$n <- nrow(m$x)
  m
}
# =============================================================================.
#' ClonalGaussian
# -----------------------------------------------------------------------------.
#' @param n number
#' @param mu mean
#' @param sigma strandard deviation
#'
#' @return list
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ClonalGaussian <- function(n, mu, sigma) {
  mu <- as.matrix(mu)
  d <- ncol(mu)
  k <- nrow(mu)
  if(length(sigma) == 1) sigma <- diag(sigma, d, d)
  x <- matrix(0, n * k, d)
  for(i in 1:k) {
    g <- 1:n + (i - 1) * n
    x[g, ] <- rmvnorm(n, mu = mu[i, ], sigma = sigma)
  }
  s <- rep(1:k, each = n)
  list(n = n * k, d = d, k = k, src = s, x = x)
}
# =============================================================================.
#' ClonalUniform
# -----------------------------------------------------------------------------.
#' @param n number
#' @param mu mean
#' @param delta half width of the range
#'
#' @return list
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ClonalUniform <- function(n, mu, delta) {
  mu <- as.matrix(mu)
  d <- ncol(mu)
  k <- nrow(mu)
  if(length(delta) == 1) delta <- matrix(delta, k, d)
  x <- matrix(0, n * k, d)
  for(i in 1:k) {
    g <- 1:n + (i - 1) * n
    for(j in 1:d) {
      x[g, j] <- runif(
        n, min = mu[i, j] - delta[i, j], max = mu[i, j] + delta[i, j]
      )
    }
  }
  s <- rep(1:k, each = n)
  list(n = n * k, d = d, k = k, src = s, x = x)
}

# COUNTS #######################################################################

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
