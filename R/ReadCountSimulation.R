# =============================================================================.
#' Define different populations for ChIP-seq simulations
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{MakeSimulation}
# -----------------------------------------------------------------------------.
#' @param chip
#' number of ChIP profiles (default = 3).
#'
#' @param replicate
#' number of replicates per ChIP profile (default = 1, no replicate).
#'
#' @param input
#' number of Input profiles (default = 2).
#'
#' @param enrichment
#' numeric vector defining the range of ChIP enrichements (default = c(1, 3)).
#'
#' @param patterns
#' character vector, made only of symbols "+", "-", "^" and "v".
#'
#' @return
#' \code{DefineSimulation} returns a \code{matrix} of parameters.
# -----------------------------------------------------------------------------.
#' @export
DefineSimulation <- function(
  chip = 3, replicate = 1, input = 2, enrichment = c(1.0, 3.0), patterns = "^"
) {

  # Invariable population
  x <- rep(1.0, chip)

  # Variable populations
  for(p in patterns) {
    if(p == "+") s <- (1:chip - 1) / (chip - 1)
    if(p == "-") s <- (chip:1 - 1) / (chip - 1)
    if(p == "^") s <- 1 - 2 * abs((1:chip - 1) / (chip - 1) - 0.5)
    if(p == "v") s <- 0 + 2 * abs((1:chip - 1) / (chip - 1) - 0.5)
    s <- s * diff(enrichment) + enrichment[1]
    x <- rbind(x, s)
  }
  rownames(x) <- NULL

  # Replicates
  x <- x[, rep(1:chip, replicate)]

  # Name ChIP columns
  a <- rep(paste0("ChIP_", LETTERS[1:chip]), replicate)
  if(replicate > 1) a <- paste0(a, gl(replicate, chip))
  colnames(x) <- a

  # Inputs columns
  m <- matrix(1.0, nrow(x), input)
  colnames(m) <- paste0("Input_", 1:input)
  m <- cbind(x, m)

  m
}

# =============================================================================.
#' Simulate ChIP-seq read count matrixes
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{DefineSimulation}
# -----------------------------------------------------------------------------.
#' @description
#' The \code{MakeSimulation} function generates random read count matrixes
#' presenting distinct populations accross simulated ChIP and Input samples.
#' By design the first population is globally invariable, such that the true
#' value of normalization factors between the simulated conditions is always
#' equal to 1.0. The simulated variable populations can be defined using the
#' \link{DefineSimulation} function.
#'
#' @param p
#' integer vector of population sizes, the first one corresponding to the
#' invariable population.
#'
#' @param m
#' matrix of parameters for the variable populations, which can be produced by
#' calling the \link{DefineSimulation} function.
#'
#' @param f
#' list of functions allowing to use custom generators for simulated
#' measurements (default = NULL, recommended).
#'
#' @return
#' \code{MakeSimulation} returns a \code{list} with the following elements:
#' \item{data}{matrix of simulated ChIP-seq read counts.}
#' \item{group}{population memberships.}
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
MakeSimulation <- function(p = NULL, m = NULL, f = NULL) {

  noise <- function(n) { 0.0 + rnbinom(n, size = 10, mu =  5) }
  alpha <- function(n) { 0.0 + rnbinom(n, size =  1, mu = 15) }
  gamma <- function(n) { 1.0 + rnbinom(n, size =  5, mu =  5) }

  n <- max(1, length(p), length(f), nrow(m))

  if(is.null(p)) p <- rep(10000, n)
  if(is.null(m)) m <- matrix(1.0, n, 2)
  if(is.null(f)) f <- list(alpha, gamma)

  if(length(p) != n | nrow(m) != n) stop("inconsistent arguments")

  k <- c()
  r <- c()
  for(i in 1:n) {
    j <- 1 + i > 0                # f[[1]] used once, then f[[2]]
    j <- 1 + (i - 1) %% length(f) # Cyclic use of f values
    j <- min(i, length(f))        # Acyclic use of f values

    a <- m[i, ]
    b <- f[[j]](p[i])
    r <- rbind(r, b %*% t(a))
    k <- c(k, rep(i, p[i]))
  }
  r <- r + matrix(noise(length(r)), nrow(r), ncol(r))

  if(! is.null(colnames(m))) colnames(r) <- colnames(m)

  list(data = round(r), group = k)
}

# HIDDEN #######################################################################

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
#' logical (default = FALSE)
#'
#' @return
#' \code{RandomCounts} returns a list.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
RandomCounts <- function(n, d, k = NULL, p = NULL, v = NULL, extended = FALSE) {

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
    sample(1:k, size = n * d, replace = TRUE, prob = p),
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
    z <- matrixStats::rowMeans2(apply(m, 2, function(i) p[i]))
    z <- z / sum(z)
    # Value substitution and result
    m <- matrix(v[m], k^2, 2)
    r <- c(r, list(m = m, d = z))
  }

  r
}
