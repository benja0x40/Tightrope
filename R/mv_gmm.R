# FUNCTIONS | MV MODEL #########################################################

# =============================================================================.
#' Multivariate normal distribution: maximum likelihood estimation
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{mv_pdf},
#'   \link{mv_rsg},
#'   \link{mv_gmm}
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param w
#' numeric weight (i.e. probability) associated to observations
#'
#' @param ML
#' logical (reserved, default = TRUE)
#'
#' @return
#' \code{mv_mle} returns a list with the following elements:
#' \item{mu}{mean vector}
#' \item{sigma}{covariance matrix}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
mv_mle <- function(x, w = NULL, ML = TRUE) {

  x <- as.matrix(x)
  n <- nrow(x)

  # default weights = uniform
  if(is.null(w)) w <- rep(1, n)
  w <- w / sum(w)

  # center vector
  mu <- colSums(w * x, na.rm = TRUE)

  # covariance matrix
  x <- sqrt(w) * (x - matrix(1, n) %*% mu)
  sigma <- crossprod(x) # equivalent to: t(x) %*% (x)

  # unbiased estimator (to be consistent with builtin functions: cov, cov.wt)
  if(! ML) sigma <- sigma / (1 - sum(w^2))

  list(mu = mu, sigma = sigma)
}

# =============================================================================.
#' Multivariate normal distribution: probability density function
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{mv_mle},
#'   \link{mv_rsg},
#'   \link{mv_gmm}
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param theta
#' list of normal distribution parameters
#'
#' @return
#' \code{mv_pdf} returns a vector of probabilities.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
mv_pdf <- function(x, theta) {
  mixtools::dmvnorm(x, theta[[1]], theta[[2]])
}

# =============================================================================.
#' Multivariate normal distribution: random sample generator
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{mv_pdf},
#'   \link{mv_mle},
#'   \link{mv_gmm}
# -----------------------------------------------------------------------------.
#' @param n
#' number of observations
#'
#' @param theta
#' list of normal distribution parameters
#'
#' @return
#' \code{mv_rsg} returns a numeric matrix representing pseudo random data
#' where rows = observations and columns = samples or conditions.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
mv_rsg <- function(n, theta) {
  ns <- length(theta)                     # number of sources
  nv <- length(theta[[1]][[2]])           # number of variables (or dimensions)
  nx <- round(n * sapply(theta, "[[", 1)) # number of observations
  spl <- data.frame(group = 0, matrix(0, sum(nx), nv), stringsAsFactors = FALSE)
  a <- b <- 0
  for(i in 1:ns) {
    a <- b + 1
    b <- b + nx[i]
    spl[a:b, 1] <- i
    spl[a:b, 1+1:nv] <- mixtools::rmvnorm(
      n = nx[i], theta[[i]][[2]], theta[[i]][[3]]
    )
  }
  spl
}

# =============================================================================.
#' Multivariate normal distribution: parameter initialization
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{mv_gmm},
#'   \link{mv_emloop}
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param ns
#' number of normal distributions to be initialised
#'
#' @return
#' \code{mv_init_param} returns a list of normal distribution parameters.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
mv_init_param <- function(x, ns) {
  x <- as.matrix(x)
  nv <- ncol(x)       # number of variables (or dimensions)
  nx <- nrow(x)       # number of observations
  r <- matrix(x[sample(1:nx, size = ns, replace = FALSE), ], nrow = ns)
  a <- mv_mle(x)
  theta <- vector("list", ns)
  for(i in 1:ns) {
    theta[[i]] <- list(alpha = 1 / ns, mu = r[i, ], sigma = a$sigma / ns)
  }
  theta
}

# FUNCTIONS | MV OPTIMIZATION ##################################################

# =============================================================================.
#' Multivariate expectation
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{mv_maximization},
#'   \link{mv_emloop}
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param theta
#' list of multivariate distribution parameters
#'
#' @param p_fun
#' probability density function (default = \link{mv_pdf})
#'
#' @return
#' \code{mv_expectation} returns a vector of probabilities.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
mv_expectation <- function(x, theta, p_fun = mv_pdf) {
  x <- as.matrix(x)
  ns <- length(theta) # number of sources
  nx <- nrow(x)       # number of observations

  p <- matrix(0, nx, ns)
  for(i in 1:ns) {
    p[, i] <- theta[[i]][[1]] * p_fun(x, theta[[i]][-1])
  }
  p <- p / rowSums(p)
  p
}

# =============================================================================.
#' Multivariate maximization
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{mv_expectation},
#'   \link{mv_emloop}
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param theta
#' list of multivariate distribution parameters
#'
#' @param p
#' vector of probabilities
#'
#' @param mle_fun
#' maximum likelyhood estimator function (default = \link{mv_mle})
#'
#' @return
#' \code{mv_maximization} returns a list of multivariate distribution
#' parameters.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
mv_maximization <- function(x, theta, p, mle_fun = mv_mle) {
  x <- as.matrix(x)
  ns <- length(theta) # number of sources
  nx <- nrow(x)       # number of observations

  k <- colSums(p)
  for(i in 1:ns) {
    theta[[i]][[1]] <- k[i] / sum(k)
    theta[[i]][-1]  <- mle_fun(x, w = p[, i])
  }
  theta
}

# =============================================================================.
#' Multivariate expectation-maximization
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{mv_maximization},
#'   \link{mv_expectation},
#'   \link{mv_gmm}
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param theta
#' list of multivariate distribution parameters
#'
#' @param p_fun
#' probability density function (default = \link{mv_pdf})
#'
#' @param mle_fun
#' maximum likelyhood estimator function (default = \link{mv_mle})
#'
#' @param max_iter
#' maximum number of iterations (default = 100)
#'
#' @return
#' \code{mv_emloop} returns a list with the following elements:
#' \item{theta}{list of multivariate distribution parameters}
#' \item{p}{vector of probabilities}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
mv_emloop <- function(
  x, theta, p_fun = mv_pdf, mle_fun = mv_mle, max_iter = 100
) {
  for(loop in 1:max_iter) {
    # expectation step
    p <- mv_expectation(x, theta, p_fun)
    if(any(is.na(p))) stop("NaN")
    # maximization step
    theta <- mv_maximization(x, theta, p, mle_fun)
  }
  list(theta = theta, p = p)
}

# =============================================================================.
#' Gaussian mixture model fitting
# -----------------------------------------------------------------------------.
#' @param x
#' numeric matrix representing multivariate data where rows = observations
#' and columns = samples or conditions.
#'
#' @param ns
#' number of models to be fitted
#'
#' @param epsilon
#' to be implemented
#'
#' @param max_iter
#' maximum number of iterations (default = 100)
#'
#' @return
#' \code{mv_gmm} returns a list with the following elements:
#' \item{theta}{list of multivariate distribution parameters}
#' \item{p}{vector of probabilities}
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
mv_gmm <- function(x, ns, epsilon = 1E-3, max_iter = 100) {

  theta <- mv_init_param(x, ns)
  mv_emloop(x, theta, p_fun = mv_pdf, mle_fun = mv_mle, max_iter = max_iter)
}
