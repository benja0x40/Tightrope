% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mv_gmm.R
\name{mv_gmm}
\alias{mv_gmm}
\title{Gaussian mixture model fitting}
\usage{
mv_gmm(x, ns, epsilon = 0.001, max_iter = 100)
}
\arguments{
\item{x}{numeric matrix representing multivariate data where rows = observations
and columns = samples or conditions.}

\item{ns}{number of models to be fitted}

\item{epsilon}{to be implemented}

\item{max_iter}{maximum number of iterations (default = 100)}
}
\value{
\code{mv_gmm} returns a list with the following elements:
\item{theta}{list of multivariate distribution parameters}
\item{p}{vector of probabilities}
}
\description{
Gaussian mixture model fitting
}
\keyword{internal}
