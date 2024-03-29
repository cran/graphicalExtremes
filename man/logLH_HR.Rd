% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation_helpers.R
\name{logLH_HR}
\alias{logLH_HR}
\title{Full censored log-likelihood of HR model}
\usage{
logLH_HR(data, Gamma, cens = FALSE)
}
\arguments{
\item{data}{Numeric \nxd matrix, containing
observations following a multivariate HR Pareto distribution.}

\item{Gamma}{Numeric \dxd matrix, representing a variogram matrix \eGamma.}

\item{cens}{Boolean. If true, then censored log-likelihood is computed.
By default, \code{cens = FALSE}.}
}
\value{
Numeric. The full censored log-likelihood of HR model.
}
\description{
Computes the full (censored) log-likelihood of HR model.
}
\keyword{internal}
