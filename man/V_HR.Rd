% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation_helpers.R
\name{V_HR}
\alias{V_HR}
\title{Compute exponent measure}
\usage{
V_HR(x, Gamma = NULL, Theta = NULL)
}
\arguments{
\item{x}{Numeric vector with \code{d} positive elements
where the exponent measure is to be evaluated.}

\item{Gamma}{d x d variogram matrix or numeric vector with d(d-1)/2 elements,
containing the upper triangular part of a variogram matrix.}

\item{Theta}{d x d precision matrix or numeric vector with d(d-1)/2 elements,
containing the upper triangular part of a precision matrix.}
}
\value{
Numeric. The exponent measure of the HR distribution.
}
\description{
Computes the exponent measure of HR distribution.
}
\details{
Only \code{Gamma} is needed for the computation. \code{Theta} is only used to
compute \code{Gamma} if necessary.
}
\keyword{internal}
