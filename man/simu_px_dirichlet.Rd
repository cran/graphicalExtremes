% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation_helpers.R
\name{simu_px_dirichlet}
\alias{simu_px_dirichlet}
\title{Simulate Dirichlet extremal functions}
\usage{
simu_px_dirichlet(n, d, idx, alpha)
}
\arguments{
\item{n}{Number of simulations.}

\item{d}{Dimension of the multivariate Pareto distribution.
In some cases this can be \code{NULL} and will be inferred from \code{par}.}

\item{idx}{Integer or numeric vector with \code{n} elements. Index/indices from
1 to \code{d}, that determine which extremal function to simulate.}

\item{alpha}{Numeric vector of size \code{d}.}
}
\value{
Numeric \nxd matrix. Simulated data.
}
\description{
Simulates Dirichlet extremal functions
}
\keyword{internal}
