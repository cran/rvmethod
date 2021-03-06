% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gaussfit.R
\name{gauss1func}
\alias{gauss1func}
\title{A Single Gaussian Absorption Feature}
\usage{
gauss1func(x, a1, mu1, sig1)
}
\arguments{
\item{x}{the vector of values at which to evaluate}

\item{a1}{the amplitude of the feature}

\item{mu1}{the center of the feature}

\item{sig1}{the spread of the feature (must be greater than 0)}
}
\value{
vector of values of the specified inverted Gaussian
}
\description{
This function returns a Gaussian absorption feature with continuum 1.0 and a
specified amplitude, center, and spread.
}
\examples{
x = seq(5000, 5003, length.out=200)
y = gauss1func(x, 0.3, 5001.5, 0.1)
plot(x, y)

}
