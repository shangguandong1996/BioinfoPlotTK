% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scale.R
\name{scale_zscale}
\alias{scale_zscale}
\title{scale_zscale}
\usage{
scale_zscale(data)
}
\arguments{
\item{data}{}
}
\value{

}
\description{
scale_zscale
}
\examples{
mt <- t(matrix(c(rep(0, 5), rnorm(5), rep(2,5)), nrow = 5))
scale_zscale(mt)

}
