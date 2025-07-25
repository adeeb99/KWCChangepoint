% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FKWC.R
\name{FKWC}
\alias{FKWC}
\title{Find changepoints using functional Kruskall-Wallis for covariance algorithm}
\usage{
FKWC(funcdata, depth = "FM_depth", beta = 10)
}
\arguments{
\item{funcdata}{Data in fdata form}

\item{depth}{Depth function of choice}

\item{beta}{Numeric penalty constant passed to PELT}
}
\value{
A list of changepoints
}
\description{
Find changepoints using functional Kruskall-Wallis for covariance algorithm
}
