% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics_s3.R
\name{req}
\alias{req}
\title{Generic method for ensuring object contains required elements}
\usage{
req(m, r)
}
\arguments{
\item{m}{object}

\item{r}{vector of elements to test for}
}
\value{
throws error if required element is not present
}
\description{
Intended to be called internally by package
}
