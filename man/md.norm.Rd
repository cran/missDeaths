% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/library.R
\name{md.norm}
\alias{md.norm}
\title{md.norm}
\usage{
md.norm(name, mean = 0, sd = 1)
}
\arguments{
\item{name}{name of the covariate}

\item{mean, sd}{mean and standard deviation}
}
\description{
Creates information of a normally distributed numeric covariate with the specified mean and standard deviation.
This function call must be added to the \code{\link{md.simparams}} object.
}
\examples{

\dontrun{
library(missDeaths)

sim = md.simparams() +
   md.norm("X", 0, 1) 
}

}
