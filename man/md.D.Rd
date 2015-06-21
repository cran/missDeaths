\name{md.D}
\alias{md.D}
\title{
Correctly analyse disease recurrence with missing at risk information using population mortality
}
\description{
Utility function that returns a data.frame containing basic demographic information compatible with the \code{\link{md.survnp}}, \code{\link{md.survcox}} and \code{\link{md.impute}} functions. 
}
\usage{
md.D(age, sex, year)	
}
\arguments{
  \item{age}{age of the patients (number of days).}
  \item{sex}{contains 1 for males and 2 for females}
  \item{year}{year of entry into the study (number of days from 1/1/1960).}
}
\seealso{
  \code{\link{md.survcox}}, \code{\link{md.survnp}}
}