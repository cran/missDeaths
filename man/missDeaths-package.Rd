\name{missDeaths-package}
\alias{missDeaths-package}
\alias{missDeaths}
\docType{package}
\title{
Correctly analyse disease recurrence with missing at risk information using population mortality
}
\description{
When analysing time to disease recurrence, we sometimes need to deal with data where all the recurrences are recorded, but no information is available on the possible deaths. This may happen when studying diseases of benign nature where patients are only seen at recurrences or in poorly designed national registries with insufficient patient identifiers to obtain their dead/alive status. 

When the average time to disease recurrence is long enough in comparison to the expected survival of the patients, the statistical analysis of such data may be significantly biased.

Under the assumption that the expected survival of an individual is not influenced by the disease itself, general population mortality tables may be used to remove this bias.  We show why the intuitive solution of simply imputing the patient's expected survival time does not give unbiased estimates of the usual quantities of interest in survival analysis and explain that the cumulative incidence function analysis does not require additional assumptions on  general population mortality.

Two new approaches are implemented in the missdeaths package: an iterative imputation method  \code{\link{md.survcox}} and a mortality adjusted at risk function  \code{\link{md.survnp}}.
}
\details{
\tabular{ll}{
Package: \tab missDeaths\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-01\cr
License: \tab GPL (>=2)\cr
}
}
\author{
Tomaz Stupnik and Maja Pohar Perme

Maintainer: Tomaz Stupnik <tomaz.stupnik@gmail.com>
}
\references{
Stupnik T., Pohar Perme M. (2015) "Analysing disease recurrence with missing at risk information." Statistics in Medicine.
    
Pohar Perme M., Stupnik T. (2014) "Analysing disease recurrence with missing at risk information." In: 35th Annual Conference of the International Society for Clinical Biostatistics, 24 - 28 August 2014, Vienna. Abstracts book, p82.
     
Stupnik T., Pohar Perme M. (2013). "Analysing disease recurrence with missing at risk information." In: International Conference Applied Statistics 2013, September 22-25, 2013, Ribno (Bled), Slovenia. Abstracts and program, p47.
}
\keyword{ missdeaths }
