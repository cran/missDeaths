\name{md.survcox}
\alias{md.survcox}
\title{
Correctly analyse disease recurrence with missing at risk information using population mortality
}
\description{
The intuitive idea to `recover' the complete data set from the observed data set by imputing expected survival times obtained from the population mortality tables, does not give unbiased estimates, since the right censored indivuduals are not a random subsample of the patients. 

An iterative approach is used in this method to estimate the conditional distribution required to correctly impute the times of deaths.
}
\usage{
md.survcox(data, f, maxtime, D, ratetable, iterations = 4, R = 50)	
}
\arguments{
  \item{data}{a data.frame in which to interpret the variables named in the formula.}
  \item{f}{a formula object, with the response on the left of a ~ operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}
  \item{maxtime}{maximum potential observation time (number of days). 
  
  where \code{status}=0 equals \code{time}. 
  
  where \code{status}=1 equals potential time of right censoring if no event would be observed.
  }
  \item{D}{demographic information compatible with \code{ratetable}, see \code{\link{md.D}}.}
  \item{ratetable}{a population mortality table, default is \code{slopop}}
  \item{iterations}{the number of iteration steps to be performed, default is 4}
  \item{R}{the number of multiple imputations performed to adjust the estimated variance of estimates, default is 50.}
}
\value{
if \code{R} equals 1 then an object of class \code{\link[survival]{coxph.object}} representing the fit. 

if \code{R} > 1 then the result of the \code{\link[mitools]{MIcombine}} of the \code{coxph} objects.
}
\references{
Stupnik T., Pohar Perme M. (2015) "Analysing disease recurrence with missing at risk information." Statistics in Medicine.
    
Pohar Perme M., Stupnik T. (2014) "Analysing disease recurrence with missing at risk information." In: 35th Annual Conference of the International Society for Clinical Biostatistics, 24 - 28 August 2014, Vienna. Abstracts book, p82.
     
Stupnik T., Pohar Perme M. (2013). "Analysing disease recurrence with missing at risk information." In: International Conference Applied Statistics 2013, September 22-25, 2013, Ribno (Bled), Slovenia. Abstracts and program, p47.
}
\seealso{
  \code{\link{md.impute}}, \code{\link[mitools]{MIcombine}}
}
\examples{
\dontrun{
library(missDeaths)
data(slopop)

data(observed)
observed$time = observed$time*365.2425
D = md.D(age=observed$age*365.2425, sex=observed$sex, year=(observed$year - 1960)*365.2425)

#fit a cox model (NOTE: estimated std error is slightly underestimated!)
cox = md.survcox(observed, Surv(time, status) ~ age + sex + iq + elevation, 
  observed$maxtime*365.2425, D, slopop, iterations=4, R=1)

#multiple imputations to correct the stimated std error
mi = md.survcox(observed, Surv(time, status) ~ age + sex + iq + elevation, 
  observed$maxtime*365.2425, D, slopop, iterations=4, R=50)
}
}
