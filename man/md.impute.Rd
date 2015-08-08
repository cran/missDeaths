\name{md.impute}
\alias{md.impute}
\title{
Correctly analyse disease recurrence with missing at risk information using population mortality
}
\description{
The intuitive idea to `recover' the complete data set from the observed data set by imputing expected survival times obtained from the population mortality tables, does not give unbiased estimates, since the right censored indivuduals are not a random subsample of the patients. 

An iterative approach is used in this method to estimate the conditional distribution required to correctly impute the times of deaths.
}
\usage{
md.impute(data, f, maxtime, D, ratetable, iterations = 4)  
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
}
\value{
an array of times with imputed times of death 
}
\references{
Stupnik T., Pohar Perme M. (2015) "Analysing disease recurrence with missing at risk information." Statistics in Medicine.
    
Pohar Perme M., Stupnik T. (2014) "Analysing disease recurrence with missing at risk information." In: 35th Annual Conference of the International Society for Clinical Biostatistics, 24 - 28 August 2014, Vienna. Abstracts book, p82.
     
Stupnik T., Pohar Perme M. (2013). "Analysing disease recurrence with missing at risk information." In: International Conference Applied Statistics 2013, September 22-25, 2013, Ribno (Bled), Slovenia. Abstracts and program, p47.
}
\seealso{
  \code{\link{md.survcox}}
}
\examples{
library(missDeaths)
data(slopop)

data(observed)
observed$time = observed$time*365.2425
D = md.D(age=observed$age*365.2425, sex=observed$sex, year=(observed$year - 1960)*365.2425)
newtimes = md.impute(observed, Surv(time, status) ~ age + sex + iq + elevation, 
  observed$maxtime*365.2425, D, slopop, iterations=4)

#Cumulative incidence function
cif = survfit(Surv(observed$time, observed$status)~1)
cif$surv = 1 - cif$surv
cif$upper = 1 - cif$upper
cif$lower = 1 - cif$lower
plot(cif)

#Net survival (NOTE: std error is slightly underestimated!)
surv.net = survfit(Surv(newtimes, observed$status)~1)
summary(surv.net, times=c(3,9)*365.2425)
plot(surv.net)

#Event free survival (NOTE: std error is slightly underestimated!)
surv.efs = survfit(Surv(newtimes, 1 * (observed$status | (newtimes != observed$time)))~1)
summary(surv.efs, times=c(3,9)*365.2425)
plot(surv.efs)
}
