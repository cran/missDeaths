\name{md.survnp}
\alias{md.survnp}
\title{
Correctly analyse disease recurrence with missing at risk information using population mortality
}
\description{
Estimates the Net and Event free survial using a is non-parametric approach that aims to correct all individuals using the unconditional survival time distribution obtained from the population mortality table. 

The idea comes from realizing that the number of observed events in the data equals the number which would be observed in case of a complete data set, but the number of patients at risk does not. Hence, this method adjusts the observed number at risk to mimic the one we would get if the data was complete.
}
\usage{
md.survnp(time, status, maxtime, D, ratetable, conf.int = 0.95)
}
\arguments{
  \item{time}{the time to event (number of days)
  }
  \item{status}{the status indicator, 0=right censored, 1=event at \code{time}
  }
  \item{maxtime}{maximum potential observation time (number of days). 
  
  where \code{status}=0 equals \code{time}. 
  
  where \code{status}=1 equals potential time of right censoring if no event would be observed.
  }
  \item{D}{demographic information compatible with \code{ratetable}, see \code{\link{md.D}}.}
  \item{ratetable}{a population mortality table, default is \code{slopop}}
  \item{conf.int}{desired coverage of the estimated confidence interval}
}
\value{
A list with components giving the estimates of net and event free survival.

  \item{time}{the times where the estimates are calculated (number of days)}
  \item{surv.net}{the estimated Net survival}
  \item{std.err.net}{the estimated standard error of Net survival estimates}
  \item{surv.efs}{the estimated Event free survival}
  \item{std.err.efs}{the estimated standard error of Event free survival estimates}
}
\examples{
\dontrun{
library(missDeaths)
library(cmprsk)
data(slopop)

data(observed)
D = md.D(age=observed$age*365.2425, sex=observed$sex, year=(observed$year - 1960)*365.2425)
np = md.survnp(observed$time*365.2425, observed$status, observed$maxtime*365.2425, D, slopop)

#calculate net survival at 3 and 9 years
w = list(list(time=np$time, est=np$surv.net, var=(np$std.err.net)^2))
timepoints(w, times=c(3,9)*365.2425)

#plot the net and event free survival curves
plot(np$time, np$surv.net)
plot(np$time, np$surv.efs)
}
}
