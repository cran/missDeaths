# @useDynLib missDeaths 


#' Initializes the missDeaths population table cache
#' 
#' @keywords internal 
#' 
#' 
#' @param poptable 
#' @export md.init
md.init <- function(poptable) {
  SurvExpInit(poptable)
}



#' Calculates individual expected survival beyond the observation time
#' 
#' @keywords internal
#' 
#' @param D Demographic information
#' @param time Observation time
#' @export md.expprep
md.expprep <- function(D, time) {
  SurvExpPrep(D, time)
}


#' Initializes the missDeaths population table cache
#' 
#' @keywords internal 
#' 
#' 
#' @param year year
#' @param age  age
#' @param prob prob
#' @param sex sex
#' @export md.survtime
md.survtime <- function(year, age, prob, sex) {
  SurvTime(year, age, prob, sex)
}


#' Initializes the missDeaths population table cache
#' 
#' @keywords internal 
#' 
#' 
#' @param year year
#' @param age  age
#' @param time time
#' @param sex sex
#' @export md.survprob
md.survprob <- function(year, age, time, sex) {
  SurvProbability(year, age, time, sex)
}


#' Initializes the missDeaths population table cache
#' 
#' @keywords internal 
#' 
#' 
#' @param year year
#' @param sex sex
#' @export md.survdump
#md.survdump <- function(year, sex) {
#  SurvDump(year, sex)
#}



#' Prepare compatible demographic information
#' 
#' Utility function that returns a data.frame containing basic demographic
#' information compatible with the \code{\link{md.survnp}},
#' \code{\link{md.survcox}} and \code{\link{md.impute}} functions.
#' 
#' 
#' @param age vector of patient ages specified as number of days or number of years.
#' @param sex vector containing 1 for males and 2 for females
#' @param year vector of years of entry into the study can either be supplied 
#' as vector of start dates or as vector of years specified in number of days from origin (1-1-1970).
#' @seealso \code{\link{md.survcox}}, \code{\link{md.survnp}}
#' @export md.D
md.D <- function(age, sex, year)
{
  if (max(age) < 150)
    age = age * 365.2425
  
  if (!is.numeric(year))
    year = as.numeric(year - as.Date("0-1-1")) - 1970*365.2425
    
  D = data.frame(age=age, sex=sex, year=year)
  
  sexlevels = levels(as.factor(D$sex))
  for (i in 1:length(sexlevels))
    if (!(sexlevels[i] %in% c("1", "2")))
      stop ("column 'sex' can only contain values 1 and 2")
  
  if ((min(D$age) < 0) || (max(D$age) > 130 * 365.2425))
    stop ("values in column 'age' out of bounds, should be within [0, 130 * 365] days")
  
  return (D)
}

md.fixsample <- function(observed)
{
  md = observed
  md$year = (md$year + md$age - 1970) * 365.2425
  md$age = md$age * 365.2425
  md$time = round(md$time * 365.2425)
  md$maxtime = round(md$maxtime * 365.2425)
  return (md)
}



#' Correctly impute missing information of possible deaths using population mortality
#' 
#' An iterative approach is used in this method to estimate the conditional
#' distribution required to correctly impute the times of deaths using
#' population mortality tables.\cr\cr
#'   Note, that simply imputing expected survival times may seem intuitive, 
#' but does not give unbiased estimates, since the right censored individuals 
#' are not a random subsample of the patients.
#' 
#' 
#' @param data a data.frame in which to interpret the variables named in the
#' formula.
#' @param f a formula object, with the response on the left of a ~ operator,
#' and the terms on the right. The response must be a survival object as
#' returned by the \code{Surv} function.
#' @param maxtime maximum potential observation time (number of days).
#' 
#' where \code{status}=0 equals \code{time}.
#' 
#' where \code{status}=1 equals potential time of right censoring if no event
#' would be observed.
#' @param D demographic information compatible with \code{md.survcox}, \code{md.impute} 
#' and \code{md.survnp}, see \code{\link{md.D}}.
#' @param ratetable a population mortality table, default is \code{slopop}
#' @param iterations the number of iteration steps to be performed, default is
#' 4
#' @return an array of times with imputed times of death that can be used instead of the
#' unavailable complete data set to get unbiased estimates, ie. in \code{\link[survival]{coxph}}. 
#' @seealso \code{\link{md.survcox}}
#' @references Stupnik T., Pohar Perme M. (2015) "Analysing disease recurrence
#' with missing at risk information." Statistics in Medicine 35. p1130-43.
#' \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6766}

#' @examples
#' library(missDeaths)
#' data(slopop)
#' 
#' data(observed)
#' observed$time = observed$time*365.2425
#' D = md.D(age=observed$age*365.2425, sex=observed$sex, year=(observed$year - 1970)*365.2425)
#' newtimes = md.impute(observed, Surv(time, status) ~ age + sex + iq + elevation, 
#'   observed$maxtime*365.2425, D, slopop, iterations=4)
#' 
#' #Cumulative incidence function
#' cif = survfit(Surv(observed$time, observed$status)~1)
#' cif$surv = 1 - cif$surv
#' cif$upper = 1 - cif$upper
#' cif$lower = 1 - cif$lower
#' plot(cif)
#' 
#' #Net survival (NOTE: std error is slightly underestimated!)
#' surv.net = survfit(Surv(newtimes, observed$status)~1)
#' summary(surv.net, times=c(3,9)*365.2425)
#' plot(surv.net)
#' 
#' #Event free survival (NOTE: std error is slightly underestimated!)
#' surv.efs = survfit(Surv(newtimes, 1 * (observed$status | (newtimes != observed$time)))~1)
#' summary(surv.efs, times=c(3,9)*365.2425)
#' plot(surv.efs)
#' 
#' @export md.impute
md.impute <- function(data, f, maxtime, D, ratetable, iterations=4)
{  
  #require(survival)
  #require(rms)
  
  D$year = 1970 + round((D$year - D$age) / 365.2425)
  f = deparse(f)
  f = gsub(" ", "", f, fixed = TRUE)

  md.init(ratetable)
  return (SimCensorX(data, maxtime, f, D, iterations))
}



#' Nonparametric analysis of disease recurrence with missing information of possible deaths
#' 
#' Estimates the Net and Event free survial using a is non-parametric approach
#' that aims to correct all individuals using the unconditional survival time
#' distribution obtained from the population mortality table.\cr\cr
#'   The idea comes from realizing that the number of observed events in the data
#' equals the number which would be observed in case of a complete data set,
#' but the number of patients at risk does not. Hence, this method adjusts the
#' observed number at risk to mimic the one we would get if the data was
#' complete.
#' 
#' 
#' @param time the time to event (number of days)
#' @param status the status indicator, 0=right censored, 1=event at \code{time}
#' @param maxtime maximum potential observation time (number of days).
#' 
#' where \code{status}=0 equals \code{time}.
#' 
#' where \code{status}=1 equals potential time of right censoring if no event
#' would be observed.
#' @param D demographic information compatible with \code{ratetable}, see
#' \code{\link{md.D}}.
#' @param ratetable a population mortality table, default is \code{slopop}
#' @param conf.int desired coverage of the estimated confidence interval
#' @return A list with components giving the estimates of net and event free
#' survival.
#' 
#' \item{time}{times where the estimates are calculated (number of days)}
#' \item{Y.net}{adjusted number of patients at risk at each time in a hypothetical world where patients don't die}
#' \item{Y.efs}{adjusted number of patients at risk at each time}
#' \item{surv.net}{the estimated Net survival} \item{std.err.net}{the estimated
#' standard error of Net survival estimates} \item{surv.efs}{the estimated
#' Event free survival} \item{std.err.efs}{the estimated standard error of
#' Event free survival estimates}
#' @references Stupnik T., Pohar Perme M. (2015) "Analysing disease recurrence
#' with missing at risk information." Statistics in Medicine 35. p1130-43.
#' \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6766}
#' @examples
#' 
#' \dontrun{
#' library(missDeaths)
#' library(cmprsk)
#' data(slopop)
#' 
#' data(observed)
#' D = md.D(age=observed$age*365.2425, sex=observed$sex, year=(observed$year - 1970)*365.2425)
#' np = md.survnp(observed$time*365.2425, observed$status, observed$maxtime*365.2425, D, slopop)
#' 
#' #calculate net survival at 3 and 9 years
#' w = list(list(time=np$time, est=np$surv.net, var=(np$std.err.net)^2))
#' timepoints(w, times=c(3,9)*365.2425)
#' 
#' #plot the net and event free survival curves
#' plot(np$time, np$surv.net)
#' plot(np$time, np$surv.efs)
#' }
#' 
#' @export md.survnp
md.survnp <- function(time, status, maxtime, D, ratetable, conf.int=0.95) 
{
  return (my.survnp(time, status, D, ratetable, maxtime))
}


#' Fit a proportional hazards regression model over disease recurrence data 
#' with missing information of possible deaths
#' 
#' An iterative approach is used in this method to estimate the conditional
#' distribution required to correctly impute the times of deaths using
#' population mortality tables.\cr\cr
#'   Note, that simply imputing expected survival times may seem intuitive, 
#' but does not give unbiased estimates, since the right censored individuals 
#' are not a random subsample of the patients.
#' 
#' @param data a data.frame in which to interpret the variables named in the
#' formula.
#' @param f a formula object, with the response on the left of a ~ operator,
#' and the terms on the right. The response must be a survival object as
#' returned by the \code{Surv} function.
#' @param maxtime maximum potential observation time (number of days).
#' 
#' where \code{status}=0 equals \code{time}.
#' 
#' where \code{status}=1 equals potential time of right censoring if no event
#' would be observed.
#' @param D demographic information compatible with \code{ratetable}, see
#' \code{\link{md.D}}.
#' @param ratetable a population mortality table, default is \code{slopop}
#' @param iterations the number of iteration steps to be performed, default is
#' 4
#' @param R the number of multiple imputations performed to adjust the
#' estimated variance of estimates, default is 50.
#' @return if \code{R} equals 1 then an object of class
#' \code{\link[survival]{coxph.object}} representing the fit.
#' 
#' if \code{R} > 1 then the result of the \code{\link[mitools]{MIcombine}} of
#' the \code{coxph} objects.
#' @seealso \code{\link{md.impute}}, \code{\link[mitools]{MIcombine}}
#' @references Stupnik T., Pohar Perme M. (2015) "Analysing disease recurrence
#' with missing at risk information." Statistics in Medicine 35. p1130-43.
#' \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.6766}
#' @examples
#' 
#' \dontrun{
#' library(missDeaths)
#' data(slopop)
#' 
#' data(observed)
#' observed$time = observed$time*365.2425
#' D = md.D(age=observed$age*365.2425, sex=observed$sex, year=(observed$year - 1970)*365.2425)
#' 
#' #fit a cox model (NOTE: estimated std error is slightly underestimated!)
#' md.survcox(observed, Surv(time, status) ~ age + sex + iq + elevation, 
#'   observed$maxtime*365.2425, D, slopop, iterations=4, R=1)
#' 
#' #multiple imputations to correct the stimated std error
#' md.survcox(observed, Surv(time, status) ~ age + sex + iq + elevation, 
#'   observed$maxtime*365.2425, D, slopop, iterations=4, R=50)
#' }
#' 
#' @export md.survcox
md.survcox <- function(data, f, maxtime, D, ratetable, iterations=4, R = 50)
{  
  D$year = 1970 + round((D$year - D$age) / 365.2425)
  ff = deparse(f)
  ff = gsub(" ", "", ff, fixed = TRUE)
    
  md.init(ratetable)   
  
  if (R == 1)
  {
    newdata = data
    newdata$time = SimCensorX(data, maxtime, ff, D, iterations)
    cox = survival::coxph(f, data=newdata)
    
    return (cox)
  }

  models = list()
  length(models) = R
  for (i in 1:R)
  {
    newdata = data
    newdata$time = SimCensorX(data, maxtime, ff, D, iterations)
    models[[i]] = survival::coxph(f, data=newdata)
  }
  
  return (MIcombine(models))
}
