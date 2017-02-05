#' @useDynLib missDeaths
#' @importFrom Rcpp evalCpp

#export(rcpp_Sample)
#export(rcpp_SurvExpInit)
#export(rcpp_BlindDeath)
#export(rcpp_SurvTime)
#export(rcpp_SimCensorX)
#export(md.fixsample)

#' @importFrom mitools MIcombine
#' @importFrom stats qchisq runif rbinom rnorm rexp
#' @importFrom MASS mvrnorm
#' @importFrom methods new callNextMethod is
#' @import survival rms relsurv cmprsk
NULL


#' Simulating and analyzing time to event data in the presence of population mortality 
#' 
#' In analysis of time to event data we may have a situation where we know that a 
#' certain non-negligible competing risk exists, but is not recorded in the data. 
#' Due to competing nature of the risks, ignoring such a risk may significantly 
#' impact the at-risk group and thus lead to biased estimates. \cr\cr
#' This problem can be found in several national registries of benign diseases, 
#' medical device implantations (e.g. hip, knee or heart pacemaker) etc. where law obliges 
#' physicians to report events whereas the information on patient deaths is unavailable; 
#' it is hence unclear how many devices are still in use at a given time.\cr\cr
#' Under the assumption that the survival of an individual is not influenced by the 
#' event under study, general population mortality tables can be used to obtain unbiased 
#' estimates of the measures of interest or to verify the assumption that the bias
#' introduced by the non-recorded deaths is truly negligible.\cr\cr
#' Two approaches are implemented in the missdeaths package: \cr
#'  - an iterative imputation method \code{\link{md.survcox}} and\cr 
#'  - a mortality adjusted at risk function \code{\link{md.survnp}}.\cr\cr
#' The package also includes a comprehensive set of functions to simulate data 
#' that can be used for better understanding of these methods (See \code{\link{md.simulate}}).
#' 
#' @examples 
#' \dontrun{
#' library(missDeaths)
#' ratetable = survexp.us
#' 
#' sim = md.simparams() +
#'           md.sex("sex", 0.5) + 
#'           md.binom("Z1", 0.5) +
#'           md.mvnorm(c("Z2", "Z3"), c(100, 0), matrix(c(225, 3, 2, 1), 2, 2)) +
#'           md.sample("Z4", c(1, 2, 3, 4), c(0.25, 0.25, 0.25, 0.25)) +
#'           md.dateuniform("birth", as.Date("1925-1-1"), as.Date("1950-1-1")) +
#'           md.dateuniform("start", as.Date("2000-1-1"), as.Date("2005-1-1")) +
#'           md.death("death", ratetable, "sex", "birth", "start") +
#'           md.eval("age", "as.numeric(start - birth)/365.2425", 80, FALSE) + 
#'           md.exp("event", "start", c("age", "sex", "Z1", "Z2"), 
#'              c(0.1, 2, 1, 0.01), 0.05/365.2425)
#' data = md.simulate(sim, 1000)
#'           
#' #construct a complete right censored data set
#' complete = md.censor(data, "start", as.Date("2010-1-1"), c("event", "death"))
#' 
#' #construct an observed right censored data set with non-reported deaths
#' observed = complete
#' observed$time[which(complete$status == 2)] = observed$maxtime[which(complete$status == 2)]
#' observed$status[which(complete$status == 2)] = 0
#' 
#' #estimate coefficients of the proportional hazards model
#' D = md.D(age=observed$age, sex=observed$sex, year=observed$start)
#' md.survcox(observed, Surv(time, status) ~ age + sex + Z1 + Z2, 
#'           observed$maxtime, D, ratetable, iterations=4, R=50)
#'           
#' #estimate net- and event-free survival
#' np = md.survnp(observed$time, observed$status, observed$maxtime, D, ratetable)
#' w = list(list(time=np$time, est=np$surv.net, var=(np$std.err.net)^2))
#' timepoints(w, times=c(3,9)*365.2425)
#' plot(np$time/365.2425, np$surv.net)
#' plot(np$time/365.2425, np$surv.efs)
#' }
#' @name missDeaths
#' @docType package
#' @author Tomaz Stupnik <tomaz.stupnik@@gmail.com> and Maja Pohar Perme
#' 
#' @references Stupnik T., Pohar Perme M. (2015) "Analysing disease recurrence
#' with missing at risk information." Statistics in Medicine 35. p1130-43.
#' \url{http://onlinelibrary.wiley.com/doi/10.1002/sim.6766/abstract}
NULL





#' A simulated dataset with non-recorded deaths
#' 
#' This data set is used to illustrate the missDeaths functions.
#' 
#' 
#' @name observed
#' @docType data
#' @format A \code{data.frame} containing 10000 observations.
#' @keywords datasets
NULL
