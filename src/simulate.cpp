#include <Rcpp.h>
#include <string>
#include "survexpcache.h"

#define UNIFORM_AGE

using namespace Rcpp;
using namespace std;

class GenParams
{
protected:
  double GetParam(Rcpp::List list, char* name, double defValue)
  {
    try
    {
      Rcpp::NumericVector v = list[name];
      return (v.length() > 0) ? v[0] : defValue;
    }
    catch(...) {}
    return defValue;
  }
  bool GetParam(Rcpp::List list, char* name, bool defValue)
  {
    try
    {
      Rcpp::LogicalVector v = list[name];
      return (v.length() > 0) ? v[0] : defValue;
    }
    catch(...) {}
    return defValue;
  }
};

class AgeGenParams : GenParams
{
public:
  double MeanAge;
  double SdAge;
  double MinAge;
  double MaxAge;
  double Skew;
  
  AgeGenParams()
  {
    MeanAge = MinAge = MaxAge = 50;
    SdAge = 0;
    Skew = 0;
  }
  AgeGenParams(Rcpp::List list, char* name)
  {
    MeanAge = 50;
    MinAge = 0;
    MaxAge = 100;
    SdAge = 0; 
    Skew = 0;
    
    try
    {
      Rcpp::List sublist = list[name];
      MeanAge = GetParam(sublist, (char*)"mean", MeanAge);
      SdAge = GetParam(sublist, (char*)"sd", MeanAge);
      MinAge = GetParam(sublist, (char*)"min", MinAge);
      MaxAge = GetParam(sublist, (char*)"max", MaxAge);
      Skew = GetParam(sublist, (char*)"skew", Skew);
    }
    catch(...){}
  }
};

class SampleGenParams : GenParams
{
public:
  int     N;
  double  TimeInterval;
  double StartYear;
  double  MaleFraction;  
  bool StartTimeZero;
  
  AgeGenParams Treatment0Age;
  AgeGenParams Treatment1Age;
  
  double TreatmentOneFraction;
  int HaircutClasses;
  
  double BaselineLambda;
  double TreatmentOneHR;
  double BetaIQ;
  double BetaElevation;
  double BetaAge;
  double BetaSex;
  double BetaHaircut[100];
  
  SampleGenParams() : Treatment0Age(), Treatment1Age()
  {
    N = 1000;
    TimeInterval = 10;
    MaleFraction = 0.5;
  }
  SampleGenParams(Rcpp::List list) : Treatment0Age(list, (char*)"age0"), Treatment1Age(list, (char*)"age1")
  {
    N = GetParam(list, (char*)"N", (double)1000);
    TimeInterval = GetParam(list, (char*)"T", (double)10);
    StartYear = GetParam(list, (char*)"startyear", (double)1990);
    
    MaleFraction = GetParam(list, (char*)"maleperc", 0.5);
    TreatmentOneFraction = GetParam(list, (char*)"treatment1perc", 0.5);
    HaircutClasses = GetParam(list, (char*)"nhaircuts", (double)5);
    StartTimeZero = GetParam(list, (char*)"firstint0", true);
    
    BaselineLambda = GetParam(list, (char*)"lambda", (double)0.1);
    TreatmentOneHR = GetParam(list, (char*)"HR", exp(1));
    BetaIQ = GetParam(list, (char*)"betaiq", (double)0);
    BetaElevation = GetParam(list, (char*)"betaelevation", (double)0);
    BetaSex = GetParam(list, (char*)"betasex", (double)0);
    BetaAge = GetParam(list, (char*)"betaage", (double)0);
    
    for (int i = 0; i < 100; i++)
      BetaHaircut[i] = 0;
  }
};

Rcpp::NumericVector SimulateAge(Rcpp::NumericVector treatment, Rcpp::NumericVector sex, SampleGenParams params)
{
  Rcpp::Environment stats("package:stats");
  Rcpp::Function runif = stats["runif"]; 
  Rcpp::Function rnorm = stats["rnorm"]; 
  Rcpp::NumericVector age;
  
#ifdef UNIFORM_AGE
  age = runif((int)treatment.length(), -1, 1);
  for (int i = 0; i < treatment.length(); i++)
  {
    double p = (treatment[i] == 0) ? params.Treatment0Age.Skew : params.Treatment1Age.Skew;
    p = (sex[i] == 1) ? -p : p;
    p = max(-0.9, min(p, 0.9));

    if (age[i] < p)
      age[i] = -1 + (-1 - age[i])/(-1 - p);
    else
      age[i] = 1 - (1 - age[i])/(1 - p);
  }
#else
  age = rnorm((int)treatment.length());
#endif

  for (int i = 0; i < treatment.length(); i++)
    if (treatment[i] == 0)
    {
      age[i] = params.Treatment0Age.MeanAge + age[i] * params.Treatment0Age.SdAge;
      if (age[i] < params.Treatment0Age.MinAge)
        age[i] = (double)params.Treatment0Age.MinAge;
      if (age[i] > params.Treatment0Age.MaxAge)
        age[i] = (double)params.Treatment0Age.MaxAge;
    }
    else
    {
      age[i] = params.Treatment1Age.MeanAge + age[i] * params.Treatment1Age.SdAge;
      if (age[i] < params.Treatment1Age.MinAge)
        age[i] = (double)params.Treatment1Age.MinAge;
      if (age[i] > params.Treatment1Age.MaxAge)
        age[i] = (double)params.Treatment1Age.MaxAge;
    }

  return age;
}

// [[Rcpp::export]]
SEXP Sample(Rcpp::List paramsList)
{
  try
  {
    SampleGenParams params = SampleGenParams(paramsList);
     
    Rcpp::Environment stats("package:stats");
    Rcpp::Function runif = stats["runif"]; 
    Rcpp::Function rnorm = stats["rnorm"]; 
    Rcpp::Function rexp = stats["rexp"]; 
           
    Rcpp::NumericVector sex = runif(params.N);
    for (int i = 0; i < params.N; i++)
      sex[i] = (sex[i] < params.MaleFraction) ? 1 : 2;
  
    Rcpp::NumericVector firstint(params.N, (double)0);
    if (!params.StartTimeZero)
      firstint = runif(params.N, 0, params.TimeInterval / 2);
 
    Rcpp::NumericVector treatment = runif(params.N);
    for (int i = 0; i < params.N; i++)
      treatment[i] = (treatment[i] > params.TreatmentOneFraction) ? 0 : 1;
    Rcpp::NumericVector haircut = runif(params.N);
    for (int i = 0; i < params.N; i++)
      haircut[i] = 1 + floor(haircut[i] * params.HaircutClasses);
    Rcpp::NumericVector iq = rnorm(params.N);
    for (int i = 0; i < params.N; i++)
      iq[i] = 100 + iq[i] * 15;
    Rcpp::NumericVector elevation = runif(params.N);
    Rcpp::NumericVector dummy = runif(params.N);
  
    Rcpp::NumericVector age = SimulateAge(treatment, sex, params);
    Rcpp::IntegerVector birthYear(params.N);  
    double meanAge = 0;
    for (int i = 0; i < params.N; i++)
    {
      birthYear[i] = round(params.StartYear + firstint[i] - age[i]);
      meanAge += age[i];
    }
    meanAge /= (double)params.N;
      
    Rcpp::NumericVector deathint = runif(params.N);
    for (int i = 0; i < params.N; i++)
    {
      double t = SurvTime(birthYear[i], age[i] * 365.2425, deathint[i], sex[i]);
      deathint[i] = t / 365.2425;
    }
    
    Rcpp::NumericVector lambda(params.N, params.BaselineLambda);
    for (int i= 0; i < params.N; i++)
    {
      if (treatment[i] > 0)
        lambda[i] = lambda[i] * params.TreatmentOneHR;
  
      double beta = 0;
      beta += params.BetaHaircut[(int)floor(haircut[i]) - 1];
      beta += params.BetaIQ * (iq[i] - 100);
      beta += params.BetaElevation * (elevation[i] - 0.5);
      beta += params.BetaAge * (age[i] - meanAge);
      beta += params.BetaSex * (sex[i] - (2.000 - params.MaleFraction));
    
      lambda[i] = lambda[i] * exp(beta);
    }
    Rcpp::NumericVector recurrint = rexp(params.N, lambda);
    
    Rcpp::NumericVector maxTime(params.N, (double)0);
    for (int i= 0; i < params.N; i++)
      maxTime[i] = params.TimeInterval - firstint[i];
  
    Rcpp::DataFrame data = Rcpp::DataFrame::create(Rcpp::Named("age") = age, Rcpp::Named("sex") = sex, Rcpp::Named("year") = birthYear, Rcpp::Named("event") = recurrint, Rcpp::Named("death") = deathint, Rcpp::Named("maxtime") = maxTime, Rcpp::Named("treatment") = treatment, Rcpp::Named("iq") = iq, Rcpp::Named("elevation") = elevation, Rcpp::Named("haircut") = haircut, Rcpp::Named("dummy") = dummy);
    return Rcpp::wrap(data);
  }
  catch (...){}
  throw std::range_error("Unknown Sample() Error");
}

// [[Rcpp::export]]
SEXP Resample(Rcpp::DataFrame data1)
{
  try
  {
    Rcpp::DataFrame data = Rcpp::clone(data1);
    
    Rcpp::NumericVector age = data["age"];
    Rcpp::NumericVector sex = data["sex"];
    Rcpp::IntegerVector birthYear = data["year"];
    Rcpp::NumericVector newdeathint = runif(data.nrows());
    for (int i = 0; i < data.nrows(); i++)
    {
      double t = SurvTime(birthYear[i], age[i] * 365.2425, newdeathint[i], sex[i]);
      newdeathint[i] = t / 365.2425;
    }    
    
    Rcpp::List myList;
    Rcpp::CharacterVector namevec = data.names();
    for (int i = 0; i < data.length(); i++) 
    {
      char *name = namevec[i];
      if (strcmp(name, "death") != 0)
        myList.push_back(data(i));   
      else
        myList.push_back(newdeathint);
    }
    myList.attr("names") = namevec;
    Rcpp::DataFrame resampleddata(myList);   
    
    return Rcpp::wrap(resampleddata);
  }
  catch(...){}
  throw std::range_error("Unknown Resample() Error");
}

// [[Rcpp::export]]
SEXP BlindDeath(Rcpp::DataFrame data1, bool zombies)//, double interval)
{
  try
  {
    Rcpp::DataFrame data = Rcpp::clone(data1);
    
    Rcpp::List myList;
    Rcpp::CharacterVector oldnamevec = data.names();
    Rcpp::CharacterVector namevec;
    for (int i = 0; i < data.length(); i++) 
    {
      char *name = oldnamevec[i];
      if (strcmp(name, "death") != 0)
      {
        myList.push_back(data(i));      
        if (strcmp(name, "event") == 0)
          name = (char*)"time";
        namevec.push_back(name);
      }
    }
    Rcpp::IntegerVector status(data.nrows());
    myList.push_back(status);
    namevec.push_back("status");
    
    myList.attr("names") = namevec;
    Rcpp::DataFrame blinddata(myList);
    
    status = blinddata["status"];
    Rcpp::NumericVector recurrence = blinddata["time"];
    Rcpp::NumericVector maxTime = blinddata["maxtime"];
    
    Rcpp::NumericVector recurrint = data["event"];
    Rcpp::NumericVector deathint = data["death"];
    for (int i = 0; i < data.nrows(); i++)
    {
      recurrence[i] = recurrint[i];
      status[i] = 1;
      
      if (recurrence[i] > maxTime[i])    //right censoring, end of interval
      {
        recurrence[i] = maxTime[i];
        status[i] = 0;
      }
      if (recurrence[i] > deathint[i])                 //censor out zombies
      {	
        if (!zombies)
          recurrence[i] = deathint[i];
        else
          recurrence[i] = maxTime[i];
        status[i] = 0;
      }
    }
    return Rcpp::wrap(blinddata);
  }
  catch(...){}
  throw std::range_error("Unknown BlindDeath() Error");
}
