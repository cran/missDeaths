#include <Rcpp.h>
#include <string>
#include "survexpcache.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
SEXP exSimCensor1(Rcpp::DataFrame data1)//, double interval)
{
  try
  {
    Rcpp::DataFrame data = Rcpp::clone(data1);
   
    Rcpp::NumericVector status = Rcpp::wrap(data["status"]);
    Rcpp::NumericVector recurrence = Rcpp::wrap(data["time"]);
    Rcpp::NumericVector age = Rcpp::wrap(data["age"]);
    Rcpp::NumericVector sex = Rcpp::wrap(data["sex"]);
    Rcpp::NumericVector year = Rcpp::wrap(data["year"]);
    
    for (int i = 0; i < data.nrows(); i++)
    {
      if (status[i] == 0)
      {
        double t = SurvTime(year[i], age[i] * 365.2425, 0.5, sex[i]);
        if (recurrence[i] > (t / 365.2425))
          recurrence[i] = t / 365.2425;
      }
    }  
    return Rcpp::wrap(data);
  }
  catch(...){}
  throw std::range_error("Unknown SimCensor1() Error");
}

// [[Rcpp::export]]
SEXP exSimCensor2(Rcpp::DataFrame data1)//, double interval)
{
  try
  {
    Rcpp::DataFrame data = Rcpp::clone(data1);
    Rcpp::Environment stats("package:stats");
    Rcpp::Function runif = stats["runif"]; 
    
    Rcpp::NumericVector badluck = runif(data.nrows());
    Rcpp::NumericVector status = Rcpp::wrap(data["status"]);
    Rcpp::NumericVector recurrence = Rcpp::wrap(data["time"]);
    Rcpp::NumericVector age = Rcpp::wrap(data["age"]);
    Rcpp::NumericVector sex = Rcpp::wrap(data["sex"]);
    Rcpp::NumericVector year = Rcpp::wrap(data["year"]);
    
    for (int i = 0; i < data.nrows(); i++)
    {
      if (status[i] == 0)
      {
        double t = SurvTime(year[i], age[i] * 365.2425, badluck[i], sex[i]);
        if (recurrence[i] > (t / 365.2425))
          recurrence[i] = t / 365.2425;
      }
    }  
    return Rcpp::wrap(data);
  }
  catch(...){}
  throw std::range_error("Unknown SimCensor2() Error");
}

// [[Rcpp::export]]
SEXP exSimCensorX(Rcpp::DataFrame data1, /*double interval,*/ SEXP form, int maxiter)
{
  // try
  {
    Rcpp::Environment rms("package:rms");
    Rcpp::Environment stats("package:stats");
    Rcpp::Function runif = stats["runif"]; 
    Rcpp::Function survest = rms["survest"]; 
    Rcpp::Function cph = rms["cph"]; 
    
    Rcpp::DataFrame data = exSimCensor2(data1);
    Rcpp::Formula formula(form);
    for (int l = 0; l < maxiter; l++)
    {  
      bool f = FALSE;
      bool t = TRUE;
      Rcpp::List cox = cph(formula, data, Rcpp::Named("surv") = t);
      
      Rcpp::NumericVector coxtime = cox["time"];  
      Rcpp::NumericVector times(min(100, coxtime.length()));
      for (int j = 0; j < times.length(); j++)
        times[j] = coxtime[round(((double)coxtime.length() - 1) * j / (times.length() - 1))];
      
      data = Rcpp::clone(data1);
      Rcpp::NumericVector status = data["status"];
      Rcpp::NumericVector recurrence = data["time"];
      Rcpp::NumericVector age = data["age"];
      Rcpp::NumericVector sex = data["sex"];
      Rcpp::NumericVector year = data["year"];
      Rcpp::NumericVector maxTime = data["maxtime"];
      
      Rcpp::List fit = survest(cox, data, Rcpp::Named("times") = times, Rcpp::Named("se.fit") = f, Rcpp::Named("conf.int") = f);
      Rcpp::NumericMatrix est = fit["surv"];
      
      for (int i = 0; i < data.nrows(); i++)
      if (status[i] == 0)
      {
        Rcpp::NumericVector surv = est(i,_);
        while(1)
        {
          Rcpp::NumericVector badluck = runif(2);            
          double t = SurvTime(year[i], age[i] * 365.2425, badluck[1], sex[i]);  
          
          double mdeath = t / 365.2425;
          double time = 0;
          double TT = maxTime[i];
          
          if (badluck[0] < surv[surv.length() - 1])
          {
            time = TT + 1;
          } 
          else 
          {
            for (int r = 0; r < surv.length(); r++)
            if (surv[r] < badluck[0])
            {
              double t1, t2, s1, s2;
              if (r > 0)
              {
                t1 = times[r - 1];
                t2 = times[r];
                s1 = surv[r - 1];
                s2 = surv[r];
              }
              else
              {
                t1 = 0;
                t2 = times[0];
                s1 = 1;
                s2 = surv[0];
              }
              time = t1 + (t2 - t1)*(s1 - badluck[0])/(s1 - s2);
              break;
            }
          }
          
          if ((time > TT) && (TT < mdeath)) break;
          if (min(time, TT) > mdeath)
          {
            recurrence[i] = mdeath;
            break;
          }
        }
      }
      //Rcpp::List l = Rcpp::List::create(Rcpp::Named("times") = times, Rcpp::Named("cox") = cox, Rcpp::Named("data") = data, Rcpp::Named("est") = est);
      //return Rcpp::wrap(l);
    }
    return Rcpp::wrap(data);
  }
  //  catch(...){}
  throw std::range_error("Unknown SimCensorX() Error");
}
