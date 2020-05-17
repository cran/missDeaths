#include <Rcpp.h>
#include <string>
#include "survexpcache.h"

using namespace Rcpp;
using namespace std;

SEXP SimCensor1(Rcpp::NumericVector time1, Rcpp::IntegerVector status1, Rcpp::DataFrame D1)
{
  if (SurvExpCache != NULL)
    throw std::range_error("SurvExpCache is NULL");
  
  try
  {
    Rcpp::NumericVector recurrence = Rcpp::clone(time1);
    Rcpp::IntegerVector status = Rcpp::clone(status1);
    Rcpp::DataFrame D = Rcpp::clone(D1);
 
    Rcpp::NumericVector age = Rcpp::wrap(D["age"]);
    Rcpp::NumericVector sex = Rcpp::wrap(D["sex"]);
    Rcpp::NumericVector year = Rcpp::wrap(D["year"]);
    
    for (int i = 0; i < D.nrows(); i++)
    {
      if (status[i] == 0)
      {
        double t = SurvTime(year[i], age[i], 0.5, sex[i]);
        if (recurrence[i] > t)
          recurrence[i] = t;
      }
    }  
    return Rcpp::wrap(recurrence);
  }
  catch(...){}
  throw std::range_error("Unknown SimCensor1() Error");
}

SEXP SimCensor2(Rcpp::NumericVector time1, Rcpp::IntegerVector status1, Rcpp::DataFrame D1)
{
  if (SurvExpCache != NULL)
    throw std::range_error("SurvExpCache is NULL");
  
  try
  {
    Rcpp::NumericVector recurrence = Rcpp::clone(time1);
    Rcpp::IntegerVector status = Rcpp::clone(status1);
    Rcpp::DataFrame D = Rcpp::clone(D1);
    
    Rcpp::Environment stats("package:stats");
    Rcpp::Function runif = stats["runif"]; 
    
    Rcpp::NumericVector badluck = runif(D1.nrows());
    Rcpp::NumericVector age = Rcpp::wrap(D["age"]);
    Rcpp::NumericVector sex = Rcpp::wrap(D["sex"]);
    Rcpp::NumericVector year = Rcpp::wrap(D["year"]);
    
    for (int i = 0; i < D.nrows(); i++)
    {
      if (status[i] == 0)
      {
        double t = SurvTime(year[i], age[i], badluck[i], sex[i]);
        if (recurrence[i] > t)
          recurrence[i] = t;
      }
    }
    return Rcpp::wrap(recurrence);
  }
  catch(...){}
  throw std::range_error("Unknown SimCensor2() Error");
}

// [[Rcpp::export]]
SEXP SimCensorX(Rcpp::DataFrame data1, Rcpp::NumericVector maxtime1, Rcpp::CharacterVector form1, Rcpp::DataFrame D1, int maxiter)
{ 
  if (SurvExpCache != NULL)
    throw std::range_error("SurvExpCache is NULL");
  
  //try
  {
    Rcpp::Environment rms("package:rms");
    Rcpp::Environment stats("package:stats");
    Rcpp::Function runif = stats["runif"]; 
    Rcpp::Function survest = rms["survest"]; 
    Rcpp::Function cph = rms["cph"]; 
    
    std::string colTime, colStatus;
    std::string text = Rcpp::as<std::string>(form1[0]);
    std::transform(text.begin(), text.end(), text.begin(), ::tolower);
    size_t start = text.find("surv(");
    size_t end = text.find(")~");
    if ((start != string::npos) && (end != string::npos))
    {
      std::string str = text.substr(start + 5, end - start - 5);
      size_t comma = str.find(',');
      if (comma != string::npos)
      {
        colTime = str.substr(0, comma);
        colStatus = str.substr(comma + 1, str.length() - comma);
      }
    }
    
    if (colTime.empty() || colStatus.empty())
      throw std::range_error("");
    
    Rcpp::DataFrame D = Rcpp::clone(D1);
    Rcpp::NumericVector age = Rcpp::wrap(D["age"]);
    Rcpp::NumericVector sex = Rcpp::wrap(D["sex"]);
    Rcpp::NumericVector year = Rcpp::wrap(D["year"]);
   
    Rcpp::DataFrame data = Rcpp::clone(data1);
    Rcpp::NumericVector time = Rcpp::wrap(data[colTime]);
    Rcpp::IntegerVector status = Rcpp::wrap(data[colStatus]);
    Rcpp::NumericVector recurrence = SimCensor2(time, status, D); 
    Rcpp::NumericVector maxTime = Rcpp::clone(maxtime1);
    
    Rcpp::Formula formula = Rcpp::Formula(form1);
    for (int l = 0; l < maxiter; l++)
    {
      for (int j = 0; j < time.length(); j++)
        time[j] = recurrence[j];     
       
      bool f = FALSE;
      bool t = TRUE;
      Rcpp::List cox = cph(formula, data, Rcpp::Named("surv") = t);
     
      Rcpp::NumericVector coxtime = cox["time"];  
      Rcpp::NumericVector times(min((R_xlen_t)100, (R_xlen_t)coxtime.length()));
      for (int j = 0; j < times.length(); j++)
        times[j] = coxtime[round(((double)coxtime.length() - 1) * j / (times.length() - 1))];
               
      data = Rcpp::clone(data1);
      recurrence = Rcpp::wrap(data[colTime]);
      status = Rcpp::wrap(data[colStatus]);

      Rcpp::List fit = survest(cox, data, Rcpp::Named("times") = times, Rcpp::Named("se.fit") = f, Rcpp::Named("conf.int") = f);
      Rcpp::NumericMatrix est = fit["surv"];
    
      for (int i = 0; i < data.nrows(); i++)
      if (status[i] == 0)
      {
        Rcpp::NumericVector surv = est(i,_);
        while(1)
        {
          Rcpp::NumericVector badluck = runif(2);            
          double t = SurvTime(year[i], age[i], badluck[1], sex[i]);  
          
          double mdeath = t;
          double time = 0;
          double TT = maxTime[i];
          
          if (badluck[0] < surv[surv.length() - 1])
          {
            time = TT + 1 * 365.2425;
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
    return Rcpp::wrap(recurrence);
  }
  //catch(...){}
  throw std::range_error("Unknown SimCensorX() Error");
}

// [[Rcpp::export]]
SEXP SurvExpPrep(Rcpp::DataFrame D1, double time)
{
  if (SurvExpCache != NULL)
    throw std::range_error("SurvExpCache is NULL");  
  
  try
  {
    Rcpp::DataFrame D = Rcpp::clone(D1);
    Rcpp::NumericVector age = Rcpp::wrap(D["age"]);
    Rcpp::NumericVector sex = Rcpp::wrap(D["sex"]);
    Rcpp::NumericVector year = Rcpp::wrap(D["year"]);
    
    Rcpp::NumericVector recurrence = Rcpp::clone(age);
    for (int i = 0; i < D.nrows(); i++)
        recurrence[i] = SurvProbability(year[i]/365.25 + 1960, age[i], time, sex[i]);
    return Rcpp::wrap(recurrence);
  }
  catch(...){}
  throw std::range_error("Unknown SurvExpPrep() Error");
}
