#include <Rcpp.h>
#include <string>
#include "survexpcache.h"

using namespace std;

class SurvCurve
{
//private:
public:
  Rcpp::NumericVector Times;
  Rcpp::NumericVector Curve;
  int Year;
  
  int Find(double time, int low, int high)
  {
    if (low >= high)
      return max(0, high);
    
		int mid = (low + high) / 2;
		if (time <= Times[mid])
			  return Find(time, low, mid);
			else 
        return (mid == low) ? high : Find(time, mid, high);
  }
  
public:
  SurvCurve(Rcpp::NumericVector curve, Rcpp::NumericVector times, int year)
  {
    Curve = curve;
    Times = times;
    Year = year;
  }
  double Probability(double time)
  {
    if (time < 0) time= 0;
    int i = Find(time, 0, Times.length()-1);
    
    if (time == Times[i])
      return Curve[i];
    
    double t1 = (i == 0) ? 0 : Times[i - 1];
    double t2 = Times[i];
    double s1 = (i == 0) ? 1 : Curve[i - 1];
    double s2 = Curve[i];
          
    return s1 - (time - t1) / (t2 - t1) * (s1 - s2);
  }
  double Probability2(double time)
  {
    if (time < 0) time= 0;
    for (int i = 0; i < Times.size(); i++)
    {
      if (Times[i] > time)
      {
          double t1 = (i == 0) ? 0 : Times[i - 1];
          double t2 = Times[i];
          double s1 = (i == 0) ? 1 : Curve[i - 1];
          double s2 = Curve[i];
          
          return s1 - (time - t1) / (t2 - t1) * (s1 - s2);
      }
    }
    return -1;  
  }
  double Age(double prob)
  {
    if (prob > 1) prob = 1;
    if (prob < 0) prob = 0;
    
    for (int i = 0; i < Curve.size(); i++)
      if (Curve[i] < prob)
      {
          double t1 = (i == 0) ? 0 : Times[i - 1];
          double t2 = Times[i];
          double p1 = (i == 0) ? 1 : Curve[i - 1];
          double p2 = Curve[i];
          
          return (t1 + (p1 - prob) * (t2 - t1) / (p1 - p2));
      }
    return -1;
  }
  double Time(double age, double prob)
  {
    double initprob = Probability(age);
    if (initprob >= 0)
    {
      double time = Age(prob * initprob);
      if (time >= 0)
      {
          return (time - age);
      }
    }
    return -1;
  }  
  int BirthYear() { return Year; }
};

class SurvExp
{
private:
  SurvCurve** FemaleCache;
  SurvCurve** MaleCache;
  int Length;
  
  void InitCache(int start, int end, Rcpp::NumericVector times, int sex, SurvCurve** cache, SEXP poptable)
  {
    Rcpp::Formula formula("~1");
    bool conditional = TRUE;
    for (int survyear = start; survyear <= end; survyear++)
    {
      Rcpp::NumericVector ptage = Rcpp::NumericVector(1, 0.0);
      Rcpp::IntegerVector ptsex = Rcpp::IntegerVector(1, 1);
      ptsex[0] = sex;
      Rcpp::NumericVector ptyear = Rcpp::NumericVector(1, (survyear - 1960) * 365.2425);
      Rcpp::DataFrame data = Rcpp::DataFrame::create(Rcpp::Named("age") = ptage, Rcpp::Named("sex") = ptsex, Rcpp::Named("year") = ptyear);
  
      Rcpp::Environment stats("package:survival");
      Rcpp::Function survexp = stats["survexp"]; 
      Rcpp::List list(survexp(formula, data, Rcpp::Named("times") = times, Rcpp::Named("conditional") = conditional, Rcpp::Named("ratetable") = poptable, Rcpp::Named("method") = "ederer"));
      Rcpp::NumericVector surv = list["surv"];
  
      cache[survyear - start] = new SurvCurve(surv, times, survyear);
    }
  }
  
public:
  SurvExp(int precision, int start, int end, SEXP poptable)
  { 
    if (end < start) end = start;
    Length = end - start + 1;
    FemaleCache = new SurvCurve*[Length];
    MaleCache = new SurvCurve*[Length];
  
    Rcpp::NumericVector times = Rcpp::NumericVector(150 * precision);
    for (int i = 0; i < times.size(); i++)
      times[i] = ((double)i / (double)precision) * 365.2425;
    
    InitCache(start, end, times, 1, MaleCache, poptable);
    InitCache(start, end, times, 2, FemaleCache, poptable);
  }
  SurvCurve* Get(int year, int sex)
  {
    SurvCurve** cache = (sex == 2) ? FemaleCache : MaleCache;
    
    for (int i = 0; i < Length; i++)
      if(cache[i]->BirthYear() == year)
        return cache[i];
    return NULL;
  }
};

SurvExp* SurvExpCache;

// [[Rcpp::export]]
void SurvExpInit(int precision, int start, int end, SEXP poptable)
{
  SurvExpCache = new SurvExp(precision, start, end, poptable);
}

// [[Rcpp::export]]
SEXP SurvProbability(int year, double age, int sex)
{
  if (SurvExpCache != NULL)
  {
    SurvCurve* curve = SurvExpCache->Get(year, sex);
    if (curve != NULL)
      return Rcpp::wrap(curve->Probability(age));//-curve->Probability2(age));
  }
  return Rcpp::wrap(-1);
}


// [[Rcpp::export]]
SEXP SurvTimeX(int year, double probability, int sex)
{
  if (SurvExpCache != NULL)
  {
    SurvCurve* curve = SurvExpCache->Get(year, sex);
    if (curve != NULL)
      return Rcpp::wrap(curve->Age(probability));
  }
  return Rcpp::wrap(-1);
}

double SurvTime(int year, double age, double probability, int sex)
{
  if (SurvExpCache != NULL)
  {
    SurvCurve* curve = SurvExpCache->Get(year, sex);
    if (curve != NULL)
      return curve->Time(age, probability);
  }
  return -1;
}

// [[Rcpp::export]]
SEXP ExpPrep(Rcpp::DataFrame data, Rcpp::NumericVector difft)
{
  Rcpp::NumericVector age((SEXP)data["age"]);
  Rcpp::IntegerVector sex((SEXP)data["sex"]);
  Rcpp::NumericVector year((SEXP)data["year"]);

  Rcpp::NumericMatrix expprep(age.length(), difft.length());
  if (SurvExpCache != NULL)
    for (int i = 0; i < age.length(); i++)
    {
      SurvCurve* curve = SurvExpCache->Get(year[i], sex[i]);
      if (curve != NULL)
      {
        double currentAge = age[i];
        double p = curve->Probability(currentAge);
        for (int j = 0; j < difft.length(); j++)
        {
          currentAge += difft[j];
          double p1 = curve->Probability(currentAge);
          expprep(i, j) = p1 / p;
          p = p1;
        }
      }
    }
  return Rcpp::wrap(expprep);
}
