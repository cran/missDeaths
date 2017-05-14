#include <Rcpp.h>
#include <string>
#include <time.h>
#include "survexpcache.h"

using namespace Rcpp;
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
  SEXP RateTable;
  
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
  SurvExp(SEXP poptable)
  { 
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    
    int precision = 12;
    int start = 1850;
    int end = tm.tm_year + 1900;
    
    Length = end - start + 1;
    FemaleCache = new SurvCurve*[Length];
    MaleCache = new SurvCurve*[Length];
  
    Rcpp::NumericVector times = Rcpp::NumericVector(150 * precision);
    for (int i = 0; i < times.size(); i++)
      times[i] = ((double)i / (double)precision) * 365.2425;
    
    InitCache(start, end, times, 1, MaleCache, poptable);
    InitCache(start, end, times, 2, FemaleCache, poptable);
    RateTable = poptable;
  }
  bool Check(SEXP poptable)
  {
    return (RateTable == poptable);
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

// [[Rcpp::export]]
void SurvExpInit(SEXP poptable)
{
  if (SurvExpCache != NULL)
    if (SurvExpCache->Check(poptable))
      return;
    
  SurvExpCache = new SurvExp(poptable);
}

// [[Rcpp::export]]
double SurvTime(double birthyear, double age, double probability, int sex)
{
  if (SurvExpCache != NULL)
  {
    int year = floor(birthyear);
    SurvCurve* curve1 = SurvExpCache->Get(year, sex);
    SurvCurve* curve2 = SurvExpCache->Get(year + 1, sex);
    
    if ((curve1 != NULL) && (curve2 != NULL))
    {
      double time1 = curve1->Time(age, probability);
      double time2 = curve2->Time(age, probability);
      
      return time1 + (time2 - time1) * (birthyear - (double)year);
    }
  }
  return -1;
}

// [[Rcpp::export]]
double SurvProbability(double birthyear, double age, double time, int sex)
{
  if (SurvExpCache != NULL)
  {
    int year = floor(birthyear);
    SurvCurve* curve1 = SurvExpCache->Get(year, sex);
    //SurvCurve* curve2 = SurvExpCache->Get(year + 1, sex);
    
    if ((curve1 != NULL)/* && (curve2 != NULL)*/)
    {
      double prob1 = curve1->Probability(age + time) / curve1->Probability(age);
      //double prob2 = curve2->Probability(age + time) / curve2->Probability(age);
      
      return prob1;// + (prob2 - prob1) * (birthyear - (double)year);
    }
  }
  return 1; 
}
