#include <Rcpp.h>
#include <string>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
DataFrame Sort(Rcpp::DataFrame data, Rcpp::NumericVector column)
{
  Rcpp::Language call("order", column);
  Rcpp::IntegerVector ndx = call.eval();
  
  Rcpp::DataFrame newData = Rcpp::clone(data);
  for (int i = 0; i < data.length(); i++)
  {
    {
      Rcpp::NumericVector vector = data[i];
      Rcpp::NumericVector newVector = newData[i];
      for (int j = 0; j < data.nrows(); j++)
        newVector[j] = vector[ndx[j] - 1];
    }
    {
      Rcpp::IntegerVector vector = data[i];
      Rcpp::IntegerVector newVector = newData[i];
      for (int j = 0; j < data.nrows(); j++)
        newVector[j] = vector[ndx[j] - 1];
    }
  }
  return newData;
}

// [[Rcpp::export]]
SEXP SimKM(Rcpp::DataFrame data1, double interval)
{
  try
  {
    Rcpp::DataFrame R;
    {
      Rcpp::DataFrame data = Rcpp::clone(data1);
      Rcpp::NumericVector agex((SEXP)data["age"]);
      Rcpp::NumericVector birthyear((SEXP)data["birthyear"]);
      Rcpp::NumericVector firstint((SEXP)data["firstint"]);
      Rcpp::NumericVector ttx((SEXP)data["recurrence"]);
      Rcpp::NumericVector koncx(agex.length());
  
      for (int i = 0; i < agex.length(); i++)
      {
        agex[i] = agex[i] * 365.2425;
        ttx[i] = ttx[i] * 365.2425;    
        koncx[i] = (interval - firstint[i]) * 365.2425;
      }  
      R = Rcpp::DataFrame::create(Rcpp::Named("age") = agex, Rcpp::Named("year") = birthyear, Rcpp::Named("sex") = data["sex"], Rcpp::Named("status") = data["status"], Rcpp::Named("recurrence") = ttx, Rcpp::Named("konc") = koncx);
      R = Sort(R, ttx);
    }
      
    Rcpp::NumericVector tt((SEXP)R["recurrence"]);
    Rcpp::NumericVector ttun(1, 0.0);
    {
      Rcpp::NumericVector konc = Rcpp::clone((SEXP)R["konc"]);
      std::sort(konc.begin(), konc.end());
  
      for (int i = 0, j = 0; (i < tt.length()) || (j < konc.length());)
      {
        double last = ttun[ttun.length() - 1];
        double candidate = 0;
        if(i < tt.length())
          if (j < konc.length())
            if (tt[i] < konc[j])
              candidate = tt[i++];
            else
              candidate = konc[j++];
          else
            candidate = tt[i++];
        else
          candidate = konc[j++];
          
        if (candidate > last)
          ttun.push_back(candidate);
      }
    }
    
    Rcpp::NumericVector difft(ttun.length() - 1);
    for (int i = 0; i < difft.length(); i++)
      difft[i] = ttun[i + 1] - ttun[i];
   
    Rcpp::NumericVector le(ttun.length() - 1);
    Rcpp::NumericVector Yov(ttun.length() - 1, 0.0);
    Rcpp::NumericVector Yosebni(R.nrows(), (double)1);
    Rcpp::NumericVector surv(ttun.length(), 1.0);
  
    Rcpp::NumericVector status((SEXP)R["status"]);
    Rcpp::NumericVector konc((SEXP)R["konc"]);
  
    Rcpp::Function ExpPrep = Rcpp::Environment::global_env()["ExpPrep"];
    Rcpp::NumericMatrix Spi = ExpPrep(R, difft);
    double cumsumLe = 0;
    for (int it = 1; it < ttun.length(); it++)
    {    
      Yov[it - 1] = 0;
      surv[it - 1] = exp(-cumsumLe);
      for (int j = 0; j < Yosebni.length(); j++)
      {
        Yosebni[j] *= Spi(j, it - 1);
        Yov[it - 1] += Yosebni[j];
      }
        
      int Ne = 0;
      for (int j = 0; j < tt.length(); j++)
      {
        if (tt[j] == ttun[it])
        {
          if (status[j] == 1)
            Ne++;
          Yosebni[j] -= 1;
        }
        
        if (konc[j] == ttun[it])
          Yosebni[j] = 0;
      }
   
      le[it - 1] = Ne / Yov[it - 1];
      cumsumLe += le[it - 1];
      ttun[it] /= 365.2425;
    }
    surv[surv.length() - 1] = exp(-cumsumLe);
  
    return Rcpp::List::create(Rcpp::Named("time") = ttun, Rcpp::Named("surv") = surv, Rcpp::Named("lambda") = le, Rcpp::Named("yov") = Yov, Rcpp::Named("R") = R, Rcpp::Named("difft") = difft);
  }
  catch (...){}
  throw std::range_error("Unknown SimKM() Error");
}
