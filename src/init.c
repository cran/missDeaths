#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP missDeaths_rcpp_SimCensorX(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP missDeaths_rcpp_SurvExpInit(SEXP);
extern SEXP missDeaths_rcpp_SurvTime(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"missDeaths_rcpp_SimCensorX",  (DL_FUNC) &missDeaths_rcpp_SimCensorX,  5},
  {"missDeaths_rcpp_SurvExpInit", (DL_FUNC) &missDeaths_rcpp_SurvExpInit, 1},
  {"missDeaths_rcpp_SurvTime",    (DL_FUNC) &missDeaths_rcpp_SurvTime,    4},
  {NULL, NULL, 0}
};

void R_init_missDeaths(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
