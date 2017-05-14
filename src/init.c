#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern SEXP missDeaths_SimCensorX(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP missDeaths_SurvExpPrep(SEXP, SEXP);
extern SEXP missDeaths_SurvExpInit(SEXP);
extern SEXP missDeaths_SurvTime(SEXP, SEXP, SEXP, SEXP);
extern SEXP missDeaths_SurvProbability(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"missDeaths_SimCensorX",  (DL_FUNC) &missDeaths_SimCensorX,  5},
  {"missDeaths_SurvExpPrep",  (DL_FUNC) &missDeaths_SurvExpPrep,  2},
  {"missDeaths_SurvExpInit", (DL_FUNC) &missDeaths_SurvExpInit, 1},
  {"missDeaths_SurvTime",    (DL_FUNC) &missDeaths_SurvTime,    4},
  {"missDeaths_SurvProbability",    (DL_FUNC) &missDeaths_SurvProbability,    4},
  {NULL, NULL, 0}
};

void R_init_missDeaths(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
