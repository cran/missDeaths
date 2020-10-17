double SurvTime(double year, double age, double probability, int sex);
double SurvProbability(double birthyear, double age, double time, int sex);
SEXP SurvDump(int year, int sex);

class SurvExp;
static SurvExp* SurvExpCache;
