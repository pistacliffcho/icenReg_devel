#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
*/
  
  /* .Call calls */
  extern SEXP _icenReg_computeConditional_p(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _icenReg_computeConditional_q(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _icenReg_ic_parList(SEXP);
extern SEXP _icenReg_R_ic_bayes(SEXP, SEXP, SEXP);
extern SEXP dGeneralGamma(SEXP, SEXP, SEXP, SEXP);
extern SEXP EMICM(SEXP, SEXP, SEXP, SEXP);
extern SEXP fastNumericInsert(SEXP, SEXP, SEXP);
extern SEXP findMI(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ic_sp_ch(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP optCliq(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP q_regTrans(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP qGeneralGamma(SEXP, SEXP, SEXP, SEXP);
extern SEXP s_regTrans(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_icenReg_computeConditional_p", (DL_FUNC) &_icenReg_computeConditional_p,  5},
  {"_icenReg_computeConditional_q", (DL_FUNC) &_icenReg_computeConditional_q,  5},
  {"_icenReg_ic_parList",           (DL_FUNC) &_icenReg_ic_parList,            1},
  {"_icenReg_R_ic_bayes",           (DL_FUNC) &_icenReg_R_ic_bayes,            3},
  {"dGeneralGamma",                 (DL_FUNC) &dGeneralGamma,                  4},
  {"EMICM",                         (DL_FUNC) &EMICM,                          4},
  {"fastNumericInsert",             (DL_FUNC) &fastNumericInsert,              3},
  {"findMI",                        (DL_FUNC) &findMI,                         5},
  {"ic_sp_ch",                      (DL_FUNC) &ic_sp_ch,                      11},
  {"optCliq",                       (DL_FUNC) &optCliq,                        5},
  {"q_regTrans",                    (DL_FUNC) &q_regTrans,                     5},
  {"qGeneralGamma",                 (DL_FUNC) &qGeneralGamma,                  4},
  {"s_regTrans",                    (DL_FUNC) &s_regTrans,                     5},
  {NULL, NULL, 0}
};

void R_init_icenReg(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
