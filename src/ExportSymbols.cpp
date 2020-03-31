#include "collapse.h"
#include <Rcpp.h>
using namespace Rcpp;

RcppExport void multi_yw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
RcppExport SEXP collapse_init(SEXP);
RcppExport SEXP dt_na(SEXP, SEXP);
RcppExport SEXP forder(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP frank(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP pacf1(SEXP, SEXP);
RcppExport SEXP rbindlist(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP setcolorder(SEXP, SEXP);
RcppExport SEXP subsetDT(SEXP, SEXP, SEXP);
RcppExport SEXP subsetVector(SEXP, SEXP);
RcppExport SEXP uniqlengths(SEXP, SEXP);

static const R_CMethodDef CEntries[]  = {
  {"C_multi_yw", (DL_FUNC) &multi_yw, 10},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"Cpp_BW", (DL_FUNC) &_collapse_BWCpp, 9},
  {"Cpp_BWm", (DL_FUNC) &_collapse_BWmCpp, 9},
  {"Cpp_BWl", (DL_FUNC) &_collapse_BWlCpp, 9},
  {"Cpp_TRA", (DL_FUNC) &_collapse_TRACpp, 4},
  {"Cpp_TRAm", (DL_FUNC) &_collapse_TRAmCpp, 4},
  {"Cpp_TRAl", (DL_FUNC) &_collapse_TRAlCpp, 4},
  {"Cpp_fNdistinct", (DL_FUNC) &_collapse_fNdistinctCpp, 5},
  {"Cpp_fNdistinctl", (DL_FUNC) &_collapse_fNdistinctlCpp, 6},
  {"Cpp_fNdistinctm", (DL_FUNC) &_collapse_fNdistinctmCpp, 6},
  {"Cpp_fNobs", (DL_FUNC) &_collapse_fNobsCpp, 3},
  {"Cpp_fNobsm", (DL_FUNC) &_collapse_fNobsmCpp, 4},
  {"Cpp_fNobsl", (DL_FUNC) &_collapse_fNobslCpp, 4},
  {"Cpp_fbstats", (DL_FUNC) &_collapse_fbstatsCpp, 10},
  {"Cpp_fbstatsm", (DL_FUNC) &_collapse_fbstatsmCpp, 9},
  {"Cpp_fbstatsl", (DL_FUNC) &_collapse_fbstatslCpp, 9},
  {"Cpp_fdiff", (DL_FUNC) &_collapse_fdiffCpp, 9},
  {"Cpp_fdiffm", (DL_FUNC) &_collapse_fdiffmCpp, 9},
  {"Cpp_fdiffl", (DL_FUNC) &_collapse_fdifflCpp, 9},
  {"Cpp_ffirst", (DL_FUNC) &_collapse_ffirstCpp, 4},
  {"Cpp_ffirstm", (DL_FUNC) &_collapse_ffirstmCpp, 5},
  {"Cpp_ffirstl", (DL_FUNC) &_collapse_ffirstlCpp, 4},
  {"Cpp_fgrowth", (DL_FUNC) &_collapse_fgrowthCpp, 10},
  {"Cpp_fgrowthm", (DL_FUNC) &_collapse_fgrowthmCpp, 10},
  {"Cpp_fgrowthl", (DL_FUNC) &_collapse_fgrowthlCpp, 10},
  {"Cpp_flaglead", (DL_FUNC) &_collapse_flagleadCpp, 8},
  {"Cpp_flagleadm", (DL_FUNC) &_collapse_flagleadmCpp, 8},
  {"Cpp_flagleadl", (DL_FUNC) &_collapse_flagleadlCpp, 8},
  {"Cpp_flast", (DL_FUNC) &_collapse_flastCpp, 4},
  {"Cpp_flastm", (DL_FUNC) &_collapse_flastmCpp, 5},
  {"Cpp_flastl", (DL_FUNC) &_collapse_flastlCpp, 4},
  {"Cpp_fmax", (DL_FUNC) &_collapse_fmaxCpp, 4},
  {"Cpp_fmaxm", (DL_FUNC) &_collapse_fmaxmCpp, 5},
  {"Cpp_fmaxl", (DL_FUNC) &_collapse_fmaxlCpp, 5},
  {"Cpp_fmean", (DL_FUNC) &_collapse_fmeanCpp, 6},
  {"Cpp_fmeanm", (DL_FUNC) &_collapse_fmeanmCpp, 7},
  {"Cpp_fmeanl", (DL_FUNC) &_collapse_fmeanlCpp, 7},
  {"Cpp_fmedian", (DL_FUNC) &_collapse_fmedianCpp, 5},
  {"Cpp_fmedianm", (DL_FUNC) &_collapse_fmedianmCpp, 6},
  {"Cpp_fmedianl", (DL_FUNC) &_collapse_fmedianlCpp, 6},
  {"Cpp_fmin", (DL_FUNC) &_collapse_fminCpp, 4},
  {"Cpp_fminm", (DL_FUNC) &_collapse_fminmCpp, 5},
  {"Cpp_fminl", (DL_FUNC) &_collapse_fminlCpp, 5},
  {"Cpp_fmode", (DL_FUNC) &_collapse_fmodeCpp, 6},
  {"Cpp_fmodel", (DL_FUNC) &_collapse_fmodelCpp, 6},
  {"Cpp_fmodem", (DL_FUNC) &_collapse_fmodemCpp, 7},
  {"Cpp_fprod", (DL_FUNC) &_collapse_fprodCpp, 5},
  {"Cpp_fprodm", (DL_FUNC) &_collapse_fprodmCpp, 6},
  {"Cpp_fprodl", (DL_FUNC) &_collapse_fprodlCpp, 6},
  {"Cpp_fscale", (DL_FUNC) &_collapse_fscaleCpp, 7},
  {"Cpp_fscalem", (DL_FUNC) &_collapse_fscalemCpp, 7},
  {"Cpp_fscalel", (DL_FUNC) &_collapse_fscalelCpp, 7},
  {"Cpp_fsum", (DL_FUNC) &_collapse_fsumCpp, 5},
  {"Cpp_fsumm", (DL_FUNC) &_collapse_fsummCpp, 6},
  {"Cpp_fsuml", (DL_FUNC) &_collapse_fsumlCpp, 6},
  {"Cpp_fvarsd", (DL_FUNC) &_collapse_fvarsdCpp, 8},
  {"Cpp_fvarsdm", (DL_FUNC) &_collapse_fvarsdmCpp, 9},
  {"Cpp_fvarsdl", (DL_FUNC) &_collapse_fvarsdlCpp, 9},
  {"Cpp_mrtl", (DL_FUNC) &_collapse_mrtl, 3},
  {"Cpp_mctl", (DL_FUNC) &_collapse_mctl, 3},
  // {"Cpp_na_rm", (DL_FUNC) &_collapse_na_rm, 1},
  // {"Cpp_fanyNAint", (DL_FUNC) &_collapse_fanyNAint, 1},
  {"Cpp_psmat", (DL_FUNC) &_collapse_psmatCpp, 4},
  {"Cpp_qF", (DL_FUNC) &_collapse_qFCpp, 3},
  {"Cpp_qG", (DL_FUNC) &_collapse_qGCpp, 3},
  {"Cpp_funique", (DL_FUNC) &_collapse_funique, 2},
  {"Cpp_setAttributes", (DL_FUNC) &_collapse_setAttributes, 2},
  {"Cpp_setattributes", (DL_FUNC) &_collapse_setattributes, 2},
  {"Cpp_setAttr", (DL_FUNC) &_collapse_setAttr, 3},
  {"Cpp_setattr_clp", (DL_FUNC) &_collapse_setattr_clp, 3},
  {"Cpp_duplAttributes", (DL_FUNC) &_collapse_duplAttributes, 2},
  {"Cpp_duplattributes", (DL_FUNC) &_collapse_duplattributes, 2},
  {"Cpp_cond_duplAttributes", (DL_FUNC) &_collapse_cond_duplAttributes, 2},
  {"Cpp_cond_duplattributes", (DL_FUNC) &_collapse_cond_duplattributes, 2},
  {"Cpp_groups2GRP", (DL_FUNC) &_collapse_groups2GRPCpp, 3},
  {"Cpp_lassign", (DL_FUNC) &_collapse_lassignCpp, 4},
  {"C_collapse_init", (DL_FUNC) &collapse_init, 1},
  {"C_dt_na",         (DL_FUNC) &dt_na,         2},
  {"C_forder",        (DL_FUNC) &forder,        6},
  {"C_frank",         (DL_FUNC) &frank,         4},
  {"C_pacf1",         (DL_FUNC) &pacf1,         2},
  {"C_rbindlist",     (DL_FUNC) &rbindlist,     4},
  {"C_setcolorder",   (DL_FUNC) &setcolorder,   2},
  {"C_subsetDT",      (DL_FUNC) &subsetDT,      3},
  {"C_subsetVector",  (DL_FUNC) &subsetVector,  2},
  {"C_uniqlengths",   (DL_FUNC) &uniqlengths,   2},
  {NULL, NULL, 0}
};

RcppExport void R_init_collapse(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
