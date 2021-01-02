#include "collapse.h"
#include <Rcpp.h>
using namespace Rcpp;

// prefix with RcppExport ? -> Yes, necessary !
RcppExport void multi_yw(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
RcppExport SEXP collapse_init(SEXP);
RcppExport SEXP dt_na(SEXP, SEXP);
RcppExport SEXP allNAv(SEXP, SEXP);
RcppExport SEXP Cradixsort(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP frankds(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP pacf1(SEXP, SEXP);
RcppExport SEXP rbindlist(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP setcolorder(SEXP, SEXP);
RcppExport SEXP subsetDT(SEXP, SEXP, SEXP);
RcppExport SEXP subsetVector(SEXP, SEXP);
RcppExport SEXP setAttributes(SEXP x, SEXP a);
RcppExport void setattributes(SEXP x, SEXP a);
// RcppExport SEXP CsetAttr(SEXP object, SEXP a, SEXP v); -> mot more efficeint than attr i.e. for row.names...
RcppExport void setattr(SEXP x, SEXP a, SEXP v);
RcppExport SEXP duplAttributes(SEXP x, SEXP y);
RcppExport void duplattributes(SEXP x, SEXP y);
RcppExport SEXP cond_duplAttributes(SEXP x, SEXP y);
RcppExport SEXP CsetAttrib(SEXP object, SEXP a);
RcppExport SEXP CcopyAttrib(SEXP to, SEXP from);
RcppExport SEXP CcopyMostAttrib(SEXP to, SEXP from);
RcppExport SEXP lassign(SEXP x, SEXP s, SEXP rows, SEXP fill);
RcppExport SEXP groups2GRP(SEXP x, SEXP lx, SEXP gs);
RcppExport SEXP Cna_rm(SEXP x);

static const R_CMethodDef CEntries[]  = {
  {"C_multi_yw", (DL_FUNC) &multi_yw, 10},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"Cpp_BW", (DL_FUNC) &_collapse_BWCpp, 10},
  {"Cpp_BWm", (DL_FUNC) &_collapse_BWmCpp, 10},
  {"Cpp_BWl", (DL_FUNC) &_collapse_BWlCpp, 10},
  {"Cpp_TRA", (DL_FUNC) &_collapse_TRACpp, 4},
  {"Cpp_TRAm", (DL_FUNC) &_collapse_TRAmCpp, 4},
  {"Cpp_TRAl", (DL_FUNC) &_collapse_TRAlCpp, 4},
  {"Cpp_fNdistinct", (DL_FUNC) &_collapse_fNdistinctCpp, 5},
  {"Cpp_fNdistinctl", (DL_FUNC) &_collapse_fNdistinctlCpp, 6},
  {"Cpp_fNdistinctm", (DL_FUNC) &_collapse_fNdistinctmCpp, 6},
  {"Cpp_pwNobsm", (DL_FUNC) &_collapse_pwNobsmCpp, 1},
  {"Cpp_fNobs", (DL_FUNC) &_collapse_fNobsCpp, 3},
  {"Cpp_fNobsm", (DL_FUNC) &_collapse_fNobsmCpp, 4},
  {"Cpp_fNobsl", (DL_FUNC) &_collapse_fNobslCpp, 4},
  {"Cpp_varying", (DL_FUNC) &_collapse_varyingCpp, 4},
  {"Cpp_varyingm", (DL_FUNC) &_collapse_varyingmCpp, 5},
  {"Cpp_varyingl", (DL_FUNC) &_collapse_varyinglCpp, 5},
  {"Cpp_fbstats", (DL_FUNC) &_collapse_fbstatsCpp, 10},
  {"Cpp_fbstatsm", (DL_FUNC) &_collapse_fbstatsmCpp, 9},
  {"Cpp_fbstatsl", (DL_FUNC) &_collapse_fbstatslCpp, 9},
  {"Cpp_ffirst", (DL_FUNC) &_collapse_ffirstCpp, 4},
  {"Cpp_ffirstm", (DL_FUNC) &_collapse_ffirstmCpp, 5},
  {"Cpp_ffirstl", (DL_FUNC) &_collapse_ffirstlCpp, 4},
  {"Cpp_fdiffgrowth", (DL_FUNC) &_collapse_fdiffgrowthCpp, 12},
  {"Cpp_fdiffgrowthm", (DL_FUNC) &_collapse_fdiffgrowthmCpp, 12},
  {"Cpp_fdiffgrowthl", (DL_FUNC) &_collapse_fdiffgrowthlCpp, 12},
  {"Cpp_flaglead", (DL_FUNC) &_collapse_flagleadCpp, 7},
  {"Cpp_flagleadm", (DL_FUNC) &_collapse_flagleadmCpp, 7},
  {"Cpp_flagleadl", (DL_FUNC) &_collapse_flagleadlCpp, 7},
  {"Cpp_flast", (DL_FUNC) &_collapse_flastCpp, 4},
  {"Cpp_flastm", (DL_FUNC) &_collapse_flastmCpp, 5},
  {"Cpp_flastl", (DL_FUNC) &_collapse_flastlCpp, 4},
  {"Cpp_fminmax", (DL_FUNC) &_collapse_fminmaxCpp, 5},
  {"Cpp_fminmaxm", (DL_FUNC) &_collapse_fminmaxmCpp, 6},
  {"Cpp_fminmaxl", (DL_FUNC) &_collapse_fminmaxlCpp, 6},
  {"Cpp_fmean", (DL_FUNC) &_collapse_fmeanCpp, 6},
  {"Cpp_fmeanm", (DL_FUNC) &_collapse_fmeanmCpp, 7},
  {"Cpp_fmeanl", (DL_FUNC) &_collapse_fmeanlCpp, 7},
  // {"Cpp_fmedian", (DL_FUNC) &_collapse_fmedianCpp, 6},
  // {"Cpp_fmedianm", (DL_FUNC) &_collapse_fmedianmCpp, 7},
  // {"Cpp_fmedianl", (DL_FUNC) &_collapse_fmedianlCpp, 7},
  {"Cpp_fnth", (DL_FUNC) &_collapse_fnthCpp, 8},
  {"Cpp_fnthm", (DL_FUNC) &_collapse_fnthmCpp, 9},
  {"Cpp_fnthl", (DL_FUNC) &_collapse_fnthlCpp, 9},
  {"Cpp_fmode", (DL_FUNC) &_collapse_fmodeCpp, 7},
  {"Cpp_fmodel", (DL_FUNC) &_collapse_fmodelCpp, 7},
  {"Cpp_fmodem", (DL_FUNC) &_collapse_fmodemCpp, 8},
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
  {"Cpp_psmat", (DL_FUNC) &_collapse_psmatCpp, 4},
  {"Cpp_qF", (DL_FUNC) &_collapse_qFCpp, 6},
  // {"Cpp_qG", (DL_FUNC) &_collapse_qGCpp, 5},
  {"Cpp_funique", (DL_FUNC) &_collapse_funiqueCpp, 2},
  {"Cpp_fdroplevels", (DL_FUNC) &_collapse_fdroplevelsCpp, 2},
  {"C_setAttributes", (DL_FUNC) &setAttributes, 2},
  {"C_setattributes", (DL_FUNC) &setattributes, 2},
  // {"C_setAttr", (DL_FUNC) &CsetAttr, 3},
  {"C_setattr", (DL_FUNC) &setattr, 3},
  {"C_duplAttributes", (DL_FUNC) &duplAttributes, 2},
  {"C_duplattributes", (DL_FUNC) &duplattributes, 2},
  {"C_cond_duplAttributes", (DL_FUNC) &cond_duplAttributes, 2},
  // {"C_cond_duplattributes", (DL_FUNC) &cond_duplattributes, 2},
  {"C_setAttrib", (DL_FUNC) &CsetAttrib, 2},
  {"C_copyAttrib", (DL_FUNC) &CcopyAttrib, 2},
  {"C_copyMostAttrib", (DL_FUNC) &CcopyMostAttrib, 2},
  {"C_groups2GRP", (DL_FUNC) &groups2GRP, 3},
  {"C_lassign", (DL_FUNC) &lassign, 4},
  {"Cpp_seqid", (DL_FUNC) &_collapse_seqid, 7},
  {"Cpp_groupid", (DL_FUNC) &_collapse_groupid, 5},
  {"C_collapse_init", (DL_FUNC) &collapse_init, 1},
  {"C_dt_na",         (DL_FUNC) &dt_na,         2},
  {"C_allNA",         (DL_FUNC) &allNAv,        2},
  {"C_na_rm",         (DL_FUNC) &Cna_rm,        1},
  {"C_radixsort",     (DL_FUNC) &Cradixsort,    6},
  {"C_frankds",       (DL_FUNC) &frankds,       4},
  {"C_pacf1",         (DL_FUNC) &pacf1,         2},
  {"C_rbindlist",     (DL_FUNC) &rbindlist,     4},
  {"C_setcolorder",   (DL_FUNC) &setcolorder,   2},
  {"C_subsetDT",      (DL_FUNC) &subsetDT,      3},
  {"C_subsetVector",  (DL_FUNC) &subsetVector,  2},
  {NULL, NULL, 0}
};

RcppExport void R_init_collapse(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
