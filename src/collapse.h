#include <Rcpp.h>
using namespace Rcpp;

// BWCpp
RcppExport SEXP _collapse_BWCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP set_meanSEXP, SEXP BSEXP, SEXP fillSEXP);
// BWmCpp
RcppExport SEXP _collapse_BWmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP set_meanSEXP, SEXP BSEXP, SEXP fillSEXP);
// BWlCpp
RcppExport SEXP _collapse_BWlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP set_meanSEXP, SEXP BSEXP, SEXP fillSEXP);
// TRACpp
RcppExport SEXP _collapse_TRACpp(SEXP xSEXP, SEXP xAGSEXP, SEXP gSEXP, SEXP retSEXP);
// TRAmCpp
RcppExport SEXP _collapse_TRAmCpp(SEXP xSEXP, SEXP xAGSEXP, SEXP gSEXP, SEXP retSEXP);
// TRAlCpp
RcppExport SEXP _collapse_TRAlCpp(SEXP xSEXP, SEXP xAGSEXP, SEXP gSEXP, SEXP retSEXP);
// fNdistinctCpp
RcppExport SEXP _collapse_fNdistinctCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP narmSEXP);
// fNdistinctlCpp
RcppExport SEXP _collapse_fNdistinctlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fNdistinctmCpp
RcppExport SEXP _collapse_fNdistinctmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fNobsCpp
RcppExport SEXP _collapse_fNobsCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP);
// fNobsmCpp
RcppExport SEXP _collapse_fNobsmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP dropSEXP);
// fNobslCpp
RcppExport SEXP _collapse_fNobslCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP dropSEXP);
// fbstatsCpp
RcppExport SEXP _collapse_fbstatsCpp(SEXP xSEXP, SEXP extSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP npgSEXP, SEXP pgSEXP, SEXP wSEXP, SEXP arraySEXP, SEXP setnSEXP, SEXP gnSEXP);
// fbstatsmCpp
RcppExport SEXP _collapse_fbstatsmCpp(SEXP xSEXP, SEXP extSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP npgSEXP, SEXP pgSEXP, SEXP wSEXP, SEXP arraySEXP, SEXP gnSEXP);
// fbstatslCpp
RcppExport SEXP _collapse_fbstatslCpp(SEXP xSEXP, SEXP extSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP npgSEXP, SEXP pgSEXP, SEXP wSEXP, SEXP arraySEXP, SEXP gnSEXP);
// fdiffCpp
RcppExport SEXP _collapse_fdiffCpp(SEXP xSEXP, SEXP nSEXP, SEXP diffSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP namesSEXP);
// fdiffmCpp
RcppExport SEXP _collapse_fdiffmCpp(SEXP xSEXP, SEXP nSEXP, SEXP diffSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP namesSEXP);
// fdifflCpp
RcppExport SEXP _collapse_fdifflCpp(SEXP xSEXP, SEXP nSEXP, SEXP diffSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP namesSEXP);
// ffirstCpp
RcppExport SEXP _collapse_ffirstCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP);
// ffirstmCpp
RcppExport SEXP _collapse_ffirstmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP dropSEXP);
// ffirstlCpp
RcppExport SEXP _collapse_ffirstlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP);
// fgrowthCpp
RcppExport SEXP _collapse_fgrowthCpp(SEXP xSEXP, SEXP nSEXP, SEXP diffSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP logdiffSEXP, SEXP namesSEXP);
// fgrowthmCpp
RcppExport SEXP _collapse_fgrowthmCpp(SEXP xSEXP, SEXP nSEXP, SEXP diffSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP logdiffSEXP, SEXP namesSEXP);
// fgrowthlCpp
RcppExport SEXP _collapse_fgrowthlCpp(SEXP xSEXP, SEXP nSEXP, SEXP diffSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP logdiffSEXP, SEXP namesSEXP);
// flagleadCpp
RcppExport SEXP _collapse_flagleadCpp(SEXP xSEXP, SEXP nSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP namesSEXP);
// flagleadmCpp
RcppExport SEXP _collapse_flagleadmCpp(SEXP xSEXP, SEXP nSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP namesSEXP);
// flagleadlCpp
RcppExport SEXP _collapse_flagleadlCpp(SEXP xSEXP, SEXP nSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP namesSEXP);
// flastCpp
RcppExport SEXP _collapse_flastCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP);
// flastmCpp
RcppExport SEXP _collapse_flastmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP dropSEXP);
// flastlCpp
RcppExport SEXP _collapse_flastlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP);
// fmaxCpp
RcppExport SEXP _collapse_fmaxCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP);
// fmaxmCpp
RcppExport SEXP _collapse_fmaxmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fmaxlCpp
RcppExport SEXP _collapse_fmaxlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fmeanCpp
RcppExport SEXP _collapse_fmeanCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP);
// fmeanmCpp
RcppExport SEXP _collapse_fmeanmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fmeanlCpp
RcppExport SEXP _collapse_fmeanlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fmedianCpp
RcppExport SEXP _collapse_fmedianCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP narmSEXP);
// fmedianmCpp
RcppExport SEXP _collapse_fmedianmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fmedianlCpp
RcppExport SEXP _collapse_fmedianlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fminCpp
RcppExport SEXP _collapse_fminCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP);
// fminmCpp
RcppExport SEXP _collapse_fminmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fminlCpp
RcppExport SEXP _collapse_fminlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fmodeCpp
RcppExport SEXP _collapse_fmodeCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP);
// fmodelCpp
RcppExport SEXP _collapse_fmodelCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP);
// fmodemCpp
RcppExport SEXP _collapse_fmodemCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fprodCpp
RcppExport SEXP _collapse_fprodCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP);
// fprodmCpp
RcppExport SEXP _collapse_fprodmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fprodlCpp
RcppExport SEXP _collapse_fprodlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fscaleCpp
RcppExport SEXP _collapse_fscaleCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP set_meanSEXP, SEXP set_sdSEXP);
// fscalemCpp
RcppExport SEXP _collapse_fscalemCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP set_meanSEXP, SEXP set_sdSEXP);
// fscalelCpp
RcppExport SEXP _collapse_fscalelCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP set_meanSEXP, SEXP set_sdSEXP);
// fsumCpp
RcppExport SEXP _collapse_fsumCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP);
// fsummCpp
RcppExport SEXP _collapse_fsummCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fsumlCpp
RcppExport SEXP _collapse_fsumlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fvarsdCpp
RcppExport SEXP _collapse_fvarsdCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP stable_algoSEXP, SEXP sdSEXP);
// fvarsdmCpp
RcppExport SEXP _collapse_fvarsdmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP stable_algoSEXP, SEXP sdSEXP, SEXP dropSEXP);
// fvarsdlCpp
RcppExport SEXP _collapse_fvarsdlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP stable_algoSEXP, SEXP sdSEXP, SEXP dropSEXP);
// mrtl
RcppExport SEXP _collapse_mrtl(SEXP XSEXP, SEXP namesSEXP, SEXP retSEXP);
// mctl
RcppExport SEXP _collapse_mctl(SEXP XSEXP, SEXP namesSEXP, SEXP retSEXP);
// na_rm
// RcppExport SEXP _collapse_na_rm(SEXP xSEXP);
// psmatCpp
RcppExport SEXP _collapse_psmatCpp(SEXP xSEXP, SEXP gSEXP, SEXP tSEXP, SEXP transposeSEXP);
// qFCpp
RcppExport SEXP _collapse_qFCpp(SEXP xSEXP, SEXP orderedSEXP, SEXP na_excludeSEXP);
// qGCpp
RcppExport SEXP _collapse_qGCpp(SEXP xSEXP, SEXP orderedSEXP, SEXP na_excludeSEXP);
// funique
RcppExport SEXP _collapse_funique(SEXP xSEXP, SEXP orderedSEXP);
// setAttributes
RcppExport SEXP _collapse_setAttributes(SEXP xSEXP, SEXP aSEXP);
// setattributes
RcppExport SEXP _collapse_setattributes(SEXP xSEXP, SEXP aSEXP);
// setAttr
RcppExport SEXP _collapse_setAttr(SEXP xSEXP, SEXP aSEXP, SEXP vSEXP);
// setattr_clp
RcppExport SEXP _collapse_setattr_clp(SEXP xSEXP, SEXP aSEXP, SEXP vSEXP);
// duplAttributes
RcppExport SEXP _collapse_duplAttributes(SEXP xSEXP, SEXP ySEXP);
// duplattributes
RcppExport SEXP _collapse_duplattributes(SEXP xSEXP, SEXP ySEXP);
// cond_duplAttributes
RcppExport SEXP _collapse_cond_duplAttributes(SEXP xSEXP, SEXP ySEXP);
// cond_duplattributes
RcppExport SEXP _collapse_cond_duplattributes(SEXP xSEXP, SEXP ySEXP);
// lassignCpp
RcppExport SEXP _collapse_lassignCpp(SEXP xSEXP, SEXP sSEXP, SEXP rowsSEXP, SEXP fillSEXP);
// fanyNAint
// RcppExport SEXP _collapse_fanyNAint(SEXP xSEXP);
// groups2GRP
RcppExport SEXP _collapse_groups2GRPCpp(SEXP xSEXP, SEXP lxSEXP, SEXP gsSEXP);
