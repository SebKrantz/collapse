#include <Rcpp.h>
using namespace Rcpp;

// BWCpp
RcppExport SEXP _collapsedev17_BWCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP thetaSEXP, SEXP set_meanSEXP, SEXP BSEXP, SEXP fillSEXP);
// BWmCpp
RcppExport SEXP _collapsedev17_BWmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP thetaSEXP, SEXP set_meanSEXP, SEXP BSEXP, SEXP fillSEXP);
// BWlCpp
RcppExport SEXP _collapsedev17_BWlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP thetaSEXP, SEXP set_meanSEXP, SEXP BSEXP, SEXP fillSEXP);
// TRACpp
RcppExport SEXP _collapsedev17_TRACpp(SEXP xSEXP, SEXP xAGSEXP, SEXP gSEXP, SEXP retSEXP);
// TRAmCpp
RcppExport SEXP _collapsedev17_TRAmCpp(SEXP xSEXP, SEXP xAGSEXP, SEXP gSEXP, SEXP retSEXP);
// TRAlCpp
RcppExport SEXP _collapsedev17_TRAlCpp(SEXP xSEXP, SEXP xAGSEXP, SEXP gSEXP, SEXP retSEXP);
// fndistinctCpp
RcppExport SEXP _collapsedev17_fndistinctCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP narmSEXP);
// fndistinctlCpp
RcppExport SEXP _collapsedev17_fndistinctlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fndistinctmCpp
RcppExport SEXP _collapsedev17_fndistinctmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP narmSEXP, SEXP dropSEXP);
// pwnobsmCpp
RcppExport SEXP _collapsedev17_pwnobsmCpp(SEXP xSEXP);
// fnobsCpp
// RcppExport SEXP _collapsedev17_fnobsCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP);
// fnobsmCpp
// RcppExport SEXP _collapsedev17_fnobsmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP dropSEXP);
// fnobslCpp
// RcppExport SEXP _collapsedev17_fnobslCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP dropSEXP);
// varyingCpp
RcppExport SEXP _collapsedev17_varyingCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP any_groupSEXP);
// varyingmCpp
RcppExport SEXP _collapsedev17_varyingmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP any_groupSEXP, SEXP dropSEXP);
// varyinglCpp
RcppExport SEXP _collapsedev17_varyinglCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP any_groupSEXP, SEXP dropSEXP);
// fbstatsCpp
RcppExport SEXP _collapsedev17_fbstatsCpp(SEXP xSEXP, SEXP extSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP npgSEXP, SEXP pgSEXP, SEXP wSEXP, SEXP arraySEXP, SEXP setnSEXP, SEXP gnSEXP);
// fbstatsmCpp
RcppExport SEXP _collapsedev17_fbstatsmCpp(SEXP xSEXP, SEXP extSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP npgSEXP, SEXP pgSEXP, SEXP wSEXP, SEXP arraySEXP, SEXP gnSEXP);
// fbstatslCpp
RcppExport SEXP _collapsedev17_fbstatslCpp(SEXP xSEXP, SEXP extSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP npgSEXP, SEXP pgSEXP, SEXP wSEXP, SEXP arraySEXP, SEXP gnSEXP);
// ffirstCpp
// RcppExport SEXP _collapsedev17_ffirstCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP);
// ffirstmCpp
// RcppExport SEXP _collapsedev17_ffirstmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP dropSEXP);
// ffirstlCpp
// RcppExport SEXP _collapsedev17_ffirstlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP);
// fdiffgrowthCpp
RcppExport SEXP _collapsedev17_fdiffgrowthCpp(SEXP xSEXP, SEXP nSEXP, SEXP diffSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP retSEXP, SEXP rhoSEXP, SEXP namesSEXP, SEXP powerSEXP);
// fdiffgrowthmCpp
RcppExport SEXP _collapsedev17_fdiffgrowthmCpp(SEXP xSEXP, SEXP nSEXP, SEXP diffSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP retSEXP, SEXP rhoSEXP, SEXP namesSEXP, SEXP powerSEXP);
// fdiffgrowthlCpp
RcppExport SEXP _collapsedev17_fdiffgrowthlCpp(SEXP xSEXP, SEXP nSEXP, SEXP diffSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP tSEXP, SEXP retSEXP, SEXP rhoSEXP, SEXP namesSEXP, SEXP powerSEXP);
// flagleadCpp
RcppExport SEXP _collapsedev17_flagleadCpp(SEXP xSEXP, SEXP nSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP tSEXP, SEXP namesSEXP);
// flagleadmCpp
RcppExport SEXP _collapsedev17_flagleadmCpp(SEXP xSEXP, SEXP nSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP tSEXP, SEXP namesSEXP);
// flagleadlCpp
RcppExport SEXP _collapsedev17_flagleadlCpp(SEXP xSEXP, SEXP nSEXP, SEXP fillSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP tSEXP, SEXP namesSEXP);
// flastCpp
// RcppExport SEXP _collapsedev17_flastCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP);
// flastmCpp
// RcppExport SEXP _collapsedev17_flastmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP dropSEXP);
// flastlCpp
// RcppExport SEXP _collapsedev17_flastlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP);
// fminmaxCpp
// RcppExport SEXP _collapsedev17_fminmaxCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP retSEXP);
// fminmaxmCpp
// RcppExport SEXP _collapsedev17_fminmaxmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP dropSEXP, SEXP retSEXP);
// fminmaxlCpp
// RcppExport SEXP _collapsedev17_fminmaxlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP narmSEXP, SEXP dropSEXP, SEXP retSEXP);
// fmeanCpp
RcppExport SEXP _collapsedev17_fmeanCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP);
// fmeanmCpp
RcppExport SEXP _collapsedev17_fmeanmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fmeanlCpp
RcppExport SEXP _collapsedev17_fmeanlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fmedianCpp
// RcppExport SEXP _collapsedev17_fmedianCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP);
// fmedianmCpp
// RcppExport SEXP _collapsedev17_fmedianmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fmedianlCpp
// RcppExport SEXP _collapsedev17_fmedianlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fnthCpp
RcppExport SEXP _collapsedev17_fnthCpp(SEXP xSEXP, SEXP QSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP retSEXP);
// fnthmCpp
RcppExport SEXP _collapsedev17_fnthmCpp(SEXP xSEXP, SEXP QSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP, SEXP retSEXP);
// fnthlCpp
RcppExport SEXP _collapsedev17_fnthlCpp(SEXP xSEXP, SEXP QSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP, SEXP retSEXP);
// fmodeCpp
RcppExport SEXP _collapsedev17_fmodeCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP retSEXP);
// fmodelCpp
RcppExport SEXP _collapsedev17_fmodelCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP retSEXP);
// fmodemCpp
RcppExport SEXP _collapsedev17_fmodemCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP, SEXP retSEXP);
// fprodCpp
RcppExport SEXP _collapsedev17_fprodCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP);
// fprodmCpp
RcppExport SEXP _collapsedev17_fprodmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fprodlCpp
RcppExport SEXP _collapsedev17_fprodlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fscaleCpp
RcppExport SEXP _collapsedev17_fscaleCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP set_meanSEXP, SEXP set_sdSEXP);
// fscalemCpp
RcppExport SEXP _collapsedev17_fscalemCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP set_meanSEXP, SEXP set_sdSEXP);
// fscalelCpp
RcppExport SEXP _collapsedev17_fscalelCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP set_meanSEXP, SEXP set_sdSEXP);
// fsumCpp
// RcppExport SEXP _collapsedev17_fsumCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP);
// fsummCpp
// RcppExport SEXP _collapsedev17_fsummCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fsumlCpp
// RcppExport SEXP _collapsedev17_fsumlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP dropSEXP);
// fvarsdCpp
RcppExport SEXP _collapsedev17_fvarsdCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP stable_algoSEXP, SEXP sdSEXP);
// fvarsdmCpp
RcppExport SEXP _collapsedev17_fvarsdmCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP stable_algoSEXP, SEXP sdSEXP, SEXP dropSEXP);
// fvarsdlCpp
RcppExport SEXP _collapsedev17_fvarsdlCpp(SEXP xSEXP, SEXP ngSEXP, SEXP gSEXP, SEXP gsSEXP, SEXP wSEXP, SEXP narmSEXP, SEXP stable_algoSEXP, SEXP sdSEXP, SEXP dropSEXP);
// mrtl
RcppExport SEXP _collapsedev17_mrtl(SEXP XSEXP, SEXP namesSEXP, SEXP retSEXP);
// mctl
RcppExport SEXP _collapsedev17_mctl(SEXP XSEXP, SEXP namesSEXP, SEXP retSEXP);
// psmatCpp
RcppExport SEXP _collapsedev17_psmatCpp(SEXP xSEXP, SEXP gSEXP, SEXP tSEXP, SEXP transposeSEXP);
// qFCpp
RcppExport SEXP _collapsedev17_qFCpp(SEXP xSEXP, SEXP sortSEXP, SEXP orderedSEXP, SEXP na_excludeSEXP, SEXP keep_attrSEXP, SEXP retSEXP);
// qGCpp
// RcppExport SEXP _collapsedev17_qGCpp(SEXP xSEXP, SEXP sortSEXP, SEXP orderedSEXP, SEXP na_excludeSEXP, SEXP retgrpSEXP);
// funiqueCpp
RcppExport SEXP _collapsedev17_funiqueCpp(SEXP xSEXP, SEXP sortSEXP);
// fdroplevelsCpp
RcppExport SEXP _collapsedev17_fdroplevelsCpp(SEXP xSEXP, SEXP check_NASEXP);

// setAttributes
// SEXP _collapsedev17_setAttributes(SEXP xSEXP, SEXP aSEXP);
// setattributes
// SEXP _collapsedev17_setattributes(SEXP xSEXP, SEXP aSEXP);
// setAttr
// RcppExport SEXP _collapsedev17_setAttr(SEXP xSEXP, SEXP aSEXP, SEXP vSEXP);
// setattr
// SEXP _collapsedev17_setattr(SEXP xSEXP, SEXP aSEXP, SEXP vSEXP);
// duplAttributes
// SEXP _collapsedev17_duplAttributes(SEXP xSEXP, SEXP ySEXP);
// duplattributes
// SEXP _collapsedev17_duplattributes(SEXP xSEXP, SEXP ySEXP);
// cond_duplAttributes
// SEXP _collapsedev17_cond_duplAttributes(SEXP xSEXP, SEXP ySEXP);
// cond_duplattributes
// RcppExport SEXP _collapsedev17_cond_duplattributes(SEXP xSEXP, SEXP ySEXP);
// seqid
RcppExport SEXP _collapsedev17_seqid(SEXP xSEXP, SEXP oSEXP, SEXP delSEXP, SEXP startSEXP, SEXP na_skipSEXP, SEXP skip_seqSEXP, SEXP check_oSEXP);
// groupid
RcppExport SEXP _collapsedev17_groupid(SEXP xSEXP, SEXP oSEXP, SEXP startSEXP, SEXP na_skipSEXP, SEXP check_oSEXP);
