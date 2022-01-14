#include <R.h>
#include <Rinternals.h>


SEXP setAttributes(SEXP x, SEXP a) {
  SET_ATTRIB(x, Rf_coerceVector(a, LISTSXP));
  Rf_classgets(x, Rf_getAttrib(x, R_ClassSymbol)); // forcing class after attribute copy !!
  return x;
}

void setattributes(SEXP x, SEXP a) {
  SET_ATTRIB(x, Rf_coerceVector(a, LISTSXP));
  // SET_OBJECT(x, TYPEOF(x)); // if(OBJECT(a))  // This does not work with ts-matrices! could also make compatible with S4 objects !
  Rf_classgets(x, Rf_getAttrib(x, R_ClassSymbol));
}

// not used !
// SEXP setAttr(SEXP x, SEXP a, SEXP v) {
//  Rf_setAttrib(x, a, v);
//  return x;
// }

// void setattr(SEXP x, SEXP a, SEXP v) {
//  Rf_setAttrib(x, a, v);
// }

SEXP duplAttributes(SEXP x, SEXP y) { // also look at data.table's keepattributes ...
  DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
  return x;
}

// R_duplicate_attr -> deep copy only of attributes -> expensive if attributes are large !
// Rf_lazy_duplicate -> duplicate on modify -> but modifies object in global environment !
// Rf_shallow_duplicate -> only duplicate pointer? -> best !!

// No speed improvement to attr<- (same slow performance for data.frame 'row.names')
// SEXP CsetAttr(SEXP object, SEXP a, SEXP v) {
//   SEXP res = Rf_shallow_duplicate(object);
//   Rf_setAttrib(res, a, v);
//   return res;
// }

// Attribute Handling - 4 Situations:
// 1 - x is classed (factor, date, time series), xAG is not classed. i.e. vector of fnobs, fmean etc.
//    -> Sallow replacing, removing class and levels attributes from x, discard attributes of xAG (if any)
//    -> or (if type matches i.e. double for date or time series), copy attributes of x unless x is a factor
// 2 - x is not classed, xAG is classed (factor, date, time series). - an unusual situation should not occurr - copy attributes of xAG, discard attributes of x
// 3 - xAG and x are classed - same as above, keep attributes of xAG, discard attributes of x
// 4 - neither x nor xAG are classed - preserve attributes of x, discard attributes of xAG (if any)
//

// if(Rf_isObject(xAG)) DUPLICATE_ATTRIB(out, xAG);
// else if(!Rf_isObject(x) || (tx == txAG && !Rf_isFactor(x))) DUPLICATE_ATTRIB(out, x);
// else {
//   SHALLOW_DUPLICATE_ATTRIB(out, x);
//   Rf_classgets(out, R_NilValue); // OK !
//   Rf_setAttrib(out, R_LevelsSymbol, R_NilValue); // if(Rf_isFactor(x)) ? faster ?
// }

// Can think further about this! but this solution appears acceptable...

SEXP copyMostAttributes(SEXP x, SEXP y) {
  int tx = TYPEOF(x);
  // -> This is about the best we can do: unlist() does not preserve dates, and we don't want to create malformed factors
  if(tx == TYPEOF(y) && (tx != INTSXP || OBJECT(x) == OBJECT(y))) {
    Rf_copyMostAttrib(y, x);
    return x;
  }
  // In any case we can preserve variable labels..
  SEXP sym_label = PROTECT(install("label"));
  SEXP lab = Rf_getAttrib(y, sym_label);
  if(TYPEOF(lab) != NILSXP) Rf_setAttrib(x, sym_label, lab);
  UNPROTECT(1);
  return x;
}


SEXP CsetAttrib(SEXP object, SEXP a) {
  int il = TYPEOF(object) == VECSXP;
  SEXP res = il ? PROTECT(Rf_shallow_duplicate(object)) : object; // needed, otherwise error !!
  SET_ATTRIB(res, PROTECT(Rf_coerceVector(a, LISTSXP)));
  Rf_classgets(res, Rf_getAttrib(res, R_ClassSymbol));
  UNPROTECT(il+1);
  return res;
}

SEXP CcopyAttrib(SEXP to, SEXP from) {
  int il = TYPEOF(to) == VECSXP;
  SEXP res = il ? PROTECT(Rf_shallow_duplicate(to)) : to;
  DUPLICATE_ATTRIB(res, from);
  UNPROTECT(il);
  return res;
}


SEXP CcopyMostAttrib(SEXP to, SEXP from) {
  int il = TYPEOF(to) == VECSXP;
  SEXP res = il ? PROTECT(Rf_shallow_duplicate(to)) : to;
  Rf_copyMostAttrib(from, res);
  if(il && isFrame(from) && length(VECTOR_ELT(to, 0)) != length(VECTOR_ELT(from, 0))) {
    Rf_setAttrib(res, R_RowNamesSymbol, Rf_getAttrib(to, R_RowNamesSymbol));
  }
  UNPROTECT(il);
  return res;
}

// No longer needed...
// Warning message: In .Call(C_duplattributes, x, y) : converting NULL pointer to R NULL
// void duplattributes(SEXP x, SEXP y) {
//  DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
//  Rf_classgets(x, Rf_getAttrib(y, R_ClassSymbol)); // This solves the warning message !!
  // just to return R_NilValue; and the SEXP... retrns NULL anyway
// }

// No longer needed... using copyMostAttributes instead
// SEXP cond_duplAttributes(SEXP x, SEXP y) {
//  if(TYPEOF(x) == TYPEOF(y)) DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
//  return x;
// }

// not used !!
// void cond_duplattributes(SEXP x, SEXP y) {
//  if(TYPEOF(x) == TYPEOF(y)) DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
// }
