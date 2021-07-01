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

void setattr(SEXP x, SEXP a, SEXP v) {
  Rf_setAttrib(x, a, v);
}

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

SEXP copyMostAttributes(SEXP x, SEXP y) {
  Rf_copyMostAttrib(y, x);
  return x;
}


SEXP CsetAttrib(SEXP object, SEXP a) {
  int il = isList(object);
  SEXP res = il ? PROTECT(Rf_shallow_duplicate(object)) : object; // needed, otherwise error !!
  SET_ATTRIB(res, PROTECT(Rf_coerceVector(a, LISTSXP)));
  Rf_classgets(res, Rf_getAttrib(res, R_ClassSymbol));
  UNPROTECT(il+1);
  return res;
}

SEXP CcopyAttrib(SEXP to, SEXP from) {
  int il = isList(to);
  SEXP res = il ? PROTECT(Rf_shallow_duplicate(to)) : to;
  DUPLICATE_ATTRIB(res, from);
  UNPROTECT(il);
  return res;
}


SEXP CcopyMostAttrib(SEXP to, SEXP from) {
  int il = isList(to);
  SEXP res = il ? PROTECT(Rf_shallow_duplicate(to)) : to;
  Rf_copyMostAttrib(from, res);
  UNPROTECT(il);
  return res;
}

// Warning message: In .Call(C_duplattributes, x, y) : converting NULL pointer to R NULL
void duplattributes(SEXP x, SEXP y) {
  DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
  Rf_classgets(x, Rf_getAttrib(y, R_ClassSymbol)); // This solves the warning message !!
  // just to return R_NilValue; and the SEXP... retrns NULL anyway
}

SEXP cond_duplAttributes(SEXP x, SEXP y) {
  if(TYPEOF(x) == TYPEOF(y)) DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
  return x;
}

// not used !!
// void cond_duplattributes(SEXP x, SEXP y) {
//  if(TYPEOF(x) == TYPEOF(y)) DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
// }
