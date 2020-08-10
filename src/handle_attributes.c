#include <R.h>
#include <Rinternals.h>

SEXP setAttributes(SEXP x, SEXP a) {
  SET_ATTRIB(x, Rf_coerceVector(a, LISTSXP));
  Rf_classgets(x, Rf_getAttrib(x, R_ClassSymbol));
  return x;
}

void setattributes(SEXP x, SEXP a) {
  SET_ATTRIB(x, Rf_coerceVector(a, LISTSXP));
  // SET_OBJECT(x, TYPEOF(x)); // if(OBJECT(a))  // This does not work with ts-matrices! could also make compatible with S4 objects !
  Rf_classgets(x, Rf_getAttrib(x, R_ClassSymbol));
}

// not used !!
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

// Warning message: In .Call(C_duplattributes, x, y) : converting NULL pointer to R NULL
void duplattributes(SEXP x, SEXP y) {
  DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
  Rf_classgets(x, Rf_getAttrib(y, R_ClassSymbol)); // This solves the warning message !!
}

SEXP cond_duplAttributes(SEXP x, SEXP y) {
  if(TYPEOF(x) == TYPEOF(y)) DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
  return x;
}

// not used !!
// void cond_duplattributes(SEXP x, SEXP y) {
//  if(TYPEOF(x) == TYPEOF(y)) DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
// }
