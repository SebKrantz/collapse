#include <R.h>
#include <Rinternals.h>

// See https://github.com/wch/r-source/blob/079f863446b5414dd96f3c29d519e4a654146364/src/main/memory.c
// and https://github.com/wch/r-source/blob/80e410a786324e0e472a25481d5dd28db8285330/src/main/attrib.c
// https://github.com/wch/r-source/blob/b6f046826c87fc10ad08acd8858921fa1a58e488/doc/manual/R-ints.texi


SEXP setAttributes(SEXP x, SEXP a) {
  SET_ATTRIB(x, coerceVector(a, LISTSXP));
  classgets(x, getAttrib(x, R_ClassSymbol)); // forcing class after attribute copy !!
  return x;
}

SEXP setattributes(SEXP x, SEXP a) {
  SET_ATTRIB(x, coerceVector(a, LISTSXP));
  // SET_OBJECT(x, TYPEOF(x)); // if(OBJECT(a))  // This does not work with ts-matrices! could also make compatible with S4 objects !
  classgets(x, getAttrib(x, R_ClassSymbol));
  return R_NilValue;
}

// not used !
// SEXP setAttr(SEXP x, SEXP a, SEXP v) {
//  setAttrib(x, a, v);
//  return x;
// }

// void setattr(SEXP x, SEXP a, SEXP v) {
//  setAttrib(x, a, v);
// }

SEXP duplAttributes(SEXP x, SEXP y) { // also look at data.table's keepattributes ...
  SHALLOW_DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
  return x;
}

// R_duplicate_attr -> deep copy only of attributes -> expensive if attributes are large !
// lazy_duplicate -> duplicate on modify -> but modifies object in global environment !
// shallow_duplicate -> only duplicate pointer? -> best !!

// No speed improvement to attr<- (same slow performance for data.frame 'row.names')
// SEXP CsetAttr(SEXP object, SEXP a, SEXP v) {
//   SEXP res = shallow_duplicate(object);
//   setAttrib(res, a, v);
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

// if(isObject(xAG)) SHALLOW_DUPLICATE_ATTRIB(out, xAG);
// else if(!isObject(x) || (tx == txAG && !isFactor(x))) SHALLOW_DUPLICATE_ATTRIB(out, x);
// else {
//   SHALLOW_DUPLICATE_ATTRIB(out, x);
//   classgets(out, R_NilValue); // OK !
//   setAttrib(out, R_LevelsSymbol, R_NilValue); // if(isFactor(x)) ? faster ?
// }

// Can think further about this! but this solution appears acceptable...

SEXP copyMostAttributes(SEXP x, SEXP y) {
  int tx = TYPEOF(x);
  // -> This is about the best we can do: unlist() does not preserve dates, and we don't want to create malformed factors
  // if(TYPEOF(x) == TYPEOF(y) && (OBJECT(x) == OBJECT(y) || (!inherits(y, "factor") && !(length(x) != length(y) && inherits(y, "ts")))))
  if(tx == TYPEOF(y) && (OBJECT(x) == OBJECT(y) || tx != INTSXP || inherits(y, "IDate") || inherits(y, "ITime")) && !(length(x) != length(y) && inherits(y, "ts"))) {
    copyMostAttrib(y, x);
    return x;
  }
  // In any case we can preserve variable labels..
  SEXP sym_label = install("label");
  SEXP lab = getAttrib(y, sym_label);
  if(TYPEOF(lab) != NILSXP) setAttrib(x, sym_label, lab);
  return x;
}


SEXP CsetAttrib(SEXP object, SEXP a) {
  if(TYPEOF(object) == VECSXP) {
    SEXP res = PROTECT(shallow_duplicate(object));
    SET_ATTRIB(res, coerceVector(a, LISTSXP));
    classgets(res, getAttrib(res, R_ClassSymbol));
    UNPROTECT(1);
    return res;
  }
  SEXP res = object;
  SET_ATTRIB(res, coerceVector(a, LISTSXP));
  classgets(res, getAttrib(res, R_ClassSymbol));
  return res;
}

SEXP CcopyAttrib(SEXP to, SEXP from) {
  if(TYPEOF(to) == VECSXP) {
    SEXP res = PROTECT(shallow_duplicate(to));
    SHALLOW_DUPLICATE_ATTRIB(res, from);
    UNPROTECT(1);
    return res;
  }
  SEXP res = to;
  SHALLOW_DUPLICATE_ATTRIB(res, from);
  return res;
}


SEXP CcopyMostAttrib(SEXP to, SEXP from) {
  if(TYPEOF(to) == VECSXP) {
    SEXP res = PROTECT(shallow_duplicate(to));
    copyMostAttrib(from, res);
    if(isFrame(from) && length(VECTOR_ELT(to, 0)) != length(VECTOR_ELT(from, 0)))
       setAttrib(res, R_RowNamesSymbol, getAttrib(to, R_RowNamesSymbol));
    UNPROTECT(1);
    return res;
  }
  SEXP res = to;
  copyMostAttrib(from, res);
  return res;
}

// No longer needed...
// Warning message: In .Call(C_duplattributes, x, y) : converting NULL pointer to R NULL
// void duplattributes(SEXP x, SEXP y) {
//  SHALLOW_DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
//  classgets(x, getAttrib(y, R_ClassSymbol)); // This solves the warning message !!
  // just to return R_NilValue; and the SEXP... retrns NULL anyway
// }

// No longer needed... using copyMostAttributes instead
// SEXP cond_duplAttributes(SEXP x, SEXP y) {
//  if(TYPEOF(x) == TYPEOF(y)) SHALLOW_DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
//  return x;
// }

// not used !!
// void cond_duplattributes(SEXP x, SEXP y) {
//  if(TYPEOF(x) == TYPEOF(y)) SHALLOW_DUPLICATE_ATTRIB(x, y); // SET_ATTRIB(x, ATTRIB(y));
// }
