#include "collapse_c.h"

void matCopyAttr(SEXP out, SEXP x, SEXP Rdrop, int ng) {
  SEXP dn = getAttrib(x, R_DimNamesSymbol);
  SEXP cn = isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1); // PROTECT ??
  if(ng == 0 && asLogical(Rdrop)) {
    if(length(cn)) setAttrib(out, R_NamesSymbol, cn);
  } else {
    SEXP dim, dn;
    dim = PROTECT(duplicate(getAttrib(x, R_DimSymbol)));
    INTEGER(dim)[0] = ng == 0 ? 1 : ng;
    dimgets(out, dim);
    if(length(cn)) {
      setAttrib(out, R_DimNamesSymbol, dn = allocVector(VECSXP, 2)); // Protected by out..
      SET_VECTOR_ELT(dn, 0, R_NilValue);
      SET_VECTOR_ELT(dn, 1, cn);
    }
    if(!isObject(x)) copyMostAttrib(x, out);
    UNPROTECT(1);
  }
}

void DFcopyAttr(SEXP out, SEXP x, int ng) {
  DUPLICATE_ATTRIB(out, x);
  if(ng == 0) {
    setAttrib(out, R_RowNamesSymbol, ScalarInteger(1));
  } else {
    SEXP rn;
    setAttrib(out, R_RowNamesSymbol, rn = allocVector(INTSXP, 2));
    INTEGER(rn)[0] = NA_INTEGER;
    INTEGER(rn)[1] = -ng;
  }
}

// Faster than rep_len(value, n) and slightly faster than matrix(value, n) (which in turn is faster than rep_len)...
SEXP falloc(SEXP value, SEXP n) {
  int l = asInteger(n), tval = TYPEOF(value);
  if(length(value) > 1) error("Must supply a single value to alloc()");
  SEXP out = PROTECT(allocVector(tval, l));
  switch(tval) {
    case INTSXP:
    case LGLSXP: {
      int val = asInteger(value), *pout = INTEGER(out);
      if(val == 0) memset(pout, 0, l*sizeof(int));
      else for(int i = 0; i != l; ++i) pout[i] = val;
      break;
    }
    case REALSXP: {
      double val = asReal(value), *pout = REAL(out);
      if(val == 0.0) memset(pout, 0.0, l*sizeof(double));
      else for(int i = 0; i != l; ++i) pout[i] = val;
      break;
    }
    case STRSXP: {
      SEXP val = asChar(value), *pout = STRING_PTR(out);
      for(int i = 0; i != l; ++i) pout[i] = val;
      break;
    }
    case VECSXP: {
      SEXP *pout = SEXPPTR(out);
      for(int i = 0; i != l; ++i) pout[i] = value;
      break;
    }
    default: error("Not supportd SEXP Type in alloc()");
  }
  copyMostAttrib(value, out);
  UNPROTECT(1);
  return out;
}


SEXP groups2GRP(SEXP x, SEXP lx, SEXP gs) {
  int l = length(x);
  SEXP out = PROTECT(allocVector(INTSXP, asInteger(lx)));
  int *pout = INTEGER(out)-1, *pgs = INTEGER(gs);
  // SEXP *px = VECTOR_PTR(x); // -> Depreciated interface: https://github.com/hadley/r-internals/blob/ea892fa79bbffe961e78dbe9c90ce4ca3bf2d9bc/vectors.md
  // Matt Dowle Commented:
  // VECTOR_PTR does exist but returns 'not safe to return vector pointer' when USE_RINTERNALS is not defined.
  // VECTOR_DATA and LIST_POINTER exist too but call VECTOR_PTR. All are clearly not intended to be used by packages.
  // The concern is overhead inside VECTOR_ELT() biting when called repetitively in a loop like we do here. That's why
  // we take the R API (INTEGER()[i], REAL()[i], etc) outside loops for the simple types even when not parallel. For this
  // type list case (VECSXP) it might be that some items are ALTREP for example, so we really should use the heavier
  // _ELT accessor (VECTOR_ELT) inside the loop in this case.
  SEXP *px = SEXPPTR(x);

  for(int j = l; j--; ) { // This can go in any direction..
    // SEXP column = VECTOR_ELT(x, j);
    int *pcolumn = INTEGER(px[j]), jp = j+1;
    for(int i = pgs[j]; i--; ) pout[pcolumn[i]] = jp; // This can go in any direction...
  }
  UNPROTECT(1);
  return out;
}

// Note: Only supports numeric data!!!!
SEXP lassign(SEXP x, SEXP s, SEXP rows, SEXP fill) {
  int l = length(x), tr = TYPEOF(rows), ss = asInteger(s), rs = LENGTH(rows);
  SEXP out = PROTECT(allocVector(VECSXP, l));
  // SEXP *px = VECTOR_PTR(x); // -> Depreciated interface: https://github.com/hadley/r-internals/blob/ea892fa79bbffe961e78dbe9c90ce4ca3bf2d9bc/vectors.md
  SEXP *px = SEXPPTR(x);
  double dfill = asReal(fill);

  if(tr == INTSXP) {
    int *rowsv = INTEGER(rows); //, vs = ss * sizeof(double);
    for(int j = l; j--; ) {
      SEXP column = px[j]; // VECTOR_ELT(x, j);
      if(length(column) != rs) error("length(rows) must match nrow(x)");
      SEXP outj;
      SET_VECTOR_ELT(out, j, outj = allocVector(REALSXP, ss));
      double *pcolumn = REAL(column), *poutj = REAL(outj);
      // memset(poutj, dfill, vs); // cannot memset missing values... can only memset 0
      for(int i = ss; i--; ) poutj[i] = dfill;
      for(int i = 0; i != rs; ++i) poutj[rowsv[i]-1] = pcolumn[i];
      SHALLOW_DUPLICATE_ATTRIB(outj, column);
    }
  } else if(tr == LGLSXP) {
    int *rowsv = LOGICAL(rows);
    if(ss != rs) error("length(rows) must match length(s) if rows is a logical vector");
    for(int j = l; j--; ) {
      SEXP column = px[j]; // VECTOR_ELT(x, j);
      SEXP outj;
      SET_VECTOR_ELT(out, j, outj = allocVector(REALSXP, ss));
      double *pcolumn = REAL(column), *poutj = REAL(outj);
      for(int i = 0, k = 0; i != rs; ++i) poutj[i] = rowsv[i] ? pcolumn[k++] : dfill;
      SHALLOW_DUPLICATE_ATTRIB(outj, column);
    }
  } else error("rows must be positive integers or a logical vector");
  DUPLICATE_ATTRIB(out, x);
  UNPROTECT(1);
  return out;
}


SEXP Cna_rm(SEXP x) {
  const int n = LENGTH(x);
  if (n < 1) return x;
  int k = 0;
  switch(TYPEOF(x)) {
  case LGLSXP:
  case INTSXP: {
    const int *xd = INTEGER(x);
    for (int i = 0; i != n; ++i) if(xd[i] == NA_INTEGER) ++k;
    if(k == 0) return x;
    SEXP out = PROTECT(allocVector(TYPEOF(x), n - k));
    int *pout = INTEGER(out);
    k = 0;
    for (int i = 0; i != n; ++i) if(xd[i] != NA_INTEGER) pout[k++] = xd[i];
    copyMostAttrib(x, out);
    UNPROTECT(1);
    return out;
  }
  case REALSXP: { // What about integer64??
    const double *xd = REAL(x);
    for (int i = 0; i != n; ++i) if(ISNAN(xd[i])) ++k;
    if(k == 0) return x;
    SEXP out = PROTECT(allocVector(REALSXP, n - k));
    double *pout = REAL(out);
    k = 0;
    for (int i = 0; i != n; ++i) if(NISNAN(xd[i])) pout[k++] = xd[i]; // using xd[i] == xd[i] is not faster !!
    copyMostAttrib(x, out);
    UNPROTECT(1);
    return out;
  }
  case STRSXP: {
    const SEXP *xd = STRING_PTR(x);
    for (int i = 0; i != n; ++i) if(xd[i] == NA_STRING) ++k;
    if(k == 0) return x;
    SEXP out = PROTECT(allocVector(STRSXP, n - k));
    SEXP *pout = STRING_PTR(out);
    k = 0;
    for (int i = 0; i != n; ++i) if(xd[i] != NA_STRING) pout[k++] = xd[i];
    copyMostAttrib(x, out);
    UNPROTECT(1);
    return out;
  }}
  error("Unsupported type '%s' passed to na_rm()", type2char(TYPEOF(x)));
}


