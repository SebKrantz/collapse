#include <R.h>
#include <Rinternals.h>

#define SEXPPTR(x) ((SEXP *)DATAPTR(x))  // to avoid overhead of looped VECTOR_ELT

SEXP fnobsC(SEXP x, SEXP Rng, SEXP g) {
  int l = length(x), ng = asInteger(Rng);
  if (ng == 0) {
    int n = 0;
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        for(int i = 0; i != l; ++i) if(px[i] == px[i]) ++n;
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        for(int i = 0; i != l; ++i) if(px[i] != NA_INTEGER) ++n;
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x);
        for(int i = 0; i != l; ++i) if(px[i] != NA_STRING) ++n;
        break;
      }
      case VECSXP: {
        SEXP *px = SEXPPTR(x);
        for(int i = 0; i != l; ++i) if(length(px[i])) ++n;
        break;
      }
      default: error("Unsupported SEXP type");
    }
    return ScalarInteger(n);
  } else { // with groups
    if(length(g) != l) error("length(g) must match NROW(X)");
    SEXP n = PROTECT(allocVector(INTSXP, ng));
    int *pn = INTEGER(n), *pg = INTEGER(g);
    memset(pn, 0, sizeof(int) * ng); --pn;
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        for(int i = 0; i != l; ++i) if(px[i] == px[i]) ++pn[pg[i]];
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        for(int i = 0; i != l; ++i) if(px[i] != NA_INTEGER) ++pn[pg[i]];
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x);
        for(int i = 0; i != l; ++i) if(px[i] != NA_STRING) ++pn[pg[i]];
        break;
      }
      case VECSXP: {
        SEXP *px = SEXPPTR(x);
        for(int i = 0; i != l; ++i) if(length(px[i])) ++pn[pg[i]];
        break;
      }
      default: error("Unsupported SEXP type");
    }
    if(!isObject(x)) {
      copyMostAttrib(x, n); // SHALLOW_DUPLICATE_ATTRIB(n, x);
    } else {
      SEXP sym_label = PROTECT(install("label")); // PROTECT ??
      setAttrib(n, sym_label, getAttrib(x, sym_label));
      UNPROTECT(1);
    }
    UNPROTECT(1);
    return n;
  }
}


SEXP fnobsmC(SEXP x, SEXP Rng, SEXP g, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol); // protect ??
  if(isNull(dim)) error("x is not a matrix");
  int ng = asInteger(Rng), drop = asInteger(Rdrop),
    l = INTEGER(dim)[0], col = INTEGER(dim)[1];

  if (ng == 0) {
    SEXP n = PROTECT(allocVector(INTSXP, col));
    int *pn = INTEGER(n);
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        for(int j = 0; j != col; ++j) {
           int nj = 0, end = l * j + l;
           for(int i = l * j; i != end; ++i) if(px[i] == px[i]) ++nj;
           pn[j] = nj;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        for(int j = 0; j != col; ++j) {
          int nj = 0, end = l * j + l;
          for(int i = l * j; i != end; ++i) if(px[i] == NA_INTEGER) ++nj;
          pn[j] = nj;
        }
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x);
        for(int j = 0; j != col; ++j) {
          int nj = 0, end = l * j + l;
          for(int i = l * j; i != end; ++i) if(px[i] == NA_STRING) ++nj;
          pn[j] = nj;
        }
        break;
      }
      default: error("Unsupported SEXP type");
    }
    SEXP dimnames = getAttrib(x, R_DimNamesSymbol); // protect ??
    if(drop && !isNull(dimnames)) setAttrib(n, R_NamesSymbol, VECTOR_ELT(dimnames, 1));
    else {
      SEXP dim2 = PROTECT(duplicate(dim)); // fastest ??
      INTEGER(dim2)[0] = 1;
      dimgets(n, dim2);
      UNPROTECT(1);
      if(!isNull(dimnames)) {
        SEXP dn;
        setAttrib(n, R_DimNamesSymbol, dn = allocVector(VECSXP, 2)); // fastest ??
        SET_VECTOR_ELT(dn, 0, R_NilValue);
        SET_VECTOR_ELT(dn, 1, VECTOR_ELT(dimnames, 1));
      }
      if(!isObject(x)) copyMostAttrib(x, n);
    }
    UNPROTECT(1);
    return n;
  } else { // with groups
    if(length(g) != l) error("length(g) must match NROW(X)");
    SEXP n = PROTECT(allocMatrix(INTSXP, ng, col));
    memset(INTEGER(n), 0, sizeof(int) * ng * col);
    int *pg = INTEGER(g);
    switch(TYPEOF(x)) {
      case REALSXP: {
        for(int j = 0; j != col; ++j) {
          int *pn = INTEGER(n) + j * ng - 1;
          double *px = REAL(x) + j * l;
          for(int i = 0; i != l; ++i) if(px[i] == px[i]) ++pn[pg[i]];
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        for(int j = 0; j != col; ++j) {
          int *pn = INTEGER(n) + j * ng - 1, *px = INTEGER(x) + j * l;
          for(int i = 0; i != l; ++i) if(px[i] == NA_INTEGER) ++pn[pg[i]];
        }
        break;
      }
      case STRSXP: {
        for(int j = 0; j != col; ++j) {
          int *pn = INTEGER(n) + j * ng - 1;
          SEXP *px = STRING_PTR(x) + j * l;
          for(int i = 0; i != l; ++i) if(px[i] == NA_STRING) ++pn[pg[i]];
        }
        break;
      }
      default: error("Unsupported SEXP type");
    }
    SEXP dimnames = getAttrib(x, R_DimNamesSymbol); // protect ??
    if(!isNull(dimnames)) {
      SEXP dn;
      setAttrib(n, R_DimNamesSymbol, dn = allocVector(VECSXP, 2)); // fastest ??
      SET_VECTOR_ELT(dn, 0, R_NilValue);
      SET_VECTOR_ELT(dn, 1, VECTOR_ELT(dimnames, 1));
    }
    if(!isObject(x)) copyMostAttrib(x, n);
    return n;
  }
}


SEXP fnobslC(SEXP x, SEXP Rng, SEXP g, SEXP Rdrop) {
  int l = length(x), ng = asInteger(Rng);
  // if(l < 1) return x;
  if(asLogical(Rdrop) && ng == 0) {
    SEXP out = PROTECT(allocVector(INTSXP, l)), *px = SEXPPTR(x);
    int *pout = INTEGER(out);
    for(int j = 0; j != l; ++j) pout[j] = INTEGER(fnobsC(px[j], Rng, g))[0];
    namesgets(out, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return out;
  } else {
    SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out), *px = SEXPPTR(x);
    for(int j = 0; j != l; ++j) pout[j] = fnobsC(px[j], Rng, g);
    DUPLICATE_ATTRIB(out, x);
    if(ng == 0) {
      setAttrib(out, R_RowNamesSymbol, ScalarInteger(1));
    } else {
      SEXP rn;
      setAttrib(out, R_RowNamesSymbol, rn = allocVector(INTSXP, 2)); // PROTECT ?? -> protexted by out
      INTEGER(rn)[0] = NA_INTEGER;
      INTEGER(rn)[1] = -ng;
    }
    UNPROTECT(1);
    return out;
  }
}
