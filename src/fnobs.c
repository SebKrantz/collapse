#include "collapse_c.h"

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
        SEXP *px = SEXPPTR(x);
        for(int i = 0; i != l; ++i) if(px[i] != NA_STRING) ++n;
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
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
        SEXP *px = SEXPPTR(x);
        for(int i = 0; i != l; ++i) if(px[i] != NA_STRING) ++pn[pg[i]];
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
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
  int ng = asInteger(Rng), l = INTEGER(dim)[0], col = INTEGER(dim)[1];

  SEXP n = PROTECT(allocVector(INTSXP, ng == 0 ? col : ng * col));
  int *pn = INTEGER(n);

  if (ng == 0) {
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        for(int j = 0; j != col; ++j) {
           int nj = 0, end = l * j + l;
           for(int i = l * j; i != end; ++i) if(NISNAN(px[i])) ++nj;
           pn[j] = nj;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        for(int j = 0; j != col; ++j) {
          int nj = 0, end = l * j + l;
          for(int i = l * j; i != end; ++i) if(px[i] != NA_INTEGER) ++nj;
          pn[j] = nj;
        }
        break;
      }
      case STRSXP: {
        SEXP *px = SEXPPTR(x);
        for(int j = 0; j != col; ++j) {
          int nj = 0, end = l * j + l;
          for(int i = l * j; i != end; ++i) if(px[i] != NA_STRING) ++nj;
          pn[j] = nj;
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        for(int j = 0; j != col; ++j) {
          int nj = 0, end = l * j + l;
          for(int i = l * j; i != end; ++i) if(length(px[i])) ++nj;
          pn[j] = nj;
        }
        break;
      }
      default: error("Unsupported SEXP type");
    }
  } else { // with groups
    if(length(g) != l) error("length(g) must match NROW(X)");
    memset(pn, 0, sizeof(int) * ng * col);
    pn -= ng + 1;
    int *pg = INTEGER(g);
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x)-l;
        for(int j = 0; j != col; ++j) {
          pn += ng; px += l;
          for(int i = 0; i != l; ++i) if(NISNAN(px[i])) ++pn[pg[i]];
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x)-l;
        for(int j = 0; j != col; ++j) {
          pn += ng; px += l;
          for(int i = 0; i != l; ++i) if(px[i] != NA_INTEGER) ++pn[pg[i]];
        }
        break;
      }
      case STRSXP: {
        SEXP *px = SEXPPTR(x)-l;
        for(int j = 0; j != col; ++j) {
          pn += ng; px += l;
          for(int i = 0; i != l; ++i) if(px[i] != NA_STRING) ++pn[pg[i]];
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x)-l;
        for(int j = 0; j != col; ++j) {
          pn += ng; px += l;
          for(int i = 0; i != l; ++i) if(length(px[i])) ++pn[pg[i]];
        }
        break;
      }
      default: error("Unsupported SEXP type");
    }
  }
  matCopyAttr(n, x, Rdrop, ng);
  UNPROTECT(1);
  return n;
}


SEXP fnobslC(SEXP x, SEXP Rng, SEXP g, SEXP Rdrop) {
  int l = length(x), ng = asInteger(Rng);
  if(l < 1) return x;
  if(asLogical(Rdrop) && ng == 0) {
    SEXP out = PROTECT(allocVector(INTSXP, l));
    const SEXP *px = SEXPPTR_RO(x);
    int *pout = INTEGER(out);
    for(int j = 0; j != l; ++j) pout[j] = INTEGER(fnobsC(px[j], Rng, g))[0];
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return out;
  } else {
    SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out);
    const SEXP *px = SEXPPTR_RO(x);
    for(int j = 0; j != l; ++j) pout[j] = fnobsC(px[j], Rng, g);
    DFcopyAttr(out, x, ng);
    UNPROTECT(1);
    return out;
  }
}
