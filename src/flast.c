#include "collapse_c.h"


SEXP flast_impl(SEXP x, int ng, SEXP g, int narm, int *gl) {

  int l = length(x), tx = TYPEOF(x);
  if (l < 2) return x; // Prevents seqfault for numeric(0) #101
  if (ng == 0) {
    SEXP out = PROTECT(allocVector(tx, 1));
    int j = l-1;
    if(narm) {
      switch(tx) {
        case REALSXP: {
          double *px = REAL(x);
          while(ISNAN(px[j]) && j != 0) --j;
          REAL(out)[0] = px[j];
          break;
        }
        case STRSXP: {
          const SEXP *px = SEXPPTR_RO(x);
          while(px[j] == NA_STRING && j != 0) --j;
          SET_STRING_ELT(out, 0, px[j]);
          break;
        }
        case INTSXP:
        case LGLSXP: {
          int *px = INTEGER(x);
          while(px[j] == NA_INTEGER && j != 0) --j;
          INTEGER(out)[0] = px[j];
          break;
        }
        case VECSXP: {
          const SEXP *px = SEXPPTR_RO(x);
          while(length(px[j]) == 0 && j != 0) --j;
          SET_VECTOR_ELT(out, 0, px[j]);
          break;
        }
        default: error("Unsupported SEXP type!");
      }
    } else {
      switch(tx) {
        case REALSXP: REAL(out)[0] = REAL(x)[l-1];
          break;
        case STRSXP: SET_STRING_ELT(out, 0, STRING_ELT(x, l-1));
          break;
        case INTSXP:
        case LGLSXP: INTEGER(out)[0] = INTEGER(x)[l-1];
          break;
        case VECSXP: SET_VECTOR_ELT(out, 0, VECTOR_ELT(x, l-1));
          break;
        default: error("Unsupported SEXP type!");
      }
    }
    if(ANY_ATTRIB(x) && !(isObject(x) && inherits(x, "ts")))
      copyMostAttrib(x, out);
    if(!isNull(getAttrib(x, R_NamesSymbol)))
      namesgets(out, ScalarString(STRING_ELT(getAttrib(x, R_NamesSymbol), j)));
    UNPROTECT(1);
    return out;
  } else { // with groups
    if(length(g) != l) error("length(g) must match nrow(X)");
    SEXP out = PROTECT(allocVector(tx, ng));
    if(narm) {
      int ngs = 0, *pg = INTEGER(g);
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int i = ng; i--; ) pout[i] = NA_REAL;
        --pout;
        for(int i = l; i--; ) {
          if(NISNAN(px[i])) {
            if(ISNAN(pout[pg[i]])) {
              pout[pg[i]] = px[i];
              if(++ngs == ng) break;
            }
          }
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int i = ng; i--; ) pout[i] = NA_STRING;
        --pout;
        for(int i = l; i--; ) {
          if(px[i] != NA_STRING) {
            if(pout[pg[i]] == NA_STRING) {
              pout[pg[i]] = px[i];
              if(++ngs == ng) break;
            }
          }
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int i = ng; i--; ) pout[i] = NA_INTEGER;
        --pout;
        for(int i = l; i--; ) {
          if(px[i] != NA_INTEGER) {
            if(pout[pg[i]] == NA_INTEGER) {
              pout[pg[i]] = px[i];
              if(++ngs == ng) break;
            }
          }
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int i = ng; i--; ) pout[i] = R_NilValue;
        --pout;
        for(int i = l; i--; ) {
          if(length(px[i])) {
            if(pout[pg[i]] == R_NilValue) {
              pout[pg[i]] = px[i];
              if(++ngs == ng) break;
            }
          }
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    } else {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int i = ng; i--; ) pout[i] = gl[i] == NA_INTEGER ? NA_REAL : px[gl[i]];
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int i = ng; i--; ) pout[i] = gl[i] == NA_INTEGER ? NA_INTEGER : px[gl[i]];
        break;
      }
      case STRSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int i = ng; i--; ) pout[i] = gl[i] == NA_INTEGER ? NA_STRING : px[gl[i]];
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int i = ng; i--; ) pout[i] = gl[i] == NA_INTEGER ? R_NilValue : px[gl[i]];
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    if(ANY_ATTRIB(x) && !(isObject(x) && inherits(x, "ts")))
      copyMostAttrib(x, out);
    UNPROTECT(1);
    return out;
  }
}

SEXP flastC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm) {
  int *pgl, ng = asInteger(Rng), narm = asLogical(Rnarm);
  if(ng == 0 || narm) {
    pgl = &ng;
    return flast_impl(x, ng, g, narm, pgl);
  }
  SEXP gl = PROTECT(allocVector(INTSXP, ng));
  int *pg = INTEGER(g);
  pgl = INTEGER(gl);
  for(int i = ng; i--; ) pgl[i] = NA_INTEGER;
  --pgl;
  for(int i = length(g); i--; ) if(pgl[pg[i]] == NA_INTEGER) pgl[pg[i]] = i;
  SEXP res = flast_impl(x, ng, g, narm, ++pgl);
  UNPROTECT(1);
  return res;
}

SEXP flastlC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm) {
  int l = length(x), *pgl, ng = asInteger(Rng), narm = asLogical(Rnarm), nprotect = 1;
  if(ng > 0 && !narm) {
    SEXP gl = PROTECT(allocVector(INTSXP, ng)); ++nprotect;
    int *pg = INTEGER(g);
    pgl = INTEGER(gl);
    for(int i = ng; i--; ) pgl[i] = NA_INTEGER;
    --pgl;
    for(int i = length(g); i--; ) if(pgl[pg[i]] == NA_INTEGER) pgl[pg[i]] = i;
    ++pgl;
  } else pgl = &l;
  SEXP out = PROTECT(allocVector(VECSXP, l));
  const SEXP *px = SEXPPTR_RO(x);
  for(int j = 0; j != l; ++j) SET_VECTOR_ELT(out, j, flast_impl(px[j], ng, g, narm, pgl));
  DFcopyAttr(out, x, ng);
  UNPROTECT(nprotect);
  return out;
}

// For matrix writing a separate function to increase efficiency.
SEXP flastmC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), ng = asInteger(Rng), narm = asLogical(Rnarm),
    l = INTEGER(dim)[0], col = INTEGER(dim)[1];
  if (l < 2) return x;
  if (ng == 0) {
    SEXP out = PROTECT(allocVector(tx, col));
    if(narm) {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int j = 0, i = l-1; j != col; ++j) {
          while(ISNAN(px[i]) && i != 0) --i;
          pout[j] = px[i]; px += l; i = l-1;
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int j = 0, i = l-1; j != col; ++j) {
          while(px[i] == NA_STRING && i != 0) --i;
          pout[j] = px[i]; px += l; i = l-1;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int j = 0, i = l-1; j != col; ++j) {
          while(px[i] == NA_INTEGER && i != 0) --i;
          pout[j] = px[i]; px += l; i = l-1;
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int j = 0, i = l-1; j != col; ++j) {
          while(length(px[i]) == 0 && i != 0) --i;
          pout[j] = px[i]; px += l; i = l-1;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    } else {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l + l-1];
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l + l-1];
        break;
      }
      case STRSXP:
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int j = 0; j != col; ++j) pout[j] = px[j * l + l-1];
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    }
    matCopyAttr(out, x, Rdrop, ng);
    UNPROTECT(1);
    return out;
  } else { // with groups
    if(length(g) != l) error("length(g) must match nrow(X)");
    SEXP out = PROTECT(allocVector(tx, ng * col));
    int *pg = INTEGER(g);
    if(narm) {
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int i = ng * col; i--; ) pout[i] = NA_REAL;
        --pout;
        for(int j = 0; j != col; ++j) {
          for(int i = l; i--; ) if(NISNAN(px[i]) && ISNAN(pout[pg[i]])) pout[pg[i]] = px[i];
          px += l; pout += ng;
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int i = ng * col; i--; ) pout[i] = NA_STRING;
        --pout;
        for(int j = 0; j != col; ++j) {
          for(int i = l; i--; ) if(px[i] != NA_STRING && pout[pg[i]] == NA_STRING) pout[pg[i]] = px[i];
          px += l; pout += ng;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int i = ng * col; i--; ) pout[i] = NA_INTEGER;
        --pout;
        for(int j = 0; j != col; ++j) {
          for(int i = l; i--; ) if(px[i] != NA_INTEGER && pout[pg[i]] == NA_INTEGER) pout[pg[i]] = px[i];
          px += l; pout += ng;
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int i = ng * col; i--; ) pout[i] = R_NilValue;
        --pout;
        for(int j = 0; j != col; ++j) {
          for(int i = l; i--; ) if(length(px[i]) && pout[pg[i]] != R_NilValue) pout[pg[i]] = px[i];
          px += l; pout += ng;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
    } else {
      SEXP gl = PROTECT(allocVector(INTSXP, ng));
      int *pgl = INTEGER(gl);
      for(int i = ng; i--; ) pgl[i] = NA_INTEGER;
      --pgl;
      for(int i = l; i--; ) if(pgl[pg[i]] == NA_INTEGER) pgl[pg[i]] = i;
      ++pgl;
      switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pout[i] = pgl[i] == NA_INTEGER ? NA_REAL : px[pgl[i]];
          px += l; pout += ng;
        }
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x), *pout = INTEGER(out);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pout[i] = pgl[i] == NA_INTEGER ? NA_INTEGER : px[pgl[i]];
          px += l; pout += ng;
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pout[i] = pgl[i] == NA_INTEGER ? NA_STRING : px[pgl[i]];
          px += l; pout += ng;
        }
        break;
      }
      case VECSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        SEXP *pout = SEXPPTR(out);
        for(int j = 0; j != col; ++j) {
          for(int i = ng; i--; ) pout[i] = pgl[i] == NA_INTEGER ? R_NilValue : px[pgl[i]];
          px += l; pout += ng;
        }
        break;
      }
      default: error("Unsupported SEXP type!");
      }
      UNPROTECT(1);
    }
    matCopyAttr(out, x, Rdrop, ng);
    UNPROTECT(1);
    return out;
  }
}
