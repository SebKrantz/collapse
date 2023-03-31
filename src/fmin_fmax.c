#include "collapse_c.h"
// #include <R_ext/Altrep.h>

static double POS_INF = 1.0/0.0;
static double NEG_INF = -1.0/0.0;

void fmin_double_impl(double *pout, double *px, int ng, int *pg, int narm, int l) {
  if(ng == 0) {
    double min;
    if(narm) {
      int j = l-1;
      min = px[j];
      while(ISNAN(min) && j!=0) min = px[--j];
      if(j != 0) for(int i = j; i--; ) {
        if(min > px[i]) min = px[i];
      }
    } else {
      min = px[0];
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i])) {
          min = px[i];
          break;
        } else {
          if(min > px[i]) min = px[i];
        }
      }
    }
    pout[0] = min;
  } else {
    if(narm) {
      for(int i = ng; i--; ) pout[i] = NA_REAL; // Other way ?
      --pout;
      for(int i = l; i--; ) if(pout[pg[i]] > px[i] || ISNAN(pout[pg[i]])) pout[pg[i]] = px[i];  // fastest
    } else {
      for(int i = ng; i--; ) pout[i] = POS_INF;
      --pout;
      for(int i = l; i--; ) if(pout[pg[i]] > px[i] || ISNAN(px[i])) pout[pg[i]] = px[i];  // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    }
  }
}

void fmin_int_impl(int *pout, int *px, int ng, int *pg, int narm, int l) {
  if(ng == 0) {
    int min;
    if(narm) {
      int j = l-1;
      min = px[j];
      while(min == NA_INTEGER && j!=0) min = px[--j];
      if(j != 0) for(int i = j; i--; ) {
        if(min > px[i] && px[i] != NA_INTEGER) min = px[i];
      }
    } else {
      min = px[0];
      for(int i = 0; i != l; ++i) {
        if(px[i] == NA_INTEGER) {
          min = NA_INTEGER;
          break;
        } else {
          if(min > px[i]) min = px[i];
        }
      }
    }
    pout[0] = min;
  } else {
    if(narm) {
      for(int i = ng; i--; ) pout[i] = NA_INTEGER;
      --pout;
      for(int i = l; i--; ) if(px[i] != NA_INTEGER && (pout[pg[i]] > px[i] || pout[pg[i]] == NA_INTEGER)) pout[pg[i]] = px[i];  // fastest??
    } else {
      for(int i = ng; i--; ) pout[i] = INT_MAX;
      --pout;
      for(int i = l; i--; ) if(pout[pg[i]] > px[i]) pout[pg[i]] = px[i];
    }
  }
}

void fmax_double_impl(double *pout, double *px, int ng, int *pg, int narm, int l) {
  if(ng == 0) {
    double max;
    if(narm) {
      int j = l-1;
      max = px[j];
      while(ISNAN(max) && j!=0) max = px[--j];
      if(j != 0) for(int i = j; i--; ) {
        if(max < px[i]) max = px[i];
      }
    } else {
      max = px[0];
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i])) {
          max = px[i];
          break;
        } else {
          if(max < px[i]) max = px[i];
        }
      }
    }
    pout[0] = max;
  } else {
    if(narm) {
      for(int i = ng; i--; ) pout[i] = NA_REAL; // Other way ?
      --pout;
      for(int i = l; i--; ) if(pout[pg[i]] < px[i] || ISNAN(pout[pg[i]])) pout[pg[i]] = px[i];  // fastest
    } else {
      for(int i = ng; i--; ) pout[i] = NEG_INF;
      --pout;
      for(int i = l; i--; ) if(pout[pg[i]] < px[i] || ISNAN(px[i])) pout[pg[i]] = px[i];  // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    }
  }
}

void fmax_int_impl(int *pout, int *px, int ng, int *pg, int narm, int l) {
  if(ng == 0) {
    int max;
    if(narm) {
      max = NA_INTEGER; // same as INT_MIN
      for(int i = l; i--; ) if(max < px[i]) max = px[i];
    } else {
      max = px[0];
      for(int i = 0; i != l; ++i) {
        if(px[i] == NA_INTEGER) {
          max = NA_INTEGER;
          break;
        } else {
          if(max < px[i]) max = px[i];
        }
      }
    }
    pout[0] = max;
  } else {
    if(narm) {
      for(int i = ng; i--; ) pout[i] = NA_INTEGER;
      --pout;
      for(int i = l; i--; ) if(pout[pg[i]] < px[i]) pout[pg[i]] = px[i];  // fastest??
    } else {
      for(int i = ng; i--; ) pout[i] = INT_MIN + 1; // best ??
      --pout;
      for(int i = l; i--; ) if(px[i] == NA_INTEGER || (pout[pg[i]] != NA_INTEGER && pout[pg[i]] < px[i])) pout[pg[i]] = px[i];
    }
  }
}


SEXP fminC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm) {
  int l = length(x), tx = TYPEOF(x), ng = asInteger(Rng), narm = asLogical(Rnarm);
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match length(x)");
  if(tx == LGLSXP) tx = INTSXP;
  // ALTREP methods for compact sequences: not safe yet and not part of the API.
  // if(ALTREP(x) && ng == 0) {
  // if(tx == INTSXP) return ALTINTEGER_MIN(x, (Rboolean)narm);
  // if(tx == REALSXP) return ALTREAL_MIN(x, (Rboolean)narm);
  // error("ALTREP object must be integer or real typed");
  // }
  SEXP out = PROTECT(allocVector(tx, ng == 0 ? 1 : ng));
  switch(tx) {
  case REALSXP: fmin_double_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l);
    break;
  case INTSXP: fmin_int_impl(INTEGER(out), INTEGER(x), ng, INTEGER(g), narm, l);
    break;
  default: error("Unsupported SEXP type");
  }
  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")))
    copyMostAttrib(x, out);
  UNPROTECT(1);
  return out;
}

SEXP fminmC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1], *pg = INTEGER(g),
    ng = asInteger(Rng), ng1 = ng == 0 ? 1 : ng, narm = asLogical(Rnarm);
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match nrow(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(tx, ng == 0 ? col : col * ng));
  switch(tx) {
  case REALSXP: {
    double *px = REAL(x), *pout = REAL(out);
    for(int j = 0; j != col; ++j) fmin_double_impl(pout + j*ng1, px + j*l, ng, pg, narm, l);
    break;
  }
  case INTSXP: {
    int *px = INTEGER(x), *pout = INTEGER(out);
    for(int j = 0; j != col; ++j) fmin_int_impl(pout + j*ng1, px + j*l, ng, pg, narm, l);
    break;
  }
  default: error("Unsupported SEXP type");
  }
  matCopyAttr(out, x, Rdrop, ng);
  UNPROTECT(1);
  return out;
}

SEXP fminlC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop) {
  int l = length(x), ng = asInteger(Rng);
  if(l < 1) return x; // needed ??
  if(ng == 0 && asLogical(Rdrop)) {
    SEXP out = PROTECT(allocVector(REALSXP, l));
    const SEXP *px = SEXPPTR_RO(x);
    double *pout = REAL(out);
    for(int j = 0; j != l; ++j) pout[j] = asReal(fminC(px[j], Rng, g, Rnarm));
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return out;
  }
  SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out);
  const SEXP *px = SEXPPTR_RO(x);
  for(int j = 0; j != l; ++j) pout[j] = fminC(px[j], Rng, g, Rnarm);
  // if(ng == 0) for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]);
  DFcopyAttr(out, x, ng);
  UNPROTECT(1);
  return out;
}


SEXP fmaxC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm) {
  int l = length(x), tx = TYPEOF(x), ng = asInteger(Rng), narm = asLogical(Rnarm);
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match length(x)");
  if(tx == LGLSXP) tx = INTSXP;
  // ALTREP methods for compact sequences: not safe yet and not part of the API.
  // if(ALTREP(x) && ng == 0) {
  // if(tx == INTSXP) return ALTINTEGER_MAX(x, (Rboolean)narm);
  // if(tx == REALSXP) return ALTREAL_MAX(x, (Rboolean)narm);
  // error("ALTREP object must be integer or real typed");
  // }
  SEXP out = PROTECT(allocVector(tx, ng == 0 ? 1 : ng));
  switch(tx) {
  case REALSXP: fmax_double_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l);
    break;
  case INTSXP: fmax_int_impl(INTEGER(out), INTEGER(x), ng, INTEGER(g), narm, l);
    break;
  default: error("Unsupported SEXP type");
  }
  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")))
    copyMostAttrib(x, out);
  UNPROTECT(1);
  return out;
}

SEXP fmaxmC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1], *pg = INTEGER(g),
    ng = asInteger(Rng), ng1 = ng == 0 ? 1 : ng, narm = asLogical(Rnarm);
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match nrow(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(tx, ng == 0 ? col : col * ng));
  switch(tx) {
  case REALSXP: {
    double *px = REAL(x), *pout = REAL(out);
    for(int j = 0; j != col; ++j) fmax_double_impl(pout + j*ng1, px + j*l, ng, pg, narm, l);
    break;
  }
  case INTSXP: {
    int *px = INTEGER(x), *pout = INTEGER(out);
    for(int j = 0; j != col; ++j) fmax_int_impl(pout + j*ng1, px + j*l, ng, pg, narm, l);
    break;
  }
  default: error("Unsupported SEXP type");
  }
  matCopyAttr(out, x, Rdrop, ng);
  UNPROTECT(1);
  return out;
}

SEXP fmaxlC(SEXP x, SEXP Rng, SEXP g, SEXP Rnarm, SEXP Rdrop) {
  int l = length(x), ng = asInteger(Rng);
  if(l < 1) return x; // needed ??
  if(ng == 0 && asLogical(Rdrop)) {
    SEXP out = PROTECT(allocVector(REALSXP, l));
    const SEXP *px = SEXPPTR_RO(x);
    double *pout = REAL(out);
    for(int j = 0; j != l; ++j) pout[j] = asReal(fmaxC(px[j], Rng, g, Rnarm));
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return out;
  }
  SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out);
  const SEXP *px = SEXPPTR_RO(x);
  for(int j = 0; j != l; ++j) pout[j] = fmaxC(px[j], Rng, g, Rnarm);
  // if(ng == 0) for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]);
  DFcopyAttr(out, x, ng);
  UNPROTECT(1);
  return out;
}
