#include "collapse_c.h"

void fprod_double_impl(double *pout, double *px, int ng, int *pg, int narm, int l) {
  if(ng == 0) {
    long double prod;
    if(narm) {
      int j = l-1;
      prod = px[j];
      while(ISNAN(prod) && j!=0) prod = px[--j];
      if(j != 0) for(int i = j; i--; ) {
        if(NISNAN(px[i])) prod *= px[i]; // Fastest ?
      }
    } else {
      prod = 1;
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i])) {
          prod = px[i];
          break;
        } else {
          prod *= px[i];
        }
      }
    }
    pout[0] = (double)prod;
  } else {
    if(narm) {
      for(int i = ng; i--; ) pout[i] = NA_REAL; // Other way ?
      --pout;
      for(int i = l; i--; ) {
        if(NISNAN(px[i])) { // faster way to code this ? -> Not Bad at all
          if(ISNAN(pout[pg[i]])) pout[pg[i]] = px[i];
          else pout[pg[i]] *= px[i];
        }
      }
    } else {
      for(int i = ng; i--; ) pout[i] = 1.0; // Other way ?
      --pout;
      for(int i = l; i--; ) pout[pg[i]] *= px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    }
  }
}

void fprod_weights_impl(double *pout, double *px, int ng, int *pg, double *pw, int narm, int l) {
  if(ng == 0) {
    long double prod;
    if(narm) {
      int j = l-1;
      while((ISNAN(px[j]) || ISNAN(pw[j])) && j!=0) --j;
      prod = px[j] * pw[j];
      if(j != 0) for(int i = j; i--; ) {
        if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
        prod *= px[i] * pw[i];
      }
    } else {
      prod = 1;
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i]) || ISNAN(pw[i])) {
          prod = px[i] + pw[i];
          break;
        } else {
          prod *= px[i] * pw[i];
        }
      }
    }
    pout[0] = (double)prod;
  } else {
    if(narm) {
      for(int i = ng; i--; ) pout[i] = NA_REAL; // Other way ?
      --pout;
      for(int i = l; i--; ) {
        if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
        if(ISNAN(pout[pg[i]])) pout[pg[i]] = px[i] * pw[i];
        else pout[pg[i]] *= px[i] * pw[i];
      }
    } else {
      for(int i = ng; i--; ) pout[i] = 1.0; // Other way ?
      --pout;
      for(int i = l; i--; ) pout[pg[i]] *= px[i] * pw[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    }
  }
}

void fprod_int_impl(double *pout, int *px, int ng, int *pg, int narm, int l) {
  if(ng == 0) {
    double prod;
    if(narm) {
      int j = l-1;
      while(ISNAN(px[j]) && j!=0) --j;
      if(j != 0) {
        prod = (double)px[j];
        for(int i = j; i--; ) if(px[i] != NA_INTEGER) prod *= (double)px[i];
      } else prod = NA_REAL;
    } else {
      prod = 1.0;
      for(int i = 0; i != l; ++i) {
        if(px[i] == NA_INTEGER) {
          prod = NA_REAL;
          break;
        } else {
          prod *= (double)px[i];
        }
      }
    }
    pout[0] = prod;
  } else {
    if(narm) {
      for(int i = ng; i--; ) pout[i] = NA_REAL;
      --pout;
      for(int i = l; i--; ) {
        if(px[i] != NA_INTEGER) {
          if(ISNAN(pout[pg[i]])) pout[pg[i]] = (double)px[i];
          else pout[pg[i]] *= (double)px[i];
        }
      }
    } else {
      for(int i = ng; i--; ) pout[i] = 1.0;
      --pout;
      for(int i = l; i--; ) {
        if(px[i] == NA_INTEGER) pout[pg[i]] = NA_REAL;
        else pout[pg[i]] *= (double)px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
      }
    }
  }
}


SEXP fprodC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm) {
  int l = length(x), tx = TYPEOF(x), ng = asInteger(Rng),
    narm = asInteger(Rnarm), nprotect = 1;
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match length(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(REALSXP, ng == 0 ? 1 : ng));
  if(isNull(w)) {
    switch(tx) {
    case REALSXP: fprod_double_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l);
      break;
    case INTSXP: fprod_int_impl(REAL(out), INTEGER(x), ng, INTEGER(g), narm, l);
      break;
    default: error("Unsupported SEXP type");
    }
  } else {
    if(l != length(w)) error("length(w) must match length(x)");
    int tw = TYPEOF(w);
    SEXP xr, wr;
    double *px, *pw;
    if(tw != REALSXP) {
      if(tw != INTSXP && tw != LGLSXP) error("weigths must be double or integer");
      wr = PROTECT(coerceVector(w, REALSXP));
      pw = REAL(wr);
      ++nprotect;
    } else pw = REAL(w);
    if(tx != REALSXP) {
      if(tx != INTSXP) error("x must be double or integer");
      xr = PROTECT(coerceVector(x, REALSXP));
      px = REAL(xr);
      ++nprotect;
    } else px = REAL(x);
    fprod_weights_impl(REAL(out), px, ng, INTEGER(g), pw, narm, l);
  }
  if(ng && !isObject(x)) copyMostAttrib(x, out);
  UNPROTECT(nprotect);
  return out;
}

SEXP fprodmC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1], *pg = INTEGER(g),
    ng = asInteger(Rng), ng1 = ng == 0 ? 1 : ng,
    narm = asInteger(Rnarm), nprotect = 1;
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match nrow(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(REALSXP, ng == 0 ? col : col * ng));
  if(isNull(w)) {
    switch(tx) {
    case REALSXP: {
      double *px = REAL(x), *pout = REAL(out);
      for(int j = 0; j != col; ++j) fprod_double_impl(pout + j*ng1, px + j*l, ng, pg, narm, l);
      break;
    }
    case INTSXP: {
      int *px = INTEGER(x);
      double *pout = REAL(out);
      for(int j = 0; j != col; ++j) fprod_int_impl(pout + j*ng1, px + j*l, ng, pg, narm, l);
      break;
    }
    default: error("Unsupported SEXP type");
    }
  } else {
    if(l != length(w)) error("length(w) must match nrow(x)");
    int tw = TYPEOF(w);
    SEXP xr, wr;
    double *px, *pw, *pout = REAL(out);
    if(tw != REALSXP) {
      if(tw != INTSXP && tw != LGLSXP) error("weigths must be double or integer");
      wr = PROTECT(coerceVector(w, REALSXP));
      pw = REAL(wr);
      ++nprotect;
    } else pw = REAL(w);
    if(tx != REALSXP) {
      if(tx != INTSXP) error("x must be double or integer");
      xr = PROTECT(coerceVector(x, REALSXP));
      px = REAL(xr);
      ++nprotect;
    } else px = REAL(x);
    for(int j = 0; j != col; ++j) fprod_weights_impl(pout + j*ng1, px + j*l, ng, pg, pw, narm, l);
  }
  matCopyAttr(out, x, Rdrop, ng);
  UNPROTECT(nprotect);
  return out;
}

SEXP fprodlC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop) {
  int l = length(x), ng = asInteger(Rng);
  if(l < 1) return x; // needed ??
  if(ng == 0 && asLogical(Rdrop)) {
    SEXP out = PROTECT(allocVector(REALSXP, l)), *px = SEXPPTR(x);
    double *pout = REAL(out);
    for(int j = 0; j != l; ++j) pout[j] = asReal(fprodC(px[j], Rng, g, w, Rnarm));
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return out;
  }
  SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out), *px = SEXPPTR(x);
  for(int j = 0; j != l; ++j) pout[j] = fprodC(px[j], Rng, g, w, Rnarm);
  if(ng == 0) for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]);
  DFcopyAttr(out, x, ng);
  UNPROTECT(1);
  return out;
}
