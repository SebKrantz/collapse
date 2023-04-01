#include "collapse_c.h"
// #include <R_ext/Altrep.h>

void fprod_double_impl(double *pout, double *px, int ng, int *pg, int narm, int l) {
  if(ng == 0) {
    long double prod;
    if(narm) {
      int j = l-1;
      while(ISNAN(px[j]) && j!=0) --j;
      prod = (long double)px[j];
      if(j != 0) for(int i = j; i--; ) {
        if(NISNAN(px[i])) prod *= px[i]; // Fastest ?
      }
    } else {
      prod = 1.0;
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
      for(int i = ng; i--; ) pout[i] = 1.0;
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
      prod = 1.0;
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
      for(int i = ng; i--; ) pout[i] = 1.0;
      --pout;
      for(int i = l; i--; ) pout[pg[i]] *= px[i] * pw[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    }
  }
}

// using long long internally is substantially faster than using doubles !!
double fprod_int_impl(int *px, int narm, int l) {
  double prod;
  if(narm) {
    int j = l-1;
    while(px[j] == NA_INTEGER && j!=0) --j;
    prod = px[j];
    if(j == 0 && (l > 1 || px[j] == NA_INTEGER)) return NA_REAL;
    for(int i = j; i--; ) if(px[i] != NA_INTEGER) prod *= px[i];
  } else {
    prod = 1;
    for(int i = 0; i != l; ++i) {
      if(px[i] == NA_INTEGER) return NA_REAL;
      prod *= px[i];
    }
  }
  return prod;
}

void fprod_int_g_impl(double *pout, int *px, int ng, int *pg, int narm, int l) {
  if(narm) {
    for(int i = ng; i--; ) pout[i] = NA_REAL;
    for(int i = l, gi; i--; ) {
      if(px[i] != NA_INTEGER) {
        gi = pg[i]-1;
        if(ISNAN(pout[gi])) pout[gi] = (double)px[i];
        else pout[gi] *= px[i];
      }
    }
  } else {
    for(int i = ng; i--; ) pout[i] = 1.0;
    --pout;
    for(int i = l; i--; ) pout[pg[i]] *= px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
  }
}


SEXP fprodC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm) {
  int l = length(x), tx = TYPEOF(x), ng = asInteger(Rng),
    narm = asLogical(Rnarm), nprotect = 1;
  if (l < 1) return tx == REALSXP ? x : ScalarReal(asReal(x)); // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match length(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(REALSXP, ng == 0 ? 1 : ng));
  if(isNull(w)) {
    switch(tx) {
      case REALSXP: fprod_double_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l);
        break;
      case INTSXP: {
        if(ng > 0) fprod_int_g_impl(REAL(out), INTEGER(x), ng, INTEGER(g), narm, l);
        else REAL(out)[0] = fprod_int_impl(INTEGER(x), narm, l);
        break;
      }
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
  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")))
    copyMostAttrib(x, out); // For example "Units" objects...
  UNPROTECT(nprotect);
  return out;
}

SEXP fprodmC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1], *pg = INTEGER(g),
      ng = asInteger(Rng), ng1 = ng == 0 ? 1 : ng,
      narm = asLogical(Rnarm), nprotect = 1;
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match nrow(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(REALSXP, ng == 0 ? col : col * ng));
  double *pout = REAL(out);
  if(isNull(w)) {
    switch(tx) {
      case REALSXP: {
        double *px = REAL(x);
        for(int j = 0; j != col; ++j) fprod_double_impl(pout + j*ng1, px + j*l, ng, pg, narm, l);
        break;
      }
      case INTSXP: {
        int *px = INTEGER(x);
        if(ng > 0) {
          for(int j = 0; j != col; ++j) fprod_int_g_impl(pout + j*ng1, px + j*l, ng, pg, narm, l);
        } else {
          for(int j = 0; j != col; ++j) pout[j] = fprod_int_impl(px + j*l, narm, l);
        }
        break;
      }
      default: error("Unsupported SEXP type");
    }
  } else {
    if(l != length(w)) error("length(w) must match nrow(x)");
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
    SEXP out = PROTECT(allocVector(REALSXP, l));
    const SEXP *px = SEXPPTR_RO(x);
    double *pout = REAL(out);
    for(int j = 0; j != l; ++j) pout[j] = REAL(fprodC(px[j], Rng, g, w, Rnarm))[0];
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return out;
  }
  SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out);
  const SEXP *px = SEXPPTR_RO(x);
  for(int j = 0; j != l; ++j) pout[j] = fprodC(px[j], Rng, g, w, Rnarm);
  // if(ng == 0) for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]);
  DFcopyAttr(out, x, ng);
  UNPROTECT(1);
  return out;
}
