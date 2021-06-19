#include "collapse_c.h"

void fsum_double_impl(double *pout, double *px, int ng, int *pg, int narm, int l) {
  if(ng == 0) {
    double sum;
    if(narm) {
      int j = l-1;
      sum = px[j];
      while(ISNAN(sum) && j!=0) sum = px[--j];
      if(j != 0) for(int i = j; i--; ) {
        if(NISNAN(px[i])) sum += px[i]; // Fastest ?
      }
    } else {
      sum = 0;
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i])) {
          sum = px[i];
          break;
        } else {
          sum += px[i];
        }
      }
    }
    pout[0] = sum;
  } else {
    if(narm) {
      for(int i = ng; i--; ) pout[i] = NA_REAL; // Other way ?
      --pout;
      for(int i = l; i--; ) {
        if(NISNAN(px[i])) { // faster way to code this ? -> Not Bad at all
          if(ISNAN(pout[pg[i]])) pout[pg[i]] = px[i];
          else pout[pg[i]] += px[i];
        }
      }
    } else {
      memset(pout, 0.0, sizeof(double) * ng);
      --pout;
      for(int i = l; i--; ) pout[pg[i]] += px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    }
  }
}

void fsum_weights_impl(double *pout, double *px, int ng, int *pg, double *pw, int narm, int l) {
  if(ng == 0) {
    double sum;
    if(narm) {
      int j = l-1;
      while((ISNAN(px[j]) || ISNAN(pw[j])) && j!=0) --j;
      sum = px[j] * pw[j];
      if(j != 0) for(int i = j; i--; ) {
        if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
        sum += px[i] * pw[i];
      }
    } else {
      sum = 0;
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i]) || ISNAN(pw[i])) {
          sum = px[i] + pw[i];
          break;
        } else {
          sum += px[i] * pw[i];
        }
      }
    }
    pout[0] = sum;
  } else {
    if(narm) {
      for(int i = ng; i--; ) pout[i] = NA_REAL; // Other way ?
      --pout;
      for(int i = l; i--; ) {
        if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
        if(ISNAN(pout[pg[i]])) pout[pg[i]] = px[i] * pw[i];
        else pout[pg[i]] += px[i] * pw[i];
      }
    } else {
      memset(pout, 0.0, sizeof(double) * ng);
      --pout;
      for(int i = l; i--; ) pout[pg[i]] += px[i] * pw[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    }
  }
}

// using long long internally is substantially faster than using doubles !!
double fsum_int_impl(int *px, int narm, int l) {
  long long sum;
  if(narm) {
    int j = l-1;
    while(px[j] == NA_INTEGER && j!=0) --j;
    sum = (long long)px[j];
    if(j == 0 && (l > 1 || px[j] == NA_INTEGER)) return NA_REAL;
    for(int i = j; i--; ) if(px[i] != NA_INTEGER) sum += (long long)px[i];
  } else {
    sum = 0;
    for(int i = 0; i != l; ++i) {
      if(px[i] == NA_INTEGER) return NA_REAL;
      sum += (long long)px[i];
    }
  }
  return (double)sum;
}

void fsum_int_g_impl(int *pout, int *px, int ng, int *pg, int narm, int l) {
  long long ckof;
  if(narm) {
    for(int i = ng; i--; ) pout[i] = NA_INTEGER;
    --pout;
    for(int i = l, lsi; i--; ) {
      if(px[i] != NA_INTEGER) {
        lsi = pout[pg[i]];
        if(lsi == NA_INTEGER) pout[pg[i]] = px[i];
        else {
          ckof = (long long)lsi + px[i];
          if(ckof > INT_MAX || ckof <= INT_MIN) error("Integer overflow in one or more groups. Integers in R are bounded between 2,147,483,647 and -2,147,483,647. The sum within each group should be in that range.");
          pout[pg[i]] = (int)ckof;
        }
      }
    }
  } else {
    memset(pout, 0, sizeof(int) * ng);
    --pout;
    for(int i = l, lsi; i--; ) {
      if(px[i] == NA_INTEGER) {
        pout[pg[i]] = NA_INTEGER;
        continue;
      }
      lsi = pout[pg[i]];
      if(lsi != NA_INTEGER) { // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
        ckof = (long long)lsi + px[i];
        if(ckof > INT_MAX || ckof <= INT_MIN) error("Integer overflow in one or more groups. Integers in R are bounded between 2,147,483,647 and -2,147,483,647. The sum within each group should be in that range.");
        pout[pg[i]] = (int)ckof;
      }
    }
  }
}


SEXP fsumC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm) {
  int l = length(x), tx = TYPEOF(x), ng = asInteger(Rng),
    narm = asLogical(Rnarm), nprotect = 1, nwl = isNull(w);
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match length(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out;
  if(!(ng == 0 && nwl && tx == INTSXP))
    out = PROTECT(allocVector(nwl ? tx : REALSXP, ng == 0 ? 1 : ng));
  if(nwl) {
    switch(tx) {
      case REALSXP: fsum_double_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l);
        break;
      case INTSXP: {
        if(ng > 0) fsum_int_g_impl(INTEGER(out), INTEGER(x), ng, INTEGER(g), narm, l);
        else {
          double sum = fsum_int_impl(INTEGER(x), narm, l);
          if(sum > INT_MAX || sum <= INT_MIN) return ScalarReal(sum); // INT_MIN is NA_INTEGER
          return ScalarInteger((int)sum);
        }
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
    fsum_weights_impl(REAL(out), px, ng, INTEGER(g), pw, narm, l);
  }
  if(ng && !isObject(x)) copyMostAttrib(x, out);
  UNPROTECT(nprotect);
  return out;
}

SEXP fsummC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1], *pg = INTEGER(g),
      ng = asInteger(Rng), ng1 = ng == 0 ? 1 : ng,
      narm = asLogical(Rnarm), nprotect = 1, nwl = isNull(w);
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match nrow(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector((nwl && ng > 0) ? tx : REALSXP, ng == 0 ? col : col * ng));
  if(nwl) {
    switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        for(int j = 0; j != col; ++j) fsum_double_impl(pout + j*ng1, px + j*l, ng, pg, narm, l);
        break;
      }
      case INTSXP: {
        int *px = INTEGER(x);
        if(ng > 0) {
          int *pout = INTEGER(out);
          for(int j = 0; j != col; ++j) fsum_int_g_impl(pout + j*ng1, px + j*l, ng, pg, narm, l);
        } else {
          double *pout = REAL(out);
          int anyoutl = 0;
          for(int j = 0; j != col; ++j) {
            double sumj = fsum_int_impl(px + j*l, narm, l);
            if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
            pout[j] = sumj;
          }
          if(anyoutl == 0) {
            SEXP iout = PROTECT(coerceVector(out, INTSXP));
            matCopyAttr(iout, x, Rdrop, ng);
            UNPROTECT(2);
            return iout;
          }
        }
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
    for(int j = 0; j != col; ++j) fsum_weights_impl(pout + j*ng1, px + j*l, ng, pg, pw, narm, l);
  }
  matCopyAttr(out, x, Rdrop, ng);
  UNPROTECT(nprotect);
  return out;
}

SEXP fsumlC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop) {
  int l = length(x), ng = asInteger(Rng);
  if(l < 1) return x; // needed ??
  if(ng == 0 && asLogical(Rdrop)) {
    SEXP out = PROTECT(allocVector(REALSXP, l)), *px = SEXPPTR(x);
    double *pout = REAL(out);
    for(int j = 0; j != l; ++j) pout[j] = asReal(fsumC(px[j], Rng, g, w, Rnarm));
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return out;
  }
  SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out), *px = SEXPPTR(x);
  for(int j = 0; j != l; ++j) pout[j] = fsumC(px[j], Rng, g, w, Rnarm);
  if(ng == 0) for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]);
  DFcopyAttr(out, x, ng);
  UNPROTECT(1);
  return out;
}
