#include "collapse_c.h"
// #include <R_ext/Altrep.h>

void fmean_double_impl(double *pout, double *px, int ng, int *pg, int *pgs, int narm, int l) {
  if(ng == 0) {
    double mean = 0;
    if(narm) {
      int n = 0;
      for(int i = 0; i != l; ++i) {
        if(NISNAN(px[i])) {
          mean += px[i]; // Fastest ?
          ++n;
        }
      }
      pout[0] = n == 0 ? NA_REAL : mean / n;
    } else {
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i])) {
          mean = px[i];
          break;
        }
        mean += px[i];
      }
      pout[0] = mean / l;
    }
  } else {
    if(narm) {
      int *n = (int*)Calloc(ng, int);
      for(int i = ng; i--; ) pout[i] = NA_REAL; // Other way ?
      for(int i = l, gi; i--; ) {
        if(NISNAN(px[i])) { // faster way to code this ? -> Not Bad at all
          gi = pg[i]-1;
          if(ISNAN(pout[gi])) {
            pout[gi] = px[i];
            n[gi] = 1;
          } else {
            pout[gi] += px[i];
            ++n[gi];
          }
        }
      }
      for(int i = ng; i--; ) pout[i] /= n[i]; // could use R_alloc above, but what about this loop?
      Free(n);
    } else {
      memset(pout, 0.0, sizeof(double) * ng);
      --pout;
      for(int i = l; i--; ) pout[pg[i]] += px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
      ++pout;
      for(int i = ng; i--; ) pout[i] /= pgs[i];
    }
  }
}

void fmean_weights_impl(double *pout, double *px, int ng, int *pg, double *pw, int narm, int l) {
  if(ng == 0) {
    double mean, sumw;
    if(narm) {
      int j = l-1;
      while((ISNAN(px[j]) || ISNAN(pw[j])) && j!=0) --j;
      sumw = pw[j];
      mean = px[j] * sumw;
      if(j != 0) for(int i = j; i--; ) {
        if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
        mean += px[i] * pw[i];
        sumw += pw[i];
      }
    } else {
      mean = 0, sumw = 0;
      for(int i = 0; i != l; ++i) {
        if(ISNAN(px[i]) || ISNAN(pw[i])) {
          mean = px[i] + pw[i];
          break;
        } else {
          mean += px[i] * pw[i];
          sumw += pw[i];
        }
      }
    }
    pout[0] = mean / sumw;
  } else {
    double *sumw = (double*)Calloc(ng, double);
    if(narm) {
      for(int i = ng; i--; ) pout[i] = NA_REAL; // Other way ?
      for(int i = l, gi; i--; ) {
        if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
        gi = pg[i]-1;
        if(ISNAN(pout[gi])) {
          pout[gi] = px[i] * pw[i];
          sumw[gi] = pw[i];
        } else {
          pout[gi] += px[i] * pw[i];
          sumw[gi] += pw[i];
        }
      }
    } else {
      memset(pout, 0.0, sizeof(double) * ng);
      for(int i = l, gi; i--; ) {
        gi = pg[i]-1;
        pout[gi] += px[i] * pw[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
        sumw[gi] += pw[i];
      }
    }
    for(int i = ng; i--; ) pout[i] /= sumw[i];
    Free(sumw);
  }
}

// using long long internally is substantially faster than using doubles !!
double fmean_int_impl(int *px, int narm, int l) {
  long long mean;
  double dmean;
  if(narm) {
    int j = l-1, k = 1;
    while(px[j] == NA_INTEGER && j!=0) --j;
    mean = px[j];
    if(j == 0 && (l > 1 || px[j] == NA_INTEGER)) return NA_REAL;
    for(int i = j; i--; ) {
      if(px[i] != NA_INTEGER) {
        mean += px[i];
        ++k;
      }
    }
    dmean = (double)mean / k;
  } else {
    mean = 0;
    for(int i = 0; i != l; ++i) {
      if(px[i] == NA_INTEGER) return NA_REAL;
      mean += px[i];
    }
    dmean = (double)mean / l;
  }
  return dmean;
}

void fmean_int_g_impl(double *pout, int *px, int ng, int *pg, int *pgs, int narm, int l) {
  if(narm) {
    int *n = (int*)Calloc(ng, int);
    for(int i = ng; i--; ) pout[i] = NA_REAL;
    for(int i = l, gi; i--; ) {
      if(px[i] != NA_INTEGER) {
        gi = pg[i]-1;
        if(ISNAN(pout[gi])) {
          pout[gi] = (double)px[i];
          n[gi] = 1;
        } else {
          pout[gi] += px[i];
          ++n[gi];
        }
      }
    }
    for(int i = ng; i--; ) pout[i] /= n[i];
    Free(n);
  } else {
    memset(pout, 0.0, sizeof(double) * ng);
    --pout;
    for(int i = l; i--; ) pout[pg[i]] += px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    ++pout;
    for(int i = ng; i--; ) pout[i] /= pgs[i];
  }
}


SEXP fmeanC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm) {
  int l = length(x), tx = TYPEOF(x), ng = asInteger(Rng),
    narm = asLogical(Rnarm), nprotect = 1, nwl = isNull(w), *pgs = &l;
  // ALTREP methods for compact sequences: not safe yet and not part of the API.
  // if(ALTREP(x) && ng == 0 && nwl) {
  // switch(tx) {
  // case INTSXP: return ALTINTEGER_SUM(x, (Rboolean)narm);
  // case LGLSXP: return ALTLOGICAL_SUM(x, (Rboolean)narm);
  // case REALSXP: return ALTREAL_SUM(x, (Rboolean)narm);
  // default: error("ALTREP object must be integer or real typed");
  // }
  // }
  if (l < 1) return tx == REALSXP ? x : ScalarReal(asReal(x)); // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match length(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(REALSXP, ng == 0 ? 1 : ng));
  if(nwl) {
    if(ng && !narm) {
      if(length(gs) == ng) pgs = INTEGER(gs);
      else {
        SEXP gs_ = PROTECT(allocVector(INTSXP, ng)); ++nprotect;
        pgs = INTEGER(gs_);
        memset(pgs, 0, sizeof(int) * ng);
        for(int i = 0, *pg = INTEGER(g); i != l; ++i) ++pgs[pg[i]-1];
      }
    }
    switch(tx) {
      case REALSXP: fmean_double_impl(REAL(out), REAL(x), ng, INTEGER(g), pgs, narm, l);
        break;
      case INTSXP: {
        if(ng > 0) fmean_int_g_impl(REAL(out), INTEGER(x), ng, INTEGER(g), pgs, narm, l);
        else REAL(out)[0] = fmean_int_impl(INTEGER(x), narm, l);
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
    fmean_weights_impl(REAL(out), px, ng, INTEGER(g), pw, narm, l);
  }
  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")))
     copyMostAttrib(x, out); // ATTRIB(x) != R_NilValue? // For example "Units" objects...
  UNPROTECT(nprotect);
  return out;
}

SEXP fmeanmC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rdrop) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1], *pg = INTEGER(g),
    ng = asInteger(Rng), ng1 = ng == 0 ? 1 : ng,
    narm = asLogical(Rnarm), nprotect = 1, *pgs = &l;
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match nrow(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(REALSXP, ng == 0 ? col : col * ng));
  double *pout = REAL(out);
  if(isNull(w)) {
    if(ng && !narm) {
      if(length(gs) == ng) pgs = INTEGER(gs);
      else {
        SEXP gs_ = PROTECT(allocVector(INTSXP, ng)); ++nprotect;
        pgs = INTEGER(gs_);
        memset(pgs, 0, sizeof(int) * ng);
        for(int i = 0, *pg = INTEGER(g); i != l; ++i) ++pgs[pg[i]-1];
      }
    }
    switch(tx) {
      case REALSXP: {
        double *px = REAL(x);
        for(int j = 0; j != col; ++j) fmean_double_impl(pout + j*ng1, px + j*l, ng, pg, pgs, narm, l);
        break;
      }
      case INTSXP: {
        int *px = INTEGER(x);
        if(ng > 0) {
          for(int j = 0; j != col; ++j) fmean_int_g_impl(pout + j*ng1, px + j*l, ng, pg, pgs, narm, l);
        } else {
          for(int j = 0; j != col; ++j) pout[j] = fmean_int_impl(px + j*l, narm, l);
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
    for(int j = 0; j != col; ++j) fmean_weights_impl(pout + j*ng1, px + j*l, ng, pg, pw, narm, l);
  }
  matCopyAttr(out, x, Rdrop, ng);
  UNPROTECT(nprotect);
  return out;
}

SEXP fmeanlC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rdrop) {
  int l = length(x), ng = asInteger(Rng);
  if(l < 1) return x; // needed ??
  if(ng == 0 && asLogical(Rdrop)) {
    SEXP out = PROTECT(allocVector(REALSXP, l)), *px = SEXPPTR(x);
    double *pout = REAL(out);
    for(int j = 0; j != l; ++j) pout[j] = REAL(fmeanC(px[j], Rng, g, gs, w, Rnarm))[0];
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return out;
  }
  SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out), *px = SEXPPTR(x);
  for(int j = 0; j != l; ++j) pout[j] = fmeanC(px[j], Rng, g, gs, w, Rnarm);
  // if(ng == 0) for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]);
  DFcopyAttr(out, x, ng);
  UNPROTECT(1);
  return out;
}
