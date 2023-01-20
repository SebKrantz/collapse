#ifdef _OPENMP
#include <omp.h>
#endif
#include "collapse_c.h"
// #include <R_ext/Altrep.h>

// Adapted from fsum.c

double fmean_double_impl(const double *restrict px, const int narm, const int l) {
  if(narm) {
    int j = l-1, n = 1;
    double mean = px[j];
    while(ISNAN(mean) && j!=0) mean = px[--j];
    if(j != 0) for(int i = j; i--; ) {
      if(ISNAN(px[i])) continue;
      mean += px[i];
      ++n;
    }
    return  mean / n;
  }
  double mean = 0;
  for(int i = 0; i != l; ++i) {
    if(ISNAN(px[i])) {
      mean = px[i];
      break;
    }
    mean += px[i];
  }
  return mean / l;
}

double fmean_double_omp_impl(const double *restrict px, const int narm, const int l, const int nthreads) {
  double mean = 0;
  if(narm) {
    int n = 0;
    #pragma omp parallel for num_threads(nthreads) reduction(+:mean,n)
    for(int i = 0; i < l; ++i) {
      if(ISNAN(px[i])) continue;
      mean += px[i];
      ++n;
    }
    return n == 0 ? NA_REAL : mean / n;
  }
  #pragma omp parallel for num_threads(nthreads) reduction(+:mean)
  for(int i = 0; i < l; ++i) mean += px[i];
  return mean / l;
}

void fmean_double_g_impl(double *restrict pout, const double *restrict px, const int ng, const int *restrict pg, const int *restrict pgs, const int narm, const int l) {
  memset(pout, 0, sizeof(double) * ng);
  if(narm) {
    int *restrict n = (int*)Calloc(ng, int);
    for(int i = 0, gi; i != l; ++i) {
      if(ISNAN(px[i])) continue;
      gi = pg[i]-1;
      pout[gi] += px[i];
      ++n[gi];
    }
    for(int i = ng; i--; ) {
      if(n[i] == 0) pout[i] = NA_REAL;
      else pout[i] /= n[i];
    }
    Free(n);
  } else {
    --pout;
    for(int i = l; i--; ) pout[pg[i]] += px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    ++pout;
    for(int i = ng; i--; ) pout[i] /= pgs[i];
  }
}

double fmean_weights_impl(const double *restrict px, const double *restrict pw, const int narm, const int l) {
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
      }
      mean += px[i] * pw[i];
      sumw += pw[i];
    }
  }
  return mean / sumw;
}

double fmean_weights_omp_impl(const double *restrict px, const double *restrict pw, const int narm, const int l, const int nthreads) {
  double mean = 0, sumw = 0;
  if(narm) {
    #pragma omp parallel for num_threads(nthreads) reduction(+:mean,sumw)
    for(int i = 0; i < l; ++i) {
      if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
      mean += px[i] * pw[i];
      sumw += pw[i];
    }
    if(mean == 0 && sumw == 0) sumw = NA_REAL;
  } else {
    #pragma omp parallel for num_threads(nthreads) reduction(+:mean,sumw)
    for(int i = 0; i < l; ++i) {
      mean += px[i] * pw[i];
      sumw += pw[i];
    }
  }
  return mean / sumw;
}

void fmean_weights_g_impl(double *restrict pout, const double *restrict px, const int ng, const int *restrict pg, const double *restrict pw, const int narm, const int l) {
  double *restrict sumw = (double*)Calloc(ng, double);
  memset(pout, 0, sizeof(double) * ng);
  if(narm) {
    for(int i = 0, gi; i != l; ++i) {
      if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
      gi = pg[i]-1;
      pout[gi] += px[i] * pw[i];
      sumw[gi] += pw[i];
    }
    for(int i = ng; i--; ) {
      if(sumw[i] == 0) pout[i] = NA_REAL;
      else pout[i] /= sumw[i];
    }
  } else {
    for(int i = 0, gi; i != l; ++i) {
      gi = pg[i]-1;
      pout[gi] += px[i] * pw[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
      sumw[gi] += pw[i];
    }
    for(int i = ng; i--; ) pout[i] /= sumw[i];
  }
  Free(sumw);
}

double fmean_int_impl(const int *restrict px, const int narm, const int l) {
  long long mean;
  double dmean;
  if(narm) {
    int j = l-1, k = 1;
    while(px[j] == NA_INTEGER && j!=0) --j;
    mean = px[j];
    if(j == 0 && (l > 1 || px[j] == NA_INTEGER)) return NA_REAL;
    for(int i = j; i--; ) {
      if(px[i] == NA_INTEGER) continue;
      mean += px[i];
      ++k;
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

double fmean_int_omp_impl(const int *restrict px, const int narm, const int l, const int nthreads) {
  long long mean = 0;
  double dmean;
  if(narm) {
    int n = 0;
    #pragma omp parallel for num_threads(nthreads) reduction(+:mean,n)
    for(int i = 0; i < l; ++i) {
      if(px[i] == NA_INTEGER) continue;
      mean += px[i];
      ++n;
    }
    dmean = n == 0 ? NA_REAL : (double)mean / n;
  } else {
    #pragma omp parallel for num_threads(nthreads) reduction(+:mean)
    for(int i = 0; i < l; ++i) mean += px[i];
    dmean = (double)mean / l;
  }
  return dmean;
}

void fmean_int_g_impl(double *restrict pout, const int *restrict px, const int ng, const int *restrict pg, const int *restrict pgs, const int narm, const int l) {
  memset(pout, 0, sizeof(double) * ng);
  if(narm) {
    int *restrict n = (int*)Calloc(ng, int);
    for(int i = 0, gi; i != l; ++i) {
      if(px[i] == NA_INTEGER) continue;
      gi = pg[i]-1;
      pout[gi] += px[i];
      ++n[gi];
    }
    for(int i = ng; i--; ) {
      if(n[i] == 0) pout[i] = NA_REAL;
      else pout[i] /= n[i];
    }
    Free(n);
  } else {
    --pout;
    for(int i = l; i--; ) pout[pg[i]] += px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    ++pout;
    for(int i = ng; i--; ) pout[i] /= pgs[i];
  }
}


SEXP fmeanC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rnthreads) {
  const int l = length(x), ng = asInteger(Rng), narm = asLogical(Rnarm), nwl = isNull(w);
  int tx = TYPEOF(x), nthreads = asInteger(Rnthreads), nprotect = 1, *restrict pgs = &nprotect;
  // ALTREP methods for compact sequences: not safe yet and not part of the API.
  // if(ALTREP(x) && ng == 0 && nwl) {
  // switch(tx) {
  // case INTSXP: return ALTINTEGER_SUM(x, (Rboolean)narm);
  // case LGLSXP: return ALTLOGICAL_SUM(x, (Rboolean)narm);
  // case REALSXP: return ALTREAL_SUM(x, (Rboolean)narm);
  // default: error("ALTREP object must be integer or real typed");
  // }
  // }
  if(l < 1) return tx == REALSXP ? x : ScalarReal(asReal(x)); // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match length(x)");
  if(nthreads > max_threads) nthreads = max_threads;
  if(l < 100000) nthreads = 1; // No improvements from multithreading on small data.
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(REALSXP, ng == 0 ? 1 : ng));
  if(nwl) {
    if(ng && !narm) {
      if(length(gs) == ng) pgs = INTEGER(gs);
      else { // TODO: this is probably slower than narm, which requires only one loop...
        SEXP gs_ = PROTECT(allocVector(INTSXP, ng)); ++nprotect;
        pgs = INTEGER(gs_);
        memset(pgs, 0, sizeof(int) * ng);
        for(int i = 0, *restrict pg = INTEGER(g); i != l; ++i) ++pgs[pg[i]-1];
      }
    }
    switch(tx) {
      case REALSXP: {
        if(ng > 0) fmean_double_g_impl(REAL(out), REAL(x), ng, INTEGER(g), pgs, narm, l);
        else REAL(out)[0] = (nthreads <= 1) ? fmean_double_impl(REAL(x), narm, l) : fmean_double_omp_impl(REAL(x), narm, l, nthreads);
        break;
      }
      case INTSXP: {
        if(ng > 0) fmean_int_g_impl(REAL(out), INTEGER(x), ng, INTEGER(g), pgs, narm, l);
        else REAL(out)[0] = nthreads <= 1 ? fmean_int_impl(INTEGER(x), narm, l) : fmean_int_omp_impl(INTEGER(x), narm, l, nthreads);
        break;
      }
      default: error("Unsupported SEXP type: '%s'", type2char(tx));
    }
  } else {
    if(l != length(w)) error("length(w) must match length(x)");
    if(TYPEOF(w) != REALSXP) {
      if(TYPEOF(w) != INTSXP && TYPEOF(w) != LGLSXP) error("weigths must be double or integer");
      w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
    }
    if(tx != REALSXP) {
      if(tx != INTSXP) error("Unsupported SEXP type: '%s'", type2char(tx));
      x = PROTECT(coerceVector(x, REALSXP)); ++nprotect;
    }
    double *restrict px = REAL(x), *restrict pw = REAL(w);
    if(ng == 0) {
      REAL(out)[0] = (nthreads <= 1) ? fmean_weights_impl(px, pw, narm, l) :
              fmean_weights_omp_impl(px, pw, narm, l, nthreads);
    } else fmean_weights_g_impl(REAL(out), px, ng, INTEGER(g), pw, narm, l);
  }
  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")))
     copyMostAttrib(x, out); // For example "Units" objects...
  UNPROTECT(nprotect);
  return out;
}

SEXP fmeanmC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  const int l = INTEGER(dim)[0], col = INTEGER(dim)[1], *restrict pg = INTEGER(g), ng = asInteger(Rng), narm = asLogical(Rnarm);
  int tx = TYPEOF(x), nthreads = asInteger(Rnthreads), nprotect = 1, *restrict pgs = &nprotect;
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match nrow(x)");
  if(l*col < 100000) nthreads = 1; // No gains from multithreading on small data
  if(nthreads > max_threads) nthreads = max_threads;
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector(REALSXP, ng == 0 ? col : col * ng));
  double *restrict pout = REAL(out);
  if(isNull(w)) {
    if(ng && !narm) {
      if(length(gs) == ng) pgs = INTEGER(gs);
      else {
        SEXP gs_ = PROTECT(allocVector(INTSXP, ng)); ++nprotect;
        pgs = INTEGER(gs_);
        memset(pgs, 0, sizeof(int) * ng);
        for(int i = 0, *restrict pg = INTEGER(g); i != l; ++i) ++pgs[pg[i]-1];
      }
    }
    switch(tx) {
      case REALSXP: {
        const double *px = REAL(x);
        if(ng == 0) {
          if(nthreads <= 1) {
            for(int j = 0; j != col; ++j) pout[j] = fmean_double_impl(px + j*l, narm, l);
          } else if(col >= nthreads) {
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) pout[j] = fmean_double_impl(px + j*l, narm, l);
          } else {
            for(int j = 0; j != col; ++j) pout[j] = fmean_double_omp_impl(px + j*l, narm, l, nthreads);
          }
        } else {
          if(nthreads <= 1 || col == 1) {
            for(int j = 0; j != col; ++j) fmean_double_g_impl(pout + j*ng, px + j*l, ng, pg, pgs, narm, l);
          } else {
            if(nthreads > col) nthreads = col;
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) fmean_double_g_impl(pout + j*ng, px + j*l, ng, pg, pgs, narm, l);
          }
        }
        break;
      }
      case INTSXP: {
        const int *px = INTEGER(x);
        if(ng > 0) {
          if(nthreads <= 1 || col == 1) {
            for(int j = 0; j != col; ++j) fmean_int_g_impl(pout + j*ng, px + j*l, ng, pg, pgs, narm, l);
          } else {
            if(nthreads > col) nthreads = col;
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) fmean_int_g_impl(pout + j*ng, px + j*l, ng, pg, pgs, narm, l);
          }
        } else {
          if(nthreads <= 1) {
            for(int j = 0; j != col; ++j) pout[j] = fmean_int_impl(px + j*l, narm, l);
          } else if(col >= nthreads) {
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) pout[j] = fmean_int_impl(px + j*l, narm, l);
          } else {
            for(int j = 0; j != col; ++j) pout[j] = fmean_int_omp_impl(px + j*l, narm, l, nthreads);
          }
        }
        break;
      }
      default: error("Unsupported SEXP type: '%s'", type2char(tx));
    }
  } else {
    if(l != length(w)) error("length(w) must match nrow(x)");
    if(TYPEOF(w) != REALSXP) {
      if(TYPEOF(w) != INTSXP && TYPEOF(w) != LGLSXP) error("weigths must be double or integer");
      w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
    }
    if(tx != REALSXP) {
      if(tx != INTSXP) error("Unsupported SEXP type: '%s'", type2char(tx));
      x = PROTECT(coerceVector(x, REALSXP)); ++nprotect;
    }
    double *px = REAL(x), *restrict pw = REAL(w), *pout = REAL(out);

    if(ng == 0) {
      if(nthreads <= 1) {
        for(int j = 0; j != col; ++j) pout[j] = fmean_weights_impl(px + j*l, pw, narm, l);
      } else if(col >= nthreads) {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) pout[j] = fmean_weights_impl(px + j*l, pw, narm, l);
      } else {
        for(int j = 0; j != col; ++j) pout[j] = fmean_weights_omp_impl(px + j*l, pw, narm, l, nthreads);
      }
    } else {
      if(nthreads <= 1 || col == 1) {
        for(int j = 0; j != col; ++j) fmean_weights_g_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l);
      } else {
        if(nthreads > col) nthreads = col;
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) fmean_weights_g_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l);
      }
    }
  }
  matCopyAttr(out, x, Rdrop, ng);
  UNPROTECT(nprotect);
  return out;
}


// For safe multithreading across data frame columns

double fmean_impl_dbl(SEXP x, int narm, int nthreads) {
  int l = length(x);
  if(l < 1) return NA_REAL;
  if(nthreads <= 1) switch(TYPEOF(x)) {
    case REALSXP: return fmean_double_impl(REAL(x), narm, l);
    case LGLSXP:
    case INTSXP: return fmean_int_impl(INTEGER(x), narm, l);
    default: error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
  }
  switch(TYPEOF(x)) {
    case REALSXP: return fmean_double_omp_impl(REAL(x), narm, l, nthreads);
    case LGLSXP:
    case INTSXP: return fmean_int_omp_impl(INTEGER(x), narm, l, nthreads);
    default: error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
  }
}

SEXP fmean_impl_SEXP(SEXP x, int narm, int nthreads) {
  return ScalarReal(fmean_impl_dbl(x, narm, nthreads));
}

double fmean_w_impl_dbl(SEXP x, double *pw, int narm, int nthreads) {
  int l = length(x);
  if(l < 1) return NA_REAL;
  if(TYPEOF(x) != REALSXP) {
    if(TYPEOF(x) != INTSXP && TYPEOF(x) != LGLSXP) error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
    x = PROTECT(coerceVector(x, REALSXP));
    double res = (nthreads <= 1) ? fmean_weights_impl(REAL(x), pw, narm, l) :
      fmean_weights_omp_impl(REAL(x), pw, narm, l, nthreads);
    UNPROTECT(1);
    return res;
  }
  return (nthreads <= 1) ? fmean_weights_impl(REAL(x), pw, narm, l) :
    fmean_weights_omp_impl(REAL(x), pw, narm, l, nthreads);
}

SEXP fmean_w_impl_SEXP(SEXP x, double *pw, int narm, int nthreads) {
  return ScalarReal(fmean_w_impl_dbl(x, pw, narm, nthreads));
}

SEXP fmean_g_impl(SEXP x, const int ng, const int *pg, const int *pgs, int narm) {
  int l = length(x);
  if(l < 1) return ScalarReal(NA_REAL);

  SEXP res = PROTECT(allocVector(REALSXP, ng));
  switch(TYPEOF(x)) {
    case REALSXP:
      fmean_double_g_impl(REAL(res), REAL(x), ng, pg, pgs, narm, l);
      break;
    case LGLSXP:
    case INTSXP:
      fmean_int_g_impl(REAL(res), INTEGER(x), ng, pg, pgs, narm, l);
      break;
    default: error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
  }

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

void fmean_g_omp_impl(SEXP x, void *pres, const int ng, const int *pg, const int *pgs, int narm) {
  switch(TYPEOF(x)) {
    case REALSXP:
      fmean_double_g_impl(pres, REAL(x), ng, pg, pgs, narm, length(x));
      break;
    case LGLSXP:
    case INTSXP:
      fmean_int_g_impl(pres, INTEGER(x), ng, pg, pgs, narm, length(x));
      break;
    default: error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
  }
}


SEXP fmean_wg_impl(SEXP x, const int ng, const int *pg, double *pw, int narm) {
  int l = length(x), nprotect = 1;
  if(l < 1) return ScalarReal(NA_REAL);

  if(TYPEOF(x) != REALSXP) {
    if(TYPEOF(x) != INTSXP && TYPEOF(x) != LGLSXP) error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
    x = PROTECT(coerceVector(x, REALSXP)); ++nprotect;
  }

  SEXP res = PROTECT(allocVector(REALSXP, ng));
  fmean_weights_g_impl(REAL(res), REAL(x), ng, pg, pw, narm, l);

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) copyMostAttrib(x, res);
  UNPROTECT(nprotect);
  return res;
}


#undef COLWISE_FMEAN_LIST
#define COLWISE_FMEAN_LIST(FUN, WFUN)                                       \
if(nwl) {                                                                  \
  if(nthreads > 1 && l >= nthreads) {                                      \
    _Pragma("omp parallel for num_threads(nthreads)")                      \
    for(int j = 0; j < l; ++j) pout[j] = FUN(px[j], narm, 1);              \
  } else {                                                                 \
    for(int j = 0; j != l; ++j) pout[j] = FUN(px[j], narm, nthreads);      \
  }                                                                        \
} else {                                                                   \
  double *restrict pw = REAL(w);                                           \
  if(nthreads > 1 && l >= nthreads) {                                      \
    _Pragma("omp parallel for num_threads(nthreads)")                      \
    for(int j = 0; j < l; ++j) pout[j] = WFUN(px[j], pw, narm, 1);         \
  } else {                                                                 \
    for(int j = 0; j != l; ++j) pout[j] = WFUN(px[j], pw, narm, nthreads); \
  }                                                                        \
}


SEXP fmeanlC(SEXP x, SEXP Rng, SEXP g, SEXP gs, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads) {
  int l = length(x), ng = asInteger(Rng), nthreads = asInteger(Rnthreads), nwl = isNull(w),
    narm = asLogical(Rnarm), nprotect = 1;

  // TODO: Disable multithreading if overall data size is small?
  if(l < 1) return x; // needed ??
  if(nthreads > max_threads) nthreads = max_threads;

  if(!nwl) {
    if(length(VECTOR_ELT(x, 0)) != length(w)) error("length(w) must match nrow(x)");
    if(TYPEOF(w) != REALSXP) {
      if(TYPEOF(w) != INTSXP && TYPEOF(w) != LGLSXP) error("weigths must be double or integer");
      w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
    }
  }

  if(ng == 0 && asLogical(Rdrop)) {
    SEXP out = PROTECT(allocVector(REALSXP, l)), *restrict px = SEXPPTR(x);
    double *restrict pout = REAL(out);
    COLWISE_FMEAN_LIST(fmean_impl_dbl, fmean_w_impl_dbl);
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(nprotect);
    return out;
  }

  SEXP out = PROTECT(allocVector(VECSXP, l)), *restrict pout = SEXPPTR(out), *restrict px = SEXPPTR(x);

  if(ng == 0) {
    COLWISE_FMEAN_LIST(fmean_impl_SEXP, fmean_w_impl_SEXP);
    // Needed because including it in an OpenMP loop together with ScalarReal() is not thread safe
    for(int j = 0; j < l; ++j) {
      SEXP xj = px[j];
      if(ATTRIB(xj) != R_NilValue && !(isObject(xj) && inherits(xj, "ts")))
        copyMostAttrib(xj, pout[j]);
    }
  } else {
    if(length(VECTOR_ELT(x, 0)) != length(g)) error("length(g) must match length(x)");
    const int *restrict pg = INTEGER(g);
    if(nthreads > l) nthreads = l;

    if(nwl) { // no weights
      int *restrict pgs = &nprotect;
      if(!narm) {
        if(length(gs) == ng) pgs = INTEGER(gs);
        else {
          SEXP gs_ = PROTECT(allocVector(INTSXP, ng)); ++nprotect;
          pgs = INTEGER(gs_);
          memset(pgs, 0, sizeof(int) * ng);
          for(int i = 0, nrx = length(g); i != nrx; ++i) ++pgs[pg[i]-1];
        }
      }

      if(nthreads > 1 && l > 1) {
        for(int j = 0; j != l; ++j) {
          SEXP xj = px[j], outj;
          SET_VECTOR_ELT(out, j, outj = allocVector(REALSXP, ng));
          if(ATTRIB(xj) != R_NilValue && !(isObject(xj) && inherits(xj, "ts"))) copyMostAttrib(xj, outj);
        }
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < l; ++j) fmean_g_omp_impl(px[j], DATAPTR(pout[j]), ng, pg, pgs, narm);
      } else {
        for(int j = 0; j != l; ++j) pout[j] = fmean_g_impl(px[j], ng, pg, pgs, narm);
      }
    } else {
      double *restrict pw = REAL(w);
      if(nthreads > 1 && l > 1) {
        int nrx = length(g);
        for(int j = 0, dup = 0; j != l; ++j) {
          SEXP xj = px[j], outj;
          SET_VECTOR_ELT(out, j, outj = allocVector(REALSXP, ng));
          if(ATTRIB(xj) != R_NilValue && !(isObject(xj) && inherits(xj, "ts"))) copyMostAttrib(xj, outj);
          if(TYPEOF(xj) != REALSXP) {
            if(TYPEOF(xj) != INTSXP && TYPEOF(xj) != LGLSXP) error("Unsupported SEXP type: '%s'", type2char(TYPEOF(xj)));
            if(dup == 0) {x = PROTECT(shallow_duplicate(x)); ++nprotect; px = SEXPPTR(x); dup = 1;}
            SET_VECTOR_ELT(x, j, coerceVector(xj, REALSXP));
          }
        }
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < l; ++j) fmean_weights_g_impl(REAL(pout[j]), REAL(px[j]), ng, pg, pw, narm, nrx);
      } else {
        for(int j = 0; j != l; ++j) pout[j] = fmean_wg_impl(px[j], ng, pg, pw, narm);
      }
    }
  }

  DFcopyAttr(out, x, ng);
  UNPROTECT(nprotect);
  return out;
}

