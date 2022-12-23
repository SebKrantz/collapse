#ifdef _OPENMP
#include <omp.h>
#endif
#include "collapse_c.h"

// Inspired by Numerical Recipes in C and data.table's quickselect.c
// Additional inspiration taken from Rfast2, and these references for sample Quantiles:
// https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample
// https://doi.org/10.2307/2684934
// https://aakinshin.net/posts/weighted-quantiles/

// Adopted from data.table's quickselect.c
static inline void iswap(int *a, int *b)           {int     tmp=*a; *a=*b; *b=tmp;}
static inline void dswap(double *a, double *b)     {double  tmp=*a; *a=*b; *b=tmp;}

// Barebones quickselect algorithm from Numerical Recipes in C
#undef QUICKSELECT
#define QUICKSELECT(SWAP)                     \
  unsigned int ir = n-1, l = 0, lp;           \
  for(;;) {                                   \
    lp = l+1;                                 \
    if (ir <= lp) {                           \
      if (ir == lp && x[ir] < x[l]) {         \
        SWAP(x+l, x+ir);                      \
      }                                       \
      break;                                  \
    } else {                                  \
      unsigned int mid=(l+ir) >> 1;           \
      SWAP(x+mid, x+lp);                      \
      if (x[l] > x[ir]) {                     \
        SWAP(x+l, x+ir);                      \
      }                                       \
      if (x[lp] > x[ir]) {                    \
        SWAP(x+lp, x+ir);                     \
      }                                       \
      if (x[l] > x[lp]) {                     \
        SWAP(x+l, x+lp);                      \
      }                                       \
      unsigned int i=lp, j=ir;                \
      a=x[lp];                                \
      for (;;) {                              \
        do i++; while (x[i] < a);             \
        do j--; while (x[j] > a);             \
        if (j < i) break;                     \
        SWAP(x+i, x+j);                       \
      }                                       \
      x[lp]=x[j];                             \
      x[j]=a;                                 \
      if (j >= elem) ir=j-1;                  \
      if (j <= elem) l=i;                     \
    }                                         \
  }                                           \
  a = x[elem];

// Quantile method switcher
// https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample
// Need to subtract 1 from h because of 0-indexing in C
#undef RETQSWITCH
#define RETQSWITCH(n)                                          \
  switch(ret) {                                                \
  case 7:                                                      \
  case 1:                                                      \
  case 2: /* quantile type 7, average, or Lower element*/      \
    h = (n - 1)*Q;                                             \
    break;                                                     \
  case 3: /* upper element*/                                   \
    h = n*Q;                                                   \
    break;                                                     \
  case 4: /* quantile type 4*/                                 \
    h = n*Q - 1.0;                                             \
    break;                                                     \
  case 5: /* quantile type 5*/                                 \
    h = n*Q - 0.5;                                             \
    break;                                                     \
  case 6: /* quantile type 6*/                                 \
    h = (n + 1)*Q - 1.0;                                       \
    break;                                                     \
  case 8: /* quantile type 8 (best according to H&F 1986)*/    \
    h = ((double)n + 1.0/3.0)*Q - 2.0/3.0;                     \
    break;                                                     \
  case 9: /* quantile type 9*/                                 \
    h = ((double)n + 1.0/4.0)*Q - 5.0/8.0;                     \
    break;                                                     \
  }



double dquickselect(double *x, const int n, const int ret, const double Q) {
  if(n == 0) return NA_REAL;
  unsigned int elem;
  double a, b, h;
  RETQSWITCH(n);
  elem = h;
  h -= elem; // Key: need to subtract elem
  QUICKSELECT(dswap);
  if((ret == 1 && n%2 == 1) || ret == 2 || ret == 3 || h == 0.0) return a;
  b = x[elem+1];
  for(int i = elem+2; i < n; ++i) if(x[i] < b) b = x[i];
  if(ret == 1 || Q == 0.5) return (a+b)/2.0;
  return a + h*(b-a); // same as (1-h)*a + h*b
}

double iquickselect(int *x, const int n, const int ret, const double Q) {
  if(n == 0) return NA_REAL;
  unsigned int elem;
  int a, b;
  double h;
  RETQSWITCH(n);
  elem = h;
  h -= elem; // Key: need to subtract elem
  QUICKSELECT(iswap);
  if((ret == 1 && n%2 == 1) || ret == 2 || ret == 3 || h == 0.0) return (double)a;
  b = x[elem+1];
  for(int i = elem+2; i < n; ++i) if(x[i] < b) b = x[i];
  if(ret == 1 || Q == 0.5) return ((double)a+(double)b)/2.0;
  return (double)a + h*(double)(b-a); // same as (1-h)*(double)a + h*(double)b
}


// --------------------------------------------------------------------------
// First a faster quantile function
// --------------------------------------------------------------------------

// Need versions that to supply the element and h

double dquickselect_elem(double *x, const int n, const unsigned int elem, double h) {
  if(n == 0) return NA_REAL;
  double a, b;
  QUICKSELECT(dswap);
  if(h == 0.0) return a;
  b = x[elem+1];
  for(int i = elem+2; i < n; ++i) if(x[i] < b) b = x[i];
  return a + h*(b-a);
}

double iquickselect_elem(int *x, const int n, const unsigned int elem, double h) {
  if(n == 0) return NA_REAL;
  int a, b;
  QUICKSELECT(iswap);
  if(h == 0.0) return (double)a;
  b = x[elem+1];
  for(int i = elem+2; i < n; ++i) if(x[i] < b) b = x[i];
  return (double)a + h*(double)(b-a);
}

// If minimum and / or maximum requested
// TODO: could do this afterwards?? in selected area??
// exact comparison??
// TODO: Bad probabilities check.
// Also: For large N radix ordering seems to be faster than repeated quickselect
#define FQUANTILE_CORE(QFUN)                                                  \
  if(probs[0] == 0 || probs[np-1] == 1) {                                     \
    x_min = x_max = x_cc[0];                                                  \
    if(probs[0] == 0 && probs[np-1] == 1) {                                   \
      for(unsigned int i = 1; i != l; ++i) {                                  \
        if(x_cc[i] < x_min) x_min = x_cc[i];                                  \
        if(x_cc[i] > x_max) x_max = x_cc[i];                                  \
      }                                                                       \
      pres[0] = x_min; pres[np-1] = x_max;                                    \
    } else if(probs[0] == 0) {                                                \
      for(unsigned int i = 1; i != l; ++i)                                    \
        if(x_cc[i] < x_min) x_min = x_cc[i];                                  \
      pres[0] = x_min;                                                        \
    } else {                                                                  \
      for(unsigned int i = 1; i != l; ++i)                                    \
        if(x_cc[i] > x_min) x_min = x_cc[i];                                  \
      pres[np-1] = x_max;                                                     \
    }                                                                         \
  }                                                                           \
  double h, Q;                                                                \
  for(int i = 0, offset = 0, ih; i < np; ++i) {                               \
    Q = probs[i];                                                             \
    if(Q > 0 && Q < 1) {                                                      \
      RETQSWITCH(l);                                                          \
      ih = h;                                                                 \
      pres[i] = QFUN(x_cc + offset, l - offset, ih - offset, h - ih);         \
      offset = ih;                                                            \
    }                                                                         \
  }


#define FQUANTILE_ORDVEC                                    \
double a, b, h, Q;                                          \
for(int i = 0, ih; i < np; ++i) {                           \
  Q = probs[i];                                             \
  if(Q > 0 && Q < 1) {                                      \
    RETQSWITCH(l);                                          \
    ih = h; a = px[po[ih]]; b = px[po[ih+1]];               \
    pres[i] = a + (h - ih) * (b - a);                       \
  } else pres[i] = px[po[(l-1)*(int)Q]];                    \
}

// TODO: Proper quantile weighting...
// See: https://aakinshin.net/posts/weighted-quantiles/
// And: https://en.wikipedia.org/wiki/Percentile#Weighted_percentile
#define WQUANTILE_CORE                             \
for(unsigned int i = 0, k = 0; i < np; ++i) {      \
  if(probs[i] > 0 && probs[i] < 1) {               \
    double wsumQ = sumw * probs[i];                \
    while(wsum < wsumQ) wsum += pw[po[k++]];       \
    if(wsum == wsumQ) {                            \
      double out = px[po[k-1]], n = 2;             \
      while(pw[po[k]] == 0) {                      \
        out += px[po[k++]];                        \
        ++n;                                       \
      }                                            \
      pres[i] = (out + px[po[k]]) / n;             \
    } else {                                       \
      double h = (wsum - wsumQ) / pw[po[k-1]],     \
             a = px[po[k-2]], b = px[po[k-1]];     \
      pres[i] = a + h * (b - a);                   \
    }                                              \
  } else pres[i] = px[po[(l-1)*(int)probs[i]]];    \
}

// TODO: args o and check.o, to compute quantiles from existing ordering vector...
// TODO: check o ??
// potentially also internal optimization ???
SEXP fquantileC(SEXP x, SEXP Rprobs, SEXP w, SEXP o, SEXP Rnarm, SEXP Rtype, SEXP checko) {

  if(TYPEOF(Rprobs) != REALSXP) error("probs needs to be a numeric vector");
  int tx = TYPEOF(x), n = length(x), np = length(Rprobs), narm = asLogical(Rnarm), ret = asInteger(Rtype), nprotect = 1;
  if(tx != REALSXP && tx != INTSXP && tx != LGLSXP) error("x needs to be numeric");

  SEXP res = PROTECT(allocVector(REALSXP, np));
  if(np == 0) { // quantile(x, numeric(0))
    UNPROTECT(nprotect);
    return res;
  }

  double *probs = REAL(Rprobs), *pres = REAL(res);
  unsigned int l = 0;

  /* TODO: switch to radixorder based on rule e.g. length(probs) > log(n) as in Rfast??
     Appears better than building checks into the arguments in R, though also decreases flexibility...
     Need to think: anything against that? and case where this is undesirable? */
  for(int i = 0; i < np; ++i) {
    if(probs[i] < 0 || probs[i] > 1) error("probabilities need to be in in range [0, 1]");
    if(i > 0 && probs[i] < probs[i-1]) error("probabilities need to be passed in ascending order");
  }

  // TODO: What about l == 0 or 1 and narm = TRUE, also with weighted...
  if(isNull(w) && isNull(o)) { // Standard: quickselect

    if(tx == REALSXP) { // Numeric data
      double *x_cc = (double *) R_alloc(n, sizeof(double)), *px = REAL(x), x_min, x_max;
      if(narm) {
        for(unsigned int i = 0; i != n; ++i) if(NISNAN(px[i])) x_cc[l++] = px[i];
        if(l <= 1) {
          for(int i = 0; i < np; ++i) pres[i] = l == 0 ? NA_REAL : x_cc[0];
          UNPROTECT(nprotect);
          return res;
        }
      } else {
        l = n;
        memcpy(x_cc, px, sizeof(double) * n);
      }
      FQUANTILE_CORE(dquickselect_elem);
    } else { // Integers
      int *x_cc = (int *) R_alloc(n, sizeof(int)), *px = INTEGER(x), x_min, x_max;
      if(narm) {
        for(unsigned int i = 0; i != n; ++i) if(px[i] != NA_INTEGER) x_cc[l++] = px[i];
        if(l <= 1) {
          for(int i = 0, ret = l == 0 ? NA_INTEGER : x_cc[0]; i < np; ++i) pres[i] = ret;
          UNPROTECT(nprotect);
          return res;
        }
      } else {
        l = n;
        memcpy(x_cc, px, sizeof(int) * n);
      }
      FQUANTILE_CORE(iquickselect_elem);
    }

  } else { // Weighted or Ordered

    int *po = &n;
    double *pw = probs, nanw0 = 0;

    if(!isNull(o)) {
      if(length(o) != n || TYPEOF(o) != INTSXP) error("o must be a valid ordering vector, of the same length as x and type integer");
      po = INTEGER(o);
      if(asLogical(checko)) { // TODO: Better way??
        for(unsigned int i = 0; i != n; ++i)
          if(po[i] < 1 || po[i] > n) error("Some elements in o are outside of range [1, length(x)]");
      }
    } else  {
      po = (int *) R_alloc(n, sizeof(int)); // Calloc ??
      num1radixsort(po, TRUE, FALSE, x);
    }

    if(!isNull(w)) {
      if(length(w) != n) error("length(w) must match length(x)");
      if(TYPEOF(w) != REALSXP) {
        if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double)");
        SEXP wd = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
        pw = REAL(wd)-1;
      } else pw = REAL(w)-1;
      nanw0 = pw[po[0]];
    }

    l = n;

    if(narm) {
      if(tx == REALSXP) { // Numeric data
        double *px = REAL(x)-1;
        if(ISNAN(px[po[0]])) error("Found missing value at the beginning of the sample. Please use option na.last = TRUE (the default) when creasting ordering vectors for use with fquantile().");
        --po; while(ISNAN(px[po[l]]) && l != 0) --l; ++po;
        if(l <= 1) {
          double ret = (l == 0 || ISNAN(nanw0)) ? NA_REAL : px[po[0]];
          for(int i = 0; i < np; ++i) pres[i] = ret;
          UNPROTECT(nprotect);
          return res;
        }
      } else {
        int *px = INTEGER(x)-1;
        if(px[po[0]] == NA_INTEGER) error("Found missing value at the beginning of the sample. Please use option na.last = TRUE (the default) when creasting ordering vectors for use with fquantile().");
        --po; while(px[po[l]] == NA_INTEGER && l != 0) --l; ++po;
        if(l <= 1) {
          for(int i = 0, ret = (l == 0 || ISNAN(nanw0)) ? NA_INTEGER : px[po[0]]; i < np; ++i) pres[i] = ret;
          UNPROTECT(nprotect);
          return res;
        }
      }
    }

    if(isNull(w)) {

      if(tx == REALSXP) { // Numeric data
        double *px = REAL(x)-1;
        FQUANTILE_ORDVEC;
      } else {
        int *px = INTEGER(x)-1;
        FQUANTILE_ORDVEC;
      }

    } else {
      double wsum = 0, sumw = 0; // wsum is running sum, sumw is the total sum

      for(unsigned int i = 0; i != l; ++i) sumw += pw[po[i]];
      if(ISNAN(sumw)) error("Missing weights in order statistics are currently only supported if x is also missing");

      if(tx == REALSXP) { // Numeric data
        double *px = REAL(x)-1;
        WQUANTILE_CORE;
      } else {
        int *px = INTEGER(x)-1;
        WQUANTILE_CORE;
      }

    }

  }

  UNPROTECT(nprotect);
  return res;
}

// TODO: quicksort for weighted statistics??

// C-implementations for different data types, parallelizable ----------------------------------

double nth_int(const int *px, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l == 1) return sorted ? px[0] : px[po[0]-1];

  int *x_cc = (int *) Calloc(l, int), n = 0;
  if(sorted) {
    if(narm) {
      for(int i = 0; i != l; ++i) if(px[i] != NA_INTEGER) x_cc[n++] = px[i];
    } else {
      n = l;
      memcpy(x_cc, px, l * sizeof(int));
    }
  } else {
    const int *pxm = px-1; // creating offset pointer to x
    if(narm) {
      for(int i = 0; i != l; ++i) if(pxm[po[i]] != NA_INTEGER) x_cc[n++] = pxm[po[i]];
    } else {
      n = l;
      for(int i = 0; i != l; ++i) x_cc[i] = pxm[po[i]];
    }
  }

  double res = iquickselect(x_cc, n, ret, Q);
  Free(x_cc);
  return res;
}

int w_nth_int(const int *restrict px, const double *restrict pw, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l == 1) {
    if(sorted) return ISNAN(pw[0]) ? NA_INTEGER : px[0];
    return ISNAN(pw[po[0]-1]) ? NA_INTEGER : px[po[0]-1];
  }
  return 1;
}

double nth_double(const double *restrict px, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l == 1) return sorted ? px[0] : px[po[0]-1];

  double *x_cc = (double *) Calloc(l, double);
  int n = 0;

  if(sorted) {
    if(narm) {
      for(int i = 0; i != l; ++i) if(NISNAN(px[i])) x_cc[n++] = px[i];
    } else {
      n = l;
      memcpy(x_cc, px, l * sizeof(double));
    }
  } else {
    const double *pxm = px-1;
    if(narm) {
      for(int i = 0; i != l; ++i) if(NISNAN(pxm[po[i]])) x_cc[n++] = pxm[po[i]];
    } else {
      n = l;
      for(int i = 0; i != l; ++i) x_cc[i] = pxm[po[i]];
    }
  }

  double res = dquickselect(x_cc, n, ret, Q);
  Free(x_cc);
  return res;
}

double w_nth_double(const double *restrict px, const double *restrict pw, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l == 1) {
    if(sorted) return ISNAN(pw[0]) ? NA_REAL : px[0];
    return ISNAN(pw[po[0]-1]) ? NA_REAL : px[po[0]-1];
  }
  // num1radixsort(ord, TRUE, FALSE, px);

  // TODO: use quicksort by group?? potentially can sort half the array?? But compare against radixsort by group and overall radix sort
  return 1.0;
}


// Implementations for R vectors -----------------------------------------------

SEXP nth_impl(SEXP x, int narm, int ret, double Q) {
  int l = length(x);
  if(l <= 1) return x;

  SEXP res;
  switch(TYPEOF(x)) {
  case REALSXP:
    PROTECT(res = ScalarReal(nth_double(REAL(x), &l, l, 1, narm, ret, Q)));
    break;
  case INTSXP:
  case LGLSXP:
    PROTECT(res = ScalarInteger(nth_int(INTEGER(x), &l, l, 1, narm, ret, Q)));
  default: error("Not Supported SEXP Type!");
  }

  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

SEXP w_nth_impl(SEXP x, double *pw, int narm, int ret, double Q) {
  int l = length(x);
  if(l <= 1) return x;

  SEXP res;
  switch(TYPEOF(x)) {
  case REALSXP:
    PROTECT(res = ScalarReal(w_nth_double(REAL(x), pw, &l, l, 1, narm, ret, Q)));
    break;
  case INTSXP:
  case LGLSXP:
    PROTECT(res = ScalarInteger(w_nth_int(INTEGER(x), pw, &l, l, 1, narm, ret, Q)));
  default: error("Not Supported SEXP Type!");
  }

  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

SEXP nth_g_impl(SEXP x, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int ret, double Q, int nthreads) {

  int l = length(x), tx = TYPEOF(x);
  if(nthreads > ng) nthreads = ng;

  SEXP res = PROTECT(allocVector(tx, ng));

  if(sorted) { // Sorted: could compute cumulative group size (= starts) on the fly... but doesn't work multithreaded...
    po = &l;
    switch(tx) {
    case REALSXP: {
      double *px = REAL(x), *pres = REAL(res);
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = pgs[gr] == 0 ? NA_REAL : nth_double(px + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
      break;
    }
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x), *pres = INTEGER(res);
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = pgs[gr] == 0 ? NA_INTEGER : nth_int(px + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
      break;
    }
    default: error("Not Supported SEXP Type!");
    }
  } else { // Not sorted. Perhaps reordering x is faster??
    switch(tx) {
    case REALSXP: {
      double *px = REAL(x), *pres = REAL(res);
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = pgs[gr] == 0 ? NA_REAL : nth_double(px, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
      break;
    }
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x), *pres = INTEGER(res);
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = pgs[gr] == 0 ? NA_INTEGER : nth_int(px, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
      break;
    }
    default: error("Not Supported SEXP Type!");
    }
  }

  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

SEXP w_nth_g_impl(SEXP x, double *pw, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int ret, double Q, int nthreads) {

  int l = length(x), tx = TYPEOF(x);
  if(nthreads > ng) nthreads = ng;

  SEXP res = PROTECT(allocVector(tx, ng));

  if(sorted) { // Sorted: could compute cumulative group size (= starts) on the fly... but doesn't work multithreaded...
    po = &l;
    switch(tx) {
    case REALSXP: {
      double *px = REAL(x), *pres = REAL(res);
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = pgs[gr] == 0 ? NA_REAL : w_nth_double(px + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
      break;
    }
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x), *pres = INTEGER(res);
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = pgs[gr] == 0 ? NA_INTEGER : w_nth_int(px + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
      break;
    }
    default: error("Not Supported SEXP Type!");
    }
  } else { // Not sorted. Perhaps reordering x is faster??
    switch(tx) {
    case REALSXP: {
      double *px = REAL(x), *pres = REAL(res);
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = pgs[gr] == 0 ? NA_REAL : w_nth_double(px, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
      break;
    }
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x), *pres = INTEGER(res);
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = pgs[gr] == 0 ? NA_INTEGER : w_nth_int(px, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
      break;
    }
    default: error("Not Supported SEXP Type!");
    }
  }

  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}


// Functions for Export --------------------------------------------------------

SEXP fnthC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rret, SEXP Rnthreads) {
  int nullg = isNull(g), nullw = isNull(w), l = length(x), narm = asLogical(Rnarm),
      ret = asInteger(Rret), nprotect = 0;
  if(l <= 1) return x;
  if(length(p) != 1) error("fnth supports only a single element / quantile. Use fquantile for multiple quantiles.");
  double Q = asReal(p);
  if(ISNAN(Q) || Q <= 0 || Q == 1) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
  if(Q > 1) {
    ret = 2; // Correct ??
    if(nullg) {
      if(Q >= l) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(l-1);
    } else {
      if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
      int ng = length(VECTOR_ELT(g, 2));
      if(Q >= l/ng) error("n needs to be between 0 and 1, or between 1 and the length(x)/ng, with ng the number of groups. Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(l/ng-1);
    }
  }
  if(nullg && nullw) return nth_impl(x, narm, ret, Q);
  double tmp = 0.0, *restrict pw = &tmp;
  if(!nullw) {
    if(length(w) != l) error("length(w) must match length(x)");
    if(TYPEOF(w) != REALSXP) {
      if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double)");
      SEXP wd = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
      pw = REAL(wd);
    } else pw = REAL(w);
  }
  if(nullg) {
    UNPROTECT(nprotect);
    return w_nth_impl(x, pw, narm, ret, Q);
  }
  if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
  const SEXP *restrict pg = SEXPPTR(g), o = pg[6];
  int sorted = LOGICAL(pg[5])[1] == 1, ng = INTEGER(pg[0])[0], *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst;
  if(l != length(pg[1])) error("length(g) must match length(x)");
  if(isNull(o)) {
    int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;
    for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
    pst = cgs + 1;
    if(sorted) po = &l;
    else {
      int *restrict count = (int *) Calloc(ng+1, int);
      po = (int *) R_alloc(l, sizeof(int)); --po;
      for(int i = 0; i != l; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
      ++po; Free(count);
    }
  } else {
    po = INTEGER(o);
    pst = INTEGER(getAttrib(o, install("starts")));
  }
  SEXP res;
  if(nullw) res = nth_g_impl(x, ng, pgs, po, pst, sorted, narm, ret, Q, asInteger(Rnthreads));
  else res = w_nth_g_impl(x, pw, ng, pgs, po, pst, sorted, narm, ret, Q, asInteger(Rnthreads));
  UNPROTECT(nprotect);
  return res;
}

// TODO: allow column-level parallelism??
SEXP fnthlC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rret, SEXP Rnthreads) {
  int nullg = isNull(g), nullw = isNull(w), l = length(x), ng = 0, nprotect = 1,
    narm = asLogical(Rnarm), drop = asLogical(Rdrop), ret = asInteger(Rret), nthreads = asInteger(Rnthreads);
  if(l < 1) return x;

  double Q = asReal(p);
  if(ISNAN(Q) || Q <= 0 || Q == 1) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
  if(Q > 1) {
    ret = 2; // Correct ??
    if(nullg) {
      if(Q >= l) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(l-1);
    } else {
      if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
      int ng = length(VECTOR_ELT(g, 2));
      if(Q >= l/ng) error("n needs to be between 0 and 1, or between 1 and the length(x)/ng, with ng the number of groups. Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(l/ng-1);
    }
  }

  SEXP out = PROTECT(allocVector(nullg && drop ? REALSXP : VECSXP, l)), *restrict px = SEXPPTR(x);
  int nrx = length(px[0]);
  double tmp = 0.0, *restrict pw = &tmp;

  if(!nullw) {
    if(length(w) != nrx) error("length(w) must match nrow(x)");
    if(TYPEOF(w) != REALSXP) {
      if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double)");
      SEXP wd = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
      pw = REAL(wd);
    } else pw = REAL(w);
  }

  if(nullg) {
    if(nthreads > l) nthreads = l;
    if(drop) {
      double *restrict pout = REAL(out);
      if(nullw) {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < l; ++j) pout[j] = REAL(nth_impl(px[j], narm, ret, Q))[0];
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < l; ++j) pout[j] = REAL(w_nth_impl(px[j], pw, narm, ret, Q))[0];
      }
      setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
      UNPROTECT(nprotect);
      return out;
    }
    SEXP *restrict pout = SEXPPTR(out);
    if(nullw) {
      #pragma omp parallel for num_threads(nthreads)
      for(int j = 0; j < l; ++j) pout[j] = nth_impl(px[j], narm, ret, Q);
    } else {
      #pragma omp parallel for num_threads(nthreads)
      for(int j = 0; j < l; ++j) pout[j] = w_nth_impl(px[j], pw, narm, ret, Q);
    }
  } else {
    if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
    const SEXP *restrict pg = SEXPPTR(g), o = pg[6];
    ng = INTEGER(pg[0])[0];
    int sorted = LOGICAL(pg[5])[1] == 1, *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst;
    if(nrx != length(pg[1])) error("length(g) must match nrow(x)");
    if(isNull(o)) {
      int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;
      for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
      pst = cgs + 1;
      if(sorted) po = &l;
      else {
        int *restrict count = (int *) Calloc(ng+1, int);
        po = (int *) R_alloc(nrx, sizeof(int)); --po;
        for(int i = 0; i != nrx; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
        ++po; Free(count);
      }
    } else {
      po = INTEGER(o);
      pst = INTEGER(getAttrib(o, install("starts")));
    }
    SEXP *restrict pout = SEXPPTR(out);
    if(nullw) { // Parallelism at sub-column level
      for(int j = 0; j < l; ++j) pout[j] = nth_g_impl(px[j], ng, pgs, po, pst, sorted, narm, ret, Q, nthreads);
    } else { // Parallelism at sub-column level
      for(int j = 0; j < l; ++j) pout[j] = w_nth_g_impl(px[j], pw, ng, pgs, po, pst, sorted, narm, ret, Q, nthreads);
    }
  }

  DFcopyAttr(out, x, ng);
  UNPROTECT(nprotect);
  return out;
}

SEXP fnthmC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rret, SEXP Rnthreads) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1],
      narm = asLogical(Rnarm), ret = asInteger(Rret), nthreads = asInteger(Rnthreads),
      nullg = isNull(g), nullw = isNull(w), nprotect = 1;
  if(l <= 1) return x; // Prevents seqfault for numeric(0) #101
  if(nthreads > col) nthreads = col;
  if(length(p) != 1) error("fnth supports only a single element / quantile. Use fquantile for multiple quantiles.");
  double Q = asReal(p);
  if(ISNAN(Q) || Q <= 0 || Q == 1) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
  if(Q > 1) {
    ret = 2; // Correct ??
    if(nullg) {
      if(Q >= l) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(l-1);
    } else {
      if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
      int ng = length(VECTOR_ELT(g, 2));
      if(Q >= l/ng) error("n needs to be between 0 and 1, or between 1 and the length(x)/ng, with ng the number of groups. Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(l/ng-1);
    }
  }

  double tmp = 0.0, *restrict pw = &tmp;
  if(!nullw) {
    if(length(w) != l) error("length(w) must match nrow(x)");
    if(TYPEOF(w) != REALSXP) {
      if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double)");
      SEXP wd = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
      pw = REAL(wd);
    } else pw = REAL(w);
  }

  if(nullg) {
    SEXP res = PROTECT(allocVector(tx, col));

    switch(tx) {
    case REALSXP: {
      double *px = REAL(x), *restrict pres = REAL(res);
      if(nullw) {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) pres[j] = nth_double(px + j*l, &l, l, 1, narm, ret, Q);
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) pres[j] = w_nth_double(px + j*l, pw, &l, l, 1, narm, ret, Q);
      }
      break;
    }
    case INTSXP:
    case LGLSXP: {  // Factor matrix not well defined object...
      int *px = INTEGER(x), *restrict pres = INTEGER(res);
      if(nullw) {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) pres[j] = nth_int(px + j*l, &l, l, 1, narm, ret, Q);
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) pres[j] = w_nth_int(px + j*l, pw, &l, l, 1, narm, ret, Q);
      }
      break;
    }
    default: error("Not Supported SEXP Type!");
    }

    matCopyAttr(res, x, Rdrop, /*ng=*/0);
    UNPROTECT(nprotect);
    return res;
  }

  // With groups
  if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
  const SEXP *restrict pg = SEXPPTR(g), o = pg[6];
  int sorted = LOGICAL(pg[5])[1] == 1, ng = INTEGER(pg[0])[0], *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst, gl = length(pg[1]);
  if(l != gl) error("length(g) must match nrow(x)");
  SEXP res = PROTECT(allocVector(tx, ng * col));

  if(isNull(o)) {
    int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;
    for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
    pst = cgs + 1;
    if(sorted) po = &l;
    else {
      int *restrict count = (int *) Calloc(ng+1, int);
      po = (int *) R_alloc(l, sizeof(int)); --po;
      for(int i = 0; i != l; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
      ++po; Free(count);
    }
  } else {
    po = INTEGER(o);
    pst = INTEGER(getAttrib(o, install("starts")));
  }

  if(sorted) { // Sorted
    switch(tx) {
    case REALSXP: {
      double *px = REAL(x), *restrict pres = REAL(res);
      if(nullw) {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) {
          int jng = j * ng;
          double *pxj = px + j * l;
          for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : nth_double(pxj + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
        }
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) {
          int jng = j * ng;
          double *pxj = px + j * l;
          for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : w_nth_double(pxj + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
        }
      }
      break;
    }
    case INTSXP:
    case LGLSXP: { // Factor matrix not well defined object...
      int *px = INTEGER(x), *restrict pres = INTEGER(res);
      if(nullw) {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) {
          int *pxj = px + j * l, jng = j * ng;
          for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : nth_int(pxj + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
        }
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) {
          int *pxj = px + j * l, jng = j * ng;
          for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : w_nth_int(pxj + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
        }
      }
      break;
    }
    default: error("Not Supported SEXP Type!");
    }
  } else { // Not sorted
    switch(tx) {
    case REALSXP: {
      double *px = REAL(x), *restrict pres = REAL(res);
      if(nullw) {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) {
          int jng = j * ng;
          double *pxj = px + j * l;
          for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : nth_double(pxj, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
        }
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) {
          int jng = j * ng;
          double *pxj = px + j * l;
          for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : w_nth_double(pxj, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
        }
      }
      break;
    }
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x), *restrict pres = INTEGER(res);
      if(nullw) {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) {
          int jng = j * ng, *pxj = px + j * l;
          for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : nth_int(pxj, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
        }
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) {
          int jng = j * ng, *pxj = px + j * l;
          for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : w_nth_int(pxj, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
        }
      }
      break;
    }
    default: error("Not Supported SEXP Type!");
    }
  }

  matCopyAttr(res, x, Rdrop, ng);
  UNPROTECT(nprotect);
  return res;
}

