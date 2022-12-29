#ifdef _OPENMP
#include <omp.h>
#endif
#include "collapse_c.h"

// Inspired by Numerical Recipes in C and data.table's quickselect.c,
// R's quantile() function, Rfast2::Quantile, and these references for sample quantiles:
// https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample
// https://doi.org/10.2307/2684934
// https://aakinshin.net/posts/weighted-quantiles/
// At large, the weighted quantile algorithm is my own cooking [(C) 2022 Sebastian Krantz]

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
/*
case 4: // quantile type 4: not exact for weighted case, and also bad properties.
  h = n*Q - 1.0;
  break;
*/

// Weighted quantiles. Idea: make h dependent on sum of weights and average weight.
#undef RETWQSWITCH
#define RETWQSWITCH(sumw, mu)                                  \
switch(ret) {                                                  \
case 7:                                                        \
case 1:                                                        \
case 2: /* quantile type 7, average, or Lower element*/        \
  h = (sumw - mu) * Q;                                         \
  break;                                                       \
case 3: /* upper element*/                                     \
  h = sumw * Q;                                                \
  break;                                                       \
case 5: /* quantile type 5*/                                   \
  h = sumw * Q - 0.5*mu;                                       \
  if(h < 0.0) h = 0.0;                                         \
  break;                                                       \
case 6: /* quantile type 6*/                                   \
  h = (sumw + mu)*Q - mu;                                      \
  if(h < 0.0) h = 0.0;                                         \
  break;                                                       \
case 8: /* quantile type 8 (best according to H&F 1986)*/      \
  h = (sumw + 1.0/3.0 * mu)*Q - 2.0/3.0 * mu;                  \
  if(h < 0.0) h = 0.0;                                         \
  break;                                                       \
case 9: /* quantile type 9*/                                   \
  h = (sumw + 1.0/4.0 * mu)*Q - 5.0/8.0 * mu;                  \
  if(h < 0.0) h = 0.0;                                         \
  break;                                                       \
}
/*
case 4: // quantile type 4: does not give exact results
  h = sumw * Q - mu;
  break;
*/
/* Redundant? -> yes!
 if(h > sumw) h = sumw;
*/

// --------------------------------------------------------------------------
// First a faster quantile function
// --------------------------------------------------------------------------

// Need versions that supply the element and h

double dquickselect_elem(double *x, const int n, const unsigned int elem, double h) {
  // if(n == 0) return NA_REAL; // done in fquantile...
  double a, b;
  QUICKSELECT(dswap);
  if(elem == n-1 || h <= 0.0) return a;
  b = x[elem+1];
  for(int i = elem+2; i < n; ++i) if(x[i] < b) b = x[i];
  return a + h*(b-a);
}

double iquickselect_elem(int *x, const int n, const unsigned int elem, double h) {
  // if(n == 0) return NA_REAL; // done in fquantile...
  int a, b;
  QUICKSELECT(iswap);
  if(elem == n-1 || h <= 0.0) return (double)a;
  b = x[elem+1];
  for(int i = elem+2; i < n; ++i) if(x[i] < b) b = x[i];
  return (double)a + h*(double)(b-a);
}

/*
 Old Solution: full initial pass to get min and max.
if(probs[0] == 0.0 || probs[np-1] == 1.0) {                                   \
  x_min = x_max = x_cc[0];                                                    \
  if(probs[0] == 0.0 && probs[np-1] == 1.0) {                                 \
    for(unsigned int i = 1; i != l; ++i) {                                    \
      if(x_cc[i] < x_min) x_min = x_cc[i];                                    \
      if(x_cc[i] > x_max) x_max = x_cc[i];                                    \
    }                                                                         \
    pres[0] = x_min; pres[np-1] = x_max;                                      \
  } else if(probs[0] == 0.0) {                                                \
    for(unsigned int i = 1; i != l; ++i)                                      \
      if(x_cc[i] < x_min) x_min = x_cc[i];                                    \
      pres[0] = x_min;                                                        \
  } else {                                                                    \
    for(unsigned int i = 1; i != l; ++i)                                      \
      if(x_cc[i] > x_min) x_min = x_cc[i];                                    \
      pres[np-1] = x_max;                                                     \
  }                                                                           \
}                                                                             \
*/
#define FQUANTILE_CORE(QFUN)                                                  \
  double h, Q;                                                                \
  int ih;                                                                     \
  for(int i = 0, offset = 0; i < np; ++i) {                                   \
    Q = probs[i];                                                             \
    if(Q > 0.0 && Q < 1.0) {                                                  \
      RETQSWITCH(l);                                                          \
      ih = h;                                                                 \
      pres[i] = QFUN(x_cc + offset, l - offset, ih - offset, h - ih);         \
      offset = ih;                                                            \
    }                                                                         \
  } /* This is much more efficient: fetching min and max ex-post */           \
  if(probs[0] == 0.0) {                                                       \
    x_min = x_cc[0];                                                          \
    for(unsigned int i = 0, end = l*probs[1]; i < end; ++i)                   \
        if(x_cc[i] < x_min) x_min = x_cc[i];                                  \
    pres[0] = (double)x_min;                                                  \
  }                                                                           \
  if(probs[np-1] == 1.0) {                                                    \
    x_max = x_cc[ih];                                                         \
    for(unsigned int i = ih+1; i < l; ++i)                                    \
        if(x_cc[i] > x_max) x_max = x_cc[i];                                  \
    pres[np-1] = (double)x_max;                                               \
  }

// If we have an ordering vector supplied as input to the function
#define FQUANTILE_ORDVEC                                    \
double a, b, h, Q;                                          \
for(int i = 0, ih; i < np; ++i) {                           \
  Q = probs[i];                                             \
  if(Q > 0.0 && Q < 1.0) {                                  \
    RETQSWITCH(l);                                          \
    ih = h; a = px[po[ih]];                                 \
    if(ih == n-1 || h <= 0.0) pres[i] = a;                  \
    else {                                                  \
      b = px[po[ih+1]];                                     \
      pres[i] = a + (h - ih) * (b - a);                     \
    }                                                       \
  } else pres[i] = px[po[(int)((l-1)*Q)]];                  \
}

// Proper quantile weighting? At least it gives the same results for equal weights of any magnitude.
// See: https://aakinshin.net/posts/weighted-quantiles/
// And: https://en.wikipedia.org/wiki/Percentile#Weighted_percentile
#define WQUANTILE_CORE                                \
double Q, h, a, wb, b;                                \
for(int i = 0, k = 0; i < np; ++i) {                  \
  Q = probs[i];                                       \
  if(Q > 0.0 && Q < 1.0) {                            \
    RETWQSWITCH(sumw, mu);                            \
    while(wsum <= h) wsum += pw[po[k++]];             \
    a = px[po[k == 0 ? 0 : k-1]];                     \
    if(k == 0 || k == l || h == 0.0) {                \
      pres[i] = a; continue;                          \
    }                                                 \
    wb = pw[po[k]];                                   \
    /* If zero weights, need to move forward*/        \
    if(wb == 0.0) {                                   \
    /* separate indices as possible: h < wsum(i-1) */ \
      int kp = k, lm = l-1;                           \
      while(kp < lm && wb == 0.0) wb = pw[po[++kp]];  \
      if(wb == 0.0) {                                 \
        k = kp; pres[i] = a; continue;                \
      }                                               \
      b = px[po[kp]];                                 \
    } else b = px[po[k]];                             \
    h = (wsum - h) / wb;                              \
    pres[i] = b + h * (a - b);                        \
  } else {  /* Since probs must be passed in order*/  \
    if(Q == 0.0) {                                    \
      while(pw[po[k]] == 0.0) ++k;                    \
    } else {                                          \
      k = l-1;                                        \
      while(pw[po[k]] == 0.0) --k;                    \
    }                                                 \
    pres[i] = px[po[k]];                              \
  }                                                   \
}


SEXP fquantileC(SEXP x, SEXP Rprobs, SEXP w, SEXP o, SEXP Rnarm, SEXP Rtype, SEXP Rnames, SEXP checko) {

  if(TYPEOF(Rprobs) != REALSXP) error("probs needs to be a numeric vector");
  int tx = TYPEOF(x), n = length(x), np = length(Rprobs), narm = asLogical(Rnarm), ret = asInteger(Rtype), nprotect = 1;
  if(tx != REALSXP && tx != INTSXP && tx != LGLSXP) error("x needs to be numeric");
  if(ret < 5 || ret > 9) error("fquantile only supports continuous quantile types 5-9. You requested type: %d", ret);

  SEXP res = PROTECT(allocVector(REALSXP, np));
  if(np == 0) { // quantile(x, numeric(0))
    UNPROTECT(nprotect);
    return res;
  }

  double *probs = REAL(Rprobs), *pres = REAL(res);
  unsigned int l = 0;

  for(int i = 0; i < np; ++i) {
    if(probs[i] < 0.0 || probs[i] > 1.0) error("probabilities need to be in in range [0, 1]");
    if(i > 0 && probs[i] < probs[i-1]) error("probabilities need to be passed in ascending order");
  }

  if(asLogical(Rnames)) {
    SEXP names = PROTECT(allocVector(STRSXP, np)); ++nprotect;
    char namei[5];
    for(int i = 0; i < np; ++i) {
      snprintf(namei, 5, "%d%%", (int)(probs[i]*100));
      SET_STRING_ELT(names, i, mkChar(namei));
    }
    namesgets(res, names);
  }

  // First the trivial case
  if(n <= 1) {

    if(!isNull(w)) {
      if(length(w) != n) error("length(w) must match length(x)");
      if(length(w) > 0) {
        double wtmp = asReal(w);
        if(wtmp == 0.0) n = 0;
        else if(ISNAN(wtmp) && NISNAN(asReal(x))) error("Missing weights in order statistics are currently only supported if x is also missing");
      }
    }

    wall0:; // If all weights are zero
    double val = n == 0 ? NA_REAL : tx == REALSXP ? REAL(x)[0] : INTEGER(x)[0] == NA_INTEGER ? NA_REAL : (double)INTEGER(x)[0];
    for(int i = 0; i < np; ++i) pres[i] = val;

  // This case: no quantile estimation, simple range
  } else if(np <= 2 && isNull(o) && (probs[0] == 0.0 || probs[0] == 1.0) && (np <= 1 || probs[1] == 1.0)) {

    // TODO: could also check weights here, but this case is presumably very rare anyway..
    SEXP rng = PROTECT(frange(x, Rnarm)); ++nprotect;
    if(TYPEOF(rng) != REALSXP) {
      rng = PROTECT(coerceVector(rng, REALSXP)); ++nprotect;
    }
    if(probs[0] == 0.0) pres[0] = REAL(rng)[0];
    else if(probs[0] == 1.0) pres[0] = REAL(rng)[1];
    if(np == 2) pres[1] = REAL(rng)[1];

  } else if(isNull(w) && isNull(o)) { // Standard: quickselect

    if(tx == REALSXP) { // Numeric data
      double *x_cc = (double *) R_alloc(n, sizeof(double)), *px = REAL(x), x_min, x_max;
      if(narm) {
        for(unsigned int i = 0; i != n; ++i) if(NISNAN(px[i])) x_cc[l++] = px[i];
        if(l <= 1) { // TODO: More elegant way to solve? Also with integers and weighted estimation ...
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
          for(int i = 0; i < np; ++i) pres[i] = l == 0 ? NA_REAL : (double)x_cc[0];
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
    double *pw = probs, nanw0 = 0.0;

    if(!isNull(o)) {
      if(length(o) != n || TYPEOF(o) != INTSXP) error("o must be a valid ordering vector, of the same length as x and type integer");
      po = INTEGER(o);
      if(asLogical(checko)) { // TODO: Better way?
        for(unsigned int i = 0; i != n; ++i)
          if(po[i] < 1 || po[i] > n) error("Some elements in o are outside of range [1, length(x)]");
      }
    } else  {
      po = (int *) R_alloc(n, sizeof(int)); // Calloc ?
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
          double val = (l == 0 || ISNAN(nanw0)) ? NA_REAL : px[po[0]];
          for(int i = 0; i < np; ++i) pres[i] = val;
          UNPROTECT(nprotect);
          return res;
        }
      } else {
        int *px = INTEGER(x)-1;
        if(px[po[0]] == NA_INTEGER) error("Found missing value at the beginning of the sample. Please use option na.last = TRUE (the default) when creasting ordering vectors for use with fquantile().");
        --po; while(px[po[l]] == NA_INTEGER && l != 0) --l; ++po;
        if(l <= 1) {
          double val = (l == 0 || ISNAN(nanw0)) ? NA_REAL : (double)px[po[0]];
          for(int i = 0; i < np; ++i) pres[i] = val;
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
      double wsum = 0.0, sumw = 0.0; // wsum is running sum, sumw is the total sum

      unsigned int nw0 = 0;
      for(unsigned int i = 0; i != l; ++i) {
        wsum = pw[po[i]];
        if(wsum == 0.0) ++nw0;
        sumw += wsum;
      }
      wsum = 0.0;
      if(ISNAN(sumw)) error("Missing weights in order statistics are currently only supported if x is also missing");
      if(sumw < 0.0) error("Weights must be positive or zero");
      if(l == nw0 || sumw == 0.0) { // error("For weighted quantile estimation, must supply at least one non-zero weight for non-NA x");
        n = 0;
        goto wall0;
      }
      double mu = sumw / (l - nw0);

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



// --------------------------------------------------------------------------
// Then: C rewrite of fnth(), now also supporting (weighted) quantiles
// --------------------------------------------------------------------------

// Without weights, we can apply quickselect at the group-level

double dquickselect(double *x, const int n, const int ret, const double Q) {
  if(n == 0) return NA_REAL;
  unsigned int elem;
  double a, b, h;
  RETQSWITCH(n);
  elem = h; h -= elem; // need to subtract elem
  QUICKSELECT(dswap);
  if((ret < 4 && (ret != 1 || n%2 == 1)) || elem == n-1 || h <= 0.0) return a;
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
  elem = h; h -= elem; // need to subtract elem
  QUICKSELECT(iswap);
  if((ret < 4 && (ret != 1 || n%2 == 1)) || elem == n-1 || h <= 0.0) return (double)a;
  b = x[elem+1];
  for(int i = elem+2; i < n; ++i) if(x[i] < b) b = x[i];
  if(ret == 1 || Q == 0.5) return ((double)a+(double)b)/2.0;
  return (double)a + h*(double)(b-a); // same as (1-h)*(double)a + h*(double)b
}

// With weights, currently radix sort of the entire vector, and then passing through by groups
// TODO: quicksort for weighted statistics at the group-level? could be a partial sort so potentially faster and parallelizable...

// TODO: Check consistency of pointer increments!!

double w_compute_h(const double *pw, const int l, const int ret, const double Q) {
  double sumw = 0.0;
  int nw0 = 0;
  for(int i = 0; i != l; ++i) {
    if(pw[i] == 0.0) ++nw0;
    sumw += pw[i];
  }
  if(ISNAN(sumw)) error("Missing weights in order statistics are currently only supported if x is also missing");
  if(sumw < 0.0) error("Weights must be positive or zero");
  if(l == nw0 || sumw == 0.0) return NA_REAL;
  double mu = sumw / (l - nw0), h;
  RETWQSWITCH(sumw, mu);
  return h;
}

// If no groups or sorted groups pxo is the ordering of x
#define WNTH_CORE                                                          \
double wsum = 0.0, wb, a; /* TODO: reuse wsum, eliminate a */              \
int k = 0;                                                                 \
while(wsum <= h) wsum += pw[pxo[k++]];                                     \
a = px[pxo[k == 0 ? 0 : k-1]];                                             \
if((ret < 4 && (ret != 1 || l%2 == 1)) || k == 0 || k == l || h == 0.0)    \
  return a; /* TODO: ret == 1 averaging */                                 \
wb = pw[pxo[k]];                                                           \
if(wb == 0.0)  /* If zero weights, need to move forward*/                  \
  while(k < l-1 && wb == 0.0) wb = pw[pxo[++k]];                           \
if(wb == 0.0) return a;                                                    \
h = (wsum - h) / wb;                                                       \
wb = px[pxo[k]];                                                           \
return wb + h * (a - wb);


double w_compute_h_grouped(const double *pw, const int *po, const int l, const int ret, const double Q) {
  double sumw = 0.0, mu, h;
  int nw0 = 0;
  for(int i = 0; i != l; ++i) {
    mu = pw[po[i]];
    if(mu == 0.0) ++nw0; // TODO: nw0 += mu == 0.0 ??
    sumw += mu;
  }
  if(ISNAN(sumw)) error("Missing weights in order statistics are currently only supported if x is also missing");
  if(sumw < 0.0) error("Weights must be positive or zero");
  if(l == nw0 || sumw == 0.0) return NA_REAL;
  mu = sumw / (l - nw0);
  RETWQSWITCH(sumw, mu);
  return h;
}

// Otherwise: double indexation, po is the ordering of the groups
#define WNTH_CORE_GROUPED                                                  \
double wsum = 0.0, wb, a; /* TODO: reuse wsum, eliminate a */              \
int k = 0;                                                                 \
while(wsum <= h) wsum += pw[pxo[po[k++]]];                                 \
a = px[pxo[po[k == 0 ? 0 : k-1]]];                                         \
if((ret < 4 && (ret != 1 || l%2 == 1)) || k == 0 || k == l || h == 0.0)    \
  return a; /* TODO: ret == 1 averaging */                                 \
wb = pw[pxo[po[k]]];                                                       \
if(wb == 0.0)  /* If zero weights, need to move forward*/                  \
  while(k < l-1 && wb == 0.0) wb = pw[pxo[po[++k]]];                       \
if(wb == 0.0) return a;                                                    \
h = (wsum - h) / wb;                                                       \
wb = px[pxo[po[k]]];                                                       \
return wb + h * (a - wb);


// Finally, in the default vector method: also provide the option to pass an ordering vector of x, even without weights

#define NTH_ORDVEC                                                          \
RETQSWITCH(l);                                                              \
double a, b;                                                                \
int ih = h; a = px[pxo[ih]];                                                \
if((ret < 4 && (ret != 1 || l%2 == 1)) || ih == l-1 || h <= 0.0) return a;  \
b = px[pxo[ih+1]];                                                          \
return (ret == 1 || Q == 0.5) ? (a+b)/2.0 : a + (h - ih) * (b - a);

#define NTH_ORDVEC_GROUPED                                                  \
RETQSWITCH(l);                                                              \
double a, b;                                                                \
int ih = h; a = px[pxo[po[ih]]];                                            \
if((ret < 4 && (ret != 1 || l%2 == 1)) || ih == l-1 || h <= 0.0) return a;  \
b = px[pxo[po[ih+1]]];                                                      \
return (ret == 1 || Q == 0.5) ? (a+b)/2.0 : a + (h - ih) * (b - a);




// C-implementations for different data types, parallelizable ----------------------------------

double nth_int(const int *px, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l == 0) return NA_REAL;
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

double nth_double(const double *restrict px, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l == 0) return NA_REAL;
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

double w_nth_int(const int *restrict px, const int *restrict pxo, const double *restrict pw, const int *restrict po, double h, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l == 0) return NA_REAL;
  if(l == 1) { // TODO: what about NA_INTEGER and NA/0 weights??
    if(sorted) return ISNAN(pw[0]) ? NA_REAL : (double)px[0];
    return ISNAN(pw[po[0]-1]) ? NA_REAL : (double)px[po[0]];
  }
  if(sorted) { // This refers to the data being sorted by groups (po), not x being sorted (pxo)
    if(narm) { // Adjusting l as necessary... do initial NA check at higher level where pxo is computed...
      while(px[pxo[l]] == NA_INTEGER && l != 0) --l;
      if(l <= 1) return (l == 0 || ISNAN(pw[pxo[, Q)) ? NA_REAL : (double)px[pxo[0]];
    }
    if(h == DBL_MIN) h = w_compute_h(pw+1, l, ret, Q);  // TODO: should only be the case if narm = TRUE, otherwise h should be passed beforehand??
    if(ISNAN(h)) return NA_REAL;
    WNTH_CORE;
  }
  if(narm) { // Adjusting l as necessary... do initial NA check at higher level where pxo is computed...
    while(px[pxo[po[l]]] == NA_INTEGER && l != 0) --l;
    if(l <= 1) return (l == 0 || ISNAN(pw[pxo[po[0]]])) ? NA_REAL : (double)px[pxo[po[0]]];
  }
  if(h = DBL_MIN) h = w_compute_h_grouped(pw, po, l, ret, Q);  // TODO: should only be the case if narm = TRUE, otherwise h should be passed beforehand??
  if(ISNAN(h)) return NA_REAL;
  WNTH_CORE_GROUPED;
}

double w_nth_double(const double *restrict px, const int *restrict pxo, const double *restrict pw, const int *restrict po, double h, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l == 0) return NA_REAL;
  if(l == 1) {
    if(sorted) return ISNAN(pw[0]) ? NA_REAL : px[0];
    return ISNAN(pw[po[0]]) ? NA_REAL : px[po[0]];
  }
  if(sorted) { // This refers to the data being sorted by groups (po), not x being sorted (pxo)
    if(narm) {
      while(ISNAN(px[pxo[l]]) && l != 0) --l;
      if(l <= 1) return (l == 0 || ISNAN(pw[pxo[0]])) ? NA_REAL : px[pxo[0]];
    }
    if(h = DBL_MIN) h = w_compute_h(pw+1, l, ret, Q);  // TODO: should only be the case if narm = TRUE, otherwise h should be passed beforehand??
    if(ISNAN(h)) return NA_REAL;
    WNTH_CORE;
  }
  if(narm) {
    while(ISNAN(px[pxo[po[l]]]) && l != 0) --l;
    if(l <= 1) return (l == 0 || ISNAN(pw[pxo[po[0]]])) ? NA_REAL : px[pxo[po[0]]];
  }
  if(h = DBL_MIN) h = w_compute_h_grouped(pw, po, l, ret, Q);  // TODO: should only be the case if narm = TRUE, otherwise h should be passed beforehand??
  if(ISNAN(h)) return NA_REAL;
  WNTH_CORE_GROUPED;
}

double ord_nth_int(const int *restrict px, const int *restrict pxo, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l == 0) return NA_REAL;
  if(l == 1) return px[sorted ? 0 : po[0]];
  if(sorted) { // This refers to the data being sorted by groups (po), not x being sorted (pxo)
    if(narm) { // Adjusting l as necessary... do initial NA check at higher level where pxo is computed...
      while(px[pxo[l]] == NA_INTEGER && l != 0) --l;
      if(l <= 1) return l == 0 ? NA_REAL : (double)px[pxo[0]];
    }
    NTH_ORDVEC;
  }
  if(narm) { // Adjusting l as necessary... do initial NA check at higher level where pxo is computed...
    while(px[pxo[po[l]]] == NA_INTEGER && l != 0) --l;
    if(l <= 1) return l == 0 ? NA_REAL : (double)px[pxo[po[0]]];
  }
  NTH_ORDVEC_GROUPED;
}

double ord_nth_double(const double *restrict px, const int *restrict pxo, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l == 0) return NA_REAL;
  if(l == 1) return px[sorted ? 0 : po[0]];
  if(sorted) { // This refers to the data being sorted by groups (po), not x being sorted (pxo)
    if(narm) { // Adjusting l as necessary... do initial NA check at higher level where pxo is computed...
      while(ISNAN(px[pxo[l]]) && l != 0) --l;
      if(l <= 1) return l == 0 ? NA_REAL : px[pxo[0]];
    }
    NTH_ORDVEC;
  }
  if(narm) { // Adjusting l as necessary... do initial NA check at higher level where pxo is computed...
    while(ISNAN(px[pxo[po[l]]]) && l != 0) --l;
    if(l <= 1) return l == 0 ? NA_REAL : px[pxo[po[0]]];
  }
  NTH_ORDVEC_GROUPED;
}


// Implementations for R vectors ---------------------------------------------------------------

SEXP nth_impl(SEXP x, int narm, int ret, double Q) {
  int l = length(x);
  if(l <= 1) return x;

  SEXP res;
  switch(TYPEOF(x)) {
    case REALSXP:
      res = ScalarReal(nth_double(REAL(x), &l, l, 1, narm, ret, Q));
      break;
    case INTSXP:
    case LGLSXP:
      res = ScalarReal(nth_int(INTEGER(x), &l, l, 1, narm, ret, Q));
      break;
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }

  if(ATTRIB(x) == R_NilValue || (isObject(x) && inherits(x, "ts"))) return res;
  // PROTECT(res); // Needed ??
  copyMostAttrib(x, res);
  // UNPROTECT(1);
  return res;
}

SEXP w_nth_impl(SEXP x, double *pxo, double *pw, int narm, int ret, double Q) { // , int nthreads
  int l = length(x);
  if(l <= 1) return x;

  // if(nthreads > 1) pxo = (int *)Calloc(l, int);
  // num1radixsort(pxo, TRUE, FALSE, x);

  SEXP res;
  switch(TYPEOF(x)) {
  case REALSXP:
    res = ScalarReal(w_nth_double(REAL(x), pxo, pw, &l, DBL_MIN, l, 1, narm, ret, Q));
    break;
  case INTSXP:
  case LGLSXP:
    res = ScalarReal(w_nth_int(INTEGER(x), pxo, pw, &l, DBL_MIN, l, 1, narm, ret, Q));
    break;
  default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }

  // if(nthreads > 1) Free(pxo);

  if(ATTRIB(x) == R_NilValue || (isObject(x) && inherits(x, "ts"))) return res;
  copyMostAttrib(x, res);
  return res;
}

SEXP ord_nth_impl(SEXP x, double *pxo, int narm, int ret, double Q) {
  int l = length(x);
  if(l <= 1) return x;

  SEXP res;
  switch(TYPEOF(x)) {
  case REALSXP:
    res = ScalarReal(ord_nth_double(REAL(x), pxo, &l, l, 1, narm, ret, Q));
    break;
  case INTSXP:
  case LGLSXP:
    res = ScalarReal(ord_nth_int(INTEGER(x), pxo, &l, l, 1, narm, ret, Q));
    break;
  default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }

  if(ATTRIB(x) == R_NilValue || (isObject(x) && inherits(x, "ts"))) return res;
  copyMostAttrib(x, res);
  return res;
}


SEXP nth_g_impl(SEXP x, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int ret, double Q, int nthreads) {

  int l = length(x);
  if(nthreads > ng) nthreads = ng;

  SEXP res = PROTECT(allocVector(REALSXP, ng));
  double *pres = REAL(res);

  if(sorted) { // Sorted: could compute cumulative group size (= starts) on the fly... but doesn't work multithreaded...
    po = &l;
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = nth_double(px + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = nth_int(px + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  } else { // Not sorted. Perhaps reordering x is faster?
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = nth_double(px, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = nth_int(px, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  }

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")) copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

SEXP w_nth_g_impl(SEXP x, double *pxo, double *pw, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int ret, double Q, int nthreads) {

  int l = length(x);
  if(nthreads > ng) nthreads = ng;

  // how to check if ordering vector is already passed ???
  // num1radixsort(pxo, TRUE, FALSE, x); // sub-column level parallelism. TODO: column-level would be faster with the radix sort... but probably does not make much sense on data frames...

  SEXP res = PROTECT(allocVector(REALSXP, ng));
  double *pres = REAL(res);

  if(sorted) { // Sorted: could compute cumulative group size (= starts) on the fly... but doesn't work multithreaded...
    po = &l;
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = w_nth_double(px + pst[gr]-1, pxo, pw, po, DBL_MIN, pgs[gr], 1, narm, ret, Q); // TODO: where to increment pointers ??
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = w_nth_int(px + pst[gr]-1, pxo, pw, po, DBL_MIN, pgs[gr], 1, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  } else { // Not sorted. Perhaps reordering x is faster?
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = w_nth_double(px, pxo, pw, po + pst[gr]-1, DBL_MIN, pgs[gr], 0, narm, ret, Q); // TODO: where to increment pointers ??
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = w_nth_int(px, pxo, pw, po + pst[gr]-1, DBL_MIN, pgs[gr], 0, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  }

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")) copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

SEXP ord_nth_g_impl(SEXP x, double *pxo, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int ret, double Q, int nthreads) {

  int l = length(x);
  if(nthreads > ng) nthreads = ng;

  SEXP res = PROTECT(allocVector(REALSXP, ng));
  double *pres = REAL(res);

  if(sorted) { // Sorted: could compute cumulative group size (= starts) on the fly... but doesn't work multithreaded...
    po = &l;
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = ord_nth_double(px + pst[gr]-1, pxo, po, pgs[gr], 1, narm, ret, Q);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = ord_nth_int(px + pst[gr]-1, pxo, po, pgs[gr], 1, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  } else { // Not sorted. Perhaps reordering x is faster?
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = ord_nth_double(px, pxo, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q); // TODO: where to increment pointers ??
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = ord_nth_int(px, pxo, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  }

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")) copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}


// // This only works for a single vector, po is the ordering of x
// SEXP w_ord_nth_impl(SEXP x, int *po, double *pw, double h, int narm, int ret, double Q) {
//   int n = length(x), l = n;
//
//   if(narm) {
//     if(tx == REALSXP) { // Numeric data
//       double *px = REAL(x)-1;
//       if(ISNAN(px[po[0]])) error("Found missing value at the beginning of the sample. Please use option na.last = TRUE (the default) when creasting ordering vectors for use with fnth().");
//       --po; while(ISNAN(px[po[l]]) && l != 0) --l; ++po;
//       if(l <= 1) return (l == 0 || ISNAN(pw[po[0]])) ? NA_REAL : px[po[0]];
//     } else {
//       int *px = INTEGER(x)-1;
//       if(px[po[0]] == NA_INTEGER) error("Found missing value at the beginning of the sample. Please use option na.last = TRUE (the default) when creasting ordering vectors for use with fnth().");
//       --po; while(px[po[l]] == NA_INTEGER && l != 0) --l; ++po;
//       if(l <= 1) return (l == 0 || ISNAN(pw[po[0]])) ? NA_REAL : (double)px[po[0]];
//     }
//   }
//
//   double res, a, b; // result and required constants
//
//   if(ISNAN(h)) { // If just an ordering vector is supplied
//
//   #define NTH_ORDVEC                                                          \
//     RETQSWITCH(l);                                                            \
//     int ih = h; a = px[po[ih]];                                               \
//     if((ret < 4 && (ret != 1 || l%2 == 1)) || ih == l-1 || h <= 0.0) res = a; \
//     else {                                                                    \
//       b = px[po[ih+1]];                                                       \
//       res == (ret == 1 || Q == 0.5) ? (a+b)/2.0 : a + (h - ih) * (b - a);     \
//     }
//
//     if(tx == REALSXP) { // Numeric data
//       double *px = REAL(x)-1;
//       NTH_ORDVEC;
//     } else {
//       int *px = INTEGER(x)-1;
//       NTH_ORDVEC;
//     }
//
//   } else { // Weighted quantile
//
//     double wsum = 0.0, wb; // wsum is running sum
//
//     if(narm || h = DBL_MIN) {
//
//       double sumw = 0.0;
//       int nw0 = 0;
//
//       for(unsigned int i = 0; i != l; ++i) {
//         wsum = pw[po[i]];
//         if(wsum == 0.0) ++nw0;
//         sumw += wsum;
//       }
//       wsum = 0.0;
//       if(ISNAN(sumw)) error("Missing weights in order statistics are currently only supported if x is also missing");
//       else if(sumw < 0.0) error("Weights must be positive or zero");
//       if(l == nw0 || sumw == 0.0) l = 0; /// TODO: what to do here ??
//       else mu = sumw / (l - nw0); // repetitive? provide as input??
//
//       RETWQSWITCH(sumw, mu);
//     }
//
//      #define WNTH_CORE                                                             \
//         while(wsum <= h) wsum += pw[po[k++]];                                      \
//         a = px[po[k == 0 ? 0 : k-1]];                                              \
//         if((ret < 4 && (ret != 1 || l%2 == 1)) || k == 0 || k == l || h == 0.0) {  \
//           res = a; /* TODO: ret == 1 averaging */                                  \
//         } else {                                                                   \
//           wb = pw[po[k]];                                                          \
//           /* If zero weights, need to move forward*/                               \
//           if(wb == 0.0)                                                            \
//              while(k < l-1 && wb == 0.0) wb = pw[po[++k]];                         \
//           if(wb == 0.0) pres = a;                                                  \
//           else {                                                                   \
//             b = px[po[k]];                                                         \
//             h = (wsum - h) / wb;                                                   \
//             res = b + h * (a - b);                                                 \
//           }                                                                        \
//         }
//       // Previous zero-weight code in fnth: forward iteration with simple averaging of adjacent order statistics.
//       // if(wsum == h) {
//       //   double out = px[po[k-1]], n = 2;
//       //   while(pw[po[k]] == 0) {
//       //     out += px[po[k++]];
//       //     ++n;
//       //   }
//       //   pres[i] = (out + px[po[k]]) / n;
//       // }
//
//     if(tx == REALSXP) { // Numeric data
//       double *px = REAL(x)-1;
//       WNTH_CORE;
//     } else {
//       int *px = INTEGER(x)-1;
//       WNTH_CORE;
//     }
//
//   }
//
//   if(ATTRIB(x) == R_NilValue || (isObject(x) && inherits(x, "ts"))) return ScalarReal(res);
//
//   SEXP out = ScalarReal(res); // Needs protection ??
//   copyMostAttrib(x, out);
//   return out;
// }
//
// // NOTE: This is too complicated: need to create lots of intermediate vectors
// SEXP w_ord_nth_g_impl(SEXP x, int *po, double *pw, int *pg, int ng, int *pgs, int sorted, double h, int narm, int ret, double Q) {
//
//   int l = length(x), tx = TYPEOF(x);
//
//   if(sorted) { // TODO: things are easy.. just pass through o
//
//   }
//
//   double *wsumQ = (double*)Calloc(ng, double)-1, *wsum = (double*)Calloc(ng, double);
//   SEXP out = PROTECT(allocVector(REALSXP, ng));
//
//   if(narm) {
//     switch(tx) {
//     case REALSXP:
//       for(int i = 0; i != l; ++i) if(NISNAN(x[i])) wsumQ[pg[i]] += pw[i];
//       break;
//     case INTSXP:
//     case LGLSXP:
//       for(int i = 0; i != l; ++i) if(x[i] != NA_INTEGER) wsumQ[pg[i]] += pw[i];
//       break;
//     default: error("Not Supported SEXP Type: '%s'", type2char(tx));
//     }
//   } else {
//     for(int i = 0; i != l; ++i) wsumQ[pg[i]] += pw[i];
//   }
//
//   for(int i = 0; i != ng; ++i) {
//     if(ISNAN(wsumQ[i])) stop("Missing weights in order statistics are currently only supported if x is also missing");
//     wsumQ[i] *= Q;
//   }
//
//
//   // wsumQ = wsumQ * Q;
//   int gi, oi;
//   if(tiesmean) {
//     std::vector<bool> seen(ng);
//     std::vector<int> n(ng, 1); // only needed for 0 weights. Could check above if there are any.
//     for(int i = 0; i != l; ++i) {
//       oi = o[i]-1;
//       gi = g[oi]-1;
//       if(seen[gi]) continue;
//       if(wsum[gi] < wsumQ[gi]) out[gi] = x[oi];
//       else {
//         if(wsum[gi] > wsumQ[gi]) {
//           seen[gi] = true;
//           continue;
//         }
//         out[gi] += (x[oi]-out[gi])/++n[gi]; // https://stackoverflow.com/questions/28820904/how-to-efficiently-compute-average-on-the-fly-moving-average
//       }
//       wsum[gi] += wg[oi];
//     }
//   } else if(lower) {
//     for(int i = 0; i != l; ++i) {
//       oi = o[i]-1;
//       gi = g[oi]-1;
//       if(wsum[gi] < wsumQ[gi]) {
//         wsum[gi] += wg[oi];
//         out[gi] = x[oi];
//       }
//     }
//   } else {
//     for(int i = 0; i != l; ++i) {
//       oi = o[i]-1;
//       gi = g[oi]-1;
//       if(wsum[gi] <= wsumQ[gi]) {
//         wsum[gi] += wg[oi];
//         out[gi] = x[oi];
//       }
//     }
//   }
//
//
//
// #define NTH_ORDVEC                                                         \
//   RETQSWITCH(l);                                                            \
//   int ih = h; a = px[po[ih]];                                               \
//   if((ret < 4 && (ret != 1 || l%2 == 1)) || ih == l-1 || h <= 0.0) res = a; \
//   else {                                                                    \
//     b = px[po[ih+1]];                                                       \
//     res == (ret == 1 || Q == 0.5) ? (a+b)/2.0 : a + (h - ih) * (b - a);     \
//   }
//
// if(tx == REALSXP) { // Numeric data
//   double *px = REAL(x)-1;
//   NTH_ORDVEC;
// } else {
//   int *px = INTEGER(x)-1;
//   NTH_ORDVEC;
// }
//
//
// }




// Functions for Export --------------------------------------------------------

// TODO: Single thread optimization: re-use array for quickselect?? should be quite easy... just check if nthreads = 1,
// otherwise assign pointer...
// Also for multiple columns with weights if na.rm = FALSE, can compute h's and supply repeatedly for each column

// Function for atomic vectors: has extra arguments o and checko for passing external ordering vector
// This is mean to speed up computation of several (grouped) quantiles on the same data
SEXP fnthC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rret, SEXP Rnthreads, SEXP o, SEXP checko) {

  int nullg = isNull(g), nullw = isNull(w), nullo = isNull(o), l = length(x), narm = asLogical(Rnarm),
      ret = asInteger(Rret), nprotect = 0;

  if(l <= 1) return x;
  if(length(p) != 1) error("fnth supports only a single element / quantile. Use fquantile for multiple quantiles.");

  double Q = asReal(p);
  if(ISNAN(Q) || Q <= 0.0 || Q == 1.0) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
  if(Q > 1.0) {
    ret = 2; // Correct ??
    if(nullg) {
      if(Q >= l) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
      Q = (Q-1.0)/(l-1);
    } else {
      if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
      int ng = INTEGER(VECTOR_ELT(g, 0))[0];
      if(Q >= (double)l/ng) error("n needs to be between 0 and 1, or between 1 and the length(x)/ng, with ng the number of groups. Use fmin and fmax for minima and maxima.");
      Q = (Q-1.0)/((double)l/ng-1.0);
    }
  }

  // First the simplest case
  if(nullg && nullw && nullo) return nth_impl(x, narm, ret, Q);

  // Creating pointers that may or may not be needed
  double *pw = &Q;
  int *pxo = &l;

  // Preprocessing o
  if(!nullo) {
    if(length(o) != l || TYPEOF(o) != INTSXP) error("o must be a valid ordering vector, of the same length as x and type integer");
    pxo = INTEGER(o);
    if(asLogical(checko)) { // TODO: Better way?
      for(unsigned int i = 0; i != l; ++i)
          if(pxo[i] < 1 || pxo[i] > l) error("Some elements in o are outside of range [1, length(x)]");
    }
  }

  // Preprocessing w, computing ordering of x if not supplied
  if(!nullw) {
    if(length(w) != l) error("length(w) must match length(x)");
    if(TYPEOF(w) != REALSXP) {
      if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double). You supplied a vector of type: '%s'", type2char(TYPEOF(w)));
      w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
    }
    pw = REAL(w);
    if(nullo) {
      nullo = 0;
      pxo = (int *) R_alloc(l, sizeof(int));
      num1radixsort(pxo, TRUE, FALSE, x);
    }
  }

  SEXP res; // result

  // If no groups, return using suitable functions
  if(nullg) {
    if(nullw) res = ord_nth_impl(x, pxo, narm, ret, Q);
    else res = w_nth_impl(x, pxo, pw, narm, ret, Q);
    UNPROTECT(nprotect);
    return res;
  }

  // Preprocessing g, computing ordering vector if not supplied
  if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
  const SEXP *restrict pg = SEXPPTR(g), ord = pg[6];
  int sorted = LOGICAL(pg[5])[1] == 1, ng = INTEGER(pg[0])[0], *restrict pgs = INTEGER(pg[2]), *restrict pord, *restrict pst;
  if(l != length(pg[1])) error("length(g) must match length(x)");
  if(isNull(ord)) {
    int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;
    for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i]; // TODO: get maxgrpn?
    pst = cgs + 1;
    if(sorted) pord = &l;
    else {
      int *restrict count = (int *) Calloc(ng+1, int);
      pord = (int *) R_alloc(l, sizeof(int)); --pord;
      for(int i = 0; i != l; ++i) pord[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
      ++pord; Free(count);
    }
  } else {
    pord = INTEGER(ord);
    pst = INTEGER(getAttrib(ord, install("starts")));
  }

  if(nullw && nullo) res = nth_g_impl(x, ng, pgs, po, pst, sorted, narm, ret, Q, asInteger(Rnthreads));
  else if(nullw) res = ord_nth_g_impl(x, pxo, ng, pgs, po, pst, sorted, narm, ret, Q, asInteger(Rnthreads));
  else res = w_nth_g_impl(x, pxo, pw, ng, pgs, po, pst, sorted, narm, ret, Q, asInteger(Rnthreads));
  UNPROTECT(nprotect);
  return res;
}

// // TODO: allow column-level parallelism??
// SEXP fnthlC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rret, SEXP Rnthreads) {
//   int nullg = isNull(g), nullw = isNull(w), l = length(x), ng = 0, nprotect = 1,
//     narm = asLogical(Rnarm), drop = asLogical(Rdrop), ret = asInteger(Rret), nthreads = asInteger(Rnthreads);
//   if(l < 1) return x;
//
//   double Q = asReal(p);
//   if(ISNAN(Q) || Q <= 0 || Q == 1) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
//   if(Q > 1) {
//     ret = 2; // Correct ??
//     if(nullg) {
//       if(Q >= l) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
//       Q = (Q-1)/(l-1);
//     } else {
//       if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
//       int ng = length(VECTOR_ELT(g, 2));
//       if(Q >= l/ng) error("n needs to be between 0 and 1, or between 1 and the length(x)/ng, with ng the number of groups. Use fmin and fmax for minima and maxima.");
//       Q = (Q-1)/(l/ng-1);
//     }
//   }
//
//   SEXP out = PROTECT(allocVector(nullg && drop ? REALSXP : VECSXP, l)), *restrict px = SEXPPTR(x);
//   int nrx = length(px[0]);
//   double tmp = 0.0, *restrict pw = &tmp;
//
//   if(!nullw) {
//     if(length(w) != nrx) error("length(w) must match nrow(x)");
//     if(TYPEOF(w) != REALSXP) {
//       if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double)");
//       SEXP wd = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
//       pw = REAL(wd);
//     } else pw = REAL(w);
//   }
//
//   if(nullg) {
//     if(nthreads > l) nthreads = l;
//     if(drop) {
//       double *restrict pout = REAL(out);
//       if(nullw) {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < l; ++j) pout[j] = REAL(nth_impl(px[j], narm, ret, Q))[0];
//       } else {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < l; ++j) pout[j] = REAL(w_nth_impl(px[j], pw, narm, ret, Q))[0];
//       }
//       setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
//       UNPROTECT(nprotect);
//       return out;
//     }
//     SEXP *restrict pout = SEXPPTR(out);
//     if(nullw) {
//       #pragma omp parallel for num_threads(nthreads)
//       for(int j = 0; j < l; ++j) pout[j] = nth_impl(px[j], narm, ret, Q);
//     } else {
//       #pragma omp parallel for num_threads(nthreads)
//       for(int j = 0; j < l; ++j) pout[j] = w_nth_impl(px[j], pw, narm, ret, Q);
//     }
//   } else {
//     if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
//     const SEXP *restrict pg = SEXPPTR(g), o = pg[6];
//     ng = INTEGER(pg[0])[0];
//     int sorted = LOGICAL(pg[5])[1] == 1, *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst;
//     if(nrx != length(pg[1])) error("length(g) must match nrow(x)");
//     if(isNull(o)) {
//       int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;
//       for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
//       pst = cgs + 1;
//       if(sorted) po = &l;
//       else {
//         int *restrict count = (int *) Calloc(ng+1, int);
//         po = (int *) R_alloc(nrx, sizeof(int)); --po;
//         for(int i = 0; i != nrx; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
//         ++po; Free(count);
//       }
//     } else {
//       po = INTEGER(o);
//       pst = INTEGER(getAttrib(o, install("starts")));
//     }
//     SEXP *restrict pout = SEXPPTR(out);
//     if(nullw) { // Parallelism at sub-column level
//       for(int j = 0; j < l; ++j) pout[j] = nth_g_impl(px[j], ng, pgs, po, pst, sorted, narm, ret, Q, nthreads);
//     } else { // Parallelism at sub-column level
//       for(int j = 0; j < l; ++j) pout[j] = w_nth_g_impl(px[j], pw, ng, pgs, po, pst, sorted, narm, ret, Q, nthreads);
//     }
//   }
//
//   DFcopyAttr(out, x, ng);
//   UNPROTECT(nprotect);
//   return out;
// }
//
// SEXP fnthmC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rret, SEXP Rnthreads) {
//   SEXP dim = getAttrib(x, R_DimSymbol);
//   if(isNull(dim)) error("x is not a matrix");
//   int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1],
//       narm = asLogical(Rnarm), ret = asInteger(Rret), nthreads = asInteger(Rnthreads),
//       nullg = isNull(g), nullw = isNull(w), nprotect = 1;
//   if(l <= 1) return x; // Prevents seqfault for numeric(0) #101
//   if(nthreads > col) nthreads = col;
//   if(length(p) != 1) error("fnth supports only a single element / quantile. Use fquantile for multiple quantiles.");
//   double Q = asReal(p);
//   if(ISNAN(Q) || Q <= 0 || Q == 1) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
//   if(Q > 1) {
//     ret = 2; // Correct ??
//     if(nullg) {
//       if(Q >= l) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
//       Q = (Q-1)/(l-1);
//     } else {
//       if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
//       int ng = length(VECTOR_ELT(g, 2));
//       if(Q >= l/ng) error("n needs to be between 0 and 1, or between 1 and the length(x)/ng, with ng the number of groups. Use fmin and fmax for minima and maxima.");
//       Q = (Q-1)/(l/ng-1);
//     }
//   }
//
//   double tmp = 0.0, *restrict pw = &tmp;
//   if(!nullw) {
//     if(length(w) != l) error("length(w) must match nrow(x)");
//     if(TYPEOF(w) != REALSXP) {
//       if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double)");
//       SEXP wd = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
//       pw = REAL(wd);
//     } else pw = REAL(w);
//   }
//
//   if(nullg) {
//     SEXP res = PROTECT(allocVector(tx, col));
//
//     switch(tx) {
//     case REALSXP: {
//       double *px = REAL(x), *restrict pres = REAL(res);
//       if(nullw) {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) pres[j] = nth_double(px + j*l, &l, l, 1, narm, ret, Q);
//       } else {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) pres[j] = w_nth_double(px + j*l, pw, &l, l, 1, narm, ret, Q);
//       }
//       break;
//     }
//     case INTSXP:
//     case LGLSXP: {  // Factor matrix not well defined object...
//       int *px = INTEGER(x), *restrict pres = INTEGER(res);
//       if(nullw) {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) pres[j] = nth_int(px + j*l, &l, l, 1, narm, ret, Q);
//       } else {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) pres[j] = w_nth_int(px + j*l, pw, &l, l, 1, narm, ret, Q);
//       }
//       break;
//     }
//     default: error("Not Supported SEXP Type!");
//     }
//
//     matCopyAttr(res, x, Rdrop, /*ng=*/0);
//     UNPROTECT(nprotect);
//     return res;
//   }
//
//   // With groups
//   if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
//   const SEXP *restrict pg = SEXPPTR(g), o = pg[6];
//   int sorted = LOGICAL(pg[5])[1] == 1, ng = INTEGER(pg[0])[0], *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst, gl = length(pg[1]);
//   if(l != gl) error("length(g) must match nrow(x)");
//   SEXP res = PROTECT(allocVector(tx, ng * col));
//
//   if(isNull(o)) {
//     int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;
//     for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
//     pst = cgs + 1;
//     if(sorted) po = &l;
//     else {
//       int *restrict count = (int *) Calloc(ng+1, int);
//       po = (int *) R_alloc(l, sizeof(int)); --po;
//       for(int i = 0; i != l; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
//       ++po; Free(count);
//     }
//   } else {
//     po = INTEGER(o);
//     pst = INTEGER(getAttrib(o, install("starts")));
//   }
//
//   if(sorted) { // Sorted
//     switch(tx) {
//     case REALSXP: {
//       double *px = REAL(x), *restrict pres = REAL(res);
//       if(nullw) {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) {
//           int jng = j * ng;
//           double *pxj = px + j * l;
//           for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : nth_double(pxj + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
//         }
//       } else {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) {
//           int jng = j * ng;
//           double *pxj = px + j * l;
//           for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : w_nth_double(pxj + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
//         }
//       }
//       break;
//     }
//     case INTSXP:
//     case LGLSXP: { // Factor matrix not well defined object...
//       int *px = INTEGER(x), *restrict pres = INTEGER(res);
//       if(nullw) {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) {
//           int *pxj = px + j * l, jng = j * ng;
//           for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : nth_int(pxj + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
//         }
//       } else {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) {
//           int *pxj = px + j * l, jng = j * ng;
//           for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : w_nth_int(pxj + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret, Q);
//         }
//       }
//       break;
//     }
//     default: error("Not Supported SEXP Type!");
//     }
//   } else { // Not sorted
//     switch(tx) {
//     case REALSXP: {
//       double *px = REAL(x), *restrict pres = REAL(res);
//       if(nullw) {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) {
//           int jng = j * ng;
//           double *pxj = px + j * l;
//           for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : nth_double(pxj, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
//         }
//       } else {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) {
//           int jng = j * ng;
//           double *pxj = px + j * l;
//           for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : w_nth_double(pxj, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
//         }
//       }
//       break;
//     }
//     case INTSXP:
//     case LGLSXP: {
//       int *px = INTEGER(x), *restrict pres = INTEGER(res);
//       if(nullw) {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) {
//           int jng = j * ng, *pxj = px + j * l;
//           for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : nth_int(pxj, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
//         }
//       } else {
//         #pragma omp parallel for num_threads(nthreads)
//         for(int j = 0; j < col; ++j) {
//           int jng = j * ng, *pxj = px + j * l;
//           for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : w_nth_int(pxj, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret, Q);
//         }
//       }
//       break;
//     }
//     default: error("Not Supported SEXP Type!");
//     }
//   }
//
//   matCopyAttr(res, x, Rdrop, ng);
//   UNPROTECT(nprotect);
//   return res;
// }
//
