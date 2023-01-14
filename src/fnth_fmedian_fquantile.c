#ifdef _OPENMP
#include <omp.h>
#endif
#include "collapse_c.h"

/*
 Inspired by Numerical Recipes in C and data.table's quickselect.c,
 R's quantile() function, Rfast2::Quantile(), and these references for sample quantiles:
   https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample
   https://doi.org/10.2307/2684934
   https://aakinshin.net/posts/weighted-quantiles/
 At large, the weighted quantile algorithm is my own cooking [(C) 2022 Sebastian Krantz]
*/

// Adopted from data.table's quickselect.c
static inline void iswap(int *a, int *b)           {int     tmp=*a; *a=*b; *b=tmp;}
static inline void dswap(double *a, double *b)     {double  tmp=*a; *a=*b; *b=tmp;}

// Barebones quickselect algorithm from Numerical Recipes in C
#undef QUICKSELECT
#define QUICKSELECT(SWAP)                                                         \
  unsigned int ir = n-1, l = 0, lp;                                               \
  for(;;) {                                                                       \
    lp = l+1;                                                                     \
    if (ir <= lp) {  /* Active partition contains 1 or 2 elements. */             \
      if (ir == lp && x[ir] < x[l]) { /* Case of 2 elements. */                   \
        SWAP(x+l, x+ir);                                                          \
      }                                                                           \
      break;                                                                      \
    } else {                                                                      \
      unsigned int mid=(l+ir) >> 1;  /* Choose median of left, center, and right elements as partitioning element a. */  \
      SWAP(x+mid, x+lp); /* Also rearrange so that arr[l] ≤ arr[l+1] ≤ arr[ir] */ \
      if (x[l] > x[ir]) {                                                         \
        SWAP(x+l, x+ir);                                                          \
      }                                                                           \
      if (x[lp] > x[ir]) {                                                        \
        SWAP(x+lp, x+ir);                                                         \
      }                                                                           \
      if (x[l] > x[lp]) {                                                         \
        SWAP(x+l, x+lp);                                                          \
      }                                                                           \
      unsigned int i=lp, j=ir; /* Initialize pointers for partitioning. */        \
      a=x[lp];    /* Partitioning element. */                                     \
      for (;;) {  /* Beginning of innermost loop. */                              \
        do i++; while (x[i] < a);  /* Scan up to find element > a. */             \
        do j--; while (x[j] > a);  /* Scan down to find element < a. */           \
        if (j < i) break;  /* Pointers crossed. Partitioning complete. */         \
        SWAP(x+i, x+j);                                                           \
      }            /* End of innermost loop. */                                   \
      x[lp]=x[j];  /* Insert partitioning element. */                             \
      x[j]=a;                                                                     \
      if (j >= elem) ir=j-1; /* if index of partitioning element j is above median index */ \
      if (j <= elem) l=i;    /* if index of partitioning element j is below median index */ \
    }                                                                             \
  }                                                                               \
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
case 7: /* quantile type 7 */                                  \
  h = (sumw - mu) * Q;                                         \
  break;                                                       \
case 1:                                                        \
case 2:                                                        \
case 3: /* average, lower or upper element (adjust algorithm)*/\
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

#undef FQUANTILE_CORE
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
// Expects px to be decremented by 1
#undef FQUANTILE_ORDVEC
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
// Expects px and pw to be decremented by 1
#undef WQUANTILE_CORE
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
        --po; while(l != 0 && ISNAN(px[po[l]])) --l; ++po;
        if(l <= 1) {
          double val = (l == 0 || ISNAN(nanw0)) ? NA_REAL : px[po[0]];
          for(int i = 0; i < np; ++i) pres[i] = val;
          UNPROTECT(nprotect);
          return res;
        }
      } else {
        int *px = INTEGER(x)-1;
        if(px[po[0]] == NA_INTEGER) error("Found missing value at the beginning of the sample. Please use option na.last = TRUE (the default) when creasting ordering vectors for use with fquantile().");
        --po; while(l != 0 && px[po[l]] == NA_INTEGER) --l; ++po;
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
  if(ret == 1) return (a+b)/2.0; //  || Q == 0.5
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
  if(ret == 1) return ((double)a+(double)b)/2.0; //  || Q == 0.5
  return (double)a + h*(double)(b-a); // same as (1-h)*(double)a + h*(double)b
}

// With weights, either radix sort of the entire vector, and then passing through by groups,
// or quicksort at the group-level

// Expects pw and po to be consistent
double w_compute_h(const double *pw, const int *po, const int l, const int sorted, const int ret, const double Q) {
  double sumw = 0.0, mu, h;
  int nw0 = 0;
  if(sorted) {
    for(int i = 0; i != l; ++i) {
      if(pw[i] == 0.0) ++nw0; // nw0 += pw[i] == 0.0 -> seems not faster...
      sumw += pw[i];
    }
  } else {
    for(int i = 0; i != l; ++i) {
      mu = pw[po[i]];
      if(mu == 0.0) ++nw0; // nw0 += mu == 0.0 -> seems not faster...
      sumw += mu;
    }
  }
  if(ISNAN(sumw)) error("Missing weights in order statistics are currently only supported if x is also missing");
  if(sumw < 0.0) error("Weights must be positive or zero");
  if(l == nw0 || sumw == 0.0) return NA_REAL;
  mu = sumw / (l - nw0);
  RETWQSWITCH(sumw, mu);
  return h;
}

// If no groups or sorted groups po is the ordering of x
// Expects pointers px and pw to be decremented by one
#undef WNTH_CORE
#define WNTH_CORE                                                          \
double wsum = pw[po[0]], wb, a;                                            \
int k = 1;                                                                 \
if(ret < 3) { /* lower (2), or average (1) element*/                       \
  while(wsum < h) wsum += pw[po[k++]];                                     \
  a = px[po[k-1]];                                                         \
  if(ret == 2 || wsum != h) return a; /* h = sumw * Q must be > 0 here */  \
  wsum = 2.0; wb = px[po[k]];                                              \
  while(pw[po[k]] == 0.0) { /* l should never be reached, I tested it */   \
    wb += px[po[++k]]; ++wsum;                                             \
  }                                                                        \
  return (a + wb) / wsum;                                                  \
}                                                                          \
while(wsum <= h) wsum += pw[po[k++]];                                      \
a = px[po[k-1]];                                                           \
if(ret == 3 || k == l || h == 0.0)                                         \
  return a;                                                                \
wb = pw[po[k]];                                                            \
if(wb == 0.0) { /* If zero weights, need to move forward*/                 \
  while(k < l-1 && wb == 0.0) wb = pw[po[++k]];                            \
  if(wb == 0.0) return a;                                                  \
}                                                                          \
h = (wsum - h) / wb;                                                       \
wb = px[po[k]];                                                            \
return wb + h * (a - wb);

// This is the same, just that the result is assigned. Needed for quicksort based implementations
// Does not require incremented pointers (depending on the content of i_cc)
#undef WNTH_CORE_QSORT
#define WNTH_CORE_QSORT                                                      \
double res, wsum = pw[i_cc[0]], wb, a;                                       \
int k = 1;                                                                   \
if(ret < 3) { /* lower (2), or average (1) element*/                         \
  while(wsum < h) wsum += pw[i_cc[k++]];                                     \
  a = x_cc[k-1];                                                             \
  if(ret == 2 || wsum != h) res = a; /* h = sumw * Q must be > 0 here */     \
  else {                                                                     \
    wsum = 2.0; wb = x_cc[k];                                                \
    while(pw[i_cc[k]] == 0.0) { /* n should never be reached, I tested it */ \
      wb += x_cc[++k]; ++wsum;                                               \
    }                                                                        \
    res = (a + wb) / wsum;                                                   \
  }                                                                          \
} else {                                                                     \
  while(wsum <= h) wsum += pw[i_cc[k++]];                                    \
  a = x_cc[k-1];                                                             \
  if(ret == 3 || k == n || h == 0.0) {                                       \
    res = a;                                                                 \
  } else {                                                                   \
    wb = pw[i_cc[k]];                                                        \
    if(wb == 0.0)  /* If zero weights, need to move forward*/                \
       while(k < n-1 && wb == 0.0) wb = pw[i_cc[++k]];                       \
    if(wb == 0.0) res = a;                                                   \
    else {                                                                   \
      h = (wsum - h) / wb;                                                   \
      wb = x_cc[k];                                                          \
      res = wb + h * (a - wb);                                               \
    }                                                                        \
  }                                                                          \
}

// Finally, in the default vector method: also provide the option to pass an ordering vector of x, even without weights
// if the groups are unsorted, po needs to be recomputed to provide the ordering within groups

// Expects pointer px to be decremented by 1
#undef NTH_ORDVEC
#define NTH_ORDVEC                                                         \
double a, b, h;                                                            \
RETQSWITCH(l);                                                             \
int ih = h; a = px[po[ih]]; h -= ih;                                       \
if((ret < 4 && (ret != 1 || l%2 == 1)) || ih == l-1 || h <= 0.0) return a; \
b = px[po[ih+1]];                                                          \
return (ret == 1) ? (a+b)/2.0 : a + h * (b - a); //  || Q == 0.5



// C-implementations for different data types, parallelizable ----------------------------------

double nth_int(const int *restrict px, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l <= 1) return l == 0 ? NA_REAL : sorted ? (double)px[0] : (double)px[po[0]-1];

  int *x_cc = (int *) Calloc(l, int), n = 0;
  if(sorted) {
    // if(narm) {
      for(int i = 0; i != l; ++i) if(px[i] != NA_INTEGER) x_cc[n++] = px[i];
    // } else {
    //   n = l;
    //   memcpy(x_cc, px, l * sizeof(int));
    // }
  } else {
    const int *pxm = px-1; // creating offset pointer to x
    // if(narm) {
      for(int i = 0; i != l; ++i) if(pxm[po[i]] != NA_INTEGER) x_cc[n++] = pxm[po[i]];
    // } else {
    //   n = l;
    //   for(int i = 0; i != l; ++i) x_cc[i] = pxm[po[i]];
    // }
  }

  double res = (narm == 0 && n != l) ? NA_REAL : iquickselect(x_cc, n, ret, Q);
  Free(x_cc);
  return res;
}

double nth_int_noalloc(const int *restrict px, const int *restrict po, int *x_cc, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l <= 1) return l == 0 ? NA_REAL : sorted ? (double)px[0] : (double)px[po[0]-1];

  int n = 0;

  if(sorted) {
    for(int i = 0; i != l; ++i) if(px[i] != NA_INTEGER) x_cc[n++] = px[i];
  } else {
    const int *pxm = px-1; // creating offset pointer to x
    for(int i = 0; i != l; ++i) if(pxm[po[i]] != NA_INTEGER) x_cc[n++] = pxm[po[i]];
  }

  return (narm == 0 && n != l) ? NA_REAL : iquickselect(x_cc, n, ret, Q);
}

double nth_double(const double *restrict px, const int *restrict po, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l <= 1) return l == 0 ? NA_REAL : sorted ? px[0] : px[po[0]-1];

  double *x_cc = (double *) Calloc(l, double);
  int n = 0;

  if(sorted) {
    // if(narm) {
      for(int i = 0; i != l; ++i) if(NISNAN(px[i])) x_cc[n++] = px[i];
    // } else {
    //   n = l;
    //   memcpy(x_cc, px, l * sizeof(double));
    // }
  } else {
    const double *pxm = px-1;
    // if(narm) {
      for(int i = 0; i != l; ++i) if(NISNAN(pxm[po[i]])) x_cc[n++] = pxm[po[i]];
    // } else {
    //   n = l;
    //   for(int i = 0; i != l; ++i) x_cc[i] = pxm[po[i]];
    // }
  }

  double res = (narm == 0 && n != l) ? NA_REAL : dquickselect(x_cc, n, ret, Q);
  Free(x_cc);
  return res;
}

double nth_double_noalloc(const double *restrict px, const int *restrict po, double *x_cc, const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l <= 1) return l == 0 ? NA_REAL : sorted ? px[0] : px[po[0]-1];

  int n = 0;

  if(sorted) {
    for(int i = 0; i != l; ++i) if(NISNAN(px[i])) x_cc[n++] = px[i];
  } else {
    const double *pxm = px-1;
    for(int i = 0; i != l; ++i) if(NISNAN(pxm[po[i]])) x_cc[n++] = pxm[po[i]];
  }

  return (narm == 0 && n != l) ? NA_REAL : dquickselect(x_cc, n, ret, Q);
}

// Expects pointer px to be decremented by 1
double nth_int_ord(const int *restrict px, const int *restrict po, int l, const int narm, const int ret, const double Q) {
  if(l <= 1) return l == 0 ? NA_REAL : (double)px[po[0]];
  if(narm) { // Adjusting l as necessary... initial NA check done in fnthC()
    while(l != 0 && px[po[l-1]] == NA_INTEGER) --l;
    if(l <= 1) return l == 0 ? NA_REAL : (double)px[po[0]];
  } else if(px[po[l-1]] == NA_INTEGER) return NA_REAL;
  NTH_ORDVEC;
}

// Expects pointer px to be decremented by 1
double nth_double_ord(const double *restrict px, const int *restrict po, int l, const int narm, const int ret, const double Q) {
  if(l <= 1) return l == 0 ? NA_REAL : px[po[0]];
  if(narm) { // Adjusting l as necessary... initial NA check done in fnthC()
    while(l != 0 && ISNAN(px[po[l-1]])) --l;
    if(l <= 1) return l == 0 ? NA_REAL : px[po[0]];
  } else if(ISNAN(px[po[l-1]])) return NA_REAL;
  NTH_ORDVEC;
}

// Expects pointers px and pw to be decremented by 1
double w_nth_int_ord(const int *restrict px, const double *restrict pw, const int *restrict po, double h, int l, const int narm, const int ret, const double Q) {
  if(l <= 1) {
    if(l == 0) return NA_REAL;
    return ISNAN(pw[po[0]]) ? NA_REAL : (double)px[po[0]];
  }
  if(narm) { // Adjusting l as necessary... initial NA check done in fnthC()
    while(l != 0 && px[po[l-1]] == NA_INTEGER) --l;
    if(l <= 1) return (l == 0 || ISNAN(pw[po[0]])) ? NA_REAL : (double)px[po[0]];
  } else if(px[po[l-1]] == NA_INTEGER) return NA_REAL;
  if(h == DBL_MIN) h = w_compute_h(pw, po, l, 0, ret, Q);
  if(ISNAN(h)) return NA_REAL;
  WNTH_CORE;
}

// Expects pointers px and pw to be decremented by 1
double w_nth_double_ord(const double *restrict px, const double *restrict pw, const int *restrict po, double h, int l, const int narm, const int ret, const double Q) {
  if(l <= 1) {
    if(l == 0) return NA_REAL;
    return ISNAN(pw[po[0]]) ? NA_REAL : px[po[0]];
  }
  if(narm) { // Adjusting l as necessary... initial NA check done in fnthC()
    while(l != 0 && ISNAN(px[po[l-1]])) --l;
    if(l <= 1) return (l == 0 || ISNAN(pw[po[0]])) ? NA_REAL : px[po[0]];
  } else if(ISNAN(px[po[l-1]])) return NA_REAL;
  if(h == DBL_MIN) h = w_compute_h(pw, po, l, 0, ret, Q);
  if(ISNAN(h)) return NA_REAL;
  WNTH_CORE;
}

// Quicksort versions: only for grouped execution (too slow on bigger vectors compared to radix sort)
// Expects pointer pw to be decremented by 1 if sorted == 0
double w_nth_int_qsort(const int *restrict px, const double *restrict pw, const int *restrict po, double h,
                       const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l <= 1) {
    if(l == 0) return NA_REAL;
    if(sorted) return ISNAN(pw[0]) ? NA_REAL : (double)px[0];
    return ISNAN(pw[po[0]]) ? NA_REAL : (double)px[po[0]-1];
  }

  int *x_cc = (int *) Calloc(l, int), *i_cc = (int *) Calloc(l, int), n = 0; // TODO: alloc i_cc afterwards if narm ??

  if(sorted) { // both the pointers to x and w need to be suitably incremented for grouped execution.
    // if(narm) {
      for(int i = 0; i != l; ++i) {
        if(px[i] != NA_INTEGER) {
          i_cc[n] = i;
          x_cc[n++] = px[i];
        }
      }
    // } else {
    //   n = l;
    //   for(int i = 0; i != l; ++i) {
    //     i_cc[i] = i;
    //     x_cc[i] = px[i];
    //   }
    // }
  } else {
    const int *pxm = px-1;
    // if(narm) {
      for(int i = 0; i != l; ++i) {
        if(pxm[po[i]] != NA_INTEGER) {
          i_cc[n] = po[i];
          x_cc[n++] = pxm[po[i]];
        }
      }
    // } else {
    //   n = l;
    //   for(int i = 0; i != l; ++i) {
    //     i_cc[i] = po[i];
    //     x_cc[i] = pxm[po[i]];
    //   }
    // }
  }

  if(narm == 0 && n != l) {
    Free(x_cc); Free(i_cc);
    return NA_REAL;
  }

  // i_cc is one-indexed
  R_qsort_int_I(x_cc, i_cc, 1, n);

  if(h == DBL_MIN) h = w_compute_h(pw, i_cc, n, 0, ret, Q);
  if(ISNAN(h)) {
    Free(x_cc); Free(i_cc);
    return NA_REAL;
  }

  WNTH_CORE_QSORT;

  Free(x_cc); Free(i_cc);
  return res;
}

// Expects pointer pw to be decremented by 1 if sorted == 0
double w_nth_double_qsort(const double *restrict px, const double *restrict pw, const int *restrict po, double h,
                          const int l, const int sorted, const int narm, const int ret, const double Q) {
  if(l <= 1) {
    if(l == 0) return NA_REAL;
    if(sorted) return ISNAN(pw[0]) ? NA_REAL : px[0];
    return ISNAN(pw[po[0]]) ? NA_REAL : px[po[0]-1];
  }

  double *x_cc = (double *) Calloc(l, double);
  int *i_cc = (int *) Calloc(l, int), n = 0; // TODO: alloc afterwards if narm ??

  if(sorted) {
    // if(narm) {
      for(int i = 0; i != l; ++i) {
        if(NISNAN(px[i])) {
          i_cc[n] = i;
          x_cc[n++] = px[i];
        }
      }
    // } else {
    //   n = l;
    //   for(int i = 0; i != l; ++i) {
    //     i_cc[i] = i;
    //     x_cc[i] = px[i];
    //   }
    // }
  } else {
    const double *pxm = px-1;
    // if(narm) {
      for(int i = 0; i != l; ++i) {
        if(NISNAN(pxm[po[i]])) {
          i_cc[n] = po[i];
          x_cc[n++] = pxm[po[i]];
        }
      }
    // } else {
    //   n = l;
    //   for(int i = 0; i != l; ++i) {
    //     i_cc[i] = po[i];
    //     x_cc[i] = pxm[po[i]];
    //   }
    // }
  }

  if(narm == 0 && n != l) {
    Free(x_cc); Free(i_cc);
    return NA_REAL;
  }

  // i_cc is one-indexed
  R_qsort_I(x_cc, i_cc, 1, n);

  if(h == DBL_MIN) h = w_compute_h(pw, i_cc, n, 0, ret, Q);

  if(ISNAN(h)) {
    Free(x_cc);
    Free(i_cc);
    return NA_REAL;
  }

  WNTH_CORE_QSORT;

  Free(x_cc);
  Free(i_cc);
  return res;
}



// Implementations for R vectors ---------------------------------------------------------------

// for safe multithreading in fnthlC()
SEXP nth_impl_plain(SEXP x, int narm, int ret, double Q) {
  int l = length(x);
  if(l <= 1) return x;

  switch(TYPEOF(x)) {
    case REALSXP: return ScalarReal(nth_double(REAL(x), &l, l, 1, narm, ret, Q));
    case INTSXP:
    case LGLSXP:  return ScalarReal(nth_int(INTEGER(x), &l, l, 1, narm, ret, Q));
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }
}

SEXP nth_impl(SEXP x, int narm, int ret, double Q) {
  if(length(x) <= 1) return x;
  if(ATTRIB(x) == R_NilValue || (isObject(x) && inherits(x, "ts")))
     return nth_impl_plain(x, narm, ret, Q);
  SEXP res = PROTECT(nth_impl_plain(x, narm, ret, Q));
  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

// for safe multithreading in fnthlC()
double nth_impl_dbl(SEXP x, int narm, int ret, double Q) {
  int l = length(x);
  if(l < 1) return NA_REAL;
  switch(TYPEOF(x)) {
    case REALSXP: return nth_double(REAL(x), &l, l, 1, narm, ret, Q);
    case INTSXP:
    case LGLSXP: return nth_int(INTEGER(x), &l, l, 1, narm, ret, Q);
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }
}

// for safe multithreading in fnthlC()
SEXP nth_impl_noalloc_plain(SEXP x, void* x_cc, int narm, int ret, double Q) {
  int l = length(x);
  if(l <= 1) return x;
  switch(TYPEOF(x)) {
    case REALSXP: return ScalarReal(nth_double_noalloc(REAL(x), &l, x_cc, l, 1, narm, ret, Q));
    case INTSXP:
    case LGLSXP: return ScalarReal(nth_int_noalloc(INTEGER(x), &l, x_cc, l, 1, narm, ret, Q));
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }
}

double nth_impl_noalloc_dbl(SEXP x, void* x_cc, int narm, int ret, double Q) {
  int l = length(x);
  if(l < 1) return NA_REAL;
  switch(TYPEOF(x)) {
    case REALSXP: return nth_double_noalloc(REAL(x), &l, x_cc, l, 1, narm, ret, Q);
    case INTSXP:
    case LGLSXP: return nth_int_noalloc(INTEGER(x), &l, x_cc, l, 1, narm, ret, Q);
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }
}

SEXP nth_ord_impl(SEXP x, int *pxo, int narm, int ret, double Q) {
  int l = length(x);
  if(l <= 1) return x;

  SEXP res;
  switch(TYPEOF(x)) {
    case REALSXP:
      res = ScalarReal(nth_double_ord(REAL(x)-1, pxo, l, narm, ret, Q));
      break;
    case INTSXP:
    case LGLSXP:
      res = ScalarReal(nth_int_ord(INTEGER(x)-1, pxo, l, narm, ret, Q));
      break;
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }

  if(ATTRIB(x) == R_NilValue || (isObject(x) && inherits(x, "ts"))) return res;
  PROTECT(res); // Needed ??
  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

// Expects pointer pw to be decremented by 1
SEXP w_nth_ord_impl_plain(SEXP x, int *pxo, double *pw, int narm, int ret, double Q, double h) {
  int l = length(x);
  if(l <= 1) return x;

  switch(TYPEOF(x)) {
    case REALSXP: return ScalarReal(w_nth_double_ord(REAL(x)-1, pw, pxo, h, l, narm, ret, Q));
    case INTSXP:
    case LGLSXP:  return ScalarReal(w_nth_int_ord(INTEGER(x)-1, pw, pxo, h, l, narm, ret, Q));
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }
}

// Expects pointer pw to be decremented by 1
SEXP w_nth_ord_impl(SEXP x, int *pxo, double *pw, int narm, int ret, double Q, double h) {
  if(length(x) <= 1) return x;
  if(ATTRIB(x) == R_NilValue || (isObject(x) && inherits(x, "ts")))
     return w_nth_ord_impl_plain(x, pxo, pw, narm, ret, Q, h);
  SEXP res = PROTECT(w_nth_ord_impl_plain(x, pxo, pw, narm, ret, Q, h));
  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

// Expects pointer pw to be decremented by 1
double w_nth_ord_impl_dbl(SEXP x, int *pxo, double *pw, int narm, int ret, double Q, double h) {
  int l = length(x);
  if(l < 1) return NA_REAL;
  switch(TYPEOF(x)) {
    case REALSXP: return w_nth_double_ord(REAL(x)-1, pw, pxo, h, l, narm, ret, Q);
    case INTSXP:
    case LGLSXP: return w_nth_int_ord(INTEGER(x)-1, pw, pxo, h, l, narm, ret, Q);
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }
}

// Expects pointer po to be decremented by 1
SEXP nth_g_impl(SEXP x, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int ret, double Q, int nthreads) {

  if(nthreads > ng) nthreads = ng;

  // TODO: if nthreads = 1, pass x_cc array of size maxgrpn repeatedly to the functions!!
  SEXP res = PROTECT(allocVector(REALSXP, ng));
  double *pres = REAL(res);

  if(sorted) { // Sorted: could compute cumulative group size (= starts) on the fly... but doesn't work multithreaded...
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x)-1;
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = nth_double(px + pst[gr], po, pgs[gr], 1, narm, ret, Q);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x)-1;
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = nth_int(px + pst[gr], po, pgs[gr], 1, narm, ret, Q);
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
          pres[gr] = nth_double(px, po + pst[gr], pgs[gr], 0, narm, ret, Q);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = nth_int(px, po + pst[gr], pgs[gr], 0, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  }

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

// Expects pointer po to be decremented by 1
SEXP nth_g_impl_noalloc(SEXP x, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int ret, double Q, void* x_cc) {

  SEXP res = PROTECT(allocVector(REALSXP, ng));
  double *pres = REAL(res);

  if(sorted) {
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x)-1;
        for(int gr = 0; gr != ng; ++gr) pres[gr] = nth_double_noalloc(px + pst[gr], po, x_cc, pgs[gr], 1, narm, ret, Q);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x)-1;
        for(int gr = 0; gr != ng; ++gr) pres[gr] = nth_int_noalloc(px + pst[gr], po, x_cc, pgs[gr], 1, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  } else {
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        for(int gr = 0; gr != ng; ++gr) pres[gr] = nth_double_noalloc(px, po + pst[gr], x_cc, pgs[gr], 0, narm, ret, Q);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        for(int gr = 0; gr != ng; ++gr) pres[gr] = nth_int_noalloc(px, po + pst[gr], x_cc, pgs[gr], 0, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  }

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

// Expects pointer po to be decremented by 1
SEXP nth_g_ord_impl(SEXP x, int ng, int *pgs, int *po, int *pst, int narm, int ret, double Q, int nthreads) {

  if(nthreads > ng) nthreads = ng;

  SEXP res = PROTECT(allocVector(REALSXP, ng));
  double *pres = REAL(res);

  switch(TYPEOF(x)) {
    case REALSXP: {
      double *px = REAL(x)-1;
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = nth_double_ord(px, po + pst[gr], pgs[gr], narm, ret, Q);
      break;
    }
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x)-1;
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = nth_int_ord(px, po + pst[gr], pgs[gr], narm, ret, Q);
      break;
    }
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

// Expects pointers pw and po to be decremented by 1
SEXP w_nth_g_ord_impl(SEXP x, double *pw, int ng, int *pgs, int *po, int *pst, int narm, int ret, double Q, int nthreads) {

  if(nthreads > ng) nthreads = ng;

  SEXP res = PROTECT(allocVector(REALSXP, ng));
  double *pres = REAL(res);

  switch(TYPEOF(x)) {
    case REALSXP: {
      double *px = REAL(x)-1;
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = w_nth_double_ord(px, pw, po + pst[gr], DBL_MIN, pgs[gr], narm, ret, Q);
      break;
    }
    case INTSXP:
    case LGLSXP: {
      int *px = INTEGER(x)-1;
      #pragma omp parallel for num_threads(nthreads)
      for(int gr = 0; gr < ng; ++gr)
        pres[gr] = w_nth_int_ord(px, pw, po + pst[gr], DBL_MIN, pgs[gr], narm, ret, Q);
      break;
    }
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

// Expects pointers pw and po to be decremented by 1
SEXP w_nth_g_qsort_impl(SEXP x, double *pw, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int ret, double Q, int nthreads) {

  if(nthreads > ng) nthreads = ng;

  SEXP res = PROTECT(allocVector(REALSXP, ng));
  double *pres = REAL(res);

  if(sorted) { // sorted by groups: need to offset both px and pw
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x)-1;
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = w_nth_double_qsort(px + pst[gr], pw + pst[gr], po, DBL_MIN, pgs[gr], 1, narm, ret, Q);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x)-1;
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = w_nth_int_qsort(px + pst[gr], pw + pst[gr], po, DBL_MIN, pgs[gr], 1, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  } else {
    switch(TYPEOF(x)) {
      case REALSXP: {
        double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = w_nth_double_qsort(px, pw, po + pst[gr], DBL_MIN, pgs[gr], 0, narm, ret, Q);
        break;
      }
      case INTSXP:
      case LGLSXP: {
        int *px = INTEGER(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = w_nth_int_qsort(px, pw, po + pst[gr], DBL_MIN, pgs[gr], 0, narm, ret, Q);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
    }
  }

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}



// Functions for Export --------------------------------------------------------

int Rties2int(SEXP x) {
  int tx = TYPEOF(x);
  if(tx == INTSXP || tx == REALSXP || tx == LGLSXP) {
    int ret = asInteger(x);
    if(ret < 1 || ret > 9 || ret == 4) error("ties must be 1, 2, 3 or 5-9, you supplied: %d", ret);
    return ret;
  }
  if(tx != STRSXP) error("ties must be integer or character");
  const char * r = CHAR(STRING_ELT(x, 0)); // translateCharUTF8()
  if(strcmp(r, "mean") == 0) return 1;
  if(strcmp(r, "min") == 0) return 2;
  if(strcmp(r, "max") == 0) return 3;
  if(strcmp(r, "q5") == 0) return 5;
  if(strcmp(r, "q6") == 0) return 6;
  if(strcmp(r, "q7") == 0) return 7;
  if(strcmp(r, "q8") == 0) return 8;
  if(strcmp(r, "q9") == 0) return 9;
  error("Unknown ties option: %s", r);
}

#undef CHECK_PROB
#define CHECK_PROB(l)                                                                                                                                                            \
  if(length(p) != 1) error("fnth supports only a single element / quantile. Use fquantile for multiple quantiles.");                                                             \
  double Q = asReal(p);                                                                                                                                                          \
  if(ISNAN(Q) || Q <= 0.0 || Q == 1.0) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");                             \
  if(Q > 1.0) {                                                                                                                                                                  \
    ret = 2; /* ties = "min" */                                                                                                                                                    \
    if(nullg) {                                                                                                                                                                  \
      if(Q >= l) error("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");                                                   \
      Q = (Q-1.0)/(l-1);                                                                                                                                                         \
    } else {                                                                                                                                                                     \
      if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");                                                                  \
      int ng = INTEGER(VECTOR_ELT(g, 0))[0];                                                                                                                                     \
      if(Q >= (double)l/ng) error("n needs to be between 0 and 1, or between 1 and the length(x)/ng, with ng the number of groups. Use fmin and fmax for minima and maxima.");   \
      Q = (Q-1.0)/((double)l/ng-1.0);                                                                                                                                            \
    }                                                                                                                                                                            \
  }

#undef CHECK_WEIGHTS
#define CHECK_WEIGHTS(l)                                         \
  if(length(w) != l) error("length(w) must match length(x)");    \
  if(TYPEOF(w) != REALSXP) {                                     \
    if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double). You supplied a vector of type: '%s'", type2char(TYPEOF(w))); \
    w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;           \
  }                                                              \
  pw = REAL(w)-1; /* All functions require decremented w pointer */


#undef CHECK_GROUPS
#define CHECK_GROUPS(nrx, cond)                                                                                \
if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");      \
const SEXP *restrict pg = SEXPPTR(g), ord = pg[6];                                                             \
ng = INTEGER(pg[0])[0];                                                                                        \
int sorted = LOGICAL(pg[5])[1] == 1, *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst, maxgrpn = 0; \
if(nrx != length(pg[1])) error("length(g) must match nrow(x)");                                                \
if(isNull(ord)) {                                                                                              \
  int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;                   \
  if(nthreads <= 1 && nullw) {                                                                                 \
    for(int i = 0; i != ng; ++i) {                                                                             \
      if(pgs[i] > maxgrpn) maxgrpn = pgs[i];                                                                   \
      cgs[i+2] = cgs[i+1] + pgs[i];                                                                            \
    }                                                                                                          \
  } else {                                                                                                     \
    for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];                                                 \
  }                                                                                                            \
  pst = cgs + 1;                                                                                               \
  if((cond)) po = &l;                                                                                          \
  else {                                                                                                       \
    int *restrict count = (int *) Calloc(ng+1, int);                                                           \
    po = (int *) R_alloc(nrx, sizeof(int)); --po;                                                              \
    for(int i = 0; i != nrx; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;                                     \
    Free(count);                                                                                               \
  }                                                                                                            \
} else {                                                                                                       \
  po = INTEGER(ord)-1;                                                                                         \
  pst = INTEGER(getAttrib(ord, install("starts")));                                                            \
  if(nthreads <= 1 && nullw) maxgrpn = asInteger(getAttrib(ord, install("maxgrpn")));                          \
}


/*
   Function for atomic vectors: has extra arguments o and checko for passing external ordering vector.
   This is meant to speed up computation of several (grouped) quantiles on the same data.
   Note that for grouped execution the ordering vector needs to take into account the grouping e.g radixorder(GRPid(), myvar).
 */
SEXP fnthC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rret, SEXP Rnthreads, SEXP o, SEXP checko) {

  int nullg = isNull(g), nullw = isNull(w), nullo = isNull(o), l = length(x), narm = asLogical(Rnarm),
      ret = Rties2int(Rret), nprotect = 0;

  CHECK_PROB(l);

  // if(l < 1) return x;
  if(l < 1 || (l == 1 && nullw)) return TYPEOF(x) == REALSXP ? x : ScalarReal(asReal(x));

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
    if((TYPEOF(x) == REALSXP && ISNAN(REAL(x)[pxo[0]-1])) || ((TYPEOF(x) == INTSXP || TYPEOF(x) == LGLSXP) && INTEGER(x)[pxo[0]-1] == NA_INTEGER))
      error("Found missing value at the beginning of the sample. Please use option na.last = TRUE (the default) when creasting ordering vectors for use with fnth().");
  }

  // Preprocessing w, computing ordering of x if not supplied
  if(!nullw) {
    CHECK_WEIGHTS(l);
    if(l == 1) {
      if(ISNAN(pw[1])) return ScalarReal(NA_REAL);
      return TYPEOF(x) == REALSXP ? x : ScalarReal(asReal(x));
    }
    if(nullo && nullg) { // for grouped execution use w_nth_g_qsort_impl() if o is not supplied.
      // nullo = 0;
      pxo = (int *) R_alloc(l, sizeof(int));
      num1radixsort(pxo, TRUE, FALSE, x);
    }
  }

  // If no groups, return using suitable functions
  if(nullg) {
    SEXP res; // result, could be put outside if() to avoid repetition below, but this seems to confuse rchk
    if(nullw) res = nth_ord_impl(x, pxo, narm, ret, Q);
    else res = w_nth_ord_impl(x, pxo, pw, narm, ret, Q, DBL_MIN);
    UNPROTECT(nprotect);
    return res;
  }

  int nthreads = asInteger(Rnthreads), ng;
  if(nthreads > max_threads) nthreads = max_threads;

  // Preprocessing g
  CHECK_GROUPS(l, sorted || !nullo);
  /*
   * Previous version: computes po if overall ordering of x is supplied to o. This is made redundant by requiring
   * the ordering o to now take into account the grouping (facilitated by R-level helper GRPid()), which provides
   * much greater speedup for repeated executions, and by the addition of w_nth_g_qsort_impl().
   *
  if((!nullw && nullo) || isNull(ord)) { // Extra case: if ordering vector supplied, need to use it to get the group elements in order
    int *restrict pgv = INTEGER(pg[1]);
    if(isNull(ord)) {
      int *cgs = (int *) R_alloc(ng+2, sizeof(int)); cgs[1] = 1;
      for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i]; // TODO: get maxgrpn?
      pst = cgs;
    } else pst = INTEGER(getAttrib(ord, install("starts")))-1;
    if(nullw && sorted) po = &l;
    else {
      int *restrict count = (int *) Calloc(ng+1, int);
      po = (int *) R_alloc(l, sizeof(int)); --po;
      if(nullw) {
        for(int i = 0; i != l; ++i) po[pst[pgv[i]] + count[pgv[i]]++] = i+1;
      } else { // This orders the elements of x within groups... e.g. starting with the first group, the indices of all elements of x in order, then the second group etc.
        --pgv;
        for(int i = 0, tmp; i != l; ++i) {
          tmp = pgv[pxo[i]];
          po[pst[tmp] + count[tmp]++] = pxo[i];
        }
      }
      Free(count);
    }
    ++pst;
  }
   */

  SEXP res; // result
  if(nullw && nullo) res = nthreads <= 1 ? nth_g_impl_noalloc(x, ng, pgs, po, pst, sorted, narm, ret, Q, R_alloc(maxgrpn, TYPEOF(x) == REALSXP ? sizeof(double) : sizeof(int))) :
                                           nth_g_impl(x, ng, pgs, po, pst, sorted, narm, ret, Q, nthreads);
  else if(nullw) res = nth_g_ord_impl(x, ng, pgs, pxo-1, pst, narm, ret, Q, nthreads);
  else if(nullo) res = w_nth_g_qsort_impl(x, pw, ng, pgs, po, pst, sorted, narm, ret, Q, nthreads);
  else res = w_nth_g_ord_impl(x, pw, ng, pgs, pxo-1, pst, narm, ret, Q, nthreads);

  UNPROTECT(nprotect);
  return res;
}



#undef COLWISE_NTH_LIST
#define COLWISE_NTH_LIST(FUN_NA, FUN, WFUN)                    \
if(nullw) {                                                    \
  if(nthreads == 1) {                                          \
    void *x_cc = Calloc(nrx, double);                          \
    for(int j = 0; j != l; ++j) pout[j] = FUN_NA(px[j], x_cc, narm, ret, Q); \
    Free(x_cc);                                                \
  } else {                                                     \
    _Pragma("omp parallel for num_threads(nthreads)")          \
    for(int j = 0; j < l; ++j) pout[j] = FUN(px[j], narm, ret, Q); \
  }                                                            \
} else { /* TODO: if narm = FALSE, can compute sumw beforehand */ \
  int *pxo = (int *) R_alloc(nrx, sizeof(int));                \
  for(int j = 0; j != l; ++j) {                                \
    num1radixsort(pxo, TRUE, FALSE, px[j]);                    \
    pout[j] = WFUN(px[j], pxo, pw, narm, ret, Q, h);           \
  }                                                            \
}
/* Multithreading: does not work with radixorder
 * } else {
   #pragma omp parallel for num_threads(nthreads)
   for(int j = 0; j < l; ++j) {
   int *pxo = (int *) Calloc(nrx, int);
   // num1radixsort(pxo, TRUE, FALSE, px[j]); // Probably cannot be parallelized, can try R_orderVector1()
   // R_orderVector1(pxo, nrx, px[j], TRUE, FALSE); // Also not thread safe, and also 0-indexed.
   // for(int i = 0; i < nrx; ++i) pxo[i] += 1;
   pout[j] = w_nth_ord_impl_dbl(px[j], pxo, pw, narm, ret, Q, h);
   Free(pxo);
   }
  }
   */

// TODO: Pre-compute weights at the group-level if narm = FALSE for list and matrix method

// Function for lists / data frames
SEXP fnthlC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rret, SEXP Rnthreads) {

  int nullg = isNull(g), nullw = isNull(w), l = length(x), ng = 0, nprotect = 1,
    narm = asLogical(Rnarm), drop = asLogical(Rdrop), ret = Rties2int(Rret), nthreads = asInteger(Rnthreads);

  if(l < 1) return x;
  if(nthreads > max_threads) nthreads = max_threads;

  SEXP out = PROTECT(allocVector(nullg && drop ? REALSXP : VECSXP, l)), *restrict px = SEXPPTR(x);
  int nrx = length(px[0]);

  CHECK_PROB(nrx);

  double *restrict pw = &Q, h = DBL_MIN;
  if(!nullw) {
    CHECK_WEIGHTS(nrx);
    if(nullg && !narm) h = w_compute_h(pw+1, &l, nrx, 1, ret, Q); // if no missing value removal, h is the same for all columns
  }

  if(nullg) { // No groups, multithreading across columns
    if(nthreads > l) nthreads = l;
    if(drop) { // drop dimensions (return vector)
      double *restrict pout = REAL(out);
      COLWISE_NTH_LIST(nth_impl_noalloc_dbl, nth_impl_dbl, w_nth_ord_impl_dbl);
      setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
      UNPROTECT(nprotect);
      return out;
    }
    // returns a list of atomic elements
    SEXP *restrict pout = SEXPPTR(out);
    COLWISE_NTH_LIST(nth_impl_noalloc_plain, nth_impl_plain, w_nth_ord_impl_plain);
    // Needed because including it in an OpenMP loop together with ScalarReal() is not thread safe
    for(int j = 0; j != l; ++j) {
      SEXP xj = px[j];
      if(ATTRIB(xj) != R_NilValue && !(isObject(xj) && inherits(xj, "ts")))
        copyMostAttrib(xj, pout[j]);
    }

  } else { // with groups: do the usual checking

    CHECK_GROUPS(nrx, sorted);

    SEXP *restrict pout = SEXPPTR(out);
    if(nullw) { // Parallelism at sub-column level
      if(nthreads <= 1) {
        void *x_cc = R_alloc(maxgrpn, sizeof(double));
        for(int j = 0; j < l; ++j) pout[j] = nth_g_impl_noalloc(px[j], ng, pgs, po, pst, sorted, narm, ret, Q, x_cc);
      } else {
        for(int j = 0; j < l; ++j) pout[j] = nth_g_impl(px[j], ng, pgs, po, pst, sorted, narm, ret, Q, nthreads);
      }
    } else { // Parallelism at sub-column level
      for(int j = 0; j < l; ++j) pout[j] = w_nth_g_qsort_impl(px[j], pw, ng, pgs, po, pst, sorted, narm, ret, Q, nthreads);
    }
  }

  DFcopyAttr(out, x, ng);
  UNPROTECT(nprotect);
  return out;
}


// Iterate over matrix columns: for integers and doubles
#undef COLWISE_NTH
#define COLWISE_NTH(tdef, FUN, FUN_NA, WFUN, ORDFUN)                                 \
  if(nullw) {                                                                        \
    if(nthreads == 1) {                                                              \
      tdef *x_cc = (tdef *) R_alloc(l, sizeof(tdef));                                \
      for(int j = 0; j < col; ++j) pres[j] = FUN_NA(px + j*l, &l, x_cc, l, 1, narm, ret, Q);  \
    } else {                                                                         \
      _Pragma("omp parallel for num_threads(nthreads)")                              \
      for(int j = 0; j < col; ++j) pres[j] = FUN(px + j*l, &l, l, 1, narm, ret, Q);  \
    }                                                                                \
  } else {                                                                           \
    /* if(nthreads == 1) { */                                                        \
      int *pxo = (int *) R_alloc(l, sizeof(int));                                    \
      for(int j = 0; j < col; ++j) {                                                 \
        ORDFUN(pxo, TRUE, FALSE, l, px + j*l);                                       \
        pres[j] = WFUN(px + j*l - 1, pw, pxo, h, l, narm, ret, Q);                   \
      }                                                                              \
  } /* else {                                                                        \
      _Pragma("omp parallel for num_threads(nthreads)")                              \
      for(int j = 0; j < col; ++j) {                                                 \
        int *pxo = (int *) Calloc(l, int);                                           \
        ORDFUN(pxo, TRUE, FALSE, l, px + j*l); // Currently cannot be parallelized   \
        pres[j] = WFUN(px + j*l - 1, pw, pxo, h, l, narm, ret, Q);                   \
        Free(pxo);                                                                   \
      }                                                                              \
    }                                                                                \
  }                                                            \
     */

// The same by groups if data already sorted by groups. px and pw should be decremented by 1
#undef COLWISE_NTH_GROUPED_SORTED
#define COLWISE_NTH_GROUPED_SORTED(tdef, FUN, FUN_NA, WFUN)                          \
if(nullw) {                                                                          \
  if(nthreads == 1) {                                                                \
    tdef *x_cc = (tdef *) R_alloc(maxgrpn, sizeof(tdef));                            \
    for(int j = 0; j != col; ++j) {                                                  \
      int jng = j * ng;                                                              \
      tdef *pxj = px + j * l;                                                        \
      for(int gr = 0; gr != ng; ++gr)                                                \
        pres[jng + gr] = FUN_NA(pxj + pst[gr], po, x_cc, pgs[gr], 1, narm, ret, Q);  \
    }                                                                                \
  } else {                                                                           \
    _Pragma("omp parallel for num_threads(nthreads)")                                \
    for(int j = 0; j < col; ++j) {                                                   \
      int jng = j * ng;                                                              \
      tdef *pxj = px + j * l;                                                        \
      for(int gr = 0; gr < ng; ++gr)                                                 \
        pres[jng + gr] = FUN(pxj + pst[gr], po, pgs[gr], 1, narm, ret, Q);           \
    }                                                                                \
  }                                                                                  \
} else {                                                                             \
  if(nthreads == 1) {                                                                \
    for(int j = 0; j != col; ++j) {                                                  \
      int jng = j * ng;                                                              \
      tdef *pxj = px + j * l;                                                        \
      for(int gr = 0; gr != ng; ++gr)                                                \
        pres[jng + gr] = WFUN(pxj + pst[gr], pw + pst[gr], po, DBL_MIN, pgs[gr], 1, narm, ret, Q); \
    }                                                                                \
  } else {                                                                           \
    _Pragma("omp parallel for num_threads(nthreads)")                                \
    for(int j = 0; j < col; ++j) {                                                   \
      int jng = j * ng;                                                              \
      tdef *pxj = px + j * l;                                                        \
      for(int gr = 0; gr < ng; ++gr)                                                 \
        pres[jng + gr] = WFUN(pxj + pst[gr], pw + pst[gr], po, DBL_MIN, pgs[gr], 1, narm, ret, Q); \
    }                                                                                \
  }                                                                                  \
}

// The more general case. po should be decremented by 1.
#undef COLWISE_NTH_GROUPED_UNSORTED
#define COLWISE_NTH_GROUPED_UNSORTED(tdef, FUN, FUN_NA, WFUN)                                      \
if(nullw) {                                                                                        \
  if(nthreads == 1) {                                                                              \
    tdef *x_cc = (tdef *) R_alloc(maxgrpn, sizeof(tdef));                                          \
    for(int j = 0; j != col; ++j) {                                                                \
      int jng = j * ng;                                                                            \
      tdef *pxj = px + j * l;                                                                      \
      for(int gr = 0; gr != ng; ++gr)                                                              \
        pres[jng + gr] = FUN_NA(pxj, po + pst[gr], x_cc, pgs[gr], 0, narm, ret, Q);                \
    }                                                                                              \
  } else {                                                                                         \
    _Pragma("omp parallel for num_threads(nthreads)")                                              \
    for(int j = 0; j < col; ++j) {                                                                 \
      int jng = j * ng;                                                                            \
      tdef *pxj = px + j * l;                                                                      \
      for(int gr = 0; gr < ng; ++gr)                                                               \
        pres[jng + gr] = FUN(pxj, po + pst[gr], pgs[gr], 0, narm, ret, Q);                         \
    }                                                                                              \
  }                                                                                                \
} else {                                                                                           \
  if(nthreads == 1) {                                                                              \
    for(int j = 0; j != col; ++j) {                                                                \
      int jng = j * ng;                                                                            \
      tdef *pxj = px + j * l;                                                                      \
      for(int gr = 0; gr != ng; ++gr)                                                              \
        pres[jng + gr] = WFUN(pxj, pw, po + pst[gr], DBL_MIN, pgs[gr], 0, narm, ret, Q);           \
    }                                                                                              \
  } else {                                                                                         \
    _Pragma("omp parallel for num_threads(nthreads)")                                              \
    for(int j = 0; j < col; ++j) {                                                                 \
      int jng = j * ng;                                                                            \
      tdef *pxj = px + j * l;                                                                      \
      for(int gr = 0; gr < ng; ++gr)                                                               \
        pres[jng + gr] = WFUN(pxj, pw, po + pst[gr], DBL_MIN, pgs[gr], 0, narm, ret, Q);           \
    }                                                                                              \
  }                                                                                                \
}

// Function for matrices: implemented at lower-level
SEXP fnthmC(SEXP x, SEXP p, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rret, SEXP Rnthreads) {

  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");

  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1],
      narm = asLogical(Rnarm), ret = Rties2int(Rret), nthreads = asInteger(Rnthreads),
      nullg = isNull(g), nullw = isNull(w), nprotect = 1;

  if(nthreads > col) nthreads = col;
  if(nthreads > max_threads) nthreads = max_threads;

  CHECK_PROB(l);

  if(l < 1 || (l == 1 && nullw)) {
    if(TYPEOF(x) == REALSXP || TYPEOF(x) == INTSXP || TYPEOF(x) == LGLSXP) return x;
    error("Unsopported SEXP type: '%s'", type2char(TYPEOF(x)));
  }

  double *restrict pw = &Q, h = DBL_MIN;
  if(!nullw) {
    CHECK_WEIGHTS(l);
    if(nullg && !narm) h = w_compute_h(pw+1, &l, l, 1, ret, Q);
  }

  if(nullg) {
    SEXP res = PROTECT(allocVector(REALSXP, col));

    switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *restrict pres = REAL(res);
        COLWISE_NTH(double, nth_double, nth_double_noalloc, w_nth_double_ord, dradixsort);
        break;
      }
      case INTSXP:
      case LGLSXP: {  // Factor matrix not well defined object...
        int *px = INTEGER(x), *restrict pres = INTEGER(res);
        COLWISE_NTH(int, nth_int, nth_int_noalloc, w_nth_int_ord, iradixsort);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(tx));
    }

    matCopyAttr(res, x, Rdrop, /*ng=*/0);
    UNPROTECT(nprotect);
    return res;
  }

  // With groups
  int ng;
  CHECK_GROUPS(l, sorted);

  SEXP res = PROTECT(allocVector(REALSXP, col * ng));

  if(sorted) { // Sorted
    switch(tx) {
      case REALSXP: {
        double *px = REAL(x)-1, *restrict pres = REAL(res);
        COLWISE_NTH_GROUPED_SORTED(double, nth_double, nth_double_noalloc, w_nth_double_qsort);
        break;
      }
      case INTSXP:
      case LGLSXP: { // Factor matrix not well defined object...
        int *px = INTEGER(x)-1, *restrict pres = INTEGER(res);
        COLWISE_NTH_GROUPED_SORTED(int, nth_int, nth_int_noalloc, w_nth_int_qsort);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(tx));
    }
  } else { // Not sorted
    switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *restrict pres = REAL(res);
        COLWISE_NTH_GROUPED_UNSORTED(double, nth_double, nth_double_noalloc, w_nth_double_qsort);
        break;
      }
      case INTSXP:
      case LGLSXP: { // Factor matrix not well defined object...
        int *px = INTEGER(x), *restrict pres = INTEGER(res);
        COLWISE_NTH_GROUPED_UNSORTED(int, nth_int, nth_int_noalloc, w_nth_int_qsort);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(tx));
    }
  }

  matCopyAttr(res, x, Rdrop, ng);
  UNPROTECT(nprotect);
  return res;
}

