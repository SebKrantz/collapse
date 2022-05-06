#include "kit.h"
#include "collapse_c.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// C-implementations for different data types ----------------------------------
// TODO: outsource and memset hash table?
// Problem: does not work in parallel, each thread needs own hash table...

int ndistinct_int(const int *px, const int *po, const int n, const int sorted, const int narm) {
  const size_t n2 = 2U * (size_t) n;
  size_t M = 256, id = 0;
  int K = 8, ndist = 0, anyNA = 0;
  while(M < n2) {
    M *= 2;
    K++;
  }
  int *h = (int*)Calloc(M, int); // Table to save the hash values, table has size M

  if(sorted) {
    for (int i = 0; i != n; ++i) {
      if(px[i] == NA_INTEGER) {
        anyNA = 1;
        continue;
      }
      id = HASH(px[i], K);
      while(h[id]) {
        if(px[h[id]-1] == px[i]) goto ibls;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      ++ndist;
      ibls:;
    }
  } else {
    for (int i = 0, xi; i != n; ++i) {
      xi = px[po[i]-1];
      if(xi == NA_INTEGER) {
        anyNA = 1;
        continue;
      }
      id = HASH(xi, K);
      while(h[id]) {
        if(px[po[h[id]-1]-1] == xi) goto ibl;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      ++ndist;
      ibl:;
    }
  }

  Free(h);
  if(narm == 0) ndist += anyNA;
  return ndist;
}

int ndistinct_fct(const int *px, const int *po, const int n, const int nlev, const int sorted, const int narm) {
  int *h = (int*)Calloc(nlev+1, int);
  int ndist = 0, anyNA = narm; // Ensures breaking works if narm = TRUE or FALSE
  if(sorted) {
    for (int i = 0, xi; i != n; ++i) {
      xi = px[i];
      if(xi == NA_INTEGER) {
        anyNA = 1;
        continue;
      }
      if(h[xi]) continue;
      ++ndist;
      if(anyNA && ndist == nlev) break;
      h[xi] = 1;
    }
  } else {
    for (int i = 0, xi; i != n; ++i) {
      xi = px[po[i]-1];
      if(xi == NA_INTEGER) {
        anyNA = 1;
        continue;
      }
      if(h[xi]) continue;
      ++ndist;
      if(anyNA && ndist == nlev) break;
      h[xi] = 1;
    }
  }
  if(narm == 0) ndist += anyNA;
  Free(h);
  return ndist;
}

int ndistinct_logi(const int *px, const int *po, const int n, const int sorted, const int narm) {
  int seenT = 0, seenF = 0, anyNA = narm; // Ensures breaking works if narm = TRUE or FALSE
  if(sorted) {
    for (int i = 0, xi; i != n; ++i) {
      xi = px[i];
      if(xi == NA_LOGICAL) {
        anyNA = 1;
      } else if(xi) {
        if(seenT) continue;
        seenT = 1;
        if(anyNA && seenF) break;
      } else {
        if(seenF) continue;
        seenF = 1;
        if(anyNA && seenT) break;
      }
    }
  } else {
    for (int i = 0, xi; i != n; ++i) {
      xi = px[po[i]-1];
      if(xi == NA_LOGICAL) {
        anyNA = 1;
      } else if(xi) {
        if(seenT) continue;
        seenT = 1;
        if(anyNA && seenF) break;
      } else {
        if(seenF) continue;
        seenF = 1;
        if(anyNA && seenT) break;
      }
    }
  }

  if(narm == 0) seenT += anyNA;
  return seenT + seenF;
}

int ndistinct_double(const double *px, const int *po, const int n, const int sorted, const int narm) {
  const size_t n2 = 2U * (size_t) n;
  size_t M = 256, id = 0;
  int K = 8, ndist = 0, anyNA = 0;
  while(M < n2) {
    M *= 2;
    K++;
  }
  int *h = (int*)Calloc(M, int); // Table to save the hash values, table has size M
  union uno tpv;
  double xi;

  if(sorted) {
    for (int i = 0; i != n; ++i) {
      if(ISNAN(px[i])) {
        anyNA = 1;
        continue;
      }
      tpv.d = px[i];
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        if(REQUAL(px[h[id]-1], px[i])) goto rbls;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      ++ndist;
      rbls:;
    }
  } else {
    for (int i = 0; i != n; ++i) {
      xi = px[po[i]-1];
      if(ISNAN(xi)) {
        anyNA = 1;
        continue;
      }
      tpv.d = xi;
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        if(REQUAL(px[po[h[id]-1]-1], xi)) goto rbl;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      ++ndist;
      rbl:;
    }
  }


  Free(h);
  if(narm == 0) ndist += anyNA;
  return ndist;
}

int ndistinct_string(const SEXP *px, const int *po, const int n, const int sorted, const int narm) {
  const size_t n2 = 2U * (size_t) n;
  size_t M = 256, id = 0;
  int K = 8, ndist = 0, anyNA = 0;
  while(M < n2) {
    M *= 2;
    K++;
  }
  int *h = (int*)Calloc(M, int); // Table to save the hash values, table has size M
  SEXP xi;

  if(sorted) {
    for (int i = 0; i != n; ++i) {
      if(px[i] == NA_STRING) {
        anyNA = 1;
        continue;
      }
      id = HASH(((intptr_t) px[i] & 0xffffffff), K);
      while(h[id]) {
        if(px[h[id]-1] == px[i]) goto sbls;
        if(++id >= M) id %= M; //++id; id %= M;
      }
      h[id] = i + 1;
      ++ndist;
      sbls:;
    }
  } else {
    for (int i = 0; i != n; ++i) {
      xi = px[po[i]-1];
      if(xi == NA_STRING) {
        anyNA = 1;
        continue;
      }
      id = HASH(((intptr_t) xi & 0xffffffff), K);
      while(h[id]) {
        if(px[po[h[id]-1]-1] == xi) goto sbl;
        if(++id >= M) id %= M; //++id; id %= M;
      }
      h[id] = i + 1;
      ++ndist;
      sbl:;
    }
  }

  Free(h);
  if(narm == 0) ndist += anyNA;
  return ndist;
}

// Implementations for R vectors -----------------------------------------------

SEXP ndistinct_impl(SEXP x, int narm) {
  int l = length(x), res;
  if(l < 1) return ScalarInteger(0);
  switch(TYPEOF(x)) {
    case REALSXP:
      res = ndistinct_double(REAL(x), &l, l, 1, narm);
      break;
  case INTSXP:  // TODO: optimize for plain integer??
      res = isFactor(x) ? ndistinct_fct(INTEGER(x), &l, l, nlevels(x), 1, narm) :
               ndistinct_int(INTEGER(x), &l, l, 1, narm);
      break;
    case LGLSXP:
      res = ndistinct_logi(INTEGER(x), &l, l, 1, narm);
      break;
    case STRSXP:
      res = ndistinct_string(STRING_PTR(x), &l, l, 1, narm);
      break;
    default: error("Not Supported SEXP Type!");
  }

  return ScalarInteger(res);
}

// TODO: Optimize grouped distinct value count for logical vectors??
SEXP ndistinct_g_impl(SEXP x, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int nthreads) {

  SEXP res = PROTECT(allocVector(INTSXP, ng));
  int l = length(x), *pres = INTEGER(res);
  if(nthreads > ng) nthreads = ng;

  if(sorted) { // Sorted: compute cumulative group size (= starts) on the fly...
    po = &l;
    int gs = 0, gsgr = 0;
    switch(TYPEOF(x)) {
      case REALSXP: {
        const double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr) {
          gsgr = pgs[gr];
          pres[gr] = gsgr == 0 ? 0 : ndistinct_double(px + gs, po, gsgr, 1, narm);
          gs += gsgr;
        }
        break;
      }
      case INTSXP: {
        const int *px = INTEGER(x);
        if(isFactor(x) && (ng == 0 || nlevels(x) < l / ng * 3)) {
          int M = nlevels(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr) {
            gsgr = pgs[gr];
            pres[gr] = gsgr == 0 ? 0 : ndistinct_fct(px + gs, po, gsgr, M, 1, narm);
            gs += gsgr;
          }
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr) {
            gsgr = pgs[gr];
            pres[gr] = gsgr == 0 ? 0 : ndistinct_int(px + gs, po, gsgr, 1, narm);
            gs += gsgr;
          }
        }
        break;
      }
      case LGLSXP: {
        const int *px = LOGICAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr) {
            gsgr = pgs[gr];
            pres[gr] = gsgr == 0 ? 0 : ndistinct_logi(px + gs, po, gsgr, 1, narm);
            gs += gsgr;
        }
        break;
      }
      case STRSXP: {
        const SEXP *px = STRING_PTR(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr) {
          gsgr = pgs[gr];
          pres[gr] = gsgr == 0 ? 0 : ndistinct_string(px + gs, po, gsgr, 1, narm);
          gs += gsgr;
        }
        break;
      }
      default: error("Not Supported SEXP Type!");
    }
  } else { // Not sorted. Perhaps reordering x is faster??
    switch(TYPEOF(x)) {
      case REALSXP: {
        const double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? 0 : ndistinct_double(px, po + pst[gr]-1, pgs[gr], 0, narm);
        break;
      }
      case INTSXP: {
        const int *px = INTEGER(x);
        if(isFactor(x) && (ng == 0 || nlevels(x) < l / ng * 3)) {
          int M = nlevels(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? 0 : ndistinct_fct(px, po + pst[gr]-1, pgs[gr], M, 0, narm);
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? 0 : ndistinct_int(px, po + pst[gr]-1, pgs[gr], 0, narm);
        }
        break;
      }
      case LGLSXP: {
        const int *px = LOGICAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? 0 : ndistinct_logi(px, po + pst[gr]-1, pgs[gr], 0, narm);
        break;
      }
      case STRSXP: {
        const SEXP *px = STRING_PTR(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? 0 : ndistinct_string(px, po + pst[gr]-1, pgs[gr], 0, narm);
        break;
      }
      default: error("Not Supported SEXP Type!");
    }
  }

  UNPROTECT(1);
  return res;
}

// Functions for Export --------------------------------------------------------

SEXP fndistinctC(SEXP x, SEXP g, SEXP Rnarm, SEXP Rnthreads) {
  if(isNull(g)) return ndistinct_impl(x, asLogical(Rnarm));
  if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
  const SEXP *pg = SEXPPTR(g), o = pg[6];
  SEXP res;
  int sorted = INTEGER(pg[5])[1], ng = INTEGER(pg[0])[0], *pgs = INTEGER(pg[2]), *po, *pst, l = length(x);
  if(l != length(pg[1])) error("length(g) must match length(x)");
  if(sorted) {
    po = pst = &l;
  } else if(isNull(o)) {
    int *count = (int *) Calloc(ng+1, int), *pgv = INTEGER(pg[1]);
    int *cgs = (int *) R_alloc(ng+2, sizeof(int)); cgs[1] = 1;
    po = (int *) R_alloc(l, sizeof(int)); --po;
    for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
    for(int i = 0; i != l; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
    pst = cgs + 1; ++po;
    Free(count);
  } else {
    po = INTEGER(o);
    pst = INTEGER(getAttrib(o, install("starts")));
  }

  PROTECT(res = ndistinct_g_impl(x, ng, pgs, po, pst, sorted, asLogical(Rnarm), asInteger(Rnthreads)));
  if(!isObject(x)) copyMostAttrib(x, res);
  else {
    SEXP sym_label = install("label");
    setAttrib(res, sym_label, getAttrib(x, sym_label));
  }
  UNPROTECT(1);
  return res;
}

SEXP fndistinctlC(SEXP x, SEXP g, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads) {
  int l = length(x), narm = asLogical(Rnarm), nthreads = asInteger(Rnthreads);
  if(l < 1) return x;
  if(isNull(g) && asLogical(Rdrop)) {
    SEXP out = PROTECT(allocVector(INTSXP, l)), *px = SEXPPTR(x);
    int *pout = INTEGER(out);
    if(nthreads > l) nthreads = l;
    #pragma omp parallel for num_threads(nthreads)
    for(int j = 0; j < l; ++j) pout[j] = INTEGER(ndistinct_impl(px[j], narm))[0];
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return out;
  } else {
    SEXP out = PROTECT(allocVector(VECSXP, l)), sym_label = PROTECT(install("label")),
      *pout = SEXPPTR(out), *px = SEXPPTR(x);
    if(isNull(g)) {
      if(nthreads > l) nthreads = l;
      #pragma omp parallel for num_threads(nthreads)
      for(int j = 0; j < l; ++j) {
        SEXP xj = px[j];
        pout[j] = ndistinct_impl(xj, narm);
        if(!isObject(xj)) copyMostAttrib(xj, pout[j]);
        else setAttrib(pout[j], sym_label, getAttrib(xj, sym_label));
      }
      DFcopyAttr(out, x, /*ng=*/0);
    } else {
      if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
      const SEXP *pg = SEXPPTR(g), o = pg[6];
      int sorted = INTEGER(pg[5])[1], ng = INTEGER(pg[0])[0], *pgs = INTEGER(pg[2]), *po, *pst, gl = length(pg[1]);

      if(sorted) {
        po = pst = &l;
      } else if(isNull(o)) {
        int *count = (int *) Calloc(ng+1, int), *pgv = INTEGER(pg[1]);
        int *cgs = (int *) R_alloc(ng+2, sizeof(int)); cgs[1] = 1;
        po = (int *) R_alloc(gl, sizeof(int)); --po;
        for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
        for(int i = 0; i != gl; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
        pst = cgs + 1; ++po;
        Free(count);
      } else {
        po = INTEGER(o);
        pst = INTEGER(getAttrib(o, install("starts")));
      }
      for(int j = 0; j != l; ++j) {
        SEXP xj = px[j];
        if(length(xj) != gl) error("length(g) must match nrow(x)");
        pout[j] = ndistinct_g_impl(xj, ng, pgs, po, pst, sorted, narm, nthreads);
        if(!isObject(xj)) copyMostAttrib(xj, pout[j]);
        else setAttrib(pout[j], sym_label, getAttrib(xj, sym_label));
      }
      DFcopyAttr(out, x, ng);
    }
    UNPROTECT(1);
    return out;
  }
}

SEXP fndistinctmC(SEXP x, SEXP g, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1],
      narm = asLogical(Rnarm), nthreads = asInteger(Rnthreads);
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  if(isNull(g)) {
    SEXP res = PROTECT(allocVector(INTSXP, col));
    int *pres = INTEGER(res);
    if(nthreads > col) nthreads = col;

    switch(tx) {
      case REALSXP: {
        double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j)
          pres[j] = ndistinct_double(px + j*l, &l, l, 1, narm);
        break;
      }
      case INTSXP: {  // Factor matrix not well defined object...
        int *px = INTEGER(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j)
          pres[j] = ndistinct_int(px + j*l, &l, l, 1, narm);
        break;
      }
      case LGLSXP: {
        int *px = INTEGER(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j)
          pres[j] = ndistinct_logi(px + j*l, &l, l, 1, narm);
        break;
      }
      case STRSXP: {
        SEXP *px = STRING_PTR(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j)
          pres[j] = ndistinct_string(px + j*l, &l, l, 1, narm);
        break;
      }
      default: error("Not Supported SEXP Type!");
    }
    matCopyAttr(res, x, Rdrop, /*ng=*/0);
    UNPROTECT(1);
    return res;
  } else { // With groups
    if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
    const SEXP *pg = SEXPPTR(g), o = pg[6];
    int sorted = INTEGER(pg[5])[1], ng = INTEGER(pg[0])[0], *pgs = INTEGER(pg[2]), *po, *pst, gl = length(pg[1]);
    if(l != gl) error("length(g) must match nrow(x)");

    SEXP res = PROTECT(allocVector(INTSXP, col * ng));
    int *pres = INTEGER(res);
    if(nthreads > col) nthreads = col; // column-level sufficient? or do sub-column level??

    if(sorted) {
      po = pst = &l;
      // Computing beforehand: does not seem faster?!
      // if(!isNull(o)) pst = INTEGER(getAttrib(o, install("starts")));
      // else {
      //   int *cgs = (int *) R_alloc(ng, sizeof(int)); cgs[0] = 1;
      //   for(int i = 1; i != ng; ++i) cgs[i] += pgs[i-1];
      //   pst = cgs;
      // }
    } else if(isNull(o)) {
      int *count = (int *) Calloc(ng+1, int), *pgv = INTEGER(pg[1]);
      int *cgs = (int *) R_alloc(ng+2, sizeof(int)); cgs[1] = 0;
      po = (int *) R_alloc(gl, sizeof(int)); --po;
      for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
      for(int i = 0; i != gl; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
      pst = cgs + 1; ++po;
      Free(count);
    } else {
      po = INTEGER(o);
      pst = INTEGER(getAttrib(o, install("starts")));
    }

    if(sorted) { // Sorted: compute cumulative group size (= starts) on the fly...
      int gsgr = 0;
      switch(TYPEOF(x)) {
        case REALSXP: {
          double *px = REAL(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            double *pxj = px + j * l;
            for(int gr = 0, gs = 0; gr < ng; ++gr) {
              gsgr = pgs[gr];
              pres[jng + gr] = gsgr == 0 ? 0 : ndistinct_double(pxj + gs, po, gsgr, 1, narm);
              gs += gsgr;
            }
          }
          break;
        }
        case INTSXP: {
          int *px = INTEGER(x);
          if(isFactor(x) && (ng == 0 || nlevels(x) < l / ng * 3)) {
            int M = nlevels(x);
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) {
              int *pxj = px + j * l, jng = j * ng;
              for(int gr = 0, gs = 0; gr < ng; ++gr) {
                gsgr = pgs[gr];
                pres[jng + gr] = gsgr == 0 ? 0 : ndistinct_fct(pxj + gs, po, gsgr, M, 1, narm);
                gs += gsgr;
              }
            }
          } else {
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) {
              int *pxj = px + j * l, jng = j * ng;
              for(int gr = 0, gs = 0; gr < ng; ++gr) {
                gsgr = pgs[gr];
                pres[jng + gr] = gsgr == 0 ? 0 : ndistinct_int(pxj + gs, po, gsgr, 1, narm);
                gs += gsgr;
              }
            }
          }
          break;
        }
        case LGLSXP: {
          int *px = LOGICAL(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int *pxj = px + j * l, jng = j * ng;
            for(int gr = 0, gs = 0; gr < ng; ++gr) {
              gsgr = pgs[gr];
              pres[jng + gr] = gsgr == 0 ? 0 : ndistinct_logi(pxj + gs, po, gsgr, 1, narm);
              gs += gsgr;
            }
          }
          break;
        }
        case STRSXP: {
          SEXP *px = STRING_PTR(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            SEXP *pxj = px + j * l;
            for(int gr = 0, gs = 0; gr < ng; ++gr) {
              gsgr = pgs[gr];
              pres[jng + gr] = gsgr == 0 ? 0 : ndistinct_string(pxj + gs, po, gsgr, 1, narm);
              gs += gsgr;
            }
          }
          break;
        }
        default: error("Not Supported SEXP Type!");
      }
    } else { // Not sorted. Perhaps reordering x is faster??
             // Todo: perhaps going first by groups, then by columns is better? saves zero group size checks...
      switch(TYPEOF(x)) {
        case REALSXP: {
          double *px = REAL(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            double *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr)
              pres[jng + gr] = pgs[gr] == 0 ? 0 : ndistinct_double(pxj, po + pst[gr]-1, pgs[gr], 0, narm);
          }
          break;
        }
        case INTSXP: {
          int *px = INTEGER(x);
          if(isFactor(x) && (ng == 0 || nlevels(x) < l / ng * 3)) {
            int M = nlevels(x);
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) {
              int jng = j * ng, *pxj = px + j * l;
              for(int gr = 0; gr < ng; ++gr)
                pres[jng + gr] = pgs[gr] == 0 ? 0 : ndistinct_fct(pxj, po + pst[gr]-1, pgs[gr], M, 0, narm);
            }
          } else {
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) {
              int jng = j * ng, *pxj = px + j * l;
              for(int gr = 0; gr < ng; ++gr)
                pres[jng + gr] = pgs[gr] == 0 ? 0 : ndistinct_int(pxj, po + pst[gr]-1, pgs[gr], 0, narm);
            }
          }
          break;
        }
        case LGLSXP: {
          int *px = LOGICAL(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng, *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr)
              pres[jng + gr] = pgs[gr] == 0 ? 0 : ndistinct_logi(pxj, po + pst[gr]-1, pgs[gr], 0, narm);
          }
          break;
        }
        case STRSXP: {
          SEXP *px = STRING_PTR(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            SEXP *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr)
              pres[jng + gr] = pgs[gr] == 0 ? 0 : ndistinct_string(pxj, po + pst[gr]-1, pgs[gr], 0, narm);
          }
          break;
        }
        default: error("Not Supported SEXP Type!");
      }
    }
    matCopyAttr(res, x, Rdrop, ng);
    UNPROTECT(1);
    return res;
  }
}

