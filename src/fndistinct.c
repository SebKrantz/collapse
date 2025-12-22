#include "collapse_c.h" // Needs to be first because includes OpenMP
#include "kit.h"

// C-implementations for different data types ----------------------------------
// TODO: outsource and memset hash table?
// Problem: does not work in parallel, each thread needs own hash table...

int ndistinct_int(const int *restrict px, const int *restrict po, const int l, const int sorted, const int narm) {
  if(l == 1) return !(narm && px[sorted ? 0 : po[0]-1] == NA_INTEGER);
  const size_t l2 = 2U * (size_t) l;
  size_t M = 256, id = 0;
  int K = 8, ndist = 0, anyNA = 0;
  while(M < l2) {
    M *= 2;
    K++;
  }
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values, table has size M

  if(sorted) {
    for (int i = 0; i != l; ++i) {
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
    for (int i = 0, xi; i != l; ++i) {
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

  R_Free(h);
  if(narm == 0) ndist += anyNA;
  return ndist;
}

int ndistinct_fct(const int *restrict px, const int *restrict po, const int l, const int nlev, const int sorted, const int narm) {
  if(l == 1) return !(narm && px[sorted ? 0 : po[0]-1] == NA_INTEGER);
  int *restrict h = (int*)R_Calloc(nlev+1, int);
  int ndist = 0, anyNA = narm; // Ensures breaking works if narm = TRUE or FALSE
  if(sorted) {
    for (int i = 0, xi; i != l; ++i) {
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
    for (int i = 0, xi; i != l; ++i) {
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
  R_Free(h);
  return ndist;
}

int ndistinct_logi(const int *restrict px, const int *restrict po, const int l, const int sorted, const int narm) {
  if(l == 1) return !(narm && px[sorted ? 0 : po[0]-1] == NA_LOGICAL);
  int seenT = 0, seenF = 0, anyNA = narm; // Ensures breaking works if narm = TRUE or FALSE
  if(sorted) {
    for (int i = 0, xi; i != l; ++i) {
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
    for (int i = 0, xi; i != l; ++i) {
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

int ndistinct_double(const double *restrict px, const int *restrict po, const int l, const int sorted, const int narm) {
  if(l == 1) return !(narm && ISNAN(px[sorted ? 0 : po[0]-1]));
  const size_t l2 = 2U * (size_t) l;
  size_t M = 256, id = 0;
  int K = 8, ndist = 0, anyNA = 0;
  while(M < l2) {
    M *= 2;
    K++;
  }
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values, table has size M
  union uno tpv;
  double xi;

  if(sorted) {
    for (int i = 0; i != l; ++i) {
      if(ISNAN(px[i])) {
        anyNA = 1;
        continue;
      }
      tpv.d = px[i] + 0.0; // to avoid -0.0 and 0.0 being different
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
    for (int i = 0; i != l; ++i) {
      xi = px[po[i]-1];
      if(ISNAN(xi)) {
        anyNA = 1;
        continue;
      }
      tpv.d = xi + 0.0;
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


  R_Free(h);
  if(narm == 0) ndist += anyNA;
  return ndist;
}

int ndistinct_string(const SEXP *restrict px, const int *restrict po, const int l, const int sorted, const int narm) {
  if(l == 1) return !(narm && px[sorted ? 0 : po[0]-1] == NA_STRING);
  const size_t l2 = 2U * (size_t) l;
  size_t M = 256, id = 0;
  int K = 8, ndist = 0, anyNA = 0;
  while(M < l2) {
    M *= 2;
    K++;
  }
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values, table has size M
  SEXP xi;

  if(sorted) {
    for (int i = 0; i != l; ++i) {
      if(px[i] == NA_STRING) {
        anyNA = 1;
        continue;
      }
      id = HASH(((uintptr_t) px[i] & 0xffffffff), K);
      while(h[id]) {
        if(px[h[id]-1] == px[i]) goto sbls;
        if(++id >= M) id %= M; //++id; id %= M;
      }
      h[id] = i + 1;
      ++ndist;
      sbls:;
    }
  } else {
    for (int i = 0; i != l; ++i) {
      xi = px[po[i]-1];
      if(xi == NA_STRING) {
        anyNA = 1;
        continue;
      }
      id = HASH(((uintptr_t) xi & 0xffffffff), K);
      while(h[id]) {
        if(px[po[h[id]-1]-1] == xi) goto sbl;
        if(++id >= M) id %= M; //++id; id %= M;
      }
      h[id] = i + 1;
      ++ndist;
      sbl:;
    }
  }

  R_Free(h);
  if(narm == 0) ndist += anyNA;
  return ndist;
}

// Implementations for R vectors -----------------------------------------------

int ndistinct_impl_int(SEXP x, int narm) {
  int l = length(x);
  if(l < 1) return 0;
  switch(TYPEOF(x)) {
    case REALSXP: return ndistinct_double(REAL(x), &l, l, 1, narm);
    case INTSXP:  // TODO: optimize for plain integer??
      return isFactor(x) ? ndistinct_fct(INTEGER(x), &l, l, nlevels(x), 1, narm) :
                           ndistinct_int(INTEGER(x), &l, l, 1, narm);
    case LGLSXP: return ndistinct_logi(LOGICAL(x), &l, l, 1, narm);
    case STRSXP: return ndistinct_string(SEXPPTR_RO(x), &l, l, 1, narm);
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }
}

SEXP ndistinct_impl(SEXP x, int narm) {
  return ScalarInteger(ndistinct_impl_int(x, narm));
}

// TODO: Optimize grouped distinct value count for logical vectors??
SEXP ndistinct_g_impl(SEXP x, const int ng, const int *restrict pgs, const int *restrict po, const int *restrict pst, const int sorted, const int narm, int nthreads) {

  SEXP res = PROTECT(allocVector(INTSXP, ng));
  int l = length(x), *restrict pres = INTEGER(res);
  if(nthreads > ng) nthreads = ng;

  if(sorted) { // Sorted: could compute cumulative group size (= starts) on the fly... but doesn't work multithreaded...
    po = &l;
    // int gs = 0, gsgr = 0; // need pst because gs += gsgr; doesn't work multithreaded...
    switch(TYPEOF(x)) {
      case REALSXP: {
        const double *px = REAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? 0 : ndistinct_double(px + pst[gr]-1, po, pgs[gr], 1, narm);
        break;
      }
      case INTSXP: {
        const int *px = INTEGER(x);
        if(isFactor(x) && nlevels(x) < l / ng * 3) {
          int M = nlevels(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? 0 : ndistinct_fct(px + pst[gr]-1, po, pgs[gr], M, 1, narm);
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? 0 : ndistinct_int(px + pst[gr]-1, po, pgs[gr], 1, narm);
        }
        break;
      }
      case LGLSXP: {
        const int *px = LOGICAL(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? 0 : ndistinct_logi(px + pst[gr]-1, po, pgs[gr], 1, narm);
        break;
      }
      case STRSXP: {
        const SEXP *px = SEXPPTR_RO(x);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? 0 : ndistinct_string(px + pst[gr]-1, po, pgs[gr], 1, narm);
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
        if(isFactor(x) && nlevels(x) < l / ng * 3) {
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
        const SEXP *px = SEXPPTR_RO(x);
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
  const SEXP *restrict pg = SEXPPTR_RO(g), o = pg[6];
  SEXP res;
  int sorted = LOGICAL(pg[5])[1] == 1, ng = INTEGER(pg[0])[0], *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst,
    l = length(x), nthreads = asInteger(Rnthreads);
  if(l != length(pg[1])) error("length(g) must match length(x)");
  if(l < 1) return ScalarInteger(0);
  if(isNull(o)) {
    int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;
    for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
    pst = cgs + 1;
    if(sorted) po = &l;
    else {
      int *restrict count = (int *) R_Calloc(ng+1, int);
      po = (int *) R_alloc(l, sizeof(int)); --po;
      for(int i = 0; i != l; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
      ++po; R_Free(count);
    }
  } else {
    po = INTEGER(o);
    pst = INTEGER(getAttrib(o, sym_starts));
  }
  if(nthreads > max_threads) nthreads = max_threads;
  PROTECT(res = ndistinct_g_impl(x, ng, pgs, po, pst, sorted, asLogical(Rnarm), nthreads));
  if(!isObject(x)) copyMostAttrib(x, res);
  else setAttrib(res, sym_label, getAttrib(x, sym_label));
  UNPROTECT(1);
  return res;
}

SEXP fndistinctlC(SEXP x, SEXP g, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads) {
  int l = length(x), narm = asLogical(Rnarm), nthreads = asInteger(Rnthreads);
  if(l < 1) return ScalarInteger(0);
  if(nthreads > max_threads) nthreads = max_threads;
  if(isNull(g) && asLogical(Rdrop)) {
    SEXP out = PROTECT(allocVector(INTSXP, l));
    const SEXP *restrict px = SEXPPTR_RO(x);
    int *restrict pout = INTEGER(out);
    if(nthreads <= 1) {
      for(int j = 0; j != l; ++j) pout[j] = ndistinct_impl_int(px[j], narm);
    } else {
      if(nthreads > l) nthreads = l;
      #pragma omp parallel for num_threads(nthreads)
      for(int j = 0; j < l; ++j) pout[j] = ndistinct_impl_int(px[j], narm);
    }
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return out;
  } else {
    SEXP out = PROTECT(allocVector(VECSXP, l)), *restrict pout = SEXPPTR(out);
    const SEXP *restrict px = SEXPPTR_RO(x);
    if(isNull(g)) {
      if(nthreads <= 1) {
        for(int j = 0; j != l; ++j) pout[j] = ndistinct_impl(px[j], narm);
      } else {
        if(nthreads > l) nthreads = l;
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < l; ++j) pout[j] = ndistinct_impl(px[j], narm);
      }
      // Not thread safe and thus taken out
      for(int j = 0; j != l; ++j) {
        SEXP xj = px[j];
        if(!isObject(xj)) copyMostAttrib(xj, pout[j]);
        else setAttrib(pout[j], sym_label, getAttrib(xj, sym_label));
      }
      DFcopyAttr(out, x, /*ng=*/0);
    } else {
      if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
      const SEXP *restrict pg = SEXPPTR_RO(g), o = pg[6];
      int sorted = LOGICAL(pg[5])[1] == 1, ng = INTEGER(pg[0])[0], *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst, gl = length(pg[1]);
      if(isNull(o)) {
        int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;
        for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
        pst = cgs + 1;
        if(sorted) po = &l;
        else {
          int *restrict count = (int *) R_Calloc(ng+1, int);
          po = (int *) R_alloc(gl, sizeof(int)); --po;
          for(int i = 0; i != gl; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
          ++po; R_Free(count);
        }
      } else {
        po = INTEGER(o);
        pst = INTEGER(getAttrib(o, sym_starts));
      }
      for(int j = 0; j != l; ++j) {
        SEXP xj = px[j];
        if(length(xj) != gl) error("length(g) must match nrow(x)");
        SET_VECTOR_ELT(out, j, ndistinct_g_impl(xj, ng, pgs, po, pst, sorted, narm, nthreads));
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
  if(l < 1) return ScalarInteger(0); // Prevents seqfault for numeric(0) #101
  if(nthreads > max_threads) nthreads = max_threads;
  if(isNull(g)) {
    SEXP res = PROTECT(allocVector(INTSXP, col));
    int *restrict pres = INTEGER(res);
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
        const SEXP *px = SEXPPTR_RO(x);
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
    const SEXP *restrict pg = SEXPPTR_RO(g), o = pg[6];
    int sorted = LOGICAL(pg[5])[1] == 1, ng = INTEGER(pg[0])[0], *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst, gl = length(pg[1]);
    if(l != gl) error("length(g) must match nrow(x)");

    SEXP res = PROTECT(allocVector(INTSXP, col * ng));
    int *restrict pres = INTEGER(res);
    if(nthreads > col) nthreads = col; // column-level sufficient? or do sub-column level??

    if(isNull(o)) {
      int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;
      for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
      pst = cgs + 1;
      if(sorted) po = &l;
      else {
        int *restrict count = (int *) R_Calloc(ng+1, int);
        po = (int *) R_alloc(l, sizeof(int)); --po;
        for(int i = 0; i != l; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
        ++po; R_Free(count);
      }
    } else {
      po = INTEGER(o);
      pst = INTEGER(getAttrib(o, sym_starts));
    }

    if(sorted) { // Sorted
      switch(TYPEOF(x)) {
        case REALSXP: {
          double *px = REAL(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            double *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr)
              pres[jng + gr] = pgs[gr] == 0 ? 0 : ndistinct_double(pxj + pst[gr]-1, po, pgs[gr], 1, narm);
          }
          break;
        }
        case INTSXP: { // Factor matrix not well defined object...
          int *px = INTEGER(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int *pxj = px + j * l, jng = j * ng;
            for(int gr = 0; gr < ng; ++gr)
              pres[jng + gr] = pgs[gr] == 0 ? 0 : ndistinct_int(pxj + pst[gr]-1, po, pgs[gr], 1, narm);
          }
          break;
        }
        case LGLSXP: {
          int *px = LOGICAL(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int *pxj = px + j * l, jng = j * ng;
            for(int gr = 0; gr < ng; ++gr)
              pres[jng + gr] = pgs[gr] == 0 ? 0 : ndistinct_logi(pxj + pst[gr]-1, po, pgs[gr], 1, narm);
          }
          break;
        }
        case STRSXP: {
          const SEXP *px = SEXPPTR_RO(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            const SEXP *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr)
              pres[jng + gr] = pgs[gr] == 0 ? 0 : ndistinct_string(pxj + pst[gr]-1, po, pgs[gr], 1, narm);
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
        case INTSXP: { // Factor matrix not well defined object...
          int *px = INTEGER(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng, *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr)
              pres[jng + gr] = pgs[gr] == 0 ? 0 : ndistinct_int(pxj, po + pst[gr]-1, pgs[gr], 0, narm);
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
          const SEXP *px = SEXPPTR_RO(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            const SEXP *pxj = px + j * l;
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

