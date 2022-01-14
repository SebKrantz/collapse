/*
 This code is adapted from the kit package: https://github.com/2005m/kit
 and licensed under a GPL-3.0 license.
*/

#include "kit.h"

// TODO: Check if division hash is not faster, or use Rcpp IndexHash
// TODO: Option to Preserva NA's ?

// ****************************************
// This function groups a single vector
// ****************************************
SEXP dupVecIndex(SEXP x) {
  const int n = length(x);
  int K, tx = TYPEOF(x);
  size_t M;
  if(n >= INT_MAX) error("Length of 'x' is too large. (Long vector not supported yet)"); // 1073741824
  if (tx == STRSXP || tx == REALSXP || tx == CPLXSXP ) {
    const size_t n2 = 2U * (size_t) n;
    M = 256;
    K = 8;
    while (M < n2) {
      M *= 2;
      K++;
    }
  } else if(tx == INTSXP) {
    if(isFactor(x)) {
      tx = 1000;
      M = (size_t)nlevels(x) + 2;
    } else M = (size_t)n;
  } else if (tx == LGLSXP) {
    M = 3;
  } else error("Type %s is not supported.", type2char(tx)); // # nocov
  int *h = (int*)Calloc(M, int); // Table to save the hash values, table has size M
  // memset(h, 0, M * sizeof(int)); // not needed??
  SEXP ans_i = PROTECT(allocVector(INTSXP, n));
  int *pans_i = INTEGER(ans_i);
  size_t id = 0, g = 0;
  switch (tx) {
  case LGLSXP:
  case 1000: // This is for factors or logical vectors where the size of the table is known
  {
    const int *px = INTEGER(x);
    for (int i = 0, j, k = (int)M-1; i != n; ++i) {
      j = (px[i] == NA_INTEGER) ? k : px[i];
      if(h[j]) pans_i[i] = h[j];
      else pans_i[i] = h[j] = ++g;
    }
  } break;
  case INTSXP: { // Faster version based on division hash...
    const int *px = INTEGER(x);
    unsigned int iid = 0, nu = (unsigned)n;
    for (int i = 0; i != n; ++i) {
      iid = (unsigned)px[i];
      if(iid >= nu) iid %= nu;
      // iid = (xi < nu) ? xi : xi % nu; // HASH(px[i], K); // get the hash value of x[i]
      while(h[iid]) { // Check if this hash value has been seen before
        if(px[h[iid]-1] == px[i]) { // Get the element of x that produced his value. if x[i] is the same, assign it the same index.
          pans_i[i] = pans_i[h[iid]-1]; // h[id];
          goto ibl;
        } // else, we move forward to the next slot, until we find an empty one... We need to keep checking against the values,
          // because if we found the same value before, we would also have put it in another slot after the initial one with the same hash value.
        if(++iid >= nu) iid %= nu; // # nocov
      } // We put the index into the empty slot.
      h[iid] = i + 1; // need + 1 because for zero the while loop gives false..
      pans_i[i] = ++g; // h[id];
      ibl:;
    }
  } break;
  case REALSXP: {
    const double *px = REAL(x);
    union uno tpv;
    for (int i = 0; i != n; ++i) {
      tpv.d = R_IsNA(px[i]) ? NA_REAL : (R_IsNaN(px[i]) ? R_NaN : px[i]);
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        if(REQUAL(px[h[id]-1], px[i])) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto rbl;
        }
        if(++id >= M) id %= M;
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      rbl:;
    }
  } break;
  case CPLXSXP: {
    const Rcomplex *px = COMPLEX(x);
    unsigned int u;
    union uno tpv;
    Rcomplex tmp;
    for (int i = 0; i != n; ++i) {
      tmp.r = (px[i].r == 0.0) ? 0.0 : px[i].r;
      tmp.i = (px[i].i == 0.0) ? 0.0 : px[i].i;
      if (C_IsNA(tmp)) {
        tmp.r = tmp.i = NA_REAL;
      } else if (C_IsNaN(tmp)) {
        tmp.r = tmp.i = R_NaN;
      }
      tpv.d = tmp.r;
      u = tpv.u[0] ^ tpv.u[1];
      tpv.d = tmp.i;
      u ^= tpv.u[0] ^ tpv.u[1];
      id = HASH(u, K);
      while(h[id]) {
        if(CEQUAL(px[h[id]-1], px[i])) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto cbl;
        }
        if(++id >= M) id %= M; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      cbl:;
    }
  } break;
  case STRSXP: {
    const SEXP *px = STRING_PTR(x);
    for (int i = 0; i != n; ++i) {
      id = HASH(((intptr_t) px[i] & 0xffffffff), K);
      while(h[id]) {
        if(px[h[id]-1] == px[i]) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto sbl;
        }
        if(++id >= M) id %= M; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      sbl:;
    }
  } break;
  }
  Free(h);
  setAttrib(ans_i, install("N.groups"), ScalarInteger(g));
  UNPROTECT(1);
  return ans_i;
}


SEXP dupVecIndexKeepNA(SEXP x) {
  const int n = length(x);
  int K, tx = TYPEOF(x);
  size_t M;
  if(n >= INT_MAX) error("Length of 'x' is too large. (Long vector not supported yet)"); // 1073741824
  if (tx == STRSXP || tx == REALSXP || tx == CPLXSXP ) {
    const size_t n2 = 2U * (size_t) n;
    M = 256;
    K = 8;
    while (M < n2) {
      M *= 2;
      K++;
    }
  } else if(tx == INTSXP) {
    if(isFactor(x)) {
      tx = 1000;
      M = (size_t)nlevels(x) + 2;
    } else M = (size_t)n;
  } else if (tx == LGLSXP) {
    M = 3;
  } else error("Type %s is not supported.", type2char(tx)); // # nocov
  int *h = (int*)Calloc(M, int); // Table to save the hash values, table has size M
  // memset(h, 0, M * sizeof(int)); // not needed??
  SEXP ans_i = PROTECT(allocVector(INTSXP, n));
  int *pans_i = INTEGER(ans_i);
  size_t id = 0, g = 0;
  switch (tx) {
  case LGLSXP:
  case 1000: // This is for factors or logical vectors where the size of the table is known
  {
    const int *px = INTEGER(x);
    for (int i = 0, j; i != n; ++i) {
      if(px[i] == NA_INTEGER) {
        pans_i[i] = NA_INTEGER;
        continue;
      }
      j = px[i];
      if(h[j]) pans_i[i] = h[j];
      else pans_i[i] = h[j] = ++g;
    }
  } break;
  case INTSXP: { // Faster version based on division hash...
    const int *px = INTEGER(x);
    unsigned int iid = 0, nu = (unsigned)n;
    for (int i = 0; i != n; ++i) {
      if(px[i] == NA_INTEGER) {
        pans_i[i] = NA_INTEGER;
        continue;
      }
      iid = (unsigned)px[i];
      if(iid >= nu) iid %= nu; // iid = (px[i] < n) ? px[i] : px[i] % n; // HASH(px[i], K); // get the hash value of x[i]
      while(h[iid]) { // Check if this hash value has been seen before
        if(px[h[iid]-1] == px[i]) { // Get the element of x that produced his value. if x[i] is the same, assign it the same index.
          pans_i[i] = pans_i[h[iid]-1]; // h[id];
          goto ibl;
        } // else, we move forward to the next slot, until we find an empty one... We need to keep checking against the values,
        // because if we found the same value before, we would also have put it in another slot after the initial one with the same hash value.
        if(++id >= nu) id %= nu; // ++iid; iid %= nu; // # nocov
      } // We put the index into the empty slot.
      h[iid] = i + 1; // need + 1 because for zero the while loop gives false..
      pans_i[i] = ++g; // h[id];
      ibl:;
    }
  } break;
  case REALSXP: {
    const double *px = REAL(x);
    union uno tpv;
    for (int i = 0; i != n; ++i) {
      if(ISNAN(px[i])) {
        pans_i[i] = NA_INTEGER;
        continue;
      }
      tpv.d = px[i];
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        if(REQUAL(px[h[id]-1], px[i])) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto rbl;
        }
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      rbl:;
    }
  } break;
  case CPLXSXP: {
    const Rcomplex *px = COMPLEX(x);
    unsigned int u;
    union uno tpv;
    Rcomplex tmp;
    for (int i = 0; i != n; ++i) {
      tmp.r = (px[i].r == 0.0) ? 0.0 : px[i].r;
      tmp.i = (px[i].i == 0.0) ? 0.0 : px[i].i;
      if(C_IsNA(tmp) || C_IsNaN(tmp)) {
        pans_i[i] = NA_INTEGER;
        continue;
      }
      tpv.d = tmp.r;
      u = tpv.u[0] ^ tpv.u[1];
      tpv.d = tmp.i;
      u ^= tpv.u[0] ^ tpv.u[1];
      id = HASH(u, K);
      while(h[id]) {
        if(CEQUAL(px[h[id]-1], px[i])) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto cbl;
        }
        if(++id >= M) id %= M; //++id; id %= M; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      cbl:;
    }
  } break;
  case STRSXP: {
    const SEXP *px = STRING_PTR(x);
    for (int i = 0; i != n; ++i) {
      if(px[i] == NA_STRING) {
        pans_i[i] = NA_INTEGER;
        continue;
      }
      id = HASH(((intptr_t) px[i] & 0xffffffff), K);
      while(h[id]) {
        if(px[h[id]-1] == px[i]) {
          pans_i[i] = pans_i[h[id]-1]; // h[id];
          goto sbl;
        }
        if(++id >= M) id %= M; //++id; id %= M;
      }
      h[id] = i + 1;
      pans_i[i] = ++g;
      sbl:;
    }
  } break;
  }
  Free(h);
  setAttrib(ans_i, install("N.groups"), ScalarInteger(g));
  UNPROTECT(1);
  return ans_i;
}

// TODO: Only one M calculation
// Think: If in the second grouping variable all entries are the same, you loop through the whole table for each value..
// TODO: Speed up for real values, i.e. system.time(group(DHSBR[1:2])) and system.time(group(wlddev)) (date), especially repeated real values appear slow !!
// --> But also integers is slow, i.e. system.time(group(DHSBR[1:2])) when DHSBR[2] is integer.

// **************************************************
// This function adds a second vector to the grouping
// **************************************************
int dupVecSecond(int *pidx, int *pans_i, SEXP x, const int n, const int ng) {

  if(length(x) != n) error("Unequal length columns");
  int K, tx = TYPEOF(x);
  size_t M;
  if (tx == INTSXP || tx == STRSXP || tx == REALSXP || tx == CPLXSXP ) {
    if(tx == INTSXP && isFactor(x) && (nlevels(x)+1) * ng <= 3 * n) {
      tx = 1000;
      M = (size_t)(nlevels(x)+1) * ng + 1;
    } else {
      const size_t n2 = 2U * (size_t) n;
      M = 256;
      K = 8;
      while (M < n2) {
        M *= 2;
        K++;
      }
      M += ng; // Here we addd the number of previous groups...
    }
  } else if (tx == LGLSXP) {
    M = (size_t)ng * 3 + 1;
  } else error("Type %s is not supported.", type2char(tx)); // # nocov
  int *h = (int*)Calloc(M, int); // Table to save the hash values, table has size M
  size_t id = 0, g = 0, hid = 0;
  switch (tx) {
  case LGLSXP:
  {
    const int *px = LOGICAL(x);
    for (int i = 0, j; i != n; ++i) {
      j = (px[i] == NA_LOGICAL) ? pidx[i] : pidx[i] + (px[i] + 1) * ng;
      if(h[j]) pans_i[i] = h[j];
      else pans_i[i] = h[j] = ++g;
    }
  } break;
  case 1000: // This is for factors if feasible...
  {
    const int *px = INTEGER(x);
    for (int i = 0, j; i != n; ++i) {
      j = (px[i] == NA_INTEGER) ? pidx[i] : pidx[i] + px[i] * ng;
      if(h[j]) pans_i[i] = h[j];
      else pans_i[i] = h[j] = ++g;
    }
  } break;
  // TODO: Think further about this! Perhaps you can also do this totally differently with a second vector capturing the unique values of idx!
  // See again what Morgan does to his matrix of single groupings...

  // Note: Combining bitwise i.e. px[i] ^ pidx[i] in all these functions seems slightly faster than multiplying (px[i] * pidx[i]) !
  case INTSXP: {
    const int *px = INTEGER(x);
    for (int i = 0; i != n; ++i) {
      id = (px[i] == NA_INTEGER) ? pidx[i] : HASH(px[i] ^ pidx[i], K) + pidx[i]; // This seems to be very fast (DHSBR test...)
      while(h[id]) {
        hid = h[id]-1;
        if(px[hid] == px[i] && pidx[hid] == pidx[i]) {
          pans_i[i] = pans_i[hid];
          goto ibl;
        }
        if(++id >= M) id %= M; // ++id; id %= M; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      ibl:;
    }
  } break;
  case REALSXP: {
    const double *px = REAL(x);
    union uno tpv;
    for (int i = 0; i != n; ++i) {
      tpv.d = R_IsNA(px[i]) ? NA_REAL : (R_IsNaN(px[i]) ? R_NaN :px[i]);
      id = HASH((tpv.u[0] + tpv.u[1]) ^ pidx[i], K) + pidx[i]; // Note: This is much faster than just adding pidx[i] to the hash value...
      while(h[id]) { // Problem: This value might be seen before, but not in combination with that pidx value...
        hid = h[id]-1; // The issue here is that REQUAL(px[hid], px[i]) could be true but pidx[hid] == pidx[i] fails, although the same combination of px and pidx could be seen earlier before...
        if(REQUAL(px[hid], px[i]) && pidx[hid] == pidx[i]) {
          pans_i[i] = pans_i[hid];
          goto rbl;
        }
        if(++id >= M) id %= M; //++id; id %= M;
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      rbl:;
    }
  } break;
  case CPLXSXP: {
    const Rcomplex *px = COMPLEX(x);
    unsigned int u;
    union uno tpv;
    Rcomplex tmp;
    for (int i = 0; i != n; ++i) {
      tmp.r = (px[i].r == 0.0) ? 0.0 : px[i].r;
      tmp.i = (px[i].i == 0.0) ? 0.0 : px[i].i;
      if (C_IsNA(tmp)) {
        tmp.r = tmp.i = NA_REAL;
      } else if (C_IsNaN(tmp)) {
        tmp.r = tmp.i = R_NaN;
      }
      tpv.d = tmp.r;
      u = tpv.u[0] ^ tpv.u[1];
      tpv.d = tmp.i;
      u ^= tpv.u[0] ^ tpv.u[1];
      id = HASH(u ^ pidx[i], K) + pidx[i];
      while(h[id]) {
        hid = h[id]-1;
        if(CEQUAL(px[hid], px[i]) && pidx[hid] == pidx[i]) {
          pans_i[i] = pans_i[hid];
          goto cbl;
        }
        if(++id >= M) id %= M; //++id; id %= M; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      cbl:;
    }
  } break;
  case STRSXP: {
    const SEXP *px = STRING_PTR(x);
    for (int i = 0; i != n; ++i) {
      id = HASH(((intptr_t) px[i] & 0xffffffff) ^ pidx[i], K) + pidx[i];
      while(h[id]) {
        hid = h[id]-1;
        if(px[hid] == px[i] && pidx[hid] == pidx[i]) {
          pans_i[i] = pans_i[hid];
          goto sbl;
        }
        if(++id >= M) id %= M; //++id; id %= M; // # nocov
      }
      h[id] = i + 1;
      pans_i[i] = ++g; // h[id];
      sbl:;
    }
  } break;
  }
  Free(h);
  // setAttrib(ans_i, sym_ng, ScalarInteger(g));
  // UNPROTECT(1);
  // return ans_i;
  return g;
}

// ************************************************************************
// This function brings everything together for vectors or lists of vectors
// ************************************************************************
SEXP groupVec(SEXP X, SEXP starts, SEXP sizes) {

  int l = length(X), islist = TYPEOF(X) == VECSXP,
    start = asLogical(starts), size = asLogical(sizes), nprotect = 0;
  // Better not exceptions to fundamental algorithms, when a couple of user-level functions return qG objects...
  // if(islist == 0 && OBJECT(X) != 0 && inherits(X, "qG") && inherits(X, "na.included")) return X; // return "qG" objects
  SEXP idx = islist ? dupVecIndex(VECTOR_ELT(X, 0)) : dupVecIndex(X);
  if(!(islist && l > 1) && start == 0 && size == 0) return idx; // l == 1 &&
  PROTECT(idx); ++nprotect;
  SEXP sym_ng = PROTECT(install("N.groups")), res; ++nprotect;
  int ng = asInteger(getAttrib(idx, sym_ng)), n = length(idx);
  if(islist && l > 1) {
    SEXP ans = PROTECT(allocVector(INTSXP, n)); ++nprotect;
    int i = 1, *pidx = INTEGER(idx), *pans = INTEGER(ans);
    for( ; i < l; ++i) {
      if(ng == n) break;
      if(i % 2) {
        ng = dupVecSecond(pidx, pans, VECTOR_ELT(X, i), n, ng);
      } else {
        ng = dupVecSecond(pans, pidx, VECTOR_ELT(X, i), n, ng);
      }
    }
    res = i % 2 ? idx : ans;
    setAttrib(res, sym_ng, ScalarInteger(ng));
  } else res = idx;
  // Cumpoting group starts and sizes attributes
  if(start || size) {
    PROTECT(res); ++nprotect;
    int *pres = INTEGER(res);
    if(start && size) { // Protect res ??
      SEXP gs, st;
      setAttrib(res, install("starts"), st = allocVector(INTSXP, ng));
      setAttrib(res, install("group.sizes"), gs = allocVector(INTSXP, ng));
      int *pgs = INTEGER(gs), *pst = INTEGER(st);
      memset(pgs, 0, sizeof(int) * ng); --pgs;
      memset(pst, 0, sizeof(int) * ng); --pst;
      for(int i = 0; i != n; ++i) {
        ++pgs[pres[i]];
        if(pst[pres[i]] == 0) pst[pres[i]] = i + 1;
      }
    } else if(start) {
      SEXP st;
      setAttrib(res, install("starts"), st = allocVector(INTSXP, ng));
      int *pst = INTEGER(st), k = 0;
      memset(pst, 0, sizeof(int) * ng); --pst;
      for(int i = 0; i != n; ++i) {
        if(pst[pres[i]] == 0) {
          pst[pres[i]] = i + 1;
          if(++k == ng) break;
        }
      }
    } else {
      SEXP gs;
      setAttrib(res, install("group.sizes"), gs = allocVector(INTSXP, ng));
      int *pgs = INTEGER(gs);
      memset(pgs, 0, sizeof(int) * ng); --pgs;
      for(int i = 0; i != n; ++i) ++pgs[pres[i]];
    }
  }
  UNPROTECT(nprotect);
  return res;
}


// This version is only for atomic vectors (factor generation)
SEXP groupAtVec(SEXP X, SEXP starts, SEXP naincl) {

  int start = asLogical(starts), nain = asLogical(naincl);
  // Note: These functions will give errors for unsupported types...
  SEXP idx = nain ? dupVecIndex(X) : dupVecIndexKeepNA(X);
  if(start == 0) return idx;
  PROTECT(idx);
  SEXP sym_ng = PROTECT(install("N.groups"));
  int ng = asInteger(getAttrib(idx, sym_ng)), n = length(idx);
  int *pidx = INTEGER(idx);
  SEXP st;
  setAttrib(idx, install("starts"), st = allocVector(INTSXP, ng));
  int *pst = INTEGER(st), k = 0;
  memset(pst, 0, sizeof(int) * ng); --pst;
  if(nain) {
    for(int i = 0; i != n; ++i) {
      if(pst[pidx[i]] == 0) {
        pst[pidx[i]] = i + 1;
        if(++k == ng) break;
      }
    }
  } else {
    for(int i = 0; i != n; ++i) {
      if(pidx[i] != NA_INTEGER && pst[pidx[i]] == 0) {
        pst[pidx[i]] = i + 1;
        if(++k == ng) break;
      }
    }
  }
  UNPROTECT(2);
  return idx;
}



// From the kit package...

/*
 *  Data.Frame
 */

// SEXP dupDataFrameR(SEXP x) { // move to matrix if possible
//
//   const SEXP *restrict px = SEXPPTR_RO(x);
//   const R_xlen_t len_x = xlength(x);
//   const R_xlen_t len_i = xlength(px[0]);
//   SEXP ans = R_NilValue;
//   SEXP mlv = PROTECT(allocMatrix(INTSXP, (int)len_i, (int)len_x));
//   for (R_xlen_t i = 0; i < len_x; ++i) {
//     memcpy(INTEGER(mlv)+i*len_i, INTEGER(PROTECT(dupVecIndexOnlyR(px[i]))), (unsigned)len_i*sizeof(int));
//   }
//   UNPROTECT((int)len_x);
//   const size_t n2 = 2U * (size_t) len_i;
//   size_t M = 256;
//   int K = 8;
//   while (M < n2) {
//     M *= 2;
//     K++;
//   }
//   R_xlen_t count = 0;
//   int *restrict h = (int*) Calloc(M, int);
//   const int *restrict v = INTEGER(mlv);
//   int *restrict pans = (int*) Calloc(len_i, int);
//   size_t id = 0;
//
//       for (R_xlen_t i = 0; i < len_i; ++i) {
//         R_xlen_t key = 0;
//         for (R_xlen_t j = 0; j < len_x; ++j) {
//           key ^= HASH(v[i+j*len_i],K)*97*(j+1);
//         }
//         id = HASH(key, K);
//         while (h[id]) {
//           for (R_xlen_t j = 0; j < len_x; ++j) {
//             if (v[h[id]-1+j*len_i] != v[i+j*len_i]) {
//               goto label1;
//             }
//           }
//           goto label2;
//           label1:;
//           id++; id %= M;
//         }
//         h[id] = (int) i + 1;
//         pans[i]++;
//         count++;
//         label2:;
//       }
//     Free(h);
//     UNPROTECT(1);
//     SEXP indx = PROTECT(allocVector(INTSXP, count));
//     int ct = 0;
//     int *restrict py = INTEGER(indx);
//     for (int i = 0; ct < count; ++i) {
//       if (pans[i]) {
//         py[ct++] = i;
//       }
//     }
//     SEXP output = PROTECT(subSetRowDataFrame(x, indx));
//     Free(pans);
//     UNPROTECT(2);
//     return output;
// }


/*
 *  Data.Frame
 */

// SEXP dupLenDataFrameR(SEXP x) {
//   const SEXP *restrict px = SEXPPTR_RO(x);
//   const R_xlen_t len_x = xlength(x);
//   // bool allT = true;
//   // const SEXPTYPE t0 = UTYPEOF(px[0]);
//   // for (int i = 1; i < len_x; ++i) {
//   //   if (UTYPEOF(px[i]) != t0) {
//   //     allT = false;
//   //     break;
//   //   }
//   // }
//   // if (allT) {
//   //   SEXP output = PROTECT(dupLenMatrixR(PROTECT(dfToMatrix(x))));
//   //   UNPROTECT(2);
//   //   return output;
//   // }
//   const R_xlen_t len_i = xlength(px[0]);
//   SEXP mlv = PROTECT(allocMatrix(INTSXP, (int)len_i, (int)len_x));
//   for (R_xlen_t i = 0; i < len_x; ++i) {
//     memcpy(INTEGER(mlv)+i*len_i, INTEGER(PROTECT(dupVecIndexOnlyR(px[i], ScalarLogical(false)))), (unsigned)len_i*sizeof(int));
//   }
//   UNPROTECT((int)len_x);
//   const size_t n2 = 2U * (size_t) len_i;
//   size_t M = 256;
//   int K = 8;
//   while (M < n2) {
//     M *= 2;
//     K++;
//   }
//   R_xlen_t count = 0;
//   int *restrict h = (int*) Calloc(M, int);
//   const int *restrict v = INTEGER(mlv);
//   size_t id = 0;
//   for (R_xlen_t i = 0; i < len_i; ++i) {
//     R_xlen_t key = 0;
//     for (R_xlen_t j = 0; j < len_x; ++j) {
//       key ^= HASH(v[i+j*len_i],K)*97*(j+1);
//     }
//     id = HASH(key, K);
//     while (h[id]) {
//       for (R_xlen_t j = 0; j < len_x; ++j) {
//         if (v[h[id]-1+j*len_i] != v[i+j*len_i]) {
//           goto label1;
//         }
//       }
//       goto label2;
//       label1:;
//       id++; id %= M;
//     }
//     h[id] = (int) i + 1;
//     count++;
//     label2:;
//   }
//   Free(h);
//   UNPROTECT(1);
//   return ScalarInteger(count);
// }
