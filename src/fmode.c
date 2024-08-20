#include "collapse_c.h" // Needs to be first because includes OpenMP
#include "kit.h"


static double NEG_INF = -1.0/0.0;

// C-implementations for different data types ----------------------------------
// TODO: outsource and memset hash table and count vector?
// Problem: does not work in parallel, each thread needs own intermediate vectors

int mode_int(const int *restrict px, const int *restrict po, const int l, const int sorted, const int narm, const int ret) {
  if(l == 1) return sorted ? px[0] : px[po[0]-1];
  const size_t l2 = 2U * (size_t) l;
  size_t M = 256, id = 0;
  int K = 8, index = 0, val, mode, max = 1, i = 0, end = l-1,
      minm = ret == 1, nfirstm = ret > 0, lastm = ret == 3;
  while(M < l2) {
    M *= 2;
    K++;
  }
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values
  int *restrict n = (int*)R_Calloc(l, int); // Table to count frequency of values

  if(sorted) {
    mode = px[0];
    if(narm) while(mode == NA_INTEGER && i < end) mode = px[++i];
    for(; i < l; ++i) {
      val = px[i];
      if(val == NA_INTEGER && narm) continue;
      id = HASH(val, K);
      while(h[id]) {
        index = h[id]-1;
        if(px[index] == val) goto ibls;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      ibls:;
      if(++n[index] >= max) {
        if(lastm || n[index] > max) {
          max = n[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  } else {
    mode = px[po[0]-1];
    if(narm) while(mode == NA_INTEGER && i < end) mode = px[po[++i]-1];
    for(; i < l; ++i) {
      val = px[po[i]-1];
      if(val == NA_INTEGER && narm) continue;
      id = HASH(val, K);
      while(h[id]) {
        index = h[id]-1;
        if(px[po[index]-1] == val) goto ibl;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      ibl:;
      if(++n[index] >= max) {
        if(lastm || n[index] > max) {
          max = n[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  }

  R_Free(h);
  R_Free(n);
  return mode;
}

int w_mode_int(const int *restrict px, const double *restrict pw, const int *restrict po, const int l, const int sorted, const int narm, const int ret) {
  if(l == 1) {
    if(sorted) return ISNAN(pw[0]) ? NA_INTEGER : px[0];
    return ISNAN(pw[po[0]-1]) ? NA_INTEGER : px[po[0]-1];
  }
  const size_t l2 = 2U * (size_t) l;
  size_t M = 256, id = 0;
  int K = 8, index = 0, val, mode, i = 0, end = l-1,
    minm = ret == 1, nfirstm = ret > 0, lastm = ret == 3;
  while(M < l2) {
    M *= 2;
    K++;
  }
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values
  double *restrict sumw = (double*)R_Calloc(l, double); // Table to save each values sum of weights
  double max = NEG_INF;

  if(sorted) {
    mode = px[0];
    if(narm) while((mode == NA_INTEGER || ISNAN(pw[i])) && i < end) mode = px[++i];
    for(; i < l; ++i) {
      val = px[i];
      if(ISNAN(pw[i]) || (val == NA_INTEGER && narm)) continue;
      id = HASH(val, K);
      while(h[id]) {
        index = h[id]-1;
        if(px[index] == val) goto ibls;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      ibls:;
      sumw[index] += pw[i];
      if(sumw[index] >= max) {
        if(lastm || sumw[index] > max) {
          max = sumw[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  } else {
    mode = px[po[0]-1];
    if(narm) while((mode == NA_INTEGER || ISNAN(pw[po[i]-1])) && i < end) mode = px[po[++i]-1];
    for(int oi; i < l; ++i) {
      oi = po[i]-1;
      val = px[oi];
      if(ISNAN(pw[oi]) || (val == NA_INTEGER && narm)) continue;
      id = HASH(val, K);
      while(h[id]) {
        index = h[id]-1;
        if(px[po[index]-1] == val) goto ibl;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      ibl:;
      sumw[index] += pw[oi];
      if(sumw[index] >= max) {
        if(lastm || sumw[index] > max) {
          max = sumw[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  }

  R_Free(h);
  R_Free(sumw);
  return mode;
}


int mode_fct_logi(const int *restrict px, const int *restrict po, const int l, const int nlev, const int sorted, const int narm, const int ret) {
  if(l == 1) return sorted ? px[0] : px[po[0]-1];
  int val, mode, max = 1, nlevp = nlev + 1, i = 0, end = l-1,
    minm = ret == 1, nfirstm = ret > 0, lastm = ret == 3;
  int *restrict n = (int*)R_Calloc(nlevp+1, int); // Table to count frequency of values

  if(sorted) {
    mode = px[0];
    if(narm) while(mode == NA_INTEGER && i < end) mode = px[++i];
    for(; i < l; ++i) {
      val = px[i];
      if(val == NA_INTEGER) {
        if(narm) continue;
        val = nlevp;
      }
      if(++n[val] >= max) {
        if(lastm || n[val] > max) {
          max = n[val];
          mode = px[i];
        } else if(nfirstm) {
          if(minm) {
            if(mode > px[i]) mode = px[i];
          } else {
            if(mode < px[i]) mode = px[i];
          }
        }
      }
    }
  } else {
    mode = px[po[0]-1];
    if(narm) while(mode == NA_INTEGER && i < end) mode = px[po[++i]-1];
    for(int xi; i < l; ++i) {
      val = xi = px[po[i]-1];
      if(val == NA_INTEGER) {
        if(narm) continue;
        val = nlevp;
      }
      if(++n[val] >= max) {
        if(lastm || n[val] > max) {
          max = n[val];
          mode = xi;
        } else if(nfirstm) {
          if(minm) {
            if(mode > xi) mode = xi;
          } else {
            if(mode < xi) mode = xi;
          }
        }
      }
    }
  }

  R_Free(n);
  return mode;
}

int w_mode_fct_logi(const int *restrict px, const double *restrict pw, const int *restrict po, const int l, const int nlev, const int sorted, const int narm, const int ret) {
  if(l == 1) {
    if(sorted) return ISNAN(pw[0]) ? NA_INTEGER : px[0];
    return ISNAN(pw[po[0]-1]) ? NA_INTEGER : px[po[0]-1];
  }
  int val, mode, nlevp = nlev + 1, i = 0, end = l-1,
    minm = ret == 1, nfirstm = ret > 0, lastm = ret == 3;
  double *restrict sumw = (double*)R_Calloc(nlevp+1, double); // Table to save each values sum of weights
  double max = NEG_INF;

  if(sorted) {
    mode = px[0];
    if(narm) while((mode == NA_INTEGER || ISNAN(pw[i])) && i < end) mode = px[++i];
    for(; i < l; ++i) {
      if(ISNAN(pw[i])) continue;
      val = px[i];
      if(val == NA_INTEGER) {
        if(narm) continue;
        val = nlevp;
      }
      sumw[val] += pw[i];
      if(sumw[val] >= max) {
        if(lastm || sumw[val] > max) {
          max = sumw[val];
          mode = px[i];
        } else if(nfirstm) {
          if(minm) {
            if(mode > px[i]) mode = px[i];
          } else {
            if(mode < px[i]) mode = px[i];
          }
        }
      }
    }
  } else {
    mode = px[po[0]-1];
    if(narm) while((mode == NA_INTEGER || ISNAN(pw[po[i]-1])) && i < end) mode = px[po[++i]-1];
    for(int oi, xoi; i < l; ++i) {
      oi = po[i]-1;
      if(ISNAN(pw[oi])) continue;
      val = xoi = px[oi];
      if(val == NA_INTEGER) {
        if(narm) continue;
        val = nlevp;
      }
      sumw[val] += pw[oi];
      if(sumw[val] >= max) {
        if(lastm || sumw[val] > max) {
          max = sumw[val];
          mode = xoi;
        } else if(nfirstm) {
          if(minm) {
            if(mode > xoi) mode = xoi;
          } else {
            if(mode < xoi) mode = xoi;
          }
        }
      }
    }
  }

  R_Free(sumw);
  return mode;
}


double mode_double(const double *restrict px, const int *restrict po, const int l, const int sorted, const int narm, const int ret) {
  if(l == 1) return sorted ? px[0] : px[po[0]-1];
  const size_t l2 = 2U * (size_t) l;
  size_t M = 256, id = 0;
  int K = 8, index = 0, max = 1, i = 0, end = l-1,
    minm = ret == 1, nfirstm = ret > 0, lastm = ret == 3;
  while(M < l2) {
    M *= 2;
    K++;
  }
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values
  int *restrict n = (int*)R_Calloc(l, int); // Table to count frequency of values
  double val, mode;
  union uno tpv;

  if(sorted) {
    mode = px[0];
    if(narm) while(ISNAN(mode) && i < end) mode = px[++i];
    for(; i < l; ++i) {
      val = px[i];
      if(ISNAN(val) && narm) continue;
      tpv.d = val;
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        index = h[id]-1;
        if(REQUAL(px[index], val)) goto rbls;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      rbls:;
      if(++n[index] >= max) {
        if(lastm || n[index] > max) {
          max = n[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  } else {
    mode = px[po[0]-1];
    if(narm) while(ISNAN(mode) && i < end) mode = px[po[++i]-1];
    for(; i < l; ++i) {
      val = px[po[i]-1];
      if(ISNAN(val) && narm) continue;
      tpv.d = val;
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        index = h[id]-1;
        if(REQUAL(px[po[index]-1], val)) goto rbl;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      rbl:;
      if(++n[index] >= max) {
        if(lastm || n[index] > max) {
          max = n[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  }

  R_Free(h);
  R_Free(n);
  return mode;
}

double w_mode_double(const double *restrict px, const double *restrict pw, const int *restrict po, const int l, const int sorted, const int narm, const int ret) {
  if(l == 1) {
    if(sorted) return ISNAN(pw[0]) ? NA_REAL : px[0];
    return ISNAN(pw[po[0]-1]) ? NA_REAL : px[po[0]-1];
  }
  const size_t l2 = 2U * (size_t) l;
  size_t M = 256, id = 0;
  int K = 8, index = 0, i = 0, end = l-1, minm = ret == 1, nfirstm = ret > 0, lastm = ret == 3;
  while(M < l2) {
    M *= 2;
    K++;
  }
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values
  double *restrict sumw = (double*)R_Calloc(l, double); // Table to save each values sum of weights
  double val, mode, max = NEG_INF;
  union uno tpv;

  if(sorted) {
    mode = px[0];
    if(narm) while((ISNAN(mode) || ISNAN(pw[i])) && i < end) mode = px[++i];
    for(; i < l; ++i) {
      val = px[i];
      if(ISNAN(pw[i]) || (ISNAN(val) && narm)) continue;
      tpv.d = val;
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        index = h[id]-1;
        if(REQUAL(px[index], val)) goto rbls;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      rbls:;
      sumw[index] += pw[i];
      if(sumw[index] >= max) {
        if(lastm || sumw[index] > max) {
          max = sumw[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  } else {
    mode = px[po[0]-1];
    if(narm) while((ISNAN(mode) || ISNAN(pw[po[i]-1])) && i < end) mode = px[po[++i]-1];
    for(int oi; i < l; ++i) {
      oi = po[i]-1;
      val = px[oi];
      if(ISNAN(pw[oi]) || (ISNAN(val) && narm)) continue;
      tpv.d = val;
      id = HASH(tpv.u[0] + tpv.u[1], K);
      while(h[id]) {
        index = h[id]-1;
        if(REQUAL(px[po[index]-1], val)) goto rbl;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      rbl:;
      sumw[index] += pw[oi];
      if(sumw[index] >= max) {
        if(lastm || sumw[index] > max) {
          max = sumw[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  }

  R_Free(h);
  R_Free(sumw);
  return mode;
}


SEXP mode_string(const SEXP *restrict px, const int *restrict po, const int l, const int sorted, const int narm, const int ret) {
  if(l == 1) return sorted ? px[0] : px[po[0]-1];
  const size_t l2 = 2U * (size_t) l;
  size_t M = 256, id = 0;
  int K = 8, index = 0, max = 1, i = 0, end = l-1,
    minm = ret == 1, nfirstm = ret > 0, lastm = ret == 3;
  while(M < l2) {
    M *= 2;
    K++;
  }
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values
  int *restrict n = (int*)R_Calloc(l, int); // Table to count frequency of values
  SEXP val, mode;

  if(sorted) {
    mode = px[0];
    if(narm) while(mode == NA_STRING && i < end) mode = px[++i];
    for(; i < l; ++i) {
      val = px[i];
      if(val == NA_STRING && narm) continue;
      id = HASH(((uintptr_t) val & 0xffffffff), K);
      while(h[id]) {
        index = h[id]-1;
        if(px[index] == val) goto sbls;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      sbls:;
      if(++n[index] >= max) {
        if(lastm || n[index] > max) {
          max = n[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  } else {
    mode = px[po[0]-1];
    if(narm) while(mode == NA_STRING && i < end) mode = px[po[++i]-1];
    for(; i < l; ++i) {
      val = px[po[i]-1];
      if(val == NA_STRING && narm) continue;
      id = HASH(((uintptr_t) val & 0xffffffff), K);
      while(h[id]) {
        index = h[id]-1;
        if(px[po[index]-1] == val) goto sbl;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      sbl:;
      if(++n[index] >= max) {
        if(lastm || n[index] > max) {
          max = n[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  }

  R_Free(h);
  R_Free(n);
  return mode;
}

SEXP w_mode_string(const SEXP *restrict px, const double *restrict pw, const int *restrict po, const int l, const int sorted, const int narm, const int ret) {
  if(l == 1) {
    if(sorted) return ISNAN(pw[0]) ? NA_STRING : px[0];
    return ISNAN(pw[po[0]-1]) ? NA_STRING : px[po[0]-1];
  }
  const size_t l2 = 2U * (size_t) l;
  size_t M = 256, id = 0;
  int K = 8, index = 0, i = 0, end = l-1, minm = ret == 1, nfirstm = ret > 0, lastm = ret == 3;
  while(M < l2) {
    M *= 2;
    K++;
  }
  int *restrict h = (int*)R_Calloc(M, int); // Table to save the hash values
  double *restrict sumw = (double*)R_Calloc(l, double); // Table to save each values sum of weights
  double max = NEG_INF;
  SEXP val, mode;

  if(sorted) {
    mode = px[0];
    if(narm) while((mode == NA_STRING || ISNAN(pw[i])) && i < end) mode = px[++i];
    for(; i != l; ++i) {
      val = px[i];
      if(ISNAN(pw[i]) || (val == NA_STRING && narm)) continue;
      id = HASH(((uintptr_t) val & 0xffffffff), K);
      while(h[id]) {
        index = h[id]-1;
        if(px[index] == val) goto sbls;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      sbls:;
      sumw[index] += pw[i];
      if(sumw[index] >= max) {
        if(lastm || sumw[index] > max) {
          max = sumw[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  } else {
    mode = px[po[0]-1];
    if(narm) while((mode == NA_STRING || ISNAN(pw[po[i]-1])) && i < end) mode = px[po[++i]-1];
    for(int oi; i != l; ++i) {
      oi = po[i]-1;
      val = px[oi];
      if(ISNAN(pw[oi]) || (val == NA_STRING && narm)) continue;
      id = HASH(((uintptr_t) val & 0xffffffff), K);
      while(h[id]) {
        index = h[id]-1;
        if(px[po[index]-1] == val) goto sbl;
        if(++id >= M) id %= M; // ++id; id %= M;
      }
      h[id] = i + 1;
      index = i; // n[i] = 1;
      sbl:;
      sumw[index] += pw[oi];
      if(sumw[index] >= max) {
        if(lastm || sumw[index] > max) {
          max = sumw[index];
          mode = val;
        } else if(nfirstm) {
          if(minm) {
            if(mode > val) mode = val;
          } else {
            if(mode < val) mode = val;
          }
        }
      }
    }
  }

  R_Free(h);
  R_Free(sumw);
  return mode;
}


// Implementations for R vectors -----------------------------------------------

// https://github.com/wch/r-source/blob/trunk/src/include/Rinlinedfuns.h
SEXP my_ScalarLogical(int x) {
  SEXP ans = allocVector(LGLSXP, 1);
  // SET_SCALAR_LVAL(ans, x); // Not part of the API
  LOGICAL(ans)[0] = x;
  return ans;
}

// Splitting this up to increase thread safety
SEXP mode_impl_plain(SEXP x, int narm, int ret) {
  int l = length(x);
  if(l <= 1) return x;

  switch(TYPEOF(x)) {
    case REALSXP: return ScalarReal(mode_double(REAL(x), &l, l, 1, narm, ret));
    case INTSXP:  return ScalarInteger(isFactor(x) ? mode_fct_logi(INTEGER(x), &l, l, nlevels(x), 1, narm, ret) :
                                    mode_int(INTEGER(x), &l, l, 1, narm, ret));
    case LGLSXP: return my_ScalarLogical(mode_fct_logi(LOGICAL(x), &l, l, 1, 1, narm, ret));
    case STRSXP: return ScalarString(mode_string(SEXPPTR(x), &l, l, 1, narm, ret));
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }
}

SEXP mode_impl(SEXP x, int narm, int ret) {
  if(length(x) <= 1) return x;
  SEXP res = PROTECT(mode_impl_plain(x, narm, ret));
  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

SEXP w_mode_impl_plain(SEXP x, double *pw, int narm, int ret) {
  int l = length(x);
  if(l <= 1) return x;

  switch(TYPEOF(x)) {
    case REALSXP: return ScalarReal(w_mode_double(REAL(x), pw, &l, l, 1, narm, ret));
    case INTSXP:  return ScalarInteger(isFactor(x) ? w_mode_fct_logi(INTEGER(x), pw, &l, l, nlevels(x), 1, narm, ret) :
                                    w_mode_int(INTEGER(x), pw, &l, l, 1, narm, ret));
    case LGLSXP:  return my_ScalarLogical(w_mode_fct_logi(LOGICAL(x), pw, &l, l, 1, 1, narm, ret));
    case STRSXP:  return ScalarString(w_mode_string(SEXPPTR(x), pw, &l, l, 1, narm, ret));
    default: error("Not Supported SEXP Type: '%s'", type2char(TYPEOF(x)));
  }
}

SEXP w_mode_impl(SEXP x, double *pw, int narm, int ret) {
  if(length(x) <= 1) return x;
  SEXP res = PROTECT(w_mode_impl_plain(x, pw, narm, ret));
  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

SEXP mode_g_impl(SEXP x, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int ret, int nthreads) {

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
          pres[gr] = pgs[gr] == 0 ? NA_REAL : mode_double(px + pst[gr]-1, po, pgs[gr], 1, narm, ret);
        break;
      }
      case INTSXP: {
        int *px = INTEGER(x), *pres = INTEGER(res);
        if(isFactor(x) && nlevels(x) < l / ng * 3) {
          int M = nlevels(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? NA_INTEGER : mode_fct_logi(px + pst[gr]-1, po, pgs[gr], M, 1, narm, ret);
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? NA_INTEGER : mode_int(px + pst[gr]-1, po, pgs[gr], 1, narm, ret);
        }
        break;
      }
      case LGLSXP: {
        int *px = LOGICAL(x), *pres = LOGICAL(res);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? NA_LOGICAL : mode_fct_logi(px + pst[gr]-1, po, pgs[gr], 1, 1, narm, ret);
        break;
      }
      case STRSXP: {
        SEXP *px = SEXPPTR(x), *pres = SEXPPTR(res);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? NA_STRING : mode_string(px + pst[gr]-1, po, pgs[gr], 1, narm, ret);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(tx));
    }
  } else { // Not sorted. Perhaps reordering x is faster??
    switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pres = REAL(res);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? NA_REAL : mode_double(px, po + pst[gr]-1, pgs[gr], 0, narm, ret);
        break;
      }
      case INTSXP: {
        int *px = INTEGER(x), *pres = INTEGER(res);
        if(isFactor(x) && nlevels(x) < l / ng * 3) {
          int M = nlevels(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? NA_INTEGER : mode_fct_logi(px, po + pst[gr]-1, pgs[gr], M, 0, narm, ret);
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? NA_INTEGER : mode_int(px, po + pst[gr]-1, pgs[gr], 0, narm, ret);
        }
        break;
      }
      case LGLSXP: {
        int *px = LOGICAL(x), *pres = LOGICAL(res);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? NA_LOGICAL : mode_fct_logi(px, po + pst[gr]-1, pgs[gr], 1, 0, narm, ret);
        break;
      }
      case STRSXP: {
        SEXP *px = SEXPPTR(x), *pres = SEXPPTR(res);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? NA_STRING : mode_string(px, po + pst[gr]-1, pgs[gr], 0, narm, ret);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(tx));
    }
  }

  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

SEXP w_mode_g_impl(SEXP x, double *pw, int ng, int *pgs, int *po, int *pst, int sorted, int narm, int ret, int nthreads) {

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
          pres[gr] = pgs[gr] == 0 ? NA_REAL : w_mode_double(px + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret);
        break;
      }
      case INTSXP: {
        int *px = INTEGER(x), *pres = INTEGER(res);
        if(isFactor(x) && nlevels(x) < l / ng * 3) {
          int M = nlevels(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? NA_INTEGER : w_mode_fct_logi(px + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], M, 1, narm, ret);
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? NA_INTEGER : w_mode_int(px + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret);
        }
        break;
      }
      case LGLSXP: {
        int *px = LOGICAL(x), *pres = LOGICAL(res);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? NA_LOGICAL : w_mode_fct_logi(px + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, 1, narm, ret);
        break;
      }
      case STRSXP: {
        SEXP *px = SEXPPTR(x), *pres = SEXPPTR(res);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? NA_STRING : w_mode_string(px + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(tx));
    }
  } else { // Not sorted. Perhaps reordering x is faster??
    switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pres = REAL(res);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? NA_REAL : w_mode_double(px, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret);
        break;
      }
      case INTSXP: {
        int *px = INTEGER(x), *pres = INTEGER(res);
        if(isFactor(x) && nlevels(x) < l / ng * 3) {
          int M = nlevels(x);
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? NA_INTEGER : w_mode_fct_logi(px, pw, po + pst[gr]-1, pgs[gr], M, 0, narm, ret);
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int gr = 0; gr < ng; ++gr)
            pres[gr] = pgs[gr] == 0 ? NA_INTEGER : w_mode_int(px, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret);
        }
        break;
      }
      case LGLSXP: {
        int *px = LOGICAL(x), *pres = LOGICAL(res);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? NA_LOGICAL : w_mode_fct_logi(px, pw, po + pst[gr]-1, pgs[gr], 1, 0, narm, ret);
        break;
      }
      case STRSXP: {
        SEXP *px = SEXPPTR(x), *pres = SEXPPTR(res);
        #pragma omp parallel for num_threads(nthreads)
        for(int gr = 0; gr < ng; ++gr)
          pres[gr] = pgs[gr] == 0 ? NA_STRING : w_mode_string(px, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret);
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(tx));
    }
  }

  copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}


// Functions for Export --------------------------------------------------------

SEXP fmodeC(SEXP x, SEXP g, SEXP w, SEXP Rnarm, SEXP Rret, SEXP Rnthreads) {
  int nullg = isNull(g), nullw = isNull(w), l = length(x), nprotect = 0;
  if(l <= 1) return x;
  if(nullg && nullw) return mode_impl(x, asLogical(Rnarm), asInteger(Rret));
  double tmp = 0.0, *restrict pw = &tmp;
  if(!nullw) {
    if(length(w) != l) error("length(w) must match length(x)");
    if(TYPEOF(w) != REALSXP) {
      if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double)");
      w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
    }
    pw = REAL(w);
  }
  if(nullg) {
    // if(TYPEOF(w) != REALSXP)
    UNPROTECT(nprotect);
    return w_mode_impl(x, pw, asLogical(Rnarm), asInteger(Rret));
  }
  if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
  const SEXP *restrict pg = SEXPPTR_RO(g), o = pg[6];
  int sorted = LOGICAL(pg[5])[1] == 1, ng = INTEGER(pg[0])[0], *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst, nthreads = asInteger(Rnthreads);
  if(l != length(pg[1])) error("length(g) must match length(x)");
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
    pst = INTEGER(getAttrib(o, install("starts")));
  }
  // if(nullw) return mode_g_impl(x, ng, pgs, po, pst, sorted, asLogical(Rnarm), asInteger(Rret), asInteger(Rnthreads));
  // if(TYPEOF(w) != REALSXP) UNPROTECT(nprotect);
  // return w_mode_g_impl(x, pw, ng, pgs, po, pst, sorted, asLogical(Rnarm), asInteger(Rret), asInteger(Rnthreads));
  // Thomas Kalibera Patch:
  if(nthreads > max_threads) nthreads = max_threads;
  SEXP res;
  if(nullw) res = mode_g_impl(x, ng, pgs, po, pst, sorted, asLogical(Rnarm), asInteger(Rret), nthreads);
  else res = w_mode_g_impl(x, pw, ng, pgs, po, pst, sorted, asLogical(Rnarm), asInteger(Rret), nthreads);
  UNPROTECT(nprotect);
  return res;
}


// TODO: allow column-level parallelism??
SEXP fmodelC(SEXP x, SEXP g, SEXP w, SEXP Rnarm, SEXP Rret, SEXP Rnthreads) {
  int nullg = isNull(g), nullw = isNull(w), l = length(x), ng = 0, nprotect = 1,
    narm = asLogical(Rnarm), ret = asInteger(Rret), nthreads = asInteger(Rnthreads);
  if(l < 1) return x;
  if(nthreads > max_threads) nthreads = max_threads;
  SEXP out = PROTECT(allocVector(VECSXP, l)), *restrict pout = SEXPPTR(out);
  const SEXP *restrict px = SEXPPTR_RO(x);
  if(nullg && nthreads > l) nthreads = l;
  if(nullg && nullw) {
    if(nthreads <= 1) {
      for(int j = 0; j != l; ++j) pout[j] = mode_impl(px[j], narm, ret);
    } else {
      #pragma omp parallel for num_threads(nthreads)
      for(int j = 0; j < l; ++j) pout[j] = mode_impl_plain(px[j], narm, ret);
      for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]); // Not thread safe and thus taken out...
    }
  } else {
    int nrx = length(px[0]);
    double tmp = 0.0, *restrict pw = &tmp;
    if(!nullw) {
      if(length(w) != nrx) error("length(w) must match nrow(x)");
      if(TYPEOF(w) != REALSXP) {
        if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double)");
        w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
      }
      pw = REAL(w);
    }
    if(nullg) {
      if(nthreads <= 1) {
        for(int j = 0; j != l; ++j) pout[j] = w_mode_impl(px[j], pw, narm, ret);
      } else {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < l; ++j) pout[j] = w_mode_impl_plain(px[j], pw, narm, ret);
        for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]); // Not thread safe and thus taken out...
      }
    } else {
      if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
      const SEXP *restrict pg = SEXPPTR_RO(g), o = pg[6];
      ng = INTEGER(pg[0])[0];
      int sorted = LOGICAL(pg[5])[1] == 1, *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst;
      if(nrx != length(pg[1])) error("length(g) must match nrow(x)");
      if(isNull(o)) {
        int *cgs = (int *) R_alloc(ng+2, sizeof(int)), *restrict pgv = INTEGER(pg[1]); cgs[1] = 1;
        for(int i = 0; i != ng; ++i) cgs[i+2] = cgs[i+1] + pgs[i];
        pst = cgs + 1;
        if(sorted) po = &l;
        else {
          int *restrict count = (int *) R_Calloc(ng+1, int);
          po = (int *) R_alloc(nrx, sizeof(int)); --po;
          for(int i = 0; i != nrx; ++i) po[cgs[pgv[i]] + count[pgv[i]]++] = i+1;
          ++po; R_Free(count);
        }
      } else {
        po = INTEGER(o);
        pst = INTEGER(getAttrib(o, install("starts")));
      }
      if(nullw) { // Parallelism at sub-column level
        for(int j = 0; j < l; ++j) pout[j] = mode_g_impl(px[j], ng, pgs, po, pst, sorted, narm, ret, nthreads);
      } else { // Parallelism at sub-column level
        for(int j = 0; j < l; ++j) pout[j] = w_mode_g_impl(px[j], pw, ng, pgs, po, pst, sorted, narm, ret, nthreads);
      }
    }
  }
  DFcopyAttr(out, x, ng);
  UNPROTECT(nprotect);
  return out;
}


SEXP fmodemC(SEXP x, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rret, SEXP Rnthreads) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1],
      narm = asLogical(Rnarm), ret = asInteger(Rret), nthreads = asInteger(Rnthreads),
      nullg = isNull(g), nullw = isNull(w), nprotect = 1;
  if(l <= 1) return x; // Prevents seqfault for numeric(0) #101
  if(nthreads > max_threads) nthreads = max_threads;
  if(nthreads > col) nthreads = col;

  double tmp = 0.0, *restrict pw = &tmp;
  if(!nullw) {
    if(length(w) != l) error("length(w) must match nrow(x)");
    if(TYPEOF(w) != REALSXP) {
      if(!(TYPEOF(w) == INTSXP || TYPEOF(w) == LGLSXP)) error("weights need to be double or integer/logical (internally coerced to double)");
      w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
    }
    pw = REAL(w);
  }

  if(nullg) {
    SEXP res = PROTECT(allocVector(tx, col));

    switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *restrict pres = REAL(res);
        if(nullw) {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) pres[j] = mode_double(px + j*l, &l, l, 1, narm, ret);
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) pres[j] = w_mode_double(px + j*l, pw, &l, l, 1, narm, ret);
        }
        break;
      }
      case INTSXP: {  // Factor matrix not well defined object...
        int *px = INTEGER(x), *restrict pres = INTEGER(res);
        if(nullw) {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) pres[j] = mode_int(px + j*l, &l, l, 1, narm, ret);
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) pres[j] = w_mode_int(px + j*l, pw, &l, l, 1, narm, ret);
        }
        break;
      }
      case LGLSXP: {
        int *px = LOGICAL(x), *restrict pres = LOGICAL(res);
        if(nullw) {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) pres[j] = mode_fct_logi(px + j*l, &l, l, 1, 1, narm, ret);
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) pres[j] = w_mode_fct_logi(px + j*l, pw, &l, l, 1, 1, narm, ret);
        }
        break;
      }
      case STRSXP: {
        SEXP *px = SEXPPTR(x), *restrict pres = SEXPPTR(res);
        if(nullw) {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) pres[j] = mode_string(px + j*l, &l, l, 1, narm, ret);
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) pres[j] = w_mode_string(px + j*l, pw, &l, l, 1, narm, ret);
        }
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(tx));
    }

    matCopyAttr(res, x, Rdrop, /*ng=*/0);
    UNPROTECT(nprotect);
    return res;
  }

  // With groups
  if(TYPEOF(g) != VECSXP || !inherits(g, "GRP")) error("g needs to be an object of class 'GRP', see ?GRP");
  const SEXP *restrict pg = SEXPPTR_RO(g), o = pg[6];
  int sorted = LOGICAL(pg[5])[1] == 1, ng = INTEGER(pg[0])[0], *restrict pgs = INTEGER(pg[2]), *restrict po, *restrict pst, gl = length(pg[1]);
  if(l != gl) error("length(g) must match nrow(x)");
  SEXP res = PROTECT(allocVector(tx, ng * col));

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
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : mode_double(pxj + pst[gr]-1, po, pgs[gr], 1, narm, ret);
          }
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            double *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : w_mode_double(pxj + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret);
          }
        }
        break;
      }
      case INTSXP: { // Factor matrix not well defined object...
        int *px = INTEGER(x), *restrict pres = INTEGER(res);
        if(nullw) {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int *pxj = px + j * l, jng = j * ng;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : mode_int(pxj + pst[gr]-1, po, pgs[gr], 1, narm, ret);
          }
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int *pxj = px + j * l, jng = j * ng;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : w_mode_int(pxj + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret);
          }
        }
        break;
      }
      case LGLSXP: {
        int *px = LOGICAL(x), *restrict pres = LOGICAL(res);
        if(nullw) {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int *pxj = px + j * l, jng = j * ng;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_LOGICAL : mode_fct_logi(pxj + pst[gr]-1, po, pgs[gr], 1, 1, narm, ret);
          }
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int *pxj = px + j * l, jng = j * ng;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_LOGICAL : w_mode_fct_logi(pxj + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, 1, narm, ret);
          }
        }
        break;
      }
      case STRSXP: {
        SEXP *px = SEXPPTR(x), *restrict pres = SEXPPTR(res);
        if(nullw) {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            SEXP *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_STRING : mode_string(pxj + pst[gr]-1, po, pgs[gr], 1, narm, ret);
          }
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            SEXP *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_STRING : w_mode_string(pxj + pst[gr]-1, pw + pst[gr]-1, po, pgs[gr], 1, narm, ret);
          }
        }
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(tx));
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
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : mode_double(pxj, po + pst[gr]-1, pgs[gr], 0, narm, ret);
          }
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            double *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_REAL : w_mode_double(pxj, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret);
          }
        }
        break;
      }
      case INTSXP: {
        int *px = INTEGER(x), *restrict pres = INTEGER(res);
        if(nullw) {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng, *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : mode_int(pxj, po + pst[gr]-1, pgs[gr], 0, narm, ret);
          }
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng, *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_INTEGER : w_mode_int(pxj, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret);
          }
        }
        break;
      }
      case LGLSXP: {
        int *px = LOGICAL(x), *restrict pres = LOGICAL(res);
        if(nullw) {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng, *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_LOGICAL : mode_fct_logi(pxj, po + pst[gr]-1, pgs[gr], 1, 0, narm, ret);
          }
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng, *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_LOGICAL : w_mode_fct_logi(pxj, pw, po + pst[gr]-1, pgs[gr], 1, 0, narm, ret);
          }
        }
        break;
      }
      case STRSXP: {
        SEXP *px = SEXPPTR(x), *restrict pres = SEXPPTR(res);
        if(nullw) {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            SEXP *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_STRING : mode_string(pxj, po + pst[gr]-1, pgs[gr], 0, narm, ret);
          }
        } else {
          #pragma omp parallel for num_threads(nthreads)
          for(int j = 0; j < col; ++j) {
            int jng = j * ng;
            SEXP *pxj = px + j * l;
            for(int gr = 0; gr < ng; ++gr) pres[jng + gr] = pgs[gr] == 0 ? NA_STRING : w_mode_string(pxj, pw, po + pst[gr]-1, pgs[gr], 0, narm, ret);
          }
        }
        break;
      }
      default: error("Not Supported SEXP Type: '%s'", type2char(tx));
    }
  }

  matCopyAttr(res, x, Rdrop, ng);
  UNPROTECT(nprotect);
  return res;
}

