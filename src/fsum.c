#ifdef _OPENMP
#include <omp.h>
#endif
#include "collapse_c.h"
// #include <R_ext/Altrep.h>

void fsum_double_impl(double *restrict pout, const double *restrict px, const int narm, const int l) {
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
}

void fsum_double_g_impl(double *restrict pout, const double *restrict px, const int ng, const int *restrict pg, const int narm, const int l) {
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

void fsum_double_omp_impl(double *restrict pout, const double *restrict px, const int narm, const int l, const int nth) {
  double sum;
  if(narm) {
    int j = 1;
    sum = px[0];
    while(ISNAN(sum) && j != l) sum = px[j++];
    if(j != l) {
      #pragma omp parallel for num_threads(nth) reduction(+:sum)
      for(int i = j; i < l; ++i) if(NISNAN(px[i])) sum += px[i]; // Fastest ?
    }
  } else {
    sum = 0;
    #pragma omp parallel for num_threads(nth) reduction(+:sum)
    for(int i = 0; i < l; ++i) sum += px[i]; // Cannot have break statements in OpenMP for loop
  }
  pout[0] = sum;
}

// This is unsafe...
// void fsum_double_g_omp_impl(double *restrict pout, double *restrict px, int ng, int *restrict pg, int narm, int l, int nth) {
//   if(narm) {
//     for(int i = ng; i--; ) pout[i] = NA_REAL;
//     #pragma omp parallel for num_threads(nth) reduction(+:pout[:ng])
//     for(int i = 0; i < l; ++i) {
//       if(!ISNAN(px[i])) {
//         if(ISNAN(pout[pg[i]-1])) pout[pg[i]-1] = px[i];
//         else pout[pg[i]-1] += px[i];
//       }
//     }
//   } else {
//     memset(pout, 0.0, sizeof(double) * ng);
//     #pragma omp parallel for num_threads(nth) reduction(+:pout[:ng]) // shared(pout)
//     for(int i = 0; i < l; ++i) {
//       // #pragma omp atomic
//       pout[pg[i]-1] += px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
//     }
//   }
// }

void fsum_weights_impl(double *restrict pout, const double *restrict px, const double *restrict pw, const int narm, const int l) {
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
}

void fsum_weights_g_impl(double *restrict pout, const double *restrict px, const int ng, const int *restrict pg, const double *restrict pw, const int narm, const int l) {
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


void fsum_weights_omp_impl(double *restrict pout, const double *restrict px, const double *restrict pw, const int narm, const int l, const int nth) {
  double sum;
  if(narm) {
    int j = 0;
    while((ISNAN(px[j]) || ISNAN(pw[j])) && j!=l) ++j;
    if(j != l) {
      sum = px[j] * pw[j];
      #pragma omp parallel for num_threads(nth) reduction(+:sum)
      for(int i = j+1; i < l; ++i) {
        if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
        sum += px[i] * pw[i];
      }
    } else sum = NA_REAL;
  } else {
    sum = 0;
    #pragma omp parallel for num_threads(nth) reduction(+:sum)
    for(int i = 0; i < l; ++i) sum += px[i] * pw[i];
  }
  pout[0] = sum;
}

// This is unsafe...
// void fsum_weights_g_omp_impl(double *restrict pout, double *restrict px, int ng, int *restrict pg, double *restrict pw, int narm, int l, int nth) {
//   if(narm) {
//     for(int i = ng; i--; ) pout[i] = NA_REAL;
//     #pragma omp parallel for num_threads(nth) reduction(+:pout[:ng])
//     for(int i = 0; i < l; ++i) {
//       if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
//       if(ISNAN(pout[pg[i]-1])) pout[pg[i]-1] = px[i] * pw[i];
//       else pout[pg[i]-1] += px[i] * pw[i];
//     }
//   } else {
//     memset(pout, 0.0, sizeof(double) * ng);
//     #pragma omp parallel for num_threads(nth) reduction(+:pout[:ng])
//     for(int i = 0; i < l; ++i) pout[pg[i]-1] += px[i] * pw[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
//   }
// }


// using long long internally is substantially faster than using doubles !!
double fsum_int_impl(const int *restrict px, const int narm, const int l) {
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

void fsum_int_g_impl(int *restrict pout, const int *restrict px, const int ng, const int *restrict pg, const int narm, const int l) {
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

double fsum_int_omp_impl(const int *restrict px, const int narm, const int l, const int nth) {
  long long sum;
  if(narm) {
    int j = 0;
    while(px[j] == NA_INTEGER && j!=l) ++j;
    if(j == l && (l > 1 || px[j-1] == NA_INTEGER)) return NA_REAL;
    sum = (long long)px[j];
    #pragma omp parallel for num_threads(nth) reduction(+:sum)
    for(int i = j+1; i < l; ++i) if(px[i] != NA_INTEGER) sum += (long long)px[i];
  } else {
    sum = 0;
    #pragma omp parallel for num_threads(nth) reduction(+:sum)
    for(int i = 0; i < l; ++i) sum += (long long)px[i];
  }
  return (double)sum;
}

// This is unsafe...
// void fsum_int_g_omp_impl(int *restrict pout, int *restrict px, int ng, int *restrict pg, int narm, int l, int nth) {
//   long long ckof;
//   if(narm) {
//     for(int i = ng; i--; ) pout[i] = NA_INTEGER;
//     int lsi;
//     #pragma omp parallel for num_threads(nth) reduction(+:pout[:ng])
//     for(int i = 0; i < l; ++i) {
//       if(px[i] != NA_INTEGER) {
//         lsi = pout[pg[i]-1];
//         if(lsi == NA_INTEGER) pout[pg[i]-1] = px[i];
//         else {
//           ckof = (long long)lsi + px[i];
//           if(ckof > INT_MAX || ckof <= INT_MIN) error("Integer overflow in one or more groups. Integers in R are bounded between 2,147,483,647 and -2,147,483,647. The sum within each group should be in that range.");
//           pout[pg[i]-1] = (int)ckof;
//         }
//       }
//     }
//   } else {
//     memset(pout, 0, sizeof(int) * ng);
//     int lsi;
//     #pragma omp parallel for num_threads(nth) reduction(+:pout[:ng])
//     for(int i = 0; i < l; ++i) {
//       if(px[i] == NA_INTEGER) {
//         pout[pg[i]-1] = NA_INTEGER;
//         continue;
//       }
//       lsi = pout[pg[i]-1];
//       if(lsi != NA_INTEGER) { // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
//         ckof = (long long)lsi + px[i];
//         if(ckof > INT_MAX || ckof <= INT_MIN) error("Integer overflow in one or more groups. Integers in R are bounded between 2,147,483,647 and -2,147,483,647. The sum within each group should be in that range.");
//         pout[pg[i]-1] = (int)ckof;
//       }
//     }
//   }
// }


SEXP fsumC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rnth) {
  int l = length(x), tx = TYPEOF(x), ng = asInteger(Rng),
    narm = asLogical(Rnarm), nth = asInteger(Rnth), nprotect = 0, nwl = isNull(w);
  // ALTREP methods for compact sequences: not safe yet and not part of the API.
  // if(ALTREP(x) && ng == 0 && nwl) {
  // switch(tx) {
  // case INTSXP: return ALTINTEGER_SUM(x, (Rboolean)narm);
  // case LGLSXP: return ALTLOGICAL_SUM(x, (Rboolean)narm);
  // case REALSXP: return ALTREAL_SUM(x, (Rboolean)narm);
  // default: error("ALTREP object must be integer or real typed");
  // }
  // }
  if (l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(ng && l != length(g)) error("length(g) must match length(x)");
  if(l < 100000) nth = 1; // No improvements from multithreading on small data.
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out;
  if(!(ng == 0 && nwl && tx == INTSXP)) {
    out = PROTECT(allocVector(nwl ? tx : REALSXP, ng == 0 ? 1 : ng));
    ++nprotect;
  }
  if(nwl) {
    switch(tx) {
      case REALSXP:
        if(ng == 0) {
          if(nth <= 1) fsum_double_impl(REAL(out), REAL(x), narm, l);
          else fsum_double_omp_impl(REAL(out), REAL(x), narm, l, nth);
        } else fsum_double_g_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l);
        // If safe sub-column-level mutithreading can be developed...
        // if(nth <= 1) {
        //   if(ng == 0) fsum_double_impl(REAL(out), REAL(x), narm, l);
        //   else fsum_double_g_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l);
        // } else {
        //   if(ng == 0) fsum_double_omp_impl(REAL(out), REAL(x), narm, l, nth);
        //   else fsum_double_g_omp_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l, nth);
        // }
        break;
      case INTSXP: {
        if(ng > 0) {
          fsum_int_g_impl(INTEGER(out), INTEGER(x), ng, INTEGER(g), narm, l);
          // If safe sub-column-level mutithreading can be developed...
          // if(nth <= 1) fsum_int_g_impl(INTEGER(out), INTEGER(x), ng, INTEGER(g), narm, l);
          // else fsum_int_g_omp_impl(INTEGER(out), INTEGER(x), ng, INTEGER(g), narm, l, nth);
        } else {
          double sum = nth <= 1 ? fsum_int_impl(INTEGER(x), narm, l) : fsum_int_omp_impl(INTEGER(x), narm, l, nth);
          if(sum > INT_MAX || sum <= INT_MIN) return ScalarReal(sum); // INT_MIN is NA_INTEGER
          return ScalarInteger(ISNAN(sum) ? NA_INTEGER : (int)sum);
        }
        break;
      }
      default: error("Unsupported SEXP type");
    }
  } else {
    if(l != length(w)) error("length(w) must match length(x)");
    int tw = TYPEOF(w);
    SEXP xr, wr;
    double *restrict px, *restrict pw;
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
    if(ng == 0) {
      if(nth <= 1) fsum_weights_impl(REAL(out), px, pw, narm, l);
      else fsum_weights_omp_impl(REAL(out), px, pw, narm, l, nth);
    } else fsum_weights_g_impl(REAL(out), px, ng, INTEGER(g), pw, narm, l);
  }
  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")))
    copyMostAttrib(x, out); // For example "Units" objects...
  UNPROTECT(nprotect);
  return out;
}

SEXP fsummC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnth) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1], *restrict pg = INTEGER(g),
      ng = asInteger(Rng), // ng1 = ng == 0 ? 1 : ng,
      narm = asLogical(Rnarm), nprotect = 1, nwl = isNull(w),
      nth = asInteger(Rnth); // , cmth = nth > 1 && col >= nth;
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(l*col < 100000) nth = 1; // No gains from multithreading on small data
  if(ng && l != length(g)) error("length(g) must match nrow(x)");
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector((nwl && ng > 0) ? tx : REALSXP, ng == 0 ? col : col * ng));
  if(nwl) {
    switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        if(ng == 0) {
          if(nth <= 1) {
            for(int j = 0; j != col; ++j) fsum_double_impl(pout + j, px + j*l, narm, l);
          } else if(col >= nth) {
            #pragma omp parallel for num_threads(nth)
            for(int j = 0; j < col; ++j) fsum_double_impl(pout + j, px + j*l, narm, l);
          } else {
            for(int j = 0; j != col; ++j) fsum_double_omp_impl(pout + j, px + j*l, narm, l, nth);
          }
        } else {
          if(nth <= 1 || col == 1) {
            for(int j = 0; j != col; ++j) fsum_double_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
          } else {
            if(nth > col) nth = col;
            #pragma omp parallel for num_threads(nth)
            for(int j = 0; j < col; ++j) fsum_double_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
          }
        }
        break;
      }
      case INTSXP: {
        int *px = INTEGER(x);
        if(ng > 0) {
          int *pout = INTEGER(out);
          if(nth <= 1 || col == 1) {
            for(int j = 0; j != col; ++j) fsum_int_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
          } else {
            if(nth > col) nth = col;
            #pragma omp parallel for num_threads(nth)
            for(int j = 0; j < col; ++j) fsum_int_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
          }
        } else {
          double *restrict pout = REAL(out);
          int anyoutl = 0;
          if(nth <= 1) {
            for(int j = 0; j != col; ++j) {
              double sumj = fsum_int_impl(px + j*l, narm, l);
              if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
              pout[j] = sumj;
            }
          } else if(col >= nth) { // If high-dimensional: column-level parallelism
            #pragma omp parallel for num_threads(nth)
            for(int j = 0; j < col; ++j) {
              double sumj = fsum_int_impl(px + j*l, narm, l);
              if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
              pout[j] = sumj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              double sumj = fsum_int_omp_impl(px + j*l, narm, l, nth);
              if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
              pout[j] = sumj;
            }
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
    double *px, *restrict pw, *pout = REAL(out);
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
    if(ng == 0) {
      if(nth <= 1) {
        for(int j = 0; j != col; ++j) fsum_weights_impl(pout + j, px + j*l, pw, narm, l);
      } else if(col >= nth) {
        #pragma omp parallel for num_threads(nth)
        for(int j = 0; j < col; ++j) fsum_weights_impl(pout + j, px + j*l, pw, narm, l);
      } else {
        for(int j = 0; j != col; ++j) fsum_weights_omp_impl(pout + j, px + j*l, pw, narm, l, nth);
      }
    } else {
      if(nth <= 1 || col == 1) {
        for(int j = 0; j != col; ++j) fsum_weights_g_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l);
      } else {
        if(nth > col) nth = col;
        #pragma omp parallel for num_threads(nth)
        for(int j = 0; j < col; ++j) fsum_weights_g_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l);
      }
    }
  }
  matCopyAttr(out, x, Rdrop, ng);
  UNPROTECT(nprotect);
  return out;
}

SEXP fsumlC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnth) {
  int l = length(x), ng = asInteger(Rng), nth = asInteger(Rnth), nprotect = 1;
  // TODO: Disable multithreading if overall data size is small?
  if(l < 1) return x; // needed ??
  if(ng == 0 && asLogical(Rdrop)) {
    SEXP out = PROTECT(allocVector(REALSXP, l)), *restrict px = SEXPPTR(x);
    double *restrict pout = REAL(out);
    if(nth > 1 && l >= nth) { // If high-dimensional: column-level parallelism
      SEXP Rnth1 = PROTECT(ScalarInteger(1)); ++nprotect;
      #pragma omp parallel for num_threads(nth)
      for(int j = 0; j < l; ++j) pout[j] = asReal(fsumC(px[j], Rng, g, w, Rnarm, Rnth1));
    } else {
      for(int j = 0; j != l; ++j) pout[j] = asReal(fsumC(px[j], Rng, g, w, Rnarm, Rnth));
    }
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(nprotect);
    return out;
  }
  SEXP out = PROTECT(allocVector(VECSXP, l)), *restrict pout = SEXPPTR(out), *restrict px = SEXPPTR(x);
  if((ng > 0 && nth > 1 && l > 1) || (ng == 0 && nth > 1 && nth >= l)) {
    if(nth > l) nth = l;
    SEXP Rnth1 = PROTECT(ScalarInteger(1)); ++nprotect; // Needed if ng == 0, otherwise double multithreading
    #pragma omp parallel for num_threads(nth)
    for(int j = 0; j < l; ++j) pout[j] = fsumC(px[j], Rng, g, w, Rnarm, Rnth1);
  } else {
    for(int j = 0; j != l; ++j) pout[j] = fsumC(px[j], Rng, g, w, Rnarm, Rnth);
  }
  // if(ng == 0) for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]);
  DFcopyAttr(out, x, ng);
  UNPROTECT(nprotect);
  return out;
}

// If effective sub-column-level multithreading can be developed...
// SEXP fsummC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnth) {
//   SEXP dim = getAttrib(x, R_DimSymbol);
//   if(isNull(dim)) error("x is not a matrix");
//   int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1], *restrict pg = INTEGER(g),
//     ng = asInteger(Rng), // ng1 = ng == 0 ? 1 : ng,
//     narm = asLogical(Rnarm), nprotect = 1, nwl = isNull(w),
//     nth = asInteger(Rnth), cmth = nth > 1 && col >= nth;
//   if (l < 1) return x; // Prevents seqfault for numeric(0) #101
//   if(nth < 100000) nth = 1; // No gains from multithreading on small data
//   if(ng && l != length(g)) error("length(g) must match nrow(x)");
//   if(tx == LGLSXP) tx = INTSXP;
//   SEXP out = PROTECT(allocVector((nwl && ng > 0) ? tx : REALSXP, ng == 0 ? col : col * ng));
//   if(nwl) {
//     switch(tx) {
//     case REALSXP: {
//       double *px = REAL(x), *pout = REAL(out);
//       if(nth <= 1) { // No multithreading
//         if(ng == 0) for(int j = 0; j != col; ++j) fsum_double_impl(pout + j, px + j*l, narm, l);
//         else for(int j = 0; j != col; ++j) fsum_double_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
//       } else { // Multithreading
//         if(ng == 0) {
//           if(cmth) { // If high-dimensional: column-level parallelism
// #pragma omp parallel for num_threads(nth)
//             for(int j = 0; j < col; ++j) fsum_double_impl(pout + j, px + j*l, narm, l);
//           } else {
//             for(int j = 0; j != col; ++j) fsum_double_omp_impl(pout + j, px + j*l, narm, l, nth);
//           }
//         } else {
//           if(cmth) { // If high-dimensional: column-level parallelism
// #pragma omp parallel for num_threads(nth)
//             for(int j = 0; j < col; ++j) fsum_double_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
//           } else {
//             for(int j = 0; j != col; ++j) fsum_double_g_omp_impl(pout + j*ng, px + j*l, ng, pg, narm, l, nth);
//           }
//         }
//       }
//       break;
//     }
//     case INTSXP: {
//       int *px = INTEGER(x);
//       if(ng > 0) {
//         int *pout = INTEGER(out);
//         if(nth <= 1) {
//           for(int j = 0; j != col; ++j) fsum_int_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
//         } else if(cmth) { // If high-dimensional: column-level parallelism
// #pragma omp parallel for num_threads(nth)
//           for(int j = 0; j < col; ++j) fsum_int_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
//         } else {
//           for(int j = 0; j != col; ++j) fsum_int_g_omp_impl(pout + j*ng, px + j*l, ng, pg, narm, l, nth);
//         }
//       } else {
//         double *pout = REAL(out);
//         int anyoutl = 0;
//         if(nth <= 1) {
//           for(int j = 0; j != col; ++j) {
//             double sumj = fsum_int_impl(px + j*l, narm, l);
//             if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
//             pout[j] = sumj;
//           }
//         } else if(cmth) { // If high-dimensional: column-level parallelism
// #pragma omp parallel for num_threads(nth)
//           for(int j = 0; j < col; ++j) {
//             double sumj = fsum_int_impl(px + j*l, narm, l);
//             if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
//             pout[j] = sumj;
//           }
//         } else {
//           for(int j = 0; j != col; ++j) {
//             double sumj = fsum_int_omp_impl(px + j*l, narm, l, nth);
//             if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
//             pout[j] = sumj;
//           }
//         }
//         if(anyoutl == 0) {
//           SEXP iout = PROTECT(coerceVector(out, INTSXP));
//           matCopyAttr(iout, x, Rdrop, ng);
//           UNPROTECT(2);
//           return iout;
//         }
//       }
//       break;
//     }
//     default: error("Unsupported SEXP type");
//     }
//   } else {
//     if(l != length(w)) error("length(w) must match nrow(x)");
//     int tw = TYPEOF(w);
//     SEXP xr, wr;
//     double *px, *pw, *pout = REAL(out);
//     if(tw != REALSXP) {
//       if(tw != INTSXP && tw != LGLSXP) error("weigths must be double or integer");
//       wr = PROTECT(coerceVector(w, REALSXP));
//       pw = REAL(wr);
//       ++nprotect;
//     } else pw = REAL(w);
//     if(tx != REALSXP) {
//       if(tx != INTSXP) error("x must be double or integer");
//       xr = PROTECT(coerceVector(x, REALSXP));
//       px = REAL(xr);
//       ++nprotect;
//     } else px = REAL(x);
//     if(nth <= 1) {
//       for(int j = 0; j != col; ++j) fsum_weights_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l);
//     } else if(cmth) {
// #pragma omp parallel for num_threads(nth)
//       for(int j = 0; j < col; ++j) fsum_weights_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l);
//     } else {
//       for(int j = 0; j != col; ++j) fsum_weights_omp_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l, nth);
//     }
//   }
//   matCopyAttr(out, x, Rdrop, ng);
//   UNPROTECT(nprotect);
//   return out;
// }

// If effective sub-column-level multithreading can be developed...
// SEXP fsumlC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnth) {
//   int l = length(x), ng = asInteger(Rng), nth = asInteger(Rnth),
//     nprotect = 1, cmth = nth > 1 && l >= nth;
//   // TODO: Disable multithreading if overall data size is small?
//   if(l < 1) return x; // needed ??
//   SEXP Rnth1;
//   if(cmth) {
//     Rnth1 = PROTECT(ScalarInteger(1));
//     ++nprotect;
//   }
//   if(ng == 0 && asLogical(Rdrop)) {
//     SEXP out = PROTECT(allocVector(REALSXP, l)), *px = SEXPPTR(x);
//     double *pout = REAL(out);
//     if(cmth) { // If high-dimensional: column-level parallelism
//       #pragma omp parallel for num_threads(nth)
//       for(int j = 0; j < l; ++j) pout[j] = asReal(fsumC(px[j], Rng, g, w, Rnarm, Rnth1));
//     } else {
//       for(int j = 0; j != l; ++j) pout[j] = asReal(fsumC(px[j], Rng, g, w, Rnarm, Rnth));
//     }
//     setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
//     UNPROTECT(nprotect);
//     return out;
//   }
//   SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out), *px = SEXPPTR(x);
//   if(cmth) {
//     #pragma omp parallel for num_threads(nth)
//     for(int j = 0; j < l; ++j) pout[j] = fsumC(px[j], Rng, g, w, Rnarm, Rnth1);
//   } else {
//     for(int j = 0; j != l; ++j) pout[j] = fsumC(px[j], Rng, g, w, Rnarm, Rnth);
//   }
//   // if(ng == 0) for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]);
//   DFcopyAttr(out, x, ng);
//   UNPROTECT(nprotect);
//   return out;
// }
