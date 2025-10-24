#include "collapse_c.h"
// #include <R_ext/Altrep.h>

double fsum_double_impl(const double *restrict px, const int narm, const int l) {
  double sum;
  if(narm == 1) {
    int j = 1;
    sum = px[0];
    while(ISNAN(sum) && j!=l) sum = px[j++];
    if(j != l) {
      #pragma omp simd reduction(+:sum)
      for(int i = j; i < l; ++i) sum += NISNAN(px[i]) ? px[i] : 0.0;
    }
  } else {
    sum = 0;
    if(narm) {
      #pragma omp simd reduction(+:sum)
      for(int i = 0; i < l; ++i) sum += NISNAN(px[i]) ? px[i] : 0.0;
    } else {
     // Should just be fast, don't stop for NA's
      #pragma omp simd reduction(+:sum)
      for(int i = 0; i < l; ++i) sum += px[i];
    }
  }
  return sum;
}

void fsum_double_g_impl(double *restrict pout, const double *restrict px, const int ng, const int *restrict pg, const int narm, const int l) {
  if(narm == 1) {
    for(int i = ng; i--; ) pout[i] = NA_REAL; // Other way ?
    --pout;
    for(int i = 0; i != l; ++i) {
      if(ISNAN(px[i])) continue; // faster way to code this ? -> Not Bad at all
      if(ISNAN(pout[pg[i]])) pout[pg[i]] = px[i];
      else pout[pg[i]] += px[i];
    }
  } else {
    memset(pout, 0, sizeof(double) * ng);
    --pout;
    if(narm == 2) {
      for(int i = 0; i != l; ++i) if(NISNAN(px[i])) pout[pg[i]] += px[i];
    } else {
      for(int i = 0; i != l; ++i) pout[pg[i]] += px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    }
  }
}

double fsum_double_omp_impl(const double *restrict px, const int narm, const int l, const int nthreads) {
  double sum;
  if(narm) {
    int j = 1;
    sum = px[0];
    while(ISNAN(sum) && j != l) sum = px[j++];
    if(j != l) {
      #pragma omp parallel for simd num_threads(nthreads) reduction(+:sum)
      for(int i = j; i < l; ++i) sum += NISNAN(px[i]) ? px[i] : 0.0;
    } else if(narm == 2) sum = 0.0;
  } else {
    sum = 0;
    #pragma omp parallel for simd num_threads(nthreads) reduction(+:sum)
    for(int i = 0; i < l; ++i) sum += px[i]; // Cannot have break statements in OpenMP for loop
  }
  return sum;
}

// This is unsafe...
// void fsum_double_g_omp_impl(double *restrict pout, double *restrict px, int ng, int *restrict pg, int narm, int l, int nthreads) {
//   if(narm) {
//     for(int i = ng; i--; ) pout[i] = NA_REAL;
//     #pragma omp parallel for num_threads(nthreads) reduction(+:pout[:ng])
//     for(int i = 0; i < l; ++i) {
//       if(!ISNAN(px[i])) {
//         if(ISNAN(pout[pg[i]-1])) pout[pg[i]-1] = px[i];
//         else pout[pg[i]-1] += px[i];
//       }
//     }
//   } else {
//     memset(pout, 0, sizeof(double) * ng);
//     #pragma omp parallel for num_threads(nthreads) reduction(+:pout[:ng]) // shared(pout)
//     for(int i = 0; i < l; ++i) {
//       // #pragma omp atomic
//       pout[pg[i]-1] += px[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
//     }
//   }
// }

double fsum_weights_impl(const double *restrict px, const double *restrict pw, const int narm, const int l) {
  double sum;
  if(narm == 1) {
    int j = 0, end = l-1;
    while((ISNAN(px[j]) || ISNAN(pw[j])) && j!=end) ++j;
    sum = px[j] * pw[j];
    if(j != end) {
      #pragma omp simd reduction(+:sum)
      for(int i = j+1; i < l; ++i) sum += (NISNAN(px[i]) && NISNAN(pw[i])) ? px[i] * pw[i] : 0.0;
    }
  } else {
    sum = 0;
    if(narm) {
      #pragma omp simd reduction(+:sum)
      for(int i = 0; i < l; ++i) sum += (NISNAN(px[i]) && NISNAN(pw[i])) ? px[i] * pw[i] : 0.0;
    } else {
      // Also here speed is key...
      #pragma omp simd reduction(+:sum)
      for(int i = 0; i < l; ++i) sum += px[i] * pw[i];
    }
  }
  return sum;
}

void fsum_weights_g_impl(double *restrict pout, const double *restrict px, const int ng, const int *restrict pg, const double *restrict pw, const int narm, const int l) {
  if(narm == 1) {
    for(int i = ng; i--; ) pout[i] = NA_REAL; // Other way ?
    --pout;
    for(int i = l; i--; ) {
      if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
      if(ISNAN(pout[pg[i]])) pout[pg[i]] = px[i] * pw[i];
      else pout[pg[i]] += px[i] * pw[i];
    }
  } else {
    memset(pout, 0, sizeof(double) * ng);
    --pout;
    if(narm == 2) {
      for(int i = l; i--; ) if(NISNAN(px[i]) && NISNAN(pw[i])) pout[pg[i]] += px[i] * pw[i];
    } else {
      for(int i = l; i--; ) pout[pg[i]] += px[i] * pw[i]; // Used to stop loop when all groups passed with NA, but probably no speed gain since groups are mostly ordered.
    }
  }
}

double fsum_weights_omp_impl(const double *restrict px, const double *restrict pw, const int narm, const int l, const int nthreads) {
  double sum;
  if(narm) {
    int j = 0;
    while(j!=l && (ISNAN(px[j]) || ISNAN(pw[j]))) ++j;
    if(j != l) {
      sum = px[j] * pw[j];
      #pragma omp parallel for simd num_threads(nthreads) reduction(+:sum)
      for(int i = j+1; i < l; ++i) sum += (NISNAN(px[i]) && NISNAN(pw[i])) ? px[i] * pw[i] : 0.0;
    } else sum = narm == 1 ? NA_REAL : 0.0;
  } else {
    sum = 0;
    #pragma omp parallel for simd num_threads(nthreads) reduction(+:sum)
    for(int i = 0; i < l; ++i) sum += px[i] * pw[i];
  }
  return sum;
}

// This is unsafe...
// void fsum_weights_g_omp_impl(double *restrict pout, double *restrict px, int ng, int *restrict pg, double *restrict pw, int narm, int l, int nthreads) {
//   if(narm) {
//     for(int i = ng; i--; ) pout[i] = NA_REAL;
//     #pragma omp parallel for num_threads(nthreads) reduction(+:pout[:ng])
//     for(int i = 0; i < l; ++i) {
//       if(ISNAN(px[i]) || ISNAN(pw[i])) continue;
//       if(ISNAN(pout[pg[i]-1])) pout[pg[i]-1] = px[i] * pw[i];
//       else pout[pg[i]-1] += px[i] * pw[i];
//     }
//   } else {
//     memset(pout, 0, sizeof(double) * ng);
//     #pragma omp parallel for num_threads(nthreads) reduction(+:pout[:ng])
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
    if(j == 0 && px[j] == NA_INTEGER) return narm == 1 ? NA_REAL : 0;
    for(int i = j; i--; ) if(px[i] != NA_INTEGER) sum += (long long)px[i];
  } else {
    sum = 0;
    for(int i = 0; i != l; ++i) {
      if(px[i] == NA_INTEGER) return NA_REAL; // Need this, otherwise result is incorrect !!
      sum += (long long)px[i];
    }
  }
  return (double)sum;
}

void fsum_int_g_impl(int *restrict pout, const int *restrict px, const int ng, const int *restrict pg, const int narm, const int l) {
  long long ckof;
  if(narm == 1) {
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
    if(narm == 2) {
      for(int i = l; i--; ) {
        if(px[i] != NA_INTEGER) {
          ckof = (long long)pout[pg[i]] + px[i];
          if(ckof > INT_MAX || ckof <= INT_MIN) error("Integer overflow in one or more groups. Integers in R are bounded between 2,147,483,647 and -2,147,483,647. The sum within each group should be in that range.");
          pout[pg[i]] = (int)ckof;
        }
      }
    } else {
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
}

double fsum_int_omp_impl(const int *restrict px, const int narm, const int l, const int nthreads) {
  long long sum;
  if(narm) {
    int j = 0;
    while(px[j] == NA_INTEGER && j!=l) ++j;
    if(j == l && px[j-1] == NA_INTEGER) return narm == 1 ? NA_REAL : 0;
    sum = (long long)px[j];
    #pragma omp parallel for simd num_threads(nthreads) reduction(+:sum)
    for(int i = j+1; i < l; ++i) sum += px[i] != NA_INTEGER ? (long long)px[i] : 0;
  } else {
    if(px[0] == NA_INTEGER || px[l-1] == NA_INTEGER) return NA_REAL;
    sum = 0;
    #pragma omp parallel for simd num_threads(nthreads) reduction(+:sum)
    for(int i = 0; i < l; ++i) sum += (long long)px[i]; // Need this, else wrong result
  }
  return (double)sum;
}

// This is unsafe...
// void fsum_int_g_omp_impl(int *restrict pout, int *restrict px, int ng, int *restrict pg, int narm, int l, int nthreads) {
//   long long ckof;
//   if(narm) {
//     for(int i = ng; i--; ) pout[i] = NA_INTEGER;
//     int lsi;
//     #pragma omp parallel for num_threads(nthreads) reduction(+:pout[:ng])
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
//     #pragma omp parallel for num_threads(nthreads) reduction(+:pout[:ng])
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


SEXP fsumC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP fill, SEXP Rnthreads) {
  int l = length(x), tx = TYPEOF(x), ng = asInteger(Rng),
    narm = asLogical(Rnarm), nthreads = asInteger(Rnthreads), nprotect = 0, nwl = isNull(w);
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
  if(l < 100000) nthreads = 1; // No improvements from multithreading on small data.
  if(narm) narm += asLogical(fill);
  if(nthreads > max_threads) nthreads = max_threads;
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
          REAL(out)[0] = (nthreads <= 1) ? fsum_double_impl(REAL(x), narm, l) :
                        fsum_double_omp_impl(REAL(x), narm, l, nthreads);
        } else fsum_double_g_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l);
        // If safe sub-column-level mutithreading can be developed...
        // if(nthreads <= 1) {
        //   if(ng == 0) fsum_double_impl(REAL(out), REAL(x), narm, l);
        //   else fsum_double_g_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l);
        // } else {
        //   if(ng == 0) fsum_double_omp_impl(REAL(out), REAL(x), narm, l, nthreads);
        //   else fsum_double_g_omp_impl(REAL(out), REAL(x), ng, INTEGER(g), narm, l, nthreads);
        // }
        break;
      case INTSXP: {
        if(ng > 0) {
          fsum_int_g_impl(INTEGER(out), INTEGER(x), ng, INTEGER(g), narm, l);
          // If safe sub-column-level mutithreading can be developed...
          // if(nthreads <= 1) fsum_int_g_impl(INTEGER(out), INTEGER(x), ng, INTEGER(g), narm, l);
          // else fsum_int_g_omp_impl(INTEGER(out), INTEGER(x), ng, INTEGER(g), narm, l, nthreads);
        } else {
          double sum = nthreads <= 1 ? fsum_int_impl(INTEGER(x), narm, l) : fsum_int_omp_impl(INTEGER(x), narm, l, nthreads);
          UNPROTECT(nprotect); // Thomas Kalibera Patch: to appease rchk.
          if(sum > INT_MAX || sum <= INT_MIN) return ScalarReal(sum); // INT_MIN is NA_INTEGER
          return ScalarInteger(ISNAN(sum) ? NA_INTEGER : (int)sum);
        }
        break;
      }
    default: error("Unsupported SEXP type: '%s'", type2char(tx));
    }
  } else {
    if(l != length(w)) error("length(w) must match length(x)");
    if(TYPEOF(w) != REALSXP) {
      if(TYPEOF(w) != INTSXP && TYPEOF(w) != LGLSXP) error("weights must be double or integer");
      w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
    }
    if(tx != REALSXP) {
      if(tx != INTSXP) error("Unsupported SEXP type: '%s'", type2char(tx));
      x = PROTECT(coerceVector(x, REALSXP)); ++nprotect;
    }
    double *restrict px = REAL(x), *restrict pw = REAL(w);
    if(ng == 0) {
      REAL(out)[0] = (nthreads <= 1) ? fsum_weights_impl(px, pw, narm, l) :
               fsum_weights_omp_impl(px, pw, narm, l, nthreads);
    } else fsum_weights_g_impl(REAL(out), px, ng, INTEGER(g), pw, narm, l);
  }
  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts")))
    copyMostAttrib(x, out); // For example "Units" objects...
  UNPROTECT(nprotect);
  return out;
}

SEXP fsummC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP fill, SEXP Rdrop, SEXP Rnthreads) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1], *restrict pg = INTEGER(g),
      ng = asInteger(Rng), // ng1 = ng == 0 ? 1 : ng,
      narm = asLogical(Rnarm), nprotect = 1, nwl = isNull(w),
      nthreads = asInteger(Rnthreads); // , cmth = nthreads > 1 && col >= nthreads;
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101
  if(l*col < 100000) nthreads = 1; // No gains from multithreading on small data
  if(ng && l != length(g)) error("length(g) must match nrow(x)");
  if(narm) narm += asLogical(fill);
  if(nthreads > max_threads) nthreads = max_threads;
  if(tx == LGLSXP) tx = INTSXP;
  SEXP out = PROTECT(allocVector((nwl && ng > 0) ? tx : REALSXP, ng == 0 ? col : col * ng));
  if(nwl) {
    switch(tx) {
      case REALSXP: {
        double *px = REAL(x), *pout = REAL(out);
        if(ng == 0) {
          if(nthreads <= 1) {
            for(int j = 0; j != col; ++j) pout[j] = fsum_double_impl(px + j*l, narm, l);
          } else if(col >= nthreads) {
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) pout[j] = fsum_double_impl(px + j*l, narm, l);
          } else {
            for(int j = 0; j != col; ++j) pout[j] = fsum_double_omp_impl(px + j*l, narm, l, nthreads);
          }
        } else {
          if(nthreads <= 1 || col == 1) {
            for(int j = 0; j != col; ++j) fsum_double_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
          } else {
            if(nthreads > col) nthreads = col;
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) fsum_double_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
          }
        }
        break;
      }
      case INTSXP: {
        int *px = INTEGER(x);
        if(ng > 0) {
          int *pout = INTEGER(out);
          if(nthreads <= 1 || col == 1) {
            for(int j = 0; j != col; ++j) fsum_int_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
          } else {
            if(nthreads > col) nthreads = col;
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) fsum_int_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
          }
        } else {
          double *restrict pout = REAL(out);
          int anyoutl = 0;
          if(nthreads <= 1) {
            for(int j = 0; j != col; ++j) {
              double sumj = fsum_int_impl(px + j*l, narm, l);
              if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
              pout[j] = sumj;
            }
          } else if(col >= nthreads) { // If high-dimensional: column-level parallelism
            #pragma omp parallel for num_threads(nthreads)
            for(int j = 0; j < col; ++j) {
              double sumj = fsum_int_impl(px + j*l, narm, l);
              if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
              pout[j] = sumj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              double sumj = fsum_int_omp_impl(px + j*l, narm, l, nthreads);
              if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
              pout[j] = sumj;
            }
          }
          if(anyoutl == 0) {
            out = PROTECT(coerceVector(out, INTSXP));
            matCopyAttr(out, x, Rdrop, ng);
            UNPROTECT(nprotect + 1);
            return out;
          }
        }
        break;
      }
      default: error("Unsupported SEXP type: '%s'", type2char(tx));
    }
  } else {
    if(l != length(w)) error("length(w) must match nrow(x)");
    if(TYPEOF(w) != REALSXP) {
      if(TYPEOF(w) != INTSXP && TYPEOF(w) != LGLSXP) error("weights must be double or integer");
      w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
    }
    if(tx != REALSXP) {
      if(tx != INTSXP) error("Unsupported SEXP type: '%s'", type2char(tx));
      x = PROTECT(coerceVector(x, REALSXP)); ++nprotect;
    }
    double *px = REAL(x), *restrict pw = REAL(w), *pout = REAL(out);

    if(ng == 0) {
      if(nthreads <= 1) {
        for(int j = 0; j != col; ++j) pout[j] = fsum_weights_impl(px + j*l, pw, narm, l);
      } else if(col >= nthreads) {
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) pout[j] = fsum_weights_impl(px + j*l, pw, narm, l);
      } else {
        for(int j = 0; j != col; ++j) pout[j] = fsum_weights_omp_impl(px + j*l, pw, narm, l, nthreads);
      }
    } else {
      if(nthreads <= 1 || col == 1) {
        for(int j = 0; j != col; ++j) fsum_weights_g_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l);
      } else {
        if(nthreads > col) nthreads = col;
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < col; ++j) fsum_weights_g_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l);
      }
    }
  }
  matCopyAttr(out, x, Rdrop, ng);
  UNPROTECT(nprotect);
  return out;
}

// For safe multithreading across data frame columns

double fsum_impl_dbl(SEXP x, int narm, int nthreads) {
  int l = length(x);
  if(l < 1) return NA_REAL;
  if(nthreads <= 1) switch(TYPEOF(x)) {
    case REALSXP: return fsum_double_impl(REAL(x), narm, l);
    case LGLSXP:
    case INTSXP: return fsum_int_impl(INTEGER(x), narm, l);
    default: error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
  }
  switch(TYPEOF(x)) {
    case REALSXP: return fsum_double_omp_impl(REAL(x), narm, l, nthreads);
    case LGLSXP:
    case INTSXP: return fsum_int_omp_impl(INTEGER(x), narm, l, nthreads);
    default: error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
  }
}

SEXP fsum_impl_SEXP(SEXP x, int narm, int nthreads) {
  return ScalarReal(fsum_impl_dbl(x, narm, nthreads));
  // This is not thread safe... need to do separate serial loop
  // SEXP res = ScalarReal(fsum_impl_dbl(x, narm, nthreads));
  // if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) {
  //   PROTECT(res);
  //   copyMostAttrib(x, res);
  //   UNPROTECT(1);
  // }
  // return res;
}

double fsum_w_impl_dbl(SEXP x, double *pw, int narm, int nthreads) {
  int l = length(x);
  if(l < 1) return NA_REAL;
  if(TYPEOF(x) != REALSXP) {
    if(TYPEOF(x) != INTSXP && TYPEOF(x) != LGLSXP) error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
    x = PROTECT(coerceVector(x, REALSXP));
    double res = (nthreads <= 1) ? fsum_weights_impl(REAL(x), pw, narm, l) :
      fsum_weights_omp_impl(REAL(x), pw, narm, l, nthreads);
    UNPROTECT(1);
    return res;
  }
  return (nthreads <= 1) ? fsum_weights_impl(REAL(x), pw, narm, l) :
    fsum_weights_omp_impl(REAL(x), pw, narm, l, nthreads);
}

SEXP fsum_w_impl_SEXP(SEXP x, double *pw, int narm, int nthreads) {
  return ScalarReal(fsum_w_impl_dbl(x, pw, narm, nthreads));
  // This is not thread safe... need to do separate serial loop
  // SEXP res = ScalarReal(fsum_w_impl_dbl(x, pw, narm, nthreads));
  // if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) {
  //   PROTECT(res);
  //   copyMostAttrib(x, res);
  //   UNPROTECT(1);
  // }
  // return res;
}

SEXP fsum_g_impl(SEXP x, const int ng, const int *pg, int narm) {
  int l = length(x);
  if(l < 1) return ScalarReal(NA_REAL);

  SEXP res;
  switch(TYPEOF(x)) {
    case REALSXP: {
      res = PROTECT(allocVector(REALSXP, ng));
      fsum_double_g_impl(REAL(res), REAL(x), ng, pg, narm, l);
      break;
    }
    case LGLSXP:
    case INTSXP:  {
      res = PROTECT(allocVector(INTSXP, ng));
      fsum_int_g_impl(INTEGER(res), INTEGER(x), ng, pg, narm, l);
      break;
    }
    default: error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
  }

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) copyMostAttrib(x, res);
  UNPROTECT(1);
  return res;
}

void fsum_g_omp_impl(SEXP x, void *pres, const int ng, const int *pg, int narm) {
  switch(TYPEOF(x)) {
    case REALSXP:
      fsum_double_g_impl(pres, REAL(x), ng, pg, narm, length(x));
      break;
    case LGLSXP:
    case INTSXP:
      fsum_int_g_impl(pres, INTEGER(x), ng, pg, narm, length(x));
      break;
    default: error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
  }
}

SEXP fsum_wg_impl(SEXP x, const int ng, const int *pg, double *pw, int narm) {
  int l = length(x), nprotect = 1;
  if(l < 1) return ScalarReal(NA_REAL);

  if(TYPEOF(x) != REALSXP) {
    if(TYPEOF(x) != INTSXP && TYPEOF(x) != LGLSXP) error("Unsupported SEXP type: '%s'", type2char(TYPEOF(x)));
    x = PROTECT(coerceVector(x, REALSXP)); ++nprotect;
  }

  SEXP res = PROTECT(allocVector(REALSXP, ng));
  fsum_weights_g_impl(REAL(res), REAL(x), ng, pg, pw, narm, l);

  if(ATTRIB(x) != R_NilValue && !(isObject(x) && inherits(x, "ts"))) copyMostAttrib(x, res);
  UNPROTECT(nprotect);
  return res;
}




#undef COLWISE_FSUM_LIST
#define COLWISE_FSUM_LIST(FUN, WFUN)                           \
if(nwl) {                                                      \
  if(nthreads > 1 && l >= nthreads) {                          \
    _Pragma("omp parallel for num_threads(nthreads)")          \
    for(int j = 0; j < l; ++j) pout[j] = FUN(px[j], narm, 1);  \
  } else {                                                     \
    for(int j = 0; j != l; ++j) pout[j] = FUN(px[j], narm, nthreads); \
  }                                                            \
} else {                                                       \
  double *restrict pw = REAL(w);                               \
  if(nthreads > 1 && l >= nthreads) {                          \
    _Pragma("omp parallel for num_threads(nthreads)")          \
    for(int j = 0; j < l; ++j) pout[j] = WFUN(px[j], pw, narm, 1); \
  } else {                                                     \
    for(int j = 0; j != l; ++j) pout[j] = WFUN(px[j], pw, narm, nthreads); \
  }                                                            \
}


SEXP fsumlC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP fill, SEXP Rdrop, SEXP Rnthreads) {
  int l = length(x), ng = asInteger(Rng), nthreads = asInteger(Rnthreads), nwl = isNull(w),
    narm = asLogical(Rnarm), nprotect = 1;

  // TODO: Disable multithreading if overall data size is small?
  if(l < 1) return x; // needed ??
  if(narm) narm += asLogical(fill);
  if(nthreads > max_threads) nthreads = max_threads;

  if(!nwl) {
    if(length(VECTOR_ELT(x, 0)) != length(w)) error("length(w) must match nrow(x)");
    if(TYPEOF(w) != REALSXP) {
      if(TYPEOF(w) != INTSXP && TYPEOF(w) != LGLSXP) error("weights must be double or integer");
      w = PROTECT(coerceVector(w, REALSXP)); ++nprotect;
    }
  }

  if(ng == 0 && asLogical(Rdrop)) {
    SEXP out = PROTECT(allocVector(REALSXP, l));
    const SEXP *restrict px = SEXPPTR_RO(x);
    double *restrict pout = REAL(out);
    COLWISE_FSUM_LIST(fsum_impl_dbl, fsum_w_impl_dbl);
    setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
    UNPROTECT(nprotect);
    return out;
  }

  SEXP out = PROTECT(allocVector(VECSXP, l)), *restrict pout = SEXPPTR(out);
  const SEXP *restrict px = SEXPPTR_RO(x);

  if(ng == 0) {
    COLWISE_FSUM_LIST(fsum_impl_SEXP, fsum_w_impl_SEXP);
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
      if(nthreads > 1 && l > 1) {
        for(int j = 0; j != l; ++j) {
          SEXP xj = px[j], outj;
          SET_VECTOR_ELT(out, j, outj = allocVector(TYPEOF(px[j]) == REALSXP ? REALSXP : INTSXP, ng));
          if(ATTRIB(xj) != R_NilValue && !(isObject(xj) && inherits(xj, "ts"))) copyMostAttrib(xj, outj);
        }
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < l; ++j) fsum_g_omp_impl(px[j], DPTR(pout[j]), ng, pg, narm);
      } else {
        for(int j = 0; j != l; ++j) SET_VECTOR_ELT(out, j, fsum_g_impl(px[j], ng, pg, narm));
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
            if(dup == 0) {x = PROTECT(shallow_duplicate(x)); ++nprotect; dup = 1;}
            SET_VECTOR_ELT(x, j, coerceVector(xj, REALSXP));
            px = SEXPPTR_RO(x); // Fix suggested by ChatGPT
          }
        }
        #pragma omp parallel for num_threads(nthreads)
        for(int j = 0; j < l; ++j) fsum_weights_g_impl(REAL(pout[j]), REAL(px[j]), ng, pg, pw, narm, nrx);
      } else {
        for(int j = 0; j != l; ++j) SET_VECTOR_ELT(out, j, fsum_wg_impl(px[j], ng, pg, pw, narm));
      }
    }
  }

  DFcopyAttr(out, x, ng);
  UNPROTECT(nprotect);
  return out;
}

// If effective sub-column-level multithreading can be developed...
// SEXP fsummC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads) {
//   SEXP dim = getAttrib(x, R_DimSymbol);
//   if(isNull(dim)) error("x is not a matrix");
//   int tx = TYPEOF(x), l = INTEGER(dim)[0], col = INTEGER(dim)[1], *restrict pg = INTEGER(g),
//     ng = asInteger(Rng), // ng1 = ng == 0 ? 1 : ng,
//     narm = asLogical(Rnarm), nprotect = 1, nwl = isNull(w),
//     nthreads = asInteger(Rnthreads), cmth = nthreads > 1 && col >= nthreads;
//   if (l < 1) return x; // Prevents seqfault for numeric(0) #101
//   if(nthreads < 100000) nthreads = 1; // No gains from multithreading on small data
//   if(ng && l != length(g)) error("length(g) must match nrow(x)");
//   if(tx == LGLSXP) tx = INTSXP;
//   SEXP out = PROTECT(allocVector((nwl && ng > 0) ? tx : REALSXP, ng == 0 ? col : col * ng));
//   if(nwl) {
//     switch(tx) {
//     case REALSXP: {
//       double *px = REAL(x), *pout = REAL(out);
//       if(nthreads <= 1) { // No multithreading
//         if(ng == 0) for(int j = 0; j != col; ++j) fsum_double_impl(pout + j, px + j*l, narm, l);
//         else for(int j = 0; j != col; ++j) fsum_double_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
//       } else { // Multithreading
//         if(ng == 0) {
//           if(cmth) { // If high-dimensional: column-level parallelism
// #pragma omp parallel for num_threads(nthreads)
//             for(int j = 0; j < col; ++j) fsum_double_impl(pout + j, px + j*l, narm, l);
//           } else {
//             for(int j = 0; j != col; ++j) fsum_double_omp_impl(pout + j, px + j*l, narm, l, nthreads);
//           }
//         } else {
//           if(cmth) { // If high-dimensional: column-level parallelism
// #pragma omp parallel for num_threads(nthreads)
//             for(int j = 0; j < col; ++j) fsum_double_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
//           } else {
//             for(int j = 0; j != col; ++j) fsum_double_g_omp_impl(pout + j*ng, px + j*l, ng, pg, narm, l, nthreads);
//           }
//         }
//       }
//       break;
//     }
//     case INTSXP: {
//       int *px = INTEGER(x);
//       if(ng > 0) {
//         int *pout = INTEGER(out);
//         if(nthreads <= 1) {
//           for(int j = 0; j != col; ++j) fsum_int_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
//         } else if(cmth) { // If high-dimensional: column-level parallelism
// #pragma omp parallel for num_threads(nthreads)
//           for(int j = 0; j < col; ++j) fsum_int_g_impl(pout + j*ng, px + j*l, ng, pg, narm, l);
//         } else {
//           for(int j = 0; j != col; ++j) fsum_int_g_omp_impl(pout + j*ng, px + j*l, ng, pg, narm, l, nthreads);
//         }
//       } else {
//         double *pout = REAL(out);
//         int anyoutl = 0;
//         if(nthreads <= 1) {
//           for(int j = 0; j != col; ++j) {
//             double sumj = fsum_int_impl(px + j*l, narm, l);
//             if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
//             pout[j] = sumj;
//           }
//         } else if(cmth) { // If high-dimensional: column-level parallelism
// #pragma omp parallel for num_threads(nthreads)
//           for(int j = 0; j < col; ++j) {
//             double sumj = fsum_int_impl(px + j*l, narm, l);
//             if(sumj > INT_MAX || sumj <= INT_MIN) anyoutl = 1;
//             pout[j] = sumj;
//           }
//         } else {
//           for(int j = 0; j != col; ++j) {
//             double sumj = fsum_int_omp_impl(px + j*l, narm, l, nthreads);
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
//       if(tw != INTSXP && tw != LGLSXP) error("weights must be double or integer");
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
//     if(nthreads <= 1) {
//       for(int j = 0; j != col; ++j) fsum_weights_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l);
//     } else if(cmth) {
// #pragma omp parallel for num_threads(nthreads)
//       for(int j = 0; j < col; ++j) fsum_weights_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l);
//     } else {
//       for(int j = 0; j != col; ++j) fsum_weights_omp_impl(pout + j*ng, px + j*l, ng, pg, pw, narm, l, nthreads);
//     }
//   }
//   matCopyAttr(out, x, Rdrop, ng);
//   UNPROTECT(nprotect);
//   return out;
// }

// If effective sub-column-level multithreading can be developed...
// SEXP fsumlC(SEXP x, SEXP Rng, SEXP g, SEXP w, SEXP Rnarm, SEXP Rdrop, SEXP Rnthreads) {
//   int l = length(x), ng = asInteger(Rng), nthreads = asInteger(Rnthreads),
//     nprotect = 1, cmth = nthreads > 1 && l >= nthreads;
//   // TODO: Disable multithreading if overall data size is small?
//   if(l < 1) return x; // needed ??
//   SEXP Rnthreads1;
//   if(cmth) {
//     Rnthreads1 = PROTECT(ScalarInteger(1));
//     ++nprotect;
//   }
//   if(ng == 0 && asLogical(Rdrop)) {
//     SEXP out = PROTECT(allocVector(REALSXP, l)), *px = SEXPPTR(x);
//     double *pout = REAL(out);
//     if(cmth) { // If high-dimensional: column-level parallelism
//       #pragma omp parallel for num_threads(nthreads)
//       for(int j = 0; j < l; ++j) pout[j] = asReal(fsumC(px[j], Rng, g, w, Rnarm, Rnthreads1));
//     } else {
//       for(int j = 0; j != l; ++j) pout[j] = asReal(fsumC(px[j], Rng, g, w, Rnarm, Rnthreads));
//     }
//     setAttrib(out, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
//     UNPROTECT(nprotect);
//     return out;
//   }
//   SEXP out = PROTECT(allocVector(VECSXP, l)), *pout = SEXPPTR(out), *px = SEXPPTR(x);
//   if(cmth) {
//     #pragma omp parallel for num_threads(nthreads)
//     for(int j = 0; j < l; ++j) pout[j] = fsumC(px[j], Rng, g, w, Rnarm, Rnthreads1);
//   } else {
//     for(int j = 0; j != l; ++j) pout[j] = fsumC(px[j], Rng, g, w, Rnarm, Rnthreads);
//   }
//   // if(ng == 0) for(int j = 0; j != l; ++j) copyMostAttrib(px[j], pout[j]);
//   DFcopyAttr(out, x, ng);
//   UNPROTECT(nprotect);
//   return out;
// }
