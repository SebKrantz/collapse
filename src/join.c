#include "collapse_c.h" // Needs to be first because includes OpenMP, to avoid namespace conflicts.
#include "kit.h"
#include "base_radixsort.h"

/* A Sort-Merge Join
 See: https://www.dcs.ed.ac.uk/home/tz/phd/thesis/node20.htm
 And: https://en.wikipedia.org/wiki/Sort-merge_join
*/

void sort_merge_join_int(const int *restrict px, const int *restrict pt, // Data pointers, decremented by 1
                         const int *restrict pox, const int *restrict pot, // Ordering vectors
                         const int nx, const int nt, int *restrict pres, // Sizes and result vector, pres should also be decremented by 1
                         const int second_pass) { // Is it the first or second pass?
  int i = 0, j = 0, tmp;
  while (i != nx && j != nt) {
    if(second_pass && pres[pox[i]] == NA_INTEGER) {
      ++i; continue; // For second passes (Need to initialize pres or add argument to function)
    }
    tmp = pt[pot[j]];
    if (px[pox[i]] == tmp) {
      pres[pox[i]] = pot[j];
      // This takes care of duplicates in x and table
      while (++i != nx && px[pox[i]] == tmp) pres[pox[i]] = pot[j];
      while (++j != nt && pt[pot[j]] == tmp);
    } else if (px[pox[i]] < tmp) { // No extra condition needed here because NA_INTEGER is the smallest integer
      pres[pox[i++]] = NA_INTEGER;
    } else ++j;
  }
  if(!second_pass) while (i < nx) pres[pox[i++]] = NA_INTEGER;
}

void sort_merge_join_double(const double *restrict px, const double *restrict pt, // Data pointers, decremented by 1
                            const int *restrict pox, const int *restrict pot, // Ordering vectors
                            const int nx, const int nt, int *restrict pres, // Sizes and result vector, pres should also be decremented by 1
                            const int second_pass) { // Is it the first or second pass?
  int i = 0, j = 0;
  double tmp;
  while (i != nx && j != nt) {
    if(second_pass && pres[pox[i]] == NA_INTEGER) {
      ++i; continue; // For second passes (Need to initialize pres or add argument to function)
    }
    tmp = pt[pot[j]];
    if (REQUAL(px[pox[i]], tmp)) {
      pres[pox[i]] = pot[j];
      // This takes care of duplicates in x and table
      while (++i != nx && REQUAL(px[pox[i]], tmp)) pres[pox[i]] = pot[j];
      while (++j != nt && REQUAL(pt[pot[j]], tmp));
    } else if (ISNAN(px[pox[i]]) || px[pox[i]] < tmp) {
      pres[pox[i++]] = NA_INTEGER;
    } else ++j;
  }
  if(!second_pass) while (i < nx) pres[pox[i++]] = NA_INTEGER;
}


void sort_merge_join_string(const SEXP *restrict px, const SEXP *restrict pt, // Data pointers, decremented by 1
                            const int *restrict pox, const int *restrict pot, // Ordering vectors
                            const int nx, const int nt, int *restrict pres, // Sizes and result vector, pres should also be decremented by 1
                            const int second_pass) { // Is it the first or second pass?
  int i = 0, j = 0;
  SEXP tmp;
  while (i != nx && j != nt) {
    if(second_pass && pres[pox[i]] == NA_INTEGER) {
      ++i; continue; // For second passes (Need to initialize pres or add argument to function)
    }
    tmp = pt[pot[j]];
    if (px[pox[i]] == tmp) {
      pres[pox[i]] = pot[j];
      // This takes care of duplicates in x and table
      while (++i != nx && px[pox[i]] == tmp) pres[pox[i]] = pot[j];
      while (++j != nt && pt[pot[j]] == tmp);
    } else if (px[pox[i]] == NA_STRING || strcmp(CHAR(px[pox[i]]), CHAR(tmp)) < 0) {
      pres[pox[i++]] = NA_INTEGER;
    } else ++j;
  }
  if(!second_pass) while (i < nx) pres[pox[i++]] = NA_INTEGER;
}


void sort_merge_join_complex(const Rcomplex *restrict px, const Rcomplex *restrict pt, // Data pointers, decremented by 1
                             const int *restrict pox, const int *restrict pot, // Ordering vectors
                             const int nx, const int nt, int *restrict pres, // Sizes and result vector, pres should also be decremented by 1
                             const int second_pass) { // Is it the first or second pass?
  int i = 0, j = 0;
  Rcomplex xi, tj;
  while (i != nx && j != nt) {
    if(second_pass && pres[pox[i]] == NA_INTEGER) {
      ++i; continue; // For second passes (Need to initialize pres or add argument to function)
    }
    tj = pt[pot[j]];
    xi = px[pox[i]];
    if (CEQUAL(xi, tj)) {
      pres[pox[i]] = pot[j];
      // This takes care of duplicates in x and table
      while (++i != nx && CEQUAL(px[pox[i]], tj)) pres[pox[i]] = pot[j];
      while (++j != nt && CEQUAL(pt[pot[j]], tj));
    } else if (ISNAN(xi.r) || ISNAN(xi.i) || xi.r < tj.r || (xi.r == tj.r && xi.i < tj.i)) { // Todo: comparison permissible ?
      pres[pox[i++]] = NA_INTEGER;
    } else ++j;
  }
  if(!second_pass) while (i < nx) pres[pox[i++]] = NA_INTEGER;
}


SEXP sort_merge_join(SEXP x, SEXP table, SEXP ox, SEXP ot) {

  if(TYPEOF(x) != VECSXP || TYPEOF(table) != VECSXP) error("x and table need to be lists");
  if(TYPEOF(ox) != INTSXP || TYPEOF(ot) != INTSXP) error("ox and ot need to be integers");
  const int nx = length(ox), nt = length(ot), *restrict pox = INTEGER(ox), *restrict pot = INTEGER(ot);
  // TODO: x and table could be atomic??
  if(length(VECTOR_ELT(x, 0)) != nx) error("nrow(x) must match length(ox)");
  if(length(VECTOR_ELT(table, 0)) != nt) error("nrow(x) must match length(ox)");

  SEXP res = PROTECT(allocVector(INTSXP, nx));
  int *restrict pres = INTEGER(res)-1;

  SEXP clist = PROTECT(coerce_to_equal_types(x, table)); // This checks that the lengths match
  const SEXP *pc = SEXPPTR_RO(clist);
  int l = length(clist);

  for (int i = 0; i < l; ++i) {
    const SEXP *pci = SEXPPTR_RO(pc[i]);
    switch(TYPEOF(pci[0])) {
      case INTSXP:
      case LGLSXP:
        sort_merge_join_int(INTEGER(pci[0])-1, INTEGER(pci[1])-1, pox, pot, nx, nt, pres, i > 0);
        break;
      case REALSXP:
        sort_merge_join_double(REAL(pci[0])-1, REAL(pci[1])-1, pox, pot, nx, nt, pres, i > 0);
        break;
      case STRSXP:
        sort_merge_join_string(SEXPPTR_RO(pci[0])-1, SEXPPTR_RO(pci[1])-1, pox, pot, nx, nt, pres, i > 0);
        break;
      case CPLXSXP:
        sort_merge_join_complex(COMPLEX(pci[0])-1, COMPLEX(pci[1])-1, pox, pot, nx, nt, pres, i > 0);
        break;
      default:
        error("Unsupported type for x/table: %s", type2char(TYPEOF(pci[0])));
    }
  }

  UNPROTECT(2);
  return res;
}

// Miscellaneous -----------------------------------------
// SECOND PASS?

// void sort_merge_join_int_second(
//     const int *restrict px, const int *restrict pt, // Data pointers, decremented by 1
//     const int *restrict pox, const int *restrict pot, // Ordering vectors
//     // int *restrict pct, // Vector to save the counts of t (needed for second pass), also decremented by 1
//     const int nx, const int nt, int *restrict pres) {
//   int i = 0, j = 0, tmp;
//   while (i != nx && j != nt) {
//     if(pres[pox[i]] == NA_INTEGER) ++i; // If missing prom previous round, move ahead
//     tmp = pt[pot[j]];
//     if (px[pox[i]] == tmp) {
//       if(pres[pox[i]] == pot[j]) ++i; // Cool, we already got this match
//       else if (pres[pox[i]] < pot[j]) { // Now this is much more interesting: different entries in previous column of the table could have been identical...
//         pres[pox[i]] = pot[j];
//         // This takes care of duplicates in x and table
//         while (i != nx && px[pox[++i]] == tmp) pres[pox[i]] = pot[j];
//         while (j != nt && pt[pot[++j]] == tmp);
//       }
//     } else if (px[pox[i]] < tmp) {
//       pres[pox[i++]] = NA_INTEGER;
//     } else j++;
//   }
//   while (i < nx) pres[i++] = NA_INTEGER;
// }
