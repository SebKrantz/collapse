#include "collapse_c.h" // Needs to be first because includes OpenMP, to avoid namespace conflicts.
#include "kit.h"
#include "base_radixsort.h"

/* A Sort-Merge Join
 See: https://www.dcs.ed.ac.uk/home/tz/phd/thesis/node20.htm
 And: https://en.wikipedia.org/wiki/Sort-merge_join

 Note: this is only used in join(..., sort = TRUE), and expects that x was
 sorted by the join columns. The default hash join used with sort = FALSE
 is implemented in match.c
*/


// FIRST PASS

void sort_merge_join_int(const int *restrict px, const int *restrict pt, // Data pointers, decremented by 1
                         int *restrict ptab, const int *restrict pot, // matches and ordering vector for table
                         const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, tmp, otj;
  while (i != nx && j != nt) {
    otj = pot[j];
    tmp = pt[otj];
    if (px[i] == tmp) {
      ptab[otj] = pres[i] = otj;
      // This takes care of duplicates in x and table
      while (++i != nx && px[i] == tmp) pres[i] = otj;
      while (++j != nt && pt[pot[j]] == tmp) ptab[pot[j]] = otj;
    } else if ((px[i] != NA_INTEGER && px[i] < tmp) || tmp == NA_INTEGER) { // NA_INTEGER is the smallest integer: assuming ordering with na.last
      pres[i++] = NA_INTEGER;
    } else ++j;
  }
  while (i < nx) pres[i++] = NA_INTEGER;
}

void sort_merge_join_double(const double *restrict px, const double *restrict pt, // Data pointers, decremented by 1
                            int *restrict ptab, const int *restrict pot, // matches and ordering vector for table
                            const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, otj;
  double tmp;
  while (i != nx && j != nt) {
    otj = pot[j];
    tmp = pt[otj];
    if (REQUAL(px[i], tmp)) {
      ptab[otj] = pres[i] = otj;
      // This takes care of duplicates in x and table
      while (++i != nx && REQUAL(px[i], tmp)) pres[i] = otj;
      while (++j != nt && REQUAL(pt[pot[j]], tmp)) ptab[pot[j]] = otj;
    } else if (px[i] < tmp || ISNAN(tmp)) { // ISNAN(px[i]) || // Ordering with na.last: allows NA's to be matched last...
      pres[i++] = NA_INTEGER;
    } else ++j;
  }
  while (i < nx) pres[i++] = NA_INTEGER;
}


void sort_merge_join_string(const SEXP *restrict px, const SEXP *restrict pt, // Data pointers, decremented by 1
                            int *restrict ptab, const int *restrict pot, // matches and ordering vector for table
                            const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, otj;
  SEXP tmp;
  while (i != nx && j != nt) {
    otj = pot[j];
    tmp = pt[otj];
    if (px[i] == tmp) {
      ptab[otj] = pres[i] = otj;
      // This takes care of duplicates in x and table
      while (++i != nx && px[i] == tmp) pres[i] = otj;
      while (++j != nt && pt[pot[j]] == tmp) ptab[pot[j]] = otj;
    } else if (tmp == NA_STRING || (px[i] != NA_STRING && strcmp(CHAR(px[i]), CHAR(tmp)) < 0)) { // Ordering with na.last
      pres[i++] = NA_INTEGER;
    } else ++j;
  }
  while (i < nx) pres[i++] = NA_INTEGER;
}


void sort_merge_join_complex(const Rcomplex *restrict px, const Rcomplex *restrict pt, // Data pointers, decremented by 1
                             int *restrict ptab, const int *restrict pot, // matches and ordering vector for table
                             const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, otj;
  Rcomplex xi, tmp;
  while (i != nx && j != nt) {
    otj = pot[j];
    tmp = pt[otj];
    xi = px[i];
    if (CEQUAL(xi, tmp)) {
      ptab[otj] = pres[i] = otj;
      // This takes care of duplicates in x and table
      while (++i != nx && CEQUAL(px[i], tmp)) pres[i] = otj;
      while (++j != nt && CEQUAL(pt[pot[j]], tmp)) ptab[pot[j]] = otj;
    } else if (xi.r < tmp.r || (xi.r == tmp.r && xi.i < tmp.i) || ISNAN(tmp.r) || ISNAN(tmp.i)) { // ISNAN(xi.r) || ISNAN(xi.i) || Ordering with na.last
      pres[i++] = NA_INTEGER;
    } else ++j;
  }
  while (i < nx) pres[i++] = NA_INTEGER;
}



// SECOND PASS

void sort_merge_join_int_second(const int *restrict px, const int *restrict pt, // Data pointers, decremented by 1
                                int *restrict ptab, const int *restrict pot, // previous matches and ordering vector for table
                                const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, tmp, grj, otj;
  while (i != nx && j != nt) {
    if (pres[i] == NA_INTEGER) {
      ++i; continue;
    }
    // if (pres[i] > pot[j]) j += pres[i] - pot[j]; // Need to keep track of previous matching, if we are in a different category 1, move ahead an appropriate number of steps
    otj = pot[j];
    grj = ptab[otj];
    tmp = pt[otj];
    if (px[i] == tmp && pres[i] == grj) {
      ptab[otj] = pres[i] = otj;
      while (++i != nx && px[i] == tmp && pres[i] == grj) pres[i] = otj;
      while (++j != nt && pt[pot[j]] == tmp && ptab[pot[j]] == grj) ptab[pot[j]] = otj;
    } else if (pres[i] < grj || (pres[i] == grj && (tmp == NA_INTEGER || (px[i] != NA_INTEGER && px[i] < tmp)))) { // Assuming ordering with na.last
      pres[i++] = NA_INTEGER;
    } else ++j;
  }
  while (i < nx) pres[i++] = NA_INTEGER;
}

void sort_merge_join_double_second(const double *restrict px, const double *restrict pt, // Data pointers, decremented by 1
                                   int *restrict ptab, const int *restrict pot, // previous matches and ordering vector for table
                                   const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, grj, otj;
  double tmp;
  while (i != nx && j != nt) {
    if (pres[i] == NA_INTEGER) {
      ++i; continue;
    }
    otj = pot[j];
    grj = ptab[otj];
    tmp = pt[otj];
    if (REQUAL(px[i], tmp) && pres[i] == grj) {
      ptab[otj] = pres[i] = otj;
      // This takes care of duplicates in x and table
      while (++i != nx && REQUAL(px[i], tmp) && pres[i] == grj) pres[i] = otj;
      while (++j != nt && REQUAL(pt[pot[j]], tmp) && ptab[pot[j]] == grj) ptab[pot[j]] = otj;
    } else if (pres[i] < grj || (pres[i] == grj && (ISNAN(tmp) || px[i] < tmp))) { // ISNAN(px[i]) || px[i] < tmp || (REQUAL(px[i], tmp) && pres[i] < grj)
      pres[i++] = NA_INTEGER;
    } else ++j;
  }
  while (i < nx) pres[i++] = NA_INTEGER;
}

void sort_merge_join_string_second(const SEXP *restrict px, const SEXP *restrict pt, // Data pointers, decremented by 1
                                   int *restrict ptab, const int *restrict pot, // previous matches and ordering vector for table
                                   const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, grj, otj;
  SEXP tmp;
  while (i != nx && j != nt) {
    if (pres[i] == NA_INTEGER) {
	    i++; continue;
    }
    otj = pot[j];
    grj = ptab[otj];
    tmp = pt[otj];
    if (px[i] == tmp && pres[i] == grj) {
      ptab[otj] = pres[i] = otj;
      while (++i != nx && px[i] == tmp && pres[i] == grj) pres[i] = otj;
      while (++j != nt && pt[pot[j]] == tmp && ptab[pot[j]] == grj) ptab[pot[j]] = otj;
    } else if (pres[i] < grj || (pres[i] == grj && (tmp == NA_STRING || (px[i] != NA_STRING && strcmp(CHAR(px[i]), CHAR(tmp)) < 0)))) { // px[i] == NA_STRING || strcmp(CHAR(px[i]), CHAR(tmp)) < 0 || (px[i] == tmp && pres[i] < grj)
      pres[i++] = NA_INTEGER;
    } else ++j;
  }
  while (i < nx) pres[i++] = NA_INTEGER;
}

void sort_merge_join_complex_second(const Rcomplex *restrict px, const Rcomplex *restrict pt, // Data pointers, decremented by 1
                                    int *restrict ptab, const int *restrict pot, // previous matches and ordering vector for table
                                    const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, grj, otj;
  Rcomplex tmp, xi;
  while (i != nx && j != nt) {
    if (pres[i] == NA_INTEGER) {
	    i++; continue;
    }
    otj = pot[j];
    grj = ptab[otj];
    tmp = pt[otj];
    xi = px[i];
    if (CEQUAL(xi, tmp) && pres[i] == grj) {
      ptab[otj] = pres[i] = otj;
      while (++i != nx && CEQUAL(px[i], tmp) && pres[i] == grj) pres[i] = otj;
      while (++j != nt && CEQUAL(pt[pot[j]], tmp) && ptab[pot[j]] == grj) ptab[pot[j]] = otj;
    } else if (pres[i] < grj || (pres[i] == grj && (xi.r < tmp.r || (xi.r == tmp.r && xi.i < tmp.i) || ISNAN(tmp.r) || ISNAN(tmp.i)))) { // ISNAN(xi.r) || ISNAN(xi.i) || xi.r < tmp.r || (xi.r == tmp.r && xi.i < tmp.i) || (CEQUAL(xi, tmp) && pres[i] < grj)
      pres[i++] = NA_INTEGER;
    } else ++j;
  }
  while (i < nx) pres[i++] = NA_INTEGER;
}



// R FUNCTION

SEXP sort_merge_join(SEXP x, SEXP table, SEXP ot, SEXP count) {

  if(TYPEOF(x) != VECSXP || TYPEOF(table) != VECSXP) error("x and table need to be lists");
  if(TYPEOF(ot) != INTSXP) error("ot needs to be integer");
  // TODO: x and table could be atomic??
  const int nx = length(VECTOR_ELT(x, 0)), nt = length(ot), *restrict pot = INTEGER(ot);
  if(length(VECTOR_ELT(table, 0)) != nt) error("nrow(table) must match length(ot)");

  SEXP res = PROTECT(allocVector(INTSXP, nx));
  int *restrict pres = INTEGER(res);
  int *ptab = (int*)Calloc(nt+1, int);

  SEXP clist = PROTECT(coerce_to_equal_types(x, table)); // This checks that the lengths match
  const SEXP *pc = SEXPPTR_RO(clist);
  int l = length(clist);

  for (int i = 0; i < l; ++i) {
    const SEXP *pci = SEXPPTR_RO(pc[i]);
    switch(TYPEOF(pci[0])) {
      case INTSXP:
      case LGLSXP:
        if(i == 0) sort_merge_join_int(INTEGER_RO(pci[0]), INTEGER_RO(pci[1])-1, ptab, pot, nx, nt, pres);
        else sort_merge_join_int_second(INTEGER_RO(pci[0]), INTEGER_RO(pci[1])-1, ptab, pot, nx, nt, pres);
        break;
      case REALSXP:
        if(i == 0) sort_merge_join_double(REAL_RO(pci[0]), REAL_RO(pci[1])-1, ptab, pot, nx, nt, pres);
        else sort_merge_join_double_second(REAL_RO(pci[0]), REAL_RO(pci[1])-1, ptab, pot, nx, nt, pres);
        break;
      case STRSXP:
        checkEncodings(pci[0]); checkEncodings(pci[1]);
        if(i == 0) sort_merge_join_string(STRING_PTR_RO(pci[0]), STRING_PTR_RO(pci[1])-1, ptab, pot, nx, nt, pres);
        else sort_merge_join_string_second(STRING_PTR_RO(pci[0]), STRING_PTR_RO(pci[1])-1, ptab, pot, nx, nt, pres);
        break;
      case CPLXSXP:
        if(i == 0) sort_merge_join_complex(COMPLEX_RO(pci[0]), COMPLEX_RO(pci[1])-1, ptab, pot, nx, nt, pres);
        else sort_merge_join_complex_second(COMPLEX_RO(pci[0]), COMPLEX_RO(pci[1])-1, ptab, pot, nx, nt, pres);
        break;
      default:
        error("Unsupported type for x/table: %s", type2char(TYPEOF(pci[0])));
    }
  }

  Free(ptab);
  if(asLogical(count)) count_match(res, nt, NA_INTEGER);
  UNPROTECT(2);
  return res;
}

