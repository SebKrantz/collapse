#include "collapse_c.h" // Needs to be first because includes OpenMP, to avoid namespace conflicts.
#include "data.table.h"
#include "kit.h"


/* A Sort-Merge Join
 See: https://www.dcs.ed.ac.uk/home/tz/phd/thesis/node20.htm
 And: https://en.wikipedia.org/wiki/Sort-merge_join

 Note: this is only used in join(..., sort = TRUE), and expects that x was
 sorted by the join columns (done at R-level).
 The default hash join used with sort = FALSE is implemented in match.c
*/


// TODO: could add any_dup condition similar to fmatch() in while loop for j, i.e. any_dup = 1;
// this would resemble the overid argument to fmatch().

// FIRST PASS

void sort_merge_join_int(const int *restrict px, const int *restrict pt, // Data pointers, decremented by 1
                         int *restrict pg, int *restrict ptab, const int *restrict pot, // matches and ordering vector for table
                         const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, g = 0, tmp, otj;
  while (i != nx && j != nt) {
    otj = pot[j];
    tmp = pt[otj];
    if (px[i] == tmp) {
      pres[i] = otj;
      pg[i] = ptab[j] = ++g;
      // This takes care of duplicates in x and table
      while (++i != nx && px[i] == tmp) {
        pres[i] = otj;
        pg[i] = g;
      }
      while (++j != nt && pt[pot[j]] == tmp) ptab[j] = g;
    } else if ((px[i] != NA_INTEGER && px[i] < tmp) || tmp == NA_INTEGER) { // NA_INTEGER is the smallest integer: assuming ordering with na.last
      pg[i] = pres[i] = NA_INTEGER; ++i;
    } else ++j;
  }
  while(i < nx) {
    pg[i] = pres[i] = NA_INTEGER; ++i;
  }
}

void sort_merge_join_double(const double *restrict px, const double *restrict pt, // Data pointers, decremented by 1
                            int *restrict pg, int *restrict ptab, const int *restrict pot, // matches and ordering vector for table
                            const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, g = 0, otj;
  double tmp;
  while (i != nx && j != nt) {
    otj = pot[j];
    tmp = pt[otj];
    if (REQUAL(px[i], tmp)) {
      pres[i] = otj;
      pg[i] = ptab[j] = ++g;
      // This takes care of duplicates in x and table
      while (++i != nx && REQUAL(px[i], tmp)) {
        pres[i] = otj;
        pg[i] = g;
      }
      while (++j != nt && REQUAL(pt[pot[j]], tmp)) ptab[j] = g;
    } else if (px[i] < tmp || ISNAN(tmp)) { // Ordering with na.last
      pg[i] = pres[i] = NA_INTEGER; ++i;
    } else ++j;
  }
  while(i < nx) {
    pg[i] = pres[i] = NA_INTEGER; ++i;
  }
}


void sort_merge_join_string(const SEXP *restrict px, const SEXP *restrict pt, // Data pointers, decremented by 1
                            int *restrict pg, int *restrict ptab, const int *restrict pot, // matches and ordering vector for table
                            const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, g = 0, otj;
  SEXP tmp;
  while (i != nx && j != nt) {
    otj = pot[j];
    tmp = pt[otj];
    if (px[i] == tmp) {
      pres[i] = otj;
      pg[i] = ptab[j] = ++g;
      // This takes care of duplicates in x and table
      while (++i != nx && px[i] == tmp) {
        pres[i] = otj;
        pg[i] = g;
      }
      while (++j != nt && pt[pot[j]] == tmp) ptab[j] = g;
    } else if (tmp == NA_STRING || (px[i] != NA_STRING && strcmp(CHAR(px[i]), CHAR(tmp)) < 0)) { // Ordering with na.last
      pg[i] = pres[i] = NA_INTEGER; ++i;
    } else ++j;
  }
  while(i < nx) {
    pg[i] = pres[i] = NA_INTEGER; ++i;
  }
}


void sort_merge_join_complex(const Rcomplex *restrict px, const Rcomplex *restrict pt, // Data pointers, decremented by 1
                             int *restrict pg, int *restrict ptab, const int *restrict pot, // matches and ordering vector for table
                             const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, g = 0, otj;
  Rcomplex xi, tmp;
  while (i != nx && j != nt) {
    otj = pot[j];
    tmp = pt[otj];
    xi = px[i];
    if (CEQUAL(xi, tmp)) {
      pres[i] = otj;
      pg[i] = ptab[j] = ++g;
      // This takes care of duplicates in x and table
      while (++i != nx && CEQUAL(px[i], tmp)) {
        pres[i] = otj;
        pg[i] = g;
      }
      while (++j != nt && CEQUAL(pt[pot[j]], tmp)) ptab[j] = g;
    } else if (xi.r < tmp.r || (xi.r == tmp.r && xi.i < tmp.i) || ISNAN(tmp.r) || ISNAN(tmp.i)) { // Ordering with na.last
      pg[i] = pres[i] = NA_INTEGER; ++i;
    } else ++j;
  }
  while(i < nx) {
    pg[i] = pres[i] = NA_INTEGER; ++i;
  }
}



// SECOND PASS

void sort_merge_join_int_second(const int *restrict px, const int *restrict pt, // Data pointers, decremented by 1
                                int *restrict pg, int *restrict ptab, const int *restrict pot, // previous matches and ordering vector for table
                                const int nx, const int nt, int *restrict pres)  // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, g = 0, tmp, grj, otj;
  while (i != nx && j != nt) {
    if (pres[i] == NA_INTEGER) {
      ++i; continue;
    }
    grj = ptab[j];
    if (grj == 0) {
      ++j; continue;
    }
    otj = pot[j];
    tmp = pt[otj];
    if (px[i] == tmp && pg[i] == grj) {
      pres[i] = otj;
      pg[i] = ptab[j] = ++g;
      // This takes care of duplicates in x and table
      while (++i != nx && px[i] == tmp && pg[i] == grj) {
        pres[i] = otj;
        pg[i] = g;
      }
      while (++j != nt && pt[pot[j]] == tmp && ptab[j] == grj) ptab[j] = g;
    } else if (pg[i] < grj || (pg[i] == grj && ((px[i] != NA_INTEGER && px[i] < tmp) || tmp == NA_INTEGER))) { // Ordering with na.last
      pg[i] = pres[i] = NA_INTEGER; ++i;
    } else ptab[j++] = 0;
  }
  while(i < nx) {
    pg[i] = pres[i] = NA_INTEGER; ++i;
  }
}

void sort_merge_join_double_second(const double *restrict px, const double *restrict pt, // Data pointers, decremented by 1
                                   int *restrict pg, int *restrict ptab, const int *restrict pot, // previous matches and ordering vector for table
                                   const int nx, const int nt, int *restrict pres) // , int pass // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, g = 0, grj, otj;
  double tmp;
  while (i != nx && j != nt) {
    if (pres[i] == NA_INTEGER) {
      ++i; continue;
    }
    grj = ptab[j];
    if (grj == 0) {
      ++j; continue;
    }
    otj = pot[j];
    tmp = pt[otj];
    if (REQUAL(px[i], tmp) && pg[i] == grj) {
      pres[i] = otj;
      pg[i] = ptab[j] = ++g;
      // This takes care of duplicates in x and table
      while (++i != nx && REQUAL(px[i], tmp) && pg[i] == grj) {
        pres[i] = otj;
        pg[i] = g;
      }
      while (++j != nt && REQUAL(pt[pot[j]], tmp) && ptab[j] == grj) ptab[j] = g;
    } else if (pg[i] < grj || (pg[i] == grj && (px[i] < tmp || ISNAN(tmp)))) { // Ordering with na.last
      pg[i] = pres[i] = NA_INTEGER; ++i;
    } else ptab[j++] = 0;
  }
  while(i < nx) {
    pg[i] = pres[i] = NA_INTEGER; ++i;
  }
}

void sort_merge_join_string_second(const SEXP *restrict px, const SEXP *restrict pt, // Data pointers, decremented by 1
                                   int *restrict pg, int *restrict ptab, const int *restrict pot, // previous matches and ordering vector for table
                                   const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, g = 0, grj, otj;
  SEXP tmp;
  while (i != nx && j != nt) {
    if (pres[i] == NA_INTEGER) {
      ++i; continue;
    }
    grj = ptab[j];
    if (grj == 0) {
      ++j; continue;
    }
    otj = pot[j];
    tmp = pt[otj];
    if (px[i] == tmp && pg[i] == grj) {
      pres[i] = otj;
      pg[i] = ptab[j] = ++g;
      // This takes care of duplicates in x and table
      while (++i != nx && px[i] == tmp && pg[i] == grj) {
        pres[i] = otj;
        pg[i] = g;
      }
      while (++j != nt && pt[pot[j]] == tmp && ptab[j] == grj) ptab[j] = g;
    } else if (pg[i] < grj || (pg[i] == grj && (tmp == NA_STRING || (px[i] != NA_STRING && strcmp(CHAR(px[i]), CHAR(tmp)) < 0)))) { // Ordering with na.last
      pg[i] = pres[i] = NA_INTEGER; ++i;
    } else ptab[j++] = 0;
  }
  while(i < nx) {
    pg[i] = pres[i] = NA_INTEGER; ++i;
  }
}

void sort_merge_join_complex_second(const Rcomplex *restrict px, const Rcomplex *restrict pt, // Data pointers, decremented by 1
                                    int *restrict pg, int *restrict ptab, const int *restrict pot, // previous matches and ordering vector for table
                                    const int nx, const int nt, int *restrict pres) // Sizes and result vector, pres should also be decremented by 1
{
  int i = 0, j = 0, g = 0, grj, otj;
  Rcomplex tmp, xi;
  while (i != nx && j != nt) {
    if (pres[i] == NA_INTEGER) {
      ++i; continue;
    }
    grj = ptab[j];
    if (grj == 0) {
      ++j; continue;
    }
    otj = pot[j];
    tmp = pt[otj];
    xi = px[i];
    if (CEQUAL(xi, tmp) && pg[i] == grj) {
      pres[i] = otj;
      pg[i] = ptab[j] = ++g;
      // This takes care of duplicates in x and table
      while (++i != nx && CEQUAL(px[i], tmp) && pg[i] == grj) {
        pres[i] = otj;
        pg[i] = g;
      }
      while (++j != nt && CEQUAL(pt[pot[j]], tmp) && ptab[j] == grj) ptab[j] = g;
    } else if (pg[i] < grj || (pg[i] == grj && (xi.r < tmp.r || (xi.r == tmp.r && xi.i < tmp.i) || ISNAN(tmp.r) || ISNAN(tmp.i)))) { // Ordering with na.last
      pg[i] = pres[i] = NA_INTEGER; ++i;
    } else ptab[j++] = 0;
  }
  while(i < nx) {
    pg[i] = pres[i] = NA_INTEGER; ++i;
  }
}


// R FUNCTION

SEXP sort_merge_join(SEXP x, SEXP table, SEXP ot, SEXP count) {

  if(TYPEOF(x) != VECSXP || TYPEOF(table) != VECSXP) error("x and table need to be lists");
  if(TYPEOF(ot) != INTSXP) error("ot needs to be integer");
  if(length(x) == 0 || length(table) == 0) error("x and table need to have a non-zero number of columns");
  // TODO: x and table could be atomic??
  const int nx = length(VECTOR_ELT(x, 0)), nt = length(ot), *restrict pot = INTEGER(ot);
  if(length(VECTOR_ELT(table, 0)) != nt) error("nrow(table) must match length(ot)");

  SEXP res = PROTECT(allocVector(INTSXP, nx));
  int *restrict pres = INTEGER(res);
  int *pg = (int*)R_Calloc(nx, int);
  int *ptab = (int*)R_Calloc(nt, int);

  SEXP clist = PROTECT(coerce_to_equal_types(x, table)); // This checks that the lengths match
  const SEXP *pc = SEXPPTR_RO(clist);
  int l = length(clist);

  for (int i = 0; i < l; ++i) {
    const SEXP *pci = SEXPPTR_RO(pc[i]);
    switch(TYPEOF(pci[0])) {
      case INTSXP:
      case LGLSXP:
        if(i == 0) sort_merge_join_int(INTEGER_RO(pci[0]), INTEGER_RO(pci[1])-1, pg, ptab, pot, nx, nt, pres);
        else sort_merge_join_int_second(INTEGER_RO(pci[0]), INTEGER_RO(pci[1])-1, pg, ptab, pot, nx, nt, pres);
        break;
      case REALSXP:
        if(i == 0) sort_merge_join_double(REAL_RO(pci[0]), REAL_RO(pci[1])-1, pg, ptab, pot, nx, nt, pres);
        else sort_merge_join_double_second(REAL_RO(pci[0]), REAL_RO(pci[1])-1, pg, ptab, pot, nx, nt, pres);
        break;
      case STRSXP:
        if(i == 0) sort_merge_join_string(SEXPPTR_RO(PROTECT(coerceUtf8IfNeeded(pci[0]))), SEXPPTR_RO(PROTECT(coerceUtf8IfNeeded(pci[1])))-1, pg, ptab, pot, nx, nt, pres);
        else sort_merge_join_string_second(SEXPPTR_RO(PROTECT(coerceUtf8IfNeeded(pci[0]))), SEXPPTR_RO(PROTECT(coerceUtf8IfNeeded(pci[1])))-1, pg, ptab, pot, nx, nt, pres);
        UNPROTECT(2);
        break;
      case CPLXSXP:
        if(i == 0) sort_merge_join_complex(COMPLEX_RO(pci[0]), COMPLEX_RO(pci[1])-1, pg, ptab, pot, nx, nt, pres);
        else sort_merge_join_complex_second(COMPLEX_RO(pci[0]), COMPLEX_RO(pci[1])-1, pg, ptab, pot, nx, nt, pres);
        break;
      default:
        error("Unsupported type for x/table: %s", type2char(TYPEOF(pci[0])));
    }
  }

  R_Free(pg);
  R_Free(ptab);
  if(asLogical(count)) count_match(res, nt, NA_INTEGER);
  UNPROTECT(2);
  return res;
}

/*
 Helper to Perform Multi-Match Join
 The input is fmatch(x, y) and group(y, group.sizes = TRUE)
*/
SEXP multi_match(SEXP m, SEXP g) {

  SEXP ng_sym = install("N.groups"), sizes_sym = install("group.sizes"), gsR = getAttrib(g, sizes_sym);
  if(isNull(gsR)) error("Internal error: g needs to be a 'qG' type vector with a 'group.sizes' attribute.");
  const int ng = asInteger(getAttrib(g, ng_sym)), ngp = ng+1;
  if(ng != length(gsR)) error("'qG' vector is invalied, 'N.groups' attribute does not match 'group.sizes' attribute");
  const int lm = length(m), l = length(g), lp = l+1,
    *gs = INTEGER(gsR)-1, *pm = INTEGER(m), *pg = INTEGER(g)-1;

  // No multiple matches, records are unique
  if(ng == l) return m;
  int n = 0;
  #pragma omp simd reduction(+:n)
  for(int i = 0; i < lm; ++i) n += pm[i] == NA_INTEGER ? 1 : gs[pg[pm[i]]];
  if(n == lm) return m;

  // This just creates an ordering vector for g, could also use radixorder on y
  int *cgs = (int*)R_alloc(ng+2, sizeof(int)); cgs[1] = 1;
  for(int i = 1; i != ngp; ++i) cgs[i+1] = cgs[i] + gs[i];
  int *restrict cnt = (int*)R_Calloc(ngp, int);
  int *po = (int*)R_alloc(l, sizeof(int)); --po;
  for(int i = 1; i != lp; ++i) po[cgs[pg[i]] + cnt[pg[i]]++] = i;
  R_Free(cnt);

  // Indices to duplicate x
  SEXP x_ind = PROTECT(allocVector(INTSXP, n));
  // Indices to duplicate y (this is the normal fmatch(x, y) vector but now accounting for multiple matches)
  SEXP y_ind = PROTECT(allocVector(INTSXP, n));
  int *px_ind = INTEGER(x_ind), *py_ind = INTEGER(y_ind);
  for(int i = 0, j = 0, q = 0, k = 0, s = 0; i != lm; ++i) {
    if(pm[i] == NA_INTEGER) {
      px_ind[j] = i+1;
      py_ind[j++] = NA_INTEGER;
      continue;
    }
    k = pg[pm[i]];
    q = cgs[k];
    s = q + gs[k];
    while(q < s) {
      px_ind[j] = i+1;
      py_ind[j++] = po[q++];
    }
  }

  if(isObject(m)) count_match(y_ind, l, NA_INTEGER);
  // SHALLOW_DUPLICATE_ATTRIB(y_ind, m);

  SEXP res = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(res, 0, x_ind);
  SET_VECTOR_ELT(res, 1, y_ind);

  UNPROTECT(3);
  return res;
}
