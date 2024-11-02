/*
 This code is adapted from the data.table package: http://r-datatable.com
 and licensed under a Mozilla Public License 2.0 (MPL-2.0) license.
*/

#include "collapse_c.h" // Needs to be first because includes OpenMP, to avoid namespace conflicts.
#include "data.table.h"
#include <Rdefines.h>
// #include <R_ext/Rdynload.h>
// #include <R_ext/Visibility.h>

static inline int imin(int a, int b) { return a < b ? a : b; }

// global constants extern in data.table.h for gcc10 -fno-common; #4091
// these are written to once here on initialization, but because of that write they can't be declared const
SEXP char_integer64;
SEXP char_nanotime;
SEXP char_factor;
SEXP char_ordered;
SEXP char_dataframe;
SEXP char_datatable;
SEXP char_sf;
SEXP sym_sorted;
SEXP sym_index;
SEXP sym_index_df;
SEXP sym_sf_column;
SEXP SelfRefSymbol;
SEXP sym_datatable_locked;
// SEXP sym_inherits;
// SEXP sym_maxgrpn;
// SEXP sym_starts;
// SEXP char_starts;
// SEXP sym_collapse_DT_alloccol;

int max_threads;

double NA_INT64_D;
long long NA_INT64_LL;
Rcomplex NA_CPLX;
size_t sizes[100];  // max appears to be FUNSXP = 99, see Rinternals.h
size_t typeorder[100];

// -> Needed for SIZEOF macro used in rbindlist Howver TYPEORDER macro and typeof is not used...
static void setSizes(void) {
  for (int i=0; i<100; ++i) { sizes[i]=0; typeorder[i]=0; }
  // only these types are currently allowed as column types :
  sizes[LGLSXP] =  sizeof(int);       typeorder[LGLSXP] =  0;
  sizes[RAWSXP] =  sizeof(Rbyte);     typeorder[RAWSXP] =  1;
  sizes[INTSXP] =  sizeof(int);       typeorder[INTSXP] =  2;   // integer and factor
  sizes[REALSXP] = sizeof(double);    typeorder[REALSXP] = 3;   // numeric and integer64
  sizes[CPLXSXP] = sizeof(Rcomplex);  typeorder[CPLXSXP] = 4;
  sizes[STRSXP] =  sizeof(SEXP *);    typeorder[STRSXP] =  5;
  sizes[VECSXP] =  sizeof(SEXP *);    typeorder[VECSXP] =  6;   // list column
  if (sizeof(char *)>8) error("Pointers are %d bytes, greater than 8. We have not tested on any architecture greater than 64bit yet.", (int)sizeof(char *));
  // One place we need the largest sizeof is the working memory malloc in reorder.c
}
// before it was SEXP  attribute_visible
SEXP collapse_init(SEXP mess) // void SEXP mess DllInfo *info
// relies on pkg/src/Makevars to mv data.table.so to datatable.so
{
  // R_registerRoutines(info, NULL, callMethods, NULL, externalMethods);
  // R_useDynamicSymbols(info, FALSE);
  setSizes();
  const char *msg = "... failed. Please forward this message to maintainer('collapse').";
  if ((int)NA_INTEGER != (int)INT_MIN) error("Checking NA_INTEGER [%d] == INT_MIN [%d] %s", NA_INTEGER, INT_MIN, msg);
  if ((int)NA_INTEGER != (int)NA_LOGICAL) error("Checking NA_INTEGER [%d] == NA_LOGICAL [%d] %s", NA_INTEGER, NA_LOGICAL, msg);
  if (sizeof(int) != 4) error("Checking sizeof(int) [%d] is 4 %s", (int)sizeof(int), msg);
  if (sizeof(double) != 8) error("Checking sizeof(double) [%d] is 8 %s", (int)sizeof(double), msg);     // 8 on both 32bit and 64bit
  // alignof not available in C99: if (alignof(double) != 8) error("Checking alignof(double) [%d] is 8 %s", alignof(double), msg);  // 8 on both 32bit and 64bit
  if (sizeof(long long) != 8) error("Checking sizeof(long long) [%d] is 8 %s", (int)sizeof(long long), msg);
  if (sizeof(char *) != 4 && sizeof(char *) != 8) error("Checking sizeof(pointer) [%d] is 4 or 8 %s", (int)sizeof(char *), msg);
  if (sizeof(SEXP) != sizeof(char *)) error("Checking sizeof(SEXP) [%d] == sizeof(pointer) [%d] %s", (int)sizeof(SEXP), (int)sizeof(char *), msg);
  if (sizeof(uint64_t) != 8) error("Checking sizeof(uint64_t) [%d] is 8 %s", (int)sizeof(uint64_t), msg);
  if (sizeof(int64_t) != 8) error("Checking sizeof(int64_t) [%d] is 8 %s", (int)sizeof(int64_t), msg);
  if (sizeof(signed char) != 1) error("Checking sizeof(signed char) [%d] is 1 %s", (int)sizeof(signed char), msg);
  if (sizeof(int8_t) != 1) error("Checking sizeof(int8_t) [%d] is 1 %s", (int)sizeof(int8_t), msg);
  if (sizeof(uint8_t) != 1) error("Checking sizeof(uint8_t) [%d] is 1 %s", (int)sizeof(uint8_t), msg);
  if (sizeof(int16_t) != 2) error("Checking sizeof(int16_t) [%d] is 2 %s", (int)sizeof(int16_t), msg);
  if (sizeof(uint16_t) != 2) error("Checking sizeof(uint16_t) [%d] is 2 %s", (int)sizeof(uint16_t), msg);

  SEXP tmp = PROTECT(allocVector(INTSXP,2));
  if (LENGTH(tmp)!=2) error("Checking LENGTH(allocVector(INTSXP,2)) [%d] is 2 %s", LENGTH(tmp), msg);
  if (TRULEN(tmp)!=0) error("Checking TRUELENGTH(allocVector(INTSXP,2)) [%d] is 0 %s", (int)TRULEN(tmp), msg);
  UNPROTECT(1);

  // According to IEEE (http://en.wikipedia.org/wiki/IEEE_754-1985#Zero) we can rely on 0.0 being all 0 bits.
  // But check here anyway just to be sure, just in case this answer is right (http://stackoverflow.com/a/2952680/403310).
  int i = 314;
  memset(&i, 0, sizeof(int));
  if (i != 0) error("Checking memset(&i,0,sizeof(int)); i == (int)0 %s", msg);
  unsigned int ui = 314;
  memset(&ui, 0, sizeof(unsigned int));
  if (ui != 0) error("Checking memset(&ui, 0, sizeof(unsigned int)); ui == (unsigned int)0 %s", msg);
  double d = 3.14;
  memset(&d, 0, sizeof(double));
  if (d != 0.0) error("Checking memset(&d, 0, sizeof(double)); d == (double)0.0 %s", msg);
  long double ld = 3.14;
  memset(&ld, 0, sizeof(long double));
  if (ld != 0.0) error("Checking memset(&ld, 0, sizeof(long double)); ld == (long double)0.0 %s", msg);

  // Variables rather than #define for NA_INT64 to ensure correct usage; i.e. not casted
  NA_INT64_LL = LLONG_MIN;
  NA_INT64_D = LLtoD(NA_INT64_LL);
  if (NA_INT64_LL != DtoLL(NA_INT64_D)) error("Conversion of NA_INT64 via double failed %lld!=%lld", NA_INT64_LL, DtoLL(NA_INT64_D));
  // LLONG_MIN when punned to double is the sign bit set and then all zeros in exponent and significand i.e. -0.0
  //   That's why we must never test for NA_INT64_D using == in double type. Must always DtoLL and compare long long types.
  //   Assigning NA_INT64_D to a REAL is ok however.
  if (NA_INT64_D != 0.0)  error("NA_INT64_D (negative -0.0) is not == 0.0.");
  if (NA_INT64_D != -0.0) error("NA_INT64_D (negative -0.0) is not ==-0.0.");
  if (ISNAN(NA_INT64_D)) error("ISNAN(NA_INT64_D) is TRUE but should not be");
  if (isnan(NA_INT64_D)) error("isnan(NA_INT64_D) is TRUE but should not be");

  NA_CPLX.r = NA_REAL;  // NA_REAL is defined as R_NaReal which is not a strict constant and thus initializer {NA_REAL, NA_REAL} can't be used in .h
  NA_CPLX.i = NA_REAL;  // https://github.com/Rdatatable/data.table/pull/3689/files#r304117234

  // create needed strings in advance for speed, same techique as R_*Symbol
  // Following R-exts 5.9.4; paragraph and example starting "Using install ..."
  // either use PRINTNAME(install()) or R_PreserveObject(mkChar()) here.
  char_integer64 = PRINTNAME(install("integer64"));
  char_nanotime =  PRINTNAME(install("nanotime"));
  // char_starts =    PRINTNAME(sym_starts = install("starts"));
  char_factor =    PRINTNAME(install("factor"));
  char_ordered =   PRINTNAME(install("ordered"));
  char_dataframe = PRINTNAME(install("data.frame"));
  char_datatable = PRINTNAME(install("data.table"));
  char_sf = PRINTNAME(install("sf"));

  if (TYPEOF(char_integer64) != CHARSXP) {
    // checking one is enough in case of any R-devel changes
    error("PRINTNAME(install(\"integer64\")) has returned %s not %s", type2char(TYPEOF(char_integer64)), type2char(CHARSXP));  // # nocov
  }

  // create commonly used symbols, same as R_*Symbol but internal to DT
  // Not really for speed but to avoid leak in situations like setAttrib(DT, install(), allocVector()) where
  // the allocVector() can happen first and then the install() could gc and free it before it is protected
  // within setAttrib. Thanks to Bill Dunlap finding and reporting. Using these symbols instead of install()
  // avoids the gc without needing an extra PROTECT and immediate UNPROTECT after the setAttrib which would
  // look odd (and devs in future might be tempted to remove them). Avoiding passing install() to API calls
  // keeps the code neat and readable. Also see grep's added to CRAN_Release.cmd to find such calls.

  sym_sorted  = install("sorted");
  sym_index   = install("index");
  sym_index_df   = install("index_df");
  sym_sf_column = install("sf_column");
  SelfRefSymbol = install(".internal.selfref");
  sym_datatable_locked = install(".data.table.locked");
  // sym_inherits = install("inherits");
  // sym_maxgrpn = install("maxgrpn");
  // sym_collapse_DT_alloccol = install("collapse_DT_alloccol");

  max_threads = OMP_NUM_PROCS;
  max_threads = imin(max_threads, OMP_THREAD_LIMIT);
  max_threads = imin(max_threads, OMP_MAX_THREADS);

  return mess;
}


inline long long DtoLL(double x) {
  union {double d; int64_t i64;} u;
  u.d = x;
  return (long long)u.i64;
}

inline double LLtoD(long long x) {
  union {double d; int64_t i64;} u;
  u.i64 = (int64_t)x;
  return u.d;
}

