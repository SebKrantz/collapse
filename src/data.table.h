/*
 This code is adapted from the data.table package: http://r-datatable.com
 and licensed under a Mozilla Public License 2.0 (MPL-2.0) license.
*/

#include <R.h>
#define USE_RINTERNALS
#include <Rinternals.h>
#include <stdint.h> // for uint64_t rather than unsigned long long
#include <stdbool.h>
// #include "types.h"

// data.table depends on R>=3.0.0 when R_xlen_t was introduced
// Before R 3.0.0, RLEN used to be switched to R_len_t as R_xlen_t wasn't available.
// We could now replace all RLEN with R_xlen_t directly. Or keep RLEN for the shorter
// name so as not to have to check closely one letter difference R_xlen_t/R_len_t. We
// might also undefine R_len_t to ensure not to use it.
typedef R_xlen_t RLEN;

// #define PRId64 "lld" // needed in rbindlist CHECK_RANGE macro -> disabled right now..
#define IS_ASCII(x) (LEVELS(x) & 64)
#define IS_TRUE(x)  (TYPEOF(x)==LGLSXP && LENGTH(x)==1 && LOGICAL(x)[0]==TRUE)
#define IS_FALSE(x) (TYPEOF(x)==LGLSXP && LENGTH(x)==1 && LOGICAL(x)[0]==FALSE)
#define IS_TRUE_OR_FALSE(x) (TYPEOF(x)==LGLSXP && LENGTH(x)==1 && LOGICAL(x)[0]!=NA_LOGICAL)
#define SIZEOF(x) sizes[TYPEOF(x)]
#define TYPEORDER(x) typeorder[x]

// for use with bit64::integer64
#define NA_INTEGER64  INT64_MIN
#define MAX_INTEGER64 INT64_MAX

// Backport macros added to R in 2017 so we don't need to update dependency from R 3.0.0
#ifndef MAYBE_REFERENCED
# define MAYBE_REFERENCED(x) ( NAMED(x) > 0 )
#endif

#ifndef ALTREP
#define ALTREP(x) 0  // for R<3.5.0, see issue #2866 and grep for "ALTREP" to see comments where it's used
#endif

#define SEXPPTR(x) ((SEXP *)DATAPTR(x))  // to avoid overhead of looped STRING_ELT and VECTOR_ELT
#define SEXPPTR_RO(x) ((const SEXP *)DATAPTR_RO(x))  // to avoid overhead of looped VECTOR_ELT

// init.c // https://stackoverflow.com/questions/1410563/what-is-the-difference-between-a-definition-and-a-declaration
extern SEXP char_integer64;
extern SEXP char_nanotime;
extern SEXP char_factor;
extern SEXP char_ordered;
extern SEXP char_dataframe;
extern SEXP char_datatable;
extern SEXP char_sf;
// not currently needed (base_radixsort uses install), but perhaps later..
extern SEXP sym_sorted;
// extern SEXP sym_maxgrpn;
// extern SEXP sym_starts;
// extern SEXP char_starts;
extern SEXP sym_index;
extern SEXP sym_index_df;
extern SEXP sym_inherits;
extern SEXP sym_sf_column;
extern SEXP SelfRefSymbol;
extern SEXP sym_datatable_locked;
// extern SEXP sym_collapse_DT_alloccol;

// data.table_init.c
SEXP collapse_init(SEXP mess);
long long DtoLL(double x);
double LLtoD(long long x);
extern double NA_INT64_D;
extern long long NA_INT64_LL;
extern Rcomplex NA_CPLX;  // initialized in init.c; see there for comments
extern size_t sizes[100];  // max appears to be FUNSXP = 99, see Rinternals.h
extern size_t typeorder[100];

// data.table_utils.c
SEXP setnames(SEXP x, SEXP nam);
bool allNA(SEXP x, bool errorForBadType);
SEXP allNAv(SEXP x, SEXP errorForBadType);
bool INHERITS(SEXP x, SEXP char_);
SEXP dt_na(SEXP x, SEXP cols, SEXP Rprop, SEXP Rcount);
SEXP frankds(SEXP xorderArg, SEXP xstartArg, SEXP xlenArg, SEXP dns);
SEXP setcolorder(SEXP x, SEXP o);

// data.table_subset.c
void setselfref(SEXP x);
SEXP Calloccol(SEXP dt);
SEXP convertNegAndZeroIdx(SEXP idx, SEXP maxArg, SEXP allowOverMax);
SEXP extendIntVec(SEXP x, int len, int val);
SEXP subsetCols(SEXP x, SEXP cols, SEXP checksf);
SEXP subsetDT(SEXP x, SEXP rows, SEXP cols, SEXP checkrows);
SEXP subsetVector(SEXP x, SEXP idx, SEXP checkidx);

// rbindlist.c
void writeNA(SEXP v, const int from, const int n);
void writeValue(SEXP target, SEXP source, const int from, const int n);
void savetl_init(void), savetl(SEXP s), savetl_end(void);
SEXP rbindlist(SEXP l, SEXP usenamesArg, SEXP fillArg, SEXP idcolArg);

