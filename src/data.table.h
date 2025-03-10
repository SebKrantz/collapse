/*
 This code is adapted from the data.table package: http://r-datatable.com
 and licensed under a Mozilla Public License 2.0 (MPL-2.0) license.
*/

#ifndef DATATABLE_H  // Check if DATATABLE_H is not defined
#define DATATABLE_H  // Define DATATABLE_H

// #define USE_RINTERNALS
#include "base_radixsort.h"
// #include <stdint.h> // for uint64_t rather than unsigned long long
#include <stdbool.h>
// #include "types.h"

#define IS_TRUE(x)  (TYPEOF(x)==LGLSXP && LENGTH(x)==1 && LOGICAL(x)[0]==TRUE)
#define IS_FALSE(x) (TYPEOF(x)==LGLSXP && LENGTH(x)==1 && LOGICAL(x)[0]==FALSE)
#define IS_TRUE_OR_FALSE(x) (TYPEOF(x)==LGLSXP && LENGTH(x)==1 && LOGICAL(x)[0]!=NA_LOGICAL)
#define SIZEOF(x) sizes[TYPEOF(x)]
#define TYPEORDER(x) typeorder[x]

// Needed for match.c and join.c
#define NEED2UTF8(s) !(IS_ASCII(s) || (s)==NA_STRING || IS_UTF8(s))
#define ENC2UTF8(s) (!NEED2UTF8(s) ? (s) : mkCharCE(translateCharUTF8(s), CE_UTF8))

// for use with bit64::integer64
#define NA_INTEGER64  INT64_MIN
#define MAX_INTEGER64 INT64_MAX

// init.c // https://stackoverflow.com/questions/1410563/what-is-the-difference-between-a-definition-and-a-declaration
extern SEXP char_integer64;
extern SEXP char_nanotime;
extern SEXP char_factor;
extern SEXP char_ordered;
extern SEXP char_dataframe;
extern SEXP char_datatable;
extern SEXP char_sf;
extern SEXP sym_sorted;
extern SEXP sym_index;
extern SEXP sym_index_df;
extern SEXP sym_sf_column;
extern SEXP SelfRefSymbol;
extern SEXP sym_datatable_locked;

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
int need2utf8(SEXP x);
SEXP coerceUtf8IfNeeded(SEXP x);
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

#endif // End of DATATABLE_H guard
