
// #include <Defn.h> // Not available in C API !!
// #include <Internal.h> // Not available in C API !!
#include <R.h>
#include <Rinternals.h>
#include <stdint.h>
// typedef uint64_t ZPOS64_T; // already defined in stdint.h
#define IS_ASCII(x) (LEVELS(x) & 64) // from data.table.h
// #define ASCII_MASK (1<<6) // evaluates to 64 !!
// # define IS_ASCII(x) ((x)->sxpinfo.gp & ASCII_MASK)
// #define IS_ASCII(x) (LEVELS(x) & ASCII_MASK)

/* It would be better to find a way to avoid abusing TRUELENGTH, but
 in the meantime replace TRUELENGTH/SET_TRUELENGTH with
 TRLEN/SET_TRLEN that cast to int to avoid warnings. */
#define TRLEN(x) ((int) TRUELENGTH(x))
#define SET_TRLEN(x, v) SET_TRUELENGTH(x, ((int) (v)))


SEXP Cradixsort(SEXP NA_last, SEXP decreasing, SEXP RETstrt, SEXP RETgs, SEXP SORTStr, SEXP args);
void Cdoubleradixsort(int *o, Rboolean NA_last, Rboolean decreasing, SEXP x);
