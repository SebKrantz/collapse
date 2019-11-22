// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// #include <R.h>
// #define USE_RINTERNALS
// #include <Rinternals.h>

NumericVector fmeanCpp(const NumericVector& x, int ng, const IntegerVector& g,
                       const SEXP& gs, const SEXP& w, bool narm);
SEXP fmeanmCpp(const NumericMatrix& x, int ng, const IntegerVector& g, const SEXP& gs,
               const SEXP& w, bool narm, bool drop);
SEXP fmeanlCpp(const List& x, int ng, const IntegerVector& g, const SEXP& gs,
               const SEXP& w, bool narm, bool drop);
