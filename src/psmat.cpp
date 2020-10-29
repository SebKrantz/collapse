// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;


template <int RTYPE>
Matrix<RTYPE> psmatCppImpl(Vector<RTYPE> x, IntegerVector g, SEXP t, bool transpose) {
  int l = x.size(), gss = g.size();
  if(gss != l) stop("length(g) must match length(x)");
  CharacterVector glevs = Rf_getAttrib(g, R_LevelsSymbol);
  int *pg = INTEGER(g);
  int ng = glevs.size(), gs = l/ng, ngp = ng+1;
  if(Rf_isNull(t)) {
    if(l%ng != 0) stop("length(x) must be a multiple of length(levels(g))");
    std::vector<int> seen(ngp);
    Matrix<RTYPE> out = transpose ? no_init_matrix(gs, ng) : no_init_matrix(ng, gs);
    if(transpose) {
      for(int i = 0; i != l; ++i) {
        if(seen[pg[i]] == gs) stop("Panel not Balanced: Need to supply timevar");
        out(seen[pg[i]]++, pg[i]-1) = x[i]; // out[(g[i]-1)*gs + seen[g[i]]++] = x[i]; not really faster...
      }
    } else {
      for(int i = 0; i != l; ++i) {
        if(seen[pg[i]] == gs) stop("Panel not Balanced: Need to supply timevar");
        out(pg[i]-1, seen[pg[i]]++) = x[i]; // out[(seen[g[i]]++)*ng + g[i]-1] = x[i]; not really faster...
      }
    }
    Rf_dimnamesgets(out, transpose ? List::create(seq_len(gs), glevs) : List::create(glevs, seq_len(gs)));
    Rf_setAttrib(out, Rf_install("transpose"), Rf_ScalarLogical(transpose));
    Rf_classgets(out, CharacterVector::create("psmat", "matrix"));
    return out;
  } else {
    int *pt = INTEGER(t);
    if(l != Rf_length(t)) stop("length(t) must match length(x)");
    // int maxt = max(t); // needed ? // check whether t.levels is same size as maxt ?
    CharacterVector tlevs = Rf_getAttrib(t, R_LevelsSymbol);
    int nt = tlevs.size();
    Matrix<RTYPE> out = transpose ? no_init_matrix(nt, ng) : no_init_matrix(ng, nt); // best way to do this ? Stable ? -> Could conditionally create vector and the coerce to matrix -> faster init ?
    if(nt != gs) std::fill(out.begin(), out.end(), Vector<RTYPE>::get_na());
    if(transpose) {
      for(int i = 0; i != l; ++i) out[(pg[i]-1)*nt + pt[i]-1] = x[i]; // out(tt[i]-1, g[i]-1) = x[i]; // tiny bit faster
    } else {
      for(int i = 0; i != l; ++i) out[(pt[i]-1)*ng + pg[i]-1] = x[i]; // out(g[i]-1, tt[i]-1) = x[i]; // tiny bit faster
    }
    Rf_dimnamesgets(out, transpose ? List::create(tlevs, glevs) : List::create(glevs, tlevs));
    Rf_setAttrib(out, Rf_install("transpose"), Rf_ScalarLogical(transpose));
    Rf_classgets(out, CharacterVector::create("psmat", "matrix"));
    return out;
  }
}

template <>
Matrix<VECSXP> psmatCppImpl(Vector<VECSXP> x, IntegerVector g, SEXP t, bool transpose) {
  stop("Not supported SEXP type!");
}

template <>
Matrix<RAWSXP> psmatCppImpl(Vector<RAWSXP> x, IntegerVector g, SEXP t, bool transpose) {
  stop("Not supported SEXP type!");
}

template <>
Matrix<EXPRSXP> psmatCppImpl(Vector<EXPRSXP> x, IntegerVector g, SEXP t, bool transpose) {
  stop("Not supported SEXP type!");
}


// [[Rcpp::export]]
SEXP psmatCpp(const SEXP& x, const IntegerVector& g, const SEXP& t = R_NilValue, bool transpose = false) {
  RCPP_RETURN_VECTOR(psmatCppImpl, x, g, t, transpose);
}


// Only Numeric Version:
// // [[Rcpp::export]]
// SEXP psmatCpp(NumericVector x, IntegerVector g, SEXP t = R_NilValue, bool transpose = false) {
//   int l = x.size(), gss = g.size();
//   if(gss != l) stop("length(g) must match length(x)");
//   CharacterVector glevs = g.attr("levels");
//   int ng = glevs.size(), gs = l/ng, ngp = ng+1;
//   if(Rf_isNull(t)) {
//     if(l%ng != 0) stop("length(x) must be a multiple of length(levels(g))");
//     IntegerVector seen(ngp);
//     NumericMatrix out = transpose ? no_init_matrix(gs, ng) : no_init_matrix(ng, gs);
//     if(transpose) {
//       for(int i = 0; i != l; ++i) {
//         if(seen[g[i]] == gs) stop("Panel not Balanced: Need to supply timevar");
//         out(seen[g[i]]++, g[i]-1) = x[i];
//       }
//     } else {
//       for(int i = 0; i != l; ++i) {
//         if(seen[g[i]] == gs) stop("Panel not Balanced: Need to supply timevar");
//         out(g[i]-1, seen[g[i]]++) = x[i];
//       }
//     }
//     out.attr("dimnames") = transpose ? List::create(seq_len(gs), glevs) : List::create(glevs, seq_len(gs));
//     return out;
//   } else {
//     IntegerVector tt = t;
//     if(l != tt.size()) stop("length(t) must match length(x)");
//     // int maxt = max(tt); // needed ?? // check whether t.levels is same size as maxt ??
//     CharacterVector tlevs = tt.attr("levels");
//     int nt = tlevs.size();
//     NumericMatrix out = transpose ? no_init_matrix(nt, ng) : no_init_matrix(ng, nt); // best way to do this ?? Stable ?? -> Could conditionally create vector and the coerce to matrix -> faster init ??
//     if(nt != gs) std::fill(out.begin(), out.end(), NA_REAL);  // memset(out, NA_REAL, sizeof(double)*ng*maxt); -> unstable !! // else balanced panel !!
//     if(transpose) {
//       for(int i = 0; i != l; ++i) out(tt[i]-1, g[i]-1) = x[i];
//     } else {
//       for(int i = 0; i != l; ++i) out(g[i]-1, tt[i]-1) = x[i];
//     }
//     out.attr("dimnames") = transpose ? List::create(tlevs, glevs) : List::create(glevs, tlevs);
//     return out;
//   }
// }
