// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;


template <int RTYPE>
Matrix<RTYPE> psmatCppImpl(Vector<RTYPE> x, IntegerVector g, SEXP t, bool transpose) {
  int l = x.size(), gss = g.size();
  if(gss != l) stop("length(g) must match length(x)");
  CharacterVector glevs = g.attr("levels");
  int ng = glevs.size(), gs = l/ng, ngp = ng+1;
  if(Rf_isNull(t)) {
    if(l%ng != 0) stop("length(x) must be a multiple of length(levels(g))");
    IntegerVector seen(ngp);
    Matrix<RTYPE> out = transpose ? no_init_matrix(gs, ng) : no_init_matrix(ng, gs);
    // List dim = ATTRIB(out);
    // SHALLOW_DUPLICATE_ATTRIB(out, x); // good ?
    // Rf_setAttrib(out, R_DimSymbol, dim[0]);
    if(transpose) {
      for(int i = 0; i != l; ++i) {
        if(seen[g[i]] == gs) stop("Panel not Balanced: Need to supply timevar");
        out(seen[g[i]]++, g[i]-1) = x[i]; // out[(g[i]-1)*gs + seen[g[i]]++] = x[i]; not really faster...
      }
    } else {
      for(int i = 0; i != l; ++i) {
        if(seen[g[i]] == gs) stop("Panel not Balanced: Need to supply timevar");
        out(g[i]-1, seen[g[i]]++) = x[i]; // out[(seen[g[i]]++)*ng + g[i]-1] = x[i]; not really faster...
      }
    }
    out.attr("dimnames") = transpose ? List::create(seq_len(gs), glevs) : List::create(glevs, seq_len(gs));
    out.attr("transpose") = transpose;
    out.attr("class") = CharacterVector::create("psmat","matrix");
    return out;
  } else {
    IntegerVector tt = t;
    if(l != tt.size()) stop("length(t) must match length(x)");
    // int maxt = max(tt); // needed ? // check whether t.levels is same size as maxt ?
    CharacterVector tlevs = tt.attr("levels");
    int nt = tlevs.size();
    Matrix<RTYPE> out = transpose ? no_init_matrix(nt, ng) : no_init_matrix(ng, nt); // best way to do this ? Stable ? -> Could conditionally create vector and the coerce to matrix -> faster init ?
    if(nt != gs) std::fill(out.begin(), out.end(), Vector<RTYPE>::get_na());  // memset(out, NA_REAL, sizeof(double)*ng*maxt); -> unstable ! // else balanced panel !
    // List dim = ATTRIB(out);
    // SHALLOW_DUPLICATE_ATTRIB(out, x); // good ?
    // Rf_setAttrib(out, R_DimSymbol, dim[0]);
    if(transpose) {
      for(int i = 0; i != l; ++i) out[(g[i]-1)*nt + tt[i]-1] = x[i]; // out(tt[i]-1, g[i]-1) = x[i]; // tiny bit faster
    } else {
      for(int i = 0; i != l; ++i) out[(tt[i]-1)*ng + g[i]-1] = x[i]; // out(g[i]-1, tt[i]-1) = x[i]; // tiny bit faster
    }
    out.attr("dimnames") = transpose ? List::create(tlevs, glevs) : List::create(glevs, tlevs);
    out.attr("transpose") = transpose;
    out.attr("class") = CharacterVector::create("psmat","matrix");
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
SEXP psmatCpp(SEXP x, IntegerVector g, SEXP t = R_NilValue, bool transpose = false) {
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
