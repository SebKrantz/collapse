// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Implemented smarter copy names !!

template <int RTYPE>
inline bool isnaNUM(typename Rcpp::traits::storage_type<RTYPE>::type x) {
  return x != x;
}

template <int RTYPE>
inline bool isnaOTH(typename Rcpp::traits::storage_type<RTYPE>::type x) {
  return x == Vector<RTYPE>::get_na();
}

template <int RTYPE>
Vector<RTYPE> ffirstCppImpl(const Vector<RTYPE>& x, int ng, const IntegerVector& g, bool narm) {
  int l = x.size(), end = l-1;
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;
  if (ng == 0) {
    if(narm) {
      int j = 0;
      auto first = x[j];
      while(isnanT(first) && j!=end) first = x[++j];
      Vector<RTYPE> out(1, first); // faster using create ??
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        out.attr("names") = CharacterVector::create(names[j]);
      }
      return out;
    } else {
      Vector<RTYPE> out(1, x[0]);
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        out.attr("names") = CharacterVector::create(names[0]);
      }
      return out;
    }
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    int ngs = 0;
    Vector<RTYPE> first = no_init_vector(ng);
    DUPLICATE_ATTRIB(first, x);
    if(narm) {
      std::fill(first.begin(), first.end(), Vector<RTYPE>::get_na());
      if(Rf_getAttrib(x, R_NamesSymbol) == R_NilValue) {
        for(int i = 0; i != l; ++i) {
          if(!isnanT(x[i])) {
            if(isnanT(first[g[i]-1])) {
              first[g[i]-1] = x[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
      } else {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        if(names.size() != l) stop("x has a names attribute of length != length(x)");
        CharacterVector newnames = no_init_vector(ng);
        for(int i = 0; i != l; ++i) {
          if(!isnanT(x[i])) {
            if(isnanT(first[g[i]-1])) {
              first[g[i]-1] = x[i];
              newnames[g[i]-1] = names[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
        first.attr("names") = newnames;
      }
    } else {
      LogicalVector gl(ng, true); // std::vector<bool> glj(ng, true);? -> Nope, not faster (see matrix method)
      if(Rf_getAttrib(x, R_NamesSymbol) == R_NilValue) {
        for(int i = 0; i != l; ++i) {
          if(gl[g[i]-1]) {
            gl[g[i]-1] = false;
            first[g[i]-1] = x[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
      } else {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        if(names.size() != l) stop("x has a names attribute of length != length(x)");
        CharacterVector newnames = no_init_vector(ng);
        for(int i = 0; i != l; ++i) {
          if(gl[g[i]-1]) {
            gl[g[i]-1] = false;
            first[g[i]-1] = x[i];
            newnames[g[i]-1] = names[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
        first.attr("names") = newnames;
      }
    }
    return first;
  }
}

template <>
Vector<CPLXSXP> ffirstCppImpl(const Vector<CPLXSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<VECSXP> ffirstCppImpl(const Vector<VECSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<RAWSXP> ffirstCppImpl(const Vector<RAWSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<EXPRSXP> ffirstCppImpl(const Vector<EXPRSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}


// [[Rcpp::export]]
SEXP ffirstCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true){
  RCPP_RETURN_VECTOR(ffirstCppImpl, x, ng, g, narm);
}



// with better handling row.names !!

template <int RTYPE>
SEXP ffirstmCppImpl(const Matrix<RTYPE>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  int l = x.nrow(), col = x.ncol();
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;

  if(ng == 0) {
    Vector<RTYPE> first = no_init_vector(col);
    if(narm) {
      int end = l-1;
      for(int j = col; j--; ) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        int k = 0;
        auto firstj = column[k];
        while(isnanT(firstj) && k!=end) firstj = column[++k];
        first[j] = firstj;
      }
    } else {
      first = x(0, _);
    }
    if(drop) first.attr("names") = colnames(x);
    else {
      first.attr("dim") = Dimension(1, col);
      colnames(first) = colnames(x);
    }
    return first;
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    Matrix<RTYPE> first = no_init_matrix(ng, col);
    if(narm) {
      std::fill(first.begin(), first.end(), Vector<RTYPE>::get_na());
      for(int j = col; j--; ) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        MatrixColumn<RTYPE> firstj = first( _ , j);
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(!isnanT(column[i])) {
            if(isnanT(firstj[g[i]-1])) {
              firstj[g[i]-1] = column[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
      }
      colnames(first) = colnames(x);
    } else {
      List dn = x.attr("dimnames");
      if(dn[0] != R_NilValue) {
        CharacterVector rn = dn[0];
        CharacterVector newrn = no_init_vector(ng);
        LogicalVector glj(ng, true); // using std::vector<bool> here is more memory efficient but not faster !!
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(glj[g[i]-1]) {
            glj[g[i]-1] = false;
            first(g[i]-1, _) = x(i, _);
            newrn[g[i]-1] = rn[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
        first.attr("dimnames") = List::create(newrn, dn[1]); // best way !!
      } else {
        LogicalVector glj(ng, true); // using std::vector<bool> here is more memory efficient but not faster !!
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(glj[g[i]-1]) {
            glj[g[i]-1] = false;
            first(g[i]-1, _) = x(i, _);
            ++ngs;
            if(ngs == ng) break;
          }
        }
        colnames(first) = colnames(x);
      }
    }
    return first;
  }
}

template <>
SEXP ffirstmCppImpl(const Matrix<CPLXSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP ffirstmCppImpl(const Matrix<VECSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP ffirstmCppImpl(const Matrix<RAWSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP ffirstmCppImpl(const Matrix<EXPRSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP ffirstmCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true, bool drop = true){
  RCPP_RETURN_MATRIX(ffirstmCppImpl, x, ng, g, narm, drop);
}




// [[Rcpp::export]]
SEXP ffirstlCpp(const List& x, int ng = 0, const IntegerVector& g = 0, bool narm = true) { // , bool drop = true

  int l = x.size();
  List first(l);

  if (ng == 0) {
    if(narm) {
      for(int j = l; j--; ) {
        int k = 0;
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          int row = column.size(), end = row-1;
          while(std::isnan(column[k]) && k!=end) ++k;
          NumericVector out(1, column[k]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          first[j] = out;
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          int row = column.size(), end = row-1;
          while(column[k] == NA_INTEGER && k!=end) ++k;
          IntegerVector out(1, column[k]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          first[j] = out;
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          int row = column.size(), end = row-1;
          while(column[k] == NA_STRING && k!=end) ++k;
          CharacterVector out(1, column[k]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          first[j] = out;
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          int row = column.size(), end = row-1;
          while(column[k] == NA_LOGICAL && k!=end) ++k;
          LogicalVector out(1, column[k]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          first[j] = out;
          break;
        }
        default:
          stop("incompatible SEXP encountered;");
          break;
        }
      }
    } else {
      for(int j = l; j--; ) {
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          NumericVector out(1, column[0]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          first[j] = out;
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          IntegerVector out(1, column[0]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          first[j] = out;
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          CharacterVector out(1, column[0]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          first[j] = out;
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          LogicalVector out(1, column[0]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          first[j] = out;
          break;
        }
        default:
          stop("incompatible SEXP encountered;");
          break;
        }
      }
    }
    // if(drop) first.attr("names") = x.attr("names");
    DUPLICATE_ATTRIB(first, x);
    first.attr("row.names") = 1;
    return first;
  } else { // With groups !!
    int gss = g.size();
    if(narm) {
      for(int j = l; j--; ) {
        int ngs = 0;
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          NumericVector firstj(ng, NA_REAL);
          for(int i = 0; i != gss; ++i) {
            if(!std::isnan(column[i])) {
              if(std::isnan(firstj[g[i]-1])) {
                firstj[g[i]-1] = column[i];
                ++ngs;
                if(ngs == ng) break;
              }
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(firstj, column);
          first[j] = firstj;
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          IntegerVector firstj(ng, NA_INTEGER);
          for(int i = 0; i != gss; ++i) {
            if(column[i] != NA_INTEGER) {
              if(firstj[g[i]-1] == NA_INTEGER) {
                firstj[g[i]-1] = column[i];
                ++ngs;
                if(ngs == ng) break;
              }
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(firstj, column);
          first[j] = firstj;
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          CharacterVector firstj(ng, NA_STRING);
          for(int i = 0; i != gss; ++i) {
            if(column[i] != NA_STRING) {
              if(firstj[g[i]-1] == NA_STRING) {
                firstj[g[i]-1] = column[i];
                ++ngs;
                if(ngs == ng) break;
              }
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(firstj, column);
          first[j] = firstj;
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          LogicalVector firstj(ng, NA_LOGICAL);
          for(int i = 0; i != gss; ++i) {
            if(column[i] != NA_LOGICAL) {
              if(firstj[g[i]-1] == NA_LOGICAL) {
                firstj[g[i]-1] = column[i];
                ++ngs;
                if(ngs == ng) break;
              }
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(firstj, column);
          first[j] = firstj;
          break;
        }
        default:
          stop("incompatible SEXP encountered;");
          break;
        }
      }
      DUPLICATE_ATTRIB(first, x);
      first.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng); // NumericVector::create(NA_REAL, -ng);
    } else {
      LogicalVector glj(ng, true); //  Much faster method !! (precomputing indices and then going through data)
      IntegerVector firstindex = no_init_vector(ng);
      int ngs = 0;
      for(int i = 0; i != gss; ++i) {
        if(glj[g[i]-1]) {
          glj[g[i]-1] = false;
          firstindex[g[i]-1] = i;
          ++ngs;
          if(ngs == ng) break;
        }
      }
      for(int j = l; j--; ) {
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          first[j] = column[firstindex];
          SHALLOW_DUPLICATE_ATTRIB(first[j], column);
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          first[j] = column[firstindex];
          SHALLOW_DUPLICATE_ATTRIB(first[j], column);
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          first[j] = column[firstindex];
          SHALLOW_DUPLICATE_ATTRIB(first[j], column);
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          first[j] = column[firstindex];
          SHALLOW_DUPLICATE_ATTRIB(first[j], column);
          break;
        }
        default:
          stop("incompatible SEXP encountered;");
          break;
        }
      }
      DUPLICATE_ATTRIB(first, x);
      if(Rf_getAttrib(x, R_RowNamesSymbol) != R_NilValue) {
        const CharacterVector& rn = Rf_getAttrib(x, R_RowNamesSymbol); // const doesn't really make a difference !!
        Rf_setAttrib(first, R_RowNamesSymbol, rn[firstindex]); //  first.attr("row.names") = rn[firstindex]; // Other sloghtly faster, but no big deal !!
      } else {
        first.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng); // NumericVector::create(NA_REAL, -ng);
      }
    }
    return first;
  }
}
