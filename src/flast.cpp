// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// TODO: Implemented smarter copy names ?!

template <int RTYPE>
inline bool isnaNUM(typename Rcpp::traits::storage_type<RTYPE>::type x) {
  return x != x;
}

template <int RTYPE>
inline bool isnaOTH(typename Rcpp::traits::storage_type<RTYPE>::type x) {
  return x == Vector<RTYPE>::get_na();
}

template <int RTYPE>
Vector<RTYPE> flastCppImpl(const Vector<RTYPE>& x, int ng, const IntegerVector& g, bool narm) {
  int l = x.size();
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;
  if (ng == 0) {
    if(narm) {
      int j = l-1;
      auto last = x[j];
      while(isnanT(last) && j!=0) last = x[--j];
      Vector<RTYPE> out(1, last); // faster using create ?
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        Rf_namesgets(out, Rf_ScalarString(names[j]));
      }
      return out;
    } else {
      Vector<RTYPE> out(1, x[l-1]);
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        Rf_namesgets(out, Rf_ScalarString(names[l-1]));
      }
      return out;
    }
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    int ngs = 0;
    Vector<RTYPE> last = no_init_vector(ng);
    DUPLICATE_ATTRIB(last, x);
    if(narm) {
      std::fill(last.begin(), last.end(), Vector<RTYPE>::get_na());
      if(Rf_getAttrib(x, R_NamesSymbol) == R_NilValue) {
        for(int i = l; i--; ) {
          if(!isnanT(x[i])) {
            if(isnanT(last[g[i]-1])) {
              last[g[i]-1] = x[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
      } else {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        if(names.size() != l) stop("x has a names attribute of length != length(x)");
        CharacterVector newnames = no_init_vector(ng);
        for(int i = l; i--; ) {
          if(!isnanT(x[i])) {
            if(isnanT(last[g[i]-1])) {
              last[g[i]-1] = x[i];
              newnames[g[i]-1] = names[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
        Rf_namesgets(last, newnames);
      }
    } else {
      LogicalVector gl(ng, true); // std::vector<bool> glj(ng, true);? -> Nope, not faster (see matrix method)
      if(Rf_getAttrib(x, R_NamesSymbol) == R_NilValue) {
        for(int i = l; i--; ) {
          if(gl[g[i]-1]) {
            gl[g[i]-1] = false;
            last[g[i]-1] = x[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
      } else {
        CharacterVector names = Rf_getAttrib(x, R_NamesSymbol);
        if(names.size() != l) stop("x has a names attribute of length != length(x)");
        CharacterVector newnames = no_init_vector(ng);
        for(int i = l; i--; ) {
          if(gl[g[i]-1]) {
            gl[g[i]-1] = false;
            last[g[i]-1] = x[i];
            newnames[g[i]-1] = names[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
        Rf_namesgets(last, newnames);
      }
    }
    return last;
  }
}

template <>
Vector<CPLXSXP> flastCppImpl(const Vector<CPLXSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<VECSXP> flastCppImpl(const Vector<VECSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<RAWSXP> flastCppImpl(const Vector<RAWSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<EXPRSXP> flastCppImpl(const Vector<EXPRSXP>& x, int ng, const IntegerVector& g, bool narm) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP flastCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true){
  RCPP_RETURN_VECTOR(flastCppImpl, x, ng, g, narm);
}





template <int RTYPE>
SEXP flastmCppImpl(const Matrix<RTYPE>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  int l = x.nrow(), col = x.ncol();
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;

  if(ng == 0) {
    Vector<RTYPE> last = no_init_vector(col);
    if(narm) {
      for(int j = col; j--; ) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        int k = l-1;
        auto lastj = column[k];
        while(isnanT(lastj) && k!=0) lastj = column[--k];
        last[j] = lastj;
      }
    } else {
      last = x(l-1, _);
    }
    if(drop) Rf_setAttrib(last, R_NamesSymbol, colnames(x));
    else {
      Rf_dimgets(last, Dimension(1, col));
      colnames(last) = colnames(x);
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, last);
    }
    return last;
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    Matrix<RTYPE> last = no_init_matrix(ng, col);
    if(narm) {
      std::fill(last.begin(), last.end(), Vector<RTYPE>::get_na());
      for(int j = col; j--; ) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        MatrixColumn<RTYPE> lastj = last( _ , j);
        int ngs = 0;
        for(int i = l; i--; ) {
          if(!isnanT(column[i])) {
            if(isnanT(lastj[g[i]-1])) {
              lastj[g[i]-1] = column[i];
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
      }
      colnames(last) = colnames(x);
    } else {
      List dn = Rf_getAttrib(x, R_DimNamesSymbol);
      if(dn[0] != R_NilValue) {
        CharacterVector rn = dn[0];
        CharacterVector newrn = no_init_vector(ng);
        LogicalVector glj(ng, true); // using std::vector<bool> here is more memory efficient but not faster
        int ngs = 0;
        for(int i = l; i--; ) {
          if(glj[g[i]-1]) {
            glj[g[i]-1] = false;
            last(g[i]-1, _) = x(i, _);
            newrn[g[i]-1] = rn[i];
            ++ngs;
            if(ngs == ng) break;
          }
        }
        Rf_dimnamesgets(last, List::create(newrn, dn[1])); // best way
      } else {
        LogicalVector glj(ng, true); // using std::vector<bool> here is more memory efficient but not faster
        int ngs = 0;
        for(int i = l; i--; ) {
          if(glj[g[i]-1]) {
            glj[g[i]-1] = false;
            last(g[i]-1, _) = x(i, _);
            ++ngs;
            if(ngs == ng) break;
          }
        }
        colnames(last) = colnames(x);
      }
    }
    if(!Rf_isObject(x)) Rf_copyMostAttrib(x, last);
    return last;
  }
}

template <>
SEXP flastmCppImpl(const Matrix<CPLXSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP flastmCppImpl(const Matrix<VECSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP flastmCppImpl(const Matrix<RAWSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP flastmCppImpl(const Matrix<EXPRSXP>& x, int ng, const IntegerVector& g, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP flastmCpp(SEXP x, int ng = 0, IntegerVector g = 0, bool narm = true, bool drop = true){
  RCPP_RETURN_MATRIX(flastmCppImpl, x, ng, g, narm, drop);
}





// [[Rcpp::export]]
SEXP flastlCpp(const List& x, int ng = 0, const IntegerVector& g = 0, bool narm = true) { // , bool drop = true
  int l = x.size();
  List last(l);

  if (ng == 0) {
    if(narm) {
      for(int j = l; j--; ) {
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          int k = column.size()-1;
          while(std::isnan(column[k]) && k!=0) --k;
          NumericVector out(1, column[k]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          last[j] = out;
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          int k = column.size()-1;
          while(column[k] == NA_INTEGER && k!=0) --k;
          IntegerVector out(1, column[k]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          last[j] = out;
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          int k = column.size()-1;
          while(column[k] == NA_STRING && k!=0) --k;
          SEXP out = Rf_ScalarString(column[k]); // CharacterVector out(1, column[k]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          last[j] = out;
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          int k = column.size()-1;
          while(column[k] == NA_LOGICAL && k!=0) --k;
          LogicalVector out(1, column[k]);
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          last[j] = out;
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
          NumericVector out(1, *(column.end()-1));
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          last[j] = out;
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          IntegerVector out(1, *(column.end()-1));
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          last[j] = out;
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          SEXP out = Rf_ScalarString(*(column.end()-1)); // CharacterVector out(1, String(*(column.end()-1)));
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          last[j] = out;
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          LogicalVector out(1, *(column.end()-1));
          SHALLOW_DUPLICATE_ATTRIB(out, column);
          last[j] = out;
          break;
        }
        default:
          stop("incompatible SEXP encountered;");
          break;
        }
      }
    }
    DUPLICATE_ATTRIB(last, x);
    Rf_setAttrib(last, R_RowNamesSymbol, Rf_ScalarInteger(1));
    return last;
  } else { // With groups
    int gss = g.size();
    if(narm) {
      for(int j = l; j--; ) {
        int ngs = 0;
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          NumericVector lastj(ng, NA_REAL);
          for(int i = gss; i--; ) {
            if(!std::isnan(column[i])) {
              if(std::isnan(lastj[g[i]-1])) {
                lastj[g[i]-1] = column[i];
                ++ngs;
                if(ngs == ng) break;
              }
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(lastj, column);
          last[j] = lastj;
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          IntegerVector lastj(ng, NA_INTEGER);
          for(int i = gss; i--; ) {
            if(column[i] != NA_INTEGER) {
              if(lastj[g[i]-1] == NA_INTEGER) {
                lastj[g[i]-1] = column[i];
                ++ngs;
                if(ngs == ng) break;
              }
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(lastj, column);
          last[j] = lastj;
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          CharacterVector lastj(ng, NA_STRING);
          for(int i = gss; i--; ) {
            if(column[i] != NA_STRING) {
              if(lastj[g[i]-1] == NA_STRING) {
                lastj[g[i]-1] = column[i];
                ++ngs;
                if(ngs == ng) break;
              }
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(lastj, column);
          last[j] = lastj;
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          LogicalVector lastj(ng, NA_LOGICAL);
          for(int i = gss; i--; ) {
            if(column[i] != NA_LOGICAL) {
              if(lastj[g[i]-1] == NA_LOGICAL) {
                lastj[g[i]-1] = column[i];
                ++ngs;
                if(ngs == ng) break;
              }
            }
          }
          SHALLOW_DUPLICATE_ATTRIB(lastj, column);
          last[j] = lastj;
          break;
        }
        default:
          stop("incompatible SEXP encountered;");
          break;
        }
      }
      DUPLICATE_ATTRIB(last, x);
      Rf_setAttrib(last, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -ng));
    } else {
      LogicalVector glj(ng, true); //  Much faster method (precomputing indices and then going through data)
      IntegerVector lastindex = no_init_vector(ng);
      int ngs = 0;
      for(int i = gss; i--; ) {
        if(glj[g[i]-1]) {
          glj[g[i]-1] = false;
          lastindex[g[i]-1] = i;
          ++ngs;
          if(ngs == ng) break;
        }
      }
      for(int j = l; j--; ) {
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          last[j] = column[lastindex];
          SHALLOW_DUPLICATE_ATTRIB(last[j], column);
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          last[j] = column[lastindex];
          SHALLOW_DUPLICATE_ATTRIB(last[j], column);
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          last[j] = column[lastindex];
          SHALLOW_DUPLICATE_ATTRIB(last[j], column);
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          last[j] = column[lastindex];
          SHALLOW_DUPLICATE_ATTRIB(last[j], column);
          break;
        }
        default:
          stop("incompatible SEXP encountered;");
          break;
        }
      }
      DUPLICATE_ATTRIB(last, x);
      if(Rf_getAttrib(x, R_RowNamesSymbol) != R_NilValue) {
        const CharacterVector& rn = Rf_getAttrib(x, R_RowNamesSymbol); // const doesn't really make a difference
        Rf_setAttrib(last, R_RowNamesSymbol, rn[lastindex]); //  last.attr("row.names") = rn[lastindex]; // Other slightly faster, but no big deal
      } else {
        Rf_setAttrib(last, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -ng));
      }
    }
    return last;
  }
}
