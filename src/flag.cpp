// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;


// 7th version: Irregular time series and panels supported !
template <int RTYPE>
Vector<RTYPE> flagleadCppImpl(const Vector<RTYPE>& x, const IntegerVector& n, const SEXP& fill,
                              int ng, const IntegerVector& g, const SEXP& t, bool names) {

  // typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
  // storage_t fil;
  Vector<RTYPE> fil(1);
  if(Rf_isNull(fill)) { // fill != fill // Not necessary !!
    fil = Vector<RTYPE>::get_na();
  } else {
    fil = as<Vector<RTYPE> >(fill); //as<storage_t>(fill); -> doesn't work for Character vector fill !!
  }
  auto ff = fil[0];

  int l = x.size(), ns = n.size(), prev = INT_MAX;
  IntegerVector absn = no_init_vector(ns);
  for(int i = 0; i != ns; ++i) {
    if(n[i] == prev) stop("duplicated values in n detected"); // because one might mistakenly pass a factor to the n-slot !!
    prev = n[i];
    if(prev < 0) {
      if(prev == NA_INTEGER) stop("NA in n");
      absn[i] = -prev;
    } else absn[i] = prev;
  }
  if(ns == 1) names = false;
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING;  // NumericVector(abs(n))
  CharacterVector colnam = names ? no_init_vector(ns) : no_init_vector(1);
  Matrix<RTYPE> out = no_init_matrix(l, ns);
  if(ng == 0) { // No groups
    if(Rf_isNull(t)) { // Ordered data
      for(int p = ns; p--; ) {
        int np = n[p];
        if(absn[p] > l) stop("lag-length exceeds length of vector");
        MatrixColumn<RTYPE> outp = out( _ , p);
        if(np>0) {
          if(names) colnam[p] = "L" + nc[p];
          int i = 0;
          while(i != np) outp[i++] = ff;
          for( ; i != l; ++i) outp[i] = x[i - np];
        } else if(np<0) {
          if(names) colnam[p] = "F" + nc[p];
          int i = l, st = l+np;
          while(i != st) outp[--i] = ff;
          for( ; i--; ) outp[i] = x[i - np];
        } else {
          if(names) colnam[p] = "--";
          outp = x;
        }
      }
    } else { // Unordered data: Timevar provided
      IntegerVector ord = t;
      if(l != ord.size()) stop("length(x) must match length(t)");
      int min = INT_MAX, max = INT_MIN, osize, temp;
      for(int i = 0; i != l; ++i) {
        if(ord[i] < min) min = ord[i];
        if(ord[i] > max) max = ord[i];
      }
      if(min == NA_INTEGER) stop("Timevar contains missing values");
      osize = max-min+1;
      if(osize > 3 * l) warning("Your time series is very irregular. Need to create an internal ordering vector of length %s to represent it.", osize);
      IntegerVector omap(osize), ord2 = no_init_vector(l);
      for(int i = 0; i != l; ++i) {
        temp = ord[i] - min; // Best ? Or direct assign to ord2[i] ? Also check for panel version..
        if(omap[temp]) stop("Repeated values in timevar");
        omap[temp] = i+1;
        ord2[i] = temp;
      }
      // return as<Vector<RTYPE> >(omap);
      for(int p = ns; p--; ) {
        int np = n[p];
        if(absn[p] > l) stop("lag-length exceeds length of vector");
        MatrixColumn<RTYPE> outp = out( _ , p);
        if(np>0) {
          if(names) colnam[p] = "L" + nc[p];
          for(int i = 0; i != l; ++i) { // Smarter solution using while ???
            if(ord2[i] >= np && (temp = omap[ord2[i] - np])) {
              outp[i] = x[temp-1];
            } else {
              outp[i] = ff;
            }
          }
        } else if(np<0) {
          if(names) colnam[p] = "F" + nc[p];
          for(int i = 0, osnp = osize+np; i != l; ++i) { // Smarter solution using while ???
            if(ord2[i] < osnp && (temp = omap[ord2[i] - np])) {
              outp[i] = x[temp-1];
            } else {
              outp[i] = ff;
            }
          }
        } else {
          if(names) colnam[p] = "--";
          outp = x;
        }
      }
    }
  } else { // With groups
    if(l != g.size()) stop("length(x) must match length(g)");
    int ags = l/ng, ngp = ng+1;
    if(Rf_isNull(t)) { // Ordered data
      // int seen[ngp], memsize = sizeof(int)*ngp;
      for(int p = ns; p--; ) {
        int np = n[p];
        if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
        MatrixColumn<RTYPE> outp = out( _ , p);
        if(np>0) {
          if(names) colnam[p] = "L" + nc[p];
          std::vector<int> seen(ngp); // memset(seen, 0, memsize);
          for(int i = 0; i != l; ++i) {
            if(seen[g[i]] == np) {
              outp[i] = x[i-np];
            } else {
              outp[i] = ff;
              ++seen[g[i]];
            }
          }
        } else if(np<0) {
          std::vector<int> seen(ngp); // memset(seen, 0, memsize);
          if(names) colnam[p] = "F" + nc[p];
          for(int i = l; i--; ) { // good??
            if(seen[g[i]] == np) {
              outp[i] = x[i-np];
            } else {
              outp[i] = ff;
              --seen[g[i]];
            }
          }
        } else {
          if(names) colnam[p] = "--";
          outp = x;
        }
      }
    } else { // Unordered data: Timevar provided
      IntegerVector ord = t;
      int temp;
      if(l != ord.size()) stop("length(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX), max(ngp, INT_MIN), cgs = no_init_vector(ngp);
      for(int i = 0; i != l; ++i) {
        temp = g[i];
        if(ord[i] < min[temp]) min[temp] = ord[i];
        if(ord[i] > max[temp]) max[temp] = ord[i];
      }
      temp = 0;
      for(int i = 1; i != ngp; ++i) {
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values");
        if(min[i] == INT_MAX) continue; // Needed in case of unused factor levels (group vector too large)
        cgs[i] = temp; // This needs to b here (for unused factor levels case...)
        max[i] -= min[i] - 1; // need max[i] which stores the complete group sizes only if p<0 e.g. if computing leads..
        temp += max[i];
      }
      // omap provides the ordering to order the vector (needed to find previous / next values)
      if(temp > 3 * l) warning("Your panel is very irregular. Need to create an internal ordering vector of length %s to represent it.", temp);
      IntegerVector omap(temp), ord2 = no_init_vector(l);
      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]];
        temp = cgs[g[i]] + ord2[i];
        if(omap[temp]) stop("Repeated values of timevar within one or more groups");
        omap[temp] = i+1; // needed to add 1 to distinguish between 0 and gap
      }
      for(int p = ns; p--; ) {
        int np = n[p];
        if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
        MatrixColumn<RTYPE> outp = out( _ , p);
        if(np>0) {
          if(names) colnam[p] = "L" + nc[p];
          for(int i = 0; i != l; ++i) {
            if(ord2[i] >= np && (temp = omap[cgs[g[i]]+ord2[i]-np])) {
              outp[i] = x[temp-1];
            } else {
              outp[i] = ff;
            }
          }
        } else if(np<0) {
          if(names) colnam[p] = "F" + nc[p];
          for(int i = 0; i != l; ++i) {
            if(ord2[i] < max[g[i]]+np && (temp = omap[cgs[g[i]]+ord2[i]-np])) {
              outp[i] = x[temp-1];
            } else {
              outp[i] = ff;
            }
          }
        } else {
          if(names) colnam[p] = "--";
          outp = x;
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  if(ns != 1) {
    Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
    Rf_dimgets(out, Dimension(l, ns));
    if(Rf_isObject(x)) {
      CharacterVector classes = Rf_getAttrib(out, R_ClassSymbol);
      classes.push_back("matrix");
      Rf_classgets(out, classes);
    } else {
      Rf_classgets(out, Rf_mkString("matrix"));
    }
    if(names) Rf_dimnamesgets(out, List::create(Rf_getAttrib(x, R_NamesSymbol), colnam));
    // out.attr("class") = CharacterVector::create(x.attr("class"),"matrix");
  }
  return out;
}

template <>
Vector<CPLXSXP> flagleadCppImpl(const Vector<CPLXSXP>& x, const IntegerVector& n, const SEXP& fill,
                               int ng, const IntegerVector& g, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Vector<VECSXP> flagleadCppImpl(const Vector<VECSXP>& x, const IntegerVector& n, const SEXP& fill,
                               int ng, const IntegerVector& g, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Vector<RAWSXP> flagleadCppImpl(const Vector<RAWSXP>& x, const IntegerVector& n, const SEXP& fill,
                               int ng, const IntegerVector& g, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Vector<EXPRSXP> flagleadCppImpl(const Vector<EXPRSXP>& x, const IntegerVector& n, const SEXP& fill,
                                int ng, const IntegerVector& g, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP flagleadCpp(SEXP x, IntegerVector n = 1, SEXP fill = R_NilValue,
                 int ng = 0, IntegerVector g = 0, SEXP t = R_NilValue, bool names = true){
  RCPP_RETURN_VECTOR(flagleadCppImpl, x, n, fill, ng, g, t, names);
}



inline SEXP coln_check(SEXP x) {
  if(Rf_isNull(x)) return NA_STRING;
  else return x; // Rf_coerceVector(x, STRSXP);
}

template <int RTYPE>
Matrix<RTYPE> flagleadmCppImpl(const Matrix<RTYPE>& x, const IntegerVector& n, const SEXP& fill,
                               int ng, const IntegerVector& g, const SEXP& t, bool names) {

  Vector<RTYPE> fil(1);
  if(Rf_isNull(fill)) { //  || fill != fill not necessary !!
    fil = Vector<RTYPE>::get_na();
  } else {
    fil = as<Vector<RTYPE> >(fill);
  }
  auto ff = fil[0];

  int l = x.nrow(), col = x.ncol(), ns = n.size(), pos = INT_MAX;
  IntegerVector absn = no_init_vector(ns);
  for(int i = 0; i != ns; ++i) {
    if(n[i] == pos) stop("duplicated values in n detected"); // because one might mistakenly pass a factor to the n-slot !!
    pos = n[i];
    if(pos < 0) {
      if(pos == NA_INTEGER) stop("NA in n");
      absn[i] = -pos;
    } else absn[i] = pos;
  }
  pos = 0;
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING; // NumericVector(abs(n))
  CharacterVector colnam = names ? no_init_vector(col*ns) : no_init_vector(1); // what if no names ??
  CharacterVector coln = names ? coln_check(colnames(x)) : NA_STRING;
  if(names && coln[0] == NA_STRING) names = false;

  Matrix<RTYPE> out = no_init_matrix(l, col*ns);

  if(ng == 0) { // No groups
    if(Rf_isNull(t)) { // Ordered data
      for(int j = 0; j != col; ++j) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p] > l) stop("lag-length exceeds length of vector");
          MatrixColumn<RTYPE> outj = out( _ , pos);
          if(np>0) {
            if(names) colnam[pos] = "L" + nc[p] + "." + coln[j];
            int i = 0;
            while(i != np) outj[i++] = ff;
            for( ; i != l; ++i) outj[i] = column[i - np];
          } else if(np<0) {
            if(names) colnam[pos] = "F" + nc[p] + "." + coln[j];
            int i = l, st = l+np;
            while(i != st) outj[--i] = ff;
            for( ; i--; ) outj[i] = column[i - np];
          } else {
            if(names) colnam[pos] = coln[j];
            outj = column;
          }
          ++pos;
        }
      }
    } else { // Unordered data: Timevar Provided
      IntegerVector ord = t;
      if(l != ord.size()) stop("length(x) must match length(t)");
      int min = INT_MAX, max = INT_MIN, osize, temp;
      for(int i = 0; i != l; ++i) {
        if(ord[i] < min) min = ord[i];
        if(ord[i] > max) max = ord[i];
      }
      if(min == NA_INTEGER) stop("Timevar contains missing values");
      osize = max-min+1;
      if(osize > 3 * l) warning("Your time series is very irregular. Need to create an internal ordering vector of length %s to represent it.", osize);
      IntegerVector omap(osize), ord2 = no_init_vector(l);
      for(int i = 0; i != l; ++i) {
        temp = ord[i] - min; // Best ? Or direct assign to ord2[i] ? Also check for panel version..
        if(omap[temp]) stop("Repeated values in timevar");
        omap[temp] = i+1;
        ord2[i] = temp;
      }
      for(int j = 0; j != col; ++j) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p] > l) stop("lag-length exceeds length of vector");
          MatrixColumn<RTYPE> outj = out( _ , pos);
          if(np>0) {
            if(names) colnam[pos] = "L" + nc[p] + "." + coln[j];
            for(int i = 0; i != l; ++i) { // Smarter solution using while ???
              if(ord2[i] >= np && (temp = omap[ord2[i] - np])) {
                outj[i] = column[temp-1];
              } else {
                outj[i] = ff;
              }
            }
          } else if(np<0) {
            if(names) colnam[pos] = "F" + nc[p] + "." + coln[j];
            for(int i = 0, osnp = osize+np; i != l; ++i) { // Smarter solution using while ???
              if(ord2[i] < osnp && (temp = omap[ord2[i] - np])) {
                outj[i] = column[temp-1];
              } else {
                outj[i] = ff;
              }
            }
          } else {
            if(names) colnam[pos] = coln[j];
            outj = column;
          }
          ++pos;
        }
      }
    }
  } else { // With groups
    if(l != g.size()) stop("length(x) must match length(g)");
    int ags = l/ng, ngp = ng+1;
    if(Rf_isNull(t)) { // Ordered data
      // int seen[ngp], memsize = sizeof(int)*ngp;
      for(int j = 0; j != col; ++j) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
          MatrixColumn<RTYPE> outj = out( _ , pos);
          if(np>0) {
            if(names) colnam[pos] = "L" + nc[p] + "." + coln[j];
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            for(int i = 0; i != l; ++i) {
              if(seen[g[i]] == np) {
                outj[i] = column[i-np];
              } else {
                outj[i] = ff;
                ++seen[g[i]];
              }
            }
          } else if(np<0) {
            if(names) colnam[pos] = "F" + nc[p] + "." + coln[j];
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            for(int i = l; i--; ) { // good??
              if(seen[g[i]] == np) {
                outj[i] = column[i-np];
              } else {
                outj[i] = ff;
                --seen[g[i]];
              }
            }
          } else {
            if(names) colnam[pos] = coln[j];
            outj = column;
          }
          ++pos;
        }
      }
    } else { // Unordered data: Timevar provided
      IntegerVector ord = t;
      int temp;
      if(l != ord.size()) stop("length(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX), max(ngp, INT_MIN), cgs = no_init_vector(ngp);
      for(int i = 0; i != l; ++i) {
        temp = g[i];
        if(ord[i] < min[temp]) min[temp] = ord[i];
        if(ord[i] > max[temp]) max[temp] = ord[i];
      }
      temp = 0;
      for(int i = 1; i != ngp; ++i) {
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values");
        if(min[i] == INT_MAX) continue; // Needed in case of unused factor levels (group vector too large)
        cgs[i] = temp; // This needs to b here (for unused factor levels case...)
        max[i] -= min[i] - 1; // need max[i] which stores the complete group sizes only if p<0 e.g. if computing leads..
        temp += max[i];
      }
      // omap provides the ordering to order the vector (needed to find previous / next values)
      if(temp > 3 * l) warning("Your panel is very irregular. Need to create an internal ordering vector of length %s to represent it.", temp);
      IntegerVector omap(temp), ord2 = no_init_vector(l), index = no_init_vector(l);
      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]];
        index[i] = cgs[g[i]] + ord2[i];
        if(omap[index[i]]) stop("Repeated values of timevar within one or more groups");
        omap[index[i]] = i+1; // needed to add 1 to distinguish between 0 and gap
      }
      for(int j = 0; j != col; ++j) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
          MatrixColumn<RTYPE> outj = out( _ , pos);
          if(np>0) {
            if(names) colnam[pos] = "L" + nc[p] + "." + coln[j];
            for(int i = 0; i != l; ++i) {
              if(ord2[i] >= np && (temp = omap[index[i]-np])) {
                outj[i] = column[temp-1];
              } else {
                outj[i] = ff;
              }
            }
          } else if(np<0) {
            if(names) colnam[pos] = "F" + nc[p] + "." + coln[j];
            for(int i = 0; i != l; ++i) { // best loop ??
              if(ord2[i] < max[g[i]]+np && (temp = omap[index[i]-np])) {
                outj[i] = column[temp-1];
              } else {
                outj[i] = ff;
              }
            }
          } else {
            if(names) colnam[pos] = coln[j];
            outj = column;
          }
          ++pos;
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  if(ns != 1) Rf_dimgets(out, Dimension(l, col*ns));
  if(names) {
    Rf_dimnamesgets(out, List::create(rownames(x), colnam)); // colnames(out) = colnam deletes row names !
  } else if(ns != 1) {
    Rf_setAttrib(out, R_DimNamesSymbol, R_NilValue);
  }
  return out;
}



template <>
Matrix<CPLXSXP> flagleadmCppImpl(const Matrix<CPLXSXP>& x, const IntegerVector& n, const SEXP& fill,
                                int ng, const IntegerVector& g, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Matrix<VECSXP> flagleadmCppImpl(const Matrix<VECSXP>& x, const IntegerVector& n, const SEXP& fill,
                                int ng, const IntegerVector& g, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Matrix<RAWSXP> flagleadmCppImpl(const Matrix<RAWSXP>& x, const IntegerVector& n, const SEXP& fill,
                                int ng, const IntegerVector& g, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Matrix<EXPRSXP> flagleadmCppImpl(const Matrix<EXPRSXP>& x, const IntegerVector& n, const SEXP& fill,
                                 int ng, const IntegerVector& g, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP flagleadmCpp(SEXP x, IntegerVector n = 1, SEXP fill = R_NilValue,
                  int ng = 0, IntegerVector g = 0, SEXP t = R_NilValue, bool names = true){
  RCPP_RETURN_MATRIX(flagleadmCppImpl, x, n, fill, ng, g, t, names);
}



// [[Rcpp::export]]
List flagleadlCpp(const List& x, const IntegerVector& n = 1, const SEXP& fill = R_NilValue,
                  int ng = 0, const IntegerVector& g = 0, const SEXP& t = R_NilValue, bool names = true) {

  bool lfill = Rf_isNull(fill);
  if(!lfill && TYPEOF(fill) == LGLSXP) lfill = Rf_asLogical(fill) == NA_LOGICAL;
  int l = x.size(), ns = n.size(), pos = INT_MAX;
  List out(l * ns);
  IntegerVector absn = no_init_vector(ns);
  for(int i = 0; i != ns; ++i) {
    if(n[i] == pos) stop("duplicated values in n detected"); // because one might mistakenly pass a factor to the n-slot !!
    pos = n[i];
    if(pos < 0) {
      if(pos == NA_INTEGER) stop("NA in n");
      absn[i] = -pos;
    } else absn[i] = pos;
  }
  pos = 0;
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING; // NumericVector(abs(n))
  CharacterVector nam = names ? no_init_vector(l*ns) : no_init_vector(1); // what if no names ??
  CharacterVector na = names ? coln_check(Rf_getAttrib(x, R_NamesSymbol)) : NA_STRING;
  if(names && na[0] == NA_STRING) names = false;

  if(ng == 0) { // No groups
    if(Rf_isNull(t)) { // Ordered data
      for(int j = 0; j != l; ++j) {
        int txj = TYPEOF(x[j]);
        switch(txj) {
        case REALSXP: {
          NumericVector column = x[j];
          int row = column.size();
          double ff = lfill ? NA_REAL : Rf_asReal(fill); // as<double>()
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > row) stop("lag-length exceeds length of vector");
            if(np>0) {
              NumericVector outjp = no_init_vector(row);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[i++] = ff;
              for( ; i != row; ++i) outjp[i] = column[i - np];
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              NumericVector outjp = no_init_vector(row);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int i = row, st = row+np;
              while(i != st) outjp[--i] = ff;
              for( ; i--; ) outjp[i] = column[i - np];
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case LGLSXP:
        case INTSXP: {
          IntegerVector column = x[j];
          int row = column.size();
          int ff = lfill ? NA_INTEGER : Rf_asInteger(fill); // as<int>()
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > row) stop("lag-length exceeds length of vector");
            if(np>0) {
              IntegerVector outjp = no_init_vector(row);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[i++] = ff;
              for( ; i != row; ++i) outjp[i] = column[i - np];
              DUPLICATE_ATTRIB(outjp, column);
              if(txj == LGLSXP) SET_TYPEOF(outjp, LGLSXP);
              out[pos] = outjp;
            } else if(np<0) {
              IntegerVector outjp = no_init_vector(row);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int i = row, st = row+np;
              while(i != st) outjp[--i] = ff;
              for( ; i--; ) outjp[i] = column[i - np];
              DUPLICATE_ATTRIB(outjp, column);
              if(txj == LGLSXP) SET_TYPEOF(outjp, LGLSXP);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = x[j];
            }
            ++pos;
          }
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          int row = column.size();
          // String ff = lfill ? NA_STRING : as<String>(fill); // String
          SEXP ff = lfill ? NA_STRING : Rf_asChar(fill);
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > row) stop("lag-length exceeds length of vector");
            if(np>0) {
              CharacterVector outjp = no_init_vector(row);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[i++] = ff;
              for( ; i != row; ++i) outjp[i] = column[i - np];
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              CharacterVector outjp = no_init_vector(row);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int i = row, st = row+np;
              while(i != st) outjp[--i] = ff;
              for( ; i--; ) outjp[i] = column[i - np];
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        default: stop("Not supported SEXP type!");
        }
      }
    } else { // Unordered data: Timevar Provided
      IntegerVector ord = t;
      int min = INT_MAX, max = INT_MIN, osize, temp, os = ord.size();
      if(Rf_length(x[0]) != os) stop("nrow(x) must match length(t)");
      for(int i = 0; i != os; ++i) {
        if(ord[i] < min) min = ord[i];
        if(ord[i] > max) max = ord[i];
      }
      if(min == NA_INTEGER) stop("Timevar contains missing values");
      osize = max-min+1;
      if(osize > 3 * os) warning("Your time series is very irregular. Need to create an internal ordering vector of length %s to represent it.", osize);
      IntegerVector omap(osize), ord2 = no_init_vector(os);
      for(int i = 0; i != os; ++i) {
        temp = ord[i] - min; // Best ? Or direct assign to ord2[i] ? Also check for panel version..
        if(omap[temp]) stop("Repeated values in timevar");
        omap[temp] = i+1;
        ord2[i] = temp;
      }
      for(int j = 0; j != l; ++j) {
        int txj = TYPEOF(x[j]);
        switch(txj) {
        case REALSXP: {
          NumericVector column = x[j];
          if(os != column.size()) stop("nrow(x) must match length(t)");
          double ff = lfill ? NA_REAL : Rf_asReal(fill); // as<double>(
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > os) stop("lag-length exceeds length of vector");
            if(np>0) {
              NumericVector outjp = no_init_vector(os);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              for(int i = 0; i != os; ++i) { // Smarter solution using while ???
                if(ord2[i] >= np && (temp = omap[ord2[i] - np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              NumericVector outjp = no_init_vector(os);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              for(int i = 0, osnp = osize+np; i != os; ++i) { // Smarter solution using while ???
                if(ord2[i] < osnp && (temp = omap[ord2[i] - np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case LGLSXP:
        case INTSXP: {
          IntegerVector column = x[j];
          if(os != column.size()) stop("length(x) must match length(t)");
          int ff = lfill ? NA_INTEGER : Rf_asInteger(fill); // as<int>(
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > os) stop("lag-length exceeds length of vector");
            if(np>0) {
              IntegerVector outjp = no_init_vector(os);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              for(int i = 0; i != os; ++i) { // Smarter solution using while ???
                if(ord2[i] >= np && (temp = omap[ord2[i] - np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              if(txj == LGLSXP) SET_TYPEOF(outjp, LGLSXP);
              out[pos] = outjp;
            } else if(np<0) {
              IntegerVector outjp = no_init_vector(os);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              for(int i = 0, osnp = osize+np; i != os; ++i) { // Smarter solution using while ???
                if(ord2[i] < osnp && (temp = omap[ord2[i] - np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              if(txj == LGLSXP) SET_TYPEOF(outjp, LGLSXP);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = x[j];
            }
            ++pos;
          }
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          if(os != column.size()) stop("length(x) must match length(t)");
          // String ff = lfill ? NA_STRING : as<String>(fill); // String ??
          SEXP ff = lfill ? NA_STRING : Rf_asChar(fill);
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > os) stop("lag-length exceeds length of vector");
            if(np>0) {
              CharacterVector outjp = no_init_vector(os);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              for(int i = 0; i != os; ++i) { // Smarter solution using while ???
                if(ord2[i] >= np && (temp = omap[ord2[i] - np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              CharacterVector outjp = no_init_vector(os);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              for(int i = 0, osnp = osize+np; i != os; ++i) { // Smarter solution using while ???
                if(ord2[i] < osnp && (temp = omap[ord2[i] - np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        default: stop("Not supported SEXP type!");
        }
      }
    }
  } else { // With groups
    int gss = g.size(), ags = gss/ng, ngp = ng+1, temp = 0;
    if(Rf_isNull(t)) { // Ordered data
      std::vector<int> seen(ngp); // int seen[ngp], memsize = sizeof(int)*ngp;
      for(int j = 0; j != l; ++j) {
        int txj = TYPEOF(x[j]);
        switch(txj) {
        case REALSXP: {
          NumericVector column = x[j];
          double ff = lfill ? NA_REAL : Rf_asReal(fill); // as<double>()
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
            if(np>0) {
              NumericVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              seen.assign(ngp, 0); // std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = 0; i != gss; ++i) {
                if(seen[g[i]] == np) {
                  outjp[i] = column[i-np];
                } else {
                  outjp[i] = ff;
                  ++seen[g[i]];
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              NumericVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              seen.assign(ngp, 0); //std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = gss; i--; ) { // good??
                if(seen[g[i]] == np) {
                  outjp[i] = column[i-np];
                } else {
                  outjp[i] = ff;
                  --seen[g[i]];
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case LGLSXP:
        case INTSXP: {
          IntegerVector column = x[j];
          int ff = lfill ? NA_INTEGER : Rf_asInteger(fill); // as<int>()
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
            if(np>0) {
              IntegerVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              seen.assign(ngp, 0); //std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = 0; i != gss; ++i) {
                if(seen[g[i]] == np) {
                  outjp[i] = column[i-np];
                } else {
                  outjp[i] = ff;
                  ++seen[g[i]];
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              if(txj == LGLSXP) SET_TYPEOF(outjp, LGLSXP);
              out[pos] = outjp;
            } else if(np<0) {
              IntegerVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              seen.assign(ngp, 0); //std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = gss; i--; ) { // good??
                if(seen[g[i]] == np) {
                  outjp[i] = column[i-np];
                } else {
                  outjp[i] = ff;
                  --seen[g[i]];
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              if(txj == LGLSXP) SET_TYPEOF(outjp, LGLSXP);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = x[j];
            }
            ++pos;
          }
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          // String ff = lfill ? NA_STRING : as<String>(fill); // String ??
          SEXP ff = lfill ? NA_STRING : Rf_asChar(fill);
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
            if(np>0) {
              CharacterVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              seen.assign(ngp, 0); //std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = 0; i != gss; ++i) {
                if(seen[g[i]] == np) {
                  outjp[i] = column[i-np];
                } else {
                  outjp[i] = ff;
                  ++seen[g[i]];
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              CharacterVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              seen.assign(ngp, 0); //std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = gss; i--; ) { // good??
                if(seen[g[i]] == np) {
                  outjp[i] = column[i-np];
                } else {
                  outjp[i] = ff;
                  --seen[g[i]];
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        default: stop("Not supported SEXP type!");
        }
      }
    } else { // Unordered data: Timevar provided
      IntegerVector ord = t;
      if(gss != ord.size()) stop("length(g) must match length(t)");
      IntegerVector min(ngp, INT_MAX), max(ngp, INT_MIN), cgs = no_init_vector(ngp);
      // return List::create(min, max); // Note: INT_MIN is the same as NA_INTEGER
      for(int i = 0; i != gss; ++i) {
        temp = g[i];
        if(ord[i] < min[temp]) min[temp] = ord[i];
        if(ord[i] > max[temp]) max[temp] = ord[i];
      }
      temp = 0;
      for(int i = 1; i != ngp; ++i) {
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values");
        if(min[i] == INT_MAX) continue; // Needed in case of unused factor levels (group vector too large)
        cgs[i] = temp; // This needs to b here (for unused factor levels case...)
        max[i] -= min[i] - 1; // need max[i] which stores the complete group sizes only if p<0 e.g. if computing leads..
        temp += max[i]; // + max[i] - min[i] + 1;
      }
      // if(min[ng] == NA_INTEGER) stop("Timevar contains missing values");
      // if(min[ng] != INT_MAX) {
      //   max[ng] -= min[ng] - 1;
      //   temp += max[ng];
      // }
      // return List::create(cgs, min, max);
      // index stores the position of the current observation in the ordered vector
      // omap provides the ordering to order the vector (needed to find previous / next values)
      if(temp > 3 * gss) warning("Your panel is very irregular. Need to create an internal ordering vector of length %s to represent it.", temp);
      IntegerVector omap(temp), ord2 = no_init_vector(gss), index = no_init_vector(gss);
      for(int i = 0; i != gss; ++i) {
        ord2[i] = ord[i] - min[g[i]]; // Need ord2 can get rid of any part ?? ??
        // if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        index[i] = cgs[g[i]] + ord2[i];
        if(omap[index[i]]) stop("Repeated values of timevar within one or more groups");
        omap[index[i]] = i+1; // needed to add 1 to distinguish between 0 and gap
      }
      // return List::create(cgs, min, max, ord2, index, omap);
      for(int j = 0; j != l; ++j) {
        int txj = TYPEOF(x[j]);
        switch(txj) {
        case REALSXP: {
          NumericVector column = x[j];
          double ff = lfill ? NA_REAL : Rf_asReal(fill); // as<double>()
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
            if(np>0) {
              NumericVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) {
                if(ord2[i] >= np && (temp = omap[index[i]-np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              NumericVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) { // best loop ??
                if(ord2[i] < max[g[i]]+np && (temp = omap[index[i]-np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case LGLSXP:
        case INTSXP: {
          IntegerVector column = x[j];
          int ff = lfill ? NA_INTEGER : Rf_asInteger(fill); // as<int>
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
            if(np>0) {
              IntegerVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) {
                if(ord2[i] >= np && (temp = omap[index[i]-np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              if(txj == LGLSXP) SET_TYPEOF(outjp, LGLSXP);
              out[pos] = outjp;
            } else if(np<0) {
              IntegerVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) { // best loop ??
                if(ord2[i] < max[g[i]]+np && (temp = omap[index[i]-np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              if(txj == LGLSXP) SET_TYPEOF(outjp, LGLSXP);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = x[j];
            }
            ++pos;
          }
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          // String ff = lfill ? NA_STRING : as<String>(fill);
          SEXP ff = lfill ? NA_STRING : Rf_asChar(fill);
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
            if(np>0) {
              CharacterVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) {
                if(ord2[i] >= np && (temp = omap[index[i]-np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              CharacterVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) { // best loop ??
                if(ord2[i] < max[g[i]]+np && (temp = omap[index[i]-np])) {
                  outjp[i] = column[temp-1];
                } else {
                  outjp[i] = ff;
                }
              }
              DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        default: stop("Not supported SEXP type!");
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  if(names) { // best way to code this ??
    Rf_namesgets(out, nam);
  } else {
    if(ns != 1) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
  }
  return out;
}

