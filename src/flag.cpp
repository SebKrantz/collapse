// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;


// 6th version: Type Dispatch and names argument !!
template <int RTYPE>
Vector<RTYPE> flagleadCppImpl(const Vector<RTYPE>& x, const IntegerVector& n, const SEXP& fill,
                              int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {

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
      LogicalVector ocheck(l, true);
      IntegerVector omap = no_init_vector(l);
      for(int i = 0; i != l; ++i) {
        if(ord[i] > l || ord[i] < 1) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i; // Note: omap is the same as order(ord) !!
        } else {
          stop("Repeated values in timevar");
        }
      }
      for(int p = ns; p--; ) {
        int np = n[p];
        if(absn[p] > l) stop("lag-length exceeds length of vector");
        MatrixColumn<RTYPE> outp = out( _ , p);
        if(np>0) {
          if(names) colnam[p] = "L" + nc[p];
          int i = 0;
          while(i != np) outp[omap[i++]] = ff;
          for( ; i != l; ++i) outp[omap[i]] = x[omap[i - np]];
        } else if(np<0) {
          if(names) colnam[p] = "F" + nc[p];
          int st = l+np, i = l;
          while(i != st) outp[omap[--i]] = ff;
          for( ; i--; ) outp[omap[i]] = x[omap[i - np]];
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
        if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
      if(l != ord.size()) stop("length(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX); // INFINITY gives bug !!
      IntegerVector gsv = Rf_isNull(gs) ? IntegerVector(ng) : as<IntegerVector>(gs); // no_init_vector(ng); // No real improvements here by using C++ arrays !!
      IntegerVector ord2 = no_init_vector(l); // use array ??
      if(Rf_isNull(gs)) {
        // std::fill(gsv.begin(), gsv.end(), 0);
        // gsv = IntegerVector(ng);
        for(int i = 0; i != l; ++i) {
          ++gsv[g[i]-1];
          if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
        }
      } else {
        // gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
      }
      IntegerVector omap(l), cgs = no_init_vector(ngp);
      // int cgs[ngp];
      // cgs[1] = 0;
      // for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; // or get "starts from forderv"
      cgs[1] = 0;
      for(int i = 1; i != ng; ++i) {
        cgs[i+1] = cgs[i] + gsv[i-1]; // or get "starts from forderv"
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?
      }
      if(min[ng] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?
      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]]; // still room for speed improvement ?? -> could get rid of the first if condition, if there is a gap, there will be a repeated value error later on !!
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        if(omap[cgs[g[i]]+ord2[i]] == 0) omap[cgs[g[i]]+ord2[i]] = i;
        else stop("Repeated values of timevar within one or more groups");
      }
      for(int p = ns; p--; ) {
        int np = n[p];
        if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
        MatrixColumn<RTYPE> outp = out( _ , p);
        if(np>0) {
          if(names) colnam[p] = "L" + nc[p];
          for(int i = 0; i != l; ++i) {
            if(ord2[i] >= np) {
              outp[i] = x[omap[cgs[g[i]]+ord2[i]-np]];
            } else {
              outp[i] = ff;
            }
          }
        } else if(np<0) {
          if(names) colnam[p] = "F" + nc[p];
          for(int i = 0; i != l; ++i) {
            if(ord2[i] < gsv[g[i]-1]+np) {
              outp[i] = x[omap[cgs[g[i]]+ord2[i]-np]];
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
    Rf_setAttrib(out, R_NamesSymbol, R_NilValue); // if(x.hasAttribute("names")) out.attr("names") = R_NilValue; // fastest ?? Rf_setAttrib(x, R_NamesSymbol, R_NilValue);
    out.attr("dim") = Dimension(l, ns); // Rf_dimgets(Dimension(l, ns));
    if(Rf_isObject(x)) { // out.attr("class") = CharacterVector::create(out.attr("class"),"matrix"); // Rf_dimgets(Dimension(l, ns));
      CharacterVector classes = out.attr("class");
      classes.push_back("matrix");
      out.attr("class") = classes;
    } else {
      out.attr("class") = "matrix"; // Rf_classgets(Rf_installChar("matrix"));
    }
    if(names) out.attr("dimnames") = List::create(x.attr("names"), colnam); // Rf_dimnamesgets(List::create(x.attr("names"), colnam));
    // out.attr("class") = CharacterVector::create(x.attr("class"),"matrix");
  }
  return out;
}

// template <>
// Vector<CPLXSXP> flagleadCppImpl(Vector<CPLXSXP> x, IntegerVector n, SEXP fill,
//                               int ng, IntegerVector g, SEXP gs, SEXP t) {
//   stop("Not supported SEXP type!");
// }

template <>
Vector<VECSXP> flagleadCppImpl(const Vector<VECSXP>& x, const IntegerVector& n, const SEXP& fill,
                               int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Vector<RAWSXP> flagleadCppImpl(const Vector<RAWSXP>& x, const IntegerVector& n, const SEXP& fill,
                               int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Vector<EXPRSXP> flagleadCppImpl(const Vector<EXPRSXP>& x, const IntegerVector& n, const SEXP& fill,
                                int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP flagleadCpp(SEXP x, IntegerVector n = 1, SEXP fill = R_NilValue,
                 int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue, bool names = true){
  RCPP_RETURN_VECTOR(flagleadCppImpl, x, n, fill, ng, g, gs, t, names);
}



inline SEXP coln_check(SEXP x) {
  if(Rf_isNull(x)) return NA_STRING;
  else return x; // Rf_coerceVector(x, STRSXP);
}

template <int RTYPE>
Matrix<RTYPE> flagleadmCppImpl(const Matrix<RTYPE>& x, const IntegerVector& n, const SEXP& fill,
                               int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {

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
      LogicalVector ocheck(l, true);
      IntegerVector omap = no_init_vector(l);
      for(int i = 0; i != l; ++i) {
        if(ord[i] > l || ord[i] < 1) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i; // Note: omap is the same as order(ord) !!
        } else {
          stop("Repeated values in timevar");
        }
      }
      for(int j = 0; j != col; ++j) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p] > l) stop("lag-length exceeds length of vector");
          MatrixColumn<RTYPE> outj = out( _ , pos);
          if(np>0) {
            if(names) colnam[pos] = "L" + nc[p] + "." + coln[j];
            int i = 0;
            while(i != np) outj[omap[i++]] = ff;
            for( ; i != l; ++i) outj[omap[i]] = column[omap[i - np]];
          } else if(np<0) {
            if(names) colnam[pos] = "F" + nc[p] + "." + coln[j];
            int st = l+np, i = l;
            while(i != st) outj[omap[--i]] = ff;
            for( ; i--; ) outj[omap[i]] = column[omap[i - np]];
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
          if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
      if(l != ord.size()) stop("length(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX);
      IntegerVector gsv = Rf_isNull(gs) ? IntegerVector(ng) : as<IntegerVector>(gs); // no_init_vector(ng); // NULL; gives compiler warning
      IntegerVector ord2 = no_init_vector(l); // See flag.cpp for any improvements on this code !!
      if(Rf_isNull(gs)) {
        // gsv = IntegerVector(ng);
        // std::fill(gsv.begin(), gsv.end(), 0);
        for(int i = 0; i != l; ++i) {
          ++gsv[g[i]-1];
          if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
        }
      } else {
        // gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
      }
      IntegerVector omap(l), cgs = no_init_vector(ngp), index = no_init_vector(l);
      // int cgs[ngp], index[l]; // See above for any improvements on this code !!
      cgs[1] = 0;
      for(int i = 1; i != ng; ++i) {
        cgs[i+1] = cgs[i] + gsv[i-1]; // or get "starts from forderv"
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?
      }
      if(min[ng] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?

      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]];
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        index[i] = cgs[g[i]]+ord2[i];
        if(omap[index[i]] == 0) omap[index[i]] = i;
        else stop("Repeated values of timevar within one or more groups");
      }
      for(int j = 0; j != col; ++j) {
        ConstMatrixColumn<RTYPE> column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
          MatrixColumn<RTYPE> outj = out( _ , pos);
          if(np>0) {
            if(names) colnam[pos] = "L" + nc[p] + "." + coln[j];
            for(int i = 0; i != l; ++i) {
              if(ord2[i] >= np) {
                outj[i] = column[omap[index[i]-np]];
              } else {
                outj[i] = ff;
              }
            }
          } else if(np<0) {
            if(names) colnam[pos] = "F" + nc[p] + "." + coln[j];
            for(int i = 0; i != l; ++i) { // best loop ??
              if(ord2[i] < gsv[g[i]-1]+np) {
                outj[i] = column[omap[index[i]-np]];
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
  if(ns != 1) out.attr("dim") = Dimension(l, col*ns);
  if(names) {
    out.attr("dimnames") = List::create(rownames(x), colnam); // colnames(out) = colnam deletes row names !!!
  } else if(ns != 1) {
    out.attr("dimnames") = R_NilValue;
  }
  return out;
}


// template <>
// Vector<CPLXSXP> flagleadCppImpl(Vector<CPLXSXP> x, IntegerVector n, SEXP fill,
//                               int ng, IntegerVector g, SEXP gs, SEXP t) {
//   stop("Not supported SEXP type!");
// }

template <>
Matrix<VECSXP> flagleadmCppImpl(const Matrix<VECSXP>& x, const IntegerVector& n, const SEXP& fill,
                                int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Matrix<RAWSXP> flagleadmCppImpl(const Matrix<RAWSXP>& x, const IntegerVector& n, const SEXP& fill,
                                int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

template <>
Matrix<EXPRSXP> flagleadmCppImpl(const Matrix<EXPRSXP>& x, const IntegerVector& n, const SEXP& fill,
                                 int ng, const IntegerVector& g, const SEXP& gs, const SEXP& t, bool names) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP flagleadmCpp(SEXP x, IntegerVector n = 1, SEXP fill = R_NilValue,
                  int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue, bool names = true){
  RCPP_RETURN_MATRIX(flagleadmCppImpl, x, n, fill, ng, g, gs, t, names);
}





// [[Rcpp::export]]
List flagleadlCpp(const List& x, const IntegerVector& n = 1, const SEXP& fill = R_NilValue,
                  int ng = 0, const IntegerVector& g = 0, const SEXP& gs = R_NilValue,
                  const SEXP& t = R_NilValue, bool names = true) {

  bool lfill = fill == R_NilValue;
  if(!lfill && TYPEOF(fill) == LGLSXP) {
    LogicalVector f = fill;
    lfill = f[0] == NA_LOGICAL;
  }
  int l = x.size(), ns = n.size(), pos = INT_MAX;
  List out(l*ns);
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
  CharacterVector na = names ? coln_check(x.attr("names")) : NA_STRING;
  if(names && na[0] == NA_STRING) names = false;

  if(ng == 0) { // No groups
    if(Rf_isNull(t)) { // Ordered data
      for(int j = 0; j != l; ++j) {
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          int row = column.size();
          double ff = lfill ? NA_REAL : as<double>(fill);
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > row) stop("lag-length exceeds length of vector");
            if(np>0) {
              NumericVector outjp = no_init_vector(row);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[i++] = ff;
              for( ; i != row; ++i) outjp[i] = column[i - np];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              NumericVector outjp = no_init_vector(row);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int i = row, st = row+np;
              while(i != st) outjp[--i] = ff;
              for( ; i--; ) outjp[i] = column[i - np];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          int row = column.size();
          int ff = lfill ? NA_INTEGER : as<int>(fill);
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > row) stop("lag-length exceeds length of vector");
            if(np>0) {
              IntegerVector outjp = no_init_vector(row);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[i++] = ff;
              for( ; i != row; ++i) outjp[i] = column[i - np];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              IntegerVector outjp = no_init_vector(row);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int i = row, st = row+np;
              while(i != st) outjp[--i] = ff;
              for( ; i--; ) outjp[i] = column[i - np];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          int row = column.size();
          String ff = lfill ? NA_STRING : as<String>(fill); // String ??
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > row) stop("lag-length exceeds length of vector");
            if(np>0) {
              CharacterVector outjp = no_init_vector(row);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[i++] = ff;
              for( ; i != row; ++i) outjp[i] = column[i - np];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              CharacterVector outjp = no_init_vector(row);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int i = row, st = row+np;
              while(i != st) outjp[--i] = ff;
              for( ; i--; ) outjp[i] = column[i - np];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          int row = column.size();
          auto ff = lfill ? NA_LOGICAL : as<bool>(fill);
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > row) stop("lag-length exceeds length of vector");
            if(np>0) {
              LogicalVector outjp = no_init_vector(row);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[i++] = ff;
              for( ; i != row; ++i) outjp[i] = column[i - np];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              LogicalVector outjp = no_init_vector(row);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int i = row, st = row+np;
              while(i != st) outjp[--i] = ff;
              for( ; i--; ) outjp[i] = column[i - np];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
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
      int os = ord.size();
      LogicalVector ocheck(os, true);
      IntegerVector omap = no_init_vector(os);
      for(int i = 0; i != os; ++i) {
        if(ord[i] > os || ord[i] < 1) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i; // Note: omap is the same as order(ord) !!
        } else {
          stop("Repeated values in timevar");
        }
      }
      for(int j = 0; j != l; ++j) {
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          if(os != column.size()) stop("length(x) must match length(t)");
          double ff = lfill ? NA_REAL : as<double>(fill);
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > os) stop("lag-length exceeds length of vector");
            if(np>0) {
              NumericVector outjp = no_init_vector(os);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[omap[i++]] = ff;
              for( ; i != os; ++i) outjp[omap[i]] = column[omap[i - np]];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              NumericVector outjp = no_init_vector(os);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int st = os+np, i = os;
              while(i != st) outjp[omap[--i]] = ff;
              for( ; i--; ) outjp[omap[i]] = column[omap[i - np]];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          if(os != column.size()) stop("length(x) must match length(t)");
          int ff = lfill ? NA_INTEGER : as<int>(fill);
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > os) stop("lag-length exceeds length of vector");
            if(np>0) {
              IntegerVector outjp = no_init_vector(os);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[omap[i++]] = ff;
              for( ; i != os; ++i) outjp[omap[i]] = column[omap[i - np]];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              IntegerVector outjp = no_init_vector(os);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int st = os+np, i = os;
              while(i != st) outjp[omap[--i]] = ff;
              for( ; i--; ) outjp[omap[i]] = column[omap[i - np]];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          if(os != column.size()) stop("length(x) must match length(t)");
          String ff = lfill ? NA_STRING : as<String>(fill); // String ??
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > os) stop("lag-length exceeds length of vector");
            if(np>0) {
              CharacterVector outjp = no_init_vector(os);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[omap[i++]] = ff;
              for( ; i != os; ++i) outjp[omap[i]] = column[omap[i - np]];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              CharacterVector outjp = no_init_vector(os);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int st = os+np, i = os;
              while(i != st) outjp[omap[--i]] = ff;
              for( ; i--; ) outjp[omap[i]] = column[omap[i - np]];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          if(os != column.size()) stop("length(x) must match length(t)");
          auto ff = lfill ? NA_LOGICAL : as<bool>(fill);
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > os) stop("lag-length exceeds length of vector");
            if(np>0) {
              LogicalVector outjp = no_init_vector(os);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              int i = 0;
              while(i != np) outjp[omap[i++]] = ff;
              for( ; i != os; ++i) outjp[omap[i]] = column[omap[i - np]];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              LogicalVector outjp = no_init_vector(os);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              int st = os+np, i = os;
              while(i != st) outjp[omap[--i]] = ff;
              for( ; i--; ) outjp[omap[i]] = column[omap[i - np]];
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
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
    int gss = g.size(), ags = gss/ng, ngp = ng+1;
    if(Rf_isNull(t)) { // Ordered data
      std::vector<int> seen(ngp); // int seen[ngp], memsize = sizeof(int)*ngp;
      for(int j = 0; j != l; ++j) {
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          double ff = lfill ? NA_REAL : as<double>(fill);
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
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
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          int ff = lfill ? NA_INTEGER : as<int>(fill);
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
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
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          String ff = lfill ? NA_STRING : as<String>(fill); // String ??
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
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
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          auto ff = lfill ? NA_LOGICAL : as<bool>(fill);
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
            if(np>0) {
              LogicalVector outjp = no_init_vector(gss);
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
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              LogicalVector outjp = no_init_vector(gss);
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
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
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
      IntegerVector min(ngp, INT_MAX); // Necessary !!!
      IntegerVector gsv = Rf_isNull(gs) ? IntegerVector(ng) : as<IntegerVector>(gs); // = no_init_vector(ng); // NULL; gives compiler warning
      IntegerVector ord2 = no_init_vector(gss); // See flag.cpp for any improvements on this code !!
      if(Rf_isNull(gs)) {
        // gsv = IntegerVector(ng);
        // std::fill(gsv.begin(), gsv.end(), 0);
        for(int i = 0; i != gss; ++i) {
          ++gsv[g[i]-1];
          if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
        }
      } else {
        // gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != gss; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
      }
      IntegerVector omap(gss), cgs = no_init_vector(ngp), index = no_init_vector(gss);
      // int cgs[ngp], index[gss]; // See flag.cpp for any improvements on this code !!
      cgs[1] = 0;
      for(int i = 1; i != ng; ++i) {
        cgs[i+1] = cgs[i] + gsv[i-1]; // or get "starts from forderv"
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?
      }
      if(min[ng] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?

      for(int i = 0; i != gss; ++i) {
        ord2[i] = ord[i] - min[g[i]];
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        index[i] = cgs[g[i]]+ord2[i];
        if(omap[index[i]] == 0) omap[index[i]] = i;
        else stop("Repeated values of timevar within one or more groups");
      }
      for(int j = 0; j != l; ++j) {
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          double ff = lfill ? NA_REAL : as<double>(fill);
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
            if(np>0) {
              NumericVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) {
                if(ord2[i] >= np) {
                  outjp[i] = column[omap[index[i]-np]];
                } else {
                  outjp[i] = ff;
                }
              }
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              NumericVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) { // best loop ??
                if(ord2[i] < gsv[g[i]-1]+np) {
                  outjp[i] = column[omap[index[i]-np]];
                } else {
                  outjp[i] = ff;
                }
              }
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          int ff = lfill ? NA_INTEGER : as<int>(fill);
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
            if(np>0) {
              IntegerVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) {
                if(ord2[i] >= np) {
                  outjp[i] = column[omap[index[i]-np]];
                } else {
                  outjp[i] = ff;
                }
              }
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              IntegerVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) { // best loop ??
                if(ord2[i] < gsv[g[i]-1]+np) {
                  outjp[i] = column[omap[index[i]-np]];
                } else {
                  outjp[i] = ff;
                }
              }
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          String ff = lfill ? NA_STRING : as<String>(fill); // String ??
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
            if(np>0) {
              CharacterVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) {
                if(ord2[i] >= np) {
                  outjp[i] = column[omap[index[i]-np]];
                } else {
                  outjp[i] = ff;
                }
              }
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              CharacterVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) { // best loop ??
                if(ord2[i] < gsv[g[i]-1]+np) {
                  outjp[i] = column[omap[index[i]-np]];
                } else {
                  outjp[i] = ff;
                }
              }
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else {
              if(names) nam[pos] = na[j];
              out[pos] = column;
            }
            ++pos;
          }
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          auto ff = lfill ? NA_LOGICAL : as<bool>(fill);
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) warning("lag-length exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
            if(np>0) {
              LogicalVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) {
                if(ord2[i] >= np) {
                  outjp[i] = column[omap[index[i]-np]];
                } else {
                  outjp[i] = ff;
                }
              }
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
              out[pos] = outjp;
            } else if(np<0) {
              LogicalVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "F" + nc[p] + "." + na[j];
              for(int i = 0; i != gss; ++i) { // best loop ??
                if(ord2[i] < gsv[g[i]-1]+np) {
                  outjp[i] = column[omap[index[i]-np]];
                } else {
                  outjp[i] = ff;
                }
              }
              SHALLOW_DUPLICATE_ATTRIB(outjp, column);
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
    out.attr("names") = nam;
  } else {
    if(ns != 1) out.attr("names") = R_NilValue;
  }
  return out;
}

