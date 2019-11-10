// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Note: for comments and innovations see fdiff.cpp, fgrowth.cpp, flaga.cpp or fdiffa.cpp !!

inline SEXP coln_check(SEXP x) { 
  if(Rf_isNull(x)) return NA_STRING;     
  else return x; 
} 

inline double do_growth(double y, double x) {
  return (y-x)*(100/x);
}

inline double do_logdiff(double y, double x) {
  return (log(y)-log(x))*100;
}


// [[Rcpp::export]]
NumericMatrix fgrowthmCpp(const NumericMatrix& x, const IntegerVector& n = 1, const IntegerVector& diff = 1, 
                          double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0, 
                          const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue, 
                          bool logdiff = false, bool names = true) { 
  
  int l = x.nrow(), col = x.ncol(), ns = n.size(), ds = diff.size(), zeros = 0, pos = 0;
  IntegerVector absn = no_init_vector(ns); 
  for(int i = ns; i--; ) {
    if(n[i] == 0) ++zeros;
    if(n[i] < 0) absn[i] = -n[i];
    else absn[i] = n[i];
  }
  auto FUN = logdiff ? do_logdiff : do_growth;
  std::string stub, stub2; 
  if(names) {
    if(logdiff) {
      stub = "Dlog";
      stub2 = "FDlog";
    } else {
      stub = "G";
      stub2 = "FG";
    }
  } 
  int ncol = ((ns-zeros)*ds+zeros)*col;
  NumericMatrix out = no_init_matrix(l, ncol); 
  CharacterVector colnam = names ? no_init_vector(ncol) : no_init_vector(1);  
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING; 
  CharacterVector diffc = names ? Rf_coerceVector(diff, STRSXP) : NA_STRING;
  CharacterVector coln = names ? coln_check(colnames(x)) : NA_STRING;
  if(names && coln[0] == NA_STRING) names = false;
  
  if(ng == 0) { // No groups 
    if(Rf_isNull(t)) { // Ordered data
      for(int j = 0; j != col; ++j) {
        NumericMatrix::ConstColumn column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0], end = np*d1; 
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            if(end >= l) stop("n * diff needs to be < nrow(x)");
            NumericMatrix::Column outp = out( _ , pos); 
            if(names) {
              if(L1) colnam[pos] = stub + diffc[0] + "." + coln[j]; 
              else colnam[pos] = "L" + nc[p] + stub + diffc[0] + "." + coln[j]; 
            }
            ++pos;
            for(int i = np; i != l; ++i) outp[i] = FUN(column[i], column[i - np]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1)-1; 
              for(int i = l-1; i != start; --i) outp[i] = FUN(outp[i], outp[i - np]); 
            }
            for(int i = end; i--; ) outp[i] = fill; 
            if(ds > 1) {
              NumericVector outtemp = outp; 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = np*dq;
                if(end >= l) stop("n * diff needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1)-1; 
                  for(int i = l-1; i != start; --i) outtemp[i] = FUN(outtemp[i], outtemp[i - np]); 
                }
                for(int i = np*L_dq; i != end; ++i) outtemp[i] = fill; 
                out( _ , pos) = outtemp; 
                if(names) {
                  if(L1) colnam[pos] = stub + diffc[q] + "." + coln[j]; 
                  else colnam[pos] = "L" + nc[p] + stub + diffc[q] + "." + coln[j]; 
                }
                ++pos;
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0], end = l+np*d1; 
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
            NumericMatrix::Column outp = out( _ , pos); 
            if(names) {
              if(F1) colnam[pos] = stub2 + diffc[0] + "." + coln[j]; 
              else colnam[pos] = "F" + nc[p] + stub + diffc[0] + "." + coln[j];  
            }
            ++pos;
            for(int i = l+np; i--; ) outp[i] = FUN(column[i], column[i - np]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int final = l+np*(k+1); 
              for(int i = 0; i != final; ++i) outp[i] = FUN(outp[i], outp[i - np]); 
            }
            for(int i = end; i != l; ++i) outp[i] = fill; 
            if(ds > 1) {
              NumericVector outtemp = outp; 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
                if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int final = l+np*(k+1); 
                  for(int i = 0; i != final; ++i) outtemp[i] = FUN(outtemp[i], outtemp[i - np]); 
                }
                for(int i = end; i != start; ++i) outtemp[i] = fill; 
                out( _ , pos) = outtemp; 
                if(names) {
                  if(F1) colnam[pos] = stub2 + diffc[q] + "." + coln[j]; 
                  else colnam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + coln[j]; 
                }
                ++pos;
              }
            }
          } else {
            out( _ , pos) = column;
            if(names) colnam[pos] = coln[j];
            ++pos;
          }
        }
      }
    } else { // Unordered data: Timevar provided 
      IntegerVector ord = t;
      if(l != ord.size()) stop("nrow(x) must match length(t)");
      LogicalVector ocheck(l, true); 
      int omap[l]; 
      for(int i = 0; i != l; ++i) { 
        if(ord[i] > l) stop("t needs to be a factor or integer vector of time-periods between 1 and nrow(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i; 
        } else {
          stop("Repeated values in timevar");
        }
      }
      for(int j = 0; j != col; ++j) {
        NumericMatrix::ConstColumn column = x( _ , j);
        for(int p = 0; p != ns; ++p) { 
          int np = n[p];
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0], end = np*d1; 
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            if(end >= l) stop("n * diff needs to be < nrow(x)");
            NumericMatrix::Column outp = out( _ , pos); 
            if(names) {
              if(L1) colnam[pos] = stub + diffc[0] + "." + coln[j]; 
              else colnam[pos] = "L" + nc[p] + stub + diffc[0] + "." + coln[j]; 
            }
            ++pos;
            for(int i = np; i != l; ++i) outp[omap[i]] = FUN(column[omap[i]], column[omap[i - np]]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1)-1; 
              for(int i = l-1; i != start; --i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]); 
            }
            for(int i = end; i--; ) outp[omap[i]] = fill; 
            if(ds > 1) {
              NumericVector outtemp = outp; 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = np*dq;
                if(end >= l) stop("n * diff needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1)-1; 
                  for(int i = l-1; i != start; --i) outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]); 
                }
                for(int i = np*L_dq; i != end; ++i) outtemp[omap[i]] = fill; 
                out( _ , pos) = outtemp; 
                if(names) {
                  if(L1) colnam[pos] = stub + diffc[q] + "." + coln[j]; 
                  else colnam[pos] = "L" + nc[p] + stub + diffc[q] + "." + coln[j]; 
                }
                ++pos;
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0], end = l+np*d1; 
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
            NumericMatrix::Column outp = out( _ , pos); 
            if(names) {
              if(F1) colnam[pos] = stub2 + diffc[0] + "." + coln[j]; 
              else colnam[pos] = "F" + nc[p] + stub + diffc[0] + "." + coln[j];  
            }
            ++pos;
            for(int i = l+np; i--; ) outp[omap[i]] = FUN(column[omap[i]], column[omap[i - np]]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int final = l+np*(k+1); 
              for(int i = 0; i != final; ++i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]); 
            }
            for(int i = end; i != l; ++i) outp[omap[i]] = fill; 
            if(ds > 1) {
              NumericVector outtemp = outp; 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
                if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int final = l+np*(k+1); 
                  for(int i = 0; i != final; ++i) outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]); 
                }
                for(int i = end; i != start; ++i) outtemp[omap[i]] = fill; 
                out( _ , pos) = outtemp; 
                if(names) {
                  if(F1) colnam[pos] = stub2 + diffc[q] + "." + coln[j]; 
                  else colnam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + coln[j]; 
                }
                ++pos;
              }
            }
          } else {
            out( _ , pos) = column;
            if(names) colnam[pos] = coln[j];
            ++pos;
          }
        }
      }
    }
  } else { // With groups
    if(l != g.size()) stop("nrow(x) must match length(g)");
    int ags = l/ng, ngp = ng+1, maxdiff = max(diff); 
    IntegerVector gsv = NULL; 
    if(Rf_isNull(t)) { // Ordered data
      if(maxdiff != 1) { 
        if(Rf_isNull(gs)) { 
          gsv = IntegerVector(ng);
          for(int i = 0; i != l; ++i) ++gsv[g[i]-1];
        } else {
          gsv = gs;
          if(ng != gsv.size()) stop("ng must match length(gs)"); 
        }
      }
      int seen[ngp], memsize = sizeof(int)*(ngp); 
      for(int j = 0; j != col; ++j) {
        NumericMatrix::ConstColumn column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0]; 
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            NumericMatrix::Column outp = out( _ , pos); 
            memset(seen, 0, memsize); 
            if(names) {
              if(L1) colnam[pos] = stub + diffc[0] + "." + coln[j]; 
              else colnam[pos] = "L" + nc[p] + stub + diffc[0] + "." + coln[j]; 
            }
            ++pos;
            for(int i = 0; i != l; ++i) { 
              if(seen[g[i]] == np) outp[i] = FUN(column[i], column[i - np]);
              else {
                outp[i] = fill;
                ++seen[g[i]];
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1); 
              memset(seen, 0, memsize); 
              for(int i = l; i--; ) { 
                if(seen[g[i]] == gsv[g[i]-1]-start) outp[i] = fill;
                else {
                  outp[i] = FUN(outp[i], outp[i - np]);
                  ++seen[g[i]];
                }
              }
            }
            if(ds > 1) {
              NumericVector outtemp = outp; 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1); 
                  memset(seen, 0, memsize); 
                  for(int i = l; i--; ) {
                    if(seen[g[i]] == gsv[g[i]-1]-start) outtemp[i] = fill;
                    else {
                      outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                      ++seen[g[i]];
                    }
                  }
                }
                out( _ , pos) = outtemp; 
                if(names) {
                  if(L1) colnam[pos] = stub + diffc[q] + "." + coln[j]; 
                  else colnam[pos] = "L" + nc[p] + stub + diffc[q] + "." + coln[j]; 
                }
                ++pos;
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0]; 
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            NumericMatrix::Column outp = out( _ , pos); 
            memset(seen, 0, memsize);
            if(names) {
              if(F1) colnam[pos] = stub2 + diffc[0] + "." + coln[j]; 
              else colnam[pos] = "F" + nc[p] + stub + diffc[0] + "." + coln[j];  
            }
            ++pos;
            for(int i = l; i--; ) { 
              if(seen[g[i]] == np) outp[i] = FUN(column[i], column[i - np]);
              else {
                outp[i] = fill;
                --seen[g[i]];
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1); 
              memset(seen, 0, memsize); 
              for(int i = 0; i != l; ++i) {
                if(seen[g[i]] == gsv[g[i]-1]+start) outp[i] = fill; 
                else {
                  outp[i] = FUN(outp[i], outp[i - np]);
                  ++seen[g[i]];
                }
              }
            }
            if(ds > 1) {
              NumericVector outtemp = outp; 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1); 
                  memset(seen, 0, memsize); 
                  for(int i = 0; i != l; ++i) {
                    if(seen[g[i]] == gsv[g[i]-1]+start) outtemp[i] = fill; 
                    else {
                      outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                      ++seen[g[i]];
                    }
                  }
                }
                out( _ , pos) = outtemp; 
                if(names) {
                  if(F1) colnam[pos] = stub2 + diffc[q] + "." + coln[j]; 
                  else colnam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + coln[j]; 
                }
                ++pos;
              }
            }
          } else {
            out( _ , pos) = column;
            if(names) colnam[pos] = coln[j];
            ++pos;
          }
        }
      }
    } else { // Unordered data: Timevar Provided
      IntegerVector ord = t; 
      if(l != ord.size()) stop("nrow(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX); 
      IntegerVector ord2 = no_init_vector(l); 
      if(Rf_isNull(gs)) { 
        gsv = IntegerVector(ng);
        for(int i = 0; i != l; ++i) { 
          ++gsv[g[i]-1];
          if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
        }
      } else {
        gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
      }
      IntegerVector omap(l); 
      int cgs[ngp], seen[ngp], index[l], memsize = sizeof(int)*(ngp); 
      cgs[1] = 0;
      for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; 
      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]]; 
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        index[i] = cgs[g[i]]+ord2[i]; // index ??
        if(omap[index[i]] == 0) omap[index[i]] = i;  
        else stop("Repeated values of timevar within one or more groups"); 
      }
      for(int j = 0; j != col; ++j) {
        NumericMatrix::ConstColumn column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0]; 
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            NumericMatrix::Column outp = out( _ , pos); 
            if(names) {
              if(L1) colnam[pos] = stub + diffc[0] + "." + coln[j]; 
              else colnam[pos] = "L" + nc[p] + stub + diffc[0] + "." + coln[j]; 
            }
            ++pos;
            for(int i = 0; i != l; ++i) { 
              if(ord2[i] >= np) {
                outp[i] = FUN(column[i], column[omap[index[i]-np]]);
              } else {
                outp[i] = fill;
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1); 
              memset(seen, 0, memsize); 
              for(int i = l; i--; ) { 
                if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outp[omap[i]] = fill;
                else {
                  outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
                  ++seen[g[omap[i]]];
                }
              }
            }
            if(ds > 1) {
              NumericVector outtemp = outp; 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1); 
                  memset(seen, 0, memsize); 
                  for(int i = l; i--; ) {
                    if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outtemp[omap[i]] = fill;
                    else {
                      outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                      ++seen[g[omap[i]]];
                    }
                  }
                }
                out( _ , pos) = outtemp; 
                if(names) {
                  if(L1) colnam[pos] = stub + diffc[q] + "." + coln[j]; 
                  else colnam[pos] = "L" + nc[p] + stub + diffc[q] + "." + coln[j]; 
                }
                ++pos;
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0]; 
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            NumericMatrix::Column outp = out( _ , pos); 
            if(names) {
              if(F1) colnam[pos] = stub2 + diffc[0] + "." + coln[j]; 
              else colnam[pos] = "F" + nc[p] + stub + diffc[0] + "." + coln[j];  
            }
            ++pos;
            for(int i = 0; i != l; ++i) { 
              if(ord2[i] < gsv[g[i]-1]+np) {
                outp[i] = FUN(column[i], column[omap[index[i]-np]]); 
              } else {
                outp[i] = fill;
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1); 
              memset(seen, 0, memsize); 
              for(int i = 0; i != l; ++i) {
                if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outp[omap[i]] = fill; 
                else {
                  outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
                  ++seen[g[omap[i]]];
                }
              }
            }
            if(ds > 1) {
              NumericVector outtemp = outp; 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1); 
                  memset(seen, 0, memsize); 
                  for(int i = 0; i != l; ++i) {
                    if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outtemp[omap[i]] = fill; 
                    else {
                      outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                      ++seen[g[omap[i]]];
                    }
                  }
                }
                out( _ , pos) = outtemp; 
                if(names) {
                  if(F1) colnam[pos] = stub2 + diffc[q] + "." + coln[j]; 
                  else colnam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + coln[j]; 
                }
                ++pos;
              }
            }
          } else {
            out( _ , pos) = column;
            if(names) colnam[pos] = coln[j];
            ++pos;
          }
        }
      }
    }
  }
  // Previous Solution:  
  // if(names) {
  //   out.attr("dimnames") = List::create(rownames(x), colnam); 
  // } else {
  //   if(ns*ds == 1) DUPLICATE_ATTRIB(out, x);
  //   // else rownames(out) = rownames(x); // redundant !! 
  // }
  
  DUPLICATE_ATTRIB(out, x);
  if(ncol != col) out.attr("dim") = Dimension(l, ncol);
  if(names) {
    out.attr("dimnames") = List::create(rownames(x), colnam); //colnames(out) = colnam; also deletes rownames !!
  } else if(ncol != col) {
    out.attr("dimnames") = R_NilValue;
  }
  
  return out;
}
