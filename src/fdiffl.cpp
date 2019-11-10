// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Note: for comments and innovations see fdiff.cpp or flaga.cpp !!

inline SEXP names_check(SEXP x) { 
  if(Rf_isNull(x)) return NA_STRING;     
  else return x; 
} 

// [[Rcpp::export]]
List fdifflCpp(const List& x, const IntegerVector& n = 1, const IntegerVector& diff = 1, 
               double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0, 
               const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue, bool names = true) { 
  
  int l = x.size(), ns = n.size(), ds = diff.size(), zeros = 0, pos = 0;
  IntegerVector absn = no_init_vector(ns); 
  for(int i = ns; i--; ) {
    if(n[i] == 0) ++zeros;
    if(n[i] < 0) absn[i] = -n[i];
    else absn[i] = n[i];
  }
  int ncol = ((ns-zeros)*ds+zeros)*l;
  List out(ncol); 
  CharacterVector nam = names ? no_init_vector(ncol) : no_init_vector(1);  
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING; 
  CharacterVector diffc = names ? Rf_coerceVector(diff, STRSXP) : NA_STRING;
  CharacterVector na = names ? names_check(x.attr("names")) : NA_STRING;
  if(names && na[0] == NA_STRING) names = false;
  
  if(ng == 0) { // No groups 
    if(Rf_isNull(t)) { // Ordered data
      for(int j = 0; j != l; ++j) {
        NumericVector column = x[j];
        int row = column.size();
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0], end = np*d1; 
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            if(end >= row) stop("n * diff needs to be < nrow(x)");
            NumericVector outjp = no_init_vector(row); 
            if(names) {
              if(L1) nam[pos] = "D" + diffc[0] + "." + na[j]; 
              else nam[pos] = "L" + nc[p] + "D" + diffc[0] + "." + na[j]; 
            }
            for(int i = np; i != row; ++i) outjp[i] = column[i] - column[i - np];
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1)-1; 
              for(int i = row-1; i != start; --i) outjp[i] -= outjp[i - np]; 
            }
            for(int i = end; i--; ) outjp[i] = fill;
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp);
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = np*dq;
                if(end >= row) stop("n * diff needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1)-1; 
                  for(int i = row-1; i != start; --i) outtemp[i] -= outtemp[i - np]; 
                }
                for(int i = np*L_dq; i != end; ++i) outtemp[i] = fill;
                if(names) {
                  if(L1) nam[pos] = "D" + diffc[q] + "." + na[j]; 
                  else nam[pos] = "L" + nc[p] + "D" + diffc[q] + "." + na[j]; 
                }
                out[pos++] = Rf_shallow_duplicate(outtemp); // or Rf_copyVector, Rf_shallow_duplicate, Rf_lazy_duplicate  http://mtweb.cs.ucl.ac.uk/mus/bin/install_R/R-3.1.1/src/main/duplicate.c
              } // https://rlang.r-lib.org/reference/duplicate.html
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0], end = row+np*d1; 
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
            NumericVector outjp = no_init_vector(row); 
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            if(names) {
              if(F1) nam[pos] = "FD" + diffc[0] + "." + na[j]; 
              else nam[pos] = "F" + nc[p] + "D" + diffc[0] + "." + na[j];  
            }
            for(int i = row+np; i--; ) outjp[i] = column[i] - column[i - np];
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int final = row+np*(k+1); 
              for(int i = 0; i != final; ++i) outjp[i] -= outjp[i - np]; 
            }
            for(int i = end; i != row; ++i) outjp[i] = fill; 
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp); 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = row+np*dq, start = row+np*L_dq;
                if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int final = row+np*(k+1); 
                  for(int i = 0; i != final; ++i) outtemp[i] -= outtemp[i - np]; 
                }
                for(int i = end; i != start; ++i) outtemp[i] = fill; 
                if(names) {
                  if(F1) nam[pos] = "FD" + diffc[q] + "." + na[j]; 
                  else nam[pos] = "F"+ nc[p] + "D" + diffc[q] + "." + na[j]; 
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else {
            if(names) nam[pos] = na[j];
            out[pos++] = column;
          }
        }
      }
    } else { // Unordered data: Timevar provided 
      IntegerVector ord = t;
      int os = ord.size(), omap[os]; 
      LogicalVector ocheck(os, true); 
      for(int i = 0; i != os; ++i) { 
        if(ord[i] > os) stop("t needs to be a factor or integer vector of time-periods between 1 and nrow(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i; 
        } else {
          stop("Repeated values in timevar");
        }
      }
      for(int j = 0; j != l; ++j) {
        NumericVector column = x[j];
        if(os != column.size()) stop("nrow(x) must match length(t)");
        for(int p = 0; p != ns; ++p) { 
          int np = n[p];
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0], end = np*d1; 
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            if(end >= os) stop("n * diff needs to be < nrow(x)");
            NumericVector outjp = no_init_vector(os); 
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            if(names) {
              if(L1) nam[pos] = "D" + diffc[0] + "." + na[j]; 
              else nam[pos] = "L" + nc[p] + "D" + diffc[0] + "." + na[j]; 
            }
            for(int i = np; i != os; ++i) outjp[omap[i]] = column[omap[i]] - column[omap[i - np]];
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1)-1; 
              for(int i = os-1; i != start; --i) outjp[omap[i]] -= outjp[omap[i - np]]; 
            }
            for(int i = end; i--; ) outjp[omap[i]] = fill; 
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp); 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = np*dq;
                if(end >= os) stop("n * diff needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1)-1; 
                  for(int i = os-1; i != start; --i) outtemp[omap[i]] -= outtemp[omap[i - np]]; 
                }
                for(int i = np*L_dq; i != end; ++i) outtemp[omap[i]] = fill; 
                if(names) {
                  if(L1) nam[pos] = "D" + diffc[q] + "." + na[j]; 
                  else nam[pos] = "L" + nc[p] + "D" + diffc[q] + "." + na[j]; 
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0], end = os+np*d1; 
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
            NumericVector outjp = no_init_vector(os); 
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            if(names) {
              if(F1) nam[pos] = "FD" + diffc[0] + "." + na[j]; 
              else nam[pos] = "F" + nc[p] + "D" + diffc[0] + "." + na[j];  
            }
            for(int i = os+np; i--; ) outjp[omap[i]] = column[omap[i]] - column[omap[i - np]];
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int final = os+np*(k+1); 
              for(int i = 0; i != final; ++i) outjp[omap[i]] -= outjp[omap[i - np]]; 
            }
            for(int i = end; i != os; ++i) outjp[omap[i]] = fill;
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp); 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = os+np*dq, start = os+np*L_dq;
                if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int final = os+np*(k+1); 
                  for(int i = 0; i != final; ++i) outtemp[omap[i]] -= outtemp[omap[i - np]]; 
                }
                for(int i = end; i != start; ++i) outtemp[omap[i]] = fill; 
                if(names) {
                  if(F1) nam[pos] = "FD" + diffc[q] + "." + na[j]; 
                  else nam[pos] = "F"+ nc[p] + "D" + diffc[q] + "." + na[j]; 
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else {
            if(names) nam[pos] = na[j];
            out[pos++] = column;
          }
        }
      }
    }
  } else { // With groups
    int gss = g.size(), ags = gss/ng, ngp = ng+1, maxdiff = max(diff); 
    IntegerVector gsv = NULL; 
    if(Rf_isNull(t)) { // Ordered data
      if(maxdiff != 1) { 
        if(Rf_isNull(gs)) { 
          gsv = IntegerVector(ng);
          for(int i = 0; i != gss; ++i) ++gsv[g[i]-1];
        } else {
          gsv = gs;
          if(ng != gsv.size()) stop("ng must match length(gs)"); 
        }
      }
      int seen[ngp], memsize = sizeof(int)*(ngp); 
      for(int j = 0; j != l; ++j) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("nrow(x) must match length(g)");
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0]; 
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            NumericVector outjp = no_init_vector(gss);
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            memset(seen, 0, memsize); 
            if(names) {
              if(L1) nam[pos] = "D" + diffc[0] + "." + na[j]; 
              else nam[pos] = "L" + nc[p] + "D" + diffc[0] + "." + na[j]; 
            }
            for(int i = 0; i != gss; ++i) { 
              if(seen[g[i]] == np) outjp[i] = column[i] - column[i - np];
              else {
                outjp[i] = fill;
                ++seen[g[i]];
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1); 
              memset(seen, 0, memsize); 
              for(int i = gss; i--; ) { 
                if(seen[g[i]] == gsv[g[i]-1]-start) outjp[i] = fill;
                else {
                  outjp[i] -= outjp[i - np];
                  ++seen[g[i]];
                }
              }
            }
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp); 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1); 
                  memset(seen, 0, memsize); 
                  for(int i = gss; i--; ) {
                    if(seen[g[i]] == gsv[g[i]-1]-start) outtemp[i] = fill;
                    else {
                      outtemp[i] -= outtemp[i - np];
                      ++seen[g[i]];
                    }
                  }
                }
                if(names) {
                  if(L1) nam[pos] = "D" + diffc[q] + "." + na[j]; 
                  else nam[pos] = "L" + nc[p] + "D" + diffc[q] + "." + na[j]; 
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0]; 
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            NumericVector outjp = no_init_vector(gss); 
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            memset(seen, 0, memsize);
            if(names) {
              if(F1) nam[pos] = "FD" + diffc[0] + "." + na[j]; 
              else nam[pos] = "F" + nc[p] + "D" + diffc[0] + "." + na[j];  
            }
            for(int i = gss; i--; ) { 
              if(seen[g[i]] == np) outjp[i] = column[i] - column[i - np];
              else {
                outjp[i] = fill;
                --seen[g[i]];
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1); 
              memset(seen, 0, memsize); 
              for(int i = 0; i != gss; ++i) {
                if(seen[g[i]] == gsv[g[i]-1]+start) outjp[i] = fill; 
                else {
                  outjp[i] -= outjp[i - np];
                  ++seen[g[i]];
                }
              }
            }
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp); 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1); 
                  memset(seen, 0, memsize); 
                  for(int i = 0; i != gss; ++i) {
                    if(seen[g[i]] == gsv[g[i]-1]+start) outtemp[i] = fill; 
                    else {
                      outtemp[i] -= outtemp[i - np];
                      ++seen[g[i]];
                    }
                  }
                }
                if(names) {
                  if(F1) nam[pos] = "FD" + diffc[q] + "." + na[j]; 
                  else nam[pos] = "F"+ nc[p] + "D" + diffc[q] + "." + na[j]; 
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else {
            if(names) nam[pos] = na[j];
            out[pos++] = column;
          }
        }
      }
    } else { // Unordered data: Timevar Provided
      IntegerVector ord = t; 
      if(gss != ord.size()) stop("length(g) must match length(t)");
      IntegerVector min(ngp, INT_MAX); 
      IntegerVector ord2 = no_init_vector(gss); 
      if(Rf_isNull(gs)) { 
        gsv = IntegerVector(ng);
        for(int i = 0; i != gss; ++i) { 
          ++gsv[g[i]-1];
          if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
        }
      } else {
        gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != gss; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
      }
      IntegerVector omap(gss); 
      int cgs[ngp], seen[ngp], index[gss], memsize = sizeof(int)*(ngp); 
      cgs[1] = 0;
      for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; 
      for(int i = 0; i != gss; ++i) {
        ord2[i] = ord[i] - min[g[i]]; 
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        index[i] = cgs[g[i]]+ord2[i]; // index ??
        if(omap[index[i]] == 0) omap[index[i]] = i;  
        else stop("Repeated values of timevar within one or more groups"); 
      }
      for(int j = 0; j != l; ++j) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("nrow(x) must match length(g)");
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0]; 
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            NumericVector outjp = no_init_vector(gss);
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            if(names) {
              if(L1) nam[pos] = "D" + diffc[0] + "." + na[j]; 
              else nam[pos] = "L" + nc[p] + "D" + diffc[0] + "." + na[j]; 
            }
            for(int i = 0; i != gss; ++i) { 
              if(ord2[i] >= np) {
                outjp[i] = column[i] - column[omap[index[i]-np]];
              } else {
                outjp[i] = fill;
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1); 
              memset(seen, 0, memsize); 
              for(int i = gss; i--; ) { 
                if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outjp[omap[i]] = fill;
                else {
                  outjp[omap[i]] -= outjp[omap[i - np]];
                  ++seen[g[omap[i]]];
                }
              }
            }
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp); 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1); 
                  memset(seen, 0, memsize); 
                  for(int i = gss; i--; ) {
                    if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outtemp[omap[i]] = fill;
                    else {
                      outtemp[omap[i]] -= outtemp[omap[i - np]];
                      ++seen[g[omap[i]]];
                    }
                  }
                }
                if(names) {
                  if(L1) nam[pos] = "D" + diffc[q] + "." + na[j]; 
                  else nam[pos] = "L" + nc[p] + "D" + diffc[q] + "." + na[j]; 
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0]; 
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0"); 
            NumericVector outjp = no_init_vector(gss); 
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            if(names) {
              if(F1) nam[pos] = "FD" + diffc[0] + "." + na[j]; 
              else nam[pos] = "F" + nc[p] + "D" + diffc[0] + "." + na[j];  
            }
            for(int i = 0; i != gss; ++i) { 
              if(ord2[i] < gsv[g[i]-1]+np) {
                outjp[i] = column[i] - column[omap[index[i]-np]]; 
              } else {
                outjp[i] = fill;
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1); 
              memset(seen, 0, memsize); 
              for(int i = 0; i != gss; ++i) {
                if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outjp[omap[i]] = fill; 
                else {
                  outjp[omap[i]] -= outjp[omap[i - np]];
                  ++seen[g[omap[i]]];
                }
              }
            }
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp); 
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1); 
                  memset(seen, 0, memsize); 
                  for(int i = 0; i != gss; ++i) {
                    if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outtemp[omap[i]] = fill; 
                    else {
                      outtemp[omap[i]] -= outtemp[omap[i - np]];
                      ++seen[g[omap[i]]];
                    }
                  }
                }
                if(names) {
                  if(F1) nam[pos] = "FD" + diffc[q] + "." + na[j]; 
                  else nam[pos] = "F"+ nc[p] + "D" + diffc[q] + "." + na[j]; 
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else {
            if(names) nam[pos] = na[j];
            out[pos++] = column;
          }
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  if(names) { // best way to code this ?? 
    out.attr("names") = nam; 
  } else if(ncol != l) {
    out.attr("names") = R_NilValue;
  }
  return out;
}

