// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// Still do: What if fill = NULL -> Delete !! -> nah, very difficult with multiple lags !!

inline SEXP names_check(SEXP x) { 
  if(Rf_isNull(x)) return NA_STRING;     
  else return x; // Rf_coerceVector(x, STRSXP);             
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
  int l = x.size(), ns = n.size(), pos = 0;
  List out(l*ns);
  IntegerVector absn = no_init_vector(ns);
  for(int i = 0; i != ns; ++i) {
    if(n[i]<0) absn[i] = -n[i];
    else absn[i] = n[i];
  }
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING; // NumericVector(abs(n))
  CharacterVector nam = names ? no_init_vector(l*ns) : no_init_vector(1); // what if no names ?? 
  CharacterVector na = names ? names_check(x.attr("names")) : NA_STRING;
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
        if(ord[i] > os) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
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
      int seen[ngp], memsize = sizeof(int)*ngp;
      for(int j = 0; j != l; ++j) {
        switch(TYPEOF(x[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          double ff = lfill ? NA_REAL : as<double>(fill);
          if(gss != column.size()) stop("length(x) must match length(g)");
          for(int p = 0; p != ns; ++p) {
            int np = n[p];
            if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
            if(np>0) {
              NumericVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              memset(seen, 0, memsize); 
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
              memset(seen, 0, memsize);
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
            if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
            if(np>0) {
              IntegerVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              memset(seen, 0, memsize); 
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
              memset(seen, 0, memsize);
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
            if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
            if(np>0) {
              CharacterVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              memset(seen, 0, memsize); 
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
              memset(seen, 0, memsize);
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
            if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
            if(np>0) {
              LogicalVector outjp = no_init_vector(gss);
              if(names) nam[pos] = "L" + nc[p] + "." + na[j];
              memset(seen, 0, memsize); 
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
              memset(seen, 0, memsize);
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
      IntegerVector gsv = NULL; 
      IntegerVector ord2 = no_init_vector(gss); // See flag.cpp for any improvements on this code !!
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
      int cgs[ngp], index[gss]; // See flag.cpp for any improvements on this code !!
      cgs[1] = 0;
      for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; // or get "starts from forderv"
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
            if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
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
            if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
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
            if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
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
            if(absn[p] > ags) stop("lag-length exceeds average group-size (%i)", ags);
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


// // [[Rcpp::export]] // List
// List flagleadlCpp(List x, IntegerVector n = 1, double fill = NA_REAL, 
//                   int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) { 
//   int l = x.size(), ns = n.size();
//   
//   List out(l*ns);
//   if(ng == 0) { // No groups 
//     if(Rf_isNull(t)) { // Ordered data
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(row);
//           int np = n[p];
//           if(np>0) {
//             int i = 0;
//             while(i != np) outjp[i++] = fill;
//             for( ; i != row; ++i) outjp[i] = column[i - np];
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             int i = row, st = row+np;
//             while(i != st) outjp[--i] = fill;
//             for( ; i--; ) outjp[i] = column[i - np];
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     } else { // Unordered data: Timevar Provided
//       IntegerVector ord = t; 
//       int os = ord.size();
//       LogicalVector ocheck(os, true); // check needed ???
//       for(int i = 0; i != os; ++i) {
//         if(ocheck[ord[i]-1]) ocheck[ord[i]-1] = false;
//         else stop("Repeated values in timevar");
//       }
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         if(column.size() != os) stop("length(x) must match length(t)");
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(os);
//           int np = n[p];
//           if(np>0) {
//             for(int i = 0; i != os; ++i) {
//               if(ord[i] > np) {
//                 outjp[i] = column[ord[i]-np-1]; 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             int st = os+np;
//             for(int i = 0; i != os; ++i) { // best loop ?? (i.e. fastest ??)
//               if(ord[i] <= st) {
//                 outjp[i] = column[ord[i]-np-1]; 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     }
//   } else { // With groups
//     int gss = g.size();
//     if(Rf_isNull(t)) { // Ordered data
//       int seen[ng], memsize = sizeof(int)*ng;
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         if(column.size() != gss) stop("nrow(x) must match length(g)");
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(gss);
//           int np = n[p];
//           if(np>0) {
//             memset(seen, 0, memsize); // fastest !!
//             for(int i = 0; i != gss; ++i) {  
//               if(seen[g[i]-1] == np) {
//                 outjp[i] = column[i-np];
//               } else {
//                 outjp[i] = fill;
//                 ++seen[g[i]-1];
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             memset(seen, 0, memsize);
//             for(int i = gss; i--; ) { // good?? 
//               if(seen[g[i]-1] == np) {
//                 outjp[i] = column[i-np];
//               } else {
//                 outjp[i] = fill;
//                 --seen[g[i]-1];
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t; 
//       if(gss != ord.size()) stop("length(x) must match length(t)");
//       IntegerVector min(ng, INT_MAX);
//       IntegerVector gsv = no_init_vector(ng); // No real improvements here by using C++ arrays !!
//       IntegerVector ord2 = no_init_vector(gss); 
//       IntegerVector g2 = no_init_vector(gss);  
//       if(Rf_isNull(gs)) {
//         std::fill(gsv.begin(), gsv.end(), 0); 
//         for(int i = 0; i != gss; ++i) {
//           g2[i] = g[i]-1; 
//           ++gsv[g2[i]];
//           if(ord[i] < min[g2[i]]) min[g2[i]] = ord[i]; 
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)"); 
//         for(int i = 0; i != gss; ++i) {
//           g2[i] = g[i]-1; 
//           if(ord[i] < min[g2[i]]) min[g2[i]] = ord[i];
//         }
//       }
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) omap[i] = new int[gsv[i]]{}; // or () // better using vector of vectors (faster initialization) ??
//       for(int i = 0; i != gss; ++i) {
//         ord2[i] = ord[i] - min[g2[i]]; 
//         if(ord2[i] >= gsv[g2[i]]) stop("Gaps in timevar within one or more groups");
//         if(omap[g2[i]][ord2[i]] == 0) omap[g2[i]][ord2[i]] = i; // fastest ??
//         else stop("Repeated values of timevar within one or more groups"); 
//       }
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         if(column.size() != gss) stop("nrow(x) must match length(g)");
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(gss);
//           int np = n[p];
//           if(np>0) {
//             for(int i = 0; i != gss; ++i) {
//               if(ord2[i] >= np) {
//                 outjp[i] = column[omap[g2[i]][ord2[i]-np]]; 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             for(int i = 0; i != gss; ++i) { // best loop ??
//               if(ord2[i] < gsv[g2[i]]+np) { 
//                 outjp[i] = column[omap[g2[i]][ord2[i]-np]]; 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     }
//   }
//   return out;
// }

// Second Numeric Version without comments !! -> 2D map
// // [[Rcpp::export]] // List
// List flagleadlCpp(List x, IntegerVector n = 1, double fill = NA_REAL, 
//                int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) { 
//   int l = x.size(), ns = n.size();
//   
//   List out(l*ns);
//   if(ng == 0) { // No groups 
//     if(Rf_isNull(t)) { // Ordered data
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(row);
//           int np = n[p];
//           if(np>0) {
//             int i = 0;
//             while(i != np) outjp[i++] = fill;
//             for( ; i != row; ++i) outjp[i] = column[i - np];
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             int i = row, st = row+np;
//             while(i != st) outjp[--i] = fill;
//             for( ; i--; ) outjp[i] = column[i - np];
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     } else { // Unordered data: Timevar Provided
//       IntegerVector ord = t; 
//       int os = ord.size();
//       LogicalVector ocheck(os, true); // check needed ???
//       for(int i = 0; i != os; ++i) {
//         if(ocheck[ord[i]-1]) ocheck[ord[i]-1] = false;
//         else stop("Repeated values in timevar");
//       }
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         if(column.size() != os) stop("length(x) must match length(t)");
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(os);
//           int np = n[p];
//           if(np>0) {
//             for(int i = 0; i != os; ++i) {
//               if(ord[i] > np) {
//                 outjp[i] = column[ord[i]-np-1]; 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             int st = os+np;
//             for(int i = 0; i != os; ++i) { // best loop ?? (i.e. fastest ??)
//               if(ord[i] <= st) {
//                 outjp[i] = column[ord[i]-np-1]; 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     }
//   } else { // With groups
//     int gss = g.size();
//     if(Rf_isNull(t)) { // Ordered data
//       int seen[ng], memsize = sizeof(int)*ng;
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         if(column.size() != gss) stop("nrow(x) must match length(g)");
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(gss);
//           int np = n[p];
//           if(np>0) {
//             memset(seen, 0, memsize); // fastest !!
//             for(int i = 0; i != gss; ++i) {  
//               if(seen[g[i]-1] == np) {
//                 outjp[i] = column[i-np];
//               } else {
//                 outjp[i] = fill;
//                 ++seen[g[i]-1];
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             memset(seen, 0, memsize);
//             for(int i = gss; i--; ) { // good?? 
//               if(seen[g[i]-1] == np) {
//                 outjp[i] = column[i-np];
//               } else {
//                 outjp[i] = fill;
//                 --seen[g[i]-1];
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t; 
//       if(gss != ord.size()) stop("length(x) must match length(t)");
//       IntegerVector min(ng, INT_MAX);
//       IntegerVector gsv = no_init_vector(ng); // No real improvements here by using C++ arrays !!
//       IntegerVector ord2 = no_init_vector(gss); 
//       IntegerVector g2 = no_init_vector(gss);  
//       if(Rf_isNull(gs)) {
//         std::fill(gsv.begin(), gsv.end(), 0); 
//         for(int i = 0; i != gss; ++i) {
//           g2[i] = g[i]-1; 
//           ++gsv[g2[i]];
//           if(ord[i] < min[g2[i]]) min[g2[i]] = ord[i]; 
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)"); 
//         for(int i = 0; i != gss; ++i) {
//           g2[i] = g[i]-1; 
//           if(ord[i] < min[g2[i]]) min[g2[i]] = ord[i];
//         }
//       }
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) omap[i] = new int[gsv[i]]{}; // or () // better using vector of vectors (faster initialization) ??
//       for(int i = 0; i != gss; ++i) {
//         ord2[i] = ord[i] - min[g2[i]]; 
//         if(ord2[i] >= gsv[g2[i]]) stop("Gaps in timevar within one or more groups");
//         if(omap[g2[i]][ord2[i]] == 0) omap[g2[i]][ord2[i]] = i; // fastest ??
//         else stop("Repeated values of timevar within one or more groups"); 
//       }
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         if(column.size() != gss) stop("nrow(x) must match length(g)");
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(gss);
//           int np = n[p];
//           if(np>0) {
//             for(int i = 0; i != gss; ++i) {
//               if(ord2[i] >= np) {
//                 outjp[i] = column[omap[g2[i]][ord2[i]-np]]; 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             for(int i = 0; i != gss; ++i) { // best loop ??
//               if(ord2[i] < gsv[g2[i]]+np) { 
//                 outjp[i] = column[omap[g2[i]][ord2[i]-np]]; 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     }
//   }
//   return out;
// }


// Second Numeric Version with lots of comments !! -> 2D map
// // [[Rcpp::export]] // List
// SEXP flaglCpp(List x, IntegerVector n = 1, double fill = NA_REAL, 
//                int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) { 
//   int l = x.size(), ns = n.size();
//   
//   List out(l*ns);
//   if(ng == 0) { // No groups 
//     if(Rf_isNull(t)) { // Ordered data
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(row);
//           int np = n[p];
//           if(np>0) {
//             int i = 0;
//             while(i != np) outjp[i++] = fill;
//             for( ; i != row; ++i) outjp[i] = column[i - np];
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             int i = row, st = row+np;
//             while(i != st) outjp[--i] = fill;
//             for( ; i--; ) outjp[i] = column[i - np];
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     } else { // Unordered data: Timevar Provided
//       IntegerVector ord = t; 
//       int os = ord.size();
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(row != os) stop("length(x) must match length(t)");
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(row);
//           int np = n[p];
//           if(np>0) {
//             for(int i = 0; i != l; ++i) {
//               if(ord[i] > np) {
//                 outjp[i] = column[ord[i]-np-1]; 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             int st = row+np;
//             for(int i = 0; i != l; ++i) { // best loop ?? (i.e. fastest ??)
//               if(ord[i] <= st) {
//                 outjp[i] = column[ord[i]-np-1]; 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     }
//   } else { // With groups
//     int gss = g.size();
//     if(Rf_isNull(t)) { // Ordered data
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(row != gss) stop("nrow(x) must match length(g)");
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(row);
//           int np = n[p];
//           if(np>0) {
//             IntegerVector seen(ng); // faster using array ???
//             for(int i = 0; i != row; ++i) {  
//               if(seen[g[i]-1] == np) {
//                 outjp[i] = column[i-np];
//               } else {
//                 outjp[i] = fill;
//                 ++seen[g[i]-1];
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             IntegerVector seen(ng); // faster using array ???
//             for(int i = row; i--; ) { // good?? 
//               if(seen[g[i]-1] == np) {
//                 outjp[i] = column[i-np];
//               } else {
//                 outjp[i] = fill;
//                 --seen[g[i]-1];
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t; // clone(t); // good ??
//       if(gss != ord.size()) stop("length(x) must match length(t)");
//       // IntegerVector g2 = clone(g); // faster ?? 
//       // if(ng != gs.size()) stop("ng must match length(gs)");
//       // 1st solution: Just tape full (sparse) array
//         // int mo = max(ord); // fastest?? -> Op. 
//         // int omap[ng][mo]; // But still problems in the code below !!! -> need min(ord) and max(ord) for each group !!!
//       // 2nd solution: Perhaps that is the solution-> get min and max by group ( as ID arrys) in the second loop, then allocate it using that. 
//       
//       // int min[ng], max[ng], gs_s = gs.size(); // Good?? could also mape one 2d array!!
//       IntegerVector min(ng, INT_MAX); //, max(ng); // Or largest int value !!, max culd be 1 ??
//       IntegerVector gsv = no_init_vector(ng);
//       IntegerVector ord2 = no_init_vector(gss); //  
//       IntegerVector g2 = no_init_vector(gss); //  
//       if(Rf_isNull(gs)) {
//         // int gs[ng]; // good?? same name ??
//         std::fill(gsv.begin(), gsv.end(), 0); // gsv = 0; // good?? same name ?? // std:: fill ???
//         for(int i = 0; i != gss; ++i) {
//           g2[i] = g[i]-1; // faster ??
//           ++gsv[g2[i]];
//           if(ord[i] < min[g2[i]]) min[g2[i]] = ord[i]; // Could also do only min and checp with gsv later !!
//           // if(ord[i] > max[g2[i]]) max[g2[i]] = ord[i];
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)");
//         for(int i = 0; i != gss; ++i) {
//           g2[i] = g[i]-1; // faster ??
//           if(ord[i] < min[g2[i]]) min[g2[i]] = ord[i];
//           // if(ord[i] > max[g2[i]]) max[g2[i]] = ord[i];
//         }
//       }
//      // List ooo(6);
//      // ooo[0] = gsv;
//      // ooo[1] = min;
//      // ooo[2] = max;
//      // NumericVector t1(3), t2(3);
//     // LogicalVector t3 = no_init_vector(ng);
//      // for(int i = 0; i != ng; ++i) { 
//      //   t1[i] = gsv[i];
//      //   t2[i] = max[i]-min[i]+1;
//      //   t3[i] = gsv[i] == max[i]-min[i]+1;
//      // }
//      // ooo[3] = t1;
//      // ooo[4] = t2;
//      // ooo[5] = t3; // t1 == t2;
//      // return ooo; // Worps up to here !!
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) {
//      //   t3[i] = gsv[i] == max[i]-min[i]+1;
//       //  if(t3[i]) {
//           omap[i] = new int[gsv[i]]; // better using vector of vectors (faster initialization) ??
//      //   } else {
//       //    stop("Gaps in the timevar within a group"); // +1 // also try to mape error for repeated timevar within a group !!
//      //   }
//       }
//       //return t3;
//       for(int i = 0; i != gss; ++i) {
//         ord2[i] = ord[i] - min[g2[i]]; // helpful for all later computations !! -> faster !!
//         if(ord2[i] >= gsv[g2[i]]) stop("Gaps in the timevar within one or more groups");
//         // if(omap[g2[i]][ord[i]-min[i]] != 0) stop("Repeated values of timevar within a group"); // fast ?? , right ??
//         omap[g2[i]][ord2[i]] = i;
//       }
//       // return ord2;
//       
//       // Old solution: 
//       // int **omap = new int*[ng];
//       // for(int i = 0; i != gss; ++i) omap[g2[i]][ord[i]-1] = i;
//       
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(row != gss) stop("nrow(x) must match length(g)");
//         for(int p = ns; p--; ) {
//           NumericVector outjp = no_init_vector(row);
//           int np = n[p];
//           if(np>0) {
//             for(int i = 0; i != row; ++i) {
//               if(ord2[i] >= np) {
//                 outjp[i] = column[omap[g2[i]][ord2[i]-np]]; // good ??? 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else if(np<0) {
//             for(int i = 0; i != row; ++i) { // best loop ??
//               if(ord2[i] < gsv[g2[i]]+np) { // easier way, rather a counting vector by group ???
//                 outjp[i] = column[omap[g2[i]][ord2[i]-np]]; // good ??? -> Yes, great and efficient !!!
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           } else out[j*ns+p] = column; // x[j];
//         }
//       }
//     }
//   }
//   return out;
// }


// Previous Versions: Only Multiple lags, but first loop through lags -> slower than first looping through columns !!
// // [[Rcpp::export]]
// List flaglCpp(List x, IntegerVector n = 1, double fill = NA_REAL, 
//                int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, SEXP t = R_NilValue) { 
//   int l = x.size(), ns = n.size();
//   
//   List out(l*ns);
//   if(ng == 0) { // No groups !!
//     // TODO: FIND out what is the right order with different datasets: First loop over lags or first over columns!!
//     for(int p = ns; p--; ) {
//       int np = n[p];
//       if(np>0) {
//         for(int j = 0; j != l; ++j) {
//           NumericVector column = x[j];
//           int row = column.size(), i = 0;
//           NumericVector outjp = no_init_vector(row);
//           while(i != np) outjp[i++] = fill;
//           for( ; i != row; ++i) outjp[i] = column[i - np];
//           out[j*ns+p] = outjp;
//         }
//       } else if(np<0) {
//         for(int j = 0; j != l; ++j) {
//           NumericVector column = x[j];
//           int i = column.size(), st = i+np;
//           NumericVector outjp = no_init_vector(i);
//           while(i != st) outjp[--i] = fill;
//           for( ; i--; ) outjp[i] = column[i - np];
//           out[j*ns+p] = outjp;
//         }
//       } else {
//         for(int j = 0; j != l; ++j) out[j*ns+p] = x[j];
//       }
//     }
//   } else {
//     int gss = g.size();
//     if (Rf_isNull(t)) { // Ordered data
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         if(np>0) {
//           for(int j = 0; j != l; ++j) {
//             NumericVector column = x[j];
//             int row = column.size();
//             if(row != gss) stop("nrow(x) must match length(g)");
//             NumericVector outjp = no_init_vector(row);
//             IntegerVector seen(ng);
//             for(int i = 0; i != row; ++i) {  
//               if(seen[g[i]-1] == np) {
//                 outjp[i] = column[i-np];
//               } else {
//                 outjp[i] = fill;
//                 ++seen[g[i]-1];
//               }
//             }
//             out[j*ns+p] = outjp;
//           }
//         } else if(np<0) {
//           for(int j = 0; j != l; ++j) {
//             NumericVector column = x[j];
//             int row = column.size();
//             if(row != gss) stop("nrow(x) must match length(g)");
//             NumericVector outjp = no_init_vector(row);
//             IntegerVector seen(ng);
//             for(int i = row; i--; ) { // good?? 
//               if(seen[g[i]-1] == np) {
//                 outjp[i] = column[i-np];
//               } else {
//                 outjp[i] = fill;
//                 --seen[g[i]-1];
//               }
//             }
//             out[j*ns+p] = outjp;
//           }
//         } else {
//           for(int j = 0; j != l; ++j) out[j*ns+p] = x[j];
//         }
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t; 
//       if(gss != ord.size()) stop("length(x) must match length(t)");
//       if(ng != gs.size()) stop("ng must match length(gs)");
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//       for(int i = 0; i != gss; ++i) omap[g[i]-1][ord[i]-1] = i;
//       for(int p = ns; p--; ) {
//         int np = n[p];
//         if(np>0) {
//           for(int j = 0; j != l; ++j) {
//             NumericVector column = x[j];
//             int row = column.size(); // This is actually redundant, could put in if condition and use gss !!
//             if(row != gss) stop("nrow(x) must match length(g)");
//             NumericVector outjp = no_init_vector(row);
//             for(int i = 0; i != row; ++i) {
//               if(ord[i] > np) {
//                 outjp[i] = column[omap[g[i]-1][ord[i]-np-1]]; // good ??? 
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           }
//         } else if(np<0) {
//           for(int j = 0; j != l; ++j) {
//             NumericVector column = x[j];
//             int row = column.size(); // This is actually redundant, could put in if condition and use gss !!
//             if(row != gss) stop("nrow(x) must match length(g)");
//             NumericVector outjp = no_init_vector(row);
//             for(int i = 0; i != row; ++i) { // best loop ??
//               if(ord[i] <= gs[g[i]-1]+np) { // easier way, rather a counting vector by group ???
//                 outjp[i] = column[omap[g[i]-1][ord[i]-np-1]]; // good ??? -> Yes, great and efficient !!!
//               } else {
//                 outjp[i] = fill;
//               }
//             }
//             out[j*ns+p] = outjp;
//           }
//         } else {
//           for(int j = 0; j != l; ++j) out[j*ns+p] = x[j];
//         }
//       }
//     }
//   }
//   return out;
// }


// Previous Version: Not necessary to nistinguish between one and more lags -> no performance gain !!
// // [[Rcpp::export]]
// List flaglCpp(List x, IntegerVector n = 1, double fill = NA_REAL, 
//               int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, SEXP t = R_NilValue) { 
//   int l = x.size(), ns = n.size();
//   
//   if(ns == 1) { // Only one lag !! (Need extra method ??) -> could also do other and then attr("dim") = R_NilValue in C++
//     List out(l);
//     int nn = n[0];
//   if(ng == 0) { // No groups !!
//       if(nn>0) {
//         for(int j = 0; j != l; ++j) {
//           NumericVector column = x[j];
//           int row = column.size(), i = 0;
//           NumericVector outjp = no_init_vector(row);
//           while(i != nn) outjp[i++] = fill;
//           for( ; i != row; ++i) outjp[i] = column[i - nn];
//           out[j] = outjp;
//         }
//       } else if(nn<0) {
//         for(int j = 0; j != l; ++j) {
//           NumericVector column = x[j];
//           int i = column.size(), st = i+nn;
//           NumericVector outjp = no_init_vector(i);
//           while(i != st) outjp[--i] = fill;
//           for( ; i--; ) outjp[i] = column[i - nn];
//           out[j] = outjp;
//         }
//       } else return x;
//   } else {
//     int gss = g.size();
//     if (Rf_isNull(t)) { // Ordered data
//       if(nn>0) {
//         for(int j = 0; j != l; ++j) {
//           NumericVector column = x[j];
//           int row = column.size();
//           if(row != gss) stop("nrow(x) must match length(g)");
//           NumericVector outjp = no_init_vector(row);
//           IntegerVector seen(ng);
//           for(int i = 0; i != row; ++i) {  
//             if(seen[g[i]-1] == nn) {
//               outjp[i] = column[i-nn];
//             } else {
//               outjp[i] = fill;
//               ++seen[g[i]-1];
//             }
//           }
//           out[j] = outjp;
//         }
//       } else if(nn<0) {
//         for(int j = 0; j != l; ++j) {
//           NumericVector column = x[j];
//           int row = column.size();
//           if(row != gss) stop("nrow(x) must match length(g)");
//           NumericVector outjp = no_init_vector(row);
//           IntegerVector seen(ng);
//           for(int i = row; i--; ) { // good?? 
//             if(seen[g[i]-1] == nn) {
//               outjp[i] = column[i-nn];
//             } else {
//               outjp[i] = fill;
//               --seen[g[i]-1];
//             }
//           }
//           out[j] = outjp;
//         }
//       } else return x;
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t; 
//       if(gss != ord.size()) stop("length(x) must match length(t)");
//       if(ng != gs.size()) stop("ng must match length(gs)");
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//       for(int i = 0; i != gss; ++i) omap[g[i]-1][ord[i]-1] = i;
//       if(nn>0) {
//         for(int j = 0; j != l; ++j) {
//           NumericVector column = x[j];
//           int row = column.size(); // This is actually redundant, could put in if condition and use gss !!
//           if(row != gss) stop("nrow(x) must match length(g)");
//           NumericVector outjp = no_init_vector(row);
//           for(int i = 0; i != row; ++i) {
//             if(ord[i] > nn) {
//               outjp[i] = column[omap[g[i]-1][ord[i]-nn-1]]; // good ??? 
//             } else {
//               outjp[i] = fill;
//             }
//           }
//           out[j] = outjp;
//         }
//       } else if(nn<0) {
//         for(int j = 0; j != l; ++j) {
//           NumericVector column = x[j];
//           int row = column.size(); // This is actually redundant, could put in if condition and use gss !!
//           if(row != gss) stop("nrow(x) must match length(g)");
//           NumericVector outjp = no_init_vector(row);
//           for(int i = 0; i != row; ++i) { // best loop ??
//             if(ord[i] <= gs[g[i]-1]+nn) { // easier way, rather a counting vector by group ???
//               outjp[i] = column[omap[g[i]-1][ord[i]-nn-1]]; // good ??? -> Yes, great and efficient !!!
//             } else {
//               outjp[i] = fill;
//             }
//           }
//           out[j] = outjp;
//         }
//       } else return x;
//     }
//   }
//   return out;
//   
//   } else { // Multiple lags !!
//     List out(l*ns);
//     if(ng == 0) { // No groups !!
//    // TODO: FIND out what is the right order with different datasets: First loop over lags or first over columns!!
//       for(int p = ns; p--; ) {
//           int np = n[p];
//         if(np>0) {
//           for(int j = 0; j != l; ++j) {
//             NumericVector column = x[j];
//             int row = column.size(), i = 0;
//             NumericVector outjp = no_init_vector(row);
//             while(i != np) outjp[i++] = fill;
//             for( ; i != row; ++i) outjp[i] = column[i - np];
//             out[j*ns+p] = outjp;
//           }
//         } else if(np<0) {
//           for(int j = 0; j != l; ++j) {
//             NumericVector column = x[j];
//             int i = column.size(), st = i+np;
//             NumericVector outjp = no_init_vector(i);
//             while(i != st) outjp[--i] = fill;
//             for( ; i--; ) outjp[i] = column[i - np];
//             out[j*ns+p] = outjp;
//           }
//         } else {
//           for(int j = 0; j != l; ++j) out[j*ns+p] = x[j];
//         }
//       }
//     } else {
//       int gss = g.size();
//       if (Rf_isNull(t)) { // Ordered data
//         for(int p = ns; p--; ) {
//           int np = n[p];
//           if(np>0) {
//             for(int j = 0; j != l; ++j) {
//               NumericVector column = x[j];
//               int row = column.size();
//               if(row != gss) stop("nrow(x) must match length(g)");
//               NumericVector outjp = no_init_vector(row);
//               IntegerVector seen(ng);
//               for(int i = 0; i != row; ++i) {  
//                 if(seen[g[i]-1] == np) {
//                   outjp[i] = column[i-np];
//                 } else {
//                   outjp[i] = fill;
//                   ++seen[g[i]-1];
//                 }
//               }
//               out[j*ns+p] = outjp;
//             }
//           } else if(np<0) {
//             for(int j = 0; j != l; ++j) {
//               NumericVector column = x[j];
//               int row = column.size();
//               if(row != gss) stop("nrow(x) must match length(g)");
//               NumericVector outjp = no_init_vector(row);
//               IntegerVector seen(ng);
//               for(int i = row; i--; ) { // good?? 
//                 if(seen[g[i]-1] == np) {
//                   outjp[i] = column[i-np];
//                 } else {
//                   outjp[i] = fill;
//                   --seen[g[i]-1];
//                 }
//               }
//               out[j*ns+p] = outjp;
//             }
//           } else {
//             for(int j = 0; j != l; ++j) out[j*ns+p] = x[j];
//           }
//         }
//       } else { // Unordered data: Timevar provided
//         IntegerVector ord = t; 
//         if(gss != ord.size()) stop("length(x) must match length(t)");
//         if(ng != gs.size()) stop("ng must match length(gs)");
//         int **omap = new int*[ng];
//         for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//         for(int i = 0; i != gss; ++i) omap[g[i]-1][ord[i]-1] = i;
//         for(int p = ns; p--; ) {
//           int np = n[p];
//           if(np>0) {
//             for(int j = 0; j != l; ++j) {
//               NumericVector column = x[j];
//               int row = column.size(); // This is actually redundant, could put in if condition and use gss !!
//               if(row != gss) stop("nrow(x) must match length(g)");
//               NumericVector outjp = no_init_vector(row);
//               for(int i = 0; i != row; ++i) {
//                 if(ord[i] > np) {
//                   outjp[i] = column[omap[g[i]-1][ord[i]-np-1]]; // good ??? 
//                 } else {
//                   outjp[i] = fill;
//                 }
//               }
//               out[j*ns+p] = outjp;
//             }
//           } else if(np<0) {
//             for(int j = 0; j != l; ++j) {
//               NumericVector column = x[j];
//               int row = column.size(); // This is actually redundant, could put in if condition and use gss !!
//               if(row != gss) stop("nrow(x) must match length(g)");
//               NumericVector outjp = no_init_vector(row);
//               for(int i = 0; i != row; ++i) { // best loop ??
//                 if(ord[i] <= gs[g[i]-1]+np) { // easier way, rather a counting vector by group ???
//                   outjp[i] = column[omap[g[i]-1][ord[i]-np-1]]; // good ??? -> Yes, great and efficient !!!
//                 } else {
//                   outjp[i] = fill;
//                 }
//               }
//               out[j*ns+p] = outjp;
//             }
//           } else {
//             for(int j = 0; j != l; ++j) out[j*ns+p] = x[j];
//           }
//         }
//       }
//     }
//     return out;
//   }
// }



// Previous version: Only one lag!
// // [[Rcpp::export]]
// List flaglCpp(List x, int n = 1, double fill = NA_REAL, 
//              int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, SEXP t = R_NilValue) { 
//   int l = x.size();
//   List out(l);
//   
//   if(ng == 0) { // No groups !!
//     for(int j = 0; j != l; ++j) {
//       NumericVector column = x[j];
//       int row = column.size();
//       NumericVector outjp = no_init_vector(row);
//       int i = 0;
//       while(i != n) outjp[i++] = fill;
//       for( ; i != row; ++i) outjp[i] = column[i - n];
//       out[j] = outjp;
//     }
//   } else {
//     int gss = g.size();
//     if (Rf_isNull(t)) { // Ordered data
//       // warning("Panel-lag computed without timevar: Assuming ordered data"); -> Do in main function 
//       for(int j = 0; j != l; ++j) {
//         NumericVector column = x[j];
//         int row = column.size();
//         if(row != gss) stop("nrow(x) must match length(g)");
//         NumericVector outjp = no_init_vector(row);
//         IntegerVector seenj(ng);
//         for(int i = 0; i != row; ++i) { 
//           if(seenj[g[i]-1] == n) {
//             outjp[i] = column[i-n];
//           } else {
//             outjp[i] = fill;
//             ++seenj[g[i]-1];
//           }
//         }
//         out[j] = outjp;
//       }
//     } else { // Unordered data: Timevar provided
//       IntegerVector ord = t; 
//       if(gss != ord.size()) stop("length(x) must match length(t)");
//       if(ng != gs.size()) stop("ng must match length(gs)");
//       int **omap = new int*[ng];
//       for(int i = 0; i != ng; ++i) omap[i] = new int[gs[i]]; // better using vector of vectors (faster initialization) ??
//       for(int i = 0; i != gss; ++i) omap[g[i]-1][ord[i]-1] = i;
//       for(int j = 0; j != l; ++j) {
//         NumericVector column = x[j];
//         int row = column.size(); // This is actually redundant, could put in if condition and use gss !!
//         if(row != gss) stop("nrow(x) must match length(g)");
//         NumericVector outjp = no_init_vector(row);
//         for(int i = 0; i != row; ++i) {
//           if(ord[i] > n) {
//             outjp[i] = column[omap[g[i]-1][ord[i]-n-1]]; // good ??? 
//           } else {
//             outjp[i] = fill;
//           }
//         }
//         out[j] = outjp;
//       }
//     }
//   }
//   return out;
// }
