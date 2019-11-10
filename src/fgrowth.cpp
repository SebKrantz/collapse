// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// For major comments and changes see fdiff.cpp code !!
// inline double do_diff(double y, double x) { // slower than fdiffCpp !!!!!!!!!!!
//   return y-x;
// }

inline double do_growth(double y, double x) {
  return (y-x)*(100/x);
}

inline double do_logdiff(double y, double x) {
  return (log(y)-log(x))*100;
}


// [[Rcpp::export]]
NumericVector fgrowthCpp(const NumericVector& x, const IntegerVector& n = 1, const IntegerVector& diff = 1, 
                         double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0, 
                         const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue, 
                         bool logdiff = false, bool names = true) { 
  
  int l = x.size(), ns = n.size(), ds = diff.size(), zeros = 0, pos = 0;
  IntegerVector absn = no_init_vector(ns); 
  for(int i = ns; i--; ) {
    if(n[i] == 0) ++zeros;
    if(n[i] < 0) absn[i] = -n[i];
    else absn[i] = n[i];
  }
  auto FUN = logdiff ? do_logdiff : do_growth;
  // double (*FUN)(double, double); // https://stackoverflow.com/questions/5582869/how-do-i-store-a-function-to-a-variable
  std::string stub, stub2, stub3; // String -> error !!
  if(names) {
    if(logdiff) {
      //    FUN = do_logdiff;
      stub = ".Dlog";
      stub2 = "Dlog";
      stub3 = ".FDlog";
    } else {
      //    FUN = do_growth;
      stub = ".G";
      stub2 = "G";
      stub3 = ".FG";
    }
  } // else {
  // FUN = logdiff ? do_logdiff : do_growth;
  // }
  int ncol = (ns-zeros)*ds+zeros;
  if(ncol == 1) names = false;
  NumericMatrix out = no_init_matrix(l, ncol); 
  CharacterVector colnam = names ? no_init_vector(ncol) : no_init_vector(1); 
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING; 
  CharacterVector diffc = names ? Rf_coerceVector(diff, STRSXP) : NA_STRING; 
  if(ng == 0) { // No groups 
    if(Rf_isNull(t)) { // Ordered data 
      for(int p = 0; p != ns; ++p) {
        int np = n[p];
        if(np>0) { // Positive lagged and iterated differences
          int d1 = diff[0], end = np*d1; 
          bool L1 = np == 1;
          if(d1 < 1) stop("diff must be a vector of integers > 0"); 
          if(end >= l) stop("n * diff needs to be < length(x)");
          NumericMatrix::Column outp = out( _ , pos); 
          if(names) {
            if(L1) colnam[pos] = stub + diffc[0]; 
            else colnam[pos] = ".L" + nc[p] + stub2 + diffc[0]; 
          }
          ++pos;
          for(int i = np; i != l; ++i) outp[i] = FUN(x[i], x[i - np]);
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1)-1; 
            for(int i = l-1; i != start; --i) outp[i] = FUN(outp[i], outp[i - np]); 
          }
          for(int i = end; i--; ) outp[i] = fill; 
          if(ds > 1) {
            NumericVector outtemp = outp; 
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1], end = np*dq;
              if(end >= l) stop("n * diff needs to be < length(x)");
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int start = np*(k+1)-1; 
                for(int i = l-1; i != start; --i) outtemp[i] = FUN(outtemp[i], outtemp[i - np]); 
              }
              for(int i = np*L_dq; i != end; ++i) outtemp[i] = fill; 
              out( _ , pos) = outtemp; 
              if(names) {
                if(L1) colnam[pos] = stub + diffc[q]; 
                else colnam[pos] = ".L" + nc[p] + stub2 + diffc[q]; 
              }
              ++pos;
            }
          }
        } else if(np<0) { // (Negative) leaded and iterated differences
          int d1 = diff[0], end = l+np*d1; 
          bool F1 = np == -1;
          if(d1 < 1) stop("diff must be a vector of integers > 0"); 
          if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
          NumericMatrix::Column outp = out( _ , pos); 
          if(names) {
            if(F1) colnam[pos] = stub3 + diffc[0]; 
            else colnam[pos] = ".F" + nc[p] + stub2 + diffc[0];  
          }
          ++pos;
          for(int i = l+np; i--; ) outp[i] = FUN(x[i], x[i - np]);
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int final = l+np*(k+1); 
            for(int i = 0; i != final; ++i) outp[i] = FUN(outp[i], outp[i - np]); 
          }
          for(int i = end; i != l; ++i) outp[i] = fill; 
          if(ds > 1) {
            NumericVector outtemp = outp; 
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
              if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int final = l+np*(k+1); 
                for(int i = 0; i != final; ++i) outtemp[i] = FUN(outtemp[i], outtemp[i - np]); 
              }
              for(int i = end; i != start; ++i) outtemp[i] = fill; 
              out( _ , pos) = outtemp; 
              if(names) {
                if(F1) colnam[pos] = stub3 + diffc[q]; 
                else colnam[pos] = ".F"+ nc[p] + stub2 + diffc[q]; 
              }
              ++pos;
            }
          }
        } else {
          out( _ , pos) = x;
          if(names) colnam[pos] = ".--";
          ++pos;
        }
      }
    } else { // Unordered data: Timevar provided 
      IntegerVector ord = t;
      if(l != ord.size()) stop("length(x) must match length(t)");
      LogicalVector ocheck(l, true); 
      int omap[l]; 
      for(int i = 0; i != l; ++i) { 
        if(ord[i] > l) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i; 
        } else {
          stop("Repeated values in timevar");
        }
      }
      for(int p = 0; p != ns; ++p) { 
        int np = n[p];
        if(np>0) { // Positive lagged and iterated differences
          int d1 = diff[0], end = np*d1; 
          bool L1 = np == 1;
          if(d1 < 1) stop("diff must be a vector of integers > 0"); 
          if(end >= l) stop("n * diff needs to be < length(x)");
          NumericMatrix::Column outp = out( _ , pos); 
          if(names) {
            if(L1) colnam[pos] = stub + diffc[0]; 
            else colnam[pos] = ".L" + nc[p] + stub2 + diffc[0]; 
          }
          ++pos;
          for(int i = np; i != l; ++i) outp[omap[i]] = FUN(x[omap[i]], x[omap[i - np]]);
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1)-1; 
            for(int i = l-1; i != start; --i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]); 
          }
          for(int i = end; i--; ) outp[omap[i]] = fill; 
          if(ds > 1) {
            NumericVector outtemp = outp; 
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1], end = np*dq;
              if(end >= l) stop("n * diff needs to be < length(x)");
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int start = np*(k+1)-1; 
                for(int i = l-1; i != start; --i) outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]); 
              }
              for(int i = np*L_dq; i != end; ++i) outtemp[omap[i]] = fill; 
              out( _ , pos) = outtemp; 
              if(names) {
                if(L1) colnam[pos] = stub + diffc[q]; 
                else colnam[pos] = ".L" + nc[p] + stub2 + diffc[q]; 
              }
              ++pos;
            }
          }
        } else if(np<0) { // (Negative) leaded and iterated differences
          int d1 = diff[0], end = l+np*d1; 
          bool F1 = np == -1;
          if(d1 < 1) stop("diff must be a vector of integers > 0"); 
          if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
          NumericMatrix::Column outp = out( _ , pos); 
          if(names) {
            if(F1) colnam[pos] = stub3 + diffc[0]; 
            else colnam[pos] = ".F" + nc[p] + stub2 + diffc[0]; 
          }
          ++pos;
          for(int i = l+np; i--; ) outp[omap[i]] = FUN(x[omap[i]], x[omap[i - np]]);
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int final = l+np*(k+1); 
            for(int i = 0; i != final; ++i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]); 
          }
          for(int i = end; i != l; ++i) outp[omap[i]] = fill; 
          if(ds > 1) {
            NumericVector outtemp = outp; 
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
              if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int final = l+np*(k+1); 
                for(int i = 0; i != final; ++i) outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]); 
              }
              for(int i = end; i != start; ++i) outtemp[omap[i]] = fill; 
              out( _ , pos) = outtemp; 
              if(names) {
                if(F1) colnam[pos] = stub3 + diffc[q]; 
                else colnam[pos] = ".F"+ nc[p] + stub2 + diffc[q]; 
              }
              ++pos;
            }
          }
        } else {
          out( _ , pos) = x;
          if(names) colnam[pos] = ".--";
          ++pos;
        }
      }
    }
  } else {
    if(l != g.size()) stop("length(x) must match length(g)");
    int ags = l/ng, ngp = ng+1, maxdiff = max(diff);  
    IntegerVector gsv = NULL; 
    if(Rf_isNull(t)) { 
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
            if(L1) colnam[pos] = stub + diffc[0]; 
            else colnam[pos] = ".L" + nc[p] + stub2 + diffc[0]; 
          }
          ++pos;
          for(int i = 0; i != l; ++i) { 
            if(seen[g[i]] == np) outp[i] = FUN(x[i], x[i - np]);
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
                int start = np*(k+1); // Right ?? -> seems so!!
                memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
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
                if(L1) colnam[pos] = stub + diffc[q]; 
                else colnam[pos] = ".L" + nc[p] + stub2 + diffc[q]; 
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
            if(F1) colnam[pos] = stub3 + diffc[0]; 
            else colnam[pos] = ".F" + nc[p] + stub2 + diffc[0]; 
          }
          ++pos;
          for(int i = l; i--; ) {  
            if(seen[g[i]] == np) outp[i] = FUN(x[i], x[i - np]);
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
                if(F1) colnam[pos] = stub3 + diffc[q]; 
                else colnam[pos] = ".F"+ nc[p] + stub2 + diffc[q]; 
              }
              ++pos;
            }
          }
        } else {
          out( _ , pos) = x;
          if(names) colnam[pos] = ".--";
          ++pos;
        }
      }
    } else { // Unordered data: Timevar Provided
      IntegerVector ord = t; 
      if(l != ord.size()) stop("length(x) must match length(t)");
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
      int cgs[ngp], seen[ngp], memsize = sizeof(int)*(ngp); 
      cgs[1] = 0;
      for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; 
      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]]; 
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        if(omap[cgs[g[i]]+ord2[i]] == 0) omap[cgs[g[i]]+ord2[i]] = i; 
        else stop("Repeated values of timevar within one or more groups"); 
      }
      for(int p = 0; p != ns; ++p) {
        int np = n[p];
        if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
        if(np>0) { // Positive lagged and iterated differences
          int d1 = diff[0]; 
          bool L1 = np == 1;
          if(d1 < 1) stop("diff must be a vector of integers > 0"); 
          NumericMatrix::Column outp = out( _ , pos);
          if(names) {
            if(L1) colnam[pos] = stub + diffc[0]; 
            else colnam[pos] = ".L" + nc[p] + stub2 + diffc[0]; 
          }
          ++pos;
          for(int i = 0; i != l; ++i) { 
            if(ord2[i] >= np) {
              outp[i] = FUN(x[i], x[omap[cgs[g[i]]+ord2[i]-np]]);
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
                if(L1) colnam[pos] = stub + diffc[q]; 
                else colnam[pos] = ".L" + nc[p] + stub2 + diffc[q]; 
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
            if(F1) colnam[pos] = stub3 + diffc[0]; 
            else colnam[pos] = ".F" + nc[p] + stub2 + diffc[0]; 
          }
          ++pos;
          for(int i = 0; i != l; ++i) { 
            if(ord2[i] < gsv[g[i]-1]+np) {
              outp[i] = FUN(x[i], x[omap[cgs[g[i]]+ord2[i]-np]]); 
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
                if(F1) colnam[pos] = stub3 + diffc[q]; 
                else colnam[pos] = ".F"+ nc[p] + stub2 + diffc[q]; 
              }
              ++pos;
            }
          }
        } else {
          out( _ , pos) = x;
          if(names) colnam[pos] = ".--";
          ++pos;
        }
      }
    }
  }
  // Previous Version  
  // if(ncol == 1) DUPLICATE_ATTRIB(out, x);
  // else if(names) out.attr("dimnames") = List::create(x.attr("names"), colnam);
  
  DUPLICATE_ATTRIB(out, x); 
  if(ncol != 1) { 
    if(x.hasAttribute("names")) out.attr("names") = R_NilValue; 
    out.attr("dim") = Dimension(l, ncol);
    if(names) out.attr("dimnames") = List::create(x.attr("names"), colnam);
  }
  
  return out;
}




// Second Version: two choices: growth and logdiff, but no names argument !! (note: viod update function doesn't work and is not faster !!)
// inline double growth1(double y, double x) {
//   return (y-x)*(100/x);
// }
// 
// inline double growth2(double y, double x) {
//   return (log(y)-log(x))*100;
// }
// 
// // [[Rcpp::export]]
// NumericVector fgrowthCpp(NumericVector x, IntegerVector n = 1, IntegerVector diff = 1, double fill = NA_REAL, 
//                        int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue, bool logdiff = false) { 
//   
//   int l = x.size(), ns = n.size(), ds = diff.size(), zeros = 0, pos = 0;
//   IntegerVector absn = no_init_vector(ns); 
//   for(int i = ns; i--; ) {
//     if(n[i] == 0) ++zeros;
//     if(n[i] < 0) absn[i] = -n[i];
//     else absn[i] = n[i];
//   }
//   
//   // #if logdiff
//   // #define growth(y,x) (std::log(y)-std::log(x))*(100) 
//   // #else
//   // #define growth(y,x) ((y)-(x))*((100)/(x)) 
//   // #endif
//   
//   // if(logdiff) {
//   //   #define growth(y,x) (std::log(y)-std::log(x))*(100) 
//   // } else {
//   //   #define growth(y,x) ((y)-(x))*((100)/(x)) 
//   // }
//   
//   auto growth = logdiff ? growth2 : growth1;
//   
//   int ncol = (ns-zeros)*ds+zeros;
//   NumericMatrix out = no_init_matrix(l, ncol); 
//   CharacterVector colnam = no_init_vector(ncol); 
//   CharacterVector nc = Rf_coerceVector(absn, STRSXP); 
//   CharacterVector diffc = Rf_coerceVector(diff, STRSXP); 
//   if(ng == 0) { // No groups 
//     if(Rf_isNull(t)) { // Ordered data 
//       for(int p = 0; p != ns; ++p) {
//         int np = n[p];
//         if(np>0) { // Positive lagged and iterated differences
//           int d1 = diff[0], end = np*d1; 
//           bool L1 = np == 1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           if(end >= l) stop("n * diff needs to be < length(x)");
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(L1) colnam[pos++] = ".G" + diffc[0]; 
//           else colnam[pos++] = ".L" + nc[p] + "G" + diffc[0]; 
//           for(int i = np; i != l; ++i) outp[i] = growth(x[i], x[i - np]);
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1)-1; 
//             for(int i = l-1; i != start; --i) outp[i] = growth(outp[i], outp[i - np]); 
//           }
//           for(int i = end; i--; ) outp[i] = fill; 
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1], end = np*dq;
//               if(end >= l) stop("n * diff needs to be < length(x)");
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1)-1; 
//                 for(int i = l-1; i != start; --i) outtemp[i] = growth(outtemp[i], outtemp[i - np]); 
//               }
//               for(int i = np*L_dq; i != end; ++i) outtemp[i] = fill; 
//               out( _ , pos) = outtemp; 
//               if(L1) colnam[pos++] = ".G" + diffc[q]; 
//               else colnam[pos++] = ".L" + nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else if(np<0) { // (Negative) leaded and iterated differences
//           int d1 = diff[0], end = l+np*d1; 
//           bool F1 = np == -1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(F1) colnam[pos++] = ".FG" + diffc[0]; 
//           else colnam[pos++] = ".F" + nc[p] + "G" + diffc[0];  
//           for(int i = l+np; i--; ) outp[i] = growth(x[i], x[i - np]);
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int final = l+np*(k+1); 
//             for(int i = 0; i != final; ++i) outp[i] = growth(outp[i], outp[i - np]); 
//           }
//           for(int i = end; i != l; ++i) outp[i] = fill; 
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
//               if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int final = l+np*(k+1); 
//                 for(int i = 0; i != final; ++i) outtemp[i] = growth(outtemp[i], outtemp[i - np]); 
//               }
//               for(int i = end; i != start; ++i) outtemp[i] = fill; 
//               out( _ , pos) = outtemp; 
//               if(F1) colnam[pos++] = ".FG" + diffc[q]; 
//               else colnam[pos++] = ".F"+ nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else {
//           out( _ , pos) = x;
//           colnam[pos++] = ".--";
//         }
//       }
//     } else { // Unordered data: Timevar provided 
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       LogicalVector ocheck(l, true); 
//       int omap[l]; 
//       for(int i = 0; i != l; ++i) { 
//         if(ord[i] > l) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
//         if(ocheck[ord[i]-1]) {
//           ocheck[ord[i]-1] = false;
//           omap[ord[i]-1] = i; 
//         } else {
//           stop("Repeated values in timevar");
//         }
//       }
//       for(int p = 0; p != ns; ++p) { 
//         int np = n[p];
//         if(np>0) { // Positive lagged and iterated differences
//           int d1 = diff[0], end = np*d1; 
//           bool L1 = np == 1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           if(end >= l) stop("n * diff needs to be < length(x)");
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(L1) colnam[pos++] = ".G" + diffc[0]; 
//           else colnam[pos++] = ".L" + nc[p] + "G" + diffc[0]; 
//           for(int i = np; i != l; ++i) outp[omap[i]] = growth(x[omap[i]], x[omap[i - np]]);
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1)-1; 
//             for(int i = l-1; i != start; --i) outp[omap[i]] = growth(outp[omap[i]], outp[omap[i - np]]); 
//           }
//           for(int i = end; i--; ) outp[omap[i]] = fill; 
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1], end = np*dq;
//               if(end >= l) stop("n * diff needs to be < length(x)");
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1)-1; 
//                 for(int i = l-1; i != start; --i) outtemp[omap[i]] = growth(outtemp[omap[i]], outtemp[omap[i - np]]); 
//               }
//               for(int i = np*L_dq; i != end; ++i) outtemp[omap[i]] = fill; 
//               out( _ , pos) = outtemp; 
//               if(L1) colnam[pos++] = ".G" + diffc[q]; 
//               else colnam[pos++] = ".L" + nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else if(np<0) { // (Negative) leaded and iterated differences
//           int d1 = diff[0], end = l+np*d1; 
//           bool F1 = np == -1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(F1) colnam[pos++] = ".FG" + diffc[0]; 
//           else colnam[pos++] = ".F" + nc[p] + "G" + diffc[0]; 
//           for(int i = l+np; i--; ) outp[omap[i]] = growth(x[omap[i]], x[omap[i - np]]);
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int final = l+np*(k+1); 
//             for(int i = 0; i != final; ++i) outp[omap[i]] = growth(outp[omap[i]], outp[omap[i - np]]); 
//           }
//           for(int i = end; i != l; ++i) outp[omap[i]] = fill; 
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
//               if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int final = l+np*(k+1); 
//                 for(int i = 0; i != final; ++i) outtemp[omap[i]] = growth(outtemp[omap[i]], outtemp[omap[i - np]]); 
//               }
//               for(int i = end; i != start; ++i) outtemp[omap[i]] = fill; 
//               out( _ , pos) = outtemp; 
//               if(F1) colnam[pos++] = ".FG" + diffc[q]; 
//               else colnam[pos++] = ".F"+ nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else {
//           out( _ , pos) = x;
//           colnam[pos++] = ".--";
//         }
//       }
//     }
//   } else {
//     if(l != g.size()) stop("length(x) must match length(g)");
//     int ags = l/ng, ngp = ng+1, maxdiff = max(diff);  
//     IntegerVector gsv = NULL; 
//     if(Rf_isNull(t)) { 
//       if(maxdiff != 1) { 
//         if(Rf_isNull(gs)) { 
//           gsv = IntegerVector(ng);
//           for(int i = 0; i != l; ++i) ++gsv[g[i]-1];
//         } else {
//           gsv = gs;
//           if(ng != gsv.size()) stop("ng must match length(gs)"); 
//         }
//       }
//       int seen[ngp], memsize = sizeof(int)*(ngp); 
//       for(int p = 0; p != ns; ++p) {
//         int np = n[p];
//         if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
//         if(np>0) { // Positive lagged and iterated differences
//           int d1 = diff[0]; 
//           bool L1 = np == 1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           NumericMatrix::Column outp = out( _ , pos); 
//           memset(seen, 0, memsize); 
//           if(L1) colnam[pos++] = ".G" + diffc[0]; 
//           else colnam[pos++] = ".L" + nc[p] + "G" + diffc[0]; 
//           for(int i = 0; i != l; ++i) { 
//             if(seen[g[i]] == np) outp[i] = growth(x[i], x[i - np]);
//             else {
//               outp[i] = fill;
//               ++seen[g[i]];
//             }
//           }
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1); 
//             memset(seen, 0, memsize); 
//             for(int i = l; i--; ) { 
//               if(seen[g[i]] == gsv[g[i]-1]-start) outp[i] = fill;
//               else {
//                 outp[i] = growth(outp[i], outp[i - np]);
//                 ++seen[g[i]];
//               }
//             }
//           }
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1];
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1); // Right ?? -> seems so!!
//                 memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
//                 for(int i = l; i--; ) {
//                   if(seen[g[i]] == gsv[g[i]-1]-start) outtemp[i] = fill;
//                   else {
//                     outtemp[i] = growth(outtemp[i], outtemp[i - np]);
//                     ++seen[g[i]];
//                   }
//                 }
//               }
//               out( _ , pos) = outtemp; 
//               if(L1) colnam[pos++] = ".G" + diffc[q]; 
//               else colnam[pos++] = ".L" + nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else if(np<0) { // (Negative) leaded and iterated differences
//           int d1 = diff[0]; 
//           bool F1 = np == -1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           NumericMatrix::Column outp = out( _ , pos); 
//           memset(seen, 0, memsize);
//           if(F1) colnam[pos++] = ".FG" + diffc[0]; 
//           else colnam[pos++] = ".F" + nc[p] + "G" + diffc[0]; 
//           for(int i = l; i--; ) {  
//             if(seen[g[i]] == np) outp[i] = growth(x[i], x[i - np]);
//             else {
//               outp[i] = fill;
//               --seen[g[i]];
//             }
//           }
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1); 
//             memset(seen, 0, memsize); 
//             for(int i = 0; i != l; ++i) {
//               if(seen[g[i]] == gsv[g[i]-1]+start) outp[i] = fill; 
//               else {
//                 outp[i] = growth(outp[i], outp[i - np]);
//                 ++seen[g[i]];
//               }
//             }
//           }
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1];
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1); 
//                 memset(seen, 0, memsize); 
//                 for(int i = 0; i != l; ++i) {
//                   if(seen[g[i]] == gsv[g[i]-1]+start) outtemp[i] = fill; 
//                   else {
//                     outtemp[i] = growth(outtemp[i], outtemp[i - np]);
//                     ++seen[g[i]];
//                   }
//                 }
//               }
//               out( _ , pos) = outtemp; 
//               if(F1) colnam[pos++] = ".FG" + diffc[q]; 
//               else colnam[pos++] = ".F"+ nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else {
//           out( _ , pos) = x;
//           colnam[pos++] = ".--";
//         }
//       }
//     } else { // Unordered data: Timevar Provided
//       IntegerVector ord = t; 
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       IntegerVector min(ngp, INT_MAX); 
//       IntegerVector ord2 = no_init_vector(l); 
//       if(Rf_isNull(gs)) { 
//         gsv = IntegerVector(ng);
//         for(int i = 0; i != l; ++i) { 
//           ++gsv[g[i]-1];
//           if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)");
//         for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//       }
//       IntegerVector omap(l); 
//       int cgs[ngp], seen[ngp], memsize = sizeof(int)*(ngp); 
//       cgs[1] = 0;
//       for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; 
//       for(int i = 0; i != l; ++i) {
//         ord2[i] = ord[i] - min[g[i]]; 
//         if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//         if(omap[cgs[g[i]]+ord2[i]] == 0) omap[cgs[g[i]]+ord2[i]] = i; 
//         else stop("Repeated values of timevar within one or more groups"); 
//       }
//       for(int p = 0; p != ns; ++p) {
//         int np = n[p];
//         if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
//         if(np>0) { // Positive lagged and iterated differences
//           int d1 = diff[0]; 
//           bool L1 = np == 1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(L1) colnam[pos++] = ".G" + diffc[0]; 
//           else colnam[pos++] = ".L" + nc[p] + "G" + diffc[0]; 
//           for(int i = 0; i != l; ++i) { 
//             if(ord2[i] >= np) {
//               outp[i] = growth(x[i], x[omap[cgs[g[i]]+ord2[i]-np]]);
//             } else {
//               outp[i] = fill;
//             }
//           }
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1); 
//             memset(seen, 0, memsize); 
//             for(int i = l; i--; ) { 
//               if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outp[omap[i]] = fill;
//               else {
//                 outp[omap[i]] = growth(outp[omap[i]], outp[omap[i - np]]);
//                 ++seen[g[omap[i]]];
//               }
//             }
//           }
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1];
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1); 
//                 memset(seen, 0, memsize); 
//                 for(int i = l; i--; ) {
//                   if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outtemp[omap[i]] = fill;
//                   else {
//                     outtemp[omap[i]] = growth(outtemp[omap[i]], outtemp[omap[i - np]]);
//                     ++seen[g[omap[i]]];
//                   }
//                 }
//               }
//               out( _ , pos) = outtemp; 
//               if(L1) colnam[pos++] = ".G" + diffc[q]; 
//               else colnam[pos++] = ".L" + nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else if(np<0) { // (Negative) leaded and iterated differences
//           int d1 = diff[0]; 
//           bool F1 = np == -1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(F1) colnam[pos++] = ".FG" + diffc[0]; 
//           else colnam[pos++] = ".F" + nc[p] + "G" + diffc[0]; 
//           for(int i = 0; i != l; ++i) { 
//             if(ord2[i] < gsv[g[i]-1]+np) {
//               outp[i] = growth(x[i], x[omap[cgs[g[i]]+ord2[i]-np]]); 
//             } else {
//               outp[i] = fill;
//             }
//           }
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1); 
//             memset(seen, 0, memsize); 
//             for(int i = 0; i != l; ++i) {
//               if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outp[omap[i]] = fill; 
//               else {
//                 outp[omap[i]] = growth(outp[omap[i]], outp[omap[i - np]]);
//                 ++seen[g[omap[i]]];
//               }
//             }
//           }
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1];
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1); 
//                 memset(seen, 0, memsize); 
//                 for(int i = 0; i != l; ++i) {
//                   if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outtemp[omap[i]] = fill; 
//                   else {
//                     outtemp[omap[i]] = growth(outtemp[omap[i]], outtemp[omap[i - np]]);
//                     ++seen[g[omap[i]]];
//                   }
//                 }
//               }
//               out( _ , pos) = outtemp; 
//               if(F1) colnam[pos++] = ".FG" + diffc[q]; 
//               else colnam[pos++] = ".F"+ nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else {
//           out( _ , pos) = x;
//           colnam[pos++] = ".--";
//         }
//       }
//     }
//   }
//   if(ncol == 1) DUPLICATE_ATTRIB(out, x);
//   else out.attr("dimnames") = List::create(x.attr("names"), colnam);
//   return out;
// }




// First Version: Only Growth, manually coded -> Note faster than function substitution !!
// // [[Rcpp::export]]
// NumericVector fgrowthCpp(NumericVector x, IntegerVector n = 1, IntegerVector diff = 1, double fill = NA_REAL, 
//                          int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue, bool logdiff = false) { 
//   
//   int l = x.size(), ns = n.size(), ds = diff.size(), zeros = 0, pos = 0;
//   IntegerVector absn = no_init_vector(ns); 
//   for(int i = ns; i--; ) {
//     if(n[i] == 0) ++zeros;
//     if(n[i] < 0) absn[i] = -n[i];
//     else absn[i] = n[i];
//   }
//   int ncol = (ns-zeros)*ds+zeros;
//   NumericMatrix out = no_init_matrix(l, ncol); 
//   CharacterVector colnam = no_init_vector(ncol); 
//   CharacterVector nc = Rf_coerceVector(absn, STRSXP); 
//   CharacterVector diffc = Rf_coerceVector(diff, STRSXP); 
//   if(ng == 0) { // No groups 
//     if(Rf_isNull(t)) { // Ordered data 
//       for(int p = 0; p != ns; ++p) {
//         int np = n[p];
//         if(np>0) { // Positive lagged and iterated differences
//           int d1 = diff[0], end = np*d1; 
//           bool L1 = np == 1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           if(end >= l) stop("n * diff needs to be < length(x)");
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(L1) colnam[pos++] = ".G" + diffc[0]; 
//           else colnam[pos++] = ".L" + nc[p] + "G" + diffc[0]; 
//           for(int i = np; i != l; ++i) outp[i] = (x[i]-x[i-np])*(100/x[i-np]);
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1)-1; 
//             for(int i = l-1; i != start; --i) outp[i] = (outp[i]-outp[i-np])*(100/outp[i-np]); 
//           }
//           for(int i = end; i--; ) outp[i] = fill; 
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1], end = np*dq;
//               if(end >= l) stop("n * diff needs to be < length(x)");
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1)-1; 
//                 for(int i = l-1; i != start; --i) outtemp[i] = (outtemp[i]-outtemp[i-np])*(100/outtemp[i-np]); 
//               }
//               for(int i = np*L_dq; i != end; ++i) outtemp[i] = fill; 
//               out( _ , pos) = outtemp; 
//               if(L1) colnam[pos++] = ".G" + diffc[q]; 
//               else colnam[pos++] = ".L" + nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else if(np<0) { // (Negative) leaded and iterated differences
//           int d1 = diff[0], end = l+np*d1; 
//           bool F1 = np == -1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(F1) colnam[pos++] = ".FG" + diffc[0]; 
//           else colnam[pos++] = ".F" + nc[p] + "G" + diffc[0];  
//           for(int i = l+np; i--; ) outp[i] = (x[i]-x[i-np])*(100/x[i-np]);
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int final = l+np*(k+1); 
//             for(int i = 0; i != final; ++i) outp[i] = (outp[i]-outp[i-np])*(100/outp[i-np]); 
//           }
//           for(int i = end; i != l; ++i) outp[i] = fill; 
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
//               if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int final = l+np*(k+1); 
//                 for(int i = 0; i != final; ++i) outtemp[i] = (outtemp[i]-outtemp[i-np])*(100/outtemp[i-np]); 
//               }
//               for(int i = end; i != start; ++i) outtemp[i] = fill; 
//               out( _ , pos) = outtemp; 
//               if(F1) colnam[pos++] = ".FG" + diffc[q]; 
//               else colnam[pos++] = ".F"+ nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else {
//           out( _ , pos) = x;
//           colnam[pos++] = ".--";
//         }
//       }
//     } else { // Unordered data: Timevar provided 
//       IntegerVector ord = t;
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       LogicalVector ocheck(l, true); 
//       int omap[l]; 
//       for(int i = 0; i != l; ++i) { 
//         if(ord[i] > l) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
//         if(ocheck[ord[i]-1]) {
//           ocheck[ord[i]-1] = false;
//           omap[ord[i]-1] = i; 
//         } else {
//           stop("Repeated values in timevar");
//         }
//       }
//       for(int p = 0; p != ns; ++p) { 
//         int np = n[p];
//         if(np>0) { // Positive lagged and iterated differences
//           int d1 = diff[0], end = np*d1; 
//           bool L1 = np == 1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           if(end >= l) stop("n * diff needs to be < length(x)");
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(L1) colnam[pos++] = ".G" + diffc[0]; 
//           else colnam[pos++] = ".L" + nc[p] + "G" + diffc[0]; 
//           for(int i = np; i != l; ++i) outp[omap[i]] = (x[omap[i]]-x[omap[i-np]])*(100/x[omap[i-np]]);
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1)-1; 
//             for(int i = l-1; i != start; --i) outp[omap[i]] = (outp[omap[i]]-outp[omap[i-np]])*(100/outp[omap[i-np]]); 
//           }
//           for(int i = end; i--; ) outp[omap[i]] = fill; 
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1], end = np*dq;
//               if(end >= l) stop("n * diff needs to be < length(x)");
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1)-1; 
//                 for(int i = l-1; i != start; --i) outtemp[omap[i]] = (outtemp[omap[i]]-outtemp[omap[i-np]])*(100/outtemp[omap[i-np]]); 
//               }
//               for(int i = np*L_dq; i != end; ++i) outtemp[omap[i]] = fill; 
//               out( _ , pos) = outtemp; 
//               if(L1) colnam[pos++] = ".G" + diffc[q]; 
//               else colnam[pos++] = ".L" + nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else if(np<0) { // (Negative) leaded and iterated differences
//           int d1 = diff[0], end = l+np*d1; 
//           bool F1 = np == -1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(F1) colnam[pos++] = ".FG" + diffc[0]; 
//           else colnam[pos++] = ".F" + nc[p] + "G" + diffc[0]; 
//           for(int i = l+np; i--; ) outp[omap[i]] = (x[omap[i]]-x[omap[i-np]])*(100/x[omap[i-np]]);
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int final = l+np*(k+1); 
//             for(int i = 0; i != final; ++i) outp[omap[i]] = (outp[omap[i]]-outp[omap[i-np]])*(100/outp[omap[i-np]]); 
//           }
//           for(int i = end; i != l; ++i) outp[omap[i]] = fill; 
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
//               if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int final = l+np*(k+1); 
//                 for(int i = 0; i != final; ++i) outtemp[omap[i]] = (outtemp[omap[i]]-outtemp[omap[i-np]])*(100/outtemp[omap[i-np]]); 
//               }
//               for(int i = end; i != start; ++i) outtemp[omap[i]] = fill; 
//               out( _ , pos) = outtemp; 
//               if(F1) colnam[pos++] = ".FG" + diffc[q]; 
//               else colnam[pos++] = ".F"+ nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else {
//           out( _ , pos) = x;
//           colnam[pos++] = ".--";
//         }
//       }
//     }
//   } else {
//     if(l != g.size()) stop("length(x) must match length(g)");
//     int ags = l/ng, ngp = ng+1, maxdiff = max(diff);  
//     IntegerVector gsv = NULL; 
//     if(Rf_isNull(t)) { 
//       if(maxdiff != 1) { 
//         if(Rf_isNull(gs)) { 
//           gsv = IntegerVector(ng);
//           for(int i = 0; i != l; ++i) ++gsv[g[i]-1];
//         } else {
//           gsv = gs;
//           if(ng != gsv.size()) stop("ng must match length(gs)"); 
//         }
//       }
//       int seen[ngp], memsize = sizeof(int)*(ngp); 
//       for(int p = 0; p != ns; ++p) {
//         int np = n[p];
//         if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
//         if(np>0) { // Positive lagged and iterated differences
//           int d1 = diff[0]; 
//           bool L1 = np == 1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           NumericMatrix::Column outp = out( _ , pos); 
//           memset(seen, 0, memsize); 
//           if(L1) colnam[pos++] = ".G" + diffc[0]; 
//           else colnam[pos++] = ".L" + nc[p] + "G" + diffc[0]; 
//           for(int i = 0; i != l; ++i) { 
//             if(seen[g[i]] == np) outp[i] = (x[i]-x[i-np])*(100/x[i-np]);
//             else {
//               outp[i] = fill;
//               ++seen[g[i]];
//             }
//           }
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1); 
//             memset(seen, 0, memsize); 
//             for(int i = l; i--; ) { 
//               if(seen[g[i]] == gsv[g[i]-1]-start) outp[i] = fill;
//               else {
//                 outp[i] = (outp[i]-outp[i-np])*(100/outp[i-np]);
//                 ++seen[g[i]];
//               }
//             }
//           }
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1];
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1); // Right ?? -> seems so!!
//                 memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
//                 for(int i = l; i--; ) {
//                   if(seen[g[i]] == gsv[g[i]-1]-start) outtemp[i] = fill;
//                   else {
//                     outtemp[i] = (outtemp[i]-outtemp[i-np])*(100/outtemp[i-np]);
//                     ++seen[g[i]];
//                   }
//                 }
//               }
//               out( _ , pos) = outtemp; 
//               if(L1) colnam[pos++] = ".G" + diffc[q]; 
//               else colnam[pos++] = ".L" + nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else if(np<0) { // (Negative) leaded and iterated differences
//           int d1 = diff[0]; 
//           bool F1 = np == -1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           NumericMatrix::Column outp = out( _ , pos); 
//           memset(seen, 0, memsize);
//           if(F1) colnam[pos++] = ".FG" + diffc[0]; 
//           else colnam[pos++] = ".F" + nc[p] + "G" + diffc[0]; 
//           for(int i = l; i--; ) {  
//             if(seen[g[i]] == np) outp[i] = (x[i]-x[i-np])*(100/x[i-np]);
//             else {
//               outp[i] = fill;
//               --seen[g[i]];
//             }
//           }
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1); 
//             memset(seen, 0, memsize); 
//             for(int i = 0; i != l; ++i) {
//               if(seen[g[i]] == gsv[g[i]-1]+start) outp[i] = fill; 
//               else {
//                 outp[i] = (outp[i]-outp[i-np])*(100/outp[i-np]);
//                 ++seen[g[i]];
//               }
//             }
//           }
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1];
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1); 
//                 memset(seen, 0, memsize); 
//                 for(int i = 0; i != l; ++i) {
//                   if(seen[g[i]] == gsv[g[i]-1]+start) outtemp[i] = fill; 
//                   else {
//                     outtemp[i] = (outtemp[i]-outtemp[i-np])*(100/outtemp[i-np]);
//                     ++seen[g[i]];
//                   }
//                 }
//               }
//               out( _ , pos) = outtemp; 
//               if(F1) colnam[pos++] = ".FG" + diffc[q]; 
//               else colnam[pos++] = ".F"+ nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else {
//           out( _ , pos) = x;
//           colnam[pos++] = ".--";
//         }
//       }
//     } else { // Unordered data: Timevar Provided
//       IntegerVector ord = t; 
//       if(l != ord.size()) stop("length(x) must match length(t)");
//       IntegerVector min(ngp, INT_MAX); 
//       IntegerVector ord2 = no_init_vector(l); 
//       if(Rf_isNull(gs)) { 
//         gsv = IntegerVector(ng);
//         for(int i = 0; i != l; ++i) { 
//           ++gsv[g[i]-1];
//           if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//         }
//       } else {
//         gsv = gs;
//         if(ng != gsv.size()) stop("ng must match length(gs)");
//         for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//       }
//       IntegerVector omap(l); 
//       int cgs[ngp], seen[ngp], memsize = sizeof(int)*(ngp); 
//       cgs[1] = 0;
//       for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; 
//       for(int i = 0; i != l; ++i) {
//         ord2[i] = ord[i] - min[g[i]]; 
//         if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//         if(omap[cgs[g[i]]+ord2[i]] == 0) omap[cgs[g[i]]+ord2[i]] = i; 
//         else stop("Repeated values of timevar within one or more groups"); 
//       }
//       for(int p = 0; p != ns; ++p) {
//         int np = n[p];
//         if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
//         if(np>0) { // Positive lagged and iterated differences
//           int d1 = diff[0]; 
//           bool L1 = np == 1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(L1) colnam[pos++] = ".G" + diffc[0]; 
//           else colnam[pos++] = ".L" + nc[p] + "G" + diffc[0]; 
//           for(int i = 0; i != l; ++i) { 
//             if(ord2[i] >= np) {
//               outp[i] = (x[i]-x[omap[cgs[g[i]]+ord2[i]-np]])*(100/x[omap[cgs[g[i]]+ord2[i]-np]]);
//             } else {
//               outp[i] = fill;
//             }
//           }
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1); 
//             memset(seen, 0, memsize); 
//             for(int i = l; i--; ) { 
//               if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outp[omap[i]] = fill;
//               else {
//                 outp[omap[i]] = (outp[omap[i]]-outp[omap[i-np]])*(100/outp[omap[i-np]]);
//                 ++seen[g[omap[i]]];
//               }
//             }
//           }
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1];
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1); 
//                 memset(seen, 0, memsize); 
//                 for(int i = l; i--; ) {
//                   if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outtemp[omap[i]] = fill;
//                   else {
//                     outtemp[omap[i]] = (outtemp[omap[i]]-outtemp[omap[i-np]])*(100/outtemp[omap[i-np]]);
//                     ++seen[g[omap[i]]];
//                   }
//                 }
//               }
//               out( _ , pos) = outtemp; 
//               if(L1) colnam[pos++] = ".G" + diffc[q]; 
//               else colnam[pos++] = ".L" + nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else if(np<0) { // (Negative) leaded and iterated differences
//           int d1 = diff[0]; 
//           bool F1 = np == -1;
//           if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//           NumericMatrix::Column outp = out( _ , pos); 
//           if(F1) colnam[pos++] = ".FG" + diffc[0]; 
//           else colnam[pos++] = ".F" + nc[p] + "G" + diffc[0]; 
//           for(int i = 0; i != l; ++i) { 
//             if(ord2[i] < gsv[g[i]-1]+np) {
//               outp[i] = (x[i]-x[omap[cgs[g[i]]+ord2[i]-np]])*(100/x[omap[cgs[g[i]]+ord2[i]-np]]); 
//             } else {
//               outp[i] = fill;
//             }
//           }
//           if(d1 > 1) for(int k = 1; k != d1; ++k) {
//             int start = np*(k+1); 
//             memset(seen, 0, memsize); 
//             for(int i = 0; i != l; ++i) {
//               if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outp[omap[i]] = fill; 
//               else {
//                 outp[omap[i]] = (outp[omap[i]]-outp[omap[i-np]])*(100/outp[omap[i-np]]);
//                 ++seen[g[omap[i]]];
//               }
//             }
//           }
//           if(ds > 1) {
//             NumericVector outtemp = outp; 
//             for(int q = 1; q != ds; ++q) {
//               int dq = diff[q], L_dq = diff[q-1];
//               if(dq <= L_dq) stop("differences must be passed in ascending order");
//               for(int k = L_dq; k != dq; ++k) {
//                 int start = np*(k+1); 
//                 memset(seen, 0, memsize); 
//                 for(int i = 0; i != l; ++i) {
//                   if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outtemp[omap[i]] = fill; 
//                   else {
//                     outtemp[omap[i]] = (outtemp[omap[i]]-outtemp[omap[i-np]])*(100/outtemp[omap[i-np]]);
//                     ++seen[g[omap[i]]];
//                   }
//                 }
//               }
//               out( _ , pos) = outtemp; 
//               if(F1) colnam[pos++] = ".FG" + diffc[q]; 
//               else colnam[pos++] = ".F"+ nc[p] + "G" + diffc[q]; 
//             }
//           }
//         } else {
//           out( _ , pos) = x;
//           colnam[pos++] = ".--";
//         }
//       }
//     }
//   }
//   if(ncol == 1) DUPLICATE_ATTRIB(out, x);
//   else out.attr("dimnames") = List::create(x.attr("names"), colnam);
//   return out;
// }

