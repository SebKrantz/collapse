// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
// #include <Rmath.h>
// #include "C:\Users\Sebastian Krantz\Documents\R\win-library\3.4\fmt\include\fmt\format.h" // doesn't work !!
using namespace Rcpp;

// Todo: Add names argument ?? 

// Final version with comments, for trials see example codes below !!
// Note: NumericVector is faster than SEXP !!
// [[Rcpp::export]]
NumericVector fdiffCpp(const NumericVector& x, const IntegerVector& n = 1, const IntegerVector& diff = 1, 
                       double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0, 
                       const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue, bool names = true) { 
  
  int l = x.size(), ns = n.size(), ds = diff.size(), zeros = 0, pos = 0;
  IntegerVector absn = no_init_vector(ns); // better way ?? 
  for(int i = ns; i--; ) {
    if(n[i] == 0) ++zeros;
    if(n[i] < 0) absn[i] = -n[i];
    else absn[i] = n[i];
  }
  int ncol = (ns-zeros)*ds+zeros;
  if(ncol == 1) names = false;
  NumericMatrix out = no_init_matrix(l, ncol); // fastest ??
  CharacterVector colnam = names ? no_init_vector(ncol) : no_init_vector(1); 
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING; //  as<CharacterVector>(absn); // fastest ?? 
  CharacterVector diffc = names ? Rf_coerceVector(diff, STRSXP) : NA_STRING; // as<CharacterVector>(diff);
  if(ng == 0) { // No groups 
    if(Rf_isNull(t)) { // Ordered data -> Perfect !!
      for(int p = 0; p != ns; ++p) {
        int np = n[p];
        if(np>0) { // Positive lagged and iterated differences
          int d1 = diff[0], end = np*d1; 
          bool L1 = np == 1;
          if(d1 < 1) stop("diff must be a vector of integers > 0"); 
          if(end >= l) stop("n * diff needs to be < length(x)");
          NumericMatrix::Column outp = out( _ , pos); // p*ds this gave D1L1, D2L1, D1L2, D2L2, but I want D1L1, D1L2, D2L1, D2L2
          if(names) {
            if(L1) colnam[pos] = ".D" + diffc[0]; 
            else colnam[pos] = ".L" + nc[p] + "D" + diffc[0]; //fmt::format("{}{}", "D.", np);
          }
          ++pos;
          for(int i = np; i != l; ++i) outp[i] = x[i] - x[i - np];
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1)-1; //  This is right !!
            for(int i = l-1; i != start; --i) outp[i] -= outp[i - np]; // needed for correct iteration 
          }
          for(int i = end; i--; ) outp[i] = fill; 
          if(ds > 1) {
            NumericVector outtemp = outp; 
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1], end = np*dq;
              if(end >= l) stop("n * diff needs to be < length(x)");
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int start = np*(k+1)-1; // good ?? -> seems good !!
                for(int i = l-1; i != start; --i) outtemp[i] -= outtemp[i - np]; // needed for correct iteration
              }
              for(int i = np*L_dq; i != end; ++i) outtemp[i] = fill; // good ?? -> yes !!
              out( _ , pos) = outtemp; // p*ds+q this gave D1L1, D2L1, D1L2, D2L2, but I want D1L1, D1L2, D2L1, D2L2
              if(names) {
                if(L1) colnam[pos] = ".D" + diffc[q]; // fmt::format("{}{}", "D", q+1); 
                else colnam[pos] = ".L" + nc[p] + "D" + diffc[q]; //fmt::format("{}{}{}{}", "D", q+1, ".", np); 
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
            if(F1) colnam[pos] = ".FD" + diffc[0]; 
            else colnam[pos] = ".F" + nc[p] + "D" + diffc[0]; // fmt::format("{}{}", "FD.", -np); 
          }
          ++pos;
          for(int i = l+np; i--; ) outp[i] = x[i] - x[i - np];
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int final = l+np*(k+1); // good ?? -> Yes !!
            for(int i = 0; i != final; ++i) outp[i] -= outp[i - np]; // gives correct iteration ?? -> yes !!
          }
          for(int i = end; i != l; ++i) outp[i] = fill; // good ?? -> yes !!
          if(ds > 1) {
            NumericVector outtemp = outp; 
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
              if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int final = l+np*(k+1); // good ?? -> yes !!
                for(int i = 0; i != final; ++i) outtemp[i] -= outtemp[i - np]; // gives correct iteration ?? -> yes !!
              }
              for(int i = end; i != start; ++i) outtemp[i] = fill; // good ?? -> yes !!
              out( _ , pos) = outtemp; 
              if(names) {
                if(F1) colnam[pos] = ".FD" + diffc[q]; //  fmt::format("{}{}", "FD", q+1); 
                else colnam[pos] = ".F"+ nc[p] + "D" + diffc[q]; // fmt::format("{}{}{}{}", "FD", q+1, ".", -np);
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
    } else { // Unordered data: Timevar provided -> works ??? ALL CORRECT AND EFFICIENT ??? -> Seems to have introduced some instability
      IntegerVector ord = t;
      if(l != ord.size()) stop("length(x) must match length(t)");
      LogicalVector ocheck(l, true); // Could do more efficiently if only one value -> but you don't use this anyway!!
      int omap[l]; // otherwise always ord[i]-1 -> efficiency gain !! -> check also for panel-lag !!!!!! // Rcpp vector faster ??
      for(int i = 0; i != l; ++i) { // integrated below !!
        if(ord[i] > l) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i; // Note: omap is the same as order(ord) !!
        } else {
          stop("Repeated values in timevar");
        }
      }
      for(int p = 0; p != ns; ++p) { // This is the same code as above, just doing everything through omap !!
        int np = n[p];
        if(np>0) { // Positive lagged and iterated differences
          int d1 = diff[0], end = np*d1; 
          bool L1 = np == 1;
          if(d1 < 1) stop("diff must be a vector of integers > 0"); 
          if(end >= l) stop("n * diff needs to be < length(x)");
          NumericMatrix::Column outp = out( _ , pos); 
          if(names) {
            if(L1) colnam[pos] = ".D" + diffc[0]; 
            else colnam[pos] = ".L" + nc[p] + "D" + diffc[0]; 
          }
          ++pos;
          for(int i = np; i != l; ++i) outp[omap[i]] = x[omap[i]] - x[omap[i - np]];
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1)-1; 
            for(int i = l-1; i != start; --i) outp[omap[i]] -= outp[omap[i - np]]; 
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
                for(int i = l-1; i != start; --i) outtemp[omap[i]] -= outtemp[omap[i - np]]; 
              }
              for(int i = np*L_dq; i != end; ++i) outtemp[omap[i]] = fill; 
              out( _ , pos) = outtemp; 
              if(names) {
                if(L1) colnam[pos] = ".D" + diffc[q]; 
                else colnam[pos] = ".L" + nc[p] + "D" + diffc[q]; 
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
            if(F1) colnam[pos] = ".FD" + diffc[0]; 
            else colnam[pos] = ".F" + nc[p] + "D" + diffc[0]; 
          }
          ++pos;
          for(int i = l+np; i--; ) outp[omap[i]] = x[omap[i]] - x[omap[i - np]];
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int final = l+np*(k+1); 
            for(int i = 0; i != final; ++i) outp[omap[i]] -= outp[omap[i - np]]; 
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
                for(int i = 0; i != final; ++i) outtemp[omap[i]] -= outtemp[omap[i - np]]; 
              }
              for(int i = end; i != start; ++i) outtemp[omap[i]] = fill; 
              out( _ , pos) = outtemp; 
              if(names) {
                if(F1) colnam[pos] = ".FD" + diffc[q]; 
                else colnam[pos] = ".F"+ nc[p] + "D" + diffc[q]; 
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
    int ags = l/ng, ngp = ng+1, maxdiff = max(diff); //  diff[ds-1]; // maybe wrong order ?? 
    // IntegerVector gsv = (Rf_isNull(gs)) ? IntegerVector(ng) : gs; //  no_init_vector(ng); // stable ?? -> Nope !!
    // http://www.cplusplus.com/articles/1AUq5Di1/ // using the ternary operator
    IntegerVector gsv = NULL; // https://stackoverflow.com/questions/1793807/declaring-a-variable-in-an-if-else-block-in-c
    if(Rf_isNull(t)) { // Ordered data -> Solution seems to work well !!!
      if(maxdiff != 1) { // The gsv solution above faster than std::fill ????????????????????
        if(Rf_isNull(gs)) { // needed for proper iteration 
          gsv = IntegerVector(ng);
          for(int i = 0; i != l; ++i) ++gsv[g[i]-1];
        } else {
          gsv = gs;
          if(ng != gsv.size()) stop("ng must match length(gs)"); 
        }
      }
      int seen[ngp], memsize = sizeof(int)*(ngp); // adding 1 here guards agaist subtracting 1 everywhere else !! -> check for lag also !!
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
            if(L1) colnam[pos] = ".D" + diffc[0]; 
            else colnam[pos] = ".L" + nc[p] + "D" + diffc[0]; 
          }
          ++pos;
          for(int i = 0; i != l; ++i) { // this loop ??
            if(seen[g[i]] == np) outp[i] = x[i] - x[i - np];
            else {
              outp[i] = fill;
              ++seen[g[i]];
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1); // right ?? -> seems so!! 
            memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
            for(int i = l; i--; ) { // correct iteration ?? -> Yes !!
              if(seen[g[i]] == gsv[g[i]-1]-start) outp[i] = fill;
              else {
                outp[i] -= outp[i - np];
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
                    outtemp[i] -= outtemp[i - np];
                    ++seen[g[i]];
                  }
                }
              }
              out( _ , pos) = outtemp; 
              if(names) {
                if(L1) colnam[pos] = ".D" + diffc[q]; 
                else colnam[pos] = ".L" + nc[p] + "D" + diffc[q]; 
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
            if(F1) colnam[pos] = ".FD" + diffc[0]; 
            else colnam[pos] = ".F" + nc[p] + "D" + diffc[0]; 
          }
          ++pos;
          for(int i = l; i--; ) { // good?? -> yes !! 
            if(seen[g[i]] == np) outp[i] = x[i] - x[i - np];
            else {
              outp[i] = fill;
              --seen[g[i]];
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1); // good ?? -> seems right !!
            memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
            for(int i = 0; i != l; ++i) {
              if(seen[g[i]] == gsv[g[i]-1]+start) outp[i] = fill; 
              else {
                outp[i] -= outp[i - np];
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
                int start = np*(k+1); // Right ?? -> seems right !!
                memset(seen, 0, memsize); // missing a problem ??, why does it work without ??????????????????????????ß
                for(int i = 0; i != l; ++i) {
                  if(seen[g[i]] == gsv[g[i]-1]+start) outtemp[i] = fill; 
                  else {
                    outtemp[i] -= outtemp[i - np];
                    ++seen[g[i]];
                  }
                }
              }
              out( _ , pos) = outtemp; 
              if(names) {
                if(F1) colnam[pos] = ".FD" + diffc[q]; 
                else colnam[pos] = ".F"+ nc[p] + "D" + diffc[q]; 
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
      // Optimal implementation for massively unordered data, formerly fdiffCpppuo3 !!
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
      int cgs[ngp], seen[ngp], memsize = sizeof(int)*(ngp); // Note: -> attempts to otimize memory use of seen[ngs] are futile !!
      cgs[1] = 0;
      for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; // or get "starts from forderv"
      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]]; 
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        if(omap[cgs[g[i]]+ord2[i]] == 0) omap[cgs[g[i]]+ord2[i]] = i; // index ??
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
            if(L1) colnam[pos] = ".D" + diffc[0]; 
            else colnam[pos] = ".L" + nc[p] + "D" + diffc[0]; 
          }
          ++pos;
          for(int i = 0; i != l; ++i) { // ordinary computation on first iteration
            if(ord2[i] >= np) {
              outp[i] = x[i] - x[omap[cgs[g[i]]+ord2[i]-np]];
            } else {
              outp[i] = fill;
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1); 
            memset(seen, 0, memsize); 
            for(int i = l; i--; ) { // correct iteration ?? optimize through i2 = omap[i]?????????????
              if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outp[omap[i]] = fill;
              else {
                outp[omap[i]] -= outp[omap[i - np]];
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
                int start = np*(k+1); // Right ?? 
                memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
                for(int i = l; i--; ) {
                  if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outtemp[omap[i]] = fill;
                  else {
                    outtemp[omap[i]] -= outtemp[omap[i - np]];
                    ++seen[g[omap[i]]];
                  }
                }
              }
              out( _ , pos) = outtemp; 
              if(names) {
                if(L1) colnam[pos] = ".D" + diffc[q]; 
                else colnam[pos] = ".L" + nc[p] + "D" + diffc[q]; 
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
            if(F1) colnam[pos] = ".FD" + diffc[0]; 
            else colnam[pos] = ".F" + nc[p] + "D" + diffc[0]; 
          }
          ++pos;
          for(int i = 0; i != l; ++i) { // Ordinary computation on first iteration 
            if(ord2[i] < gsv[g[i]-1]+np) {
              outp[i] = x[i] - x[omap[cgs[g[i]]+ord2[i]-np]]; 
            } else {
              outp[i] = fill;
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1); // good ?? 
            memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
            for(int i = 0; i != l; ++i) {
              if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outp[omap[i]] = fill; // was -start -> error !!
              else {
                outp[omap[i]] -= outp[omap[i - np]];
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
                int start = np*(k+1); // Right ?? 
                memset(seen, 0, memsize); // missing a problem ??, why does it work without ??
                for(int i = 0; i != l; ++i) {
                  if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outtemp[omap[i]] = fill; // was -start -> error !!
                  else {
                    outtemp[omap[i]] -= outtemp[omap[i - np]];
                    ++seen[g[omap[i]]];
                  }
                }
              }
              out( _ , pos) = outtemp; 
              if(names) {
                if(F1) colnam[pos] = ".FD" + diffc[q]; 
                else colnam[pos] = ".F"+ nc[p] + "D" + diffc[q]; 
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
        // Rf_setAttrib(out, "dimnames", x.attr("names"));
  // Previous: 
  // if(ncol == 1) DUPLICATE_ATTRIB(out, x);
  // else if(names) out.attr("dimnames") = List::create(x.attr("names"), colnam);
        // else out.attr("dimnames") = List::create(x.attr("names"), R_NilValue); // Redundat, no names means nothing !!
  
  DUPLICATE_ATTRIB(out, x); 
  if(ncol != 1) { 
    if(x.hasAttribute("names")) out.attr("names") = R_NilValue; 
    out.attr("dim") = Dimension(l, ncol);
    if(names) out.attr("dimnames") = List::create(x.attr("names"), colnam);
  }
  
  return out;
}


// Unordered panel lag other implementation 3: Implemented Above !!
// // [[Rcpp::export]]
// SEXP fdiffCpppuo3(NumericVector x, IntegerVector n = 1, IntegerVector diff = 1, double fill = NA_REAL, 
//                   int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) { 
//   
//   int l = x.size(), ns = n.size(), ds = diff.size(), zeros = 0, pos = 0;
//   IntegerVector absn = no_init_vector(ns); // better way ?? 
//   for(int i = ns; i--; ) {
//     if(n[i] == 0) ++zeros;
//     if(n[i] < 0) absn[i] = -n[i];
//     else absn[i] = n[i];
//   }
//   int ncol = (ns-zeros)*ds+zeros;
//   NumericMatrix out = no_init_matrix(l, ncol); // fastest ??
//   CharacterVector colnam = no_init_vector(ncol); 
//   CharacterVector nc = Rf_coerceVector(absn, STRSXP);  
//   CharacterVector diffc = Rf_coerceVector(diff, STRSXP); 
//   if(l != g.size()) stop("length(x) must match length(g)");
//   int ags = l/ng, maxdiff = diff[ds-1], ngp = ng+1; 
//   
//   IntegerVector ord = t; 
//   if(l != ord.size()) stop("length(x) must match length(t)");
//   IntegerVector min(ngp, INT_MAX); 
//   IntegerVector gsv = NULL; 
//   IntegerVector ord2 = no_init_vector(l); 
//   if(Rf_isNull(gs)) { 
//     gsv = IntegerVector(ng);
//     for(int i = 0; i != l; ++i) { 
//       ++gsv[g[i]-1];
//       if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//     }
//   } else {
//     gsv = gs;
//     if(ng != gsv.size()) stop("ng must match length(gs)");
//     for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//   }
//   IntegerVector omap(l); 
//   int cgs[ngp], seen[ngp], memsize = sizeof(int)*(ngp); // i2 = 0;
//   cgs[1] = 0;
//   for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; // or get "starts from forderv"
//   for(int i = 0; i != l; ++i) {
//     ord2[i] = ord[i] - min[g[i]]; 
//     if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//     if(omap[cgs[g[i]]+ord2[i]] == 0) omap[cgs[g[i]]+ord2[i]] = i; // fastest ??
//     else stop("Repeated values of timevar within one or more groups"); 
//   }
//   // return omap;
//   for(int p = 0; p != ns; ++p) {
//     int np = n[p];
//     if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
//     if(np>0) { // Positive lagged and iterated differences
//       int d1 = diff[0]; 
//       bool L1 = np == 1;
//       if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//       NumericMatrix::Column outp = out( _ , pos); 
//       if(L1) colnam[pos++] = ".D" + diffc[0]; 
//       else colnam[pos++] = ".L" + nc[p] + "D" + diffc[0]; 
//       for(int i = 0; i != l; ++i) { // ordinary computation on first iteration
//         if(ord2[i] >= np) {
//           outp[i] = x[i] - x[omap[cgs[g[i]]+ord2[i]-np]];
//         } else {
//           outp[i] = fill;
//         }
//       }
//       if(d1 > 1) for(int k = 1; k != d1; ++k) {
//         int start = np*(k+1); 
//         memset(seen, 0, memsize); 
//         for(int i = l; i--; ) { // correct iteration ?? 
//           if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outp[omap[i]] = fill;
//           else {
//             outp[omap[i]] -= outp[omap[i - np]];
//             ++seen[g[omap[i]]];
//           }
//         }
//       }
//       if(ds > 1) {
//         NumericVector outtemp = outp; 
//         for(int q = 1; q != ds; ++q) {
//           int dq = diff[q], L_dq = diff[q-1];
//           if(dq <= L_dq) stop("differences must be passed in ascending order");
//           for(int k = L_dq; k != dq; ++k) {
//             int start = np*(k+1); // Right ?? 
//             memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
//             for(int i = l; i--; ) {
//               if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outtemp[omap[i]] = fill;
//               else {
//                 outtemp[omap[i]] -= outtemp[omap[i - np]];
//                 ++seen[g[omap[i]]];
//               }
//             }
//           }
//           out( _ , pos) = outtemp; 
//           if(L1) colnam[pos++] = ".D" + diffc[q]; 
//           else colnam[pos++] = ".L" + nc[p] + "D" + diffc[q]; 
//         }
//       }
//     } else if(np<0) { // (Negative) leaded and iterated differences
//       int d1 = diff[0]; 
//       bool F1 = np == -1;
//       if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//       NumericMatrix::Column outp = out( _ , pos); 
//       if(F1) colnam[pos++] = ".FD" + diffc[0]; 
//       else colnam[pos++] = ".F" + nc[p] + "D" + diffc[0]; 
//       for(int i = 0; i != l; ++i) { // Ordinary computation on first iteration 
//         if(ord2[i] < gsv[g[i]-1]+np) {
//           outp[i] = x[i] - x[omap[cgs[g[i]]+ord2[i]-np]]; 
//         } else {
//           outp[i] = fill;
//         }
//       }
//       if(d1 > 1) for(int k = 1; k != d1; ++k) {
//         int start = np*(k+1); // good ?? 
//         memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
//         for(int i = 0; i != l; ++i) {
//           if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outp[omap[i]] = fill; // was -start -> error !!
//           else {
//             outp[omap[i]] -= outp[omap[i - np]];
//             ++seen[g[omap[i]]];
//           }
//         }
//       }
//       if(ds > 1) {
//         NumericVector outtemp = outp; 
//         for(int q = 1; q != ds; ++q) {
//           int dq = diff[q], L_dq = diff[q-1];
//           if(dq <= L_dq) stop("differences must be passed in ascending order");
//           for(int k = L_dq; k != dq; ++k) {
//             int start = np*(k+1); // Right ?? 
//             memset(seen, 0, memsize); // missing !!!!!!!!!!!!!!!!! a problem ??, why does it work without ??
//             for(int i = 0; i != l; ++i) {
//               if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outtemp[omap[i]] = fill; // was -start -> error !!
//               else {
//                 outtemp[omap[i]] -= outtemp[omap[i - np]];
//                 ++seen[g[omap[i]]];
//               }
//             }
//           }
//           out( _ , pos) = outtemp; 
//           if(F1) colnam[pos++] = ".FD" + diffc[q]; 
//           else colnam[pos++] = ".F"+ nc[p] + "D" + diffc[q]; 
//         }
//       }
//     } else {
//       out( _ , pos) = x;
//       colnam[pos++] = ".--";
//     }
//   }
//   if(ncol == 1) SHALLOW_DUPLICATE_ATTRIB(out, x);
//   else out.attr("dimnames") = List::create(x.attr("names"), colnam);
//   return out;
// }


// Unordered panel lag other implementation 2: Work with 1D array any only access thrugh that -> same as before, need counter: most momory efficient and fastest if relatively ordered !!  -> slowest !!!!!!!!!
// Note: the optimal version implemented above is a bot less memory efficient (needs ord2), but fastest if the data is really unordered !! 
// // [[Rcpp::export]]
// SEXP fdiffCpppuo2(NumericVector x, IntegerVector n = 1, IntegerVector diff = 1, double fill = NA_REAL, 
//                   int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) { 
//   
//   int l = x.size(), ns = n.size(), ds = diff.size(), zeros = 0, pos = 0;
//   IntegerVector absn = no_init_vector(ns); // better way ?? 
//   for(int i = ns; i--; ) {
//     if(n[i] == 0) ++zeros;
//     if(n[i] < 0) absn[i] = -n[i];
//     else absn[i] = n[i];
//   }
//   int ncol = (ns-zeros)*ds+zeros;
//   NumericMatrix out = no_init_matrix(l, ncol); // fastest ??
//   CharacterVector colnam = no_init_vector(ncol); 
//   CharacterVector nc = Rf_coerceVector(absn, STRSXP);  
//   CharacterVector diffc = Rf_coerceVector(diff, STRSXP); 
//   if(l != g.size()) stop("length(x) must match length(g)");
//   int ags = l/ng, maxdiff = diff[ds-1], ngp = ng+1; 
//   
//   IntegerVector ord = t; 
//   if(l != ord.size()) stop("length(x) must match length(t)");
//   IntegerVector min(ngp, INT_MAX); // Best solution ?? array ?? 
//   IntegerVector gsv = NULL; // The gsv solution above faster than std::fill ????????????????????
//   if(Rf_isNull(gs)) { // needed for proper iteration 
//     gsv = IntegerVector(ng);
//     for(int i = 0; i != l; ++i) { 
//       ++gsv[g[i]-1];
//       if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//     }
//   } else {
//     gsv = gs;
//     if(ng != gsv.size()) stop("ng must match length(gs)");
//     for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//   }
//   // int *omap = new int[l]{}; // slower than integer vector !!
//   IntegerVector omap(l); // best ??????????????????????????? -> yes !!!! // also see if using R's order is faster ?? !!-> but then no integrity checks !!!
//   // int omap[l]; //  omaps = sizeof(int)*l; // all these solutions are unstable, for some reason!!
//   // memset(omap, 0, omaps); 
//   // std::fill_n(omap, l, 0);
//   int cgs[ngp], seen[ngp], memsize = sizeof(int)*(ngp), ord2 = 0; // i2 = 0;
//   cgs[1] = 0;
//   for(int i = 2; i != ngp; ++i) cgs[i] = cgs[i-1] + gsv[i-2]; // or get "starts from forderv"
//   for(int i = 0; i != l; ++i) {
//     ord2 = ord[i] - min[g[i]]; // alternative idea: create 1D omap or access it as 1D, and then do the same as above !!
//     if(ord2 >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//     if(omap[cgs[g[i]]+ord2] == 0) omap[cgs[g[i]]+ord2] = i; // fastest ??
//     else stop("Repeated values of timevar within one or more groups"); 
//   }
//   // return omap;
//   for(int p = 0; p != ns; ++p) {
//     int np = n[p];
//     if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
//     if(np>0) { // Positive lagged and iterated differences
//       int d1 = diff[0]; 
//       bool L1 = np == 1;
//       if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//       NumericMatrix::Column outp = out( _ , pos); 
//       memset(seen, 0, memsize); 
//       if(L1) colnam[pos++] = ".D" + diffc[0]; 
//       else colnam[pos++] = ".L" + nc[p] + "D" + diffc[0]; 
//       for(int i = 0; i != l; ++i) { // Faster way, i.e. by using i2 = omap[i]; ?? 
//         if(seen[g[omap[i]]] == np) outp[omap[i]] = x[omap[i]] - x[omap[i - np]];
//         else {
//           outp[omap[i]] = fill;
//           ++seen[g[omap[i]]];
//         }
//       }
//       if(d1 > 1) for(int k = 1; k != d1; ++k) {
//         int start = np*(k+1); 
//         memset(seen, 0, memsize); 
//         for(int i = l; i--; ) { // correct iteration ?? 
//           if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outp[omap[i]] = fill;
//           else {
//             outp[omap[i]] -= outp[omap[i - np]];
//             ++seen[g[omap[i]]];
//           }
//         }
//       }
//       if(ds > 1) {
//         NumericVector outtemp = outp; 
//         for(int q = 1; q != ds; ++q) {
//           int dq = diff[q], L_dq = diff[q-1];
//           if(dq <= L_dq) stop("differences must be passed in ascending order");
//           for(int k = L_dq; k != dq; ++k) {
//             int start = np*(k+1); // Right ?? 
//             memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
//             for(int i = l; i--; ) {
//               if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outtemp[omap[i]] = fill;
//               else {
//                 outtemp[omap[i]] -= outtemp[omap[i - np]];
//                 ++seen[g[omap[i]]];
//               }
//             }
//           }
//           out( _ , pos) = outtemp; 
//           if(L1) colnam[pos++] = ".D" + diffc[q]; 
//           else colnam[pos++] = ".L" + nc[p] + "D" + diffc[q]; 
//         }
//       }
//     } else if(np<0) { // (Negative) leaded and iterated differences
//       int d1 = diff[0]; 
//       bool F1 = np == -1;
//       if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//       NumericMatrix::Column outp = out( _ , pos); 
//       memset(seen, 0, memsize);
//       if(F1) colnam[pos++] = ".FD" + diffc[0]; 
//       else colnam[pos++] = ".F" + nc[p] + "D" + diffc[0]; 
//       for(int i = l; i--; ) { // good??
//         if(seen[g[omap[i]]] == np) outp[omap[i]] = x[omap[i]] - x[omap[i - np]];
//         else {
//           outp[omap[i]] = fill;
//           --seen[g[omap[i]]];
//         }
//       }
//       if(d1 > 1) for(int k = 1; k != d1; ++k) {
//         int start = np*(k+1); // good ?? 
//         memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
//         for(int i = 0; i != l; ++i) {
//           if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outp[omap[i]] = fill; // was -start -> error !!
//           else {
//             outp[omap[i]] -= outp[omap[i - np]];
//             ++seen[g[omap[i]]];
//           }
//         }
//       }
//       if(ds > 1) {
//         NumericVector outtemp = outp; 
//         for(int q = 1; q != ds; ++q) {
//           int dq = diff[q], L_dq = diff[q-1];
//           if(dq <= L_dq) stop("differences must be passed in ascending order");
//           for(int k = L_dq; k != dq; ++k) {
//             int start = np*(k+1); // Right ?? 
//             memset(seen, 0, memsize); // missing !!!!!!!!!!!!!!!!! a problem ??, why does it work without ??
//             for(int i = 0; i != l; ++i) {
//               if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outtemp[omap[i]] = fill; // was -start -> error !!
//               else {
//                 outtemp[omap[i]] -= outtemp[omap[i - np]];
//                 ++seen[g[omap[i]]];
//               }
//             }
//           }
//           out( _ , pos) = outtemp; 
//           if(F1) colnam[pos++] = ".FD" + diffc[q]; 
//           else colnam[pos++] = ".F"+ nc[p] + "D" + diffc[q]; 
//         }
//       }
//     } else {
//       out( _ , pos) = x;
//       colnam[pos++] = ".--";
//     }
//   }
//   if(ncol == 1) SHALLOW_DUPLICATE_ATTRIB(out, x);
//   else out.attr("dimnames") = List::create(x.attr("names"), colnam);
//   return out;
// }



// Unordered panel lag other implementation 1: Make a copy of the vector for each iteration + 2d map  -> slowest !!!!!!!!!
// // [[Rcpp::export]]
// SEXP fdiffCpppuo1(NumericVector x, IntegerVector n = 1, IntegerVector diff = 1, double fill = NA_REAL, 
//                   int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP t = R_NilValue) { 
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
//   if(l != g.size()) stop("length(x) must match length(g)");
//   int ags = l/ng, maxdiff = diff[ds-1], ngp = ng+1; 
//   IntegerVector gsv = NULL; // The gsv solution above faster than std::fill ????????????????????
//   IntegerVector ord = t; 
//   if(l != ord.size()) stop("length(x) must match length(t)");
//   IntegerVector min(ngp, INT_MAX); 
//   IntegerVector ord2 = no_init_vector(l); 
//   if(Rf_isNull(gs)) { // needed for proper iteration 
//     gsv = IntegerVector(ng);
//     for(int i = 0; i != l; ++i) { 
//       ++gsv[g[i]-1];
//       if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//     }
//   } else {
//     gsv = gs;
//     if(ng != gsv.size()) stop("ng must match length(gs)");
//     for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//   }
//   int **omap = new int*[ngp]; 
//   for(int i = 0; i != ng; ++i) omap[i+1] = new int[gsv[i]]{}; 
//   for(int i = 0; i != l; ++i) {
//     ord2[i] = ord[i] - min[g[i]]; 
//     if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//     if(omap[g[i]][ord2[i]] == 0) omap[g[i]][ord2[i]] = i; 
//     else stop("Repeated values of timevar within one or more groups"); 
//   }
//   for(int p = 0; p != ns; ++p) {
//     int np = n[p];
//     if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
//     if(np>0) { // Positive lagged and iterated differences
//       int d1 = diff[0]; 
//       bool L1 = np == 1;
//       if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//       NumericMatrix::Column outp = out( _ , pos); 
//       if(L1) colnam[pos++] = ".D" + diffc[0]; 
//       else colnam[pos++] = ".L" + nc[p] + "D" + diffc[0]; 
//       for(int i = 0; i != l; ++i) { 
//         if(ord2[i] >= np) outp[i] = x[i] - x[omap[g[i]][ord2[i]-np]]; 
//         else outp[i] = fill;
//       }
//       if(d1 > 1) for(int k = 1; k != d1; ++k) {
//         int start = np*(k+1); 
//         NumericVector outpcopy = no_init_vector(l);
//         outpcopy = outp;
//         for(int i = l; i--; ) { 
//           if(ord2[i] >= start) outp[i] -= outpcopy[omap[g[i]][ord2[i]-np]]; 
//           else outp[i] = fill;
//         }
//       }
//       if(ds > 1) {
//         NumericVector outtemp = outp; 
//         for(int q = 1; q != ds; ++q) {
//           int dq = diff[q], L_dq = diff[q-1];
//           if(dq <= L_dq) stop("differences must be passed in ascending order");
//           for(int k = L_dq; k != dq; ++k) {
//             int start = np*(k+1); 
//             NumericVector outpcopy = no_init_vector(l);
//             outpcopy = outtemp;
//             for(int i = l; i--; ) {
//               if(ord2[i] >= start) outtemp[i] -= outpcopy[omap[g[i]][ord2[i]-np]]; 
//               else outtemp[i] = fill;
//             }
//           }
//           out( _ , pos) = outtemp; 
//           if(L1) colnam[pos++] = ".D" + diffc[q]; 
//           else colnam[pos++] = ".L" + nc[p] + "D" + diffc[q]; 
//         }
//       }
//     } else if(np<0) { // (Negative) leaded and iterated differences
//       int d1 = diff[0]; 
//       bool F1 = np == -1;
//       if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//       NumericMatrix::Column outp = out( _ , pos); 
//       if(F1) colnam[pos++] = ".FD" + diffc[0]; 
//       else colnam[pos++] = ".F" + nc[p] + "D" + diffc[0]; 
//       for(int i = l; i--; ) { // good??
//         if(ord2[i] < gsv[g[i]-1]+np) outp[i] = x[i] - x[omap[g[i]][ord2[i]-np]]; // probably fastest, but best ??
//         else outp[i] = fill;
//       }
//       if(d1 > 1) for(int k = 1; k != d1; ++k) {
//         int start = np*(k+1); // good ?? 
//         NumericVector outpcopy = no_init_vector(l);
//         outpcopy = outp;
//         for(int i = 0; i != l; ++i) {
//           if(ord2[i] < gsv[g[i]-1]+start) outp[i] -= outpcopy[omap[g[i]][ord2[i]-np]]; // probably fastest, but best ??
//           else outp[i] = fill;
//         }
//       }
//       if(ds > 1) {
//         NumericVector outtemp = outp; 
//         for(int q = 1; q != ds; ++q) {
//           int dq = diff[q], L_dq = diff[q-1];
//           if(dq <= L_dq) stop("differences must be passed in ascending order");
//           for(int k = L_dq; k != dq; ++k) {
//             int start = np*(k+1); // Right ?? 
//             NumericVector outpcopy = no_init_vector(l);
//             outpcopy = outtemp;
//             for(int i = 0; i != l; ++i) {
//               if(ord2[i] < gsv[g[i]-1]+start) outtemp[i] -= outpcopy[omap[g[i]][ord2[i]-np]]; // probably fastest, but best ??
//               else outtemp[i] = fill;
//             }
//           }
//           out( _ , pos) = outtemp; 
//           if(F1) colnam[pos++] = ".FD" + diffc[q]; 
//           else colnam[pos++] = ".F"+ nc[p] + "D" + diffc[q]; 
//         }
//       }
//     } else {
//       out( _ , pos) = x;
//       colnam[pos++] = ".--";
//     }
//   }
//   if(ncol == 1) SHALLOW_DUPLICATE_ATTRIB(out, x);
//   else out.attr("dimnames") = List::create(x.attr("names"), colnam);
//   return out;
// }



// Earlier attempt of unordered panel-lag: wrong iteration for both positive and negative differences, but a lot of trying around !!
// IntegerVector ord = t; 
// if(l != ord.size()) stop("length(x) must match length(t)");
// IntegerVector min(ngp, INT_MAX); // Best solution ?? array ?? 
// IntegerVector ord2 = no_init_vector(l); 
// if(Rf_isNull(gs)) { // needed for proper iteration 
//   gsv = IntegerVector(ng);
//   for(int i = 0; i != l; ++i) { 
//     ++gsv[g[i]-1];
//     if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
//   }
// } else {
//   gsv = gs;
//   if(ng != gsv.size()) stop("ng must match length(gs)");
//   for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
// }
// int **omap = new int*[ngp]; // best ?? 
// for(int i = 0; i != ng; ++i) omap[i+1] = new int[gsv[i]]{}; // best ?? // or () // better using vector of vectors (faster initialization) ??
// for(int i = 0; i != l; ++i) {
//   ord2[i] = ord[i] - min[g[i]]; // alternative idea: create 1D omap or access it as 1D, and then do the same as above !!
//   if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
//   if(omap[g[i]][ord2[i]] == 0) omap[g[i]][ord2[i]] = i; // fastest ??
//   else stop("Repeated values of timevar within one or more groups"); 
// }
// for(int p = 0; p != ns; ++p) {
//   int np = n[p];
//   if(absn[p]*maxdiff > ags) stop("abs(n * diff) exceeds average group size: %i", ags);
//   if(np>0) { // Positive lagged and iterated differences
//     int d1 = diff[0]; 
//     bool L1 = np == 1;
//     if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//     NumericMatrix::Column outp = out( _ , pos); 
//     if(L1) colnam[pos++] = ".D" + diffc[0]; 
//     else colnam[pos++] = ".L" + nc[p] + "D" + diffc[0]; 
//     // This only works once !! -> 4 OPTIONS: 
//     // - Make a copy of the vector for each iteration
//     // - Work with 1D array any only access thrugh that -> same ad before, need counter
//     // - for further iterations access the 2d array as 1D !!
//     // - Double loop: Groups and elements !!
//     for(int i = 0; i != l; ++i) { 
//       if(ord2[i] >= np) outp[i] = x[i] - x[omap[g[i]][ord2[i]-np]]; // probably fastest, but best ??
//       else outp[i] = fill;
//     }
//     if(d1 > 1) for(int k = 1; k != d1; ++k) {
//       int start = np*(k+1); // right ????
//       for(int i = l; i--; ) { // all needs to be in order !! just repeating the above won't work !!
//         if(ord2[i] >= start) outp[i] -= outp[omap[g[i]][ord2[i]-np]]; // probably fastest, but best ??
//         else outp[i] = fill;
//         // if(seen[g[i]] == gsv[g[i]-1]-start) outp[i] = fill;
//         // else {
//         //   outp[i] -= outp[i - np];
//         //   ++seen[g[i]];
//         // }
//       }
//     }
//     if(ds > 1) {
//       NumericVector outtemp = outp; 
//       for(int q = 1; q != ds; ++q) {
//         int dq = diff[q], L_dq = diff[q-1];
//         if(dq <= L_dq) stop("differences must be passed in ascending order");
//         for(int k = L_dq; k != dq; ++k) {
//           int start = np*(k+1); // Right ?? 
//           for(int i = l; i--; ) {
//             if(ord2[i] >= start) outtemp[i] -= outtemp[omap[g[i]][ord2[i]-np]]; // probably fastest, but best ??
//             else outtemp[i] = fill;
//             // if(seen[g[i]] == gsv[g[i]-1]-start) outtemp[i] = fill;
//             // else {
//             //   outtemp[i] -= outtemp[i - np];
//             //   ++seen[g[i]];
//             // }
//           }
//         }
//         out( _ , pos) = outtemp; 
//         if(L1) colnam[pos++] = ".D" + diffc[q]; 
//         else colnam[pos++] = ".L" + nc[p] + "D" + diffc[q]; 
//       }
//     }
//   } else if(np<0) { // (Negative) leaded and iterated differences
//     int d1 = diff[0]; 
//     bool F1 = np == -1;
//     if(d1 < 1) stop("diff must be a vector of integers > 0"); 
//     NumericMatrix::Column outp = out( _ , pos); 
//     if(F1) colnam[pos++] = ".FD" + diffc[0]; 
//     else colnam[pos++] = ".F" + nc[p] + "D" + diffc[0]; 
//     for(int i = l; i--; ) { // good??
//       if(ord2[i] < gsv[g[i]-1]+np) outp[i] = x[i] - x[omap[g[i]][ord2[i]-np]]; // probably fastest, but best ??
//       else outp[i] = fill;
//       // if(seen[g[i]] == np) outp[i] = x[i] - x[i - np];
//       // else {
//       //   outp[i] = fill;
//       //   --seen[g[i]];
//       // }
//     }
//     if(d1 > 1) for(int k = 1; k != d1; ++k) {
//       int start = np*(k+1); // good ?? 
//       for(int i = 0; i != l; ++i) {
//         if(ord2[i] < gsv[g[i]-1]+start) outp[i] -= outp[omap[g[i]][ord2[i]-np]]; // probably fastest, but best ??
//         else outp[i] = fill;
//         // if(seen[g[i]] == gsv[g[i]-1]-start) outp[i] = fill;
//         // else {
//         //   outp[i] -= outp[i - np];
//         //   ++seen[g[i]];
//         // }
//       }
//     }
//     if(ds > 1) {
//       NumericVector outtemp = outp; 
//       for(int q = 1; q != ds; ++q) {
//         int dq = diff[q], L_dq = diff[q-1];
//         if(dq <= L_dq) stop("differences must be passed in ascending order");
//         for(int k = L_dq; k != dq; ++k) {
//           int start = np*(k+1); // Right ?? 
//           for(int i = 0; i != l; ++i) {
//             if(ord2[i] < gsv[g[i]-1]+start) outtemp[i] -= outtemp[omap[g[i]][ord2[i]-np]]; // probably fastest, but best ??
//             else outtemp[i] = fill;
//             // if(seen[g[i]] == gsv[g[i]-1]-start) outtemp[i] = fill;
//             // else {
//             //   outtemp[i] -= outtemp[i - np];
//             //   ++seen[g[i]];
//             // }
//           }
//         }
//         out( _ , pos) = outtemp; 
//         if(F1) colnam[pos++] = ".FD" + diffc[q]; 
//         else colnam[pos++] = ".F"+ nc[p] + "D" + diffc[q]; 
//       }
//     }
//   } else {
//     out( _ , pos) = x;
//     colnam[pos++] = ".--";
//   }
// }


// Earlier attempt of ordered panel-lag: Wrong iteration for both positive and negative difference !! -> Need to run through from the other side!!
// if(Rf_isNull(t)) { // Ordered data
//   int seen[ng], memsize = sizeof(int)*ng;
//   for(int p = 0; p != ns; ++p) {
//     int np = n[p];
//     if(np>0) { // Positive lagged and iterated differences
//       int d1 = diff[0];
//       bool L1 = np == 1;
//       if(d1 < 1) stop("diff must be a vector of integers > 0");
//       NumericMatrix::Column outp = out( _ , pos);
//       memset(seen, 0, memsize);
//       if(L1) colnam[pos++] = ".D" + diffc[0];
//       else colnam[pos++] = ".L" + nc[p] + "D" + diffc[0];
//       for(int i = 0; i != l; ++i) { // this loop ??
//         if(seen[g[i]-1] == np) outp[i] = x[i] - x[i-np];
//         else {
//           outp[i] = fill;
//           ++seen[g[i]-1];
//         }
//       }
//       if(d1 > 1) for(int k = 1; k != d1; ++k) {
//         int start = np*(k+1); // right ????
//         memset(seen, 0, memsize); // Needed, because it loops from the beginning !!
//         for(int i = l; i--; ) { // correct iteration ??
//           if(seen[g[i]-1] == start) outp[i] -= outp[i-np];
//           else {
//             outp[i] = fill;
//             ++seen[g[i]-1];
//           }
//         }
//       }
//       if(ds > 1) {
//         NumericVector outtemp = outp;
//         for(int q = 1; q != ds; ++q) {
//           int dq = diff[q];
//           for(int k = diff[q-1]; k != dq; ++k) {
//             int start = np*k-1; // Right ?? Also no new memset needed ??
//             for(int i = l; i--; ) {
//               if(seen[g[i]-1] == start) outtemp[i] -= outtemp[i-np];
//               else {
//                 outp[i] = fill;
//                 ++seen[g[i]-1];
//               }
//             }
//           }
//           out( _ , pos) = outtemp;
//           if(L1) colnam[pos++] = ".D" + diffc[q];
//           else colnam[pos++] = ".L" + nc[p] + "D" + diffc[q];
//         }
//       }
//     } else if(np<0) { // (Negative) leaded and iterated differences
//       int d1 = diff[0];
//       bool F1 = np == -1;
//       if(d1 < 1) stop("diff must be a vector of integers > 0");
//       NumericMatrix::Column outp = out( _ , pos);
//       memset(seen, 0, memsize);
//       if(F1) colnam[pos++] = ".FD" + diffc[0];
//       else colnam[pos++] = ".F" + nc[p] + "D" + diffc[0];
//       for(int i = l; i--; ) { // good??
//         if(seen[g[i]-1] == np) outp[i] = x[i] - x[i-np];
//         else {
//           outp[i] = fill;
//           --seen[g[i]-1];
//         }
//       }
//       if(d1 > 1) for(int k = 1; k != d1; ++k) {
//         int final = l+np*k; // good ??
//         for(int i = 0; i != l; ++i) { // good?? correct iteration ??
//           if(seen[g[i]-1] == final) outp[i] -= outp[i-np];
//           else {
//             outp[i] = fill;
//             --seen[g[i]-1];
//           }
//         }
//       }
//       if(ds > 1) {
//         NumericVector outtemp = outp;
//         for(int q = 1; q != ds; ++q) {
//           int dq = diff[q];
//           for(int k = diff[q-1]; k != dq; ++k) {
//             int final = l+np*k; // good ??
//             for(int i = 0; i != l; ++i) { // good?? correct iteration ??
//               if(seen[g[i]-1] == final) outtemp[i] -= outtemp[i-np];
//               else {
//                 outtemp[i] = fill;
//                 --seen[g[i]-1];
//               }
//             }
//           }
//           out( _ , pos) = outtemp;
//           if(F1) colnam[pos++] = ".FD" + diffc[q];
//           else colnam[pos++] = ".F"+ nc[p] + "D" + diffc[q];
//         }
//       }
//     } else {
//       out( _ , pos) = x;
//       colnam[pos++] = ".--";
//     }
//   }
// }


// Earlier plaxing around with unordered differences (non-panel)!! -> not working well, but many approaches and ideas !!
// IntegerVector ord = t;
// if(l != ord.size()) stop("length(x) must match length(t)");
// LogicalVector ocheck(l, true); // Could do more efficiently if only one value -> but you don't use this anyway!!
// int omap[l+1]; // otherwise always ord[i]-1 -> efficiency gain !! -> check also for panel-lag !!!!!!
// for(int i = 0; i != l; ++i) { // integrated below !!
//   if(ord[i] > l) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
//   if(ocheck[ord[i]-1]) ocheck[ord[i]-1] = false;
//   else stop("Repeated values in timevar");
//   omap[ord[i]] = i; // Note: omap is the same as order(ord) !!
// }
// // return omap;
// for(int p = 0; p != ns; ++p) {
//   int np = n[p];
//   if(np>0) { // Positive lagged and iterated differences
//     int d1 = diff[0];
//     bool L1 = np == 1;
//     if(d1 < 1) stop("diff must be a vector of integers > 0");
//     if(np*d1 >= l) stop("n * diff needs to be < length(x)");
//     NumericMatrix::Column outp = out( _ , pos);
//     if(L1) colnam[pos++] = ".D" + diffc[0];
//     else colnam[pos++] = ".L" + nc[p] + "D" + diffc[0];
//     for(int i = 0; i != l; ++i) {
//       if(ord[i] > np) outp[i] = x[i] - x[omap[ord[i]-np]];
//       else outp[i] = fill;
//     }
//     if(d1 > 1) for(int k = 1; k != d1; ++k) {
//       int start = np*(k+1); // -1 right?? // Note: culd just do for(i in 1:10) if(i>1) outp[omap[i]] = y[omap[i]] - y[omap[i-1]] else outp[omap[i]] = NA and continue as usual
//       for(int i = l; i != 0; --i) { // Needed for correct iteration (omap runs from 1 to l)
//         if(ord[omap[i]] > start) outp[omap[i]] -= outp[omap[i-np]]; // omap must not shift !!
//         else outp[omap[i]] = fill;
//       }
//       // for(int i = l; i--; ) { // Needed for correct iteration
//       //   if(ord[omap[i+1+k]] > start) outp[i] -= outp[omap[ord[omap[i+1+k]]-np+k]]; // +k needed because upon iteration also the omap must shift ??
//       //   else outp[omap[i+1+k]] = fill;
//       // }
//       // for(int i = l; i--; ) { // Needed for correct iteration
//       //   if(ord[i] > start) outp[i] -= outp[omap[ord[i]-np+k]]; // +k needed because upon iteration also the omap must shift ??
//       //   else outp[i] = fill;
//       // }
//     }
//     if(ds > 1) {
//       NumericVector outtemp = outp;
//       for(int q = 1; q != ds; ++q) {
//         int dq = diff[q], L_dq = diff[q-1];
//         if(np*dq >= l) stop("n * diff needs to be < length(x)");
//         if(dq <= L_dq) stop("differences must be passed in ascending order");
//         for(int k = L_dq; k != dq; ++k) {
//           int start = np*(k+1); // -1 right??
//           for(int i = l; i != 0; --i) { // Needed for correct iteration (omap runs from 1 to l)
//             if(ord[omap[i]] > start) outtemp[omap[i]] -= outtemp[omap[i-np]]; // omap must not shift !!
//             else outp[omap[i]] = fill;
//           }
//           // for(int i = l; i--; ) { // Needed for correct iteration
//           //   if(ord[i] > start) outtemp[i] -= outtemp[omap[ord[i]-np+k]]; // +k needed because upon iteration also the omap must shift ??
//           //   else outtemp[i] = fill;
//           // }
//         }
//         out( _ , pos) = outtemp;
//         if(L1) colnam[pos++] = ".D" + diffc[q];
//         else colnam[pos++] = ".L" + nc[p] + "D" + diffc[q];
//       }
//     }
//   } else if(np<0) { // (Negative) leaded and iterated differences
//     int d1 = diff[0], start = l+np, end = l+np*d1;
//     bool F1 = np == -1;
//     if(d1 < 1) stop("diff must be a vector of integers > 0");
//     if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
//     NumericMatrix::Column outp = out( _ , pos);
//     if(F1) colnam[pos++] = ".FD" + diffc[0];
//     else colnam[pos++] = ".F" + nc[p] + "D" + diffc[0];
//     for(int i = 0; i != l; ++i) { // Like this or nested loop ??
//       if(ord[i] <= start) outp[i] = x[i] - x[omap[ord[i]-np]];
//       else outp[i] = fill;
//     }
//     if(d1 > 1) for(int k = 1; k != d1; ++k) {
//       int final = l+np*(k+1); // good ??
//       for(int i = 0; i != l; ++i) {
//         if(ord[i] <= final) outp[i] -= outp[omap[ord[i]-np]];
//         else outp[i] = fill;
//       }
//     }
//     if(ds > 1) {
//       NumericVector outtemp = outp;
//       for(int q = 1; q != ds; ++q) {
//         int dq = diff[q], L_dq = diff[q-1], end = l+np*dq;
//         if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
//         if(dq <= L_dq) stop("differences must be passed in ascending order");
//         for(int k = L_dq; k != dq; ++k) {
//           int final = l+np*(k+1); // good ??
//           for(int i = 0; i != l; ++i) {
//             if(ord[i] <= final) outtemp[i] -= outtemp[omap[ord[i]-np]];
//             else outp[i] = fill;
//           }
//         }
//         out( _ , pos) = outtemp;
//         if(F1) colnam[pos++] = ".FD" + diffc[q];
//         else colnam[pos++] = ".F"+ nc[p] + "D" + diffc[q];
//       }
//     }
//   } else {
//     out( _ , pos) = x;
//     colnam[pos++] = ".--";
//   }
// }
