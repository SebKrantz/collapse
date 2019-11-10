// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

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
          CharacterVector out(1, column[k]);
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
          CharacterVector out(1, *(column.end()-1));
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
    last.attr("row.names") = 1;
    return last;
  } else { // With groups !!
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
      last.attr("row.names") = NumericVector::create(NA_REAL, -ng);
    } else {
      LogicalVector glj(ng, true); //  Much faster method !! (precomputing indices and then going through data)
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
        const CharacterVector& rn = Rf_getAttrib(x, R_RowNamesSymbol); // const doesn't really make a difference !!
        Rf_setAttrib(last, R_RowNamesSymbol, rn[lastindex]); //  last.attr("row.names") = rn[lastindex]; // Other sloghtly faster, but no big deal !!
      } else {
        last.attr("row.names") = NumericVector::create(NA_REAL, -ng);
      }
    }
    return last;
  }
}


// Previous Version: Before efficient attribute copy !!
// // [[Rcpp::export]]
// SEXP flastlCpp(List x, int ng = 0, IntegerVector g = 0, bool narm = true, bool drop = true) {
//   int l = x.size();
//   
//   if (ng == 0) {
//     List last(l); 
//     if(narm) {
//       for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) { // Faster than using iterator ?? // https://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
//         case REALSXP: {
//           NumericVector column = x[j];
//           int k = column.size()-1;
//           double lasti = column[k];
//           while(std::isnan(lasti) && k!=0) lasti = column[--k];
//           if(column.hasAttribute("class")) {
//             NumericVector r(1, lasti);
//             CharacterVector an = wrap(column.attributeNames());
//             for(int i = 0; i != an.size(); ++i) {
//               String ss = an[i];
//               r.attr(ss) = column.attr(ss);
//             }
//             last[j] = r;
//           } else { 
//             last[j] = lasti; 
//           }
//           break;    
//         }
//         case INTSXP: {
//           IntegerVector column = x[j];
//           int k = column.size()-1, lasti = column[k];
//           while(lasti == NA_INTEGER && k!=0) lasti = column[--k];
//           if(Rf_isFactor(column)) { // or if(column.hasAttribute("levels")){ ?? 
//             IntegerVector r(1, lasti);
//             r.attr("class") = column.attr("class");
//             r.attr("levels") = column.attr("levels");
//             last[j] = r;
//           } else {
//             last[j] = lasti;
//           }
//           break;
//         }
//         case STRSXP: {
//           CharacterVector column = x[j];
//           int k = column.size()-1;
//           String lasti = column[k];
//           while(lasti == NA_STRING && k!=0) lasti = column[--k];
//           last[j] = lasti;
//           break;
//         }
//         case LGLSXP: {
//           LogicalVector column = x[j];
//           int k = column.size()-1;
//           bool lasti = column[k];
//           while(lasti == NA_LOGICAL && k!=0) lasti = column[--k];
//           last[j] = lasti;
//           break;
//         }
//         default: {
//           stop("incompatible SEXP encountered;");
//         }
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) { // Faster than using iterator ?? // https://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
//         case REALSXP: {
//           NumericVector column = x[j];
//           if(column.hasAttribute("class")) {
//             NumericVector r(1, *(column.end()-1));
//             CharacterVector an = wrap(column.attributeNames());
//             for(int i = 0; i != an.size(); ++i) {
//               String s(an[i]);
//               r.attr(s) = column.attr(s);
//             }
//             last[j] = r;
//           } else {
//             last[j] = *(column.end()-1);
//           }
//           break;    
//         }
//         case INTSXP: {
//           IntegerVector column = x[j];
//           if(Rf_isFactor(column)) { // or if(column.hasAttribute("levels")){ ?? 
//             IntegerVector r(1, *(column.end()-1));
//             r.attr("class") = column.attr("class");
//             r.attr("levels") = column.attr("levels");
//             last[j] = r;
//           } else {
//             last[j] = *(column.end()-1); 
//           }
//           break;
//         }
//         case STRSXP: {
//           CharacterVector column = x[j];
//           String s(*(column.end()-1));
//           last[j] = s; 
//           break;
//         }
//         case LGLSXP: {
//           LogicalVector column = x[j];
//           last[j] = LogicalVector::create(*(column.end()-1)); 
//           break;
//         }
//         default: {
//           stop("incompatible SEXP encountered;");
//         }
//         }
//       }
//     }
//     if(drop) last.attr("names") = x.attr("names");
//     return last;
//   } else { // With groups !!
//     int row = 0;
//     List last(l);
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) { // Faster than using iterator ?? // https://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
//         case REALSXP: {
//           NumericVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           NumericVector lastj(ng, NA_REAL);
//           int ngs = 0;
//           for(int i = row; i--; ) {
//             if(!std::isnan(column[i])) {
//               if(std::isnan(lastj[g[i]-1])) {
//                 lastj[g[i]-1] = column[i];
//                 ++ngs;
//                 if(ngs == ng) break;
//               }
//             } 
//           }
//           if(column.hasAttribute("class")) {
//             CharacterVector an = wrap(column.attributeNames());
//             for(int i = 0; i != an.size(); ++i) {
//               String s(an[i]);
//               lastj.attr(s) = column.attr(s);
//             }
//           } 
//           last[j] = lastj;
//           break;    
//         }
//         case INTSXP: {
//           IntegerVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           IntegerVector lastj(ng, NA_INTEGER);
//           int ngs = 0;
//           for(int i = row; i--; ) {
//             if(column[i] != NA_INTEGER) {
//               if(lastj[g[i]-1] == NA_INTEGER) {
//                 lastj[g[i]-1] = column[i];
//                 ++ngs;
//                 if(ngs == ng) break;
//               }
//             } 
//           }
//           if(Rf_isFactor(column)) { // or if(column.hasAttribute("levels")){ ?? 
//             lastj.attr("class") = column.attr("class");
//             lastj.attr("levels") = column.attr("levels");
//           } 
//           last[j] = lastj;
//           break;    
//         }
//         case STRSXP: {
//           CharacterVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           CharacterVector lastj(ng, NA_STRING);
//           int ngs = 0;
//           for(int i = row; i--; ) {
//             if(column[i] != NA_STRING) {
//               if(lastj[g[i]-1] == NA_STRING) {
//                 lastj[g[i]-1] = column[i];
//                 ++ngs;
//                 if(ngs == ng) break;
//               }
//             } 
//           }
//           last[j] = lastj;
//           break;    
//         }
//         case LGLSXP: {
//           LogicalVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           LogicalVector lastj(ng, NA_LOGICAL);
//           int ngs = 0;
//           for(int i = row; i--; ) {
//             if(column[i] != NA_LOGICAL) { 
//               if(lastj[g[i]-1] == NA_LOGICAL) { 
//                 lastj[g[i]-1] = column[i];
//                 ++ngs;
//                 if(ngs == ng) break;
//               }
//             } 
//           }
//           last[j] = lastj;
//           break;    
//         }
//         default: {
//           stop("incompatible SEXP encountered;");
//         }
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) { // Faster than using iterator ?? // https://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
//         case REALSXP: {
//           NumericVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           NumericVector lastj = no_init_vector(ng); 
//           LogicalVector glj(ng, true); 
//           int ngs = 0;
//           for(int i = row; i--; ) {
//             if(glj[g[i]-1]) {
//               glj[g[i]-1] = false;
//               lastj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }          
//           if(column.hasAttribute("class")) {
//             CharacterVector an = wrap(column.attributeNames());
//             for(int i = 0; i != an.size(); ++i) {
//               String s(an[i]);
//               lastj.attr(s) = column.attr(s);
//             }
//           } 
//           last[j] = lastj;
//           break;    
//         }
//         case INTSXP: {
//           IntegerVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           IntegerVector lastj = no_init_vector(ng); 
//           LogicalVector glj(ng, true); 
//           int ngs = 0;
//           for(int i = row; i--; ) {
//             if(glj[g[i]-1]) {
//               glj[g[i]-1] = false;
//               lastj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }          
//           if(Rf_isFactor(column)) { // or if(column.hasAttribute("levels")){ ?? 
//             lastj.attr("class") = column.attr("class");
//             lastj.attr("levels") = column.attr("levels");
//           } 
//           last[j] = lastj;
//           break;    
//         }
//         case STRSXP: {
//           CharacterVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           CharacterVector lastj = no_init_vector(ng); 
//           LogicalVector glj(ng, true); 
//           int ngs = 0;
//           for(int i = row; i--; ) {
//             if(glj[g[i]-1]) {
//               glj[g[i]-1] = false;
//               lastj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }          
//           last[j] = lastj;
//           break;    
//         }
//         case LGLSXP: {
//           LogicalVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           LogicalVector lastj = no_init_vector(ng); 
//           LogicalVector glj(ng, true); 
//           int ngs = 0;
//           for(int i = row; i--; ) {
//             if(glj[g[i]-1]) {
//               glj[g[i]-1] = false;
//               lastj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }          
//           last[j] = lastj;
//           break;    
//         }
//         default: {
//           stop("incompatible SEXP encountered;");
//         }
//         }
//       }
//     }
//     return last;
//   }
// }


// Only Numeric Version
// // [[Rcpp::export]]
// SEXP flastlCpp(List x, int ng = 0, IntegerVector g = 0,
//                 bool narm = true, bool drop = true) {
//   int l = x.size(), row = 0;
//   
//   if (ng == 0) {
//     NumericVector last = no_init_vector(l); // Good and faster !! -> Only bad for summation commands !!
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         int k = row-1;
//         double lasti = column[k];
//         while(std::isnan(lasti) && k!=0) lasti = column[--k];
//         last[j] = lasti;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         last[j] = column[row-1];
//       }
//     }
//     if(drop) {
//       last.attr("names") = x.attr("names");
//       return last;
//     } else {
//       List out(l);
//       for(int j = l; j--; ) out[j] = last[j];
//       return out;
//     }
//   } else { // With groups !!
//     List last(l);
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector lastj(ng, NA_REAL);
//         int ngs = 0;
//         for(int i = row; i--; ) {
//           if(!std::isnan(column[i])) {
//             if(std::isnan(lastj[g[i]-1])) {
//               lastj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           } 
//         }
//         last[j] = lastj;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector lastj = no_init_vector(ng); // stable? -> Yes !! (only bad for summation commands)
//         LogicalVector glj(ng, true); // Other way around ?? -> Nope, this is faster !!
//         int ngs = 0;
//         for(int i = row; i--; ) {
//           if(glj[g[i]-1]) {
//             glj[g[i]-1] = false;
//             lastj[g[i]-1] = column[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         }
//         last[j] = lastj;
//       }
//     }
//     return last;
//   }
// }
