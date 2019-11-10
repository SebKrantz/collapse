// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

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
      first.attr("row.names") = NumericVector::create(NA_REAL, -ng);
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
        first.attr("row.names") = NumericVector::create(NA_REAL, -ng);
      }
    }
    return first;
  }
}

// previous groupd na.rm = FALSE version: (new solution is much faster !!!)

// for(int j = l; j--; ) {
//   int ngs = 0;
//   switch(TYPEOF(x[j])) { // useing std::vector<bool> instead here is more memory efficient but not faster !!!!!!!!!!!!!!!
//   case REALSXP: {
//     NumericVector column = x[j];
//     if(gss != column.size()) stop("length(g) must match nrow(X)");
//     NumericVector firstj = no_init_vector(ng); // stable? -> Yes !! (only bad for summation commands)
//     LogicalVector glj(ng, true); // Other way around ?? -> Nope, this is faster !!
//     for(int i = 0; i != gss; ++i) {
//       if(glj[g[i]-1]) {
//         glj[g[i]-1] = false;
//         firstj[g[i]-1] = column[i];
//         ++ngs;
//         if(ngs == ng) break;
//       }
//     }          
//     SHALLOW_DUPLICATE_ATTRIB(firstj, column);
//     first[j] = firstj;
//     break;    
//   }
//   case INTSXP: {
//     IntegerVector column = x[j];
//     if(gss != column.size()) stop("length(g) must match nrow(X)");
//     IntegerVector firstj = no_init_vector(ng); 
//     LogicalVector glj(ng, true); 
//     for(int i = 0; i != gss; ++i) {
//       if(glj[g[i]-1]) {
//         glj[g[i]-1] = false;
//         firstj[g[i]-1] = column[i];
//         ++ngs;
//         if(ngs == ng) break;
//       }
//     }          
//     SHALLOW_DUPLICATE_ATTRIB(firstj, column);
//     first[j] = firstj;
//     break;    
//   }
//   case STRSXP: {
//     CharacterVector column = x[j];
//     if(gss != column.size()) stop("length(g) must match nrow(X)");
//     CharacterVector firstj = no_init_vector(ng); 
//     LogicalVector glj(ng, true); 
//     for(int i = 0; i != gss; ++i) {
//       if(glj[g[i]-1]) {
//         glj[g[i]-1] = false;
//         firstj[g[i]-1] = column[i];
//         ++ngs;
//         if(ngs == ng) break;
//       }
//     }          
//     SHALLOW_DUPLICATE_ATTRIB(firstj, column);
//     first[j] = firstj;
//     break;    
//   }
//   case LGLSXP: {
//     LogicalVector column = x[j];
//     if(gss != column.size()) stop("length(g) must match nrow(X)");
//     LogicalVector firstj = no_init_vector(ng); 
//     LogicalVector glj(ng, true); 
//     for(int i = 0; i != gss; ++i) {
//       if(glj[g[i]-1]) {
//         glj[g[i]-1] = false;
//         firstj[g[i]-1] = column[i];
//         ++ngs;
//         if(ngs == ng) break;
//       }
//     }          
//     SHALLOW_DUPLICATE_ATTRIB(firstj, column);
//     first[j] = firstj;
//     break;    
//   }
//   default: 
//     stop("incompatible SEXP encountered;");
//     break;
//   }
// }

// Previous Version: Before efficient attribute copy !!
// // [[Rcpp::export]]
// SEXP ffirstlCpp(List x, int ng = 0, IntegerVector g = 0, bool narm = true, bool drop = true) {
//   int l = x.size(), row = 0;
//   
//   if (ng == 0) {
//     List first(l); // Good and faster !! -> Only bad for summation commands !!
//     if(narm) {
//       for(int j = l; j--; ) {
//         int k = 0;
//         switch(TYPEOF(x[j])) { // Faster than using iterator ?? // https://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
//           case REALSXP: {
//             NumericVector column = x[j];
//             row = column.size();
//             double firsti = column[0];
//             while(std::isnan(firsti) && k!=row-1) firsti = column[++k];
//             if(column.hasAttribute("class")) {
//               // if(column.inherits("Date")) {
//                 NumericVector r(1, firsti); // Yes !! fully covers data and datetime // fastest ?? 
//                 CharacterVector an = wrap(column.attributeNames());
//                 for(int i = 0; i != an.size(); ++i) {
//                   String ss = an[i];
//                   r.attr(ss) = column.attr(ss);
//                 }
//                 // if(column.hasAttribute("tzone")) r.attr("tzone") = column.attr("tzone");
//                 first[j] = r;
//                 //   } else if (column.inherits("POSIXlt")) {
//                 //   NumericVector r(1, firsti);
//                 //   r.attr("class") = column.attr("class");
//                 //   r.attr("names") = column.attr("names");
//                 //   first[j] = r;
//                 // } else {
//                 // first[j] = firsti; 
//                 // }
//               // first[j] = DateVector::create(column[k]); // not working !!
//             } else { // Could also do datetime here !!! 
//               first[j] = firsti; 
//             }
//             break;    
//           }
//           case INTSXP: {
//             IntegerVector column = x[j];
//             row = column.size();
//             int firsti = column[0];
//             while(firsti == NA_INTEGER && k!=row-1) firsti = column[++k];
//             if(Rf_isFactor(column)) { // or if(column.hasAttribute("levels")){ ?? 
//               IntegerVector r(1, firsti);
//               r.attr("class") = column.attr("class");
//               r.attr("levels") = column.attr("levels");
//               first[j] = r;
//             } else {
//               first[j] = firsti;
//             }
//             break;
//           }
//           case STRSXP: {
//             CharacterVector column = x[j];
//             row = column.size();
//             String firsti = column[0];
//             while(firsti == NA_STRING && k!=row-1) firsti = column[++k];
//             first[j] = firsti;
//             break;
//           }
//           case LGLSXP: {
//             LogicalVector column = x[j];
//             row = column.size();
//             bool firsti = column[0];
//             while(firsti == NA_LOGICAL && k!=row-1) firsti = column[++k];
//             first[j] = firsti;
//             break;
//           }
//          // case VECSXP: { // only for date times  !!
//         //    DatetimeVector column = x[j];
//         //    row = column.size();
//         //    Datetime firsti = column[0];
//         //    while(firsti != firsti && k!=row-1) firsti = column[++k];
//         //    // List r = List::create(firsti);
//         //    // r.attr("class") = column.attr("class");
//         //    // r.attr("names") = column.attr("names");
//         //    // first[j] = r;
//         //    // first[j] = column[k];
//         //    first[j] = firsti;
//         //    break;
//         //  }
//           default: {
//             stop("incompatible SEXP encountered;");
//           }
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) { // Faster than using iterator ?? // https://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
//           case REALSXP: {
//             NumericVector column = x[j];
//             if(column.hasAttribute("class")) {
//               NumericVector r(1, column[0]);
//               CharacterVector an = wrap(column.attributeNames());
//               for(int i = 0; i != an.size(); ++i) {
//                 String s(an[i]);
//                 r.attr(s) = column.attr(s);
//               }
//               first[j] = r;
//             } else {
//               first[j] = column[0];
//             }
//             break;    
//           }
//           case INTSXP: {
//             IntegerVector column = x[j];
//             if(Rf_isFactor(column)) { // or if(column.hasAttribute("levels")){ ?? 
//               IntegerVector r(1, column[0]);
//               r.attr("class") = column.attr("class");
//               r.attr("levels") = column.attr("levels");
//               first[j] = r;
//             } else {
//               first[j] = column[0]; 
//             }
//             break;
//           }
//           case STRSXP: {
//             CharacterVector column = x[j];
//             String s(column[0]);
//             first[j] = s; // good ?? -> Yes, this is necessary !!
//             break;
//           }
//           case LGLSXP: {
//             LogicalVector column = x[j];
//             // bool s(column[0]); // good ?? or LogicalVector::create
//             first[j] = LogicalVector::create(column[0]); // Better, because NA_LOGICAL is not standard C++ type !!
//             break;
//           }
//           default: {
//             stop("incompatible SEXP encountered;");
//           }
//         }
//       }
//     }
//     if(drop) first.attr("names") = x.attr("names");
//     return first;
//   } else { // With groups !!
//     List first(l);
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) { // Faster than using iterator ?? // https://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
//           case REALSXP: {
//             NumericVector column = x[j];
//             row = column.size();
//             if(gss != row) stop("length(g) must match nrow(X)");
//             NumericVector firstj(ng, NA_REAL);
//             int ngs = 0;
//             for(int i = 0; i != row; ++i) {
//               if(!std::isnan(column[i])) {
//                 if(std::isnan(firstj[g[i]-1])) {
//                   firstj[g[i]-1] = column[i];
//                   ++ngs;
//                   if(ngs == ng) break;
//                 }
//               } 
//             }
//             if(column.hasAttribute("class")) {
//               CharacterVector an = wrap(column.attributeNames());
//               for(int i = 0; i != an.size(); ++i) {
//                 String s(an[i]);
//                 firstj.attr(s) = column.attr(s);
//               }
//             } 
//             first[j] = firstj;
//             break;    
//           }
//           case INTSXP: {
//             IntegerVector column = x[j];
//             row = column.size();
//             if(gss != row) stop("length(g) must match nrow(X)");
//             IntegerVector firstj(ng, NA_INTEGER);
//             int ngs = 0;
//             for(int i = 0; i != row; ++i) {
//               if(column[i] != NA_INTEGER) { // https://teuder.github.io/rcpp4everyone_en/240_na_nan_inf.html
//                 if(firstj[g[i]-1] == NA_INTEGER) {
//                   firstj[g[i]-1] = column[i];
//                   ++ngs;
//                   if(ngs == ng) break;
//                 }
//               } 
//             }
//             if(Rf_isFactor(column)) { // or if(column.hasAttribute("levels")){ ?? 
//               firstj.attr("class") = column.attr("class");
//               firstj.attr("levels") = column.attr("levels");
//             } 
//             first[j] = firstj;
//             break;    
//           }
//           case STRSXP: {
//             CharacterVector column = x[j];
//             row = column.size();
//             if(gss != row) stop("length(g) must match nrow(X)");
//             CharacterVector firstj(ng, NA_STRING);
//             int ngs = 0;
//             for(int i = 0; i != row; ++i) {
//               if(column[i] != NA_STRING) {
//                 if(firstj[g[i]-1] == NA_STRING) {
//                   firstj[g[i]-1] = column[i];
//                   ++ngs;
//                   if(ngs == ng) break;
//                 }
//               } 
//             }
//             first[j] = firstj;
//             break;    
//           }
//           case LGLSXP: {
//             LogicalVector column = x[j];
//             row = column.size();
//             if(gss != row) stop("length(g) must match nrow(X)");
//             LogicalVector firstj(ng, NA_LOGICAL);
//             int ngs = 0;
//             for(int i = 0; i != row; ++i) {
//               if(column[i] != NA_LOGICAL) { // https://teuder.github.io/rcpp4everyone_en/240_na_nan_inf.html
//                 if(firstj[g[i]-1] == NA_LOGICAL) { // Works ??? fast ?? 
//                   firstj[g[i]-1] = column[i];
//                   ++ngs;
//                   if(ngs == ng) break;
//                 }
//               } 
//             }
//             first[j] = firstj;
//             break;    
//           }
//           default: {
//             stop("incompatible SEXP encountered;");
//           }
//         }
//       }
//     } else {
//       for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) { // Faster than using iterator ?? // https://gallery.rcpp.org/articles/rcpp-wrap-and-recurse/
//         case REALSXP: {
//           NumericVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           NumericVector firstj = no_init_vector(ng); // stable? -> Yes !! (only bad for summation commands)
//           LogicalVector glj(ng, true); // Other way around ?? -> Nope, this is faster !!
//           int ngs = 0;
//           for(int i = 0; i != row; ++i) {
//             if(glj[g[i]-1]) {
//               glj[g[i]-1] = false;
//               firstj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }          
//           if(column.hasAttribute("class")) {
//             CharacterVector an = wrap(column.attributeNames());
//             for(int i = 0; i != an.size(); ++i) {
//               String s(an[i]);
//               firstj.attr(s) = column.attr(s);
//             }
//           } 
//           first[j] = firstj;
//           break;    
//         }
//         case INTSXP: {
//           IntegerVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           IntegerVector firstj = no_init_vector(ng); // stable? -> Yes !! (only bad for summation commands)
//           LogicalVector glj(ng, true); // Other way around ?? -> Nope, this is faster !!
//           int ngs = 0;
//           for(int i = 0; i != row; ++i) {
//             if(glj[g[i]-1]) {
//               glj[g[i]-1] = false;
//               firstj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }          
//           if(Rf_isFactor(column)) { // or if(column.hasAttribute("levels")){ ?? 
//             firstj.attr("class") = column.attr("class");
//             firstj.attr("levels") = column.attr("levels");
//           } 
//           first[j] = firstj;
//           break;    
//         }
//         case STRSXP: {
//           CharacterVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           CharacterVector firstj = no_init_vector(ng); // stable? -> Yes !! (only bad for summation commands)
//           LogicalVector glj(ng, true); // Other way around ?? -> Nope, this is faster !!
//           int ngs = 0;
//           for(int i = 0; i != row; ++i) {
//             if(glj[g[i]-1]) {
//               glj[g[i]-1] = false;
//               firstj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }          
//           first[j] = firstj;
//           break;    
//         }
//         case LGLSXP: {
//           LogicalVector column = x[j];
//           row = column.size();
//           if(gss != row) stop("length(g) must match nrow(X)");
//           LogicalVector firstj = no_init_vector(ng); // stable? -> Yes !! (only bad for summation commands)
//           LogicalVector glj(ng, true); // Other way around ?? -> Nope, this is faster !!
//           int ngs = 0;
//           for(int i = 0; i != row; ++i) {
//             if(glj[g[i]-1]) {
//               glj[g[i]-1] = false;
//               firstj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           }          
//           first[j] = firstj;
//           break;    
//         }
//         default: {
//           stop("incompatible SEXP encountered;");
//         }
//         }
//       }
//     }
//     return first;
//   }
// }

// // [[Rcpp::export]]
// SEXP ffirstlCpp(List x, int ng = 0, IntegerVector g = 0, bool narm = true, bool drop = true){
//   RCPP_RETURN_VECTOR(ffirstlCppImpl, x, ng, g, narm, drop);
// }


// Only for Numeric Vectors
// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // [[Rcpp::export]]
// SEXP ffirstlCpp(List x, int ng = 0, IntegerVector g = 0,
//                bool narm = true, bool drop = true) {
//   int l = x.size(), row = 0;
//   
//   if (ng == 0) {
//     NumericVector first = no_init_vector(l); // Good and faster !! -> Only bad for summation commands !!
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         int k = 0;
//         double firsti = column[k];
//         while(std::isnan(firsti) && k!=row-1) firsti = column[++k];
//         first[j] = firsti;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         first[j] = column[0];
//       }
//     }
//     if(drop) {
//       first.attr("names") = x.attr("names");
//       return first;
//     } else {
//       List out(l);
//       for(int j = l; j--; ) out[j] = first[j];
//       return out;
//     }
//   } else { // With groups !!
//     List first(l);
//     int gss = g.size();
//     if(narm) {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector firstj(ng, NA_REAL);
//         int ngs = 0;
//         for(int i = 0; i != row; ++i) {
//           if(!std::isnan(column[i])) {
//             if(std::isnan(firstj[g[i]-1])) {
//               firstj[g[i]-1] = column[i];
//               ++ngs;
//               if(ngs == ng) break;
//             }
//           } 
//         }
//         first[j] = firstj;
//       }
//     } else {
//       for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(gss != row) stop("length(g) must match nrow(X)");
//         NumericVector firstj = no_init_vector(ng); // stable? -> Yes !! (only bad for summation commands)
//         LogicalVector glj(ng, true); // Other way around ?? -> Nope, this is faster !!
//         int ngs = 0;
//         for(int i = 0; i != row; ++i) {
//           if(glj[g[i]-1]) {
//             glj[g[i]-1] = false;
//             firstj[g[i]-1] = column[i];
//             ++ngs;
//             if(ngs == ng) break;
//           }
//         }
//         first[j] = firstj;
//       }
//     }
//     return first;
//   }
// }
