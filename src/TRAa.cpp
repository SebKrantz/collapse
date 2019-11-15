#include <Rcpp.h>
using namespace Rcpp;

// Cases:
// 1- replace
// 2- replace with NA rm
// 3- demean
// 4- demean with global mean added
// 5- Proportion
// 6- Percentages
// 7- Add
// 8- Multiply

// Todo:: Check !!

// [[Rcpp::export]]
SEXP TRAmCpp(const SEXP& x, const SEXP& xAG, const IntegerVector& g = 0, int ret = 1) {
  int gs = g.size();

  if(ret <= 2) {
    switch(TYPEOF(x)) {
    case REALSXP:
    {
      NumericMatrix xx = x;
      int l = xx.nrow(), col = xx.ncol();
      NumericMatrix out = no_init_matrix(l, col);
      if(gs == 1) {
        NumericVector AG = xAG;
        if(AG.size() != col) stop("If g = NULL, length(STATS) must match ncol(x)");
        switch(ret) {
        case 1:
          for(int j = col; j--; ) out(_,j) = rep(AG[j], l);
          break;
        case 2:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            double sumj = AG[j];
            for(int i = l; i--; ) {
              if(std::isnan(column[i])) colo[i] = column[i];
              else colo[i] = sumj;
            }
          }
          break;
        default: stop("Unknown Transformation");
        }
      } else {
        NumericMatrix AG = xAG;
        if(AG.ncol() != col) stop("ncol(STATS) must match ncol(x)");
        if(gs != l) stop("length(g) must match nrow(x)");
        switch(ret) {
        case 1:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) colo[i] = sumj[g[i]-1];
          }
          break;
        case 2:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            NumericMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) {
              if(std::isnan(column[i])) colo[i] = column[i];
              else colo[i] = sumj[g[i]-1];
            }
          }
          break;
        default: stop("Unknown Transformation");
        }
      }
      DUPLICATE_ATTRIB(out, xx); // or x ?? which is faster ??
      return out;
    }
    case INTSXP:
    {
      IntegerMatrix xx = x;
      int l = xx.nrow(), col = xx.ncol();
      IntegerMatrix out = no_init_matrix(l, col);
      if(gs == 1) {
        IntegerVector AG = xAG;
        if(AG.size() != col) stop("If g = NULL, length(STATS) must match ncol(x)");
        switch(ret) {
        case 1:
          for(int j = col; j--; ) out(_,j) = rep(AG[j], l);
          break;
        case 2:
          for(int j = col; j--; ) {
            IntegerMatrix::Column colo = out( _ , j);
            IntegerMatrix::Column column = xx( _ , j);
            int sumj = AG[j];
            for(int i = l; i--; ) {
              if(column[i] == NA_INTEGER) colo[i] = NA_INTEGER;
              else colo[i] = sumj;
            }
          }
          break;
        default: stop("Unknown Transformation");
        }
      } else {
        IntegerMatrix AG = xAG;
        if(AG.ncol() != col) stop("ncol(STATS) must match ncol(x)");
        if(gs != l) stop("length(g) must match nrow(x)");
        switch(ret) {
        case 1:
          for(int j = col; j--; ) {
            IntegerMatrix::Column colo = out( _ , j);
            IntegerMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) colo[i] = sumj[g[i]-1];
          }
          break;
        case 2:
          for(int j = col; j--; ) {
            IntegerMatrix::Column colo = out( _ , j);
            IntegerMatrix::Column column = xx( _ , j);
            IntegerMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) {
              if(column[i] == NA_INTEGER) colo[i] = NA_INTEGER;
              else colo[i] = sumj[g[i]-1];
            }
          }
          break;
        default: stop("Unknown Transformation");
        }
      }
      DUPLICATE_ATTRIB(out, xx); // or x ?? which is faster ??
      return out;
    }
    case STRSXP:
    {
      CharacterMatrix xx = x;
      int l = xx.nrow(), col = xx.ncol();
      CharacterMatrix out = no_init_matrix(l, col);
      if(gs == 1) {
        CharacterVector AG = xAG;
        if(AG.size() != col) stop("If g = NULL, length(STATS) must match ncol(x)");
        switch(ret) {
        case 1:
          for(int j = col; j--; ) {
            CharacterMatrix::Column outj = out(_,j);
            std::fill(outj.begin(), outj.end(), AG[j]);
          }
          break;
        case 2:
          for(int j = col; j--; ) {
            CharacterMatrix::Column colo = out( _ , j);
            CharacterMatrix::Column column = xx( _ , j);
            String sumj = AG[j];
            for(int i = l; i--; ) {
              if(column[i] == NA_STRING) colo[i] = NA_STRING;
              else colo[i] = sumj;
            }
          }
          break;
        default: stop("Unknown Transformation");
        }
      } else {
        CharacterMatrix AG = xAG;
        if(AG.ncol() != col) stop("ncol(STATS) must match ncol(x)");
        if(gs != l) stop("length(g) must match nrow(x)");
        switch(ret) {
        case 1:
          for(int j = col; j--; ) {
            CharacterMatrix::Column colo = out( _ , j);
            CharacterMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) colo[i] = sumj[g[i]-1];
          }
          break;
        case 2:
          for(int j = col; j--; ) {
            CharacterMatrix::Column colo = out( _ , j);
            CharacterMatrix::Column column = xx( _ , j);
            CharacterMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) {
              if(column[i] == NA_STRING) colo[i] = NA_STRING;
              else colo[i] = sumj[g[i]-1];
            }
          }
          break;
        default: stop("Unknown Transformation");
        }
      }
      DUPLICATE_ATTRIB(out, xx); // or x ?? which is faster ??
      return out;
    }
    case LGLSXP:
    {
      LogicalMatrix xx = x;
      int l = xx.nrow(), col = xx.ncol();
      LogicalMatrix out = no_init_matrix(l, col);
      if(gs == 1) {
        LogicalVector AG = xAG;
        if(AG.size() != col) stop("length(STATS) must match ncol(x)");
        switch(ret) {
        case 1:
          for(int j = col; j--; ) {
            LogicalMatrix::Column outj = out(_,j);
            std::fill(outj.begin(), outj.end(), AG[j]);
          }
          break;
        case 2:
          for(int j = col; j--; ) {
            LogicalMatrix::Column colo = out( _ , j);
            LogicalMatrix::Column column = xx( _ , j);
            bool sumj = AG[j];
            for(int i = l; i--; ) {
              if(column[i] == NA_LOGICAL) colo[i] = NA_LOGICAL;
              else colo[i] = sumj;
            }
          }
          break;
        default: stop("Unknown Transformation");
        }
      } else {
        LogicalMatrix AG = xAG;
        if(AG.ncol() != col) stop("ncol(STATS) must match ncol(x)");
        if(gs != l) stop("length(g) must match nrow(x)");
        switch(ret) {
        case 1:
          for(int j = col; j--; ) {
            LogicalMatrix::Column colo = out( _ , j);
            LogicalMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) colo[i] = sumj[g[i]-1];
          }
          break;
        case 2:
          for(int j = col; j--; ) {
            LogicalMatrix::Column colo = out( _ , j);
            LogicalMatrix::Column column = xx( _ , j);
            LogicalMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) {
              if(column[i] == NA_LOGICAL) colo[i] = NA_LOGICAL;
              else colo[i] = sumj[g[i]-1];
            }
          }
          break;
        default: stop("Unknown Transformation");
        }
      }
      DUPLICATE_ATTRIB(out, xx); // or x ?? which is faster ??
      return out;
    }
    default:
      stop("Not supported SEXP type!");
    }
  } else {
    switch(TYPEOF(x)) {
    case REALSXP:
    case INTSXP:
    {
      NumericMatrix xx = x;
      int l = xx.nrow(), col = xx.ncol();
      NumericMatrix out = no_init_matrix(l, col);
      if(gs == 1) {
        NumericVector AG = xAG;
        if(AG.size() != col) stop("length(STATS) must match ncol(x)");
        switch(ret) {
        case 3:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            colo = column - AG[j];
          }
          break;
        case 4: stop("This transformation can only be computed with groups!");
        case 5:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            colo = column * (1/AG[j]); // faster ??
          }
          break;
        case 6:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            colo = column * (100/AG[j]);
          }
          break;
        case 7:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            colo = column + AG[j];
          }
          break;
        case 8:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            colo = column * AG[j];
          }
          break;
        default: stop("Unknown Transformation");
        }
      } else {
        NumericMatrix AG = xAG;
        if(AG.ncol() != col) stop("ncol(STATS) must match ncol(x)");
        if(gs != l) stop("length(g) must match nrow(x)");
        switch(ret) {
        case 3:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            NumericMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) colo[i] = column[i] - sumj[g[i]-1];
          }
          break;
        case 4:
          { // Needed ??
            for(int j = col; j--; ) {
            NumericMatrix::Column column = xx( _ , j);
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column sumj = AG( _ , j);
            double OM = 0;
            int n = 0;
            for(int i = l; i--; ) { // Faster way ??
              if(std::isnan(column[i])) colo[i] = column[i];
              else { // Problem: if one sumj remained NA, oM becomes NA !!!
                colo[i] = column[i] - sumj[g[i]-1];
                if(std::isnan(sumj[g[i]-1])) continue; // solves the issue !!
                OM += sumj[g[i]-1]; // column[i]; good ??
                ++n;
              }
            }
            OM = OM / n;
            colo = colo + OM; // Fastest ??
          }
            break;
          }
        case 5:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            NumericMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) colo[i] = column[i] * (1/sumj[g[i]-1]); // fastest ??
          }
          break;
        case 6:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            NumericMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) colo[i] = column[i] * (100/sumj[g[i]-1]);
          }
          break;
        case 7:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            NumericMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) colo[i] = column[i] + sumj[g[i]-1];
          }
          break;
        case 8:
          for(int j = col; j--; ) {
            NumericMatrix::Column colo = out( _ , j);
            NumericMatrix::Column column = xx( _ , j);
            NumericMatrix::Column sumj = AG( _ , j);
            for(int i = l; i--; ) colo[i] = column[i] * sumj[g[i]-1];
          }
          break;
        default: stop("Unknown Transformation");
        }
      }
      DUPLICATE_ATTRIB(out, xx); // or x ?? which is faster ??
      return out;
    }
    case STRSXP: stop("The requested transformation is not possible with strings");
    case LGLSXP: stop("The requested transformation is not possible with logical data");
    default:
      stop("Not supported SEXP type!");
    }
  }
}

// Old Version: Only Numeric Data
// // [[Rcpp::export]]
// NumericMatrix TRAmCpp(NumericMatrix x, SEXP xAG, IntegerVector g = 0, int ret = 0) {
//   int l = x.nrow(), col = x.ncol(), gs = g.size();
//   NumericMatrix out = no_init_matrix(l, col);
//
//   if (gs == 1) { // ng redundant !! -> still do speed check!!
//     NumericVector AG = xAG; // initialize better ?? -> Nope, good !!
//     if(AG.size() != col) stop("length(AG) must match ncol(x)");
//     switch(ret) {
//       case 1:
//         for(int j = col; j--; ) out(_,j) = rep(AG[j], l);
//         break;
//       case 2: // works without extra brackets ?? -> Same !!
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           double sumj = AG[j];
//           for(int i = l; i--; ) {
//             if(std::isnan(column[i])) colo[i] = column[i];
//             else colo[i] = sumj;
//           }
//         }
//         break;
//       case 3:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           colo = column - AG[j];
//         }
//         break;
//       case 4: return x;
//       case 5:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           colo = column / AG[j];
//         }
//         break;
//       case 6:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           colo = column * (100/AG[j]);
//         }
//         break;
//       case 7:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           colo = column + AG[j];
//         }
//         break;
//       case 8:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           colo = column * AG[j];
//         }
//         break;
//     }
//
//   } else { // With groups !!!
//     if(gs != l) stop("g must match nrow(x)");
//     NumericMatrix AG = xAG; // initialize better ??
//     if(AG.ncol() != col) stop("ncol(AG) must match ncol(x)");
//
//     switch(ret) {
//     case 1: // Possible without brackets ?? -> Yes!!
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column sumj = AG( _ , j);
//           for(int i = l; i--; ) colo[i] = sumj[g[i]-1];
//         }
//         break;
//       case 2:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column sumj = AG( _ , j);
//           for(int i = l; i--; ) {
//             if(std::isnan(column[i])) colo[i] = column[i];
//             else colo[i] = sumj[g[i]-1];
//           }
//         }
//         break;
//       case 3:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column sumj = AG( _ , j);
//           for(int i = l; i--; ) colo[i] = column[i] - sumj[g[i]-1];
//         }
//         break;
//       case 4:
//       { // Needed ??
//           for(int j = col; j--; ) {
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column sumj = AG( _ , j);
//           double OM = 0;
//           int n = 0;
//           for(int i = l; i--; ) { // Faster way ??
//             if(std::isnan(column[i])) colo[i] = column[i];
//             else {
//               colo[i] = column[i] - sumj[g[i]-1];
//               OM += sumj[g[i]-1]; // column[i]; good ??
//               ++n;
//             }
//           }
//           OM = OM / n;
//           colo = colo + OM; // Fastest ??
//         }
//           break;
//       }
//       case 5:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column sumj = AG( _ , j);
//           for(int i = l; i--; ) colo[i] = column[i] / sumj[g[i]-1];
//         }
//         break;
//       case 6:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column sumj = AG( _ , j);
//           for(int i = l; i--; ) colo[i] = column[i] * (100/sumj[g[i]-1]);
//         }
//         break;
//       case 7:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column sumj = AG( _ , j);
//           for(int i = l; i--; ) colo[i] = column[i] + sumj[g[i]-1];
//         }
//         break;
//       case 8:
//         for(int j = col; j--; ) {
//           NumericMatrix::Column colo = out( _ , j);
//           NumericMatrix::Column column = x( _ , j);
//           NumericMatrix::Column sumj = AG( _ , j);
//           for(int i = l; i--; ) colo[i] = column[i] * sumj[g[i]-1];
//         }
//         break;
//     }
//   }
//   out.attr("dimnames") = x.attr("dimnames"); // NEW, faster ?? -> Not working !!
//   return out;
// }
