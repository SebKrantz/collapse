// [[Rcpp::plugins(cpp11)]]
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

// Todo: Make list input working if g = 0 !!

// [[Rcpp::export]]
List TRAlCpp(const List& x, const SEXP& xAG, const IntegerVector& g = 0, int ret = 1) {
  int l = x.size(), gs = g.size();
  List out(l);

  if(gs == 1) { // ng redundant !! -> still do speed check!!
    if(ret <= 2) {
      switch(TYPEOF(xAG)) {
      case VECSXP: {
        List AG = xAG;
        if(AG.size() != l) stop("length(STATS) must match length(x)");
        switch(ret) {
        case 1: {
          for(int j = l; j--; ) {
          switch(TYPEOF(x[j])) {
          case REALSXP: {
            NumericVector column = x[j];
            NumericVector outj(column.size(), AG[j]);
            SHALLOW_DUPLICATE_ATTRIB(outj, column); // Here or before filling ??
            out[j] = outj;
            break;
          }
          case INTSXP: {
            IntegerVector column = x[j];
            IntegerVector outj(column.size(), AG[j]);
            SHALLOW_DUPLICATE_ATTRIB(outj, column); // Here or before filling ??
            out[j] = outj;
            break;
          }
          case STRSXP: {
            CharacterVector column = x[j];
            CharacterVector outj(column.size(), wrap(AG[j]));
            SHALLOW_DUPLICATE_ATTRIB(outj, column); // Here or before filling ??
            out[j] = outj;
            break;
          }
          case LGLSXP: {
            LogicalVector column = x[j];
            LogicalVector outj(column.size(), AG[j]);
            SHALLOW_DUPLICATE_ATTRIB(outj, column); // Here or before filling ??
            out[j] = outj;
            break;
          }
          default: stop("list x element of unsupported type");
          }
        }
          break;
        }
        case 2: {
          for(int j = l; j--; ) {
          switch(TYPEOF(x[j])) {
          case REALSXP: {
            NumericVector column = x[j];
            int row = column.size();
            NumericVector sgj = no_init_vector(row);
            double sumj = AG[j];
            for(int i = row; i--; ) {
              if(std::isnan(column[i])) sgj[i] = column[i];
              else sgj[i] = sumj;
            }
            SHALLOW_DUPLICATE_ATTRIB(sgj, column); // Here or before filling ??
            out[j] = sgj;
            break;
          }
          case INTSXP: {
            IntegerVector column = x[j];
            int row = column.size();
            IntegerVector sgj = no_init_vector(row);
            int sumj = AG[j];
            for(int i = row; i--; ) {
              if(column[i] == NA_INTEGER) sgj[i] = NA_INTEGER;
              else sgj[i] = sumj;
            }
            SHALLOW_DUPLICATE_ATTRIB(sgj, column); // Here or before filling ??
            out[j] = sgj;
            break;
          }
          case STRSXP: {
            CharacterVector column = x[j];
            int row = column.size();
            CharacterVector sgj = no_init_vector(row);
            String sumj = AG[j];
            for(int i = row; i--; ) {
              if(column[i] == NA_STRING) sgj[i] = NA_STRING;
              else sgj[i] = sumj;
            }
            SHALLOW_DUPLICATE_ATTRIB(sgj, column); // Here or before filling ??
            out[j] = sgj;
            break;
          }
          case LGLSXP: {
            LogicalVector column = x[j];
            int row = column.size();
            LogicalVector sgj = no_init_vector(row);
            bool sumj = AG[j];
            for(int i = row; i--; ) {
              if(column[i] == NA_LOGICAL) sgj[i] = NA_LOGICAL;
              else sgj[i] = sumj;
            }
            SHALLOW_DUPLICATE_ATTRIB(sgj, column); // Here or before filling ??
            out[j] = sgj;
            break;
          }
          default: stop("list x element of unsupported type");
          }
        }
          break;
        }
        }
        break;
      }
      case REALSXP: {
        NumericVector AG = xAG;
        if(AG.size() != l) stop("length(STATS) must match length(x)");
        switch(ret) {
        case 1: {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            NumericVector outj(column.size(), AG[j]);
            SHALLOW_DUPLICATE_ATTRIB(outj, column); // Here or before filling ??
            out[j] = outj;
          }
          break;
        }
        case 2: {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            int row = column.size();
            NumericVector sgj = no_init_vector(row);
            double sumj = AG[j];
            for(int i = row; i--; ) {
              if(std::isnan(column[i])) sgj[i] = column[i];
              else sgj[i] = sumj;
          }
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // Here or before filling ??
          out[j] = sgj;
        }
          break;
        }
        }
        break;
      }
      case INTSXP: {
        IntegerVector AG = xAG;
        if(AG.size() != l) stop("length(STATS) must match length(x)");
        switch(ret) {
        case 1: {
          for(int j = l; j--; ) {
          IntegerVector column = x[j];
          IntegerVector outj(column.size(), AG[j]);
          SHALLOW_DUPLICATE_ATTRIB(outj, column); // Here or before filling ??
          out[j] = outj;
        }
          break;
        }
        case 2: {
          for(int j = l; j--; ) {
          IntegerVector column = x[j];
          int row = column.size();
          IntegerVector sgj = no_init_vector(row);
          int sumj = AG[j];
          for(int i = row; i--; ) {
            if(column[i] == NA_INTEGER) sgj[i] = NA_INTEGER;
            else sgj[i] = sumj;
          }
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // Here or before filling ??
          out[j] = sgj;
        }
          break;
        }
        }
        break;
      }
      case STRSXP: {
        CharacterVector AG = xAG;
        if(AG.size() != l) stop("length(STATS) must match length(x)");
        switch(ret) {
        case 1: {
          for(int j = l; j--; ) {
          CharacterVector column = x[j];
          CharacterVector outj(column.size(), AG[j]);
          SHALLOW_DUPLICATE_ATTRIB(outj, column); // Here or before filling ??
          out[j] = outj;
        }
          break;
        }
        case 2: {
          for(int j = l; j--; ) {
          CharacterVector column = x[j];
          int row = column.size();
          CharacterVector sgj = no_init_vector(row);
          String sumj = AG[j];
          for(int i = row; i--; ) {
            if(column[i] == NA_STRING) sgj[i] = NA_STRING;
            else sgj[i] = sumj;
          }
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // Here or before filling ??
          out[j] = sgj;
        }
          break;
        }
        }
        break;
      }
      case LGLSXP: {
        LogicalVector AG = xAG;
        if(AG.size() != l) stop("length(STATS) must match length(x)");
        switch(ret) {
        case 1: {
          for(int j = l; j--; ) {
          LogicalVector column = x[j];
          LogicalVector outj(column.size(), AG[j]);
          SHALLOW_DUPLICATE_ATTRIB(outj, column);
          out[j] = outj;
        }
          break;
        }
        case 2: {
          for(int j = l; j--; ) {
          LogicalVector column = x[j];
          int row = column.size();
          LogicalVector sgj = no_init_vector(row);
          bool sumj = AG[j];
          for(int i = row; i--; ) {
            if(column[i] == NA_LOGICAL) sgj[i] = NA_LOGICAL;
            else sgj[i] = sumj;
          }
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // Here or before filling ??
          out[j] = sgj;
        }
          break;
        }
        }
        break;
      }
      default: stop("Not supported SEXP type!");
      } // Faster way ?? Switch statements other way around ???
    } else {
      NumericVector AG = NULL;
      if(TYPEOF(xAG) == VECSXP) {
        AG = NumericVector(l);
        List temp = xAG;
        if(temp.size() != l) stop("length(STATS) must match length(x)");
        for(int i = l; i--; ) AG[i] = temp[i];
      } else {
        AG = xAG;
        if(AG.size() != l) stop("length(STATS) must match length(x)");
      }
       // Works for Lists ??
      switch(ret) {
      case 3: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        NumericVector outj = column - AG[j];
        SHALLOW_DUPLICATE_ATTRIB(outj, column);
        out[j] = outj;
      }
        break;
      }
      case 4: stop("This transformation can only be performed with groups!");
      case 5: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        NumericVector outj = column * (1/AG[j]);
        SHALLOW_DUPLICATE_ATTRIB(outj, column);
        out[j] = outj;
      }
        break;
      }
      case 6: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        NumericVector outj = column * (100/AG[j]);
        SHALLOW_DUPLICATE_ATTRIB(outj, column);
        out[j] = outj;
      }
        break;
      }
      case 7: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        NumericVector outj = column + AG[j];
        SHALLOW_DUPLICATE_ATTRIB(outj, column);
        out[j] = outj;
      }
        break;
      }
      case 8: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        NumericVector outj = column * AG[j];
        SHALLOW_DUPLICATE_ATTRIB(outj, column);
        out[j] = outj;
      }
        break;
      }
      default: stop("Unknown Transformation");
      }
    }
  } else {
    List AG = xAG; // initialize better ??
    if(AG.size() != l) stop("length(STATS) must match length(x)");
    if(ret <= 2) {
      switch(ret) {
      case 1: {
      for(int j = l; j--; ) {
       switch(TYPEOF(AG[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          if(column.size() != gs) stop("length(g) must match nrow(x)");
          NumericVector sgj = no_init_vector(gs);
          NumericVector sumj = AG[j];
          for(int i = gs; i--; ) sgj[i] = sumj[g[i]-1];
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
          out[j] = sgj;
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          if(column.size() != gs) stop("length(g) must match nrow(x)");
          IntegerVector sgj = no_init_vector(gs);
          IntegerVector sumj = AG[j];
          for(int i = gs; i--; ) sgj[i] = sumj[g[i]-1];
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
          out[j] = sgj;
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          if(column.size() != gs) stop("length(g) must match nrow(x)");
          CharacterVector sgj = no_init_vector(gs);
          CharacterVector sumj = AG[j];
          for(int i = gs; i--; ) sgj[i] = sumj[g[i]-1];
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
          out[j] = sgj;
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          if(column.size() != gs) stop("length(g) must match nrow(x)");
          LogicalVector sgj = no_init_vector(gs);
          LogicalVector sumj = AG[j];
          for(int i = gs; i--; ) sgj[i] = sumj[g[i]-1];
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
          out[j] = sgj;
          break;
        }
        default: stop("list x element of unsupported type");
       }
      }
      break;
      }
      case 2: {
        for(int j = l; j--; ) {
        switch(TYPEOF(AG[j])) {
        case REALSXP: {
          NumericVector column = x[j];
          if(column.size() != gs) stop("length(g) must match nrow(x)");
          NumericVector sgj = no_init_vector(gs);
          NumericVector sumj = AG[j];
          for(int i = gs; i--; ) {
            if(std::isnan(column[i])) sgj[i] = column[i];
            else sgj[i] = sumj[g[i]-1];
          }
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
          out[j] = sgj;
          break;
        }
        case INTSXP: {
          IntegerVector column = x[j];
          if(column.size() != gs) stop("length(g) must match nrow(x)");
          IntegerVector sgj = no_init_vector(gs);
          IntegerVector sumj = AG[j];
          for(int i = gs; i--; ) {
            if(column[i] == NA_INTEGER) sgj[i] = NA_INTEGER;
            else sgj[i] = sumj[g[i]-1];
          }
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
          out[j] = sgj;
          break;
        }
        case STRSXP: {
          CharacterVector column = x[j];
          if(column.size() != gs) stop("length(g) must match nrow(x)");
          CharacterVector sgj = no_init_vector(gs);
          CharacterVector sumj = AG[j];
          for(int i = gs; i--; ) {
            if(column[i] == NA_STRING) sgj[i] = NA_STRING;
            else sgj[i] = sumj[g[i]-1];
          }
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
          out[j] = sgj;
          break;
        }
        case LGLSXP: {
          LogicalVector column = x[j];
          if(column.size() != gs) stop("length(g) must match nrow(x)");
          LogicalVector sgj = no_init_vector(gs);
          LogicalVector sumj = AG[j];
          for(int i = gs; i--; ) {
            if(column[i] == NA_LOGICAL) sgj[i] = NA_LOGICAL;
            else sgj[i] = sumj[g[i]-1];
          }
          SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
          out[j] = sgj;
          break;
        }
        default: stop("list x element of unsupported type");
        }
        }
      break;
      }
      }
    } else {
      switch(ret) {
      case 3: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(column.size() != gs) stop("length(g) must match nrow(x)");
        NumericVector sgj = no_init_vector(gs);
        NumericVector sumj = AG[j];
        for(int i = gs; i--; ) sgj[i] = column[i] - sumj[g[i]-1];
        SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
        out[j] = sgj;
      }
        break;
      }
      case 4: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(column.size() != gs) stop("length(g) must match nrow(x)");
        NumericVector sgj = no_init_vector(gs);
        NumericVector sumj = AG[j];
        double OM = 0;
        int n = 0;
        for(int i = gs; i--; ) { // Faster way ??
          if(std::isnan(column[i])) sgj[i] = column[i];
          else {
            sgj[i] = column[i] - sumj[g[i]-1];
            if(std::isnan(sumj[g[i]-1])) continue;
            OM += sumj[g[i]-1]; // column[i]; // good??
            ++n;
          }
        }
        OM = OM / n;
        sgj = sgj + OM; // Fastest !!
        SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
        out[j] = sgj;
      }
        break;
      }
      case 5: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(column.size() != gs) stop("length(g) must match nrow(x)");
        NumericVector sgj = no_init_vector(gs);
        NumericVector sumj = AG[j];
        for(int i = gs; i--; ) sgj[i] = column[i] * (1/sumj[g[i]-1]); // fastest ??
        SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
        out[j] = sgj;
      }
        break;
      }
      case 6: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(column.size() != gs) stop("length(g) must match nrow(x)");
        NumericVector sgj = no_init_vector(gs);
        NumericVector sumj = AG[j];
        for(int i = gs; i--; ) sgj[i] = column[i] * (100 / sumj[g[i]-1]);
        SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
        out[j] = sgj;
      }
        break;
      }
      case 7: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(column.size() != gs) stop("length(g) must match nrow(x)");
        NumericVector sgj = no_init_vector(gs);
        NumericVector sumj = AG[j];
        for(int i = gs; i--; ) sgj[i] = column[i] + sumj[g[i]-1];
        SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
        out[j] = sgj;
      }
        break;
      }
      case 8: {
        for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(column.size() != gs) stop("length(g) must match nrow(x)");
        NumericVector sgj = no_init_vector(gs);
        NumericVector sumj = AG[j];
        for(int i = gs; i--; ) sgj[i] = column[i] * sumj[g[i]-1];
        SHALLOW_DUPLICATE_ATTRIB(sgj, column); // here or before filling??
        out[j] = sgj;
      }
        break;
      }
      default: stop("Unknown transformation!");
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  return out;
}

//Draft: Too complex !!
// // [[Rcpp::export]]
// List TRAlCpp(List x, SEXP xAG, IntegerVector g = 0, int ret = 1) {
//   int l = x.size(), row = 0, gs = g.size();
//   List out(l);
//
//   if (gs == 1) { // ng redundant !! -> still do speed check!!
//     switch(TYPEOF(xAG)) {
//     case VECSXP: {
//       List AG = xAG;
//       if(AG.size() != l) stop("length(AG) must match length(x)");
//       switch(ret) {
//       case 1: {
//         for(int i = l; i--; ) {
//         switch(TYPEOF(x[i])) {
//         case REALSXP: {
//           NumericVector column = x[i];
//           row = column.size();
//           NumericVector AGi = AG[i];
//           out[i] = rep(AGi[0], row);
//           break;
//         }
//         case INTSXP: {
//           IntegerVector column = x[i];
//           row = column.size();
//           IntegerVector AGi = AG[i];
//           out[i] = rep(AGi[0], row);
//           break;
//         }
//         case STRSXP: {
//           CharacterVector column = x[i];
//           row = column.size();
//           CharacterVector AGi = AG[i];
//           out[i] = rep(AGi[0], row);
//           break;
//         }
//         case LGLSXP: {
//           break;
//         }
//         default: stop("list x element of unsupported type");
//         }
//       }
//         break;
//       }
//       case 2: {
//         for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) {
//         case REALSXP: {
//           break;
//         }
//         case INTSXP: {
//           break;
//         }
//         case STRSXP: {
//           break;
//         }
//         case LGLSXP: {
//           break;
//         }
//         default: stop("list x element of unsupported type");
//         }
//         NumericVector column = x[j];
//         row = column.size();
//         NumericVector sgj = no_init_vector(row);
//         double sumj = AG[j];
//         for(int i = row; i--; ) {
//           if(std::isnan(column[i])) sgj[i] = column[i];
//           else sgj[i] = sumj;
//         }
//         out[j] = sgj;
//       }
//         break;
//       }
//       case 3: {
//         for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) {
//         case REALSXP: {
//           break;
//         }
//         case INTSXP: {
//           break;
//         }
//         case STRSXP: {
//           break;
//         }
//         case LGLSXP: {
//           break;
//         }
//         default: stop("list x element of unsupported type");
//         }
//         NumericVector column = x[j];
//         out[j] = column - AG[j];
//       }
//         break;
//       }
//       case 4: return x;
//       case 5: {
//         for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) {
//         case REALSXP: {
//           break;
//         }
//         case INTSXP: {
//           break;
//         }
//         case STRSXP: {
//           break;
//         }
//         case LGLSXP: {
//           break;
//         }
//         default: stop("list x element of unsupported type");
//         }
//         NumericVector column = x[j];
//         out[j] = column / AG[j];
//       }
//         break;
//       }
//       case 6: {
//         for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) {
//         case REALSXP: {
//           break;
//         }
//         case INTSXP: {
//           break;
//         }
//         case STRSXP: {
//           break;
//         }
//         case LGLSXP: {
//           break;
//         }
//         default: stop("list x element of unsupported type");
//         }
//         NumericVector column = x[j];
//         out[j] = column * (100/AG[j]);
//       }
//         break;
//       }
//       case 7: {
//         for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) {
//         case REALSXP: {
//           break;
//         }
//         case INTSXP: {
//           break;
//         }
//         case STRSXP: {
//           break;
//         }
//         case LGLSXP: {
//           break;
//         }
//         default: stop("list x element of unsupported type");
//         }
//         NumericVector column = x[j];
//         out[j] = column + AG[j];
//       }
//         break;
//       }
//       case 8: {
//         for(int j = l; j--; ) {
//         switch(TYPEOF(x[j])) {
//         case REALSXP: {
//           break;
//         }
//         case INTSXP: {
//           break;
//         }
//         case STRSXP: {
//           break;
//         }
//         case LGLSXP: {
//           break;
//         }
//         default: stop("list x element of unsupported type");
//         }
//         NumericVector column = x[j];
//         out[j] = column * AG[j];
//       }
//         break;
//       }
//       }
//       break;
//     }
//     case REALSXP: {
//       NumericVector AG = xAG;
//       if(AG.size() != l) stop("length(AG) must match length(x)");
//       switch(ret) {
//       case 1: {
//         for(int i = l; i--; ) {
//           NumericVector column = x[i];
//           row = column.size();
//           out[i] = rep(AG[i], row);
//         }
//         break;
//       }
//       case 2: {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           row = column.size();
//           NumericVector sgj = no_init_vector(row);
//           double sumj = AG[j];
//           for(int i = row; i--; ) {
//             if(std::isnan(column[i])) sgj[i] = column[i];
//             else sgj[i] = sumj;
//         }
//         out[j] = sgj;
//       }
//         break;
//       }
//       case 3: {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           out[j] = column - AG[j];
//         }
//         break;
//       }
//       case 4: return x;
//       case 5: {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           out[j] = column / AG[j];
//         }
//         break;
//       }
//       case 6: {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           out[j] = column * (100/AG[j]);
//         }
//         break;
//       }
//       case 7: {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           out[j] = column + AG[j];
//         }
//         break;
//       }
//       case 8: {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           out[j] = column * AG[j];
//         }
//         break;
//       }
//       }
//       break;
//     }
//     case INTSXP: {
//       IntegerVector AG = xAG;
//       if(AG.size() != l) stop("length(AG) must match length(x)");
//       switch(ret) {
//       case 1: {
//         for(int i = l; i--; ) {
//         IntegerVector column = x[i];
//         row = column.size();
//         out[i] = rep(AG[i], row);
//       }
//         break;
//       }
//       case 2: {
//         for(int j = l; j--; ) {
//         IntegerVector column = x[j];
//         row = column.size();
//         IntegerVector sgj = no_init_vector(row);
//         int sumj = AG[j];
//         for(int i = row; i--; ) {
//           if(column[i] == NA_INTEGER) sgj[i] = column[i];
//           else sgj[i] = sumj;
//         }
//         out[j] = sgj;
//       }
//         break;
//       }
//       case 3: {
//         for(int j = l; j--; ) {
//         IntegerVector column = x[j];
//         out[j] = column - AG[j];
//       }
//         break;
//       }
//       case 4: return x;
//       case 5: {
//         for(int j = l; j--; ) {
//         IntegerVector column = x[j];
//         out[j] = column / AG[j];
//       }
//         break;
//       }
//       case 6: {
//         for(int j = l; j--; ) {
//         IntegerVector column = x[j];
//         out[j] = column * (100/AG[j]);
//       }
//         break;
//       }
//       case 7: {
//         for(int j = l; j--; ) {
//         IntegerVector column = x[j];
//         out[j] = column + AG[j];
//       }
//         break;
//       }
//       case 8: {
//         for(int j = l; j--; ) {
//         IntegerVector column = x[j];
//         out[j] = column * AG[j];
//       }
//         break;
//       }
//       }
//       break;
//     }
//     case STRSXP: {
//       CharacterVector AG = xAG;
//       if(AG.size() != l) stop("length(AG) must match length(x)");
//       switch(ret) {
//       case 1: {
//         for(int i = l; i--; ) {
//         CharacterVector column = x[i];
//         row = column.size();
//         CharacterVector outi = no_init_vector(row);
//         std::fill(outi.begin(), outi.end(), AG[i]);
//         out[i] = outi;
//       }
//         break;
//       }
//       case 2: {
//         for(int j = l; j--; ) {
//         CharacterVector column = x[j];
//         row = column.size();
//         CharacterVector sgj = no_init_vector(row);
//         String sumj = AG[j];
//         for(int i = row; i--; ) {
//           if(column[i] == NA_STRING) sgj[i] = column[i];
//           else sgj[i] = sumj;
//         }
//         out[j] = sgj;
//       }
//         break;
//       }
//       default: stop("The requested transformation is not possible with strings");
//       }
//       break;
//     }
//     case LGLSXP: {
//       LogicalVector AG = xAG;
//       if(AG.size() != l) stop("length(AG) must match length(x)");
//       switch(ret) {
//       case 1: {
//         for(int i = l; i--; ) {
//         LogicalVector column = x[i];
//         row = column.size();
//         LogicalVector outi = no_init_vector(row);
//         std::fill(outi.begin(), outi.end(), AG[i]);
//         out[i] = outi;
//       }
//         break;
//       }
//       case 2: {
//         for(int j = l; j--; ) {
//         LogicalVector column = x[j];
//         row = column.size();
//         LogicalVector sgj = no_init_vector(row);
//         bool sumj = AG[j];
//         for(int i = row; i--; ) {
//           if(column[i] == NA_LOGICAL) sgj[i] = column[i];
//           else sgj[i] = sumj;
//         }
//         out[j] = sgj;
//       }
//         break;
//       }
//       default: stop("The requested transformation is not possible with logical values");
//       }
//       break;
//     }
//     default: stop("Not supported SEXP type!");
//     }
//   } else {
//     List AG = xAG; // initialize better ??
//     if(AG.size() != l) stop("length(AG) must match length(x)");
//     switch(ret) {
//     case 1: {
//       for(int j = l; j--; ) {
//       NumericVector column = x[j];
//       row = column.size();
//       if(row != gs) stop("g must match nrow(x)");
//       NumericVector sgj = no_init_vector(row);
//       NumericVector sumj = AG[j];
//       for(int i = row; i--; ) sgj[i] = sumj[g[i]-1];
//       out[j] = sgj;
//     }
//       break;
//     }
//     case 2: {
//       for(int j = l; j--; ) {
//       NumericVector column = x[j];
//       row = column.size();
//       if(row != gs) stop("g must match nrow(x)");
//       NumericVector sgj = no_init_vector(row);
//       NumericVector sumj = AG[j];
//       for(int i = row; i--; ) {
//         if(std::isnan(column[i])) sgj[i] = column[i];
//         else sgj[i] = sumj[g[i]-1];
//       }
//       out[j] = sgj;
//     }
//       break;
//     }
//     case 3: {
//       for(int j = l; j--; ) {
//       NumericVector column = x[j];
//       row = column.size();
//       if(row != gs) stop("g must match nrow(x)");
//       NumericVector sgj = no_init_vector(row);
//       NumericVector sumj = AG[j];
//       for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1];
//       out[j] = sgj;
//     }
//       break;
//     }
//     case 4: {
//       for(int j = l; j--; ) {
//       NumericVector column = x[j];
//       row = column.size();
//       if(row != gs) stop("g must match nrow(x)");
//       NumericVector sgj = no_init_vector(row);
//       NumericVector sumj = AG[j];
//       double OM = 0;
//       int n = 0;
//       for(int i = row; i--; ) { // Faster way ??
//         if(std::isnan(column[i])) sgj[i] = column[i];
//         else {
//           sgj[i] = column[i] - sumj[g[i]-1];
//           OM += sumj[g[i]-1]; // column[i]; // good??
//           ++n;
//         }
//       }
//       OM = OM / n;
//       sgj = sgj + OM; // Fastest !!
//       out[j] = sgj;
//     }
//       break;
//     }
//     case 5: {
//       for(int j = l; j--; ) {
//       NumericVector column = x[j];
//       row = column.size();
//       if(row != gs) stop("g must match nrow(x)");
//       NumericVector sgj = no_init_vector(row);
//       NumericVector sumj = AG[j];
//       for(int i = row; i--; ) sgj[i] = column[i] / sumj[g[i]-1];
//       out[j] = sgj;
//     }
//       break;
//     }
//     case 6: {
//       for(int j = l; j--; ) {
//       NumericVector column = x[j];
//       row = column.size();
//       if(row != gs) stop("g must match nrow(x)");
//       NumericVector sgj = no_init_vector(row);
//       NumericVector sumj = AG[j];
//       for(int i = row; i--; ) sgj[i] = column[i] * (100 / sumj[g[i]-1]);
//       out[j] = sgj;
//     }
//       break;
//     }
//     case 7: {
//       for(int j = l; j--; ) {
//       NumericVector column = x[j];
//       row = column.size();
//       if(row != gs) stop("g must match nrow(x)");
//       NumericVector sgj = no_init_vector(row);
//       NumericVector sumj = AG[j];
//       for(int i = row; i--; ) sgj[i] = column[i] + sumj[g[i]-1];
//       out[j] = sgj;
//     }
//       break;
//     }
//     case 8: {
//       for(int j = l; j--; ) {
//       NumericVector column = x[j];
//       row = column.size();
//       if(row != gs) stop("g must match nrow(x)");
//       NumericVector sgj = no_init_vector(row);
//       NumericVector sumj = AG[j];
//       for(int i = row; i--; ) sgj[i] = column[i] * sumj[g[i]-1];
//       out[j] = sgj;
//     }
//       break;
//     }
//     }
//   }
//   return out;
// }


// Previous Version: Only Numeric Data
// // [[Rcpp::export]]
// List TRAlCpp(List x, SEXP xAG, IntegerVector g = 0, int ret = 0) {
//   int l = x.size(), row = 0, gs = g.size();
//   List out(l);
//
//   if (gs == 1) { // ng redundant !! -> still do speed check!!
//     NumericVector AG = xAG; // initialize better ?? -> Nope, good !!
//     if(AG.size() != l) stop("length(AG) must match length(x)");
//     switch(ret) {
//       case 1: {
//         for(int i = l; i--; ) {
//           NumericVector column = x[i];
//           row = column.size();
//           out[i] = rep(AG[i], row);
//         }
//         break;
//       }
//       case 2: {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           row = column.size();
//           NumericVector sgj = no_init_vector(row);
//             double sumj = AG[j];
//             for(int i = row; i--; ) {
//               if(std::isnan(column[i])) sgj[i] = column[i];
//               else sgj[i] = sumj;
//             }
//             out[j] = sgj;
//         }
//         break;
//       }
//       case 3: {
//         for(int j = l; j--; ) {
//             NumericVector column = x[j];
//             out[j] = column - AG[j];
//         }
//         break;
//       }
//       case 4: return x;
//       case 5: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         out[j] = column / AG[j];
//       }
//         break;
//       }
//       case 6: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         out[j] = column * (100/AG[j]);
//       }
//         break;
//       }
//       case 7: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         out[j] = column + AG[j];
//       }
//         break;
//       }
//       case 8: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         out[j] = column * AG[j];
//       }
//         break;
//       }
//     }
//   } else {
//     List AG = xAG; // initialize better ??
//     if(AG.size() != l) stop("length(AG) must match length(x)");
//     switch(ret) {
//       case 1: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(row != gs) stop("g must match nrow(x)");
//         NumericVector sgj = no_init_vector(row);
//         NumericVector sumj = AG[j];
//         for(int i = row; i--; ) sgj[i] = sumj[g[i]-1];
//         out[j] = sgj;
//         }
//         break;
//       }
//       case 2: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(row != gs) stop("g must match nrow(x)");
//         NumericVector sgj = no_init_vector(row);
//         NumericVector sumj = AG[j];
//         for(int i = row; i--; ) {
//           if(std::isnan(column[i])) sgj[i] = column[i];
//           else sgj[i] = sumj[g[i]-1];
//         }
//         out[j] = sgj;
//         }
//         break;
//       }
//       case 3: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(row != gs) stop("g must match nrow(x)");
//         NumericVector sgj = no_init_vector(row);
//         NumericVector sumj = AG[j];
//         for(int i = row; i--; ) sgj[i] = column[i] - sumj[g[i]-1];
//         out[j] = sgj;
//         }
//         break;
//       }
//       case 4: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(row != gs) stop("g must match nrow(x)");
//         NumericVector sgj = no_init_vector(row);
//         NumericVector sumj = AG[j];
//         double OM = 0;
//         int n = 0;
//         for(int i = row; i--; ) { // Faster way ??
//           if(std::isnan(column[i])) sgj[i] = column[i];
//           else {
//           sgj[i] = column[i] - sumj[g[i]-1];
//           OM += sumj[g[i]-1]; // column[i]; // good??
//           ++n;
//           }
//         }
//         OM = OM / n;
//         sgj = sgj + OM; // Fastest !!
//         out[j] = sgj;
//       }
//         break;
//       }
//       case 5: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(row != gs) stop("g must match nrow(x)");
//         NumericVector sgj = no_init_vector(row);
//         NumericVector sumj = AG[j];
//         for(int i = row; i--; ) sgj[i] = column[i] / sumj[g[i]-1];
//         out[j] = sgj;
//       }
//         break;
//       }
//       case 6: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(row != gs) stop("g must match nrow(x)");
//         NumericVector sgj = no_init_vector(row);
//         NumericVector sumj = AG[j];
//         for(int i = row; i--; ) sgj[i] = column[i] * (100 / sumj[g[i]-1]);
//         out[j] = sgj;
//       }
//         break;
//       }
//       case 7: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(row != gs) stop("g must match nrow(x)");
//         NumericVector sgj = no_init_vector(row);
//         NumericVector sumj = AG[j];
//         for(int i = row; i--; ) sgj[i] = column[i] + sumj[g[i]-1];
//         out[j] = sgj;
//       }
//         break;
//       }
//       case 8: {
//         for(int j = l; j--; ) {
//         NumericVector column = x[j];
//         row = column.size();
//         if(row != gs) stop("g must match nrow(x)");
//         NumericVector sgj = no_init_vector(row);
//         NumericVector sumj = AG[j];
//         for(int i = row; i--; ) sgj[i] = column[i] * sumj[g[i]-1];
//         out[j] = sgj;
//       }
//         break;
//       }
//     }
//   }
//   return out;
// }
