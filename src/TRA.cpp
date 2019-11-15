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

// Todo: Checks

// [[Rcpp::export]]
SEXP TRACpp(const SEXP& x, const SEXP& xAG, const IntegerVector& g = 0, int ret = 1) {
  int gs = g.size();
  if(ret <= 2) {
    switch(TYPEOF(x)) {
    case REALSXP:
    {
      NumericVector xx = x;
      NumericVector AG = xAG;
      int l = xx.size();
      NumericVector out = no_init_vector(l);
      if(gs == 1) {
        if(AG.size() != 1) stop("If g = NULL, STATS needs to be an atomic element!");
        double AGx = AG[0];
        switch(ret) {
        case 1:
          std::fill(out.begin(), out.end(), AGx);
          break;
        case 2:
          for(int i = l; i--; ) {
            if(std::isnan(xx[i])) out[i] = xx[i];
            else out[i] = AGx;
          }
          break;
        default: stop("Unknown Transformation");
        }
      } else {
        if(gs != l) stop("g must match nrow(x)");
        switch(ret) {
        case 1:
          for(int i = l; i--; ) out[i] = AG[g[i]-1];
          break;
        case 2:
          for(int i = l; i--; ) {
            if(std::isnan(xx[i])) out[i] = xx[i];
            else out[i] = AG[g[i]-1];
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
      IntegerVector xx = x;
      IntegerVector AG = xAG;
      int l = xx.size();
      IntegerVector out = no_init_vector(l);
      if(gs == 1) {
        if(AG.size() != 1) stop("If g = NULL, STATS needs to be an atomic element!");
        int AGx = AG[0];
        switch(ret) {
        case 1:
          std::fill(out.begin(), out.end(), AGx);
          break;
        case 2:
          for(int i = l; i--; ) {
            if(xx[i] == NA_INTEGER) out[i] = xx[i];
            else out[i] = AGx;
          }
          break;
        default: stop("Unknown Transformation");
        }
      } else {
        if(gs != l) stop("g must match nrow(x)");
        switch(ret) {
        case 1:
          for(int i = l; i--; ) out[i] = AG[g[i]-1];
          break;
        case 2:
          for(int i = l; i--; ) {
            if(xx[i] == NA_INTEGER) out[i] = xx[i];
            else out[i] = AG[g[i]-1];
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
      CharacterVector xx = x;
      CharacterVector AG = xAG;
      int l = xx.size();
      CharacterVector out = no_init_vector(l);
      // if(ret > 2) stop("The requested transformation is not possible with strings");
      if(gs == 1) {
        if(AG.size() != 1) stop("If g = NULL, STATS needs to be an atomic element!");
        String AGx = AG[0];
        switch(ret) {
        case 1:
          std::fill(out.begin(), out.end(), AGx);
          break;
        case 2:
          for(int i = l; i--; ) {
            if(xx[i] == NA_STRING) out[i] = xx[i];
            else out[i] = AGx;
          }
          break;
        default: stop("Unknown Transformation");
        }
      } else {
        if(gs != l) stop("g must match nrow(x)");
        switch(ret) {
        case 1:
          for(int i = l; i--; ) out[i] = AG[g[i]-1];
          break;
        case 2:
          for(int i = l; i--; ) {
            if(xx[i] == NA_STRING) out[i] = xx[i];
            else out[i] = AG[g[i]-1];
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
      LogicalVector xx = x;
      LogicalVector AG = xAG;
      int l = xx.size();
      LogicalVector out = no_init_vector(l);
      // if(ret > 2) stop("The requested transformation is not possible with strings");
      if(gs == 1) {
        if(AG.size() != 1) stop("If g = NULL, STATS needs to be an atomic element!");
        bool AGx = AG[0];
        switch(ret) {
        case 1:
          std::fill(out.begin(), out.end(), AGx);
          break;
        case 2:
          for(int i = l; i--; ) {
            if(xx[i] == NA_LOGICAL) out[i] = xx[i];
            else out[i] = AGx;
          }
          break;
        default: stop("Unknown Transformation");
        }
      } else {
        if(gs != l) stop("g must match nrow(x)");
        switch(ret) {
        case 1:
          for(int i = l; i--; ) out[i] = AG[g[i]-1];
          break;
        case 2:
          for(int i = l; i--; ) {
            if(xx[i] == NA_LOGICAL) out[i] = xx[i];
            else out[i] = AG[g[i]-1];
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
      NumericVector xx = x;
      NumericVector AG = xAG;
      int l = xx.size();
      NumericVector out = no_init_vector(l);
      if(gs == 1) {
        if(AG.size() != 1) stop("If g = NULL, STATS needs to be an atomic element!");
        double AGx = AG[0];
        switch(ret) {
        case 3:
          out = xx - AGx;
          break;
        case 4: stop("This transformation can only be performed with groups!");
        case 5:
          out = xx / AGx;
          break;
        case 6:
          out = xx * (100 / AGx);
          break;
        case 7:
          out = xx + AGx;
          break;
        case 8:
          out = xx * AGx;
          break;
        default: stop("Unknown Transformation");
        }
      } else {
        if(gs != l) stop("g must match nrow(x)");
        switch(ret) {
        case 3:
          for(int i = l; i--; ) out[i] = xx[i] - AG[g[i]-1];
          break;
        case 4:
          {
            double OM = 0;
            int n = 0;
            for(int i = l; i--; ) { // Faster way ??
              if(std::isnan(xx[i])) out[i] = xx[i];
              else { // Problem: if one AG remained NA, oM becomes NA !!!
                out[i] = xx[i] - AG[g[i]-1];
                if(std::isnan(AG[g[i]-1])) continue; // solves the issue !!
                OM += AG[g[i]-1]; // x[i]; // we want the overall average stat, not x
                ++n;
              }
            }
            OM = OM / n;
            out = out + OM; // Fastest ??
            break;
          }
        case 5:
          for(int i = l; i--; ) out[i] = xx[i] / AG[g[i]-1];
          break;
        case 6:
          for(int i = l; i--; ) out[i] = xx[i] * (100/AG[g[i]-1]);
          break;
        case 7:
          for(int i = l; i--; ) out[i] = xx[i] + AG[g[i]-1];
          break;
        case 8:
          for(int i = l; i--; ) out[i] = xx[i] * AG[g[i]-1];
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

// Previous Version: Only for Numeric Vectors
// // [[Rcpp::export]]
// NumericVector TRACpp(NumericVector x, NumericVector xAG, IntegerVector g = 0, int ret = 0) {
//   int l = x.size(), gs = g.size();
//
//   if (gs == 1) { // ng redundant !! -> still do speed check!!
//     double AG = xAG[0]; // initialize better ?? -> Nope, good !!
//     switch(ret) {
//       case 1: return rep(AG, l);
//       case 2:
//       {
//         NumericVector out = no_init_vector(l); // Initialize ??
//         for(int i = l; i--; ) {
//           if(std::isnan(x[i])) out[i] = x[i];
//           else out[i] = AG;
//         }
//         return out;
//       }
//       case 3: return x - AG;
//       case 4: return x;
//       case 5: return x / AG;
//       case 6: return x * (100/AG);
//       case 7: return x + AG;
//       case 8: return x * AG;
//     }
//   } else {
//     if(gs != l) stop("g must match nrow(x)");
//     NumericVector AG = xAG; // initialize better ??
//     NumericVector out = no_init_vector(l); // Initialize ??
//     switch(ret) {
//       case 1:
//         for(int i = l; i--; ) out[i] = AG[g[i]-1];
//         return out;
//       case 2:
//         for(int i = l; i--; ) {
//           if(std::isnan(x[i])) out[i] = x[i];
//           else out[i] = AG[g[i]-1];
//         }
//         return out;
//       case 3:
//         for(int i = l; i--; ) out[i] = x[i] - AG[g[i]-1];
//         return out;
//       case 4:
//       {
//         double OM = 0;
//         int n = 0;
//         for(int i = l; i--; ) { // Faster way ??
//           if(std::isnan(x[i])) out[i] = x[i];
//           else {
//             out[i] = x[i] - AG[g[i]-1];
//             OM += AG[g[i]-1]; // x[i]; // we want the overall average stat, not x
//             ++n;
//           }
//         }
//         OM = OM / n;
//         out = out + OM; // Fastest ??
//         return out;
//       }
//       case 5:
//         for(int i = l; i--; ) out[i] = x[i] / AG[g[i]-1];
//         return out;
//       case 6:
//         for(int i = l; i--; ) out[i] = x[i] * (100/AG[g[i]-1]);
//         return out;
//       case 7:
//         for(int i = l; i--; ) out[i] = x[i] + AG[g[i]-1];
//         return out;
//       case 8:
//         for(int i = l; i--; ) out[i] = x[i] * AG[g[i]-1];
//         return out;
//     }
//   }
//   return R_NilValue;
// }
