// // [[Rcpp::plugins(cpp11)]]
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
// 9- Modulus
// 10- Subtract Modulus

// int(x * (1/y)) -> This gave the UBSAN error if NaN !!!
inline double mymod(double x, double y) {
  double z(x * (1/y));
  return (z == z) ? x - (int)(z) * y : z; // faster than x - (int)(x/y) * y; // also C-style conversions seem to be faster ?
}

// #define mymod(x, y) (x - ((int)(x/y) * y)) // Macro: not faster !

// template<typename T, typename U>
// constexpr double mymod (T x, U mod)
// {
//   return !mod ? x : x - mod * static_cast<long long>(x / mod);
// }

// int(x * (1/y)) -> This gave the UBSAN error if NaN !!!
inline double myremain(double x, double y) {
  double z(x * (1/y));
  return (z == z) ? (int)(z) * y : z; //   (int)(x * (1/y)) * y; <- This would be enough, but doesn't keep missing values in x!
}


SEXP ret1(const SEXP& x, const SEXP& xAG, const SEXP& g) {
  int tx = TYPEOF(x), txAG = TYPEOF(xAG), l = Rf_length(x), gs = Rf_length(g);
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  int *pg;
  bool nog = gs == 1;
  if(nog) {
    if(Rf_length(xAG) != 1) stop("If g = NULL, NROW(STATS) needs to be 1");
  } else {
    if(gs != l) stop("length(g) must match NROW(x)");
    pg = INTEGER(g);
  }
  SEXP out = PROTECT(Rf_allocVector(txAG, l));

  switch(txAG) {
    case REALSXP:
    {
      double *pout = REAL(out);
      if(nog) { // memset(pout, Rf_asReal(xAG), l * sizeof(double)); memset only works with 0 !!
        double AG = Rf_asReal(xAG);
        for(int i = l; i--; ) pout[i] = AG;
      } else {
        double *AG = REAL(xAG)-1;
        for(int i = l; i--; ) pout[i] = AG[pg[i]];
      }
      break;
    }
    case INTSXP:
    {
      int *pout = INTEGER(out);
      if(nog) {
        int AG = Rf_asInteger(xAG);
        for(int i = l; i--; ) pout[i] = AG;
      } else {
        int *AG = INTEGER(xAG)-1;
        for(int i = l; i--; ) pout[i] = AG[pg[i]];
      }
      break;
    }
    case STRSXP:
    {
      // CharacterVector AG = xAG;
      // if(nog) out = CharacterVector(l, String(AG[0]));
      // else {
      //   CharacterVector pout = out;
      //   for(int i = l; i--; ) pout[i] = AG[pg[i]-1];
      // }
      // break;
      SEXP *pout = STRING_PTR(out);
      if(nog) {
        SEXP AG = Rf_asChar(xAG);
        for(int i = l; i--; ) pout[i] = AG; // SET_STRING_ELT(out, i, AG); // Without pointer -> much slower!
      } else {
        SEXP *AG = STRING_PTR(xAG)-1;
        for(int i = l; i--; ) pout[i] = AG[pg[i]]; // SET_STRING_ELT(out, i, AG[pg[i]]); // Without pointer -> much slower!
      }
      break;
    }
    case LGLSXP:
    {
      int *pout = LOGICAL(out);
      if(nog) {
        int AG = Rf_asLogical(xAG);
        for(int i = l; i--; ) pout[i] = AG;
      } else {
        int *AG = LOGICAL(xAG)-1;
        for(int i = l; i--; ) pout[i] = AG[pg[i]];
      }
      break;
    }
    default:
      stop("Not supported SEXP type!");
  }

  // Attribute Handling - 4 Situations:
  // 1 - x is classed (factor, date, time series), xAG is not classed. i.e. vector of fnobs, fmean etc.
  //    -> Sallow replacing, removing class and levels attributes from x, discard attributes of xAG (if any)
  //    -> or (if type matches i.e. double for date or time series), copy attributes of x unless x is a factor
  // 2 - x is not classed, xAG is classed (factor, date, time series). - an unusual situation should not occurr - copy attributes of xAG, discard attributes of x
  // 3 - xAG and x are classed - same as above, keep attributes of xAG, discard attributes of x
  // 4 - neither x nor xAG are classed - preserve attributes of x, discard attributes of xAG (if any)
  //

  if(Rf_isObject(xAG)) SHALLOW_DUPLICATE_ATTRIB(out, xAG);
  else if(!Rf_isObject(x) || (tx == txAG && !Rf_isFactor(x))) SHALLOW_DUPLICATE_ATTRIB(out, x);
  else {
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    Rf_classgets(out, R_NilValue); // OK !
    Rf_setAttrib(out, R_LevelsSymbol, R_NilValue); // if(Rf_isFactor(x)) ? faster ?
  }

  UNPROTECT(1);
  return out;
}

SEXP ret2(const SEXP& x, const SEXP& xAG, const SEXP& g) {
  int l = Rf_length(x), gs = Rf_length(g), tx = TYPEOF(x), txAG = TYPEOF(xAG);
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  int *pg;
  bool nog = gs == 1;
  if(nog) {
    if(Rf_length(xAG) != 1) stop("If g = NULL, NROW(STATS) needs to be 1");
  } else {
    if(gs != l) stop("length(g) must match NROW(x)");
    pg = INTEGER(g); // Wmaybe uninitialized
  }

  SEXP out = PROTECT(Rf_allocVector(txAG, l));

  switch(tx) {
  case REALSXP:
  {
    double *px = REAL(x);
    switch(txAG) {
      case REALSXP: {
        double *pout = REAL(out);
        if(nog) {
          double AG = Rf_asReal(xAG);
          for(int i = l; i--; ) pout[i] = (std::isnan(px[i])) ? NA_REAL : AG;
        } else {
          double *AG = REAL(xAG)-1;
          for(int i = l; i--; ) pout[i] = (std::isnan(px[i])) ? NA_REAL : AG[pg[i]];
        }
        break;
      }
      case INTSXP: {
        int *pout = INTEGER(out);
        if(nog) {
          int AG = Rf_asInteger(xAG);
          for(int i = l; i--; ) pout[i] = (std::isnan(px[i])) ? NA_INTEGER : AG;
        } else {
          int *AG = INTEGER(xAG)-1;
          for(int i = l; i--; ) pout[i] = (std::isnan(px[i])) ? NA_INTEGER : AG[pg[i]];
        }
        break;
      }
      case STRSXP: {
        SEXP *pout = STRING_PTR(out);
        if(nog) {
          SEXP AG = Rf_asChar(xAG);
          for(int i = l; i--; ) pout[i] = (std::isnan(px[i])) ? NA_STRING : AG;
        } else {
          SEXP *AG = STRING_PTR(xAG)-1;
          for(int i = l; i--; ) pout[i] = (std::isnan(px[i])) ? NA_STRING : AG[pg[i]];
        }
        break;
      }
      case LGLSXP: {
        int *pout = LOGICAL(out);
        if(nog) {
          int AG = Rf_asLogical(xAG);
          for(int i = l; i--; ) pout[i] = (std::isnan(px[i])) ? NA_LOGICAL : AG;
        } else {
          int *AG = LOGICAL(xAG)-1;
          for(int i = l; i--; ) pout[i] = (std::isnan(px[i])) ? NA_LOGICAL : AG[pg[i]];
        }
        break;
      }
      default:
        stop("Not supported SEXP type!");
    }
    break;
  }
  case INTSXP:
  {
    int *px = INTEGER(x);
    switch(txAG) {
      case REALSXP: {
        double *pout = REAL(out);
        if(nog) {
          double AG = Rf_asReal(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_INTEGER) ? NA_REAL : AG;
        } else {
          double *AG = REAL(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_INTEGER) ? NA_REAL : AG[pg[i]];
        }
        break;
      }
      case INTSXP: {
        int *pout = INTEGER(out);
        if(nog) {
          int AG = Rf_asInteger(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_INTEGER) ? NA_INTEGER : AG;
        } else {
          int *AG = INTEGER(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_INTEGER) ? NA_INTEGER : AG[pg[i]];
        }
        break;
      }
      case STRSXP: {
        SEXP *pout = STRING_PTR(out);
        if(nog) {
          SEXP AG = Rf_asChar(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_INTEGER) ? NA_STRING : AG;
        } else {
          SEXP *AG = STRING_PTR(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_INTEGER) ? NA_STRING : AG[pg[i]];
        }
        break;
      }
      case LGLSXP: {
        int *pout = LOGICAL(out);
        if(nog) {
          int AG = Rf_asLogical(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_INTEGER) ? NA_LOGICAL : AG;
        } else {
          int *AG = LOGICAL(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_INTEGER) ? NA_LOGICAL : AG[pg[i]];
        }
        break;
      }
      default:
        stop("Not supported SEXP type!");
    }
    break;
  }
  case STRSXP:
  {
    SEXP *px = STRING_PTR(x);
    switch(txAG) {
      case REALSXP: {
        double *pout = REAL(out);
        if(nog) {
          double AG = Rf_asReal(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_STRING) ? NA_REAL : AG;
        } else {
          double *AG = REAL(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_STRING) ? NA_REAL : AG[pg[i]];
        }
        break;
      }
      case INTSXP: {
        int *pout = INTEGER(out);
        if(nog) {
          int AG = Rf_asInteger(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_STRING) ? NA_INTEGER : AG;
        } else {
          int *AG = INTEGER(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_STRING) ? NA_INTEGER : AG[pg[i]];
        }
        break;
      }
      case STRSXP: {
        SEXP *pout = STRING_PTR(out);
        if(nog) {
          SEXP AG = Rf_asChar(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_STRING) ? NA_STRING : AG;
        } else {
          SEXP *AG = STRING_PTR(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_STRING) ? NA_STRING : AG[pg[i]];
        }
        break;
      }
      case LGLSXP: {
        int *pout = LOGICAL(out);
        if(nog) {
          int AG = Rf_asLogical(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_STRING) ? NA_LOGICAL : AG;
        } else {
          int *AG = LOGICAL(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_STRING) ? NA_LOGICAL : AG[pg[i]];
        }
        break;
      }
      default:
        stop("Not supported SEXP type!");
    }
    break;
  }
  case LGLSXP:
  {
    int *px = LOGICAL(x);
    switch(txAG) {
      case REALSXP: {
        double *pout = REAL(out);
        if(nog) {
          double AG = Rf_asReal(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_LOGICAL) ? NA_REAL : AG;
        } else {
          double *AG = REAL(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_LOGICAL) ? NA_REAL : AG[pg[i]];
        }
        break;
      }
      case INTSXP: {
        int *pout = INTEGER(out);
        if(nog) {
          int AG = Rf_asInteger(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_LOGICAL) ? NA_INTEGER : AG;
        } else {
          int *AG = INTEGER(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_LOGICAL) ? NA_INTEGER : AG[pg[i]];
        }
        break;
      }
      case STRSXP: {
        SEXP *pout = STRING_PTR(out);
        if(nog) {
          SEXP AG = Rf_asChar(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_LOGICAL) ? NA_STRING : AG;
        } else {
          SEXP *AG = STRING_PTR(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_LOGICAL) ? NA_STRING : AG[pg[i]];
        }
        break;
      }
      case LGLSXP: {
        int *pout = LOGICAL(out);
        if(nog) {
          int AG = Rf_asLogical(xAG);
          for(int i = l; i--; ) pout[i] = (px[i] == NA_LOGICAL) ? NA_LOGICAL : AG;
        } else {
          int *AG = LOGICAL(xAG)-1;
          for(int i = l; i--; ) pout[i] = (px[i] == NA_LOGICAL) ? NA_LOGICAL : AG[pg[i]];
        }
        break;
      }
      default:
        stop("Not supported SEXP type!");
    }
    break;
  }
  default:
    stop("Not supported SEXP type!");
  }

  if(Rf_isObject(xAG)) SHALLOW_DUPLICATE_ATTRIB(out, xAG);
  else if(!Rf_isObject(x) || (tx == txAG && !Rf_isFactor(x))) SHALLOW_DUPLICATE_ATTRIB(out, x);
  else {
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    Rf_classgets(out, R_NilValue); // OK !
    Rf_setAttrib(out, R_LevelsSymbol, R_NilValue);
  }

  UNPROTECT(1);
  return out;
}

// TODO: allow integer input ??
SEXP retoth(const NumericVector& x, const NumericVector& xAG, const SEXP& g, int ret = 3) {
  int gs = Rf_length(g), l = x.size();
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  NumericVector out = no_init_vector(l);
    if(gs == 1) {
      if(xAG.size() != 1) stop("If g = NULL, STATS needs to be an atomic element!");
      double AGx = xAG[0];
      switch(ret) {
      case 3:
        out = x - AGx;
        break;
      case 4: stop("This transformation can only be performed with groups!");
      case 5:
        out = x * (1/AGx);
        break;
      case 6:
        out = x * (100 / AGx);
        break;
      case 7:
        out = x + AGx;
        break;
      case 8:
        out = x * AGx;
        break;
      case 9:
        for(int i = 0; i != l; ++i) out[i] = mymod(x[i], AGx);
        break;
      case 10:
        for(int i = 0; i != l; ++i) out[i] = myremain(x[i], AGx);
        break;
      default: stop("Unknown Transformation");
      }
    } else {
      if(gs != l) stop("length(g) must match nrow(x)");
      double *px = REAL(x), *pout = REAL(out), *pAG = REAL(xAG)-1;
      int *pg = INTEGER(g);
      switch(ret) {
      case 3:
        for(int i = l; i--; ) pout[i] = px[i] - pAG[pg[i]];
        break;
      case 4:
        {
          long double OM = 0; // better precision
          int n = 0;
          for(int i = l; i--; ) {
            if(std::isnan(px[i])) pout[i] = px[i];
            else {
              pout[i] = px[i] - pAG[pg[i]];
              if(std::isnan(pAG[pg[i]])) continue; // If one AG remained NA, OM becomes NA
              OM += pAG[pg[i]];
              ++n;
            }
          }
          OM = OM / n;
          out = out + (double)OM; // Fastest ?
          break;
        }
      case 5:
        for(int i = l; i--; ) pout[i] = px[i] / pAG[pg[i]]; // Fastest ?
        break;
      case 6:
        for(int i = l; i--; ) pout[i] = px[i] * (100 / pAG[pg[i]]);
        break;
      case 7:
        for(int i = l; i--; ) pout[i] = px[i] + pAG[pg[i]];
        break;
      case 8:
        for(int i = l; i--; ) pout[i] = px[i] * pAG[pg[i]];
        break;
      case 9:
        for(int i = l; i--; ) pout[i] = mymod(px[i], pAG[pg[i]]);
        break;
      case 10:
        for(int i = l; i--; ) pout[i] = myremain(px[i], pAG[pg[i]]);
        break;
      default: stop("Unknown Transformation");
      }
    }
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    return out;
}

// [[Rcpp::export]]
SEXP TRACpp(const SEXP& x, const SEXP& xAG, const IntegerVector& g = 0, int ret = 1) {
  if(ret <= 2) {
    if(ret == 1) return ret1(x, xAG, g);
    return ret2(x, xAG, g);
  }
  return retoth(x, xAG, g, ret);
}


// [[Rcpp::export]]
List TRAlCpp(const List& x, const SEXP& xAG, const IntegerVector& g = 0, int ret = 1) {
  int l = x.size();
  if(Rf_length(xAG) != l) stop("NCOL(x) must match NCOL(STATS)");
  List out(l);

  switch(TYPEOF(xAG)) {
  case VECSXP: {
    List AG = xAG;
    if(ret == 1)      for(int j = l; j--; ) out[j] = ret1(x[j], AG[j], g);
    else if(ret == 2) for(int j = l; j--; ) out[j] = ret2(x[j], AG[j], g);
    else              for(int j = l; j--; ) out[j] = retoth(x[j], AG[j], g, ret);
    break;
  }
  case REALSXP: {
    NumericVector AG = xAG;
    if(ret == 1)      for(int j = l; j--; ) out[j] = ret1(x[j], Rf_ScalarReal(AG[j]), g);
    else if(ret == 2) for(int j = l; j--; ) out[j] = ret2(x[j], Rf_ScalarReal(AG[j]), g);
    else              for(int j = l; j--; ) out[j] = retoth(x[j], Rf_ScalarReal(AG[j]), g, ret);
    break;
  }
  case INTSXP: {
    IntegerVector AG = xAG;
    if(ret == 1)      for(int j = l; j--; ) out[j] = ret1(x[j], Rf_ScalarInteger(AG[j]), g);
    else if(ret == 2) for(int j = l; j--; ) out[j] = ret2(x[j], Rf_ScalarInteger(AG[j]), g);
    else              for(int j = l; j--; ) out[j] = retoth(x[j], Rf_ScalarInteger(AG[j]), g, ret);
    break;
  }
  case STRSXP: {
    CharacterVector AG = xAG;
    if(ret == 1)      for(int j = l; j--; ) out[j] = ret1(x[j], Rf_ScalarString(AG[j]), g); // Rf_ScalarString ? -> Not really necessary, a scalar string is still a SEXP ...
    else if(ret == 2) for(int j = l; j--; ) out[j] = ret2(x[j], Rf_ScalarString(AG[j]), g);
    else              stop("The requested transformation is not possible with strings");
    break;
  }
  case LGLSXP: {
    LogicalVector AG = xAG;
    if(ret == 1)      for(int j = l; j--; ) out[j] = ret1(x[j], Rf_ScalarLogical(AG[j]), g);
    else if(ret == 2) for(int j = l; j--; ) out[j] = ret2(x[j], Rf_ScalarLogical(AG[j]), g);
    else              stop("The requested transformation is not possible with logical data");
    break;
  }
  default: stop("Not supported SEXP type!");
  }
  SHALLOW_DUPLICATE_ATTRIB(out, x);
  return out;
}

// TODO: "replace" method for matrices is a bit slower than before, but overall pretty good!

// [[Rcpp::export]]
SEXP TRAmCpp(const SEXP& x, const SEXP& xAG, const IntegerVector& g = 0, int ret = 1) {
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  if(Rf_isNull(dim)) stop("x is not a matrix");
  int tx = TYPEOF(x), txAG = TYPEOF(xAG), gs = g.size(),
    row = INTEGER(dim)[0], col = INTEGER(dim)[1], *pg = INTEGER(g), ng = 0;
  bool nog = gs == 1;
  if(nog) {
    if(Rf_length(xAG) != col) stop("If g = NULL, NROW(STATS) needs to be 1");
  } else {
    if(gs != row) stop("length(g) must match ncol(x)");
    if(Rf_ncols(xAG) != col) stop("ncol(STATS) must match ncol(x)");
    ng = Rf_nrows(xAG);
  }

  if(ret <= 2) {

    SEXP out = PROTECT(Rf_allocVector(txAG, row * col));

    if(ret == 1) {

      switch(txAG) {
        case REALSXP:
        {
          double *pout = REAL(out);
          if(nog) {
            double *pAG = REAL(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = AG[pg[i]];
            }
          }
          break;
        }
        case INTSXP:
        {
          int *pout = INTEGER(out);
          if(nog) {
            int *pAG = INTEGER(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = INTEGER(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = AG[pg[i]];
            }
          }
          break;
        }
        case STRSXP:
        {
          SEXP *pout = STRING_PTR(out);
          if(nog) {
            SEXP *pAG = STRING_PTR(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              SEXP AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              SEXP *AG = STRING_PTR(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = AG[pg[i]];
            }
          }
          break;
        }
        case LGLSXP:
        {
          int *pout = LOGICAL(out);
          if(nog) {
            int *pAG = LOGICAL(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = LOGICAL(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = AG[pg[i]];
            }
          }
          break;
        }
        default:
          stop("Not supported SEXP type!");
      }

    } else { // ret == 2

      switch(tx) {
      case REALSXP:
      {
        double *px = REAL(x);
        switch(txAG) {
        case REALSXP:
        {
          double *pout = REAL(out);
          if(nog) {
            double *pAG = REAL(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (std::isnan(px[i])) ? NA_REAL : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (std::isnan(px[i + s])) ? NA_REAL : AG[pg[i]];
            }
          }
          break;
        }
        case INTSXP:
        {
          int *pout = INTEGER(out);
          if(nog) {
            int *pAG = INTEGER(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (std::isnan(px[i])) ? NA_INTEGER : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = INTEGER(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (std::isnan(px[i + s])) ? NA_INTEGER : AG[pg[i]];
            }
          }
          break;
        }
        case STRSXP:
        {
          SEXP *pout = STRING_PTR(out);
          if(nog) {
            SEXP *pAG = STRING_PTR(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              SEXP AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (std::isnan(px[i])) ? NA_STRING : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              SEXP *AG = STRING_PTR(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (std::isnan(px[i + s])) ? NA_STRING : AG[pg[i]];
            }
          }
          break;
        }
        case LGLSXP:
        {
          int *pout = LOGICAL(out);
          if(nog) {
            int *pAG = LOGICAL(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (std::isnan(px[i])) ? NA_LOGICAL : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = LOGICAL(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (std::isnan(px[i + s])) ? NA_LOGICAL : AG[pg[i]];
            }
          }
          break;
        }
        default:
          stop("Not supported SEXP type!");
        }
        break;
      }
      case INTSXP:
      {
        int *px = INTEGER(x);
        switch(txAG) {
        case REALSXP:
        {
          double *pout = REAL(out);
          if(nog) {
            double *pAG = REAL(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_REAL : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_INTEGER) ? NA_REAL : AG[pg[i]];
            }
          }
          break;
        }
        case INTSXP:
        {
          int *pout = INTEGER(out);
          if(nog) {
            int *pAG = INTEGER(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_INTEGER : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = INTEGER(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_INTEGER) ? NA_INTEGER : AG[pg[i]];
            }
          }
          break;
        }
        case STRSXP:
        {
          SEXP *pout = STRING_PTR(out);
          if(nog) {
            SEXP *pAG = STRING_PTR(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              SEXP AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_STRING : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              SEXP *AG = STRING_PTR(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_INTEGER) ? NA_STRING : AG[pg[i]];
            }
          }
          break;
        }
        case LGLSXP:
        {
          int *pout = LOGICAL(out);
          if(nog) {
            int *pAG = LOGICAL(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_LOGICAL : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = LOGICAL(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_INTEGER) ? NA_LOGICAL : AG[pg[i]];
            }
          }
          break;
        }
        default:
          stop("Not supported SEXP type!");
        }
        break;
      }
      case STRSXP:
      {
        SEXP *px = STRING_PTR(x);
        switch(txAG) {
        case REALSXP:
        {
          double *pout = REAL(out);
          if(nog) {
            double *pAG = REAL(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_STRING) ? NA_REAL : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_STRING) ? NA_REAL : AG[pg[i]];
            }
          }
          break;
        }
        case INTSXP:
        {
          int *pout = INTEGER(out);
          if(nog) {
            int *pAG = INTEGER(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_STRING) ? NA_INTEGER : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = INTEGER(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_STRING) ? NA_INTEGER : AG[pg[i]];
            }
          }
          break;
        }
        case STRSXP:
        {
          SEXP *pout = STRING_PTR(out);
          if(nog) {
            SEXP *pAG = STRING_PTR(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              SEXP AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_STRING) ? NA_STRING : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              SEXP *AG = STRING_PTR(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_STRING) ? NA_STRING : AG[pg[i]];
            }
          }
          break;
        }
        case LGLSXP:
        {
          int *pout = LOGICAL(out);
          if(nog) {
            int *pAG = LOGICAL(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_STRING) ? NA_LOGICAL : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = LOGICAL(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_STRING) ? NA_LOGICAL : AG[pg[i]];
            }
          }
          break;
        }
        default:
          stop("Not supported SEXP type!");
        }
        break;
      }
      case LGLSXP:
      {
        int *px = LOGICAL(x);
        switch(txAG) {
        case REALSXP:
        {
          double *pout = REAL(out);
          if(nog) {
            double *pAG = REAL(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_LOGICAL) ? NA_REAL : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_LOGICAL) ? NA_REAL : AG[pg[i]];
            }
          }
          break;
        }
        case INTSXP:
        {
          int *pout = INTEGER(out);
          if(nog) {
            int *pAG = INTEGER(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_LOGICAL) ? NA_INTEGER : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = INTEGER(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_LOGICAL) ? NA_INTEGER : AG[pg[i]];
            }
          }
          break;
        }
        case STRSXP:
        {
          SEXP *pout = STRING_PTR(out);
          if(nog) {
            SEXP *pAG = STRING_PTR(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              SEXP AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_LOGICAL) ? NA_STRING : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              SEXP *AG = STRING_PTR(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_LOGICAL) ? NA_STRING : AG[pg[i]];
            }
          }
          break;
        }
        case LGLSXP:
        {
          int *pout = LOGICAL(out);
          if(nog) {
            int *pAG = LOGICAL(xAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_LOGICAL) ? NA_LOGICAL : AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = LOGICAL(xAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_LOGICAL) ? NA_LOGICAL : AG[pg[i]];
            }
          }
          break;
        }
        default:
          stop("Not supported SEXP type!");
        }
        break;
      }
      default:
        stop("Not supported SEXP type!");
      }

    }

    if(Rf_isObject(xAG)) SHALLOW_DUPLICATE_ATTRIB(out, xAG);
    else if(!Rf_isObject(x) || (tx == txAG && !Rf_isFactor(x))) SHALLOW_DUPLICATE_ATTRIB(out, x);
    else {
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      Rf_classgets(out, R_NilValue); // OK !
      Rf_setAttrib(out, R_LevelsSymbol, R_NilValue);
    }

    UNPROTECT(1);
    return out;

  } else { // ret > 2

    SEXP out = PROTECT(Rf_allocVector(REALSXP, row * col));
    double *pout = REAL(out), *px;
    NumericVector xxAG = xAG;

    switch(tx) {
    case REALSXP:
    case INTSXP:
    {
      if(tx == INTSXP) { // TODO: Better solution !
        NumericVector xx = x; // Rf_coerceVector(x, REALSXP);
        px = REAL(xx);
      } else {
        px = REAL(x);
      }

      switch(ret) {
        case 3: {
          if(nog) {
            double *pAG = REAL(xxAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = px[i] - AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xxAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = px[i + s] - AG[pg[i]];
            }
          }
          break;
        }
        case 4: {
          if(nog) stop("This transformation can only be computed with groups!");
          for(int j = 0; j != col; ++j) {
            int s = j * row, n = 0;
            long double OM = 0;
            double *AG = REAL(xxAG) + j * ng - 1;
            for(int i = 0; i != row; ++i) {
              if(std::isnan(px[i + s])) pout[i + s] = px[i + s];
              else {
                pout[i + s] = px[i + s] - AG[pg[i]];
                if(std::isnan(AG[pg[i]])) continue;
                OM += AG[pg[i]];
                ++n;
              }
            }
            double OMD = double(OM / n);
            for(int i = row; i--; ) pout[i + s] += OMD;
          }
          break;
        }
        case 5: {
          if(nog) {
            double *pAG = REAL(xxAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = 1 / pAG[j];
              for(int i = s; i != e; ++i) pout[i] = px[i] * AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xxAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = px[i + s] * (1 / AG[pg[i]]);
            }
          }
          break;
        }
        case 6: {
          if(nog) {
            double *pAG = REAL(xxAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = 100 / pAG[j];
              for(int i = s; i != e; ++i) pout[i] = px[i] * AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xxAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = px[i + s] * (100 / AG[pg[i]]);
            }
          }
          break;
        }
        case 7: {
          if(nog) {
            double *pAG = REAL(xxAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = px[i] + AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xxAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = px[i + s] + AG[pg[i]];
            }
          }
          break;
        }
        case 8: {
          if(nog) {
            double *pAG = REAL(xxAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = px[i] * AGj;
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xxAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = px[i + s] * AG[pg[i]];
            }
          }
          break;
        }
        case 9: {
          if(nog) {
            double *pAG = REAL(xxAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = mymod(px[i], AGj);
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xxAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = mymod(px[i + s], AG[pg[i]]);
            }
          }
          break;
        }
        case 10: {
          if(nog) {
            double *pAG = REAL(xxAG);
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = myremain(px[i], AGj);
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = REAL(xxAG) + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = myremain(px[i + s], AG[pg[i]]);
            }
          }
          break;
        }
        default: stop("Unknown Transformation");
      }
      break;
    }
    case STRSXP: stop("The requested transformation is not possible with strings");
    case LGLSXP: stop("The requested transformation is not possible with logical data");
    default: stop("Not supported SEXP type!");
    }

    SHALLOW_DUPLICATE_ATTRIB(out, x);

    UNPROTECT(1);
    return out;

  }
}





