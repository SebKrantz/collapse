#include "collapse_c.h"


// Cases:
// 0- replace_NA (only replace missing values)
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
static inline double modulus_impl(double x, double y) {
  double z = x * (1/y);
  return (z == z) ? x - (int)(z) * y : z; // faster than x - (int)(x/y) * y;
}

// #define modulus_impl(x, y) (x - ((int)(x/y) * y)) // Macro: not faster !

// template<typename T, typename U>
// constexpr double modulus_impl (T x, U mod)
// {
//   return !mod ? x : x - mod * static_cast<long long>(x / mod);
// }

// int(x * (1/y)) -> This gave the UBSAN error if NaN !!!
static inline double remainder_impl(double x, double y) {
  double z = x * (1/y);
  return (z == z) ? (int)(z) * y : z; //   (int)(x * (1/y)) * y; <- This would be enough, but doesn't keep missing values in x!
}


int TtI(SEXP x) {
  if(TYPEOF(x) != STRSXP) error("FUN must be integer or character");
  const char * r = CHAR(STRING_ELT(x, 0)); // translateCharUTF8()
  if(strcmp(r, "replace_NA") == 0) return 0;
  if(strcmp(r, "replace_fill") == 0) return 1;
  if(strcmp(r, "replace") == 0) return 2;
  if(strcmp(r, "-") == 0) return 3;
  if(strcmp(r, "-+") == 0) return 4;
  if(strcmp(r, "/") == 0) return 5;
  if(strcmp(r, "%") == 0) return 6;
  if(strcmp(r, "+") == 0) return 7;
  if(strcmp(r, "*") == 0) return 8;
  if(strcmp(r, "%%") == 0) return 9;
  if(strcmp(r, "-%%") == 0) return 10;
  error("Unknown transformation!");
}


SEXP ret1(SEXP x, SEXP xAG, SEXP g, int set) {
  int tx = TYPEOF(x), txAG = TYPEOF(xAG), l = length(x), gs = length(g);
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  int *pg = &l, nog = gs <= 1;
  if(nog) {
    if(length(xAG) != 1) error("If g = NULL, NROW(STATS) needs to be 1");
  } else {
    if(TYPEOF(g) != INTSXP) error("g must be integer typed, please report this as g should have been internally grouped");
    if(gs != l) error("length(g) must match NROW(x)");
    pg = INTEGER(g);
  }

  if(set && txAG != tx) error("if set = TRUE with option 'replace_fill', x and STATS need to have identical data types");

  SEXP out = set == 0 ? PROTECT(allocVector(txAG, l)) : x;

  switch(txAG) {
    case REALSXP:
    {
      double *pout = REAL(out);
      if(nog) {
        double AG = asReal(xAG);
        for(int i = 0; i != l; ++i) pout[i] = AG;
      } else {
        double *AG = REAL(xAG)-1;
        for(int i = 0; i != l; ++i) pout[i] = AG[pg[i]];
      }
      break;
    }
    case INTSXP:
    case LGLSXP:
    {
      int *pout = INTEGER(out);
      if(nog) {
        int AG = asInteger(xAG);
        for(int i = 0; i != l; ++i) pout[i] = AG;
      } else {
        int *AG = INTEGER(xAG)-1;
        for(int i = 0; i != l; ++i) pout[i] = AG[pg[i]];
      }
      break;
    }
    case CPLXSXP:
    {
      Rcomplex *pout = COMPLEX(out);
      if(nog) {
        Rcomplex AG = asComplex(xAG);
        for(int i = 0; i != l; ++i) pout[i] = AG;
      } else {
        Rcomplex *AG = COMPLEX(xAG)-1;
        for(int i = 0; i != l; ++i) pout[i] = AG[pg[i]];
      }
      break;
    }
    case STRSXP:
    {
      SEXP *pout = STRING_PTR(out);
      if(nog) {
        SEXP AG = asChar(xAG);
        for(int i = 0; i != l; ++i) pout[i] = AG;
      } else {
        SEXP *AG = STRING_PTR(xAG)-1;
        for(int i = 0; i != l; ++i) pout[i] = AG[pg[i]];
      }
      break;
    }
    case VECSXP:
    {
      SEXP *pout = SEXPPTR(out);
      if(nog) {
        for(int i = 0; i != l; ++i) pout[i] = xAG;
      } else {
        SEXP *AG = SEXPPTR(xAG)-1;
        for(int i = 0; i != l; ++i) pout[i] = AG[pg[i]];
      }
      break;
    }
    case RAWSXP:
    {
      Rbyte *pout = RAW(out);
      if(nog) {
        Rbyte AG = RAW_ELT(xAG, 0);
        for(int i = 0; i != l; ++i) pout[i] = AG;
      } else {
        Rbyte *AG = RAW(xAG)-1;
        for(int i = 0; i != l; ++i) pout[i] = AG[pg[i]];
      }
      break;
    }
    default: error("Not supported SEXP type!");
  }

  // Attribute Handling - 4 Situations:
  // 1 - x is classed (factor, date, time series), xAG is not classed. i.e. vector of fnobs, fmean etc.
  //    -> Sallow replacing, removing class and levels attributes from x, discard attributes of xAG (if any)
  //    -> or (if type matches i.e. double for date or time series), copy attributes of x unless x is a factor
  // 2 - x is not classed, xAG is classed (factor, date, time series). - an unusual situation should not occurr - copy attributes of xAG, discard attributes of x
  // 3 - xAG and x are classed - same as above, keep attributes of xAG, discard attributes of x
  // 4 - neither x nor xAG are classed - preserve attributes of x, discard attributes of xAG (if any)
  //
  if(set == 0) {
    if(isObject(xAG)) SHALLOW_DUPLICATE_ATTRIB(out, xAG);
    else if(!isObject(x) || (tx == txAG && !isFactor(x))) SHALLOW_DUPLICATE_ATTRIB(out, x);
    else {
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      classgets(out, R_NilValue); // OK !
      setAttrib(out, R_LevelsSymbol, R_NilValue); // if(isFactor(x)) ? faster ?
    }
    UNPROTECT(1);
  }
  return out;
}

SEXP ret2(SEXP x, SEXP xAG, SEXP g, int set) {
  int l = length(x), gs = length(g), tx = TYPEOF(x), txAG = TYPEOF(xAG);

  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  int *pg = &l, nog = gs <= 1;
  if(nog) {
    if(length(xAG) != 1) error("If g = NULL, NROW(STATS) needs to be 1");
  } else {
    if(TYPEOF(g) != INTSXP) error("g must be integer typed, please report this as g should have been internally grouped");
    if(gs != l) error("length(g) must match NROW(x)");
    pg = INTEGER(g); // Wmaybe uninitialized
  }

  if(set && txAG != tx) error("if set = TRUE with option 'replace', x and STATS need to have identical data types");

  SEXP out = set == 0 ? PROTECT(allocVector(txAG, l)) : x;

  switch(tx) {
  case REALSXP:
  {
    double *px = REAL(x);
    switch(txAG) {
      case REALSXP: {
        double *pout = REAL(out);
        if(nog) {
          double AG = asReal(xAG);
          for(int i = 0; i != l; ++i) pout[i] = ISNAN(px[i]) ? NA_REAL : AG;
        } else {
          double *AG = REAL(xAG)-1;
          for(int i = 0; i != l; ++i) pout[i] = ISNAN(px[i]) ? NA_REAL : AG[pg[i]];
        }
        break;
      }
      case LGLSXP:
      case INTSXP: {
        int *pout = INTEGER(out);
        if(nog) {
          int AG = asInteger(xAG);
          for(int i = 0; i != l; ++i) pout[i] = ISNAN(px[i]) ? NA_INTEGER : AG;
        } else {
          int *AG = INTEGER(xAG)-1;
          for(int i = 0; i != l; ++i) pout[i] = ISNAN(px[i]) ? NA_INTEGER : AG[pg[i]];
        }
        break;
      }
      case STRSXP: {
        SEXP *pout = STRING_PTR(out);
        if(nog) {
          SEXP AG = asChar(xAG);
          for(int i = 0; i != l; ++i) pout[i] = ISNAN(px[i]) ? NA_STRING : AG;
        } else {
          SEXP *AG = STRING_PTR(xAG)-1;
          for(int i = 0; i != l; ++i) pout[i] = ISNAN(px[i]) ? NA_STRING : AG[pg[i]];
        }
        break;
      }
      default:
        error("Not supported SEXP type!");
    }
    break;
  }
  case INTSXP:
  case LGLSXP:
  {
    int *px = INTEGER(x);
    switch(txAG) {
      case REALSXP: {
        double *pout = REAL(out);
        if(nog) {
          double AG = asReal(xAG);
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_REAL : AG;
        } else {
          double *AG = REAL(xAG)-1;
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_REAL : AG[pg[i]];
        }
        break;
      }
      case LGLSXP:
      case INTSXP: {
        int *pout = INTEGER(out);
        if(nog) {
          int AG = asInteger(xAG);
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_INTEGER : AG;
        } else {
          int *AG = INTEGER(xAG)-1;
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_INTEGER : AG[pg[i]];
        }
        break;
      }
      case STRSXP: {
        SEXP *pout = STRING_PTR(out);
        if(nog) {
          SEXP AG = asChar(xAG);
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_STRING : AG;
        } else {
          SEXP *AG = STRING_PTR(xAG)-1;
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_STRING : AG[pg[i]];
        }
        break;
      }
      default:
        error("Not supported SEXP type!");
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
          double AG = asReal(xAG);
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_STRING) ? NA_REAL : AG;
        } else {
          double *AG = REAL(xAG)-1;
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_STRING) ? NA_REAL : AG[pg[i]];
        }
        break;
      }
      case LGLSXP:
      case INTSXP: {
        int *pout = INTEGER(out);
        if(nog) {
          int AG = asInteger(xAG);
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_STRING) ? NA_INTEGER : AG;
        } else {
          int *AG = INTEGER(xAG)-1;
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_STRING) ? NA_INTEGER : AG[pg[i]];
        }
        break;
      }
      case STRSXP: {
        SEXP *pout = STRING_PTR(out);
        if(nog) {
          SEXP AG = asChar(xAG);
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_STRING) ? NA_STRING : AG;
        } else {
          SEXP *AG = STRING_PTR(xAG)-1;
          for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_STRING) ? NA_STRING : AG[pg[i]];
        }
        break;
      }
      default:
        error("Not supported SEXP type!");
    }
    break;
  }
  default:
    error("Not supported SEXP type!");
  }

  if(set == 0) {
    if(isObject(xAG)) SHALLOW_DUPLICATE_ATTRIB(out, xAG);
    else if(!isObject(x) || (tx == txAG && !isFactor(x))) SHALLOW_DUPLICATE_ATTRIB(out, x);
    else {
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      classgets(out, R_NilValue); // OK !
      setAttrib(out, R_LevelsSymbol, R_NilValue);
    }
    UNPROTECT(1);
  }
  return out;
}

// New: Option "replace_NA"
SEXP ret0(SEXP x, SEXP xAG, SEXP g, int set) {
  int l = length(x), gs = length(g), tx = TYPEOF(x), txAG = TYPEOF(xAG);
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  int *pg = &l, nog = gs <= 1;
  if(nog) {
    if(length(xAG) != 1) error("If g = NULL, NROW(STATS) needs to be 1");
  } else {
    if(TYPEOF(g) != INTSXP) error("g must be integer typed, please report this as g should have been internally grouped");
    if(gs != l) error("length(g) must match NROW(x)");
    pg = INTEGER(g); // Wmaybe uninitialized
  }

  SEXP out = set == 0 ? PROTECT(allocVector(tx, l)) : x;

  switch(tx) {
    case REALSXP:
    {
      double *px = REAL(x), *pout = REAL(out);
      if(nog) {
        if(txAG != REALSXP && txAG != INTSXP && txAG != LGLSXP) error("STATS needs to be numeric to replace NA's in numeric data!");
        double AG = asReal(xAG);
        for(int i = 0; i != l; ++i) pout[i] = ISNAN(px[i]) ? AG : px[i];
      } else {
        switch(txAG) {
          case REALSXP: {
            double *AG = REAL(xAG)-1;
            for(int i = 0; i != l; ++i) pout[i] = ISNAN(px[i]) ? AG[pg[i]] : px[i];
            break;
          }
          case LGLSXP:
          case INTSXP: {
            int *AG = INTEGER(xAG)-1;
            for(int i = 0; i != l; ++i) pout[i] = ISNAN(px[i]) ? AG[pg[i]] : px[i];
            break;
          }
          case STRSXP: error("Cannot replace missing values in double with a string");
          default: error("Not supported SEXP type!");
        }
      }
      break;
    }
    case LGLSXP:
    case INTSXP:
    {
      int *px = INTEGER(x), *pout = INTEGER(out);
      if(nog) {
        if(txAG != REALSXP && txAG != INTSXP && txAG != LGLSXP) error("STATS needs to be numeric to replace NA's in numeric data!");
        int AG = asInteger(xAG);
        for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_INTEGER) ? AG : px[i];
      } else {
        switch(txAG) {
          case REALSXP: {
            double *AG = REAL(xAG)-1;
            for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_INTEGER) ? AG[pg[i]] : px[i];
            break;
          }
          case LGLSXP:
          case INTSXP: {
            int *AG = INTEGER(xAG)-1;
            for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_INTEGER) ? AG[pg[i]] : px[i];
            break;
          }
          case STRSXP: error("Cannot replace missing values in integer with a string");
          default: error("Not supported SEXP type!");
        }
      }
      break;
    }
    case STRSXP:
    {
      SEXP *px = STRING_PTR(x), *pout = STRING_PTR(out);
      if(nog) {
        SEXP AG = asChar(xAG);
        for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_STRING) ? AG : px[i];
      } else {
        switch(txAG) {
          case REALSXP:
          case LGLSXP:
          case INTSXP: error("Cannot replace missing values in string with numeric data");
          case STRSXP: {
            SEXP *AG = STRING_PTR(xAG)-1;
            for(int i = 0; i != l; ++i) pout[i] = (px[i] == NA_STRING) ? AG[pg[i]] : px[i];
            break;
          }
          default: error("Not supported SEXP type!");
        }
      }
      break;
    }
    default:
      error("Not supported SEXP type!");
  }

  if(set == 0) {
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    UNPROTECT(1);
  }
  return out;
}

// TODO: allow integer input ??
SEXP retoth(SEXP x, SEXP xAG, SEXP g, int ret, int set) {
  int gs = length(g), l = length(x), txAG = TYPEOF(xAG);
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  SEXP out = set == 0 ? PROTECT(allocVector(REALSXP, l)) : x;

  if(gs <= 1) {
      if(length(xAG) != 1) error("If g = NULL, STATS needs to be an atomic element!");
      if(txAG != REALSXP && txAG != INTSXP && txAG != LGLSXP) error("for these transformations STATS needs to be numeric!");

  #define NOGOPLOOP                                                                    \
      switch(ret) {                                                                    \
      case 3:                                                                          \
        for(int i = 0; i != l; ++i) pout[i] = px[i] - AGx;                             \
        break;                                                                         \
      case 4: error("This transformation can only be performed with groups!");         \
      case 5: {                                                                        \
        double v = 1 / AGx;                                                            \
        for(int i = 0; i != l; ++i) pout[i] = px[i] * v;                               \
        break;                                                                         \
      }                                                                                \
      case 6: {                                                                        \
        double v = 100 / AGx;                                                          \
        for(int i = 0; i != l; ++i) pout[i] = px[i] * v;                               \
        break;                                                                         \
      }                                                                                \
      case 7:                                                                          \
        for(int i = 0; i != l; ++i) pout[i] = px[i] + AGx;                             \
        break;                                                                         \
      case 8:                                                                          \
        for(int i = 0; i != l; ++i) pout[i] = px[i] * AGx;                             \
        break;                                                                         \
      case 9:                                                                          \
        for(int i = 0; i != l; ++i) pout[i] = modulus_impl(px[i], AGx);                       \
        break;                                                                         \
      case 10:                                                                         \
        for(int i = 0; i != l; ++i) pout[i] = remainder_impl(px[i], AGx);                    \
        break;                                                                         \
      default: error("Unknown Transformation");                                        \
      }

      switch(TYPEOF(x)) {
        case REALSXP: {
          double AGx = asReal(xAG), *pout = REAL(out), *px = REAL(x);
          NOGOPLOOP
          break;
        }
        case INTSXP:
        case LGLSXP: {
          if(set) {
            int AGx = asInteger(xAG);
            int *pout = INTEGER(out), *px = INTEGER(x);
            NOGOPLOOP
          } else {
            double AGx = asReal(xAG), *pout = REAL(out);
            int *px = INTEGER(x);
            NOGOPLOOP
          }
          break;
        }
        default: error("x needs to be double or integer");
      }


  } else {
      if(TYPEOF(g) != INTSXP) error("g must be integer typed, please report this as g should have been internally grouped");
      if(gs != l) error("length(g) must match nrow(x)");
      int *pg = INTEGER(g);

    #define GOPLOOP                                                           \
      switch(ret) {                                                           \
        case 3:                                                               \
          for(int i = 0; i != l; ++i) pout[i] = px[i] - pAG[pg[i]];           \
          break;                                                              \
        case 4:                                                               \
          {                                                                   \
            long double OM = 0;                                               \
            int n = 0;                                                        \
            for(int i = 0; i != l; ++i) {                                     \
              if(ISNAN(px[i])) pout[i] = px[i];                               \
              else {                                                          \
                pout[i] = px[i] - pAG[pg[i]];                                 \
                if(ISNAN(pAG[pg[i]])) continue;                               \
                OM += pAG[pg[i]];                                             \
                ++n;                                                          \
              }                                                               \
            }                                                                 \
            OM /= n;                                                          \
            double dOM = (double)OM;                                          \
            for(int i = 0; i != l; ++i) pout[i] += dOM;                       \
            break;                                                            \
          }                                                                   \
        case 5:                                                               \
          for(int i = 0; i != l; ++i) pout[i] = px[i] / pAG[pg[i]];           \
          break;                                                              \
        case 6:                                                               \
          for(int i = 0; i != l; ++i) pout[i] = px[i] / pAG[pg[i]] * 100;     \
          break;                                                              \
        case 7:                                                               \
          for(int i = 0; i != l; ++i) pout[i] = px[i] + pAG[pg[i]];           \
          break;                                                              \
        case 8:                                                               \
          for(int i = 0; i != l; ++i) pout[i] = px[i] * pAG[pg[i]];           \
          break;                                                              \
        case 9:                                                               \
          for(int i = 0; i != l; ++i) pout[i] = modulus_impl(px[i], pAG[pg[i]]);     \
          break;                                                              \
        case 10:                                                              \
          for(int i = 0; i != l; ++i) pout[i] = remainder_impl(px[i], pAG[pg[i]]);  \
          break;                                                              \
        default: error("Unknown Transformation");                             \
      }

    #define TXAGSWITCH                                                                       \
      switch(txAG) {                                                                         \
        case REALSXP: {                                                                      \
          double *pAG = REAL(xAG)-1;                                                         \
          GOPLOOP                                                                            \
          break;                                                                             \
        }                                                                                    \
        case INTSXP:                                                                         \
        case LGLSXP: {                                                                       \
          int *pAG = INTEGER(xAG)-1;                                                         \
          GOPLOOP                                                                            \
          break;                                                                             \
        }                                                                                    \
        default: error("STATS needs to be integer or real for statistical transformations"); \
      }

      switch(TYPEOF(x)) {
        case REALSXP: {
          double *px = REAL(x), *pout = REAL(out);
          TXAGSWITCH
          break;
        }
        case INTSXP:
        case LGLSXP: {
        int *px = INTEGER(x);
        if(set) {
          int *pout = INTEGER(out);
          TXAGSWITCH
        } else {
          double *pout = REAL(out);
          TXAGSWITCH
        }
        break;
        }
        default: error("x needs to be double or integer");
      }

  }
  if(set == 0) {
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    UNPROTECT(1);
  }
  return out;
}

SEXP TRAC(SEXP x, SEXP xAG, SEXP g, SEXP Rret, SEXP Rset) {
  if(length(Rret) != 1) error("can only perform one transformation at a time");
  int ret = (TYPEOF(Rret) == STRSXP) ? TtI(Rret) : asInteger(Rret), set = asLogical(Rset);
  switch(ret) {
    case 0: return ret0(x, xAG, g, set);
    case 1: return ret1(x, xAG, g, set);
    case 2: return ret2(x, xAG, g, set);
    default: return retoth(x, xAG, g, ret, set);
  }
}

SEXP TRAlC(SEXP x, SEXP xAG, SEXP g, SEXP Rret, SEXP Rset) {
  if(length(Rret) != 1) error("can only perform one transformation at a time");
  int l = length(x), set = asLogical(Rset),
    ret = (TYPEOF(Rret) == STRSXP) ? TtI(Rret) : asInteger(Rret);

  if(length(xAG) != l) error("NCOL(x) must match NCOL(STATS)");

  // This is allocated anyway, but not returned if set = TRUE
  SEXP out = PROTECT(allocVector(VECSXP, l));
  SEXP *px = SEXPPTR(x);

  // Need SET_VECTOR_ELT here because we are allocating... (otherwise sometimes segfault)
#define RETLOOPS(v)                                                                       \
  switch(ret) {                                                                           \
  case 0:                                                                                 \
    for(int j = 0; j != l; ++j) {                                                         \
      SET_VECTOR_ELT(out, j, ret0(px[j], PROTECT(v), g, set)); UNPROTECT(1);              \
    }                                                                                     \
    break;                                                                                \
  case 1:                                                                                 \
    for(int j = 0; j != l; ++j) {                                                         \
      SET_VECTOR_ELT(out, j, ret1(px[j], PROTECT(v), g, set)); UNPROTECT(1);              \
    }                                                                                     \
    break;                                                                                \
  case 2:                                                                                 \
    for(int j = 0; j != l; ++j) {                                                         \
      SET_VECTOR_ELT(out, j, ret2(px[j], PROTECT(v), g, set)); UNPROTECT(1);              \
    }                                                                                     \
    break;                                                                                \
  default:                                                                                \
    for(int j = 0; j != l; ++j) {                                                         \
      SET_VECTOR_ELT(out, j, retoth(px[j], PROTECT(v), g, ret, set)); UNPROTECT(1);       \
    }                                                                                     \
  }

  switch(TYPEOF(xAG)) {
    case VECSXP: {
      SEXP *pAG = SEXPPTR(xAG);
      RETLOOPS(pAG[j])
      break;
    }
    case REALSXP: {
      double *pAG = REAL(xAG);
      RETLOOPS(ScalarReal(pAG[j]))
      break;
    }
    case LGLSXP:
    case INTSXP: {
      int *pAG = INTEGER(xAG);
      RETLOOPS(ScalarInteger(pAG[j]))
      break;
    }
    case CPLXSXP: {
      Rcomplex *pAG = COMPLEX(xAG);
      RETLOOPS(ScalarComplex(pAG[j]))
      break;
    }
    case RAWSXP: {
      Rbyte *pAG = RAW(xAG);
      RETLOOPS(ScalarRaw(pAG[j]))
      break;
    }
    case STRSXP: {
      SEXP *pAG = STRING_PTR(xAG);
      RETLOOPS(ScalarString(pAG[j]))
      break;
    }
    default: error("Not supported SEXP type!");
  }

  if(set == 0) SHALLOW_DUPLICATE_ATTRIB(out, x);
  UNPROTECT(1);
  return set ? x : out;
}

// TODO: "replace" method for matrices is a bit slower than before, but overall pretty good!

SEXP TRAmC(SEXP x, SEXP xAG, SEXP g, SEXP Rret, SEXP Rset) {
  SEXP dim = getAttrib(x, R_DimSymbol);
  if(isNull(dim)) error("x is not a matrix");
  if(length(Rret) != 1) error("can only perform one transformation at a time");
  int tx = TYPEOF(x), txAG = TYPEOF(xAG), gs = length(g),
    row = INTEGER(dim)[0], col = INTEGER(dim)[1], *pg = &gs, ng = 0,
    set = asLogical(Rset),
    ret = (TYPEOF(Rret) == STRSXP) ? TtI(Rret) : asInteger(Rret),
    nog = gs <= 1;

  if(nog) {
    if(length(xAG) != col) error("If g = NULL, NROW(STATS) needs to be 1");
  } else {
    if(TYPEOF(g) != INTSXP) error("g must be integer typed, please report this as g should have been internally grouped");
    if(gs != row) error("length(g) must match ncol(x)");
    if(ncols(xAG) != col) error("ncol(STATS) must match ncol(x)");
    pg = INTEGER(g);
    ng = nrows(xAG);
  }

  if(ret <= 2) {

    if(ret > 0) {

      if(set && txAG != tx) error("if set = TRUE with option 'replace_fill', x and STATS need to have identical data types");

      SEXP out = set ? x : PROTECT(allocVector(txAG, row * col));

      if(ret == 1) {

        switch(txAG) {
          case REALSXP:
          {
            double *pout = REAL(out), *pAG = REAL(xAG);
            if(nog) {
              for(int j = 0; j != col; ++j) {
                int s = j * row, e = s + row;
                double AGj = pAG[j];
                for(int i = s; i != e; ++i) pout[i] = AGj;
              }
            } else {
              for(int j = 0; j != col; ++j) {
                int s = j * row;
                double *AG = pAG + j * ng - 1;
                for(int i = 0; i != row; ++i) pout[i + s] = AG[pg[i]];
              }
            }
            break;
          }
          case INTSXP:
          case LGLSXP:
          {
            int *pout = INTEGER(out), *pAG = INTEGER(xAG);
            if(nog) {
              for(int j = 0; j != col; ++j) {
                int s = j * row, e = s + row, AGj = pAG[j];
                for(int i = s; i != e; ++i) pout[i] = AGj;
              }
            } else {
              for(int j = 0; j != col; ++j) {
                int s = j * row, *AG = pAG + j * ng - 1;
                for(int i = 0; i != row; ++i) pout[i + s] = AG[pg[i]];
              }
            }
            break;
          }
          case STRSXP:
          {
            SEXP *pout = STRING_PTR(out), *pAG = STRING_PTR(xAG);
            if(nog) {
              for(int j = 0; j != col; ++j) {
                int s = j * row, e = s + row;
                SEXP AGj = pAG[j];
                for(int i = s; i != e; ++i) pout[i] = AGj;
              }
            } else {
              for(int j = 0; j != col; ++j) {
                int s = j * row;
                SEXP *AG = pAG + j * ng - 1;
                for(int i = 0; i != row; ++i) pout[i + s] = AG[pg[i]];
              }
            }
            break;
          }
          default: error("Not supported SEXP type!");
        }

      } else {

        switch(tx) {
          case REALSXP:
          {
            double *px = REAL(x);
            switch(txAG) {
              case REALSXP:
              {
                double *pout = REAL(out), *pAG = REAL(xAG);
                if(nog) {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, e = s + row;
                    double AGj = pAG[j];
                    for(int i = s; i != e; ++i) pout[i] = (ISNAN(px[i])) ? NA_REAL : AGj;
                  }
                } else {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row;
                    double *AG = pAG + j * ng - 1;
                    for(int i = 0; i != row; ++i) pout[i + s] = (ISNAN(px[i + s])) ? NA_REAL : AG[pg[i]];
                  }
                }
                break;
              }
              case INTSXP:
              case LGLSXP:
              {
                int *pout = INTEGER(out), *pAG = INTEGER(xAG);
                if(nog) {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, e = s + row, AGj = pAG[j];
                    for(int i = s; i != e; ++i) pout[i] = (ISNAN(px[i])) ? NA_INTEGER : AGj;
                  }
                } else {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, *AG = pAG + j * ng - 1;
                    for(int i = 0; i != row; ++i) pout[i + s] = (ISNAN(px[i + s])) ? NA_INTEGER : AG[pg[i]];
                  }
                }
                break;
              }
              case STRSXP:
              {
                SEXP *pout = STRING_PTR(out), *pAG = STRING_PTR(xAG);
                if(nog) {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, e = s + row;
                    SEXP AGj = pAG[j];
                    for(int i = s; i != e; ++i) pout[i] = (ISNAN(px[i])) ? NA_STRING : AGj;
                  }
                } else {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row;
                    SEXP *AG = pAG + j * ng - 1;
                    for(int i = 0; i != row; ++i) pout[i + s] = (ISNAN(px[i + s])) ? NA_STRING : AG[pg[i]];
                  }
                }
                break;
              }
              default: error("Not supported SEXP type!");
            }
            break;
          }
          case INTSXP:
          case LGLSXP:
          {
            int *px = INTEGER(x);
            switch(txAG) {
              case REALSXP:
              {
                double *pout = REAL(out), *pAG = REAL(xAG);
                if(nog) {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, e = s + row;
                    double AGj = pAG[j];
                    for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_REAL : AGj;
                  }
                } else {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row;
                    double *AG = pAG + j * ng - 1;
                    for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_INTEGER) ? NA_REAL : AG[pg[i]];
                  }
                }
                break;
              }
              case INTSXP:
              case LGLSXP:
              {
                int *pout = INTEGER(out), *pAG = INTEGER(xAG);
                if(nog) {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, e = s + row, AGj = pAG[j];
                    for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_INTEGER : AGj;
                  }
                } else {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, *AG = pAG + j * ng - 1;
                    for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_INTEGER) ? NA_INTEGER : AG[pg[i]];
                  }
                }
                break;
              }
              case STRSXP:
              {
                SEXP *pout = STRING_PTR(out), *pAG = STRING_PTR(xAG);
                if(nog) {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, e = s + row;
                    SEXP AGj = pAG[j];
                    for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_INTEGER) ? NA_STRING : AGj;
                  }
                } else {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row;
                    SEXP *AG = pAG + j * ng - 1;
                    for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_INTEGER) ? NA_STRING : AG[pg[i]];
                  }
                }
                break;
              }
              default: error("Not supported SEXP type!");
            }
            break;
          }
          case STRSXP:
          {
            SEXP *px = STRING_PTR(x);
            switch(txAG) {
              case REALSXP:
              {
                double *pout = REAL(out), *pAG = REAL(xAG);
                if(nog) {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, e = s + row;
                    double AGj = pAG[j];
                    for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_STRING) ? NA_REAL : AGj;
                  }
                } else {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row;
                    double *AG = pAG + j * ng - 1;
                    for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_STRING) ? NA_REAL : AG[pg[i]];
                  }
                }
                break;
              }
              case INTSXP:
              case LGLSXP:
              {
                int *pout = INTEGER(out), *pAG = INTEGER(xAG);
                if(nog) {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, e = s + row, AGj = pAG[j];
                    for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_STRING) ? NA_INTEGER : AGj;
                  }
                } else {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, *AG = pAG + j * ng - 1;
                    for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_STRING) ? NA_INTEGER : AG[pg[i]];
                  }
                }
                break;
              }
              case STRSXP:
              {
                SEXP *pout = STRING_PTR(out), *pAG = STRING_PTR(xAG);
                if(nog) {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row, e = s + row;
                    SEXP AGj = pAG[j];
                    for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_STRING) ? NA_STRING : AGj;
                  }
                } else {
                  for(int j = 0; j != col; ++j) {
                    int s = j * row;
                    SEXP *AG = pAG + j * ng - 1;
                    for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_STRING) ? NA_STRING : AG[pg[i]];
                  }
                }
                break;
              }
              default: error("Not supported SEXP type!");
            }
            break;
          }
          default:
            error("Not supported SEXP type!");
        }

      }

      if(set == 0) {
        if(isObject(xAG)) SHALLOW_DUPLICATE_ATTRIB(out, xAG);
        else if(!isObject(x) || (tx == txAG && !isFactor(x))) SHALLOW_DUPLICATE_ATTRIB(out, x);
        else {
          SHALLOW_DUPLICATE_ATTRIB(out, x);
          classgets(out, R_NilValue); // OK !
          setAttrib(out, R_LevelsSymbol, R_NilValue);
        }
        UNPROTECT(1);
      }

      return out;

    } else { // ret == 0

      if(ret != 0) error("Unknown Transformation!");

      SEXP out = set ? x : PROTECT(allocVector(tx, row * col));

      switch(tx) {
      case REALSXP:
      {
        double *px = REAL(x), *pout = REAL(out);
        switch(txAG) {
        case REALSXP:
        {
          double *pAG = REAL(xAG);
          if(nog) {
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (ISNAN(px[i])) ? AGj : px[i];
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = pAG + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (ISNAN(px[i + s])) ? AG[pg[i]] : px[i + s];
            }
          }
          break;
        }
        case INTSXP:
        case LGLSXP:
        {
          int *pAG = INTEGER(xAG);
          if(nog) {
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              double AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (ISNAN(px[i])) ? AGj : px[i];
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = pAG + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (ISNAN(px[i + s])) ? AG[pg[i]] : px[i + s];
            }
          }
          break;
        }
        case STRSXP: error("Cannot replace missing values in double with a string");
        default: error("Not supported SEXP type!");
        }
        break;
      }
      case INTSXP:
      case LGLSXP:
      {
        int *px = INTEGER(x), *pout = INTEGER(out);
        switch(txAG) {
        case REALSXP:
        {
          double *pAG = REAL(xAG);
          if(nog) {
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_INTEGER) ? AGj : px[i];
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              double *AG = pAG + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_INTEGER) ? AG[pg[i]] : px[i + s];
            }
          }
          break;
        }
        case INTSXP:
        case LGLSXP:
        {
          int *pAG = INTEGER(xAG);
          if(nog) {
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row, AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_INTEGER) ? AGj : px[i];
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row, *AG = pAG + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_INTEGER) ? AG[pg[i]] : px[i + s];
            }
          }
          break;
        }
        case STRSXP: error("Cannot replace missing values in integer with a string");
        default: error("Not supported SEXP type!");
        }
        break;
      }
      case STRSXP:
      {
        SEXP *px = STRING_PTR(x), *pout = STRING_PTR(out);
        switch(txAG) {
        case REALSXP:
        case INTSXP:
        case LGLSXP: error("Cannot replace missing values in string with numeric data");
        case STRSXP:
        {
          SEXP *pAG = STRING_PTR(xAG);
          if(nog) {
            for(int j = 0; j != col; ++j) {
              int s = j * row, e = s + row;
              SEXP AGj = pAG[j];
              for(int i = s; i != e; ++i) pout[i] = (px[i] == NA_STRING) ? AGj : px[i];
            }
          } else {
            for(int j = 0; j != col; ++j) {
              int s = j * row;
              SEXP *AG = pAG + j * ng - 1;
              for(int i = 0; i != row; ++i) pout[i + s] = (px[i + s] == NA_STRING) ? AG[pg[i]] : px[i + s];
            }
          }
          break;
        }
        default: error("Not supported SEXP type!");
        }
        break;
      }
      default:
        error("Not supported SEXP type!");
      }

      if(set == 0) {
        SHALLOW_DUPLICATE_ATTRIB(out, x);
        UNPROTECT(1);
      }

      return out;

    }
  }

  // ret > 2

  int nprotect = 0;
  SEXP out = set ? x : PROTECT(allocVector(REALSXP, row * col));
  double *pAG;

  if(txAG != REALSXP) {
    if(txAG != INTSXP && txAG != LGLSXP) error("STATS needs to be double, integer or logical");
    SEXP xxAG = PROTECT(coerceVector(xAG, REALSXP)); ++nprotect;
    pAG = REAL(xxAG);
  } else pAG = REAL(xAG);

#define MATNUMTRALOOP                                                              \
  switch(ret) {                                                                    \
    case 3: {                                                                      \
      if(nog) {                                                                    \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row, e = s + row;                                              \
        double AGj = pAG[j];                                                       \
        for(int i = s; i != e; ++i) pout[i] = px[i] - AGj;                         \
      }                                                                            \
    } else {                                                                       \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row;                                                           \
        double *AG = pAG + j * ng - 1;                                             \
        for(int i = 0; i != row; ++i) pout[i + s] = px[i + s] - AG[pg[i]];         \
      }                                                                            \
    }                                                                              \
    break;                                                                         \
    }                                                                              \
    case 4: {                                                                      \
      if(nog) error("This transformation can only be computed with groups!");      \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row, n = 0;                                                    \
        long double OM = 0;                                                        \
        double *AG = pAG + j * ng - 1;                                             \
        for(int i = 0; i != row; ++i) {                                            \
          if(ISNAN(px[i + s])) pout[i + s] = px[i + s];                            \
          else {                                                                   \
            pout[i + s] = px[i + s] - AG[pg[i]];                                   \
            if(ISNAN(AG[pg[i]])) continue;                                         \
            OM += AG[pg[i]];                                                       \
            ++n;                                                                   \
          }                                                                        \
        }                                                                          \
        OM /= n;                                                                   \
        double OMD = (double)OM;                                                   \
        for(int i = row; i--; ) pout[i + s] += OMD;                                \
      }                                                                            \
      break;                                                                       \
    }                                                                              \
    case 5: {                                                                      \
      if(nog) {                                                                    \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row, e = s + row;                                              \
        double AGj = 1 / pAG[j];                                                   \
        for(int i = s; i != e; ++i) pout[i] = px[i] * AGj;                         \
      }                                                                            \
    } else {                                                                       \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row;                                                           \
        double *AG = pAG + j * ng - 1;                                             \
        for(int i = 0; i != row; ++i) pout[i + s] = px[i + s] * (1 / AG[pg[i]]);   \
      }                                                                            \
    }                                                                              \
    break;                                                                         \
    }                                                                              \
    case 6: {                                                                      \
      if(nog) {                                                                    \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row, e = s + row;                                              \
        double AGj = 100 / pAG[j];                                                 \
        for(int i = s; i != e; ++i) pout[i] = px[i] * AGj;                         \
      }                                                                            \
    } else {                                                                       \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row;                                                           \
        double *AG = pAG + j * ng - 1;                                             \
        for(int i = 0; i != row; ++i) pout[i + s] = px[i + s] * (100 / AG[pg[i]]); \
      }                                                                            \
    }                                                                              \
    break;                                                                         \
    }                                                                              \
    case 7: {                                                                      \
      if(nog) {                                                                    \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row, e = s + row;                                              \
        double AGj = pAG[j];                                                       \
        for(int i = s; i != e; ++i) pout[i] = px[i] + AGj;                         \
      }                                                                            \
    } else {                                                                       \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row;                                                           \
        double *AG = pAG + j * ng - 1;                                             \
        for(int i = 0; i != row; ++i) pout[i + s] = px[i + s] + AG[pg[i]];         \
      }                                                                            \
    }                                                                              \
    break;                                                                         \
    }                                                                              \
    case 8: {                                                                      \
      if(nog) {                                                                    \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row, e = s + row;                                              \
        double AGj = pAG[j];                                                       \
        for(int i = s; i != e; ++i) pout[i] = px[i] * AGj;                         \
      }                                                                            \
    } else {                                                                       \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row;                                                           \
        double *AG = pAG + j * ng - 1;                                             \
        for(int i = 0; i != row; ++i) pout[i + s] = px[i + s] * AG[pg[i]];         \
      }                                                                            \
    }                                                                              \
    break;                                                                         \
    }                                                                              \
    case 9: {                                                                      \
      if(nog) {                                                                    \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row, e = s + row;                                              \
        double AGj = pAG[j];                                                       \
        for(int i = s; i != e; ++i) pout[i] = modulus_impl(px[i], AGj);                   \
      }                                                                            \
    } else {                                                                       \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row;                                                           \
        double *AG = pAG + j * ng - 1;                                             \
        for(int i = 0; i != row; ++i) pout[i + s] = modulus_impl(px[i + s], AG[pg[i]]);   \
      }                                                                            \
    }                                                                              \
    break;                                                                         \
    }                                                                              \
    case 10: {                                                                     \
      if(nog) {                                                                    \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row, e = s + row;                                              \
        double AGj = pAG[j];                                                       \
        for(int i = s; i != e; ++i) pout[i] = remainder_impl(px[i], AGj);                \
      }                                                                            \
    } else {                                                                       \
      for(int j = 0; j != col; ++j) {                                              \
        int s = j * row;                                                           \
        double *AG = pAG + j * ng - 1;                                             \
        for(int i = 0; i != row; ++i) pout[i + s] = remainder_impl(px[i + s], AG[pg[i]]);\
      }                                                                            \
    }                                                                              \
    break;                                                                         \
    }                                                                              \
    default: error("Unknown Transformation");                                      \
  }

  switch(tx) {
    case REALSXP:
    {
      double *pout = REAL(out), *px = REAL(x);
      MATNUMTRALOOP
      break;
    }
    case INTSXP:
    case LGLSXP:
    {
      int *px = INTEGER(x);
      if(set) {
        int *pout = INTEGER(out);
        MATNUMTRALOOP
      } else {
        double *pout = REAL(out);
        MATNUMTRALOOP
      }
      break;
    }
    default: error("Not supported SEXP type!");
  }

  if(set == 0) {
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    UNPROTECT(nprotect + 1);
  } else if(nprotect > 0) UNPROTECT(nprotect);

  return out;
}





