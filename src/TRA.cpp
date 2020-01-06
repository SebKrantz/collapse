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
        if(gs != l) stop("length(g) must match nrow(x)");
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
        if(gs != l) stop("length(g) must match nrow(x)");
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
        if(gs != l) stop("length(g) must match nrow(x)");
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
        if(gs != l) stop("length(g) must match nrow(x)");
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
        if(gs != l) stop("length(g) must match nrow(x)");
        switch(ret) {
        case 3:
          for(int i = l; i--; ) out[i] = xx[i] - AG[g[i]-1];
          break;
        case 4:
          {
            long double OM = 0; // better precision !!
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
            out = out + (double)OM; // Fastest ??
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
            long double OM = 0; // gives better numeric precision !! (closer to W!!)
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
            colo = colo + (double)OM; // Fastest ??
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
      NumericVector AG = no_init_vector(l); // NULL; // gives compile warning !!
      if(TYPEOF(xAG) == VECSXP) {
        // AG = NumericVector(l); // stable now ??
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
        long double OM = 0;
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
        sgj = sgj + (double)OM; // Fastest !!
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
