#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fmeanCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0,
                       const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, bool narm = true) {
  int l = x.size();

  if (Rf_isNull(w)) { // No weights
    if (ng == 0) {
      if(narm) {
        int j = l-1, n = 1; // 1 because for-loop starts from 2
        // long double sum = x[j];
        double sum = x[j];
        while(std::isnan(sum) && j!=0) sum = x[--j];
        if(j != 0) for(int i = j; i--; ) {
          if(std::isnan(x[i])) continue;
          sum += x[i]; // Fastest ?
          ++n;
        }
        sum = sum/n;
        return NumericVector::create(sum); // :create((double)sum)
      } else {
        // long double sum = 0;
        double sum = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(x[i])) {
            sum = x[i];
            break;
          } else {
            sum += x[i];
          }
        }
        sum = sum/l;
        return NumericVector::create(sum); // :create((double)sum)
      }
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      if(narm) {
        NumericVector sum(ng, NA_REAL); // Other way ?
        IntegerVector n(ng, 1); // could also do no_init_vector and then add n[g[i]-1] = 1 in fir if condition... -> Nope, that is slower
        for(int i = l; i--; ) {
          if(!std::isnan(x[i])) { // faster way to code this ? -> Not Bad at all -> index for g[i]-1? -> Nope, no noticeable improvement
            if(std::isnan(sum[g[i]-1])) sum[g[i]-1] = x[i];
            else {
              sum[g[i]-1] += x[i];
              ++n[g[i]-1];
            }
          }
        }
        for(int i = ng; i--; ) sum[i] /= n[i]; // if(n[i] == 0) stop("group size of 0 encountered"); -> No check possible when initializing at 1
        DUPLICATE_ATTRIB(sum, x);
        return sum;
      } else {
        NumericVector sum(ng); // no_init_vector // good? -> yes, but not initializing is numerically unstable..
        int ngs = 0;
        if(Rf_isNull(gs)) {
          IntegerVector gsv(ng);
          if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
          for(int i = 0; i != l; ++i) {
            if(std::isnan(x[i])) {
              if(!std::isnan(sum[g[i]-1])) {
                sum[g[i]-1] = x[i];
                ++ngs;
                if(ngs == ng) break;
              }
            } else {
              sum[g[i]-1] += x[i];
              ++gsv[g[i]-1];
            }
          }
          for(int i = ng; i--; ) sum[i] /= gsv[i]; // Adding n takes twice as long,
        } else {
          IntegerVector gsv = gs;
          if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
          for(int i = 0; i != l; ++i) {
            if(std::isnan(x[i])) {
              if(!std::isnan(sum[g[i]-1])) {
                sum[g[i]-1] = x[i];
                ++ngs;
                if(ngs == ng) break;
              }
            } else {
              sum[g[i]-1] += x[i];
            }
          }
          for(int i = ng; i--; ) {
            if(gsv[i] == 0) stop("group size of 0 encountered");
            sum[i] /= gsv[i]; // This is good because adding n takes twice as long, if factor, supply gs = tabulate(f,nlevels(f))
          }
        }
        DUPLICATE_ATTRIB(sum, x);
        return sum;
      }
    }
  } else { // With weights
    NumericVector wg = w; // wg(w) Identical speed
    if(l != wg.size()) stop("length(w) must match length(x)");
    if (ng == 0) {
      if(narm) {
        int j = l-1; // 1 because for-loop starts from 2
        while((std::isnan(x[j]) || std::isnan(wg[j])) && j!=0) --j; // This does not make a difference in performance but is more parsimonious.
        // long double sum = x[j]*wg[j], sumw = wg[j];
        double sum = x[j]*wg[j], sumw = wg[j];
        if(j != 0) for(int i = j; i--; ) {
          if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
          sum += x[i]*wg[i]; // Fastest ??
          sumw += wg[i];
        }
        sum = sum/sumw;
        return NumericVector::create(sum); // :create((double)sum)
      } else {
        // long double sum = 0, sumw = 0;
        double sum = 0, sumw = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(x[i]) || std::isnan(wg[i])) { // good, check both ? -> yes
            sum = x[i]+wg[i];
            break;
          } else {
            sum += x[i]*wg[i];
            sumw += wg[i];
          }
        }
        sum = sum/sumw;
        return NumericVector::create(sum); // :create((double)sum)
      }
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      if(narm) {
        NumericVector sum(ng, NA_REAL); // Other way ? -> Nope, this is as good as it gets
        NumericVector sumw(ng); // = no_init_vector(ng); // no init works!! but gives valgrind issue
        for(int i = l; i--; ) {
          if(std::isnan(x[i]) || std::isnan(wg[i])) continue; // faster way to code this ? -> Not Bad at all -> index for g[i]-1? -> Nope, no noticeable improvement
          if(std::isnan(sum[g[i]-1])) {
            sum[g[i]-1] = x[i]*wg[i];
            sumw[g[i]-1] = wg[i];
          } else {
            sum[g[i]-1] += x[i]*wg[i];
            sumw[g[i]-1] += wg[i];
          }
        }
        sum = sum/sumw; // good ? better return sum/sumw? -> Nope, slower
        DUPLICATE_ATTRIB(sum, x);
        return sum;
      } else {
        NumericVector sum(ng), sumw(ng); // good? -> yes ! //  = no_init_vector// Not initializing numerically unstable
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(x[i]) || std::isnan(wg[i])) {
            if(!std::isnan(sum[g[i]-1])) {
              sum[g[i]-1] = sumw[g[i]-1] = x[i]+wg[i]; // or NA_REAL ? -> Nope, good
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            sum[g[i]-1] += x[i]*wg[i];
            sumw[g[i]-1] += wg[i];
          }
        }
        sum = sum/sumw;
        DUPLICATE_ATTRIB(sum, x);
        return sum;
      }
    }
  }
}





// [[Rcpp::export]]
SEXP fmeanmCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0, const SEXP& gs = R_NilValue,
               const SEXP& w = R_NilValue, bool narm = true, bool drop = true) {
  int l = x.nrow(), col = x.ncol();

  if(Rf_isNull(w)) { // No weights
    if(ng == 0) {
      NumericVector sum = no_init_vector(col); //  Initialize faster -> Nope
      if(narm) {
        for(int j = col; j--; ) { // Instead Am(j,_) you can use Am.row(j).
          NumericMatrix::ConstColumn column = x( _ , j);
          int k = l-1, nj = 1;
          // long double sumj = column[k];
          double sumj = column[k];
          while(std::isnan(sumj) && k!=0) sumj = column[--k];
          if(k != 0) for(int i = k; i--; ) {
            if(std::isnan(column[i])) continue;
            sumj += column[i];
            ++nj;
          }
          sumj = sumj/nj;
          sum[j] = sumj; // (double)sumj;
        }
      } else {
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j);
          // long double sumj = 0;
          double sumj = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(column[i])) {
              sumj = column[i];
              break;
            } else {
              sumj += column[i];
            }
          }
          sumj = sumj/l;
          sum[j] = sumj; // (double)sumj;
        }
      }
      if(drop) sum.attr("names") = colnames(x); // Slight speed loss 31 to 34 milliseconds on WDIM, but doing it in R not faster
      else {
        sum.attr("dim") = Dimension(1, col);
        // sum.attr("dimnames") = List::create(R_NilValue,colnames(x));
        colnames(sum) = colnames(x); // NEW! faster than R ? -> yes, good
      }
      return sum;
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      if(narm) {
        NumericMatrix sum = no_init_matrix(ng, col);
        std::fill(sum.begin(), sum.end(), NA_REAL); // fastest ? or create vector and declare as matrix ?
        // NumericVector sumt(ng*col, NA_REAL); // A tiny speed gain, but not much !! Same memory efficiency
        // sumt.attr("dim") = Dimension(ng, col);
        // NumericMatrix sum = as<NumericMatrix>(sumt);
        IntegerVector nj(ng); // = no_init_vector(ng); // better for valgrind
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j);
          NumericMatrix::Column sumj = sum( _ , j);
           // int nj[ng]; // Numerically stable and faster and more memory efficient than before
          for(int i = l; i--; ) {
            if(!std::isnan(column[i])) {
              if(std::isnan(sumj[g[i]-1])) {
                sumj[g[i]-1] = column[i];
                nj[g[i]-1] = 1;
              } else {
                sumj[g[i]-1] += column[i];
                ++nj[g[i]-1];
              }
            }
          }
          for(int i = ng; i--; ) sumj[i] /= nj[i]; // if(gsv[i] == 0) stop("group size of 0 encountered"); cant check when not initializing
        }
        colnames(sum) = colnames(x);  // extremely efficient
        return sum;
      } else {
        NumericMatrix sum(ng, col); // no init numerically unstable
        if(Rf_isNull(gs)) {
          // int gsv[ng], memsize = sizeof(int)*ng;
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            NumericMatrix::Column sumj = sum( _ , j);
            // memset(gsv, 0, memsize); // still a tiny bit faster than std::vector, but both have the same memory efficiency
            std::vector<int> gsv(ng);
            int ngs = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(column[i])) {
                if(!std::isnan(sumj[g[i]-1])) {
                  sumj[g[i]-1] = column[i];
                  ++ngs;
                  if(ngs == ng) break;
                }
              } else {
                sumj[g[i]-1] += column[i];
                ++gsv[g[i]-1];
              }
            }
            for(int i = ng; i--; ) sumj[i] /= gsv[i];
          }
        } else {
          IntegerVector gsv = gs;
          if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            NumericMatrix::Column sumj = sum( _ , j);
            int ngs = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(column[i])) {
                if(!std::isnan(sumj[g[i]-1])) {
                  sumj[g[i]-1] = column[i];
                  ++ngs;
                  if(ngs == ng) break;
                }
              } else {
                sumj[g[i]-1] += column[i];
              }
            }
            for(int i = ng; i--; ) {
              if(gsv[i] == 0) stop("group size of 0 encountered");
              sumj[i] /= gsv[i];
            }
          }
        }
        colnames(sum) = colnames(x);  // quite efficient
        return sum;
      }
    }
  } else { // With weights
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match nrow(X)");
    if(ng == 0) {
      NumericVector sum = no_init_vector(col); // Initialize faster -> Nope
      if(narm) {
        for(int j = col; j--; ) { // Instead Am(j,_) you can use Am.row(j).
          NumericMatrix::ConstColumn column = x( _ , j);
          int k = l-1;
          while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
          // long double sumj = column[k]*wg[k], sumwj = wg[k];
          double sumj = column[k]*wg[k], sumwj = wg[k];
          if(k != 0) for(int i = k; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
            sumj += column[i]*wg[i];
            sumwj += wg[i];
          }
          sumj = sumj/sumwj;
          sum[j] = sumj; // (double)sumj;
        }
      } else {
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j);
          // long double sumj = 0, sumwj = 0;
          double sumj = 0, sumwj = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) {
              sumj = column[i]+wg[i];
              break;
            } else {
              sumj += column[i]*wg[i];
              sumwj += wg[i];
            }
          }
          sumj = sumj/sumwj;
          sum[j] = sumj; // (double)sumj;
        }
      }
      if(drop) sum.attr("names") = colnames(x); // Slight speed loss 31 to 34 milliseconds on WDIM, but doing it in R not faster
      else {
        sum.attr("dim") = Dimension(1, col);
        // sum.attr("dimnames") = List::create(R_NilValue,colnames(x));
        colnames(sum) = colnames(x); // NEW! faster than R ? -> yes, good
      }
      return sum;
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      if(narm) {
        NumericMatrix sum = no_init_matrix(ng, col);
        std::fill(sum.begin(), sum.end(), NA_REAL);
        // NumericMatrix sumw = no_init_matrix(ng, col); // Numerically stable ? -> Yes
        NumericVector sumwj(ng); // = no_init_vector(ng);
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j);
          NumericMatrix::Column sumj = sum( _ , j);
          // NumericVector sumwj = no_init_vector(ng);
          // NumericMatrix::Column sumwj = sumw( _ , j);
          // double sumwj[ng]; // Numerically stable, Slightly faster and a lot more memory efficient (but long double is a lot slower)
          for(int i = l; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
            if(std::isnan(sumj[g[i]-1])) {
              sumj[g[i]-1] = column[i]*wg[i];
              sumwj[g[i]-1] = wg[i];
            } else {
              sumj[g[i]-1] += column[i]*wg[i];
              sumwj[g[i]-1] += wg[i];
            }
          }
          for(int i = ng; i--; ) sumj[i] /= sumwj[i];
          // sumj = sumj/sumwj; // This gives error because sumj is matrix column !
        }
        colnames(sum) = colnames(x);  // quite efficient
        return sum;
      } else {
        NumericMatrix sum(ng, col); // no init numerically unstable
        // NumericMatrix sumw(ng, col); // also here ? -> Nope
        // double sumwj[ng]; // Also a bit faster and a lot more memory efficient
        // int memsize = sizeof(double)*ng;
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j);
          NumericMatrix::Column sumj = sum( _ , j);
          // NumericMatrix::Column sumwj = sumw( _ , j);
          std::vector<double> sumwj(ng); // memset(sumwj, 0, memsize);
          int ngs = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) {
              if(!std::isnan(sumj[g[i]-1])) {
                sumj[g[i]-1] = sumwj[g[i]-1] = column[i]+wg[i]; // or NA_REAL ? -> Nope, good
                ++ngs;
                if(ngs == ng) break;
              }
            } else {
              sumj[g[i]-1] += column[i]*wg[i];
              sumwj[g[i]-1] += wg[i];
            }
          }
          for(int i = ng; i--; ) sumj[i] /= sumwj[i];
          // sumj = sumj/sumwj; // This gives erriir because sumj is matrix column
        }
        colnames(sum) = colnames(x);  // quite efficient
        return sum;
      }
    }
  }
}





// [[Rcpp::export]]
SEXP fmeanlCpp(const List& x, int ng = 0, const IntegerVector& g = 0, const SEXP& gs = R_NilValue,
               const SEXP& w = R_NilValue, bool narm = true, bool drop = true) {
  int l = x.size();

  if(Rf_isNull(w)) { // No weights
    if(ng == 0) {
      NumericVector sum(l); // not initializing not faster WIth NWDI (35 instead of 32 milliseconds)
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          int k = column.size()-1, ni = 1;
          // long double sumi = column[k]; // long double gives 45 instead of 35 milliseconds
          double sumi = column[k];
          while(std::isnan(sumi) && k!=0) sumi = column[--k];
          if(k != 0) for(int i = k; i--; ) {
            if(std::isnan(column[i])) continue;
            sumi += column[i];
            ++ni;
          }
          sumi = sumi/ni;
          sum[j] = sumi; // (double)sumi;
        }
      } else {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          // long double sumi = 0;
          double sumi = 0;
          int row = column.size();
          for(int i = 0; i != row; ++i) {
            if(std::isnan(column[i])) {
              sumi = column[i];
              break;
            } else {
              sumi += column[i];
            }
          }
          sumi = sumi/row;
          sum[j] = sumi; // (double)sumi;
        }
      }
      if(drop) {
        sum.attr("names") = x.attr("names");
        return sum;
      } else {
        List out(l);
        for(int j = l; j--; ) {
          out[j] = sum[j];
          SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
        }
        DUPLICATE_ATTRIB(out, x);
        out.attr("row.names") = 1;
        return out;
      }
    } else { // With groups
      List sum(l);
      int gss = g.size();
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          NumericVector sumj(ng, NA_REAL);
          // IntegerVector nj(ng, 1); // faster than no_init_vector ? -> good, cannot divide by interger 0, also else numerically unstable and no speed loss
          std::vector<int>  nj(ng, 1); // better memory allocation, and nearly same speed as integer array -> doesn't work because sets all byte to 1 -> https://stackoverflow.com/questions/14761015/memset-an-array-to-1
          for(int i = gss; i--; ) {
            if(!std::isnan(column[i])) { // faster way to code this ? -> Not Bad at all, 54.. millisec for WDIM
              if(std::isnan(sumj[g[i]-1])) sumj[g[i]-1] = column[i];
              else {
                sumj[g[i]-1] += column[i];
                ++nj[g[i]-1];
              }
            }
          }
          for(int i = ng; i--; ) sumj[i] /= nj[i];
          SHALLOW_DUPLICATE_ATTRIB(sumj, column);
          sum[j] = sumj;
        }
      } else {
        if(Rf_isNull(gs)) {
          // int gsv[ng], memsize = sizeof(int)*ng;
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector sumj(ng); //  = no_init_vector //  Not initializing seems to be numerically unstable
            std::vector<int> gsv(ng); // memset(gsv, 0, memsize);
            int ngs = 0;
            for(int i = 0; i != gss; ++i) {
              if(std::isnan(column[i])) {
                if(!std::isnan(sumj[g[i]-1])) {
                  sumj[g[i]-1] = column[i];
                  ++ngs;
                  if(ngs == ng) break;
                }
              } else {
                sumj[g[i]-1] += column[i];
                ++gsv[g[i]-1];
              }
            }
            for(int i = ng; i--; ) sumj[i] /= gsv[i];
            SHALLOW_DUPLICATE_ATTRIB(sumj, column);
            sum[j] = sumj;
          }
        } else {
          IntegerVector gsv = gs;
          if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector sumj(ng); //  = no_init_vector //  Not initializing seems to be numerically unstable
            int ngs = 0;
            for(int i = 0; i != gss; ++i) {
              if(std::isnan(column[i])) {
                if(!std::isnan(sumj[g[i]-1])) {
                  sumj[g[i]-1] = column[i];
                  ++ngs;
                  if(ngs == ng) break;
                }
              } else {
                sumj[g[i]-1] += column[i];
              }
            }
            for(int i = ng; i--; ) {
              if(gsv[i] == 0) stop("group size of 0 encountered");
              sumj[i] /= gsv[i];
            }
            SHALLOW_DUPLICATE_ATTRIB(sumj, column);
            sum[j] = sumj;
          }
        }
      }
      DUPLICATE_ATTRIB(sum, x);
      sum.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng);
      return sum;
    }
  } else { // With weights
    NumericVector wg = w; // wg(w) //  No speed loss ? -> Yes, and possibly even faster
    int wgs = wg.size();
    if (ng == 0) {
      NumericVector sum(l); // not initializing not faster WIth NWDI (35 instead of 32 milliseconds)
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          if(column.size() != wgs) stop("length(w) must match nrow(X)"); // Really necessary ?
          int k = wgs-1;
          while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
          // long double sumi = column[k]*wg[k], sumwi = wg[k];
          double sumi = column[k]*wg[k], sumwi = wg[k];
          if(k != 0) for(int i = k; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
            sumi += column[i]*wg[i];
            sumwi += wg[i];
          }
          sumi = sumi/sumwi;
          sum[j] = sumi; // (double)sumi;
        }
      } else {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          if(column.size() != wgs) stop("length(w) must match nrow(X)"); // Really necessary ?
          // long double sumi = 0, sumwi = 0;
          double sumi = 0, sumwi = 0;
          for(int i = 0; i != wgs; ++i) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) {
              sumi = column[i]+wg[i];
              break;
            } else {
              sumi += column[i]*wg[i];
              sumwi += wg[i];
            }
          }
          sumi = sumi/sumwi;
          sum[j] = sumi; // (double)sumi;
        }
      }
      if(drop) {
        sum.attr("names") = x.attr("names");
        return sum;
      } else {
        List out(l);
        for(int j = l; j--; ) {
          out[j] = sum[j];
          SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
        }
        DUPLICATE_ATTRIB(out, x);
        out.attr("row.names") = 1;
        return out;
      }
    } else { // With groups
      List sum(l);
      int gss = g.size();
      if(wgs != gss) stop("length(w) must match length(g)");
      if(narm) {
        NumericVector sumwj(ng);  // = no_init_vector(ng); // stable and faster ?
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          NumericVector sumj(ng, NA_REAL);
          // no_init_vector is faster and stable (you only divide by it every round)
          // double sumwj[ng];
          for(int i = gss; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
            if(std::isnan(sumj[g[i]-1])) {
              sumj[g[i]-1] = column[i]*wg[i];
              sumwj[g[i]-1] = wg[i];
            } else {
              sumj[g[i]-1] += column[i]*wg[i];
              sumwj[g[i]-1] += wg[i];
            }
          }
          sumj = sumj/sumwj;
          // for(int i = ng; i--; ) sumj[i] /= sumwj[i];
          SHALLOW_DUPLICATE_ATTRIB(sumj, column);
          sum[j] = sumj;
        }
      } else {
        // double sumwj[ng];
        // int memsize = sizeof(double)*ng;
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          if(gss != column.size()) stop("length(g) must match nrow(X)");
          NumericVector sumj(ng), sumwj(ng); //  = no_init_vector //  Not initializing seems to be numerically unstable
          // NumericVector sumwj(ng); // Also here not initializing is numerically unstable
          // memset(sumwj, 0, memsize);
          int ngs = 0;
          for(int i = 0; i != gss; ++i) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) {
              if(!std::isnan(sumj[g[i]-1])) {
                sumj[g[i]-1] = sumwj[g[i]-1] = column[i]+wg[i]; // or NA_REAL ? -> Nope, good
                ++ngs;
                if(ngs == ng) break;
              }
            } else {
              sumj[g[i]-1] += column[i]*wg[i];
              sumwj[g[i]-1] += wg[i];
            }
          }
          sumj = sumj/sumwj;
          // for(int i = ng; i--; ) sumj[i] /= sumwj[i];
          SHALLOW_DUPLICATE_ATTRIB(sumj, column);
          sum[j] = sumj;
        }
      }
      DUPLICATE_ATTRIB(sum, x);
      sum.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng);
      return sum;
    }
  }
}
