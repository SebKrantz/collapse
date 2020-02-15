#include <Rcpp.h>
using namespace Rcpp;

// Note: All comments are in fvar.cpp

// [[Rcpp::export]]
NumericVector fscaleCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0,
                        const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue,
                        bool narm = true, bool stable_algo = true) {
  int l = x.size();
  NumericVector out = no_init_vector(l);

  if(stable_algo) { // WELFORDS ONLINE METHOD ---------------------------------------------------------
    if(Rf_isNull(w)) { // No weights
      if (ng == 0) {
        if(narm) {
          int j = l-1;
          // double n = 0;
          // long double mean = 0, d1 = 0, M2 = 0;
          double n = 0, mean = 0, d1 = 0, M2 = 0;
          while(std::isnan(x[j]) && j!=0) --j;
          if(j != 0) {
            for(int i = j+1; i--; ) {
              if(std::isnan(x[i])) continue;
              d1 = x[i]-mean;
              mean += d1 * (1 / ++n);
              M2 += d1*(x[i]-mean);
            }
            M2 = 1/sqrt(M2/(n-1));
            if(std::isnan(M2)) {
              std::fill(out.begin(), out.end(), NA_REAL);
            } else {
              out = (x-mean)*M2;
              // for(int i = 0; i != l; ++i) out[i] = (x[i]-mean)*M2; // For some weird reason this is a lot slower than the above vectorized (2.5 milliseconds vs. 80 microseconds on a WDI series, the grouped version is much faster !!)
            }
          } else {
            std::fill(out.begin(), out.end(), NA_REAL);
          }
        } else {
          // double n = 0;
          // long double mean = 0, d1 = 0, M2 = 0;
          double n = 0, mean = 0, d1 = 0, M2 = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(x[i])) {
              std::fill(out.begin(), out.end(), NA_REAL);
              return out;
            } else {
              d1 = x[i]-mean;
              mean += d1*(1 / ++n);
              M2 += d1*(x[i]-mean);
            }
          }
          M2 = 1/sqrt(M2/(l-1));
          if(std::isnan(M2)) {
            std::fill(out.begin(), out.end(), NA_REAL);
          } else {
            out = (x-mean)*M2;
            // for(int i = 0; i != l; ++i) out[i] = (x[i]-mean)*M2;
          }
        }
      } else { // with groups
        if(g.size() != l) stop("length(g) must match nrow(X)");
        // long double d1 = 0;
        double d1 = 0;
        if(narm) {
          NumericVector M2(ng, NA_REAL);
          NumericVector mean = no_init_vector(ng);
          NumericVector n(ng, 1.0);
          for(int i = l; i--; ) {
            if(std::isnan(x[i])) continue;
            if(std::isnan(M2[g[i]-1])) {
              mean[g[i]-1] = x[i];
              M2[g[i]-1] = 0;
            } else {
              d1 = x[i]-mean[g[i]-1];
              mean[g[i]-1] += d1 * (1 / ++n[g[i]-1]);
              M2[g[i]-1] += d1*(x[i]-mean[g[i]-1]);
            }
          }
          for(int i = ng; i--; ) if(!std::isnan(M2[i])) M2[i] = 1/sqrt(M2[i]/(n[i]-1));
          for(int i = 0; i != l; ++i) out[i] = (x[i]-mean[g[i]-1])*M2[g[i]-1];
        } else {
          NumericVector M2(ng), mean(ng), n(ng);
          int ngs = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(M2[g[i]-1])) continue;
            if(std::isnan(x[i])) {
              M2[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) {
                std::fill(out.begin(), out.end(), NA_REAL);
                return out;
              }
            } else {
              d1 = x[i]-mean[g[i]-1];
              mean[g[i]-1] += d1 * (1 / ++n[g[i]-1]);
              M2[g[i]-1] += d1*(x[i]-mean[g[i]-1]);
            }
          }
          for(int i = ng; i--; ) if(!std::isnan(M2[i])) M2[i] = 1/sqrt(M2[i]/(n[i]-1));
          for(int i = 0; i != l; ++i) out[i] = (x[i]-mean[g[i]-1])*M2[g[i]-1];
        }
      }
    } else { // With weights
      NumericVector wg = w;
      if(l != wg.size()) stop("length(w) must match length(x)");
      if (ng == 0) {
        if(narm) {
          int j = l-1;
          // long double sumw = 0, mean = 0, M2 = 0, d1 = 0;
          double sumw = 0, mean = 0, M2 = 0, d1 = 0;
          while((std::isnan(x[j]) || std::isnan(wg[j])) && j!=0) --j;
          if(j != 0) {
            for(int i = j+1; i--; ) {
              if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
              sumw += wg[i];
              d1 = x[i] - mean;
              mean += d1 * (wg[i] / sumw);
              M2 += wg[i] * d1 * (x[i] - mean);
            }
            M2 = 1/sqrt(M2/(sumw-1));
            if(std::isnan(M2)) {
              std::fill(out.begin(), out.end(), NA_REAL);
            } else {
              out = (x-mean)*M2;
              // for(int i = 0; i != l; ++i) out[i] = (x[i]-mean)*M2;
            }
          } else {
            std::fill(out.begin(), out.end(), NA_REAL);
          }
        } else {
          // long double sumw = 0, mean = 0, M2 = 0, d1 = 0;
          double sumw = 0, mean = 0, M2 = 0, d1 = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(x[i]) || std::isnan(wg[i])) {
              std::fill(out.begin(), out.end(), NA_REAL);
              return out;
            } else {
              sumw += wg[i];
              d1 = x[i] - mean;
              mean += d1 * (wg[i] / sumw);
              M2 += wg[i] * d1 * (x[i] - mean);
            }
          }
          M2 = 1/sqrt(M2/(sumw-1));
          if(std::isnan(M2)) {
            std::fill(out.begin(), out.end(), NA_REAL);
          } else {
            out = (x-mean)*M2;
            // for(int i = 0; i != l; ++i) out[i] = (x[i]-mean)*M2;
          }
        }
      } else { // with groups
        if(g.size() != l) stop("length(g) must match nrow(X)");
        // long double d1 = 0;
        double d1 = 0;
        if(narm) {
          NumericVector M2(ng, NA_REAL);
          NumericVector sumw = no_init_vector(ng);
          NumericVector mean = no_init_vector(ng);
          for(int i = l; i--; ) {
            if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
            if(std::isnan(M2[g[i]-1])) {
              sumw[g[i]-1] = wg[i];
              mean[g[i]-1] = x[i];
              M2[g[i]-1] = 0;
            } else {
              sumw[g[i]-1] += wg[i];
              d1 = x[i] - mean[g[i]-1];
              mean[g[i]-1] += d1 * (wg[i] / sumw[g[i]-1]);
              M2[g[i]-1] += wg[i] * d1 * (x[i] - mean[g[i]-1]);
            }
          }
          for(int i = ng; i--; ) if(!std::isnan(M2[i])) M2[i] = 1/sqrt(M2[i]/(sumw[i]-1));
          for(int i = 0; i != l; ++i) out[i] = (x[i]-mean[g[i]-1])*M2[g[i]-1];
        } else {
          NumericVector M2(ng), sumw(ng), mean(ng);
          int ngs = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(M2[g[i]-1])) continue;
            if(std::isnan(x[i]) || std::isnan(wg[i])) {
              M2[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) {
                std::fill(out.begin(), out.end(), NA_REAL);
                return out;
              }
            } else {
              sumw[g[i]-1] += wg[i];
              d1 = x[i] - mean[g[i]-1];
              mean[g[i]-1] += d1 * (wg[i] / sumw[g[i]-1]);
              M2[g[i]-1] += wg[i] * d1 * (x[i] - mean[g[i]-1]);
            }
          }
          for(int i = ng; i--; ) if(!std::isnan(M2[i])) M2[i] = 1/sqrt(M2[i]/(sumw[i]-1));
          for(int i = 0; i != l; ++i) out[i] = (x[i]-mean[g[i]-1])*M2[g[i]-1];
        }
      }
    }

  } else { // ONE-PASS METHOD ---------------------------------------------------------
    if(Rf_isNull(w)) { // No weights
      if (ng == 0) {
        if(narm) {
          int j = l-1, n = 1;
          long double sum = x[j], sq_sum;
          while(std::isnan(sum) && j!=0) sum = x[--j];
          sq_sum = sum*sum;
          if(j != 0) {
            for(int i = j; i--; ) {
              if(std::isnan(x[i])) continue;
              sum += x[i];
              sq_sum += pow(x[i],2);
              ++n;
            }
            sum /= n;
            sq_sum = 1/sqrt((sq_sum-pow(sum,2)*n)/(n-1));
            if(std::isnan(sq_sum)) {
              std::fill(out.begin(), out.end(), NA_REAL);
            } else {
              out = (x-sum)*sq_sum;
              // for(int i = 0; i != l; ++i) out[i] = (x[i]-sum)*sq_sum;
            }
          } else {
            std::fill(out.begin(), out.end(), NA_REAL);
          }
        } else {
          long double sum = 0, sq_sum = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(x[i])) {
              std::fill(out.begin(), out.end(), NA_REAL);
              return out;
            } else {
              sum += x[i];
              sq_sum += pow(x[i],2);
            }
          }
          sum /= l;
          sq_sum = 1/sqrt((sq_sum-pow(sum,2)*l)/(l-1));
          if(std::isnan(sq_sum)) {
            std::fill(out.begin(), out.end(), NA_REAL);
          } else {
            out = (x-sum)*sq_sum;
            // for(int i = 0; i != l; ++i) out[i] = (x[i]-sum)*sq_sum;
          }
        }
      } else { // with groups
        if(g.size() != l) stop("length(g) must match nrow(X)");
        if(narm) {
          NumericVector sq_sum(ng, NA_REAL);
          NumericVector sum = no_init_vector(ng);
          IntegerVector n(ng, 1);
          for(int i = l; i--; ) {
            if(std::isnan(x[i])) continue;
            if(std::isnan(sq_sum[g[i]-1])) {
              sum[g[i]-1] = x[i];
              sq_sum[g[i]-1] = pow(x[i],2);
            } else {
              sum[g[i]-1] += x[i];
              sq_sum[g[i]-1] += pow(x[i],2);
              ++n[g[i]-1];
            }
          }
          for(int i = ng; i--; ) {
            sum[i] /= n[i];
            if(std::isnan(sq_sum[i])) continue;
            sq_sum[i] = 1/sqrt((sq_sum[i] - pow(sum[i],2)*n[i])/(n[i]-1));
          }
          for(int i = 0; i != l; ++i) out[i] = (x[i]-sum[g[i]-1])*sq_sum[g[i]-1];
        } else {
          NumericVector sq_sum(ng), sum(ng);
          IntegerVector gsv = no_init_vector(ng); // NULL; gives compiler warning
          int ngs = 0;
          if(Rf_isNull(gs)) {
            // gsv = IntegerVector(ng);
            std::fill(gsv.begin(), gsv.end(), 0);
            for(int i = 0; i != l; ++i) {
              if(std::isnan(x[i])) {
                if(std::isnan(sq_sum[g[i]-1])) continue;
                sq_sum[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) {
                  std::fill(out.begin(), out.end(), NA_REAL);
                  return out;
                }
              } else {
                sum[g[i]-1] += x[i];
                sq_sum[g[i]-1] += pow(x[i],2);
                ++gsv[g[i]-1];
              }
            }
          } else {
            gsv = gs;
            if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
            for(int i = 0; i != l; ++i) {
              if(std::isnan(x[i])) {
                if(std::isnan(sq_sum[g[i]-1])) continue;
                sq_sum[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) {
                  std::fill(out.begin(), out.end(), NA_REAL);
                  return out;
                }
              } else {
                sum[g[i]-1] += x[i];
                sq_sum[g[i]-1] += pow(x[i],2);
              }
            }
          }
          for(int i = ng; i--; ) {
            sum[i] /= gsv[i];
            if(std::isnan(sq_sum[i])) continue;
            sq_sum[i] = 1/sqrt((sq_sum[i] - pow(sum[i],2)*gsv[i])/(gsv[i]-1));
          }
          for(int i = 0; i != l; ++i) out[i] = (x[i]-sum[g[i]-1])*sq_sum[g[i]-1];
        }
      }
    } else { // With weights
      NumericVector wg = w;
      if(l != wg.size()) stop("length(w) must match length(x)");
      if (ng == 0) {
        if(narm) {
          int j = l-1;
          while((std::isnan(x[j]) || std::isnan(wg[j])) && j!=0) --j;
          long double sumw = wg[j], sum = x[j]*sumw, sq_sum = sum*x[j];
          if(j != 0) {
            for(int i = j; i--; ) {
              if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
              sum += x[i]*wg[i];
              sumw += wg[i];
              sq_sum += pow(x[i],2)*wg[i];
            }
            sum /= sumw;
            sq_sum = 1/sqrt((sq_sum - pow(sum,2)*sumw)/(sumw-1));
            if(std::isnan(sq_sum)) {
              std::fill(out.begin(), out.end(), NA_REAL);
            } else {
              out = (x-sum)*sq_sum;
              // for(int i = 0; i != l; ++i) out[i] = (x[i]-sum)*sq_sum;
            }
          } else {
            std::fill(out.begin(), out.end(), NA_REAL);
          }
        } else {
          long double sum = 0, sumw = 0, sq_sum = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(x[i]) || std::isnan(wg[i])) {
              std::fill(out.begin(), out.end(), NA_REAL);
              return out;
            } else {
              sum += x[i]*wg[i];
              sumw += wg[i];
              sq_sum += pow(x[i],2)*wg[i];
            }
          }
          sum /= sumw;
          sq_sum = 1/sqrt((sq_sum - pow(sum,2)*sumw)/(sumw-1));
          if(std::isnan(sq_sum)) {
            std::fill(out.begin(), out.end(), NA_REAL);
          } else {
            out = (x-sum)*sq_sum;
            // for(int i = 0; i != l; ++i) out[i] = (x[i]-sum)*sq_sum;
          }
        }
      } else { // with groups
        if(g.size() != l) stop("length(g) must match nrow(X)");
        if(narm) {
          NumericVector sq_sum(ng, NA_REAL);
          NumericVector sumw = no_init_vector(ng);
          NumericVector sum = no_init_vector(ng);
          for(int i = l; i--; ) {
            if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
            if(std::isnan(sq_sum[g[i]-1])) {
              sum[g[i]-1] = x[i]*wg[i];
              sumw[g[i]-1] = wg[i];
              sq_sum[g[i]-1] = pow(x[i],2)*wg[i];
            } else {
              sum[g[i]-1] += x[i]*wg[i];
              sumw[g[i]-1] += wg[i];
              sq_sum[g[i]-1] += pow(x[i],2)*wg[i];
            }
          }
          for(int i = ng; i--; ) {
            sum[i] /= sumw[i];
            if(std::isnan(sq_sum[i])) continue;
            sq_sum[i] = 1/sqrt((sq_sum[i] - pow(sum[i],2)*sumw[i])/(sumw[i]-1));
          }
          for(int i = 0; i != l; ++i) out[i] = (x[i]-sum[g[i]-1])*sq_sum[g[i]-1];
        } else {
          NumericVector sq_sum(ng), sumw(ng), sum(ng);
          int ngs = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(sq_sum[g[i]-1])) continue;
            if(std::isnan(x[i]) || std::isnan(wg[i])) {
              sq_sum[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) {
                std::fill(out.begin(), out.end(), NA_REAL);
                return out;
              }
            } else {
              sum[g[i]-1] += x[i]*wg[i];
              sumw[g[i]-1] += wg[i];
              sq_sum[g[i]-1] += pow(x[i],2)*wg[i];
            }
          }
          for(int i = ng; i--; ) {
            sum[i] /= sumw[i];
            if(std::isnan(sq_sum[i])) continue;
            sq_sum[i] = 1/sqrt((sq_sum[i] - pow(sum[i],2)*sumw[i])/(sumw[i]-1));
          }
          for(int i = 0; i != l; ++i) out[i] = (x[i]-sum[g[i]-1])*sq_sum[g[i]-1];
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  return out;
}



// [[Rcpp::export]]
NumericMatrix fscalemCpp(NumericMatrix x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue,
                         SEXP w = R_NilValue, bool narm = true, bool stable_algo = true) {

  int l = x.nrow(), col = x.ncol();
  NumericMatrix out = no_init_matrix(l, col);

  if(stable_algo) { // WELFORDS ONLINE METHOD -------------------------------------
    if (Rf_isNull(w)) { // No weights
      if(ng == 0) {
        if(narm) {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            int k = l-1;
            // double ni = 0;
            // long double meani = 0, d1i = 0, M2i = 0;
            double ni = 0, meani = 0, d1i = 0, M2i = 0;
            while(std::isnan(column[k]) && k!=0) --k;
            if(k != 0) {
              for(int i = k+1; i--; ) {
                if(std::isnan(column[i])) continue;
                d1i = column[i]-meani;
                meani += d1i * (1 / ++ni);
                M2i += d1i*(column[i]-meani);
              }
              M2i = 1/sqrt(M2i/(ni-1));
              if(std::isnan(M2i)) {
                std::fill(outj.begin(), outj.end(), NA_REAL);
              } else {
                outj = (column-meani)*M2i;
                // for(int i = 0; i != l; ++i) outj[i] = (column[i]-meani)*M2i; // For some weird reason this is a lot slower than the above vectorized (4.7 seconds vs. 0.14 seconds on WDIM, the grouped version is much faster !!)
              }
            } else {
              std::fill(outj.begin(), outj.end(), NA_REAL);
            }
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            // double ni = 0;
            // long double meani = 0, d1i = 0, M2i = 0;
            double ni = 0, meani = 0, d1i = 0, M2i = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(column[i])) {
                M2i = NA_REAL;
                break;
              } else {
                d1i = column[i]-meani;
                meani += d1i * (1 / ++ni);
                M2i += d1i*(column[i]-meani);
              }
            }
            M2i = 1/sqrt(M2i/(l-1));
            if(std::isnan(M2i)) {
              std::fill(outj.begin(), outj.end(), NA_REAL);
            } else {
              outj = (column-meani)*M2i;
              // for(int i = 0; i != l; ++i) outj[i] = (column[i]-meani)*M2i;
            }
          }
        }
      } else { // with groups
        if(g.size() != l) stop("length(g) must match nrow(X)");
        if(narm) {
          NumericVector meanj = no_init_vector(ng), nj = no_init_vector(ng); // stable !!
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            NumericVector M2j(ng, NA_REAL);
            // double meanj[ng], nj[ng]; // stable and faster ??
            // std::vector<double> M2j(ng, NA_REAL);
            // long double d1j = 0;
            double d1j = 0;
            for(int i = l; i--; ) {
              if(std::isnan(column[i])) continue;
              if(std::isnan(M2j[g[i]-1])) {
                meanj[g[i]-1] = column[i];
                M2j[g[i]-1] = 0;
                nj[g[i]-1] = 1;
              } else {
                d1j = column[i]-meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
                M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
              }
            }
            for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = 1/sqrt(M2j[i]/(nj[i]-1));
            for(int i = 0; i != l; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            std::vector<double> meanj(ng), M2j(ng), nj(ng); // faster using std::vector ??
            // long double d1j = 0;
            double d1j = 0;
            int ngs = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(M2j[g[i]-1])) continue;
              if(std::isnan(column[i])) {
                M2j[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) {
                  std::fill(outj.begin(), outj.end(), NA_REAL);
                  goto loopend;
                }
              } else {
                d1j = column[i]-meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
                M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
              }
            }
            for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = 1/sqrt(M2j[i]/(nj[i]-1));
            for(int i = 0; i != l; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
            loopend:;
          }
        }
      }
    } else { // With weights
      NumericVector wg = w;
      if(l != wg.size()) stop("length(w) must match nrow(X)");
      if(ng == 0) {
        if(narm) {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            int k = l-1;
            // long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
            double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
            while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
            if(k != 0) {
              for(int i = k+1; i--; ) {
                if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
                sumwi += wg[i];
                d1i = column[i] - meani;
                meani += d1i * (wg[i] / sumwi);
                M2i += wg[i] * d1i * (column[i] - meani);
              }
              M2i = 1/sqrt(M2i/(sumwi-1));
              if(std::isnan(M2i)) {
                std::fill(outj.begin(), outj.end(), NA_REAL);
              } else {
                outj = (column-meani)*M2i;
                // for(int i = 0; i != l; ++i) outj[i] = (column[i]-meani)*M2i;
              }
            } else {
              std::fill(outj.begin(), outj.end(), NA_REAL);
            }
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            // long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
            double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                M2i = NA_REAL;
                break;
              } else {
                sumwi += wg[i];
                d1i = column[i] - meani;
                meani += d1i * (wg[i] / sumwi);
                M2i += wg[i] * d1i * (column[i] - meani);
              }
            }
            M2i = 1/sqrt(M2i/(sumwi-1));
            if(std::isnan(M2i)) {
              std::fill(outj.begin(), outj.end(), NA_REAL);
            } else {
              outj = (column-meani)*M2i;
              // for(int i = 0; i != l; ++i) outj[i] = (column[i]-meani)*M2i;
            }
          }
        }
      } else { // with groups and weights
        if(g.size() != l) stop("length(g) must match nrow(X)");
        if(narm) {
          std::vector<double> meanj(ng), sumwj(ng);  // NumericVector meanj = no_init_vector(ng), sumwj = no_init_vector(ng); // stable !!
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            NumericVector M2j(ng, NA_REAL);
            // std::vector<double> M2j(ng, NA_REAL); // faster ??
            // double meanj[ng], sumwj[ng]; // stable and faster ??
            // long double d1j = 0;
            double d1j = 0;
            for(int i = l; i--; ) {
              if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
              if(std::isnan(M2j[g[i]-1])) {
                sumwj[g[i]-1] = wg[i];
                meanj[g[i]-1] = column[i];
                M2j[g[i]-1] = 0;
              } else {
                sumwj[g[i]-1] += wg[i];
                d1j = column[i] - meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (wg[i] / sumwj[g[i]-1]);
                M2j[g[i]-1] += wg[i] * d1j * (column[i] - meanj[g[i]-1]);
              }
            }
            for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = 1/sqrt(M2j[i]/(sumwj[i]-1));
            for(int i = 0; i != l; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            std::vector<double> M2j(ng), meanj(ng), sumwj(ng); // faster than NumericVector ??
            // long double d1j = 0;
            double d1j = 0;
            int ngs = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(M2j[g[i]-1])) continue;
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                M2j[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) {
                  std::fill(outj.begin(), outj.end(), NA_REAL);
                  goto loopend2;
                }
              } else {
                sumwj[g[i]-1] += wg[i];
                d1j = column[i] - meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (wg[i] / sumwj[g[i]-1]);
                M2j[g[i]-1] += wg[i] * d1j * (column[i] - meanj[g[i]-1]);
              }
            }
            for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = 1/sqrt(M2j[i]/(sumwj[i]-1));
            for(int i = 0; i != l; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
            loopend2:;
          }
        }
      }
    }

  } else { // ONE-PASS METHOD -------------------------------------

    if (Rf_isNull(w)) { // No weights
      if(ng == 0) {
        if(narm) {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            int k = l-1, nj = 1;
            long double sumj = column[k], sq_sumj = 0;
            while(std::isnan(sumj) && k!=0) sumj = column[--k];
            sq_sumj = sumj*sumj;
            if(k != 0) {
              for(int i = k; i--; ) {
                if(std::isnan(column[i])) continue;
                sumj += column[i];
                sq_sumj += pow(column[i],2);
                ++nj;
              }
              sumj /= nj;
              sq_sumj = 1/sqrt((sq_sumj-pow(sumj,2)*nj)/(nj-1));
              if(std::isnan(sq_sumj)) {
                std::fill(outj.begin(), outj.end(), NA_REAL);
              } else {
                outj = (column-sumj)*sq_sumj;
                // for(int i = 0; i != l; ++i) outj[i] = (column[i]-sumj)*sq_sumj;
              }
            } else {
              std::fill(outj.begin(), outj.end(), NA_REAL);
            }
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            long double sumj = 0, sq_sumj = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(column[i])) {
                sq_sumj = NA_REAL;
                break;
              } else {
                sumj += column[i];
                sq_sumj += pow(column[i],2);
              }
            }
            if(std::isnan(sq_sumj)) {
              std::fill(outj.begin(), outj.end(), NA_REAL);
            } else {
              sumj /= l;
              sq_sumj = 1/sqrt((sq_sumj-pow(sumj,2)*l)/(l-1));
              if(std::isnan(sq_sumj)) {
                std::fill(outj.begin(), outj.end(), NA_REAL);
              } else {
                outj = (column-sumj)*sq_sumj;
                // for(int i = 0; i != l; ++i) outj[i] = (column[i]-sumj)*sq_sumj;
              }
            }
          }
        }
      } else { // with groups
        if(g.size() != l) stop("length(g) must match nrow(X)");
        if(narm) {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            NumericVector sq_sumj(ng, NA_REAL), sumj = no_init_vector(ng);
            IntegerVector nj = no_init_vector(ng);
            // std::vector<double> sq_sumj(ng, NA_REAL); // faster ??
            // double sumj[ng]; // stable and faster ??
            // int nj[ng]; // stable and faster ??
            for(int i = l; i--; ) {
              if(std::isnan(column[i])) continue;
              if(std::isnan(sq_sumj[g[i]-1])) {
                sumj[g[i]-1] = column[i];
                sq_sumj[g[i]-1] = pow(column[i],2);
                nj[g[i]-1] = 1;
              } else {
                sumj[g[i]-1] += column[i];
                sq_sumj[g[i]-1] += pow(column[i],2);
                ++nj[g[i]-1];
              }
            }
            for(int i = ng; i--; ) {
              sumj[i] /= nj[i];
              if(std::isnan(sq_sumj[i])) continue;
              sq_sumj[i] = 1/sqrt((sq_sumj[i] - pow(sumj[i],2)*nj[i])/(nj[i]-1));
            }
            for(int i = 0; i != l; ++i) outj[i] = (column[i]-sumj[g[i]-1])*sq_sumj[g[i]-1];
          }
        } else {
          if(Rf_isNull(gs)) {
            // int gsv[ng], memsize = sizeof(int)*ng;
            for(int j = col; j--; ) {
              NumericMatrix::Column column = x( _ , j);
              NumericMatrix::Column outj = out( _ , j);
              std::vector<double> sq_sumj(ng), sumj(ng); // faster than NumericVector ??
              std::vector<int> gsv(ng); // memset(gsv, 0, memsize);
              int ngs = 0;
              for(int i = 0; i != l; ++i) {
                if(std::isnan(column[i])) {
                  if(std::isnan(sq_sumj[g[i]-1])) continue;
                  sq_sumj[g[i]-1] = column[i];
                  ++ngs;
                  if(ngs == ng) {
                    std::fill(outj.begin(), outj.end(), NA_REAL);
                    goto loopend3;
                  }
                } else {
                  sumj[g[i]-1] += column[i];
                  sq_sumj[g[i]-1] += pow(column[i],2);
                  ++gsv[g[i]-1];
                }
              }
              for(int i = ng; i--; ) {
                sumj[i] /= gsv[i];
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = 1/sqrt((sq_sumj[i] - pow(sumj[i],2)*gsv[i])/(gsv[i]-1));
              }
              for(int i = 0; i != l; ++i) outj[i] = (column[i]-sumj[g[i]-1])*sq_sumj[g[i]-1];
              loopend3:;
            }
          } else {
            IntegerVector gsv = gs;
            if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
            for(int j = col; j--; ) {
              NumericMatrix::Column column = x( _ , j);
              NumericMatrix::Column outj = out( _ , j);
              std::vector<double> sq_sumj(ng), sumj(ng); // faster than NumericVector ??
              int ngs = 0;
              for(int i = 0; i != l; ++i) {
                if(std::isnan(column[i])) {
                  if(std::isnan(sq_sumj[g[i]-1])) continue;
                  sq_sumj[g[i]-1] = column[i];
                  ++ngs;
                  if(ngs == ng) {
                    std::fill(outj.begin(), outj.end(), NA_REAL);
                    goto loopend5;
                  }
                } else {
                  sumj[g[i]-1] += column[i];
                  sq_sumj[g[i]-1] += pow(column[i],2);
                }
              }
              for(int i = ng; i--; ) {
                sumj[i] /= gsv[i];
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = 1/sqrt((sq_sumj[i] - pow(sumj[i],2)*gsv[i])/(gsv[i]-1));
              }
              for(int i = 0; i != l; ++i) outj[i] = (column[i]-sumj[g[i]-1])*sq_sumj[g[i]-1];
              loopend5:;
            }
          }
        }
      }
    } else { // With weights
      NumericVector wg = w;
      if(l != wg.size()) stop("length(w) must match nrow(X)");
      if(ng == 0) {
        if(narm) {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            int k = l-1;
            while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
            long double sumwj = wg[k], sumj = column[k]*sumwj, sq_sumj = column[k]*sumj;
            if(k != 0) {
              for(int i = k; i--; ) {
                if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
                sumj += column[i]*wg[i];
                sumwj += wg[i];
                sq_sumj += pow(column[i],2)*wg[i];
              }
              sumj /= sumwj;
              sq_sumj = 1/sqrt((sq_sumj - pow(sumj,2)*sumwj)/(sumwj-1));
              if(std::isnan(sq_sumj)) {
                std::fill(outj.begin(), outj.end(), NA_REAL);
              } else {
                outj = (column-sumj)*sq_sumj;
                // for(int i = 0; i != l; ++i) outj[i] = (column[i]-sumj)*sq_sumj;
              }
            } else {
              std::fill(outj.begin(), outj.end(), NA_REAL);
            }
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            long double sumj = 0, sumwj = 0, sq_sumj = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                sq_sumj = NA_REAL;
                break;
              } else {
                sumj += column[i]*wg[i];
                sumwj += wg[i];
                sq_sumj += pow(column[i],2)*wg[i];
              }
            }
            if(std::isnan(sq_sumj)) {
              std::fill(outj.begin(), outj.end(), NA_REAL);
            } else {
              sumj /= sumwj;
              sq_sumj = 1/sqrt((sq_sumj - pow(sumj,2)*sumwj)/(sumwj-1));
              if(std::isnan(sq_sumj)) {
                std::fill(outj.begin(), outj.end(), NA_REAL);
              } else {
                outj = (column-sumj)*sq_sumj;
                // for(int i = 0; i != l; ++i) outj[i] = (column[i]-sumj)*sq_sumj;
              }
            }
          }
        }
      } else { // with groups and weights
        if(g.size() != l) stop("length(g) must match nrow(X)");
        if(narm) {
          std::vector<double> sumj(ng), sumwj(ng); // stable !!
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            NumericVector sq_sumj(ng, NA_REAL);
            // std::vector<double> sq_sumj(ng, NA_REAL); // faster ??
            // double sumj[ng], sumwj[ng]; // stable and faster ??
            for(int i = l; i--; ) {
              if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
              if(std::isnan(sq_sumj[g[i]-1])) {
                sumj[g[i]-1] = column[i]*wg[i];
                sumwj[g[i]-1] = wg[i];
                sq_sumj[g[i]-1] = pow(column[i],2)*wg[i];
              } else {
                sumj[g[i]-1] += column[i]*wg[i];
                sumwj[g[i]-1] += wg[i];
                sq_sumj[g[i]-1] += pow(column[i],2)*wg[i];
              }
            }
            for(int i = ng; i--; ) {
              sumj[i] /= sumwj[i];
              if(std::isnan(sq_sumj[i])) continue;
              sq_sumj[i] = 1/sqrt((sq_sumj[i] - pow(sumj[i],2)*sumwj[i])/(sumwj[i]-1));
            }
            for(int i = 0; i != l; ++i) outj[i] = (column[i]-sumj[g[i]-1])*sq_sumj[g[i]-1];
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::Column column = x( _ , j);
            NumericMatrix::Column outj = out( _ , j);
            std::vector<double> sq_sumj(ng), sumj(ng), sumwj(ng); // faster than NumericVector ??
            int ngs = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(sq_sumj[g[i]-1])) continue;
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                sq_sumj[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) {
                  std::fill(outj.begin(), outj.end(), NA_REAL);
                  goto loopend4;
                }
              } else {
                sumj[g[i]-1] += column[i]*wg[i];
                sumwj[g[i]-1] += wg[i];
                sq_sumj[g[i]-1] += pow(column[i],2)*wg[i];
              }
            }
            for(int i = ng; i--; ) {
              sumj[i] /= sumwj[i]; //
              if(std::isnan(sq_sumj[i])) continue;
              sq_sumj[i] = 1/sqrt((sq_sumj[i] - pow(sumj[i],2)*sumwj[i])/(sumwj[i]-1));
            }
            for(int i = 0; i != l; ++i) outj[i] = (column[i]-sumj[g[i]-1])*sq_sumj[g[i]-1];
            loopend4:;
          }
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  return out;
}




// [[Rcpp::export]]
List fscalelCpp(List x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue,
                SEXP w = R_NilValue, bool narm = true, bool stable_algo = true) {

  int l = x.size();
  List out(l);

  if(stable_algo) { // WELFORDS ONLINE METHOD -------------------------------------
    if (Rf_isNull(w)) { // No weights
      if(ng == 0) {
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j]; // outj = NULL gives compile warning
            int row = column.size(), k = row-1;
            NumericVector outj = no_init_vector(row);
            // double ni = 0;
            // long double meani = 0, d1i = 0, M2i = 0;
            double ni = 0, meani = 0, d1i = 0, M2i = 0;
            while(std::isnan(column[k]) && k!=0) --k;
            if(k != 0) {
              for(int i = k+1; i--; ) {
                if(std::isnan(column[i])) continue;
                d1i = column[i]-meani;
                meani += d1i * (1 / ++ni);
                M2i += d1i*(column[i]-meani);
              }
              M2i = 1/sqrt(M2i/(ni-1));
              if(std::isnan(M2i)) {
                outj = rep(NA_REAL, row); // fastest option !!
              } else {
                outj = (column-meani)*M2i;
              }
            } else {
              outj = rep(NA_REAL, row);
            }
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j]; //  outj = NULL; gives compile warning
            int row = column.size();
            NumericVector outj = no_init_vector(row);
            // double ni = 0;
            // long double meani = 0, d1i = 0, M2i = 0;
            double ni = 0, meani = 0, d1i = 0, M2i = 0;
            for(int i = 0; i != row; ++i) {
              if(std::isnan(column[i])) {
                M2i = NA_REAL;
                break;
              } else {
                d1i = column[i]-meani;
                meani += d1i * (1 / ++ni);
                M2i += d1i*(column[i]-meani);
              }
            }
            M2i = 1/sqrt(M2i/(row-1));
            if(std::isnan(M2i)) {
              outj = rep(NA_REAL, row);
            } else {
              outj = (column-meani)*M2i;
            }
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        }
      } else { // with groups
        int gss = g.size();
        if(narm) {
          NumericVector meanj = no_init_vector(ng), nj = no_init_vector(ng); // stable !!
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector M2j(ng, NA_REAL);
            // std::vector<double> M2j(ng, NA_REAL); // faster and stable ??
            // double meanj[ng], nj[ng];
            // long double d1j = 0;
            double d1j = 0;
            for(int i = gss; i--; ) {
              if(std::isnan(column[i])) continue;
              if(std::isnan(M2j[g[i]-1])) {
                meanj[g[i]-1] = column[i];
                M2j[g[i]-1] = 0;
                nj[g[i]-1] = 1;
              } else {
                d1j = column[i]-meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
                M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
              }
            }
            for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = 1/sqrt(M2j[i]/(nj[i]-1));
            NumericVector outj = no_init_vector(gss);
            for(int i = 0; i != gss; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector outj = no_init_vector(gss);
            std::vector<double> meanj(ng), M2j(ng), nj(ng); // faster than NumericVector ??
            // long double d1j = 0;
            double d1j = 0;
            int ngs = 0;
            {
              for(int i = 0; i != gss; ++i) {
                if(std::isnan(M2j[g[i]-1])) continue;
                if(std::isnan(column[i])) {
                  M2j[g[i]-1] = NA_REAL;
                  ++ngs;
                  if(ngs == ng) {
                    outj = rep(NA_REAL, gss);
                    goto loopend;
                  }
                } else {
                  d1j = column[i]-meanj[g[i]-1];
                  meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
                  M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
                }
              }
              for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = 1/sqrt(M2j[i]/(nj[i]-1));
              for(int i = 0; i != gss; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
            }
            loopend:;
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        }
      }
    } else { // With weights
      NumericVector wg = w;
      int wgs = wg.size();
      if(ng == 0) {
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j]; // outj = NULL; gives compile warning
            if(wgs != column.size()) stop("length(w) must match nrow(X)");
            NumericVector outj = no_init_vector(wgs);
            int k = wgs-1;
            // long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
            double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
            while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
            if(k != 0) {
              for(int i = k+1; i--; ) {
                if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
                sumwi += wg[i];
                d1i = column[i] - meani;
                meani += d1i * (wg[i] / sumwi);
                M2i += wg[i] * d1i * (column[i] - meani);
              }
              M2i = 1/sqrt(M2i/(sumwi-1));
              if(std::isnan(M2i)) {
                outj = rep(NA_REAL, wgs);
              } else {
                outj = (column-meani)*M2i;
              }
            } else {
              outj = rep(NA_REAL, wgs);
            }
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j]; //  outj = NULL; gives compile warning
            if(wgs != column.size()) stop("length(w) must match nrow(X)");
            NumericVector outj = no_init_vector(wgs);
            // long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
            double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
            for(int i = 0; i != wgs; ++i) {
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                M2i = NA_REAL;
                break;
              } else {
                sumwi += wg[i];
                d1i = column[i] - meani;
                meani += d1i * (wg[i] / sumwi);
                M2i += wg[i] * d1i * (column[i] - meani);
              }
            }
            M2i = 1/sqrt(M2i/(sumwi-1));
            if(std::isnan(M2i)) {
              outj = rep(NA_REAL, wgs);
            } else {
              outj = (column-meani)*M2i;
            }
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        }
      } else { // with groups and weights
        int gss = g.size();
        if(gss != wgs) stop("length(w) must match length(g)");
        if(narm) {
          // NumericVector meanj = no_init_vector(ng), sumwj = no_init_vector(ng);
          std::vector<double> meanj(ng), sumwj(ng); // stable !!
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector M2j(ng, NA_REAL);
            // std::vector<double> M2j(ng, NA_REAL); // faster and stable ??
            // double meanj[ng], sumwj[ng];
            // long double d1j = 0;
            double d1j = 0;
            for(int i = gss; i--; ) {
              if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
              if(std::isnan(M2j[g[i]-1])) {
                sumwj[g[i]-1] = wg[i];
                meanj[g[i]-1] = column[i];
                M2j[g[i]-1] = 0;
              } else {
                sumwj[g[i]-1] += wg[i];
                d1j = column[i] - meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (wg[i] / sumwj[g[i]-1]);
                M2j[g[i]-1] += wg[i] * d1j * (column[i] - meanj[g[i]-1]);
              }
            }
            for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = 1/sqrt(M2j[i]/(sumwj[i]-1));
            NumericVector outj = no_init_vector(gss);
            for(int i = 0; i != gss; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector outj = no_init_vector(gss);
            std::vector<double> M2j(ng), meanj(ng), sumwj(ng); // faster than NumericVector??
            // long double d1j = 0;
            double d1j = 0;
            int ngs = 0;
            {
              for(int i = 0; i != gss; ++i) {
                if(std::isnan(M2j[g[i]-1])) continue;
                if(std::isnan(column[i]) || std::isnan(wg[i])) {
                  M2j[g[i]-1] = NA_REAL;
                  ++ngs;
                  if(ngs == ng) {
                    outj = rep(NA_REAL, gss);
                    goto loopend2;
                  }
                } else {
                  sumwj[g[i]-1] += wg[i];
                  d1j = column[i] - meanj[g[i]-1];
                  meanj[g[i]-1] += d1j * (wg[i] / sumwj[g[i]-1]);
                  M2j[g[i]-1] += wg[i] * d1j * (column[i] - meanj[g[i]-1]);
                }
              }
              for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = 1/sqrt(M2j[i]/(sumwj[i]-1));
              for(int i = 0; i != gss; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
            }
            loopend2:;
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        }
      }
    }

  } else { // ONE-PASS METHOD -------------------------------------

    if (Rf_isNull(w)) { // No weights
      if(ng == 0) {
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j]; // outj = NULL; gives compile error
            int row = column.size(), k = row-1, nj = 1;
            NumericVector outj = no_init_vector(row);
            long double sumj = column[k], sq_sumj = 0;
            while(std::isnan(sumj) && k!=0) sumj = column[--k];
            sq_sumj = sumj*sumj;
            if(k != 0) {
              for(int i = k; i--; ) {
                if(std::isnan(column[i])) continue;
                sumj += column[i];
                sq_sumj += pow(column[i],2);
                ++nj;
              }
              sumj /= nj;
              sq_sumj = 1/sqrt((sq_sumj-pow(sumj,2)*nj)/(nj-1));
              if(std::isnan(sq_sumj)) {
                outj = rep(NA_REAL, row);
              } else {
                outj = (column-sumj)*sq_sumj;
              }
            } else {
              outj = rep(NA_REAL, row);
            }
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j]; // outj = NULL; gives compile warning
            int row = column.size();
            NumericVector outj = no_init_vector(row);
            long double sumj = 0, sq_sumj = 0;
            for(int i = 0; i != row; ++i) {
              if(std::isnan(column[i])) {
                sq_sumj = NA_REAL;
                break;
              } else {
                sumj += column[i];
                sq_sumj += pow(column[i],2);
              }
            }
            if(std::isnan(sq_sumj)) {
              outj = rep(NA_REAL, row);
            } else {
              sumj /= row;
              sq_sumj = 1/sqrt((sq_sumj-pow(sumj,2)*row)/(row-1));
              if(std::isnan(sq_sumj)) {
                outj = rep(NA_REAL, row);
              } else {
                outj = (column-sumj)*sq_sumj;
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        }
      } else { // with groups
        int gss = g.size();
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector sq_sumj(ng, NA_REAL), sumj(ng); //  = no_init_vector
            IntegerVector nj(ng); //  = no_init_vector
            // std::vector<double> sq_sumj(ng, NA_REAL); // faster and stable ??
            // double sumj[ng];
            // int nj[ng];
            for(int i = gss; i--; ) {
              if(std::isnan(column[i])) continue;
              if(std::isnan(sq_sumj[g[i]-1])) {
                sumj[g[i]-1] = column[i];
                sq_sumj[g[i]-1] = pow(column[i],2);
                nj[g[i]-1] = 1;
              } else {
                sumj[g[i]-1] += column[i];
                sq_sumj[g[i]-1] += pow(column[i],2);
                ++nj[g[i]-1];
              }
            }
            for(int i = ng; i--; ) {
              sumj[i] /= nj[i];
              if(std::isnan(sq_sumj[i])) continue;
              sq_sumj[i] = 1/sqrt((sq_sumj[i] - pow(sumj[i],2)*nj[i])/(nj[i]-1));
            }
            NumericVector outj = no_init_vector(gss);
            for(int i = 0; i != gss; ++i) outj[i] = (column[i]-sumj[g[i]-1])*sq_sumj[g[i]-1];
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        } else {
          if(Rf_isNull(gs)) {
            // int gsv[ng], memsize = sizeof(int)*ng;
            for(int j = l; j--; ) {
              NumericVector column = x[j];
              if(gss != column.size()) stop("length(g) must match nrow(X)");
              int ngs = 0;
              NumericVector outj = no_init_vector(gss);
              std::vector<double> sq_sumj(ng), sumj(ng); // faster than NumericVector??
              std::vector<int> gsv(ng); // memset(gsv, 0, memsize);
              {
                for(int i = 0; i != gss; ++i) {
                  if(std::isnan(column[i])) {
                    if(std::isnan(sq_sumj[g[i]-1])) continue;
                    sq_sumj[g[i]-1] = column[i];
                    ++ngs;
                    if(ngs == ng) {
                      outj = rep(NA_REAL, gss);
                      goto loopend3;
                    }
                  } else {
                    sumj[g[i]-1] += column[i];
                    sq_sumj[g[i]-1] += pow(column[i],2);
                    ++gsv[g[i]-1];
                  }
                }
                for(int i = ng; i--; ) {
                  sumj[i] /= gsv[i];
                  if(std::isnan(sq_sumj[i])) continue;
                  sq_sumj[i] = 1/sqrt((sq_sumj[i] - pow(sumj[i],2)*gsv[i])/(gsv[i]-1));
                }
                for(int i = 0; i != gss; ++i) outj[i] = (column[i]-sumj[g[i]-1])*sq_sumj[g[i]-1];
              }
              loopend3:;
              SHALLOW_DUPLICATE_ATTRIB(outj, column);
              out[j] = outj;
            }
          } else {
            IntegerVector gsv = gs;
            if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
            for(int j = l; j--; ) {
              NumericVector column = x[j];
              if(gss != column.size()) stop("length(g) must match nrow(X)");
              int ngs = 0;
              NumericVector outj = no_init_vector(gss);
              std::vector<double> sq_sumj(ng), sumj(ng); // faster than NumericVector??
              {
                for(int i = 0; i != gss; ++i) {
                  if(std::isnan(column[i])) {
                    if(std::isnan(sq_sumj[g[i]-1])) continue;
                    sq_sumj[g[i]-1] = column[i];
                    ++ngs;
                    if(ngs == ng) {
                      outj = rep(NA_REAL, gss);
                      goto loopend5;
                    }
                  } else {
                    sumj[g[i]-1] += column[i];
                    sq_sumj[g[i]-1] += pow(column[i],2);
                  }
                }
                for(int i = ng; i--; ) {
                  sumj[i] /= gsv[i];
                  if(std::isnan(sq_sumj[i])) continue;
                  sq_sumj[i] = 1/sqrt((sq_sumj[i] - pow(sumj[i],2)*gsv[i])/(gsv[i]-1));
                }
                for(int i = 0; i != gss; ++i) outj[i] = (column[i]-sumj[g[i]-1])*sq_sumj[g[i]-1];
              }
              loopend5:;
              SHALLOW_DUPLICATE_ATTRIB(outj, column);
              out[j] = outj;
            }
          }
        }
      }
    } else { // With weights
      NumericVector wg = w;
      int wgs = wg.size();
      if(ng == 0) {
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j]; // outj = NULL; gives compile warning
            if(wgs != column.size()) stop("length(w) must match nrow(X)");
            NumericVector outj = no_init_vector(wgs);
            int k = wgs-1;
            while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
            long double sumwj = wg[k], sumj = column[k]*sumwj, sq_sumj = column[k]*sumj;
            if(k != 0) {
              for(int i = k; i--; ) {
                if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
                sumj += column[i]*wg[i];
                sumwj += wg[i];
                sq_sumj += pow(column[i],2)*wg[i];
              }
              sumj /= sumwj;
              sq_sumj = 1/sqrt((sq_sumj - pow(sumj,2)*sumwj)/(sumwj-1));
              if(std::isnan(sq_sumj)) {
                outj = rep(NA_REAL, wgs);
              } else {
                outj = (column-sumj)*sq_sumj;
              }
            } else {
              outj = rep(NA_REAL, wgs);
            }
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j]; // outj = NULL; gives compile warning
            if(wgs != column.size()) stop("length(w) must match nrow(X)");
            NumericVector outj = no_init_vector(wgs);
            long double sumj = 0, sumwj = 0, sq_sumj = 0;
            for(int i = 0; i != wgs; ++i) {
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                sq_sumj = NA_REAL;
                break;
              } else {
                sumj += column[i]*wg[i];
                sumwj += wg[i];
                sq_sumj += pow(column[i],2)*wg[i];
              }
            }
            if(std::isnan(sq_sumj)) {
              outj = rep(NA_REAL, wgs);
            } else {
              sumj /= sumwj;
              sq_sumj = 1/sqrt((sq_sumj - pow(sumj,2)*sumwj)/(sumwj-1));
              if(std::isnan(sq_sumj)) {
                outj = rep(NA_REAL, wgs);
              } else {
                outj = (column-sumj)*sq_sumj;
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        }
      } else { // with groups and weights
        int gss = g.size();
        if(gss != wgs) stop("length(w) must match length(g)");
        if(narm) {
          // NumericVector sumj = no_init_vector(ng), sumwj = no_init_vector(ng);
          std::vector<double> sumj(ng), sumwj(ng); // stable !!
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector sq_sumj(ng, NA_REAL);
            // std::vector<double> sq_sumj(ng, NA_REAL); // faster and stable ??
            // double sumj[ng], sumwj[ng];
            for(int i = gss; i--; ) {
              if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
              if(std::isnan(sq_sumj[g[i]-1])) {
                sumj[g[i]-1] = column[i]*wg[i];
                sumwj[g[i]-1] = wg[i];
                sq_sumj[g[i]-1] = pow(column[i],2)*wg[i];
              } else {
                sumj[g[i]-1] += column[i]*wg[i];
                sumwj[g[i]-1] += wg[i];
                sq_sumj[g[i]-1] += pow(column[i],2)*wg[i];
              }
            }
            for(int i = ng; i--; ) {
              sumj[i] /= sumwj[i];
              if(std::isnan(sq_sumj[i])) continue;
              sq_sumj[i] = 1/sqrt((sq_sumj[i] - pow(sumj[i],2)*sumwj[i])/(sumwj[i]-1));
            }
            NumericVector outj = no_init_vector(gss);
            for(int i = 0; i != gss; ++i) outj[i] = (column[i]-sumj[g[i]-1])*sq_sumj[g[i]-1];
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector outj = no_init_vector(gss);
            std::vector<double> sq_sumj(ng), sumj(ng), sumwj(ng); // faster than NumericVector ??
            int ngs = 0;
            {
              for(int i = 0; i != gss; ++i) {
                if(std::isnan(sq_sumj[g[i]-1])) continue;
                if(std::isnan(column[i]) || std::isnan(wg[i])) {
                  sq_sumj[g[i]-1] = NA_REAL;
                  ++ngs;
                  if(ngs == ng) {
                    outj = rep(NA_REAL, gss);
                    goto loopend4;
                  }
                } else {
                  sumj[g[i]-1] += column[i]*wg[i];
                  sumwj[g[i]-1] += wg[i];
                  sq_sumj[g[i]-1] += pow(column[i],2)*wg[i];
                }
              }
              for(int i = ng; i--; ) {
                sumj[i] /= sumwj[i]; //
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = 1/sqrt((sq_sumj[i] - pow(sumj[i],2)*sumwj[i])/(sumwj[i]-1));
              }
              for(int i = 0; i != gss; ++i) outj[i] = (column[i]-sumj[g[i]-1])*sq_sumj[g[i]-1];
            }
            loopend4:;
            SHALLOW_DUPLICATE_ATTRIB(outj, column);
            out[j] = outj;
          }
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  return out;
}
