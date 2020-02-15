#include <Rcpp.h>
using namespace Rcpp;

// Note: All comments are in fvar.cpp

// [[Rcpp::export]]
NumericVector fvarsdCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0, const SEXP& gs = R_NilValue,
                        const SEXP& w = R_NilValue, bool narm = true, bool stable_algo = true, bool sd = true) {
  int l = x.size();
  if(l == 1) return NumericVector::create(NA_REAL);

  if(stable_algo) { // WELFORDS ONLINE METHOD ---------------------------------------------------------
    if(Rf_isNull(w)) { // No weights
      if(ng == 0) {
        if(narm) {
          int j = l-1;
          // double n = 0;
          // long double mean = 0, d1 = 0, M2 = 0; // LD really necessary ?? what about speed ??
          double n = 0, mean = 0, d1 = 0, M2 = 0;
          while(std::isnan(x[j]) && j!=0) --j;
          if(j != 0) {
            for(int i = j+1; i--; ) {
              if(std::isnan(x[i])) continue;
              d1 = x[i]-mean;
              mean += d1 * (1 / ++n);
              M2 += d1*(x[i]-mean);
            }
            M2 = M2/(n-1);
            if(sd) M2 = sqrt(M2);
            if(std::isnan(M2)) M2 = NA_REAL;
            return NumericVector::create(M2); // ::create((double)M2)
          } else return NumericVector::create(NA_REAL);
        } else {
          // double n = 0;
          // long double mean = 0, d1 = 0, M2 = 0;
          double n = 0, mean = 0, d1 = 0, M2 = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(x[i])) {
              return NumericVector::create(NA_REAL);
            } else {
              d1 = x[i]-mean;
              mean += d1*(1 / ++n);
              M2 += d1*(x[i]-mean);
            }
          }
          M2 = M2/(l-1);
          if(sd) M2 = sqrt(M2);
          if(std::isnan(M2)) M2 = NA_REAL;
          return NumericVector::create(M2); // ::create((double)M2)
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
          if(sd) {
            for(int i = ng; i--; ) {
              if(std::isnan(M2[i])) M2[i] = NA_REAL;
              else {
                M2[i] = sqrt(M2[i]/(n[i]-1));
                if(std::isnan(M2[i])) M2[i] = NA_REAL;
              }
            }
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(M2[i])) M2[i] = NA_REAL;
              else {
                M2[i] /= n[i]-1;
                if(std::isnan(M2[i])) M2[i] = NA_REAL;
              }
            }
          }
          DUPLICATE_ATTRIB(M2, x);
          return M2;
        } else {
          NumericVector M2(ng), mean(ng), n(ng);
          int ngs = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(M2[g[i]-1])) continue;
            if(std::isnan(x[i])) {
              M2[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) break;
            } else {
              d1 = x[i]-mean[g[i]-1];
              mean[g[i]-1] += d1 * (1 / ++n[g[i]-1]);
              M2[g[i]-1] += d1*(x[i]-mean[g[i]-1]);
            }
          }
          if(sd) {
            for(int i = ng; i--; ) {
              if(std::isnan(M2[i])) M2[i] = NA_REAL;
              else {
                M2[i] = sqrt(M2[i]/(n[i]-1));
                if(std::isnan(M2[i])) M2[i] = NA_REAL;
              }
            }
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(M2[i])) M2[i] = NA_REAL;
              else {
                M2[i] /= n[i]-1;
                if(std::isnan(M2[i])) M2[i] = NA_REAL;
              }
            }
          }
          DUPLICATE_ATTRIB(M2, x);
          return M2;
        }
      }
    } else { // With weights
      NumericVector wg = w;
      if(l != wg.size()) stop("length(w) must match length(x)");
      if(ng == 0) {
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
            M2 /= sumw-1;
            if(sd) M2 = sqrt(M2);
            if(std::isnan(M2)) M2 = NA_REAL;
            return NumericVector::create(M2); // create((double)M2)
          } else return NumericVector::create(NA_REAL);
        } else {
          // long double sumw = 0, mean = 0, M2 = 0, d1 = 0;
          double sumw = 0, mean = 0, M2 = 0, d1 = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(x[i]) || std::isnan(wg[i])) {
              return NumericVector::create(NA_REAL);
            } else {
              sumw += wg[i];
              d1 = x[i] - mean;
              mean += d1 * (wg[i] / sumw);
              M2 += wg[i] * d1 * (x[i] - mean);
            }
          }
          M2 /= sumw-1;
          if(sd) M2 = sqrt(M2);
          if(std::isnan(M2)) M2 = NA_REAL;
          return NumericVector::create(M2); // create((double)M2)
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
          if(sd) {
            for(int i = ng; i--; ) {
              if(std::isnan(M2[i])) M2[i] = NA_REAL;
              else {
                M2[i] = sqrt(M2[i]/(sumw[i]-1));
                if(std::isnan(M2[i])) M2[i] = NA_REAL;
              }
            }
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(M2[i])) M2[i] = NA_REAL;
              else {
                M2[i] /= sumw[i]-1;
                if(std::isnan(M2[i])) M2[i] = NA_REAL;
              }
            }
          }
          DUPLICATE_ATTRIB(M2, x);
          return M2;
        } else {
          NumericVector M2(ng), sumw(ng), mean(ng);
          int ngs = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(M2[g[i]-1])) continue;
            if(std::isnan(x[i]) || std::isnan(wg[i])) {
              M2[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) break;
            } else {
              sumw[g[i]-1] += wg[i];
              d1 = x[i] - mean[g[i]-1];
              mean[g[i]-1] += d1 * (wg[i] / sumw[g[i]-1]);
              M2[g[i]-1] += wg[i] * d1 * (x[i] - mean[g[i]-1]);
            }
          }
          if(sd) {
            for(int i = ng; i--; ) {
              if(std::isnan(M2[i])) M2[i] = NA_REAL;
              else {
                M2[i] = sqrt(M2[i]/(sumw[i]-1));
                if(std::isnan(M2[i])) M2[i] = NA_REAL;
              }
            }
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(M2[i])) M2[i] = NA_REAL;
              else {
                M2[i] /= sumw[i]-1;
                if(std::isnan(M2[i])) M2[i] = NA_REAL;
              }
            }
          }
          DUPLICATE_ATTRIB(M2, x);
          return M2;
        }
      }
    }

  } else { // ONE-PASS METHOD ---------------------------------------------------------
    if(Rf_isNull(w)) { // No weights
      if(ng == 0) {
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
            sq_sum = (sq_sum - pow(sum/n,2)*n)/(n-1);
            if(sd) sq_sum = sqrt(sq_sum);
            if(std::isnan(sq_sum)) sq_sum = NA_REAL;
            return NumericVector::create((double)sq_sum);
          } else return NumericVector::create(NA_REAL);
        } else {
          long double sum = 0, sq_sum = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(x[i])) {
              return NumericVector::create(x[i]);
            } else {
              sum += x[i];
              sq_sum += pow(x[i],2);
            }
          }
          sq_sum = (sq_sum - pow(sum/l,2)*l)/(l-1);
          if(sd) sq_sum = sqrt(sq_sum);
          if(std::isnan(sq_sum)) sq_sum = NA_REAL;
          return NumericVector::create((double)sq_sum);
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
          if(sd) {
            for(int i = ng; i--; ) {
              if(std::isnan(sq_sum[i])) continue;
              sq_sum[i] = sqrt((sq_sum[i] - pow(sum[i]/n[i],2)*n[i])/(n[i]-1));
              if(std::isnan(sq_sum[i])) sq_sum[i] = NA_REAL;
            }
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(sq_sum[i])) continue;
              sq_sum[i] = (sq_sum[i] - pow(sum[i]/n[i],2)*n[i])/(n[i]-1);
              if(std::isnan(sq_sum[i])) sq_sum[i] = NA_REAL;
            }
          }
          DUPLICATE_ATTRIB(sq_sum, x);
          return sq_sum;
        } else {
          NumericVector sq_sum(ng), sum(ng);
          IntegerVector gsv = no_init_vector(ng); // NULL; gives compile warning
          int ngs = 0;
          if(Rf_isNull(gs)) {
            // gsv = IntegerVector(ng);
            std::fill(gsv.begin(), gsv.end(), 0);
            for(int i = 0; i != l; ++i) {
              if(std::isnan(x[i])) {
                if(std::isnan(sq_sum[g[i]-1])) continue;
                sq_sum[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) break;
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
                if(ngs == ng) break;
              } else {
                sum[g[i]-1] += x[i];
                sq_sum[g[i]-1] += pow(x[i],2);
              }
            }
          }
          if(sd) {
            for(int i = ng; i--; ) {
              if(std::isnan(sq_sum[i])) continue;
              sq_sum[i] = sqrt((sq_sum[i] - pow(sum[i]/gsv[i],2)*gsv[i])/(gsv[i]-1));
              if(std::isnan(sq_sum[i])) sq_sum[i] = NA_REAL;
            }
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(sq_sum[i])) continue;
              sq_sum[i] = (sq_sum[i] - pow(sum[i]/gsv[i],2)*gsv[i])/(gsv[i]-1);
              if(std::isnan(sq_sum[i])) sq_sum[i] = NA_REAL;
            }
          }
          DUPLICATE_ATTRIB(sq_sum, x);
          return sq_sum;
        }
      }
    } else { // With weights
      NumericVector wg = w;
      if(l != wg.size()) stop("length(w) must match length(x)");
      if(ng == 0) {
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
            sq_sum = (sq_sum - pow(sum/sumw,2)*sumw)/(sumw-1);
            if(sd) sq_sum = sqrt(sq_sum);
            if(std::isnan(sq_sum)) sq_sum = NA_REAL;
            return NumericVector::create((double)sq_sum);
          } else return NumericVector::create(NA_REAL);
        } else {
          long double sum = 0, sumw = 0, sq_sum = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(x[i]) || std::isnan(wg[i])) {
              return NumericVector::create(NA_REAL);
            } else {
              sum += x[i]*wg[i];
              sumw += wg[i];
              sq_sum += pow(x[i],2)*wg[i];
            }
          }
          sq_sum = (sq_sum - pow(sum/sumw,2)*sumw)/(sumw-1);
          if(sd) sq_sum = sqrt(sq_sum);
          if(std::isnan(sq_sum)) sq_sum = NA_REAL;
          return NumericVector::create((double)sq_sum);
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
          if(sd) {
            for(int i = ng; i--; ) {
              if(std::isnan(sq_sum[i])) continue;
              sq_sum[i] = sqrt((sq_sum[i] - pow(sum[i]/sumw[i],2)*sumw[i])/(sumw[i]-1));
              if(std::isnan(sq_sum[i])) sq_sum[i] = NA_REAL;
            }
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(sq_sum[i])) continue;
              sq_sum[i] = (sq_sum[i] - pow(sum[i]/sumw[i],2)*sumw[i])/(sumw[i]-1);
              if(std::isnan(sq_sum[i])) sq_sum[i] = NA_REAL;
            }
          }
          DUPLICATE_ATTRIB(sq_sum, x);
          return sq_sum;
        } else {
          NumericVector sq_sum(ng), sumw(ng), sum(ng);
          int ngs = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(sq_sum[g[i]-1])) continue;
            if(std::isnan(x[i]) || std::isnan(wg[i])) {
              sq_sum[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) break;
            } else {
              sum[g[i]-1] += x[i]*wg[i];
              sumw[g[i]-1] += wg[i];
              sq_sum[g[i]-1] += pow(x[i],2)*wg[i];
            }
          }
          if(sd) {
            for(int i = ng; i--; ) {
              if(std::isnan(sq_sum[i])) continue;
              sq_sum[i] = sqrt((sq_sum[i] - pow(sum[i]/sumw[i],2)*sumw[i])/(sumw[i]-1));
              if(std::isnan(sq_sum[i])) sq_sum[i] = NA_REAL;
            }
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(sq_sum[i])) continue;
              sq_sum[i] = (sq_sum[i] - pow(sum[i]/sumw[i],2)*sumw[i])/(sumw[i]-1);
              if(std::isnan(sq_sum[i])) sq_sum[i] = NA_REAL;
            }
          }
          DUPLICATE_ATTRIB(sq_sum, x);
          return sq_sum;
        }
      }
    }
  }
}





// [[Rcpp::export]]
SEXP fvarsdmCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0,
                const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue,
                bool narm = true, bool stable_algo = true,
                bool sd = true, bool drop = true) {
  int l = x.nrow(), col = x.ncol();

  if(stable_algo) { // WELFORDS ONLINE METHOD -------------------------------------
    if(Rf_isNull(w)) { // No weights
      if(ng == 0) {
        NumericVector out = no_init_vector(col);
        if(narm) {
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
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
              M2i /= ni-1;
              if(sd) M2i = sqrt(M2i);
              if(std::isnan(M2i)) M2i = NA_REAL;
              out[j] = M2i; // (double)M2i;
            } else out[j] = NA_REAL;
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
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
            M2i /= l-1;
            if(sd) M2i = sqrt(M2i);
            if(std::isnan(M2i)) M2i = NA_REAL;
            out[j] = M2i; // (double)M2i;
          }
        }
        if(drop) out.attr("names") = colnames(x);
        else {
          out.attr("dim") = Dimension(1, col);
          colnames(out) = colnames(x);
        }
        return out;
      } else { // with groups
        if(g.size() != l) stop("length(g) must match nrow(X)");
        if(narm) {
          NumericMatrix M2 = no_init_matrix(ng, col);
          std::fill(M2.begin(), M2.end(), NA_REAL);
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            NumericMatrix::Column M2j = M2( _ , j);
            double d1j = 0; // , meanj[ng], nj[ng];
            NumericVector meanj = no_init_vector(ng), nj = no_init_vector(ng);
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
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] = sqrt(M2j[i]/(nj[i]-1));
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] /= nj[i]-1;
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            }
          }
          colnames(M2) = colnames(x);
          return M2;
        } else {
          NumericMatrix M2(ng, col);
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            NumericMatrix::Column M2j = M2( _ , j);
            std::vector<double> meanj(ng), nj(ng);
            double d1j = 0;
            int ngs = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(M2j[g[i]-1])) continue;
              if(std::isnan(column[i])) {
                M2j[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) break;
              } else {
                d1j = column[i]-meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
                M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
              }
            }
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] = sqrt(M2j[i]/(nj[i]-1));
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] /= nj[i]-1;
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            }
          }
          colnames(M2) = colnames(x);
          return M2;
        }
      }
    } else { // With weights
      NumericVector wg = w;
      if(l != wg.size()) stop("length(w) must match nrow(X)");
      if(ng == 0) {
        NumericVector out = no_init_vector(col);
        if(narm) {
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
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
              M2i /= sumwi-1;
              if(sd) M2i = sqrt(M2i);
              if(std::isnan(M2i)) M2i = NA_REAL;
              out[j] = M2i; // (double)M2i;
            } else out[j] = NA_REAL;
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
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
            M2i /= sumwi-1;
            if(sd) M2i = sqrt(M2i);
            if(std::isnan(M2i)) M2i = NA_REAL;
            out[j] = M2i; // (double)M2i;
          }
        }
        if(drop) out.attr("names") = colnames(x);
        else {
          out.attr("dim") = Dimension(1, col);
          colnames(out) = colnames(x);
        }
        return out;
      } else { // with groups and weights
        if(g.size() != l) stop("length(g) must match nrow(X)");
        if(narm) {
          NumericMatrix M2 = no_init_matrix(ng, col);
          std::fill(M2.begin(), M2.end(), NA_REAL);
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            NumericMatrix::Column M2j = M2( _ , j);
            double d1j = 0; // meanj[ng], sumwj[ng];
            NumericVector meanj = no_init_vector(ng), sumwj = no_init_vector(ng);
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
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] = sqrt(M2j[i]/(sumwj[i]-1));
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] /= sumwj[i]-1;
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            }
          }
          colnames(M2) = colnames(x);
          return M2;
        } else {
          NumericMatrix M2(ng, col);
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            NumericMatrix::Column M2j = M2( _ , j);
            std::vector<double> meanj(ng), sumwj(ng);
            double d1j = 0;
            int ngs = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(M2j[g[i]-1])) continue;
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                M2j[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) break;
              } else {
                sumwj[g[i]-1] += wg[i];
                d1j = column[i] - meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (wg[i] / sumwj[g[i]-1]);
                M2j[g[i]-1] += wg[i] * d1j * (column[i] - meanj[g[i]-1]);
              }
            }
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] = sqrt(M2j[i]/(sumwj[i]-1));
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] /= sumwj[i]-1;
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            }
          }
          colnames(M2) = colnames(x);
          return M2;
        }
      }
    }

  } else { // ONE-PASS METHOD -------------------------------------

    if(Rf_isNull(w)) { // No weights
      if(ng == 0) {
        NumericVector out = no_init_vector(col);
        if(narm) {
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
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
              sq_sumj = (sq_sumj-pow(sumj/nj,2)*nj)/(nj-1);
              if(sd) sq_sumj = sqrt(sq_sumj);
              if(std::isnan(sq_sumj)) sq_sumj = NA_REAL;
              out[j] = (double)sq_sumj;
            } else out[j] = NA_REAL;
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
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
            if(!std::isnan(sq_sumj)) {
              sq_sumj = (sq_sumj-pow(sumj/l,2)*l)/(l-1);
              if(sd) sq_sumj = sqrt(sq_sumj);
              if(std::isnan(sq_sumj)) sq_sumj = NA_REAL;
            }
            out[j] = (double)sq_sumj;
          }
        }
        if(drop) out.attr("names") = colnames(x);
        else {
          out.attr("dim") = Dimension(1, col);
          colnames(out) = colnames(x);
        }
        return out;
      } else { // with groups
        if(g.size() != l) stop("length(g) must match nrow(X)");
        if(narm) {
          NumericMatrix sq_sum = no_init_matrix(ng, col);
          std::fill(sq_sum.begin(), sq_sum.end(), NA_REAL);
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            NumericMatrix::Column sq_sumj = sq_sum( _ , j);
            NumericVector sumj = no_init_vector(ng); // double sumj[ng];
            IntegerVector nj = no_init_vector(ng); // int nj[ng];
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
              if(std::isnan(sq_sumj[i])) continue;
              sq_sumj[i] = (sq_sumj[i] - pow(sumj[i]/nj[i],2)*nj[i])/(nj[i]-1);
              if(sd) sq_sumj[i] = sqrt(sq_sumj[i]);
              if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
            }
          }
          colnames(sq_sum) = colnames(x);
          return sq_sum;
        } else {
          NumericMatrix sq_sum(ng, col);
          if(Rf_isNull(gs)) {
            // int gsv[ng], memsize = sizeof(int)*ng;
            for(int j = col; j--; ) {
              NumericMatrix::ConstColumn column = x( _ , j);
              NumericMatrix::Column sq_sumj = sq_sum( _ , j);
              std::vector<int> gsv(ng); // memset(gsv, 0, memsize);
              std::vector<double> sumj(ng);
              int ngs = 0;
              for(int i = 0; i != l; ++i) {
                if(std::isnan(column[i])) {
                  if(std::isnan(sq_sumj[g[i]-1])) continue;
                  sq_sumj[g[i]-1] = column[i];
                  ++ngs;
                  if(ngs == ng) break;
                } else {
                  sumj[g[i]-1] += column[i];
                  sq_sumj[g[i]-1] += pow(column[i],2);
                  ++gsv[g[i]-1];
                }
              }
              for(int i = ng; i--; ) {
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = (sq_sumj[i] - pow(sumj[i]/gsv[i],2)*gsv[i])/(gsv[i]-1);
                if(sd) sq_sumj[i] = sqrt(sq_sumj[i]);
                if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
              }
            }
          } else {
            IntegerVector gsv = gs;
            if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
            for(int j = col; j--; ) {
              NumericMatrix::ConstColumn column = x( _ , j);
              NumericMatrix::Column sq_sumj = sq_sum( _ , j);
              std::vector<double> sumj(ng);
              int ngs = 0;
              for(int i = 0; i != l; ++i) {
                if(std::isnan(column[i])) {
                  if(std::isnan(sq_sumj[g[i]-1])) continue;
                  sq_sumj[g[i]-1] = column[i];
                  ++ngs;
                  if(ngs == ng) break;
                } else {
                  sumj[g[i]-1] += column[i];
                  sq_sumj[g[i]-1] += pow(column[i],2);
                }
              }
              for(int i = ng; i--; ) {
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = (sq_sumj[i] - pow(sumj[i]/gsv[i],2)*gsv[i])/(gsv[i]-1);
                if(sd) sq_sumj[i] = sqrt(sq_sumj[i]);
                if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
              }
            }
          }
          colnames(sq_sum) = colnames(x);
          return sq_sum;
        }
      }
    } else { // With weights
      NumericVector wg = w;
      if(l != wg.size()) stop("length(w) must match nrow(X)");
      if(ng == 0) {
        NumericVector out = no_init_vector(col);
        if(narm) {
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
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
              sq_sumj = (sq_sumj - pow(sumj/sumwj,2)*sumwj)/(sumwj-1);
              if(sd) sq_sumj = sqrt(sq_sumj);
              if(std::isnan(sq_sumj)) sq_sumj = NA_REAL;
              out[j] = (double)sq_sumj;
            } else out[j] = NA_REAL;
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
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
            if(!std::isnan(sq_sumj)) {
              sq_sumj = (sq_sumj - pow(sumj/sumwj,2)*sumwj)/(sumwj-1);
              if(sd) sq_sumj = sqrt(sq_sumj);
              if(std::isnan(sq_sumj)) sq_sumj = NA_REAL;
            }
            out[j] = (double)sq_sumj;
          }
        }
        if(drop) out.attr("names") = colnames(x);
        else {
          out.attr("dim") = Dimension(1, col);
          colnames(out) = colnames(x);
        }
        return out;
      } else { // with groups and weights
        if(g.size() != l) stop("length(g) must match nrow(X)");
        if(narm) {
          NumericMatrix sq_sum = no_init_matrix(ng, col);
          std::fill(sq_sum.begin(), sq_sum.end(), NA_REAL);
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            NumericMatrix::Column sq_sumj = sq_sum( _ , j);
            NumericVector sumj = no_init_vector(ng), sumwj = no_init_vector(ng); // double sumj[ng], sumwj[ng];
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
              if(std::isnan(sq_sumj[i])) continue;
              sq_sumj[i] = (sq_sumj[i] - pow(sumj[i]/sumwj[i],2)*sumwj[i])/(sumwj[i]-1);
              if(sd) sq_sumj[i] = sqrt(sq_sumj[i]);
              if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
            }
          }
          colnames(sq_sum) = colnames(x);
          return sq_sum;
        } else {
          NumericMatrix sq_sum(ng, col);
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            NumericMatrix::Column sq_sumj = sq_sum( _ , j);
            std::vector<double> sumj(ng), sumwj(ng);
            int ngs = 0;
            for(int i = 0; i != l; ++i) {
              if(std::isnan(sq_sumj[g[i]-1])) continue;
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                sq_sumj[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) break;
              } else {
                sumj[g[i]-1] += column[i]*wg[i];
                sumwj[g[i]-1] += wg[i];
                sq_sumj[g[i]-1] += pow(column[i],2)*wg[i];
              }
            }
            for(int i = ng; i--; ) {
              if(std::isnan(sq_sumj[i])) continue;
              sq_sumj[i] = (sq_sumj[i] - pow(sumj[i]/sumwj[i],2)*sumwj[i])/(sumwj[i]-1);
              if(sd) sq_sumj[i] = sqrt(sq_sumj[i]);
              if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
            }
          }
          colnames(sq_sum) = colnames(x);
          return sq_sum;
        }
      }
    }
  }
}





// [[Rcpp::export]]
SEXP fvarsdlCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
                const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue,
                bool narm = true, bool stable_algo = true,
                bool sd = true, bool drop = true) {
  int l = x.size();

  if(stable_algo) { // WELFORDS ONLINE METHOD -------------------------------------
    if(Rf_isNull(w)) { // No weights
      if(ng == 0) {
        NumericVector out(l);
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            int k = column.size()-1;
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
              M2i /= ni-1;
              if(sd) M2i = sqrt(M2i);
              if(std::isnan(M2i)) M2i = NA_REAL;
              out[j] = M2i; // (double)M2i;
            } else out[j] = NA_REAL;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            int row = column.size();
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
            M2i /= row-1;
            if(sd) M2i = sqrt(M2i);
            if(std::isnan(M2i)) M2i = NA_REAL;
            out[j] = M2i; // (double)M2i;
          }
        }
        if(drop) {
          out.attr("names") = x.attr("names");
          return out;
        } else {
          List res(l);
          for(int j = l; j--; ) {
            res[j] = out[j];
            SHALLOW_DUPLICATE_ATTRIB(res[j], x[j]);
          }
          DUPLICATE_ATTRIB(res, x);
          res.attr("row.names") = 1;
          return res;
        }
      } else { // With groups
        List out(l);
        int gss = g.size();
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector M2j(ng, NA_REAL), nj(ng, 1.0), meanj = no_init_vector(ng);
            double d1j = 0; // meanj[ng]
            // std::vector<double> nj(ng, 1.0);
            for(int i = gss; i--; ) {
              if(std::isnan(column[i])) continue;
              if(std::isnan(M2j[g[i]-1])) {
                meanj[g[i]-1] = column[i];
                M2j[g[i]-1] = 0;
              } else {
                d1j = column[i]-meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
                M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
              }
            }
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] = sqrt(M2j[i]/(nj[i]-1));
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] /= nj[i]-1;
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(M2j, column);
            out[j] = M2j;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector M2j(ng);
            std::vector<double> meanj(ng), nj(ng);
            double d1j = 0;
            int ngs = 0;
            for(int i = 0; i != gss; ++i) {
              if(std::isnan(M2j[g[i]-1])) continue;
              if(std::isnan(column[i])) {
                M2j[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) break;
              } else {
                d1j = column[i]-meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
                M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
              }
            }
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] = sqrt(M2j[i]/(nj[i]-1));
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] /= nj[i]-1;
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(M2j, column);
            out[j] = M2j;
          }
        }
        DUPLICATE_ATTRIB(out, x);
        out.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng); // NumericVector::create(NA_REAL, -ng);
        return out;
      }
    } else { // With weights
      NumericVector wg = w;
      int wgs = wg.size();
      if(ng == 0) {
        NumericVector out(l);
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(column.size() != wgs) stop("length(w) must match nrow(X)");
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
              M2i /= sumwi-1;
              if(sd) M2i = sqrt(M2i);
              if(std::isnan(M2i)) M2i = NA_REAL;
              out[j] = M2i; // (double)M2i;
            } else out[j] = NA_REAL;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(column.size() != wgs) stop("length(w) must match nrow(X)");
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
            M2i /= sumwi-1;
            if(sd) M2i = sqrt(M2i);
            if(std::isnan(M2i)) M2i = NA_REAL;
            out[j] = M2i; // (double)M2i;
          }
        }
        if(drop) {
          out.attr("names") = x.attr("names");
          return out;
        } else {
          List res(l);
          for(int j = l; j--; ) {
            res[j] = out[j];
            SHALLOW_DUPLICATE_ATTRIB(res[j], x[j]);
          }
          DUPLICATE_ATTRIB(res, x);
          res.attr("row.names") = 1;
          return res;
        }
      } else {
        List out(l);
        int gss = g.size();
        if(wgs != gss) stop("length(w) must match length(g)");
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector M2j(ng, NA_REAL), meanj = no_init_vector(ng), sumwj = no_init_vector(ng);
            double d1j = 0; // , sumwj[ng], meanj[ng];
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
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] = sqrt(M2j[i]/(sumwj[i]-1));
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] /= sumwj[i]-1;
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(M2j, column);
            out[j] = M2j;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector M2j(ng);
            std::vector<double> sumwj(ng), meanj(ng);
            double d1j = 0;
            int ngs = 0;
            for(int i = 0; i != gss; ++i) {
              if(std::isnan(M2j[g[i]-1])) continue;
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                M2j[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) break;
              } else {
                sumwj[g[i]-1] += wg[i];
                d1j = column[i] - meanj[g[i]-1];
                meanj[g[i]-1] += d1j * (wg[i] / sumwj[g[i]-1]);
                M2j[g[i]-1] += wg[i] * d1j * (column[i] - meanj[g[i]-1]);
              }
            }
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] = sqrt(M2j[i]/(sumwj[i]-1));
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                else {
                  M2j[i] /= sumwj[i]-1;
                  if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
                }
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(M2j, column);
            out[j] = M2j;
          }
        }
        DUPLICATE_ATTRIB(out, x);
        out.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng); // NumericVector::create(NA_REAL, -ng);
        return out;
      }
    }


  } else { // ONE-PASS METHOD -------------------------------------

    if(Rf_isNull(w)) { // No weights
      if(ng == 0) {
        NumericVector out(l);
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            int k = column.size()-1, ni = 1;
            long double sumi = column[k], sq_sumi = 0;
            while(std::isnan(sumi) && k!=0) sumi = column[--k];
            sq_sumi = sumi*sumi;
            if(k != 0) {
              for(int i = k; i--; ) {
                if(std::isnan(column[i])) continue;
                sumi += column[i];
                sq_sumi += pow(column[i],2);
                ++ni;
              }
              sq_sumi = (sq_sumi-pow(sumi/ni,2)*ni)/(ni-1);
              if(sd) sq_sumi = sqrt(sq_sumi);
              if(std::isnan(sq_sumi)) sq_sumi = NA_REAL;
              out[j] = (double)sq_sumi;
            } else out[j] = NA_REAL;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            long double sumi = 0, sq_sumi = 0;
            int row = column.size();
            for(int i = 0; i != row; ++i) {
              if(std::isnan(column[i])) {
                sq_sumi = NA_REAL;
                break;
              } else {
                sumi += column[i];
                sq_sumi += pow(column[i],2);
              }
            }
            if(!std::isnan(sq_sumi)) {
              sq_sumi = (sq_sumi - pow(sumi/row,2)*row)/(row-1);
              if(sd) sq_sumi = sqrt(sq_sumi);
              if(std::isnan(sq_sumi)) sq_sumi = NA_REAL;
            }
            out[j] = (double)sq_sumi;
          }
        }
        if(drop) {
          out.attr("names") = x.attr("names");
          return out;
        } else {
          List res(l);
          for(int j = l; j--; ) {
            res[j] = out[j];
            SHALLOW_DUPLICATE_ATTRIB(res[j], x[j]);
          }
          DUPLICATE_ATTRIB(res, x);
          res.attr("row.names") = 1;
          return res;
        }
      } else { // With groups
        List out(l);
        int gss = g.size();
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector sq_sumj(ng, NA_REAL), sumj = no_init_vector(ng);
            // double sumj[ng];
            std::vector<int> nj(ng, 1);
            for(int i = gss; i--; ) {
              if(std::isnan(column[i])) continue;
              if(std::isnan(sq_sumj[g[i]-1])) {
                sumj[g[i]-1] = column[i];
                sq_sumj[g[i]-1] = pow(column[i],2);
              } else {
                sumj[g[i]-1] += column[i];
                sq_sumj[g[i]-1] += pow(column[i],2);
                ++nj[g[i]-1];
              }
            }
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = sqrt((sq_sumj[i] - pow(sumj[i]/nj[i],2)*nj[i])/(nj[i]-1));
                if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = (sq_sumj[i] - pow(sumj[i]/nj[i],2)*nj[i])/(nj[i]-1);
                if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(sq_sumj, column);
            out[j] = sq_sumj;
          }
        } else {
          if(Rf_isNull(gs)) {
            // int gsv[ng], memsize = sizeof(int)*ng;
            for(int j = l; j--; ) {
              NumericVector column = x[j];
              if(gss != column.size()) stop("length(g) must match nrow(X)");
              NumericVector sq_sumj(ng), sumj(ng);
              std::vector<int> gsv(ng); // memset(gsv, 0, memsize);
              int ngs = 0;
              for(int i = 0; i != gss; ++i) {
                if(std::isnan(column[i])) {
                  if(std::isnan(sq_sumj[g[i]-1])) continue;
                  sq_sumj[g[i]-1] = column[i];
                  ++ngs;
                  if(ngs == ng) break;
                } else {
                  sumj[g[i]-1] += column[i];
                  sq_sumj[g[i]-1] += pow(column[i],2);
                  ++gsv[g[i]-1];
                }
              }
              if(sd) {
                for(int i = ng; i--; ) {
                  if(std::isnan(sq_sumj[i])) continue;
                  sq_sumj[i] = sqrt((sq_sumj[i] - pow(sumj[i]/gsv[i],2)*gsv[i])/(gsv[i]-1));
                  if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
                }
              } else {
                for(int i = ng; i--; ) {
                  if(std::isnan(sq_sumj[i])) continue;
                  sq_sumj[i] = (sq_sumj[i] - pow(sumj[i]/gsv[i],2)*gsv[i])/(gsv[i]-1);
                  if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
                }
              }
              SHALLOW_DUPLICATE_ATTRIB(sq_sumj, column);
              out[j] = sq_sumj;
            }
          } else {
            IntegerVector gsv = gs;
            if(gsv.size() != ng) stop("Vector of group-sizes must match number of groups");
            for(int j = l; j--; ) {
              NumericVector column = x[j];
              if(gss != column.size()) stop("length(g) must match nrow(X)");
              NumericVector sq_sumj(ng);
              std::vector<double> sumj(ng);
              int ngs = 0;
              for(int i = 0; i != gss; ++i) {
                if(std::isnan(column[i])) {
                  if(std::isnan(sq_sumj[g[i]-1])) continue;
                  sq_sumj[g[i]-1] = column[i];
                  ++ngs;
                  if(ngs == ng) break;
                } else {
                  sumj[g[i]-1] += column[i];
                  sq_sumj[g[i]-1] += pow(column[i],2);
                }
              }
              if(sd) {
                for(int i = ng; i--; ) {
                  if(std::isnan(sq_sumj[i])) continue;
                  sq_sumj[i] = sqrt((sq_sumj[i] - pow(sumj[i]/gsv[i],2)*gsv[i])/(gsv[i]-1));
                  if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
                }
              } else {
                for(int i = ng; i--; ) {
                  if(std::isnan(sq_sumj[i])) continue;
                  sq_sumj[i] = (sq_sumj[i] - pow(sumj[i]/gsv[i],2)*gsv[i])/(gsv[i]-1);
                  if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
                }
              }
              SHALLOW_DUPLICATE_ATTRIB(sq_sumj, column);
              out[j] = sq_sumj;
            }
          }
        }
        DUPLICATE_ATTRIB(out, x);
        out.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng); // NumericVector::create(NA_REAL, -ng);
        return out;
      }
    } else { // With weights
      NumericVector wg = w;
      int wgs = wg.size();
      if(ng == 0) {
        NumericVector out(l);
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(column.size() != wgs) stop("length(w) must match nrow(X)");
            int k = wgs-1;
            while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
            long double sumwi = wg[k], sumi = column[k]*sumwi, sq_sumi = column[k]*sumi;
            if(k != 0) {
              for(int i = k; i--; ) {
                if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
                sumi += column[i]*wg[i];
                sumwi += wg[i];
                sq_sumi += pow(column[i],2)*wg[i];
              }
              sq_sumi = (sq_sumi - pow(sumi/sumwi,2)*sumwi)/(sumwi-1);
              if(sd) sq_sumi = sqrt(sq_sumi);
              if(std::isnan(sq_sumi)) sq_sumi = NA_REAL;
              out[j] = (double)sq_sumi;
            } else out[j] = NA_REAL;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(column.size() != wgs) stop("length(w) must match nrow(X)");
            long double sumi = 0, sumwi = 0, sq_sumi = 0;
            for(int i = 0; i != wgs; ++i) {
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                sq_sumi = NA_REAL;
                break;
              } else {
                sumi += column[i]*wg[i];
                sumwi += wg[i];
                sq_sumi += pow(column[i],2)*wg[i];
              }
            }
            if(!std::isnan(sq_sumi)) {
              sq_sumi = (sq_sumi - pow(sumi/sumwi,2)*sumwi)/(sumwi-1);
              if(sd) sq_sumi = sqrt(sq_sumi);
              if(std::isnan(sq_sumi)) sq_sumi = NA_REAL;
            }
            out[j] = (double)sq_sumi;
          }
        }
        if(drop) {
          out.attr("names") = x.attr("names");
          return out;
        } else {
          List res(l);
          for(int j = l; j--; ) {
            res[j] = out[j];
            SHALLOW_DUPLICATE_ATTRIB(res[j], x[j]);
          }
          DUPLICATE_ATTRIB(res, x);
          res.attr("row.names") = 1;
          return res;
        }
      } else { // With groups and weights
        List out(l);
        int gss = g.size();
        if(wgs != gss) stop("length(w) must match length(g)");
        if(narm) {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector sq_sumj(ng, NA_REAL), sumj = no_init_vector(ng), sumwj = no_init_vector(ng);
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
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = sqrt((sq_sumj[i] - pow(sumj[i]/sumwj[i],2)*sumwj[i])/(sumwj[i]-1));
                if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = (sq_sumj[i] - pow(sumj[i]/sumwj[i],2)*sumwj[i])/(sumwj[i]-1);
                if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(sq_sumj, column);
            out[j] = sq_sumj;
          }
        } else {
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            NumericVector sq_sumj(ng);
            std::vector<double> sumwj(ng), sumj(ng);
            int ngs = 0;
            for(int i = 0; i != gss; ++i) {
              if(std::isnan(sq_sumj[g[i]-1])) continue;
              if(std::isnan(column[i]) || std::isnan(wg[i])) {
                sq_sumj[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) break;
              } else {
                sumj[g[i]-1] += column[i]*wg[i];
                sumwj[g[i]-1] += wg[i];
                sq_sumj[g[i]-1] += pow(column[i],2)*wg[i];
              }
            }
            if(sd) {
              for(int i = ng; i--; ) {
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = sqrt((sq_sumj[i] - pow(sumj[i]/sumwj[i],2)*sumwj[i])/(sumwj[i]-1));
                if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
              }
            } else {
              for(int i = ng; i--; ) {
                if(std::isnan(sq_sumj[i])) continue;
                sq_sumj[i] = (sq_sumj[i] - pow(sumj[i]/sumwj[i],2)*sumwj[i])/(sumwj[i]-1);
                if(std::isnan(sq_sumj[i])) sq_sumj[i] = NA_REAL;
              }
            }
            SHALLOW_DUPLICATE_ATTRIB(sq_sumj, column);
            out[j] = sq_sumj;
          }
        }
        DUPLICATE_ATTRIB(out, x);
        out.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng); // NumericVector::create(NA_REAL, -ng);
        return out;
      }
    }
  }
}
