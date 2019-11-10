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
          double n = 0; 
          long double mean = 0, d1 = 0, M2 = 0; // LD really necessary ?? what about speed ?? 
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
            return NumericVector::create((double)M2);
          } else return NumericVector::create(NA_REAL); 
        } else {
          double n = 0; 
          long double mean = 0, d1 = 0, M2 = 0;
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
          return NumericVector::create((double)M2);
        }
      } else { // with groups
        if(g.size() != l) stop("length(g) must match nrow(X)");
        long double d1 = 0;
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
          long double sumw = 0, mean = 0, M2 = 0, d1 = 0;
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
            return NumericVector::create((double)M2);
          } else return NumericVector::create(NA_REAL); 
        } else { 
          long double sumw = 0, mean = 0, M2 = 0, d1 = 0;
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
          return NumericVector::create((double)M2);
        }
      } else { // with groups
        if(g.size() != l) stop("length(g) must match nrow(X)");
        long double d1 = 0;
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
          IntegerVector gsv = NULL;
          int ngs = 0;
          if(Rf_isNull(gs)) {
            gsv = IntegerVector(ng);
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
