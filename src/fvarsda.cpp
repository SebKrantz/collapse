#include <Rcpp.h>
using namespace Rcpp;

// Note: All comments are in fvara.cpp

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
            double ni = 0;
            long double meani = 0, d1i = 0, M2i = 0;
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
              out[j] = (double)M2i;            
            } else out[j] = NA_REAL; 
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            double ni = 0;
            long double meani = 0, d1i = 0, M2i = 0;
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
            out[j] = (double)M2i;
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
            double meanj[ng], nj[ng], d1j = 0;
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
            long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
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
              out[j] = (double)M2i;
            } else out[j] = NA_REAL;
          }
        } else {
          for(int j = col; j--; ) {
            NumericMatrix::ConstColumn column = x( _ , j);
            long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
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
            out[j] = (double)M2i;
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
            double meanj[ng], sumwj[ng], d1j = 0;
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
            double sumj[ng]; 
            int nj[ng]; 
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
            int gsv[ng], memsize = sizeof(int)*ng;
            for(int j = col; j--; ) {
              NumericMatrix::ConstColumn column = x( _ , j); 
              NumericMatrix::Column sq_sumj = sq_sum( _ , j);
              memset(gsv, 0, memsize);
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
            double sumj[ng], sumwj[ng];
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
