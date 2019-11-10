#include <Rcpp.h>
using namespace Rcpp;

// Note: All comments are in fvara.cpp
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
            for(int j = col; j--; ) { 
              NumericMatrix::Column column = x( _ , j);
              NumericMatrix::Column outj = out( _ , j);
              // NumericVector meanj = no_init_vector(ng);
              // NumericVector nj = no_init_vector(ng); 
              double meanj[ng], nj[ng]; // stable and faster ?? 
              // NumericVector M2j(ng, NA_REAL);
              std::vector<double> M2j(ng, NA_REAL);
              long double d1j = 0;
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
              long double d1j = 0;
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
            for(int j = col; j--; ) { 
              NumericMatrix::Column column = x( _ , j); 
              NumericMatrix::Column outj = out( _ , j);
              // NumericVector M2j(ng, NA_REAL);
              // NumericVector meanj = no_init_vector(ng); 
              // NumericVector sumwj = no_init_vector(ng); 
              std::vector<double> M2j(ng, NA_REAL); // faster ?? 
              double meanj[ng], sumwj[ng]; // stable and faster ?? 
              long double d1j = 0;
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
              long double d1j = 0;
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
              // NumericVector sq_sumj(ng, NA_REAL);
              // NumericVector sumj = no_init_vector(ng); 
              // IntegerVector nj = no_init_vector(ng); 
              std::vector<double> sq_sumj(ng, NA_REAL); // faster ?? 
              double sumj[ng]; // stable and faster ?? 
              int nj[ng]; // stable and faster ?? 
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
              int gsv[ng], memsize = sizeof(int)*ng;
              for(int j = col; j--; ) {
                NumericMatrix::Column column = x( _ , j); 
                NumericMatrix::Column outj = out( _ , j);
                std::vector<double> sq_sumj(ng), sumj(ng); // faster than NumericVector ?? 
                memset(gsv, 0, memsize);
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
            for(int j = col; j--; ) { 
              NumericMatrix::Column column = x( _ , j); 
              NumericMatrix::Column outj = out( _ , j);
              // NumericVector sq_sumj(ng, NA_REAL);
              // NumericVector sumj = no_init_vector(ng); 
              // NumericVector sumwj = no_init_vector(ng);
              std::vector<double> sq_sumj(ng, NA_REAL); // faster ?? 
              double sumj[ng], sumwj[ng]; // stable and faster ?? 
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
