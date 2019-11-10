#include <Rcpp.h>
using namespace Rcpp;

// Note: All comments are in fvarl.cpp

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
            NumericVector column = x[j], outj = NULL;
            int row = column.size(), k = row-1;
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
            NumericVector column = x[j], outj = NULL;
            int row = column.size();
            double ni = 0;
            long double meani = 0, d1i = 0, M2i = 0;
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
          for(int j = l; j--; ) { 
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            // NumericVector meanj = no_init_vector(ng);
            // NumericVector nj = no_init_vector(ng); 
            // NumericVector M2j(ng, NA_REAL);
            std::vector<double> M2j(ng, NA_REAL); // faster and stable ?? 
            double meanj[ng], nj[ng]; 
            long double d1j = 0;
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
            long double d1j = 0;
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
            NumericVector column = x[j], outj = NULL;
            if(wgs != column.size()) stop("length(w) must match nrow(X)");
            int k = wgs-1;
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
            NumericVector column = x[j], outj = NULL;
            if(wgs != column.size()) stop("length(w) must match nrow(X)");
            long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
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
        if(narm) {
          for(int j = l; j--; ) { 
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            // NumericVector M2j(ng, NA_REAL);
            // NumericVector meanj = no_init_vector(ng); 
            // NumericVector sumwj = no_init_vector(ng); 
            std::vector<double> M2j(ng, NA_REAL); // faster and stable ?? 
            double meanj[ng], sumwj[ng]; 
            long double d1j = 0;
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
            long double d1j = 0;
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
            NumericVector column = x[j], outj = NULL;
            int row = column.size(), k = row-1, nj = 1;
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
            NumericVector column = x[j], outj = NULL;
            int row = column.size();
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
            // NumericVector sq_sumj(ng, NA_REAL);
            // NumericVector sumj = no_init_vector(ng); 
            // IntegerVector nj = no_init_vector(ng);
            std::vector<double> sq_sumj(ng, NA_REAL); // faster and stable ?? 
            double sumj[ng];
            int nj[ng]; 
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
            int gsv[ng], memsize = sizeof(int)*ng;
            for(int j = l; j--; ) {
              NumericVector column = x[j];
              if(gss != column.size()) stop("length(g) must match nrow(X)");
              int ngs = 0;
              NumericVector outj = no_init_vector(gss);
              std::vector<double> sq_sumj(ng), sumj(ng); // faster than NumericVector??
              memset(gsv, 0, memsize);
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
            NumericVector column = x[j], outj = NULL;
            if(wgs != column.size()) stop("length(w) must match nrow(X)");
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
            NumericVector column = x[j], outj = NULL;
            if(wgs != column.size()) stop("length(w) must match nrow(X)");
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
          for(int j = l; j--; ) { 
            NumericVector column = x[j];
            if(gss != column.size()) stop("length(g) must match nrow(X)");
            // NumericVector sq_sumj(ng, NA_REAL);
            // NumericVector sumj = no_init_vector(ng); 
            // NumericVector sumwj = no_init_vector(ng);
            std::vector<double> sq_sumj(ng, NA_REAL); // faster and stable ?? 
            double sumj[ng], sumwj[ng];
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
