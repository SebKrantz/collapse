#include <Rcpp.h>
using namespace Rcpp;

// Note: All comments are in fvarl.cpp

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
          for(int j = l; j--; ) {
            NumericVector column = x[j];
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
            M2i /= row-1;
            if(sd) M2i = sqrt(M2i);
            if(std::isnan(M2i)) M2i = NA_REAL; 
            out[j] = (double)M2i;
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
            NumericVector M2j(ng, NA_REAL);
            double meanj[ng], d1j = 0;
            std::vector<double> nj(ng, 1.0);
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
        out.attr("row.names") = NumericVector::create(NA_REAL, -ng);
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
          for(int j = l; j--; ) {
            NumericVector column = x[j];
            if(column.size() != wgs) stop("length(w) must match nrow(X)"); 
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
            M2i /= sumwi-1;
            if(sd) M2i = sqrt(M2i);
            if(std::isnan(M2i)) M2i = NA_REAL;
            out[j] = (double)M2i;
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
            NumericVector M2j(ng, NA_REAL);
            double sumwj[ng], meanj[ng], d1j = 0;
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
        out.attr("row.names") = NumericVector::create(NA_REAL, -ng);
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
            NumericVector sq_sumj(ng, NA_REAL);
            double sumj[ng];
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
            int gsv[ng], memsize = sizeof(int)*ng;
            for(int j = l; j--; ) {
              NumericVector column = x[j];
              if(gss != column.size()) stop("length(g) must match nrow(X)");
              NumericVector sq_sumj(ng);
              memset(gsv, 0, memsize);
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
        out.attr("row.names") = NumericVector::create(NA_REAL, -ng);
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
            NumericVector sq_sumj(ng, NA_REAL);
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
        out.attr("row.names") = NumericVector::create(NA_REAL, -ng);
        return out;
      }
    }
  }
}
