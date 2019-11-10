#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector BWCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0, 
                    const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, 
                    bool narm = true, bool option = false, bool B = false) {
  int l = x.size();
  NumericVector out = no_init_vector(l);
  
  if (Rf_isNull(w)) { // No weights !!
    if (ng == 0) {
      if(!B && option) stop("For this return option a grouping vector needs to be supplied");
      if(narm) {
        int j = l-1, n = 1; // 1 because for-loop starts from 2
        long double sum = x[j]; 
        while(std::isnan(sum) && j!=0) sum = x[--j]; 
        if(j != 0) for(int i = j; i--; ) {
          if(std::isnan(x[i])) continue;
          sum += x[i]; // Fastest ??
          ++n;
        } 
        sum = sum/n;
        if(B) {
          if(option) std::fill(out.begin(), out.end(), (double)sum); // fastest ?? -> yes !!
          else {
            for(int i = 0; i != l; ++i) {
              if(std::isnan(x[i])) out[i] = x[i];
              else out[i] = sum; // double conversion -> nope, slower !!
            }
          }
        } else {
          out = x - sum; // conversion to double not necessary !!
        }
      } else {
        long double sum = 0;
        for(int i = 0; i != l; ++i) { 
          if(std::isnan(x[i])) {
            sum = x[i]; 
            break;
          } else { 
            sum += x[i];
          }
        }
        sum = sum/l;
        if(B) {
          std::fill(out.begin(), out.end(), (double)sum); // fastes ?? 
        } else {
          out = x - sum; // conversion to double not necessary !!
        }
      }
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      if(narm) {
        NumericVector sum(ng, NA_REAL); // Other way ??
        IntegerVector n(ng, 1); // could also do no_init_vector and then add n[g[i]-1] = 1 in fir if condition... -> Nope, that is slower !!!
        for(int i = l; i--; ) { 
          if(!std::isnan(x[i])) { // faster way to code this ??? -> Not Bad at all -> index for g[i]-1?? -> Nope, no noticeable improvement !!
            if(std::isnan(sum[g[i]-1])) sum[g[i]-1] = x[i];
            else {
              sum[g[i]-1] += x[i]; 
              ++n[g[i]-1];
            }
          }
        }
        if(B) {
          for(int i = ng; i--; ) sum[i] /= n[i];
          if(option) {
            for(int i = 0; i != l; ++i) out[i] = sum[g[i]-1];
          } else {
            for(int i = 0; i != l; ++i) {
              if(std::isnan(x[i])) out[i] = x[i];
              else out[i] = sum[g[i]-1]; 
            }
          }
        } else {
          if(!option) {
            for(int i = ng; i--; ) sum[i] /= n[i]; // faster using two loops? or combine ? -> two loos (this solution) is a lot faster !!!!!!!
            for(int i = 0; i != l; ++i) out[i] = x[i] - sum[g[i]-1]; // best loop ?? -> just as fast as the other one !!
          } else {
            int on = 0;
            double osum = 0;
            for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
              if(std::isnan(sum[i])) continue; // solves the issue !!
              osum += sum[i];
              on += n[i];
              sum[i] /= n[i]; // fastest ?? 
            }
            osum = osum/on;
            for(int i = 0; i != l; ++i) out[i] = x[i] - sum[g[i]-1] + osum;
          }
        }
      } else {
        NumericVector sum(ng); // no_init_vector // good?? -> yes, but not initializing is numerically unstable.. 
        IntegerVector gsv = NULL;
        int ngs = 0;
        if(Rf_isNull(gs)) {
          gsv = IntegerVector(ng);
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
        } else {
          gsv = gs;
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
        }
        if(B) {
          for(int i = ng; i--; ) sum[i] /= gsv[i];
          for(int i = 0; i != l; ++i) out[i] = sum[g[i]-1];
        } else {
          if(!option) {
            for(int i = ng; i--; ) sum[i] /= gsv[i]; 
            for(int i = 0; i != l; ++i) out[i] = x[i] - sum[g[i]-1]; 
          } else {
            int on = 0;
            double osum = 0;
            for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
              if(std::isnan(sum[i])) continue; // solves the issue !!
              osum += sum[i];
              on += gsv[i];
              sum[i] /= gsv[i]; // fastest ?? 
            }
            osum = osum/on;
            for(int i = 0; i != l; ++i) out[i] = x[i] - sum[g[i]-1] + osum;
          }
        }
      }
    }
  } else { // With weights
    NumericVector wg = w; // wg(w) Identical speed
    if(l != wg.size()) stop("length(w) must match length(x)");
    if (ng == 0) {
      if(!B && option) stop("For this return option a grouping vector needs to be supplied");
      if(narm) {
        int j = l-1; // 1 because for-loop starts from 2
        while((std::isnan(x[j]) || std::isnan(wg[j])) && j!=0) --j; // This does not make a difference in performance but is more parsimonious. 
        long double sum = x[j]*wg[j], sumw = wg[j]; 
        if(j != 0) for(int i = j; i--; ) {
          if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
          sum += x[i]*wg[i]; // Fastest ??
          sumw += wg[i];
        } 
        sum = sum/sumw;
        if(B) {
          if(option) std::fill(out.begin(), out.end(), (double)sum); // fastes ?? 
          else {
            for(int i = 0; i != l; ++i) {
              if(std::isnan(x[i])) out[i] = x[i];
              else out[i] = sum; // double conversion ?? 
            }
          }
        } else {
          out = x - sum; // conversion to double not necessary !!
        }
      } else {
        long double sum = 0, sumw = 0;
        for(int i = 0; i != l; ++i) { 
          if(std::isnan(x[i]) || std::isnan(wg[i])) { // good, check both ?? -> yes!!
            sum = x[i]+wg[i]; 
            break;
          } else { 
            sum += x[i]*wg[i];
            sumw += wg[i];
          }
        }
        sum = sum/sumw;
        if(B) {
          std::fill(out.begin(), out.end(), (double)sum); // fastes ?? 
        } else {
          out = x - sum; // conversion to double not necessary !!
        }
      }
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      if(narm) {
        NumericVector sum(ng, NA_REAL); // Other way ?? -> Nope, this is as good as it gets !!
        NumericVector sumw = no_init_vector(ng); // what if only NA ?? -> Works !! for some reason no problem !!, and faster !!
        for(int i = l; i--; ) { 
          if(std::isnan(x[i]) || std::isnan(wg[i])) continue; // faster way to code this ??? -> Not Bad at all -> index for g[i]-1?? -> Nope, no noticeable improvement !!
          if(std::isnan(sum[g[i]-1])) {
            sum[g[i]-1] = x[i]*wg[i];
            sumw[g[i]-1] = wg[i];
          } else {
            sum[g[i]-1] += x[i]*wg[i]; 
            sumw[g[i]-1] += wg[i];
          }
        }
        if(B) {
          sum = sum/sumw;
          if(option) {
            for(int i = 0; i != l; ++i) out[i] = sum[g[i]-1];
          } else {
            for(int i = 0; i != l; ++i) {
              if(std::isnan(x[i])) out[i] = x[i];
              else out[i] = sum[g[i]-1]; 
            }
          }
        } else {
          if(!option) {
            sum = sum/sumw; 
            for(int i = 0; i != l; ++i) out[i] = x[i] - sum[g[i]-1]; 
          } else {
            double osum = 0, osumw = 0;
            for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
              if(std::isnan(sum[i])) continue; // solves the issue !!
              osum += sum[i];
              osumw += sumw[i];
              sum[i] /= sumw[i]; // fastest ?? 
            }
            osum = osum/osumw;
            for(int i = 0; i != l; ++i) out[i] = x[i] - sum[g[i]-1] + osum;
          }
        }
      } else {
        NumericVector sum(ng); // good?? -> yes !! //  = no_init_vector// Not initializing numerically unstable !!
        NumericVector sumw(ng); // Necessary !!
        int ngs = 0;
        for(int i = 0; i != l; ++i) { 
          if(std::isnan(x[i]) || std::isnan(wg[i])) { 
            if(!std::isnan(sum[g[i]-1])) {
              sum[g[i]-1] = sumw[g[i]-1] = x[i]+wg[i]; // or NA_REAL ?? -> Nope, good !!
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            sum[g[i]-1] += x[i]*wg[i];
            sumw[g[i]-1] += wg[i];
          }
        }
        if(B) {
          sum = sum/sumw;
          for(int i = 0; i != l; ++i) out[i] = sum[g[i]-1];
        } else {
          if(!option) {
            sum = sum/sumw; 
            for(int i = 0; i != l; ++i) out[i] = x[i] - sum[g[i]-1]; 
          } else {
            double osum = 0, osumw = 0;
            for(int i = ng; i--; ) { // Problem: if one sum remained NA, osum becomes NA !!!
              if(std::isnan(sum[i])) continue; // solves the issue !!
              osum += sum[i];
              osumw += sumw[i];
              sum[i] /= sumw[i]; // fastest ?? 
            }
            osum = osum/osumw;
            for(int i = 0; i != l; ++i) out[i] = x[i] - sum[g[i]-1] + osum;
          }
        }
      }
    }
  }
  DUPLICATE_ATTRIB(out, x);
  return out;
}
