// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
#include "base_radixsort.h"
}

auto isnan2 = [](double x) { return x != x; };
auto nisnan = [](double x) { return x == x; };

// Note: For changes see fmedian.cpp

//[[Rcpp::export]]
NumericVector fnthCpp(const NumericVector& x, double Q = 0.5, int ng = 0, const IntegerVector& g = 0,
                      const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, bool narm = true, bool upper = false) {
  int l = x.size();
  if(Q <= 0 || Q == 1) stop("Q needs to be positive, and please use fmin or fmax instead of Q = 0 or Q = 1.");
  if(Q > 1) Q = Q/l;
  bool lower = !upper;

  if(Rf_isNull(w)) {
    if(ng == 0) {
      NumericVector out(1);
      if(narm) {
        NumericVector xd = no_init_vector(l);
        auto pend = std::remove_copy_if(x.begin(), x.end(), xd.begin(), isnan2);
        int sz = pend - xd.begin(), nth = lower ? (sz-1)*Q : sz*Q;
        if(sz == 0) return NumericVector::create(x[0]);
        std::nth_element(xd.begin(), xd.begin()+nth, pend);
        out = xd[nth];
      } else {
        for(int i = 0; i != l; ++i) if(isnan2(x[i])) return NumericVector::create(x[i]);
        int nth = lower ? (l-1)*Q : l*Q;
        NumericVector xd = Rf_duplicate(x);
        std::nth_element(xd.begin(), xd.begin()+nth, xd.end());
        out = xd[nth];
      }
      return out;
    } else { // with groups
      if(l != g.size()) stop("length(g) must match length(x)");
      std::vector<std::vector<double> > gmap(ng);
      int ngp = ng+1;
      std::vector<int> gcount(ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != l; ++i) ++gcount[g[i]];
        for(int i = 0; i != ng; ++i) {
          if(gcount[i+1] == 0) stop("group size of 0 encountered");
          gmap[i] = std::vector<double> (gcount[i+1]);
          gcount[i+1] = 0;
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("group size of 0 encountered");
          gmap[i] = std::vector<double>(gsv[i]);
        }
      }

      if(narm) {
        NumericVector out(ng, NA_REAL);
        for(int i = 0; i != l; ++i) if(nisnan(x[i])) gmap[g[i]-1][gcount[g[i]]++] = x[i];
        for(int i = 0; i != ng; ++i) {
          if(gcount[i+1] != 0) {
            int n = gcount[i+1], nth = lower ? (n-1)*Q : n*Q;
            auto begin = gmap[i].begin(), mid = begin + nth;
              std::nth_element(begin, mid, begin + n);
              out[i] = *(mid);
          }
        }
        DUPLICATE_ATTRIB(out, x);
        return out;
      } else {
        NumericVector out(ng);
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(isnan2(x[i])) {
            if(nisnan(out[g[i]-1])) {
              out[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            gmap[g[i]-1][gcount[g[i]]++] = x[i];
          }
        }
        for(int i = 0; i != ng; ++i) {
          if(isnan2(out[i])) continue;
          int n = gcount[i+1], nth = lower ? (n-1)*Q : n*Q;
          auto begin = gmap[i].begin(), mid = begin + nth;
          std::nth_element(begin, mid, begin + n);
          out[i] = *(mid);
        }
        DUPLICATE_ATTRIB(out, x);
        return out;
      }
    }
  } else {
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match length(x)");
    IntegerVector o = no_init_vector(l);
    int *ord = INTEGER(o);
    Cdoubleradixsort(ord, TRUE, FALSE, x); // starts from 1....

    if(ng == 0) {
      double wsumh = 0, sumw = wg[o[0]-1];
      int k = 1;
      if(narm) {
        for(int i = 0; i != l; ++i) if(nisnan(x[i]) && nisnan(wg[i])) wsumh += wg[i];
        wsumh *= Q;
      } else {
        if(isnan2(x[o[l-1]-1])) return NumericVector::create(NA_REAL);
        wsumh = std::accumulate(wg.begin(), wg.end(), 0.0) * Q;
        if(isnan2(wsumh)) return NumericVector::create(NA_REAL);
      }
      if(lower) {
        while(sumw < wsumh) sumw += wg[o[k++]-1];
      } else {
        while(sumw <= wsumh) sumw += wg[o[k++]-1];
      }
      return NumericVector::create(x[o[k-1]-1]);
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
      NumericVector wsumh(ng), wsum(ng), out(ng, NA_REAL);
      if(narm) {
        for(int i = 0; i != l; ++i) if(nisnan(x[i]) && nisnan(wg[i])) wsumh[g[i]-1] += wg[i];
      } else {
        for(int i = 0; i != l; ++i) wsumh[g[i]-1] += wg[i];
      }
      wsumh = wsumh * Q;
      int gi, oi;
      if(lower) {
        for(IntegerVector::const_iterator it = o.begin(), end = o.end(); it != end; ++it) { // const iterator still faster ?
          oi = *it-1;
          gi = g[oi]-1;
          if(wsum[gi] < wsumh[gi]) {
            wsum[gi] += wg[oi];
            out[gi] = x[oi];
          }
        }
      } else {
        for(IntegerVector::const_iterator it = o.begin(), end = o.end(); it != end; ++it) { // const iterator still faster ?
          oi = *it-1;
          gi = g[oi]-1;
          if(wsum[gi] <= wsumh[gi]) {
            wsum[gi] += wg[oi];
            out[gi] = x[oi];
          }
        }
      }
      return out;
    }
  }
}

