// [[Rcpp::plugins(cpp11)]]
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
#include "base_radixsort.h"
}

auto isnan2 = [](double x) { return x != x; };
auto nisnan = [](double x) { return x == x; };


//[[Rcpp::export]]
NumericVector fnthCpp(const NumericVector& x, double Q = 0.5, int ng = 0, const IntegerVector& g = 0,
                      const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, bool narm = true, int ret = 1) {
  int l = x.size();
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  bool tiesmean, lower;
  if(Q <= 0 || Q == 1) stop("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
  if(Q > 1) {
    tiesmean = false;
    lower = true;
    if(ng == 0) {
      if(Q >= l) stop("n needs to be between 0 and 1, or between 1 and length(x). Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(l-1);
    } else {
      if(Q >= l/ng) stop("n needs to be between 0 and 1, or between 1 and the length(x)/ng, with ng the number of groups. Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(l/ng-1);
    }
  } else {
    tiesmean = ret == 1;
    lower = ret != 3;
  }

  if(Rf_isNull(w)) {
    if(ng == 0) {
      NumericVector out(1);
      if(narm) {
        NumericVector xd = no_init_vector(l);
        auto pend = std::remove_copy_if(x.begin(), x.end(), xd.begin(), isnan2);
        int sz = pend - xd.begin(), nth = lower ? (sz-1)*Q : sz*Q; // return NumericVector::create(sz, (sz-1)*Q, sz*Q-1, sz*Q);
        if(sz == 0) return Rf_ScalarReal(x[0]);
        std::nth_element(xd.begin(), xd.begin()+nth, pend);
        out = (tiesmean && sz%2 == 0) ? (xd[nth] + *(std::min_element(xd.begin()+nth+1, pend)))*0.5 : xd[nth];
      } else {
        for(int i = 0; i != l; ++i) if(isnan2(x[i])) return Rf_ScalarReal(x[i]);
        NumericVector xd = Rf_duplicate(x);
        int nth = lower ? (l-1)*Q : l*Q;
        std::nth_element(xd.begin(), xd.begin()+nth, xd.end());
        out = (tiesmean && l%2 == 0) ? (xd[nth] + *(std::min_element(xd.begin()+nth+1, xd.end())))*0.5 : xd[nth];
      }
      return out;
    } else { // with groups
      if(l != g.size()) stop("length(g) must match length(x)");
      int ngp = ng+1;
      std::vector<std::vector<double> > gmap(ngp);
      std::vector<int> gcount(ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != l; ++i) ++gcount[g[i]];
        for(int i = 1; i != ngp; ++i) {
          if(gcount[i] == 0) stop("Group size of 0 encountered. This is probably because of unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i] = std::vector<double> (gcount[i]);
          gcount[i] = 0;
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("Group size of 0 encountered. This is probably because of unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i+1] = std::vector<double> (gsv[i]);
        }
      }

      if(narm) {
        NumericVector out(ng, NA_REAL);
        for(int i = 0; i != l; ++i) if(nisnan(x[i])) gmap[g[i]][gcount[g[i]]++] = x[i];
        for(int i = 1; i != ngp; ++i) {
          if(gcount[i] != 0) {
            int n = gcount[i], nth = lower ? (n-1)*Q : n*Q;
            auto begin = gmap[i].begin(), mid = begin + nth, end = begin + n;
            std::nth_element(begin, mid, end);
            out[i-1] = (tiesmean && n%2 == 0) ? (*(mid) + *(std::min_element(mid+1, end)))*0.5 : *(mid);
          }
        }
        if(!Rf_isObject(x)) Rf_copyMostAttrib(x, out);
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
            gmap[g[i]][gcount[g[i]]++] = x[i];
          }
        }
        for(int i = 0; i != ng; ++i) {
          if(isnan2(out[i])) continue;
          int n = gcount[i+1], nth = lower ? (n-1)*Q : n*Q;
          auto begin = gmap[i+1].begin(), mid = begin + nth, end = begin + n;
          std::nth_element(begin, mid, end);
          out[i] = (tiesmean && n%2 == 0) ? (*(mid) + *(std::min_element(mid+1, end)))*0.5 : *(mid);
        }
        if(!Rf_isObject(x)) Rf_copyMostAttrib(x, out);
        return out;
      }
    }
  } else { // with weights
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match length(x)");
    IntegerVector o = no_init_vector(l);
    int *ord = INTEGER(o);
    Cdoubleradixsort(ord, TRUE, FALSE, x); // starts from 1

    if(ng == 0) {
      double wsumQ = 0, wsum = wg[o[0]-1];
      int k = 1;
      if(narm) {
        for(int i = 0; i != l; ++i) if(nisnan(x[i])) wsumQ += wg[i]; //  && nisnan(wg[i])
        if(wsumQ == 0) return Rf_ScalarReal(NA_REAL);
        wsumQ *= Q;
      } else {
        if(isnan2(x[o[l-1]-1])) return Rf_ScalarReal(NA_REAL);
        wsumQ = std::accumulate(wg.begin(), wg.end(), 0.0) * Q;
      }
      if(isnan2(wsumQ)) stop("Missing weights in order statistics are currently only supported if x is also missing"); // return Rf_ScalarReal(NA_REAL);
      if(lower) {
        while(wsum < wsumQ) wsum += wg[o[k++]-1];
        if(tiesmean && wsum == wsumQ) {
          double out = x[o[k-1]-1], n = 2;
          while(wg[o[k]-1] == 0) {
            out += x[o[k++]-1];
            ++n;
          }
          return Rf_ScalarReal((out + x[o[k]-1]) / n);
        }
      } else {
        while(wsum <= wsumQ) wsum += wg[o[k++]-1];
      }
      return Rf_ScalarReal(x[o[k-1]-1]);
    } else { // with groups and weights
      if(l != g.size()) stop("length(g) must match length(x)");
      NumericVector wsumQ(ng), wsum(ng), out(ng, NA_REAL);
      if(narm) {
        for(int i = 0; i != l; ++i) if(nisnan(x[i])) wsumQ[g[i]-1] += wg[i]; //  && nisnan(wg[i])
        for(int i = ng; i--; ) {
          if(isnan2(wsumQ[i])) stop("Missing weights in order statistics are currently only supported if x is also missing");
          wsumQ[i] *= Q;
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(nisnan(x[i])) {
            wsumQ[g[i]-1] += wg[i];
          } else {
            wsum[g[i]-1] = DBL_MAX; // OK ? // needed ? // = wsumQ[g[i]-1]
          }
        }
        for(int i = ng; i--; ) {
          if(isnan2(wsumQ[i]) && wsum[i] != DBL_MAX) stop("Missing weights in order statistics are currently only supported if x is also missing");
          wsumQ[i] *= Q;
        }
      }
      // wsumQ = wsumQ * Q;
      int gi, oi;
      if(tiesmean) {
        std::vector<bool> seen(ng);
        std::vector<int> n(ng, 1); // only needed for 0 weights. Could check above if there are any.
        for(int i = 0; i != l; ++i) {
          oi = o[i]-1;
          gi = g[oi]-1;
          if(seen[gi]) continue;
          if(wsum[gi] < wsumQ[gi]) out[gi] = x[oi];
          else {
            if(wsum[gi] > wsumQ[gi]) {
              seen[gi] = true;
              continue;
            }
            out[gi] += (x[oi]-out[gi])/++n[gi]; // https://stackoverflow.com/questions/28820904/how-to-efficiently-compute-average-on-the-fly-moving-average
          }
          wsum[gi] += wg[oi];
        }
      } else if(lower) {
        for(int i = 0; i != l; ++i) {
          oi = o[i]-1;
          gi = g[oi]-1;
          if(wsum[gi] < wsumQ[gi]) {
            wsum[gi] += wg[oi];
            out[gi] = x[oi];
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          oi = o[i]-1;
          gi = g[oi]-1;
          if(wsum[gi] <= wsumQ[gi]) {
            wsum[gi] += wg[oi];
            out[gi] = x[oi];
          }
        }
      }
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, out);
      return out;
    }
  }
}





//[[Rcpp::export]]
SEXP fnthmCpp(const NumericMatrix& x, double Q = 0.5, int ng = 0, const IntegerVector& g = 0,
              const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, bool narm = true,
              bool drop = true, int ret = 1) {
  int l = x.nrow(), col = x.ncol();
  bool tiesmean, lower;
  if(Q <= 0 || Q == 1) stop("n needs to be between 0 and 1, or between 1 and nrow(x). Use fmin and fmax for minima and maxima.");
  if(Q > 1) {
    tiesmean = false;
    lower = true;
    if(ng == 0) {
      if(Q >= l) stop("n needs to be between 0 and 1, or between 1 and nrow(x). Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(l-1);
    } else {
      if(Q >= l/ng) stop("n needs to be between 0 and 1, or between 1 and the nrow(x)/ng, with ng the number of groups. Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(l/ng-1);
    }
  } else {
    tiesmean = ret == 1;
    lower = ret != 3;
  }

  if(Rf_isNull(w)) {
    if(ng == 0) {
      NumericVector out = no_init_vector(col);
      if(narm) {
        NumericVector column = no_init_vector(l);
        auto begin = column.begin();
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn colj = x( _ , j);
          auto pend = std::remove_copy_if(colj.begin(), colj.end(), begin, isnan2);
          int sz = pend - begin, nth = lower ? (sz-1)*Q : sz*Q;
          if(sz == 0) {
            out[j] = colj[0];
            continue;
          }
          std::nth_element(begin, begin+nth, pend);
          out[j] = (tiesmean && sz%2 == 0) ? (column[nth] + *(std::min_element(begin+nth+1, pend)))*0.5 : column[nth];
        }
      } else {
        int nth = lower ? (l-1)*Q : l*Q;
        bool tm = tiesmean && l%2 == 0;
        for(int j = col; j--; ) {
          {
            NumericMatrix::ConstColumn colj = x( _ , j);
            for(int i = 0; i != l; ++i) {
              if(isnan2(colj[i])) {
                out[j] = colj[i];
                goto endloop;
              }
            }
            NumericVector column = Rf_duplicate(wrap(colj)); // best ?
            std::nth_element(column.begin(), column.begin()+nth, column.end());
            out[j] = tm ? (column[nth] + *(std::min_element(column.begin()+nth+1, column.end())))*0.5 : column[nth];
          }
          endloop:;
        }
      }
      if(drop) Rf_setAttrib(out, R_NamesSymbol, colnames(x));
      else {
        Rf_dimgets(out, Dimension(1, col));
        colnames(out) = colnames(x);
        if(!Rf_isObject(x)) Rf_copyMostAttrib(x, out);
      }
      return out;
    } else { // with groups
      if(l != g.size()) stop("length(g) must match nrow(x)");
      int ngp = ng+1;
      std::vector<std::vector<double> > gmap(ngp);
      std::vector<int> gcount(ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != l; ++i) ++gcount[g[i]];
        for(int i = 1; i != ngp; ++i) {
          if(gcount[i] == 0) stop("Group size of 0 encountered. This is probably because of unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i] = std::vector<double> (gcount[i]);
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("Group size of 0 encountered. This is probably because of unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i+1] = std::vector<double> (gsv[i]);
        }
      }

      if(narm) {
        NumericMatrix out = no_init_matrix(ng, col);
        std::fill(out.begin(), out.end(), NA_REAL);
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j);
          NumericMatrix::Column nthj = out( _ , j);
          gcount.assign(ngp, 0);
          for(int i = 0; i != l; ++i) if(nisnan(column[i])) gmap[g[i]][gcount[g[i]]++] = column[i];
          for(int i = 1; i != ngp; ++i) {
            if(gcount[i] != 0) {
              int n = gcount[i], nth = lower ? (n-1)*Q : n*Q;
              auto begin = gmap[i].begin(), mid = begin + nth, end = begin + n;
              std::nth_element(begin, mid, end);
              nthj[i-1] = (tiesmean && n%2 == 0) ? (*(mid) + *(std::min_element(mid+1, end)))*0.5 : *(mid);
            }
          }
        }
        colnames(out) = colnames(x);
        if(!Rf_isObject(x)) Rf_copyMostAttrib(x, out);
        return out;
      } else {
        NumericMatrix out(ng, col); // no init numerically unstable
        for(int j = col; j--; ) {
          NumericMatrix::ConstColumn column = x( _ , j);
          NumericMatrix::Column nthj = out( _ , j);
          gcount.assign(ngp, 0);
          int ngs = 0;
          for(int i = 0; i != l; ++i) {
            if(isnan2(column[i])) {
              if(nisnan(nthj[g[i]-1])) {
                nthj[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) break;
              }
            } else {
              gmap[g[i]][gcount[g[i]]++] = column[i];
            }
          }
          for(int i = 0; i != ng; ++i) {
            if(isnan2(nthj[i])) continue;
            int n = gcount[i+1], nth = lower ? (n-1)*Q : n*Q;
            auto begin = gmap[i+1].begin(), mid = begin + nth, end = begin + n;
            std::nth_element(begin, mid, end);
            nthj[i] = (tiesmean && n%2 == 0) ? (*(mid) + *(std::min_element(mid+1, end)))*0.5 : *(mid);
          }
        }
        colnames(out) = colnames(x);
        if(!Rf_isObject(x)) Rf_copyMostAttrib(x, out);
        return out;
      }
    }
  } else { // with weights
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match nrow(x)");
    IntegerVector o = no_init_vector(l);
    int *ord = INTEGER(o);

    if(ng == 0) {
      NumericVector out = no_init_vector(col);
      double wsumQ = 0;
      if(!narm) {
        wsumQ = std::accumulate(wg.begin(), wg.end(), 0.0) * Q;
        if(isnan2(wsumQ)) {
          stop("Missing weights in order statistics are currently only supported if x is also missing");
          // std::fill(out.begin(), out.end(), NA_REAL);
          // goto outnth;
        }
      }
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        Cdoubleradixsort(ord, TRUE, FALSE, wrap(column)); // starts from 1....
        if(narm) {
          wsumQ = 0;
          for(int i = 0; i != l; ++i) if(nisnan(column[i])) wsumQ += wg[i]; //  && nisnan(wg[i])
          if(wsumQ == 0) {
            out[j] = NA_REAL;
            continue;
          }
          if(isnan2(wsumQ)) stop("Missing weights in order statistics are currently only supported if x is also missing");
          wsumQ *= Q;
        } else {
          if(isnan2(column[o[l-1]-1])) {
            out[j] = NA_REAL;
            continue;
          }
        }
        double wsum = wg[o[0]-1];
        int k = 1; // what about all missing ? -> gives dsort error...

        if(lower) {
          while(wsum < wsumQ) wsum += wg[o[k++]-1];
          if(tiesmean && wsum == wsumQ) {
            double outtmp = column[o[k-1]-1], n = 2;
            while(wg[o[k]-1] == 0) {
              outtmp += column[o[k++]-1];
              ++n;
            }
            out[j] = (outtmp + column[o[k]-1]) / n;
            continue;
          }
        } else {
          while(wsum <= wsumQ) wsum += wg[o[k++]-1];
        }
        out[j] = column[o[k-1]-1];
      }
      // outnth:
        if(drop) Rf_setAttrib(out, R_NamesSymbol, colnames(x));
        else {
          Rf_dimgets(out, Dimension(1, col));
          colnames(out) = colnames(x);
          if(!Rf_isObject(x)) Rf_copyMostAttrib(x, out);
        }
        return out;
    } else { // with groups and weights
      if(l != g.size()) stop("length(g) must match nrow(x)");
      int gi, oi;
      NumericMatrix out = no_init_matrix(ng, col);
      std::fill(out.begin(), out.end(), NA_REAL);
      NumericVector wsumQ(ng);
      if(!narm) {
        for(int i = 0; i != l; ++i) wsumQ[g[i]-1] += wg[i];
        wsumQ = wsumQ * Q;
      }
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column nthj = out( _ , j);
        Cdoubleradixsort(ord, TRUE, FALSE, wrap(column));
        NumericVector wsum(ng);
        if(narm) {
          std::fill(wsumQ.begin(), wsumQ.end(), 0.0);
          for(int i = 0; i != l; ++i) if(nisnan(column[i])) wsumQ[g[i]-1] += wg[i]; //  && nisnan(wg[i])
          for(int i = ng; i--; ) {
            if(isnan2(wsumQ[i])) stop("Missing weights in order statistics are currently only supported if x is also missing");
            wsumQ[i] *= Q;
          }
        } else {
          for(int i = 0; i != l; ++i) if(isnan2(column[i])) wsum[g[i]-1] = DBL_MAX; // OK ?? // needed??
          for(int i = ng; i--; ) if(isnan2(wsumQ[i]) && wsum[i] != DBL_MAX) stop("Missing weights in order statistics are currently only supported if x is also missing");
        }

        if(tiesmean) {
          std::vector<bool> seen(ng);
          std::vector<int> n(ng, 1); // only needed for 0 weights. Could check above if there are any.
          for(int i = 0; i != l; ++i) {
            oi = o[i]-1;
            gi = g[oi]-1;
            if(seen[gi]) continue;
            if(wsum[gi] < wsumQ[gi]) nthj[gi] = column[oi];
            else {
              if(wsum[gi] > wsumQ[gi]) {
                seen[gi] = true;
                continue;
              }
              nthj[gi] += (column[oi]-nthj[gi])/++n[gi];
            }
            wsum[gi] += wg[oi];
          }
        } else if(lower) {
          for(int i = 0; i != l; ++i) {
            oi = o[i]-1;
            gi = g[oi]-1;
            if(wsum[gi] < wsumQ[gi]) {
              wsum[gi] += wg[oi];
              nthj[gi] = column[oi];
            }
          }
        } else {
          for(int i = 0; i != l; ++i) {
            oi = o[i]-1;
            gi = g[oi]-1;
            if(wsum[gi] <= wsumQ[gi]) {
              wsum[gi] += wg[oi];
              nthj[gi] = column[oi];
            }
          }
        }
      }
      colnames(out) = colnames(x);
      if(!Rf_isObject(x)) Rf_copyMostAttrib(x, out);
      return out;
    }
  }
}





//[[Rcpp::export]]
SEXP fnthlCpp(const List& x, double Q = 0.5, int ng = 0, const IntegerVector& g = 0,
              const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, bool narm = true,
              bool drop = true, int ret = 1) {
  int l = x.size(), lx1 = Rf_length(x[0]);
  bool tiesmean, lower;
  if(Q <= 0 || Q == 1) stop("n needs to be between 0 and 1, or between 1 and nrow(x). Use fmin and fmax for minima and maxima.");
  if(Q > 1) {
    tiesmean = false;
    lower = true;
    if(ng == 0) {
      if(Q >= lx1) stop("n needs to be between 0 and 1, or between 1 and nrow(x). Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(lx1-1);
    } else {
      if(Q >= lx1/ng) stop("n needs to be between 0 and 1, or between 1 and the nrow(x)/ng, with ng the number of groups. Use fmin and fmax for minima and maxima.");
      Q = (Q-1)/(lx1/ng-1);
    }
  } else {
    tiesmean = ret == 1;
    lower = ret != 3;
  }

  if(Rf_isNull(w)) {
    if(ng == 0) {
      NumericVector out = no_init_vector(l);
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector colj = x[j];
          NumericVector column = no_init_vector(colj.size());
          auto begin = column.begin();
          auto pend = std::remove_copy_if(colj.begin(), colj.end(), begin, isnan2);
          int sz = pend - begin, nth = lower ? (sz-1)*Q : sz*Q;
          if(sz == 0) {
            out[j] = colj[0];
            continue;
          }
          std::nth_element(begin, begin+nth, pend);
          out[j] = (tiesmean && sz%2 == 0) ? (column[nth] + *(std::min_element(begin+nth+1, pend)))*0.5 : column[nth];
        }
      } else {
        for(int j = l; j--; ) {
          {
            NumericVector colj = x[j];
            int row = colj.size(), nth = lower ? (row-1)*Q : row*Q;
            for(int i = 0; i != row; ++i) {
              if(isnan2(colj[i])) {
                out[j] = colj[i];
                goto endloop;
              }
            }
            NumericVector column = Rf_duplicate(colj);
            auto begin = column.begin();
            std::nth_element(begin, begin+nth, column.end());
            out[j] = (tiesmean && row%2 == 0) ? (column[nth] + *(std::min_element(begin+nth+1, column.end())))*0.5 : column[nth];
          }
          endloop:;
        }
      }
      if(drop) {
        Rf_setAttrib(out, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
        return out;
      } else {
        List out(l);
        for(int j = l; j--; ) {
          out[j] = out[j];
          DUPLICATE_ATTRIB(out[j], x[j]);
        }
        DUPLICATE_ATTRIB(out, x);
        Rf_setAttrib(out, R_RowNamesSymbol, Rf_ScalarInteger(1));
        return out;
      }

    } else { // with groups
      if(lx1 != g.size()) stop("length(g) must match nrow(x)");
      List out(l);
      int ngp = ng+1;
      std::vector<std::vector<double> > gmap(ngp);
      std::vector<int> gcount(ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != lx1; ++i) ++gcount[g[i]];
        for(int i = 1; i != ngp; ++i) {
          if(gcount[i] == 0) stop("Group size of 0 encountered. This is probably because of unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i] = std::vector<double> (gcount[i]);
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("Group size of 0 encountered. This is probably because of unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i+1] = std::vector<double> (gsv[i]);
        }
      }
      if(narm) {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          NumericVector nthj(ng, NA_REAL);
          if(lx1 != column.size()) stop("length(g) must match nrow(x)");
          gcount.assign(ngp, 0);
          for(int i = 0; i != lx1; ++i) if(nisnan(column[i])) gmap[g[i]][gcount[g[i]]++] = column[i];
          for(int i = 1; i != ngp; ++i) {
            if(gcount[i] != 0) {
              int n = gcount[i], nth = lower ? (n-1)*Q : n*Q;
              auto begin = gmap[i].begin(), mid = begin + nth, end = begin + n;
              std::nth_element(begin, mid, end);
              nthj[i-1] = (tiesmean && n%2 == 0) ? (*(mid) + *(std::min_element(mid+1, end)))*0.5 : *(mid);
            }
          }
          DUPLICATE_ATTRIB(nthj, column);
          out[j] = nthj;
        }
      } else {
        for(int j = l; j--; ) {
          NumericVector column = x[j];
          NumericVector nthj(ng); // no init numerically unstable
          if(lx1 != column.size()) stop("length(g) must match nrow(x)");
          gcount.assign(ngp, 0);
          int ngs = 0;
          for(int i = 0; i != lx1; ++i) {
            if(isnan2(column[i])) {
              if(nisnan(nthj[g[i]-1])) {
                nthj[g[i]-1] = NA_REAL;
                ++ngs;
                if(ngs == ng) break;
              }
            } else {
              gmap[g[i]][gcount[g[i]]++] = column[i];
            }
          }
          for(int i = 0; i != ng; ++i) {
            if(isnan2(nthj[i])) continue;
            int n = gcount[i+1], nth = lower ? (n-1)*Q : n*Q;
            auto begin = gmap[i+1].begin(), mid = begin + nth, end = begin + n;
            std::nth_element(begin, mid, end);
            nthj[i] = (tiesmean && n%2 == 0) ? (*(mid) + *(std::min_element(mid+1, end)))*0.5 : *(mid);
          }
          DUPLICATE_ATTRIB(nthj, column);
          out[j] = nthj;
        }
      }
      DUPLICATE_ATTRIB(out, x);
      Rf_setAttrib(out, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -ng));
      return out;
    }
  } else { // with weights
    NumericVector wg = w;
    if(lx1 != wg.size()) stop("length(w) must match nrow(x)");
    IntegerVector o = no_init_vector(lx1);
    int *ord = INTEGER(o);

    if(ng == 0) {
      NumericVector out = no_init_vector(l);
      double wsumQ = 0;
      if(!narm) {
        wsumQ = std::accumulate(wg.begin(), wg.end(), 0.0) * Q;
        if(isnan2(wsumQ)) {
          stop("Missing weights in order statistics are currently only supported if x is also missing");
          // std::fill(out.begin(), out.end(), NA_REAL);
          // goto outnth;
        }
      }
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(lx1 != column.size()) stop("length(w) must match nrow(x)");
        Cdoubleradixsort(ord, TRUE, FALSE, column); // starts from 1
        if(narm) {
          wsumQ = 0;
          for(int i = 0; i != lx1; ++i) if(nisnan(column[i])) wsumQ += wg[i]; //  && nisnan(wg[i])
          if(wsumQ == 0) {
            out[j] = NA_REAL;
            continue;
          }
          if(isnan2(wsumQ)) stop("Missing weights in order statistics are currently only supported if x is also missing");
          wsumQ *= Q;
        } else {
          if(isnan2(column[o[lx1-1]-1])) {
            out[j] = NA_REAL;
            continue;
          }
        }
        double wsum = wg[o[0]-1];
        int k = 1; // what about all missing ? -> gives dsort error...

        if(lower) {
          while(wsum < wsumQ) wsum += wg[o[k++]-1];
          if(tiesmean && wsum == wsumQ) {
            double outtmp = column[o[k-1]-1], n = 2;
            while(wg[o[k]-1] == 0) {
              outtmp += column[o[k++]-1];
              ++n;
            }
            out[j] = (outtmp + column[o[k]-1]) / n;
            continue;
          }
        } else {
          while(wsum <= wsumQ) wsum += wg[o[k++]-1];
        }
        out[j] = column[o[k-1]-1];
      }
      // outnth:
        if(drop) {
          Rf_setAttrib(out, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
          return out;
        } else {
          List out(l);
          for(int j = l; j--; ) {
            out[j] = out[j];
            DUPLICATE_ATTRIB(out[j], x[j]);
          }
          DUPLICATE_ATTRIB(out, x);
          Rf_setAttrib(out, R_RowNamesSymbol, Rf_ScalarInteger(1));
          return out;
        }
    } else { // with groups and weights
      if(lx1 != g.size()) stop("length(w) must match length(g)");
      int gi, oi;
      List out(l);
      NumericVector wsumQ(ng);
      if(!narm) {
        for(int i = 0; i != lx1; ++i) wsumQ[g[i]-1] += wg[i];
        wsumQ = wsumQ * Q;
      }
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(lx1 != column.size()) stop("length(w) must match nrow(x)");
        Cdoubleradixsort(ord, TRUE, FALSE, column);
        NumericVector wsum(ng), nthj(ng, NA_REAL);
        if(narm) {
          std::fill(wsumQ.begin(), wsumQ.end(), 0.0);
          for(int i = 0; i != lx1; ++i) if(nisnan(column[i])) wsumQ[g[i]-1] += wg[i]; //  && nisnan(wg[i])
          for(int i = ng; i--; ) {
            if(isnan2(wsumQ[i])) stop("Missing weights in order statistics are currently only supported if x is also missing");
            wsumQ[i] *= Q;
          }
        } else {
          for(int i = 0; i != lx1; ++i) if(isnan2(column[i])) wsum[g[i]-1] = DBL_MAX; // OK ?? // needed??
          for(int i = ng; i--; ) if(isnan2(wsumQ[i]) && wsum[i] != DBL_MAX) stop("Missing weights in order statistics are currently only supported if x is also missing");
        }

        if(tiesmean) {
          std::vector<bool> seen(ng);
          std::vector<int> n(ng, 1); // only needed for 0 weights. Could check above if there are any.
          for(int i = 0; i != lx1; ++i) {
            oi = o[i]-1;
            gi = g[oi]-1;
            if(seen[gi]) continue;
            if(wsum[gi] < wsumQ[gi]) nthj[gi] = column[oi];
            else {
              if(wsum[gi] > wsumQ[gi]) {
                seen[gi] = true;
                continue;
              }
              nthj[gi] += (column[oi]-nthj[gi])/++n[gi];
            }
            wsum[gi] += wg[oi];
          }
        } else if(lower) {
          for(int i = 0; i != lx1; ++i) {
            oi = o[i]-1;
            gi = g[oi]-1;
            if(wsum[gi] < wsumQ[gi]) {
              wsum[gi] += wg[oi];
              nthj[gi] = column[oi];
            }
          }
        } else {
          for(int i = 0; i != lx1; ++i) {
            oi = o[i]-1;
            gi = g[oi]-1;
            if(wsum[gi] <= wsumQ[gi]) {
              wsum[gi] += wg[oi];
              nthj[gi] = column[oi];
            }
          }
        }
        DUPLICATE_ATTRIB(nthj, column);
        out[j] = nthj;
      }
      DUPLICATE_ATTRIB(out, x);
      Rf_setAttrib(out, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -ng));
      return out;
    }
  }
}
