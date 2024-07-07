#include <Rcpp.h>
using namespace Rcpp;

// Notes:
// for mean there are 2 options: "overall.mean" = R_NegInf adds the overall mean. default is centering on 0, or centering on a mean provided, or FALSE = R_PosInf -> no centering, scaling preserves mean
// for sd there is "within.sd" = R_NegInf, scaling by the frequency weighted within-group sd, default is 1, or scaling by a sd provided.
// All other comments are in fvar.cpp (in C++ folder, not on GitHub)

// [[Rcpp::export]]
NumericVector fscaleCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0, const SEXP& w = R_NilValue,
                        bool narm = true, double set_mean = 0, double set_sd = 1) { // could set mean and sd with SEXP, but complicated...
  int l = x.size();
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  NumericVector out = no_init_vector(l);
  //   SHALLOW_DUPLICATE_ATTRIB(out, x); // Any speed loss or overwriting attributes ?
  if (Rf_isNull(w)) { // No weights
    if (ng == 0) {
      if(set_sd == R_NegInf) stop("within.sd can only be calculated when a grouping vector is supplied");
      if(set_mean == R_NegInf) stop("without groups, centering on the overall mean amounts to scaling without centering, so use mean = FALSE instead, or supply a grouping vector to subtract out group means.");
      double n = 0, mean = 0, d1 = 0, M2 = 0;
      if(narm) {
        int j = l-1;
        while(std::isnan(x[j]) && j!=0) --j;
        if(j != 0) {
          for(int i = j+1; i--; ) {
            if(std::isnan(x[i])) continue;
            d1 = x[i]-mean;
            mean += d1 * (1 / ++n);
            M2 += d1*(x[i]-mean);
          }
          M2 = set_sd/sqrt(M2/(n-1)); // good ? -> Yes, works !
        } else { // use goto to make code simpler ?
          std::fill(out.begin(), out.end(), NA_REAL);
          SHALLOW_DUPLICATE_ATTRIB(out, x);
          return out;
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(std::isnan(x[i])) {
            std::fill(out.begin(), out.end(), NA_REAL);
            SHALLOW_DUPLICATE_ATTRIB(out, x);
            return out;
          } else {
            d1 = x[i]-mean;
            mean += d1*(1 / ++n);
            M2 += d1*(x[i]-mean);
          }
        }
        M2 = set_sd/sqrt(M2/(l-1));
      }
      if(std::isnan(M2)) {
        std::fill(out.begin(), out.end(), NA_REAL);
      } else {
        if(set_mean == 0) out = (x-mean)*M2;
        else if(set_mean == R_PosInf) out = (x-mean)*M2 + mean; // best ? // !R_FINITE(set_mean)
        else out = (x-mean)*M2 + set_mean; // best ?
      }
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      double d1 = 0, gl_mean = 0; // Best way of doing this ? How can you declare variables in global scope ?
      // NumericVector mean =  narm ? no_init_vector(ng) : NumericVector(ng); // works but valgrind issue
      // NumericVector M2 = narm ? NumericVector(ng, NA_REAL) : NumericVector(ng);
      // NumericVector n = narm ? NumericVector(ng, 1.0) : NumericVector(ng);
      NumericVector mean(ng), n(ng, (narm) ? 1.0 : 0.0), M2(ng, (narm) ? NA_REAL : 0.0);
      if(narm) {
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
      } else {
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(M2[g[i]-1])) continue;
          if(std::isnan(x[i])) {
            M2[g[i]-1] = NA_REAL;
            ++ngs;
            if(ngs == ng) {
              std::fill(out.begin(), out.end(), NA_REAL);
              SHALLOW_DUPLICATE_ATTRIB(out, x);
              return out;
            }
          } else {
            d1 = x[i]-mean[g[i]-1];
            mean[g[i]-1] += d1 * (1 / ++n[g[i]-1]);
            M2[g[i]-1] += d1*(x[i]-mean[g[i]-1]);
          }
        }
      }
      if(set_sd == R_NegInf) {
        double within_sd = 0;
        int sum_n = 0;
        if(set_mean == R_NegInf) {
          for(int i = ng; i--; ) {
            if(std::isnan(M2[i])) continue;
            within_sd += M2[i];
            M2[i] = 1/sqrt(M2[i]/(n[i]-1));
            gl_mean += mean[i]*n[i];
            sum_n += n[i];
          }
          gl_mean /= sum_n;
        } else {
          for(int i = ng; i--; ) {
            if(std::isnan(M2[i])) continue;
            within_sd += M2[i];
            M2[i] = 1/sqrt(M2[i]/(n[i]-1));
            sum_n += n[i];
          }
          gl_mean = set_mean;
        }
        within_sd = sqrt(within_sd/(sum_n-1));
        M2 = M2 * within_sd; // fastest ?
      } else {
        if(set_mean == R_NegInf) {
          int sum_n = 0;
          for(int i = ng; i--; ) {
            if(std::isnan(M2[i])) continue;
            M2[i] = set_sd/sqrt(M2[i]/(n[i]-1));
            gl_mean += mean[i]*n[i];
            sum_n += n[i];
          }
          gl_mean /= sum_n;
        } else {
          gl_mean = set_mean;
          for(int i = ng; i--; ) if(!std::isnan(M2[i])) M2[i] = set_sd/sqrt(M2[i]/(n[i]-1));
        }
      }
      if(set_mean == 0) {
        for(int i = 0; i != l; ++i) out[i] = (x[i]-mean[g[i]-1])*M2[g[i]-1];
      } else if(set_mean == R_PosInf) {
        for(int i = 0; i != l; ++i) out[i] = (x[i]-mean[g[i]-1])*M2[g[i]-1] + mean[g[i]-1]; // best ?
      } else {
        for(int i = 0; i != l; ++i) out[i] = (x[i]-mean[g[i]-1])*M2[g[i]-1] + gl_mean; // best ?
      }
    }
  } else { // With weights
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match length(x)");
    if (ng == 0) {
      if(set_sd == R_NegInf) stop("within.sd can only be calculated when a grouping vector is supplied");
      if(set_mean == R_NegInf) stop("without groups, centering on the overall mean amounts to scaling without centering, so use mean = FALSE instead, or supply a grouping vector to subtract out group means.");
      double sumw = 0, mean = 0, M2 = 0, d1 = 0;
      if(narm) {
        int j = l-1;
        while((std::isnan(x[j]) || std::isnan(wg[j]) || wg[j] == 0) && j!=0) --j;
        if(j != 0) {
          for(int i = j+1; i--; ) {
            if(std::isnan(x[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
            sumw += wg[i];
            d1 = x[i] - mean;
            mean += d1 * (wg[i] / sumw);
            M2 += wg[i] * d1 * (x[i] - mean);
          }
        } else {
          std::fill(out.begin(), out.end(), NA_REAL);
          SHALLOW_DUPLICATE_ATTRIB(out, x);
          return out;
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(std::isnan(x[i]) || std::isnan(wg[i])) {
            std::fill(out.begin(), out.end(), NA_REAL);
            SHALLOW_DUPLICATE_ATTRIB(out, x);
            return out;
          } else {
            if(wg[i] == 0) continue;
            sumw += wg[i];
            d1 = x[i] - mean;
            mean += d1 * (wg[i] / sumw);
            M2 += wg[i] * d1 * (x[i] - mean);
          }
        }
      }
      M2 = set_sd/sqrt(M2/(sumw-1));
      if(std::isnan(M2)) {
        std::fill(out.begin(), out.end(), NA_REAL);
      } else {
        if(set_mean == 0) out = (x-mean)*M2;
        else if(set_mean == R_PosInf) out = (x-mean)*M2 + mean; // best ?
        else out = (x-mean)*M2 + set_mean; // best ?
      }
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      double d1 = 0, gl_mean = 0; // Best way of doing this ? How can you declare variables in overall scope ?
      // NumericVector M2 = narm ? NumericVector(ng, NA_REAL) : NumericVector(ng);
      NumericVector M2(ng, (narm) ? NA_REAL : 0.0), mean(ng), sumw(ng); // = narm ? no_init_vector(ng) : NumericVector(ng); // works but valgrind issues
      // NumericVector sumw = narm ? no_init_vector(ng) : NumericVector(ng);
      if(narm) {
        for(int i = l; i--; ) {
          if(std::isnan(x[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
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
      } else {
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(M2[g[i]-1])) continue;
          if(std::isnan(x[i]) || std::isnan(wg[i])) {
            M2[g[i]-1] = NA_REAL;
            ++ngs;
            if(ngs == ng) {
              std::fill(out.begin(), out.end(), NA_REAL);
              SHALLOW_DUPLICATE_ATTRIB(out, x);
              return out;
            }
          } else {
            if(wg[i] == 0) continue;
            sumw[g[i]-1] += wg[i];
            d1 = x[i] - mean[g[i]-1];
            mean[g[i]-1] += d1 * (wg[i] / sumw[g[i]-1]);
            M2[g[i]-1] += wg[i] * d1 * (x[i] - mean[g[i]-1]);
          }
        }
      }
      if(set_sd == R_NegInf) {
        double within_sd = 0, sum_sumw = 0;
        if(set_mean == R_NegInf) {
          for(int i = ng; i--; ) {
            if(std::isnan(M2[i])) continue;
            within_sd += M2[i];
            M2[i] = 1/sqrt(M2[i]/(sumw[i]-1));
            gl_mean += mean[i]*sumw[i];
            sum_sumw += sumw[i];
          }
          gl_mean /= sum_sumw;
        } else {
          for(int i = ng; i--; ) {
            if(std::isnan(M2[i])) continue;
            within_sd += M2[i];
            M2[i] = 1/sqrt(M2[i]/(sumw[i]-1));
            sum_sumw += sumw[i];
          }
          gl_mean = set_mean;
        }
        within_sd = sqrt(within_sd/(sum_sumw-1));
        M2 = M2 * within_sd; // fastest ?
      } else {
        if(set_mean == R_NegInf) {
          double sum_sumw = 0;
          for(int i = ng; i--; ) {
            if(std::isnan(M2[i])) continue;
            M2[i] = set_sd/sqrt(M2[i]/(sumw[i]-1));
            gl_mean += mean[i]*sumw[i];
            sum_sumw += sumw[i];
          }
          gl_mean /= sum_sumw;
        } else {
          gl_mean = set_mean;
          for(int i = ng; i--; ) if(!std::isnan(M2[i])) M2[i] = set_sd/sqrt(M2[i]/(sumw[i]-1));
        }
      }
      if(set_mean == 0) {
        for(int i = 0; i != l; ++i) out[i] = (x[i]-mean[g[i]-1])*M2[g[i]-1];
      } else if(set_mean == R_PosInf) {
        for(int i = 0; i != l; ++i) out[i] = (x[i]-mean[g[i]-1])*M2[g[i]-1] + mean[g[i]-1]; // best ?
      } else {
        for(int i = 0; i != l; ++i) out[i] = (x[i]-mean[g[i]-1])*M2[g[i]-1] + gl_mean; // best ?
      }
    }
  }
  SHALLOW_DUPLICATE_ATTRIB(out, x);
  return out;
}


// [[Rcpp::export]]
NumericMatrix fscalemCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0, const SEXP& w = R_NilValue,
                         bool narm = true, double set_mean = 0, double set_sd = 1) {

  int l = x.nrow(), col = x.ncol();
  NumericMatrix out = no_init_matrix(l, col);

  if (Rf_isNull(w)) { // No weights
    if(ng == 0) {
      if(set_sd == R_NegInf) stop("within.sd can only be calculated when a grouping vector is supplied");
      if(set_mean == R_NegInf) stop("without groups, centering on the overall mean amounts to scaling without centering, so use mean = FALSE instead, or supply a grouping vector to subtract out group means.");
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column outj = out( _ , j);
        double nj = 0, meanj = 0, d1 = 0, M2j = 0;
        if(narm) { // faster using 2 loops over columns ?
          int k = l-1;
          while(std::isnan(column[k]) && k!=0) --k;
          if(k != 0) {
            for(int i = k+1; i--; ) {
              if(std::isnan(column[i])) continue;
              d1 = column[i]-meanj;
              meanj += d1 * (1 / ++nj);
              M2j += d1*(column[i]-meanj);
            }
            M2j = set_sd/sqrt(M2j/(nj-1));
          } else {
            std::fill(outj.begin(), outj.end(), NA_REAL);
            continue; // Necessary
          }
        } else {
          for(int i = 0; i != l; ++i) {
            if(std::isnan(column[i])) {
              M2j = NA_REAL;
              break;
            } else {
              d1 = column[i]-meanj;
              meanj += d1 * (1 / ++nj);
              M2j += d1*(column[i]-meanj);
            }
          }
          M2j = set_sd/sqrt(M2j/(l-1));
        }
        if(std::isnan(M2j)) {
          std::fill(outj.begin(), outj.end(), NA_REAL);
        } else {
          if(set_mean == 0) outj = (column-meanj)*M2j;
          else if(set_mean == R_PosInf) outj = (column-meanj)*M2j + meanj; // best ?
          else outj = (column-meanj)*M2j + set_mean; // best ?
        }
      }
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      // Better way ?
      NumericVector meanj(ng), nj(ng), M2j(ng);
      // NumericVector meanj = no_init_vector(ng), nj = no_init_vector(ng), M2j = no_init_vector(ng); // Works but valgrind issue
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column outj = out( _ , j);
        double d1 = 0, gl_meanj = 0;
        if(narm) { // better do two loops ??
          std::fill(M2j.begin(), M2j.end(), NA_REAL);
          for(int i = l; i--; ) {
            if(std::isnan(column[i])) continue;
            if(std::isnan(M2j[g[i]-1])) {
              meanj[g[i]-1] = column[i];
              M2j[g[i]-1] = 0;
              nj[g[i]-1] = 1;
            } else {
              d1 = column[i]-meanj[g[i]-1];
              meanj[g[i]-1] += d1 * (1 / ++nj[g[i]-1]);
              M2j[g[i]-1] += d1*(column[i]-meanj[g[i]-1]);
            }
          }
        } else {
          for(int i = ng; i--; ) meanj[i] = M2j[i] = nj[i] = 0;
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
              d1 = column[i]-meanj[g[i]-1];
              meanj[g[i]-1] += d1 * (1 / ++nj[g[i]-1]);
              M2j[g[i]-1] += d1*(column[i]-meanj[g[i]-1]);
            }
          }
        }
        if(set_sd == R_NegInf) { // best way of coding ? Goes through all the if conditions for every column...
          double within_sdj = 0;
          int sum_nj = 0;
          if(set_mean == R_NegInf) {
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              within_sdj += M2j[i];
              M2j[i] = 1/sqrt(M2j[i]/(nj[i]-1));
              gl_meanj += meanj[i]*nj[i];
              sum_nj += nj[i];
            }
            gl_meanj /= sum_nj;
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              within_sdj += M2j[i];
              M2j[i] = 1/sqrt(M2j[i]/(nj[i]-1));
              sum_nj += nj[i];
            }
            gl_meanj = set_mean;
          }
          within_sdj = sqrt(within_sdj/(sum_nj-1));
          M2j = M2j * within_sdj; // fastest ?
        } else {
          if(set_mean == R_NegInf) {
            int sum_nj = 0;
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              M2j[i] = set_sd/sqrt(M2j[i]/(nj[i]-1));
              gl_meanj += meanj[i]*nj[i];
              sum_nj += nj[i];
            }
            gl_meanj /= sum_nj;
          } else {
            gl_meanj = set_mean;
            for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = set_sd/sqrt(M2j[i]/(nj[i]-1));
          }
        }
        if(set_mean == 0) {
          for(int i = 0; i != l; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
        } else if(set_mean == R_PosInf) {
          for(int i = 0; i != l; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1] + meanj[g[i]-1]; // best ?
        } else {
          for(int i = 0; i != l; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1] + gl_meanj; // best ?
        }
        loopend:;
      }
    }
  } else { // With weights
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match nrow(X)");
    if(ng == 0) {
      if(set_sd == R_NegInf) stop("within.sd can only be calculated when a grouping vector is supplied");
      if(set_mean == R_NegInf) stop("without groups, centering on the overall mean amounts to scaling without centering, so use mean = FALSE instead, or supply a grouping vector to subtract out group means.");
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column outj = out( _ , j);
        double sumwj = 0, meanj = 0, M2j = 0, d1 = 0;
        if(narm) {
          int k = l-1;
          while((std::isnan(column[k]) || std::isnan(wg[k]) || wg[k] == 0) && k!=0) --k;
          if(k != 0) {
            for(int i = k+1; i--; ) {
              if(std::isnan(column[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
              sumwj += wg[i];
              d1 = column[i] - meanj;
              meanj += d1 * (wg[i] / sumwj);
              M2j += wg[i] * d1 * (column[i] - meanj);
            }
          } else {
            std::fill(outj.begin(), outj.end(), NA_REAL);
            continue; // Necessary
          }
        } else {
          for(int i = 0; i != l; ++i) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) {
              M2j = NA_REAL;
              break;
            } else {
              if(wg[i] == 0) continue;
              sumwj += wg[i];
              d1 = column[i] - meanj;
              meanj += d1 * (wg[i] / sumwj);
              M2j += wg[i] * d1 * (column[i] - meanj);
            }
          }
        }
        M2j = set_sd/sqrt(M2j/(sumwj-1));
        if(std::isnan(M2j)) {
          std::fill(outj.begin(), outj.end(), NA_REAL);
        } else {
          if(set_mean == 0) outj = (column-meanj)*M2j;
          else if(set_mean == R_PosInf) outj = (column-meanj)*M2j + meanj; // best ?
          else outj = (column-meanj)*M2j + set_mean; // best ?
        }
      }
    } else { // with groups and weights
      if(g.size() != l) stop("length(g) must match nrow(X)");
      // Works but valgrind issue
      // NumericVector meanj = no_init_vector(ng), sumwj = no_init_vector(ng), M2j = no_init_vector(ng);
      NumericVector meanj(ng), sumwj(ng), M2j(ng); // better for valgrind
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column outj = out( _ , j);
        double d1 = 0, gl_meanj = 0;
        if(narm) {
          std::fill(M2j.begin(), M2j.end(), NA_REAL);
          for(int i = l; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
            if(std::isnan(M2j[g[i]-1])) {
              sumwj[g[i]-1] = wg[i];
              meanj[g[i]-1] = column[i];
              M2j[g[i]-1] = 0;
            } else {
              sumwj[g[i]-1] += wg[i];
              d1 = column[i] - meanj[g[i]-1];
              meanj[g[i]-1] += d1 * (wg[i] / sumwj[g[i]-1]);
              M2j[g[i]-1] += wg[i] * d1 * (column[i] - meanj[g[i]-1]);
            }
          }
        } else {
          for(int i = ng; i--; ) meanj[i] = M2j[i] = sumwj[i] = 0;
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
              if(wg[i] == 0) continue;
              sumwj[g[i]-1] += wg[i];
              d1 = column[i] - meanj[g[i]-1];
              meanj[g[i]-1] += d1 * (wg[i] / sumwj[g[i]-1]);
              M2j[g[i]-1] += wg[i] * d1 * (column[i] - meanj[g[i]-1]);
            }
          }
        }
        if(set_sd == R_NegInf) { // best way of coding ? Goes through all the if conditions for every column...
          double within_sdj = 0, sum_sumwj = 0;
          if(set_mean == R_NegInf) {
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              within_sdj += M2j[i];
              M2j[i] = 1/sqrt(M2j[i]/(sumwj[i]-1));
              gl_meanj += meanj[i]*sumwj[i];
              sum_sumwj += sumwj[i];
            }
            gl_meanj /= sum_sumwj;
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              within_sdj += M2j[i];
              M2j[i] = 1/sqrt(M2j[i]/(sumwj[i]-1));
              sum_sumwj += sumwj[i];
            }
            gl_meanj = set_mean;
          }
          within_sdj = sqrt(within_sdj/(sum_sumwj-1));
          M2j = M2j * within_sdj; // fastest ?
        } else {
          if(set_mean == R_NegInf) {
            double sum_sumwj = 0;
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              M2j[i] = set_sd/sqrt(M2j[i]/(sumwj[i]-1));
              gl_meanj += meanj[i]*sumwj[i];
              sum_sumwj += sumwj[i];
            }
            gl_meanj /= sum_sumwj;
          } else {
            gl_meanj = set_mean;
            for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = set_sd/sqrt(M2j[i]/(sumwj[i]-1));
          }
        }
        if(set_mean == 0) {
          for(int i = 0; i != l; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
        } else if(set_mean == R_PosInf) {
          for(int i = 0; i != l; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1] + meanj[g[i]-1]; // best ?
        } else {
          for(int i = 0; i != l; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1] + gl_meanj; // best ?
        }
        loopend2:;
      }
    }
  }
  SHALLOW_DUPLICATE_ATTRIB(out, x);
  return out;
}

// [[Rcpp::export]]
List fscalelCpp(const List& x, int ng = 0, const IntegerVector& g = 0, const SEXP& w = R_NilValue,
                bool narm = true, double set_mean = 0, double set_sd = 1) {

  int l = x.size();
  List out(l);

  if (Rf_isNull(w)) { // No weights
    if(ng == 0) {
      if(set_sd == R_NegInf) stop("within.sd can only be calculated when a grouping vector is supplied");
      if(set_mean == R_NegInf) stop("without groups, centering on the overall mean amounts to scaling without centering, so use mean = FALSE instead, or supply a grouping vector to subtract out group means.");
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        int row = column.size();
        NumericVector outj = no_init_vector(row);
        double nj = 0, meanj = 0, d1 = 0, M2j = 0;
        if(narm) {
          int k = row-1;
          while(std::isnan(column[k]) && k!=0) --k;
          if(k != 0) {
            for(int i = k+1; i--; ) {
              if(std::isnan(column[i])) continue;
              d1 = column[i]-meanj;
              meanj += d1 * (1 / ++nj);
              M2j += d1*(column[i]-meanj);
            }
            M2j = set_sd/sqrt(M2j/(nj-1));
          } else {
            std::fill(outj.begin(), outj.end(), NA_REAL); // outj = rep(NA_REAL, row); // fastest option ! (faster than std::fill)
            goto loopend; // Necessary
          }
        } else {
          for(int i = 0; i != row; ++i) {
            if(std::isnan(column[i])) {
              M2j = NA_REAL;
              break;
            } else {
              d1 = column[i]-meanj;
              meanj += d1 * (1 / ++nj);
              M2j += d1*(column[i]-meanj);
            }
          }
          M2j = set_sd/sqrt(M2j/(row-1));
        }
        if(std::isnan(M2j)) {
          std::fill(outj.begin(), outj.end(), NA_REAL);
        } else {
          if(set_mean == 0) outj = (column-meanj)*M2j;
          else if(set_mean == R_PosInf) outj = (column-meanj)*M2j + meanj; // best ?
          else outj = (column-meanj)*M2j + set_mean; // best ?
        }
        loopend:;
        SHALLOW_DUPLICATE_ATTRIB(outj, column);
        out[j] = outj;
      }
    } else { // with groups
      int gss = g.size();
      // Better way ?
      NumericVector meanj(ng), nj(ng), M2j(ng);
      // NumericVector meanj = no_init_vector(ng), nj = no_init_vector(ng), M2j = no_init_vector(ng); // Works but valgrind issue
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector outj = no_init_vector(gss);
        double d1 = 0, gl_meanj = 0;
        if(narm) { // better do two loops ?
          std::fill(M2j.begin(), M2j.end(), NA_REAL);
          for(int i = gss; i--; ) {
            if(std::isnan(column[i])) continue;
            if(std::isnan(M2j[g[i]-1])) {
              meanj[g[i]-1] = column[i];
              M2j[g[i]-1] = 0;
              nj[g[i]-1] = 1;
            } else {
              d1 = column[i]-meanj[g[i]-1];
              meanj[g[i]-1] += d1 * (1 / ++nj[g[i]-1]);
              M2j[g[i]-1] += d1*(column[i]-meanj[g[i]-1]);
            }
          }
        } else {
          for(int i = ng; i--; ) meanj[i] = M2j[i] = nj[i] = 0;
          int ngs = 0;
          for(int i = 0; i != gss; ++i) {
            if(std::isnan(M2j[g[i]-1])) continue;
            if(std::isnan(column[i])) {
              M2j[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) {
                std::fill(outj.begin(), outj.end(), NA_REAL);
                goto loopend2;
              }
            } else {
              d1 = column[i]-meanj[g[i]-1];
              meanj[g[i]-1] += d1 * (1 / ++nj[g[i]-1]);
              M2j[g[i]-1] += d1*(column[i]-meanj[g[i]-1]);
            }
          }
        }
        if(set_sd == R_NegInf) { // best way of coding ? Goes through all the if conditions for every column...
          double within_sdj = 0;
          int sum_nj = 0;
          if(set_mean == R_NegInf) {
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              within_sdj += M2j[i];
              M2j[i] = 1/sqrt(M2j[i]/(nj[i]-1));
              gl_meanj += meanj[i]*nj[i];
              sum_nj += nj[i];
            }
            gl_meanj /= sum_nj;
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              within_sdj += M2j[i];
              M2j[i] = 1/sqrt(M2j[i]/(nj[i]-1));
              sum_nj += nj[i];
            }
            gl_meanj = set_mean;
          }
          within_sdj = sqrt(within_sdj/(sum_nj-1));
          M2j = M2j * within_sdj; // fastest ?
        } else {
          if(set_mean == R_NegInf) {
            int sum_nj = 0;
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              M2j[i] = set_sd/sqrt(M2j[i]/(nj[i]-1));
              gl_meanj += meanj[i]*nj[i];
              sum_nj += nj[i];
            }
            gl_meanj /= sum_nj;
          } else {
            gl_meanj = set_mean;
            for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = set_sd/sqrt(M2j[i]/(nj[i]-1));
          }
        }
        if(set_mean == 0) {
          for(int i = 0; i != gss; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
        } else if(set_mean == R_PosInf) {
          for(int i = 0; i != gss; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1] + meanj[g[i]-1]; // best ?
        } else {
          for(int i = 0; i != gss; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1] + gl_meanj; // best ?
        }
        loopend2:;
        SHALLOW_DUPLICATE_ATTRIB(outj, column);
        out[j] = outj;
      }
    }
  } else { // With weights
    NumericVector wg = w;
    int wgs = wg.size();
    if(ng == 0) {
      if(set_sd == R_NegInf) stop("within.sd can only be calculated when a grouping vector is supplied");
      if(set_mean == R_NegInf) stop("without groups, centering on the overall mean amounts to scaling without centering, so use mean = FALSE instead, or supply a grouping vector to subtract out group means.");
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(wgs != column.size()) stop("length(w) must match nrow(X)");
        NumericVector outj = no_init_vector(wgs);
        double sumwj = 0, meanj = 0, M2j = 0, d1 = 0;
        if(narm) {
          int k = wgs-1;
          while((std::isnan(column[k]) || std::isnan(wg[k]) || wg[k] == 0) && k!=0) --k;
          if(k != 0) {
            for(int i = k+1; i--; ) {
              if(std::isnan(column[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
              sumwj += wg[i];
              d1 = column[i] - meanj;
              meanj += d1 * (wg[i] / sumwj);
              M2j += wg[i] * d1 * (column[i] - meanj);
            }
          } else {
            std::fill(outj.begin(), outj.end(), NA_REAL);
            goto loopend3; // Necessary
          }
        } else {
          for(int i = 0; i != wgs; ++i) {
            if(std::isnan(column[i]) || std::isnan(wg[i])) {
              M2j = NA_REAL;
              break;
            } else {
              if(wg[i] == 0) continue;
              sumwj += wg[i];
              d1 = column[i] - meanj;
              meanj += d1 * (wg[i] / sumwj);
              M2j += wg[i] * d1 * (column[i] - meanj);
            }
          }
        }
        M2j = set_sd/sqrt(M2j/(sumwj-1));
        if(std::isnan(M2j)) {
          std::fill(outj.begin(), outj.end(), NA_REAL);
        } else {
          if(set_mean == 0) outj = (column-meanj)*M2j;
          else if(set_mean == R_PosInf) outj = (column-meanj)*M2j + meanj; // best ?
          else outj = (column-meanj)*M2j + set_mean; // best ?
        }
        loopend3:;
        SHALLOW_DUPLICATE_ATTRIB(outj, column);
        out[j] = outj;
      }
    } else { // with groups and weights
      int gss = g.size();
      if(gss != wgs) stop("length(w) must match length(g)");
      NumericVector meanj(ng), sumwj(ng), M2j(ng);
      // NumericVector meanj = no_init_vector(ng), sumwj = no_init_vector(ng), M2j = no_init_vector(ng); // Works but valgrind issue
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("length(g) must match nrow(X)");
        NumericVector outj = no_init_vector(gss);
        double d1 = 0, gl_meanj = 0;
        if(narm) {
          std::fill(M2j.begin(), M2j.end(), NA_REAL);
          for(int i = gss; i--; ) {
            if(std::isnan(column[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
            if(std::isnan(M2j[g[i]-1])) {
              sumwj[g[i]-1] = wg[i];
              meanj[g[i]-1] = column[i];
              M2j[g[i]-1] = 0;
            } else {
              sumwj[g[i]-1] += wg[i];
              d1 = column[i] - meanj[g[i]-1];
              meanj[g[i]-1] += d1 * (wg[i] / sumwj[g[i]-1]);
              M2j[g[i]-1] += wg[i] * d1 * (column[i] - meanj[g[i]-1]);
            }
          }
        } else {
          for(int i = ng; i--; ) meanj[i] = M2j[i] = sumwj[i] = 0;
          int ngs = 0;
          for(int i = 0; i != gss; ++i) {
            if(std::isnan(M2j[g[i]-1])) continue;
            if(std::isnan(column[i]) || std::isnan(wg[i])) {
              M2j[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) {
                std::fill(outj.begin(), outj.end(), NA_REAL);
                goto loopend4;
              }
            } else {
              if(wg[i] == 0) continue;
              sumwj[g[i]-1] += wg[i];
              d1 = column[i] - meanj[g[i]-1];
              meanj[g[i]-1] += d1 * (wg[i] / sumwj[g[i]-1]);
              M2j[g[i]-1] += wg[i] * d1 * (column[i] - meanj[g[i]-1]);
            }
          }
        }
        if(set_sd == R_NegInf) { // best way of coding ? Goes through all the if conditions for every column...
          double within_sdj = 0, sum_sumwj = 0;
          if(set_mean == R_NegInf) {
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              within_sdj += M2j[i];
              M2j[i] = 1/sqrt(M2j[i]/(sumwj[i]-1));
              gl_meanj += meanj[i]*sumwj[i];
              sum_sumwj += sumwj[i];
            }
            gl_meanj /= sum_sumwj;
          } else {
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              within_sdj += M2j[i];
              M2j[i] = 1/sqrt(M2j[i]/(sumwj[i]-1));
              sum_sumwj += sumwj[i];
            }
            gl_meanj = set_mean;
          }
          within_sdj = sqrt(within_sdj/(sum_sumwj-1));
          M2j = M2j * within_sdj; // fastest ?
        } else {
          if(set_mean == R_NegInf) {
            double sum_sumwj = 0;
            for(int i = ng; i--; ) {
              if(std::isnan(M2j[i])) continue;
              M2j[i] = set_sd/sqrt(M2j[i]/(sumwj[i]-1));
              gl_meanj += meanj[i]*sumwj[i];
              sum_sumwj += sumwj[i];
            }
            gl_meanj /= sum_sumwj;
          } else {
            gl_meanj = set_mean;
            for(int i = ng; i--; ) if(!std::isnan(M2j[i])) M2j[i] = set_sd/sqrt(M2j[i]/(sumwj[i]-1));
          }
        }
        if(set_mean == 0) {
          for(int i = 0; i != gss; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1];
        } else if(set_mean == R_PosInf) {
          for(int i = 0; i != gss; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1] + meanj[g[i]-1]; // best ?
        } else {
          for(int i = 0; i != gss; ++i) outj[i] = (column[i]-meanj[g[i]-1])*M2j[g[i]-1] + gl_meanj; // best ?
        }
        loopend4:;
        SHALLOW_DUPLICATE_ATTRIB(outj, column);
        out[j] = outj;
      }
    }
  }
  SHALLOW_DUPLICATE_ATTRIB(out, x);
  return out;
}



