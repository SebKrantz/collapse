#include <Rcpp/Lighter>
using namespace Rcpp;

// Return Options:
// ret = 1 - differences
// ret = 2 - log differences
// ret = 3 - log-difference growth rates
// ret = 4 - exact growth rates
// Also: if rho != 1, quasi-differencing and log differencing with rho... i.e. for Cochrane-Orcutt regression

// This Approach: currently does not support iterated differences on irregular time-series and panel data !
// TODO: Make comprehensive...

// Note: Now taking logs in R -> Faster and smaller compiled code !
// ... some systems get this wrong, possibly depends on what libs are loaded //
// static inline double R_log(double x) {
//   return x > 0 ? log(x) : x == 0 ? R_NegInf : R_NaN;
// }

template <typename F>
NumericVector fdiffgrowthCppImpl(const NumericVector& x, const IntegerVector& n = 1, const IntegerVector& diff = 1,
                                 double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0,
                                 const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue,
                                 std::string stub = "", bool names = true, F FUN = [](double y, double x) { return y-x; }) {

  int l = x.size(), ns = n.size(), ds = diff.size(), zeros = 0, pos = INT_MAX;
  IntegerVector absn = no_init_vector(ns);
  for(int i = ns; i--; ) {
    if(n[i] == pos) stop("duplicated values in n detected"); // because one might mistakenly pass a factor to the n-slot
    pos = n[i];
    if(pos == 0) ++zeros;
    if(pos < 0) {
      if(pos == NA_INTEGER) stop("NA in n");
      absn[i] = -pos;
    } else absn[i] = pos;
  }
  pos = 0;
  std::string stub2 = names ? "F" + stub : "";

  int ncol = (ns-zeros)*ds+zeros;
  if(ncol == 1) names = false;
  NumericMatrix out = no_init_matrix(l, ncol);
  CharacterVector colnam = names ? no_init_vector(ncol) : no_init_vector(1);
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING;
  CharacterVector diffc = names ? Rf_coerceVector(diff, STRSXP) : NA_STRING;
  if(ng == 0) { // No groups
    if(Rf_isNull(t)) { // Ordered data
      for(int p = 0; p != ns; ++p) {
        int np = n[p];
        if(np>0) { // Positive lagged and iterated differences
          int d1 = diff[0], end = np*d1;
          bool L1 = np == 1;
          if(d1 < 1) stop("diff must be a vector of integers > 0");
          if(end >= l) stop("n * diff needs to be < length(x)");
          NumericMatrix::Column outp = out( _ , pos);
          if(names) {
            if(L1) colnam[pos] = stub + diffc[0];
            else colnam[pos] = "L" + nc[p] + stub + diffc[0];
          }
          ++pos;
          for(int i = np; i != l; ++i) outp[i] = FUN(x[i], x[i - np]);
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1)-1;
            for(int i = l-1; i != start; --i) outp[i] = FUN(outp[i], outp[i - np]);
          }
          for(int i = end; i--; ) outp[i] = fill;
          if(ds > 1) {
            NumericVector outtemp = outp;
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1], end = np*dq;
              if(end >= l) stop("n * diff needs to be < length(x)");
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int start = np*(k+1)-1;
                for(int i = l-1; i != start; --i) outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
              }
              for(int i = np*L_dq; i != end; ++i) outtemp[i] = fill;
              out( _ , pos) = outtemp;
              if(names) {
                if(L1) colnam[pos] = stub + diffc[q];
                else colnam[pos] = "L" + nc[p] + stub + diffc[q];
              }
              ++pos;
            }
          }
        } else if(np<0) { // (Negative) leaded and iterated differences
          int d1 = diff[0], end = l+np*d1;
          bool F1 = np == -1;
          if(d1 < 1) stop("diff must be a vector of integers > 0");
          if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
          NumericMatrix::Column outp = out( _ , pos);
          if(names) {
            if(F1) colnam[pos] = stub2 + diffc[0];
            else colnam[pos] = "F" + nc[p] + stub + diffc[0];
          }
          ++pos;
          for(int i = l+np; i--; ) outp[i] = FUN(x[i], x[i - np]);
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int final = l+np*(k+1);
            for(int i = 0; i != final; ++i) outp[i] = FUN(outp[i], outp[i - np]);
          }
          for(int i = end; i != l; ++i) outp[i] = fill;
          if(ds > 1) {
            NumericVector outtemp = outp;
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
              if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int final = l+np*(k+1);
                for(int i = 0; i != final; ++i) outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
              }
              for(int i = end; i != start; ++i) outtemp[i] = fill;
              out( _ , pos) = outtemp;
              if(names) {
                if(F1) colnam[pos] = stub2 + diffc[q];
                else colnam[pos] = "F"+ nc[p] + stub + diffc[q];
              }
              ++pos;
            }
          }
        } else {
          out( _ , pos) = x;
          if(names) colnam[pos] = "--";
          ++pos;
        }
      }
    } else { // Unordered data: Timevar provided
      IntegerVector ord = t;
      if(l != ord.size()) stop("length(x) must match length(t)");
      int min = INT_MAX, max = INT_MIN, osize, temp;
      for(int i = 0; i != l; ++i) {
        if(ord[i] < min) min = ord[i];
        if(ord[i] > max) max = ord[i];
      }
      if(min == NA_INTEGER) stop("Timevar contains missing values");
      osize = max-min+1;
      bool regular = osize == l;
      IntegerVector omap(osize), ord2 = regular ? no_init_vector(1) : no_init_vector(l);
      if(!regular) { // Irregular time series
        if(osize > 10000000 && osize > 3 * l) warning("Your time series is very irregular. Need to create an internal ordering vector of length %s to represent it.", osize);
        if(Rcpp::max(diff) > 1) stop("Iterations are currently only supported for regular time series. See ?seqid to identify the regular sequences in your time series, or just apply this function multiple times.");
        for(int i = 0; i != l; ++i) {
          temp = ord[i] - min; // Best ? Or direct assign to ord2[i] ? Also check for panel version..
          if(omap[temp]) stop("Repeated values in timevar");
          omap[temp] = i+1;
          ord2[i] = temp;
        }
      } else { // Regular time series
        for(int i = 0; i != l; ++i) {
          temp = ord[i] - min;
          if(omap[temp]) stop("Repeated values in timevar");
          omap[temp] = i;
        }
      }
      for(int p = 0; p != ns; ++p) {
        int np = n[p];
        if(np>0) { // Positive lagged and iterated differences
          int d1 = diff[0], end = np*d1;
          bool L1 = np == 1;
          if(d1 < 1) stop("diff must be a vector of integers > 0");
          if(end >= l) stop("n * diff needs to be < length(x)");
          NumericMatrix::Column outp = out( _ , pos);
          if(names) {
            if(L1) colnam[pos] = stub + diffc[0];
            else colnam[pos] = "L" + nc[p] + stub + diffc[0];
          }
          ++pos;
          if(regular) {
            for(int i = np; i != l; ++i) outp[omap[i]] = FUN(x[omap[i]], x[omap[i - np]]);
          } else {
            for(int i = 0; i != l; ++i) { // Smarter solution using while ???
              if(ord2[i] >= np && (temp = omap[ord2[i] - np])) {
                outp[i] = FUN(x[i], x[temp-1]);
              } else {
                outp[i] = fill;
              }
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1)-1;
            for(int i = l-1; i != start; --i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
          }
          if(regular) for(int i = end; i--; ) outp[omap[i]] = fill;
          if(ds > 1) {
            NumericVector outtemp = outp;
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1], end = np*dq;
              if(end >= l) stop("n * diff needs to be < length(x)");
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int start = np*(k+1)-1;
                for(int i = l-1; i != start; --i) outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
              }
              for(int i = np*L_dq; i != end; ++i) outtemp[omap[i]] = fill;
              out( _ , pos) = outtemp;
              if(names) {
                if(L1) colnam[pos] = stub + diffc[q];
                else colnam[pos] = "L" + nc[p] + stub + diffc[q];
              }
              ++pos;
            }
          }
        } else if(np<0) { // (Negative) leaded and iterated differences
          int d1 = diff[0], end = l+np*d1;
          bool F1 = np == -1;
          if(d1 < 1) stop("diff must be a vector of integers > 0");
          if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
          NumericMatrix::Column outp = out( _ , pos);
          if(names) {
            if(F1) colnam[pos] = stub2 + diffc[0];
            else colnam[pos] = "F" + nc[p] + stub + diffc[0];
          }
          ++pos;
          if(regular) {
            for(int i = l+np; i--; ) outp[omap[i]] = FUN(x[omap[i]], x[omap[i - np]]);
          } else {
            for(int i = 0, osnp = osize + np; i != l; ++i) { // Smarter solution using while ???
              if(ord2[i] < osnp && (temp = omap[ord2[i] - np])) {
                outp[i] = FUN(x[i], x[temp-1]);
              } else {
                outp[i] = fill;
              }
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int final = l+np*(k+1);
            for(int i = 0; i != final; ++i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
          }
          if(regular) for(int i = end; i != l; ++i) outp[omap[i]] = fill;
          if(ds > 1) {
            NumericVector outtemp = outp;
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
              if(end <= 0) stop("abs(n * diff) needs to be < length(x)");
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int final = l+np*(k+1);
                for(int i = 0; i != final; ++i) outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
              }
              for(int i = end; i != start; ++i) outtemp[omap[i]] = fill;
              out( _ , pos) = outtemp;
              if(names) {
                if(F1) colnam[pos] = stub2 + diffc[q];
                else colnam[pos] = "F"+ nc[p] + stub + diffc[q];
              }
              ++pos;
            }
          }
        } else {
          out( _ , pos) = x;
          if(names) colnam[pos] = "--";
          ++pos;
        }
      }
    }
  } else {
    if(l != g.size()) stop("length(x) must match length(g)");
    int ags = l/ng, ngp = ng+1, maxdiff = max(diff);
    if(Rf_isNull(t)) {
      bool cond = !Rf_isNull(gs);
      IntegerVector gsv = (cond || maxdiff == 1) ? no_init_vector(1) : IntegerVector(ng);
      int *pgsv = cond ? INTEGER(gs)-1 : INTEGER(gsv)-1;
      if(maxdiff != 1) {
        if(cond) {
          if(ng != Rf_length(gs)) stop("ng must match length(gs)");
        } else {
          for(int i = 0; i != l; ++i) ++pgsv[g[i]];
        }
      }
      // int seen[ngp], memsize = sizeof(int)*(ngp);
      for(int p = 0; p != ns; ++p) {
        int np = n[p];
        if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
        if(np>0) { // Positive lagged and iterated differences
          int d1 = diff[0];
          bool L1 = np == 1;
          if(d1 < 1) stop("diff must be a vector of integers > 0");
          NumericMatrix::Column outp = out( _ , pos);
          std::vector<int> seen(ngp); // memset(seen, 0, memsize);
          if(names) {
            if(L1) colnam[pos] = stub + diffc[0];
            else colnam[pos] = "L" + nc[p] + stub + diffc[0];
          }
          ++pos;
          for(int i = 0; i != l; ++i) {
            if(seen[g[i]] == np) outp[i] = FUN(x[i], x[i - np]);
            else {
              outp[i] = fill;
              ++seen[g[i]];
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1);
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            for(int i = l; i--; ) {
              if(seen[g[i]] == pgsv[g[i]]-start) outp[i] = fill;
              else {
                outp[i] = FUN(outp[i], outp[i - np]);
                ++seen[g[i]];
              }
            }
          }
          if(ds > 1) {
            NumericVector outtemp = outp;
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1];
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int start = np*(k+1); // Right ? -> seems so
                std::vector<int> seen(ngp); // memset(seen, 0, memsize); // Needed, because it loops from the beginning
                for(int i = l; i--; ) {
                  if(seen[g[i]] == pgsv[g[i]]-start) outtemp[i] = fill;
                  else {
                    outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                    ++seen[g[i]];
                  }
                }
              }
              out( _ , pos) = outtemp;
              if(names) {
                if(L1) colnam[pos] = stub + diffc[q];
                else colnam[pos] = "L" + nc[p] + stub + diffc[q];
              }
              ++pos;
            }
          }
        } else if(np<0) { // (Negative) leaded and iterated differences
          int d1 = diff[0];
          bool F1 = np == -1;
          if(d1 < 1) stop("diff must be a vector of integers > 0");
          NumericMatrix::Column outp = out( _ , pos);
          std::vector<int> seen(ngp); // memset(seen, 0, memsize);
          if(names) {
            if(F1) colnam[pos] = stub2 + diffc[0];
            else colnam[pos] = "F" + nc[p] + stub + diffc[0];
          }
          ++pos;
          for(int i = l; i--; ) {
            if(seen[g[i]] == np) outp[i] = FUN(x[i], x[i - np]);
            else {
              outp[i] = fill;
              --seen[g[i]];
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1);
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            for(int i = 0; i != l; ++i) {
              if(seen[g[i]] == pgsv[g[i]]+start) outp[i] = fill;
              else {
                outp[i] = FUN(outp[i], outp[i - np]);
                ++seen[g[i]];
              }
            }
          }
          if(ds > 1) {
            NumericVector outtemp = outp;
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1];
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int start = np*(k+1);
                std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                for(int i = 0; i != l; ++i) {
                  if(seen[g[i]] == pgsv[g[i]]+start) outtemp[i] = fill;
                  else {
                    outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                    ++seen[g[i]];
                  }
                }
              }
              out( _ , pos) = outtemp;
              if(names) {
                if(F1) colnam[pos] = stub2 + diffc[q];
                else colnam[pos] = "F"+ nc[p] + stub + diffc[q];
              }
              ++pos;
            }
          }
        } else {
          out( _ , pos) = x;
          if(names) colnam[pos] = "--";
          ++pos;
        }
      }
    } else { // Unordered data: Timevar Provided
      IntegerVector ord = t;
      int temp;
      if(l != ord.size()) stop("length(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX), max(ngp, INT_MIN), cgs = no_init_vector(ngp);
      for(int i = 0; i != l; ++i) {
        temp = g[i];
        if(ord[i] < min[temp]) min[temp] = ord[i];
        if(ord[i] > max[temp]) max[temp] = ord[i];
      }
      temp = 0;
      for(int i = 1; i != ngp; ++i) {
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values");
        if(min[i] == INT_MAX) continue; // Needed in case of unused factor levels (group vector too large)
        cgs[i] = temp; // This needs to b here (for unused factor levels case...)
        max[i] -= min[i] - 1; // need max[i] which stores the complete group sizes only if p<0 e.g. if computing leads..
        temp += max[i];
      }
      // omap provides the ordering to order the vector (needed to find previous / next values)
      bool regular = temp == l;
      IntegerVector omap(temp), ord2 = no_init_vector(l);
      if(!regular) { // Irregular panel
        if(temp > 10000000 && temp > 3 * l) warning("Your panel is very irregular. Need to create an internal ordering vector of length %s to represent it.", temp);
        if(maxdiff > 1) stop("Iterations are currently only supported for regular panels. See ?seqid to identify the regular sequences in your panel, or just apply this function multiple times.");
        for(int i = 0; i != l; ++i) {
          ord2[i] = ord[i] - min[g[i]];
          temp = cgs[g[i]] + ord2[i];
          if(omap[temp]) stop("Repeated values of timevar within one or more groups");
          omap[temp] = i+1; // needed to add 1 to distinguish between 0 and gap
        }
      } else { // Regular panel
        for(int i = 0; i != l; ++i) {
          ord2[i] = ord[i] - min[g[i]];
          temp = cgs[g[i]] + ord2[i];
          if(omap[temp]) stop("Repeated values of timevar within one or more groups");
          omap[temp] = i;
        }
      }

      for(int p = 0; p != ns; ++p) {
        int np = n[p];
        if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
        if(np>0) { // Positive lagged and iterated differences
          int d1 = diff[0];
          bool L1 = np == 1;
          if(d1 < 1) stop("diff must be a vector of integers > 0");
          NumericMatrix::Column outp = out( _ , pos);
          if(names) {
            if(L1) colnam[pos] = stub + diffc[0];
            else colnam[pos] = "L" + nc[p] + stub + diffc[0];
          }
          ++pos;
          if(regular) {
            for(int i = 0; i != l; ++i) {
              if(ord2[i] >= np) {
                outp[i] = FUN(x[i], x[omap[cgs[g[i]]+ord2[i]-np]]);
              } else {
                outp[i] = fill;
              }
            }
          } else {
            for(int i = 0; i != l; ++i) {
              if(ord2[i] >= np && (temp = omap[cgs[g[i]]+ord2[i]-np])) {
                outp[i] = FUN(x[i], x[temp-1]);
              } else {
                outp[i] = fill;
              }
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1);
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            for(int i = l; i--; ) {
              if(seen[g[omap[i]]] == max[g[omap[i]]]-start) outp[omap[i]] = fill;
              else {
                outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
                ++seen[g[omap[i]]];
              }
            }
          }
          if(ds > 1) {
            NumericVector outtemp = outp;
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1];
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int start = np*(k+1);
                std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                for(int i = l; i--; ) {
                  if(seen[g[omap[i]]] == max[g[omap[i]]]-start) outtemp[omap[i]] = fill;
                  else {
                    outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                    ++seen[g[omap[i]]];
                  }
                }
              }
              out( _ , pos) = outtemp;
              if(names) {
                if(L1) colnam[pos] = stub + diffc[q];
                else colnam[pos] = "L" + nc[p] + stub + diffc[q];
              }
              ++pos;
            }
          }
        } else if(np<0) { // (Negative) leaded and iterated differences
          int d1 = diff[0];
          bool F1 = np == -1;
          if(d1 < 1) stop("diff must be a vector of integers > 0");
          NumericMatrix::Column outp = out( _ , pos);
          if(names) {
            if(F1) colnam[pos] = stub2 + diffc[0];
            else colnam[pos] = "F" + nc[p] + stub + diffc[0];
          }
          ++pos;
          if(regular) {
            for(int i = 0; i != l; ++i) {
              if(ord2[i] < max[g[i]]+np) {
                outp[i] = FUN(x[i], x[omap[cgs[g[i]]+ord2[i]-np]]);
              } else {
                outp[i] = fill;
              }
            }
          } else {
            for(int i = 0; i != l; ++i) { // Smarter solution using while ???
              if(ord2[i] < max[g[i]]+np && (temp = omap[cgs[g[i]]+ord2[i]-np])) {
                outp[i] = FUN(x[i], x[temp-1]);
              } else {
                outp[i] = fill;
              }
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1);
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            for(int i = 0; i != l; ++i) {
              if(seen[g[omap[i]]] == max[g[omap[i]]]+start) outp[omap[i]] = fill;
              else {
                outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
                ++seen[g[omap[i]]];
              }
            }
          }
          if(ds > 1) {
            NumericVector outtemp = outp;
            for(int q = 1; q != ds; ++q) {
              int dq = diff[q], L_dq = diff[q-1];
              if(dq <= L_dq) stop("differences must be passed in ascending order");
              for(int k = L_dq; k != dq; ++k) {
                int start = np*(k+1);
                std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                for(int i = 0; i != l; ++i) {
                  if(seen[g[omap[i]]] == max[g[omap[i]]]+start) outtemp[omap[i]] = fill;
                  else {
                    outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                    ++seen[g[omap[i]]];
                  }
                }
              }
              out( _ , pos) = outtemp;
              if(names) {
                if(F1) colnam[pos] = stub2 + diffc[q];
                else colnam[pos] = "F"+ nc[p] + stub + diffc[q];
              }
              ++pos;
            }
          }
        } else {
          out( _ , pos) = x;
          if(names) colnam[pos] = "--";
          ++pos;
        }
      }
    }
  }
  // Previous Version
  // if(ncol == 1) SHALLOW_DUPLICATE_ATTRIB(out, x);
  // else if(names) out.attr("dimnames") = List::create(x.attr("names"), colnam);

  SHALLOW_DUPLICATE_ATTRIB(out, x);
  if(ncol != 1) {
    Rf_setAttrib(out, R_NamesSymbol, R_NilValue); // if(x.hasAttribute("names")) out.attr("names") = R_NilValue;
    Rf_dimgets(out, Dimension(l, ncol));
    if(Rf_isObject(x)) { //  && !x.inherits("pseries") -> lag matrix in plm is not a pseries anymore anyway...
      CharacterVector classes = Rf_getAttrib(out, R_ClassSymbol);
      classes.push_back("matrix");
      Rf_classgets(out, classes);
    } // else {
      // Rf_classgets(out, Rf_mkString("matrix"));
      // }
    if(names) Rf_dimnamesgets(out, List::create(Rf_getAttrib(x, R_NamesSymbol), colnam));
  }

  return out;
}


// [[Rcpp::export]]
NumericVector fdiffgrowthCpp(const NumericVector& x, const IntegerVector& n = 1, const IntegerVector& diff = 1,
                             double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0,
                             const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue,
                             int ret = 1, double rho = 1, bool names = true, double power = 1) {

  std::string stub;
  if(ret < 4) {
    double rho2;
    if(ret == 3) {
      rho2 = 1;
      if(power != 1) stop("High-powered log-difference growth rates are currently not supported");
      if(names) stub = "Dlog";
    } else {
      rho2 = rho;
      if(names) stub = (ret == 1 && rho == 1) ? "D" : (ret == 1) ? "QD" : (rho == 1) ? "Dlog" : "QDlog"; // QD for quasi-differences
    }
    return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho2](double y, double x) { return y-rho2*x; }); // return y-x; same efficiency as return y-rho*x; when rho = 1 -> smart compiler !, and reduced file size !!
  } else if (ret == 4) {
    if(names) stub = "G";                   // same speed as fixing 100 ! Faster using (y/x-1)*rho or (x*(1/x)-1)*rho ?
    if(power == 1) return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return (y-x)*(rho/x); }); // definitely much faster !!
    return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho, power](double y, double x) { return (pow(y/x, power)-1)*rho; }); // without: 375 kb
  } else stop("Unknown return option!");
}

inline SEXP coln_check(SEXP x) {
  return Rf_isNull(x) ? NA_STRING : x;
}


template <typename F>
NumericMatrix fdiffgrowthmCppImpl(const NumericMatrix& x, const IntegerVector& n = 1, const IntegerVector& diff = 1,
                                  double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0,
                                  const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue,
                                  std::string stub = "", bool names = true, F FUN = [](double y, double x) { return y-x; }) {

  int l = x.nrow(), col = x.ncol(), ns = n.size(), ds = diff.size(), zeros = 0, pos = INT_MAX;
  IntegerVector absn = no_init_vector(ns);
  for(int i = ns; i--; ) {
    if(n[i] == pos) stop("duplicated values in n detected"); // because one might mistakenly pass a factor to the n-slot
    pos = n[i];
    if(pos == 0) ++zeros;
    if(pos < 0) {
      if(pos == NA_INTEGER) stop("NA in n");
      absn[i] = -pos;
    } else absn[i] = pos;
  }
  pos = 0;
  std::string stub2 = names ? "F" + stub : "";

  int ncol = ((ns-zeros)*ds+zeros)*col;
  NumericMatrix out = no_init_matrix(l, ncol);
  CharacterVector colnam = names ? no_init_vector(ncol) : no_init_vector(1);
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING;
  CharacterVector diffc = names ? Rf_coerceVector(diff, STRSXP) : NA_STRING;
  CharacterVector coln = names ? coln_check(colnames(x)) : NA_STRING;
  if(names && coln[0] == NA_STRING) names = false;

  if(ng == 0) { // No groups
    if(Rf_isNull(t)) { // Ordered data
      for(int j = 0; j != col; ++j) {
        NumericMatrix::ConstColumn column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0], end = np*d1;
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            if(end >= l) stop("n * diff needs to be < nrow(x)");
            NumericMatrix::Column outp = out( _ , pos);
            if(names) {
              if(L1) colnam[pos] = stub + diffc[0] + "." + coln[j];
              else colnam[pos] = "L" + nc[p] + stub + diffc[0] + "." + coln[j];
            }
            ++pos;
            for(int i = np; i != l; ++i) outp[i] = FUN(column[i], column[i - np]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1)-1;
              for(int i = l-1; i != start; --i) outp[i] = FUN(outp[i], outp[i - np]);
            }
            for(int i = end; i--; ) outp[i] = fill;
            if(ds > 1) {
              NumericVector outtemp = outp;
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = np*dq;
                if(end >= l) stop("n * diff needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1)-1;
                  for(int i = l-1; i != start; --i) outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                }
                for(int i = np*L_dq; i != end; ++i) outtemp[i] = fill;
                out( _ , pos) = outtemp;
                if(names) {
                  if(L1) colnam[pos] = stub + diffc[q] + "." + coln[j];
                  else colnam[pos] = "L" + nc[p] + stub + diffc[q] + "." + coln[j];
                }
                ++pos;
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0], end = l+np*d1;
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
            NumericMatrix::Column outp = out( _ , pos);
            if(names) {
              if(F1) colnam[pos] = stub2 + diffc[0] + "." + coln[j];
              else colnam[pos] = "F" + nc[p] + stub + diffc[0] + "." + coln[j];
            }
            ++pos;
            for(int i = l+np; i--; ) outp[i] = FUN(column[i], column[i - np]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int final = l+np*(k+1);
              for(int i = 0; i != final; ++i) outp[i] = FUN(outp[i], outp[i - np]);
            }
            for(int i = end; i != l; ++i) outp[i] = fill;
            if(ds > 1) {
              NumericVector outtemp = outp;
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
                if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int final = l+np*(k+1);
                  for(int i = 0; i != final; ++i) outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                }
                for(int i = end; i != start; ++i) outtemp[i] = fill;
                out( _ , pos) = outtemp;
                if(names) {
                  if(F1) colnam[pos] = stub2 + diffc[q] + "." + coln[j];
                  else colnam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + coln[j];
                }
                ++pos;
              }
            }
          } else {
            out( _ , pos) = column;
            if(names) colnam[pos] = coln[j];
            ++pos;
          }
        }
      }
    } else { // Unordered data: Timevar provided
      IntegerVector ord = t;
      if(l != ord.size()) stop("length(x) must match length(t)");
      int min = INT_MAX, max = INT_MIN, osize, temp;
      for(int i = 0; i != l; ++i) {
        if(ord[i] < min) min = ord[i];
        if(ord[i] > max) max = ord[i];
      }
      if(min == NA_INTEGER) stop("Timevar contains missing values");
      osize = max-min+1;
      bool regular = osize == l;
      IntegerVector omap(osize), ord2 = regular ? no_init_vector(1) : no_init_vector(l);
      if(!regular) { // Irregular time series
        if(osize > 10000000 && osize > 3 * l) warning("Your time series is very irregular. Need to create an internal ordering vector of length %s to represent it.", osize);
        if(Rcpp::max(diff) > 1) stop("Iterations are currently only supported for regular time series. See ?seqid to identify the regular sequences in your time series, or just apply this function multiple times.");
        for(int i = 0; i != l; ++i) {
          temp = ord[i] - min; // Best ? Or direct assign to ord2[i] ? Also check for panel version..
          if(omap[temp]) stop("Repeated values in timevar");
          omap[temp] = i+1;
          ord2[i] = temp;
        }
      } else { // Regular time series
        for(int i = 0; i != l; ++i) {
          temp = ord[i] - min;
          if(omap[temp]) stop("Repeated values in timevar");
          omap[temp] = i;
        }
      }
      for(int j = 0; j != col; ++j) {
        NumericMatrix::ConstColumn column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0], end = np*d1;
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            if(end >= l) stop("n * diff needs to be < nrow(x)");
            NumericMatrix::Column outp = out( _ , pos);
            if(names) {
              if(L1) colnam[pos] = stub + diffc[0] + "." + coln[j];
              else colnam[pos] = "L" + nc[p] + stub + diffc[0] + "." + coln[j];
            }
            ++pos;
            if(regular) {
              for(int i = np; i != l; ++i) outp[omap[i]] = FUN(column[omap[i]], column[omap[i - np]]);
            } else {
              for(int i = 0; i != l; ++i) { // Smarter solution using while ???
                if(ord2[i] >= np && (temp = omap[ord2[i] - np])) {
                  outp[i] = FUN(column[i], column[temp-1]);
                } else {
                  outp[i] = fill;
                }
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1)-1;
              for(int i = l-1; i != start; --i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
            }
            if(regular) for(int i = end; i--; ) outp[omap[i]] = fill;
            if(ds > 1) {
              NumericVector outtemp = outp;
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = np*dq;
                if(end >= l) stop("n * diff needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1)-1;
                  for(int i = l-1; i != start; --i) outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                }
                for(int i = np*L_dq; i != end; ++i) outtemp[omap[i]] = fill;
                out( _ , pos) = outtemp;
                if(names) {
                  if(L1) colnam[pos] = stub + diffc[q] + "." + coln[j];
                  else colnam[pos] = "L" + nc[p] + stub + diffc[q] + "." + coln[j];
                }
                ++pos;
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0], end = l+np*d1;
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
            NumericMatrix::Column outp = out( _ , pos);
            if(names) {
              if(F1) colnam[pos] = stub2 + diffc[0] + "." + coln[j];
              else colnam[pos] = "F" + nc[p] + stub + diffc[0] + "." + coln[j];
            }
            ++pos;
            if(regular) {
              for(int i = l+np; i--; ) outp[omap[i]] = FUN(column[omap[i]], column[omap[i - np]]);
            } else {
              for(int i = 0, osnp = osize + np; i != l; ++i) { // Smarter solution using while ???
                if(ord2[i] < osnp && (temp = omap[ord2[i] - np])) {
                  outp[i] = FUN(column[i], column[temp-1]);
                } else {
                  outp[i] = fill;
                }
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int final = l+np*(k+1);
              for(int i = 0; i != final; ++i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
            }
            if(regular) for(int i = end; i != l; ++i) outp[omap[i]] = fill;
            if(ds > 1) {
              NumericVector outtemp = outp;
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = l+np*dq, start = l+np*L_dq;
                if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int final = l+np*(k+1);
                  for(int i = 0; i != final; ++i) outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                }
                for(int i = end; i != start; ++i) outtemp[omap[i]] = fill;
                out( _ , pos) = outtemp;
                if(names) {
                  if(F1) colnam[pos] = stub2 + diffc[q] + "." + coln[j];
                  else colnam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + coln[j];
                }
                ++pos;
              }
            }
          } else {
            out( _ , pos) = column;
            if(names) colnam[pos] = coln[j];
            ++pos;
          }
        }
      }
    }
  } else { // With groups
    if(l != g.size()) stop("nrow(x) must match length(g)");
    int ags = l/ng, ngp = ng+1, maxdiff = max(diff);
    if(Rf_isNull(t)) { // Ordered data
      bool cond = !Rf_isNull(gs);
      IntegerVector gsv = (cond || maxdiff == 1) ? no_init_vector(1) : IntegerVector(ng);
      int *pgsv = cond ? INTEGER(gs)-1 : INTEGER(gsv)-1;
      if(maxdiff != 1) {
        if(cond) {
          if(ng != Rf_length(gs)) stop("ng must match length(gs)");
        } else {
          for(int i = 0; i != l; ++i) ++pgsv[g[i]];
        }
      }
      // int seen[ngp], memsize = sizeof(int)*(ngp);
      for(int j = 0; j != col; ++j) {
        NumericMatrix::ConstColumn column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0];
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            NumericMatrix::Column outp = out( _ , pos);
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            if(names) {
              if(L1) colnam[pos] = stub + diffc[0] + "." + coln[j];
              else colnam[pos] = "L" + nc[p] + stub + diffc[0] + "." + coln[j];
            }
            ++pos;
            for(int i = 0; i != l; ++i) {
              if(seen[g[i]] == np) outp[i] = FUN(column[i], column[i - np]);
              else {
                outp[i] = fill;
                ++seen[g[i]];
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = l; i--; ) {
                if(seen[g[i]] == pgsv[g[i]]-start) outp[i] = fill;
                else {
                  outp[i] = FUN(outp[i], outp[i - np]);
                  ++seen[g[i]];
                }
              }
            }
            if(ds > 1) {
              NumericVector outtemp = outp;
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1);
                  std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                  for(int i = l; i--; ) {
                    if(seen[g[i]] == pgsv[g[i]]-start) outtemp[i] = fill;
                    else {
                      outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                      ++seen[g[i]];
                    }
                  }
                }
                out( _ , pos) = outtemp;
                if(names) {
                  if(L1) colnam[pos] = stub + diffc[q] + "." + coln[j];
                  else colnam[pos] = "L" + nc[p] + stub + diffc[q] + "." + coln[j];
                }
                ++pos;
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0];
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            NumericMatrix::Column outp = out( _ , pos);
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            if(names) {
              if(F1) colnam[pos] = stub2 + diffc[0] + "." + coln[j];
              else colnam[pos] = "F" + nc[p] + stub + diffc[0] + "." + coln[j];
            }
            ++pos;
            for(int i = l; i--; ) {
              if(seen[g[i]] == np) outp[i] = FUN(column[i], column[i - np]);
              else {
                outp[i] = fill;
                --seen[g[i]];
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = 0; i != l; ++i) {
                if(seen[g[i]] == pgsv[g[i]]+start) outp[i] = fill;
                else {
                  outp[i] = FUN(outp[i], outp[i - np]);
                  ++seen[g[i]];
                }
              }
            }
            if(ds > 1) {
              NumericVector outtemp = outp;
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1);
                  std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                  for(int i = 0; i != l; ++i) {
                    if(seen[g[i]] == pgsv[g[i]]+start) outtemp[i] = fill;
                    else {
                      outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                      ++seen[g[i]];
                    }
                  }
                }
                out( _ , pos) = outtemp;
                if(names) {
                  if(F1) colnam[pos] = stub2 + diffc[q] + "." + coln[j];
                  else colnam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + coln[j];
                }
                ++pos;
              }
            }
          } else {
            out( _ , pos) = column;
            if(names) colnam[pos] = coln[j];
            ++pos;
          }
        }
      }
    } else { // Unordered data: Timevar Provided
      IntegerVector ord = t;
      int temp;
      if(l != ord.size()) stop("length(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX), max(ngp, INT_MIN), cgs = no_init_vector(ngp);
      for(int i = 0; i != l; ++i) {
        temp = g[i];
        if(ord[i] < min[temp]) min[temp] = ord[i];
        if(ord[i] > max[temp]) max[temp] = ord[i];
      }
      temp = 0;
      for(int i = 1; i != ngp; ++i) {
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values");
        if(min[i] == INT_MAX) continue; // Needed in case of unused factor levels (group vector too large)
        cgs[i] = temp; // This needs to b here (for unused factor levels case...)
        max[i] -= min[i] - 1; // need max[i] which stores the complete group sizes only if p<0 e.g. if computing leads..
        temp += max[i];
      }
      // omap provides the ordering to order the vector (needed to find previous / next values)
      bool regular = temp == l;
      IntegerVector omap(temp), ord2 = no_init_vector(l), index = no_init_vector(l);
      if(!regular) { // Irregular panel
        if(temp > 10000000 && temp > 3 * l) warning("Your panel is very irregular. Need to create an internal ordering vector of length %s to represent it.", temp);
        if(maxdiff > 1) stop("Iterations are currently only supported for regular panels. See ?seqid to identify the regular sequences in your panel, or just apply this function multiple times.");
        for(int i = 0; i != l; ++i) {
          ord2[i] = ord[i] - min[g[i]];
          index[i] = cgs[g[i]] + ord2[i];
          if(omap[index[i]]) stop("Repeated values of timevar within one or more groups");
          omap[index[i]] = i+1; // needed to add 1 to distinguish between 0 and gap
        }
      } else { // Regular panel
        for(int i = 0; i != l; ++i) {
          ord2[i] = ord[i] - min[g[i]];
          index[i] = cgs[g[i]] + ord2[i];
          if(omap[index[i]]) stop("Repeated values of timevar within one or more groups");
          omap[index[i]] = i;
        }
      }

      for(int j = 0; j != col; ++j) {
        NumericMatrix::ConstColumn column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0];
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            NumericMatrix::Column outp = out( _ , pos);
            if(names) {
              if(L1) colnam[pos] = stub + diffc[0] + "." + coln[j];
              else colnam[pos] = "L" + nc[p] + stub + diffc[0] + "." + coln[j];
            }
            ++pos;
            if(regular) {
              for(int i = 0; i != l; ++i) {
                if(ord2[i] >= np) {
                  outp[i] = FUN(column[i], column[omap[index[i]-np]]);
                } else {
                  outp[i] = fill;
                }
              }
            } else {
              for(int i = 0; i != l; ++i) {
                if(ord2[i] >= np && (temp = omap[index[i]-np])) {
                  outp[i] = FUN(column[i], column[temp-1]);
                } else {
                  outp[i] = fill;
                }
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = l; i--; ) {
                if(seen[g[omap[i]]] == max[g[omap[i]]]-start) outp[omap[i]] = fill;
                else {
                  outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
                  ++seen[g[omap[i]]];
                }
              }
            }
            if(ds > 1) {
              NumericVector outtemp = outp;
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1);
                  std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                  for(int i = l; i--; ) {
                    if(seen[g[omap[i]]] == max[g[omap[i]]]-start) outtemp[omap[i]] = fill;
                    else {
                      outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                      ++seen[g[omap[i]]];
                    }
                  }
                }
                out( _ , pos) = outtemp;
                if(names) {
                  if(L1) colnam[pos] = stub + diffc[q] + "." + coln[j];
                  else colnam[pos] = "L" + nc[p] + stub + diffc[q] + "." + coln[j];
                }
                ++pos;
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0];
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            NumericMatrix::Column outp = out( _ , pos);
            if(names) {
              if(F1) colnam[pos] = stub2 + diffc[0] + "." + coln[j];
              else colnam[pos] = "F" + nc[p] + stub + diffc[0] + "." + coln[j];
            }
            ++pos;
            if(regular) {
              for(int i = 0; i != l; ++i) {
                if(ord2[i] < max[g[i]]+np) {
                  outp[i] = FUN(column[i], column[omap[index[i]-np]]);
                } else {
                  outp[i] = fill;
                }
              }
            } else {
              for(int i = 0; i != l; ++i) { // Smarter solution using while ???
                if(ord2[i] < max[g[i]]+np && (temp = omap[index[i]-np])) {
                  outp[i] = FUN(column[i], column[temp-1]);
                } else {
                  outp[i] = fill;
                }
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = 0; i != l; ++i) {
                if(seen[g[omap[i]]] == max[g[omap[i]]]+start) outp[omap[i]] = fill;
                else {
                  outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
                  ++seen[g[omap[i]]];
                }
              }
            }
            if(ds > 1) {
              NumericVector outtemp = outp;
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1);
                  std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                  for(int i = 0; i != l; ++i) {
                    if(seen[g[omap[i]]] == max[g[omap[i]]]+start) outtemp[omap[i]] = fill;
                    else {
                      outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                      ++seen[g[omap[i]]];
                    }
                  }
                }
                out( _ , pos) = outtemp;
                if(names) {
                  if(F1) colnam[pos] = stub2 + diffc[q] + "." + coln[j];
                  else colnam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + coln[j];
                }
                ++pos;
              }
            }
          } else {
            out( _ , pos) = column;
            if(names) colnam[pos] = coln[j];
            ++pos;
          }
        }
      }
    }
  }
  // Previous Solution:
  // if(names) {
  //   out.attr("dimnames") = List::create(rownames(x), colnam);
  // } else {
  //   if(ns*ds == 1) SHALLOW_DUPLICATE_ATTRIB(out, x);
  //   // else rownames(out) = rownames(x); // redundant !!
  // }

  SHALLOW_DUPLICATE_ATTRIB(out, x);
  if(ncol != col) Rf_dimgets(out, Dimension(l, ncol));
  if(names) {
    Rf_dimnamesgets(out, List::create(rownames(x), colnam)); // colnames(out) = colnam; also deletes rownames !
  } else if(ncol != col) {
    Rf_setAttrib(out, R_DimNamesSymbol, R_NilValue);
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix fdiffgrowthmCpp(const NumericMatrix& x, const IntegerVector& n = 1, const IntegerVector& diff = 1,
                              double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0,
                              const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue,
                              int ret = 1, double rho = 1, bool names = true, double power = 1) {
  std::string stub;
  if(ret < 4) {
    double rho2;
    if(ret == 3) {
      rho2 = 1;
      if(power != 1) stop("High-powered log-difference growth rates are currently not supported");
      if(names) stub = "Dlog";
    } else {
      rho2 = rho;
      if(names) stub = (ret == 1 && rho == 1) ? "D" : (ret == 1) ? "QD" : (rho == 1) ? "Dlog" : "QDlog"; // QD for quasi-differences
    }
    return fdiffgrowthmCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho2](double y, double x) { return y-rho2*x; }); // return y-x; same efficiency as return y-rho*x; when rho = 1 -> smart compiler !, and reduced file size !!
  } else if (ret == 4) {
    if(names) stub = "G";
    if(power == 1) return fdiffgrowthmCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return (y-x)*(rho/x); }); // same speed as fixing 100 ! Faster using (y/x-1)*rho or (x*(1/x)-1)*rho ?
    return fdiffgrowthmCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho, power](double y, double x) { return (pow(y/x, power)-1)*rho; });
  } else stop("Unknown return option!");
}


template <typename F>
List fdiffgrowthlCppImpl(const List& x, const IntegerVector& n = 1, const IntegerVector& diff = 1,
                         double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0,
                         const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue,
                         std::string stub = "", bool names = true, F FUN = [](double y, double x) { return y-x; }) { // const needed for #if response...

  int l = x.size(), ns = n.size(), ds = diff.size(), zeros = 0, pos = INT_MAX;
  IntegerVector absn = no_init_vector(ns);
  for(int i = ns; i--; ) {
    if(n[i] == pos) stop("duplicated values in n detected"); // because one might mistakenly pass a factor to the n-slot
    pos = n[i];
    if(pos == 0) ++zeros;
    if(pos < 0) {
      if(pos == NA_INTEGER) stop("NA in n");
      absn[i] = -pos;
    } else absn[i] = pos;
  }
  pos = 0;
  std::string stub2 = names ? "F" + stub : "";

  int ncol = ((ns-zeros)*ds+zeros)*l;
  List out(ncol);
  CharacterVector nam = names ? no_init_vector(ncol) : no_init_vector(1);
  CharacterVector nc = names ? Rf_coerceVector(absn, STRSXP) : NA_STRING;
  CharacterVector diffc = names ? Rf_coerceVector(diff, STRSXP) : NA_STRING;
  CharacterVector na = names ? coln_check(Rf_getAttrib(x, R_NamesSymbol)) : NA_STRING;
  if(names && na[0] == NA_STRING) names = false;

  if(ng == 0) { // No groups
    if(Rf_isNull(t)) { // Ordered data
      for(int j = 0; j != l; ++j) {
        NumericVector column = x[j];
        int row = column.size();
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0], end = np*d1;
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            if(end >= row) stop("n * diff needs to be < nrow(x)");
            NumericVector outjp = no_init_vector(row);
            if(names) {
              if(L1) nam[pos] = stub + diffc[0] + "." + na[j];
              else nam[pos] = "L" + nc[p] + stub + diffc[0] + "." + na[j];
            }
            for(int i = np; i != row; ++i) outjp[i] = FUN(column[i], column[i - np]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1)-1;
              for(int i = row-1; i != start; --i) outjp[i] = FUN(outjp[i], outjp[i - np]);
            }
            for(int i = end; i--; ) outjp[i] = fill;
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp);
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = np*dq;
                if(end >= row) stop("n * diff needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1)-1;
                  for(int i = row-1; i != start; --i) outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                }
                for(int i = np*L_dq; i != end; ++i) outtemp[i] = fill;
                if(names) {
                  if(L1) nam[pos] = stub + diffc[q] + "." + na[j];
                  else nam[pos] = "L" + nc[p] + stub + diffc[q] + "." + na[j];
                }
                out[pos++] = Rf_shallow_duplicate(outtemp); // or Rf_copyVector, Rf_shallow_duplicate, Rf_lazy_duplicate  http://mtweb.cs.ucl.ac.uk/mus/bin/install_R/R-3.1.1/src/main/duplicate.c
              } // https://rlang.r-lib.org/reference/duplicate.html
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0], end = row+np*d1;
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
            NumericVector outjp = no_init_vector(row);
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            if(names) {
              if(F1) nam[pos] = stub2 + diffc[0] + "." + na[j];
              else nam[pos] = "F" + nc[p] + stub + diffc[0] + "." + na[j];
            }
            for(int i = row+np; i--; ) outjp[i] = FUN(column[i], column[i - np]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int final = row+np*(k+1);
              for(int i = 0; i != final; ++i) outjp[i] = FUN(outjp[i], outjp[i - np]);
            }
            for(int i = end; i != row; ++i) outjp[i] = fill;
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp);
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = row+np*dq, start = row+np*L_dq;
                if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int final = row+np*(k+1);
                  for(int i = 0; i != final; ++i) outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                }
                for(int i = end; i != start; ++i) outtemp[i] = fill;
                if(names) {
                  if(F1) nam[pos] = stub2 + diffc[q] + "." + na[j];
                  else nam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + na[j];
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else {
            if(names) nam[pos] = na[j];
            out[pos++] = column;
          }
        }
      }
    } else { // Unordered data: Timevar provided
      IntegerVector ord = t;
      int min = INT_MAX, max = INT_MIN, osize, temp, os = ord.size();
      if(Rf_length(x[0]) != os) stop("length(x) must match length(t)");
      for(int i = 0; i != os; ++i) {
        if(ord[i] < min) min = ord[i];
        if(ord[i] > max) max = ord[i];
      }
      if(min == NA_INTEGER) stop("Timevar contains missing values");
      osize = max-min+1;
      bool regular = osize == os;
      IntegerVector omap(osize), ord2 = regular ? no_init_vector(1) : no_init_vector(os);
      if(!regular) { // Irregular time series
        if(osize > 10000000 && osize > 3 * os) warning("Your time series is very irregular. Need to create an internal ordering vector of length %s to represent it.", osize);
        if(Rcpp::max(diff) > 1) stop("Iterations are currently only supported for regular time series. See ?seqid to identify the regular sequences in your time series, or just apply this function multiple times.");
        for(int i = 0; i != os; ++i) {
          temp = ord[i] - min; // Best ? Or direct assign to ord2[i] ? Also check for panel version..
          if(omap[temp]) stop("Repeated values in timevar");
          omap[temp] = i+1;
          ord2[i] = temp;
        }
      } else { // Regular time series
        for(int i = 0; i != os; ++i) {
          temp = ord[i] - min;
          if(omap[temp]) stop("Repeated values in timevar");
          omap[temp] = i;
        }
      }

      for(int j = 0; j != l; ++j) {
        NumericVector column = x[j];
        if(os != column.size()) stop("nrow(x) must match length(t)");
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0], end = np*d1;
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            if(end >= os) stop("n * diff needs to be < nrow(x)");
            NumericVector outjp = no_init_vector(os);
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            if(names) {
              if(L1) nam[pos] = stub + diffc[0] + "." + na[j];
              else nam[pos] = "L" + nc[p] + stub + diffc[0] + "." + na[j];
            }
            if(regular) {
              for(int i = np; i != os; ++i) outjp[omap[i]] = FUN(column[omap[i]], column[omap[i - np]]);
            } else {
              for(int i = 0; i != os; ++i) { // Smarter solution using while ???
                if(ord2[i] >= np && (temp = omap[ord2[i] - np])) {
                  outjp[i] = FUN(column[i], column[temp-1]);
                } else {
                  outjp[i] = fill;
                }
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1)-1;
              for(int i = os-1; i != start; --i) outjp[omap[i]] = FUN(outjp[omap[i]], outjp[omap[i - np]]);
            }
            if(regular) for(int i = end; i--; ) outjp[omap[i]] = fill;
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp);
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = np*dq;
                if(end >= os) stop("n * diff needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1)-1;
                  for(int i = os-1; i != start; --i) outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                }
                for(int i = np*L_dq; i != end; ++i) outtemp[omap[i]] = fill;
                if(names) {
                  if(L1) nam[pos] = stub + diffc[q] + "." + na[j];
                  else nam[pos] = "L" + nc[p] + stub + diffc[q] + "." + na[j];
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0], end = os+np*d1;
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
            NumericVector outjp = no_init_vector(os);
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            if(names) {
              if(F1) nam[pos] = stub2 + diffc[0] + "." + na[j];
              else nam[pos] = "F" + nc[p] + stub + diffc[0] + "." + na[j];
            }
            if(regular) {
              for(int i = os+np; i--; ) outjp[omap[i]] = FUN(column[omap[i]], column[omap[i - np]]);
            } else {
              for(int i = 0, osnp = osize + np; i != os; ++i) { // Smarter solution using while ???
                if(ord2[i] < osnp && (temp = omap[ord2[i] - np])) {
                  outjp[i] = FUN(column[i], column[temp-1]);
                } else {
                  outjp[i] = fill;
                }
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int final = os+np*(k+1);
              for(int i = 0; i != final; ++i) outjp[omap[i]] = FUN(outjp[omap[i]], outjp[omap[i - np]]);
            }
            if(regular) for(int i = end; i != os; ++i) outjp[omap[i]] = fill;
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp);
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1], end = os+np*dq, start = os+np*L_dq;
                if(end <= 0) stop("abs(n * diff) needs to be < nrow(x)");
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int final = os+np*(k+1);
                  for(int i = 0; i != final; ++i) outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                }
                for(int i = end; i != start; ++i) outtemp[omap[i]] = fill;
                if(names) {
                  if(F1) nam[pos] = stub2 + diffc[q] + "." + na[j];
                  else nam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + na[j];
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else {
            if(names) nam[pos] = na[j];
            out[pos++] = column;
          }
        }
      }
    }
  } else { // With groups
    int gss = g.size(), ags = gss/ng, ngp = ng+1, maxdiff = max(diff);
    if(Rf_isNull(t)) { // Ordered data
      bool cond = !Rf_isNull(gs);
      IntegerVector gsv = (cond || maxdiff == 1) ? no_init_vector(1) : IntegerVector(ng);
      int *pgsv = cond ? INTEGER(gs)-1 : INTEGER(gsv)-1;
      if(maxdiff != 1) {
        if(cond) {
          if(ng != Rf_length(gs)) stop("ng must match length(gs)");
        } else {
          for(int i = 0; i != gss; ++i) ++pgsv[g[i]];
        }
      }
      // int seen[ngp], memsize = sizeof(int)*(ngp);
      for(int j = 0; j != l; ++j) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("nrow(x) must match length(g)");
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0];
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            NumericVector outjp = no_init_vector(gss);
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            if(names) {
              if(L1) nam[pos] = stub + diffc[0] + "." + na[j];
              else nam[pos] = "L" + nc[p] + stub + diffc[0] + "." + na[j];
            }
            for(int i = 0; i != gss; ++i) {
              if(seen[g[i]] == np) outjp[i] = FUN(column[i], column[i - np]);
              else {
                outjp[i] = fill;
                ++seen[g[i]];
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = gss; i--; ) {
                if(seen[g[i]] == pgsv[g[i]]-start) outjp[i] = fill;
                else {
                  outjp[i] = FUN(outjp[i], outjp[i - np]);
                  ++seen[g[i]];
                }
              }
            }
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp);
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1);
                  std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                  for(int i = gss; i--; ) {
                    if(seen[g[i]] == pgsv[g[i]]-start) outtemp[i] = fill;
                    else {
                      outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                      ++seen[g[i]];
                    }
                  }
                }
                if(names) {
                  if(L1) nam[pos] = stub + diffc[q] + "." + na[j];
                  else nam[pos] = "L" + nc[p] + stub + diffc[q] + "." + na[j];
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0];
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            NumericVector outjp = no_init_vector(gss);
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            if(names) {
              if(F1) nam[pos] = stub2 + diffc[0] + "." + na[j];
              else nam[pos] = "F" + nc[p] + stub + diffc[0] + "." + na[j];
            }
            for(int i = gss; i--; ) {
              if(seen[g[i]] == np) outjp[i] = FUN(column[i], column[i - np]);
              else {
                outjp[i] = fill;
                --seen[g[i]];
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = 0; i != gss; ++i) {
                if(seen[g[i]] == pgsv[g[i]]+start) outjp[i] = fill;
                else {
                  outjp[i] = FUN(outjp[i], outjp[i - np]);
                  ++seen[g[i]];
                }
              }
            }
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp);
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1);
                  std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                  for(int i = 0; i != gss; ++i) {
                    if(seen[g[i]] == pgsv[g[i]]+start) outtemp[i] = fill;
                    else {
                      outtemp[i] = FUN(outtemp[i], outtemp[i - np]);
                      ++seen[g[i]];
                    }
                  }
                }
                if(names) {
                  if(F1) nam[pos] = stub2 + diffc[q] + "." + na[j];
                  else nam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + na[j];
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else {
            if(names) nam[pos] = na[j];
            out[pos++] = column;
          }
        }
      }
    } else { // Unordered data: Timevar Provided
      IntegerVector ord = t;
      int temp;
      if(gss != ord.size()) stop("length(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX), max(ngp, INT_MIN), cgs = no_init_vector(ngp);
      for(int i = 0; i != gss; ++i) {
        temp = g[i];
        if(ord[i] < min[temp]) min[temp] = ord[i];
        if(ord[i] > max[temp]) max[temp] = ord[i];
      }
      temp = 0;
      for(int i = 1; i != ngp; ++i) {
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values");
        if(min[i] == INT_MAX) continue; // Needed in case of unused factor levels (group vector too large)
        cgs[i] = temp; // This needs to b here (for unused factor levels case...)
        max[i] -= min[i] - 1; // need max[i] which stores the complete group sizes only if p<0 e.g. if computing leads..
        temp += max[i];
      }
      // omap provides the ordering to order the vector (needed to find previous / next values)
      bool regular = temp == gss;
      IntegerVector omap(temp), ord2 = no_init_vector(gss), index = no_init_vector(gss);
      if(!regular) { // Irregular panel
        if(temp > 10000000 && temp > 3 * gss) warning("Your panel is very irregular. Need to create an internal ordering vector of length %s to represent it.", temp);
        if(maxdiff > 1) stop("Iterations are currently only supported for regular panels. See ?seqid to identify the regular sequences in your panel, or just apply this function multiple times.");
        for(int i = 0; i != gss; ++i) {
          ord2[i] = ord[i] - min[g[i]];
          index[i] = cgs[g[i]] + ord2[i];
          if(omap[index[i]]) stop("Repeated values of timevar within one or more groups");
          omap[index[i]] = i+1; // needed to add 1 to distinguish between 0 and gap
        }
      } else { // Regular panel
        for(int i = 0; i != gss; ++i) {
          ord2[i] = ord[i] - min[g[i]];
          index[i] = cgs[g[i]] + ord2[i];
          if(omap[index[i]]) stop("Repeated values of timevar within one or more groups");
          omap[index[i]] = i;
        }
      }

      for(int j = 0; j != l; ++j) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("nrow(x) must match length(g)");
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25. Use fdroplevels() to remove unused factor levels from your data.", ags);
          if(np>0) { // Positive lagged and iterated differences
            int d1 = diff[0];
            bool L1 = np == 1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            NumericVector outjp = no_init_vector(gss);
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            if(names) {
              if(L1) nam[pos] = stub + diffc[0] + "." + na[j];
              else nam[pos] = "L" + nc[p] + stub + diffc[0] + "." + na[j];
            }
            if(regular) {
              for(int i = 0; i != gss; ++i) {
                if(ord2[i] >= np) {
                  outjp[i] = FUN(column[i], column[omap[index[i]-np]]);
                } else {
                  outjp[i] = fill;
                }
              }
            } else {
              for(int i = 0; i != gss; ++i) {
                if(ord2[i] >= np && (temp = omap[index[i]-np])) {
                  outjp[i] = FUN(column[i], column[temp-1]);
                } else {
                  outjp[i] = fill;
                }
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = gss; i--; ) {
                if(seen[g[omap[i]]] == max[g[omap[i]]]-start) outjp[omap[i]] = fill;
                else {
                  outjp[omap[i]] = FUN(outjp[omap[i]], outjp[omap[i - np]]);
                  ++seen[g[omap[i]]];
                }
              }
            }
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp);
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1);
                  std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                  for(int i = gss; i--; ) {
                    if(seen[g[omap[i]]] == max[g[omap[i]]]-start) outtemp[omap[i]] = fill;
                    else {
                      outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                      ++seen[g[omap[i]]];
                    }
                  }
                }
                if(names) {
                  if(L1) nam[pos] = stub + diffc[q] + "." + na[j];
                  else nam[pos] = "L" + nc[p] + stub + diffc[q] + "." + na[j];
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else if(np<0) { // (Negative) leaded and iterated differences
            int d1 = diff[0];
            bool F1 = np == -1;
            if(d1 < 1) stop("diff must be a vector of integers > 0");
            NumericVector outjp = no_init_vector(gss);
            SHALLOW_DUPLICATE_ATTRIB(outjp, column);
            if(names) {
              if(F1) nam[pos] = stub2 + diffc[0] + "." + na[j];
              else nam[pos] = "F" + nc[p] + stub + diffc[0] + "." + na[j];
            }
            if(regular) {
              for(int i = 0; i != gss; ++i) {
                if(ord2[i] < max[g[i]]+np) {
                  outjp[i] = FUN(column[i], column[omap[index[i]-np]]);
                } else {
                  outjp[i] = fill;
                }
              }
            } else {
              for(int i = 0; i != gss; ++i) { // Smarter solution using while ???
                if(ord2[i] < max[g[i]]+np && (temp = omap[index[i]-np])) {
                  outjp[i] = FUN(column[i], column[temp-1]);
                } else {
                  outjp[i] = fill;
                }
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = 0; i != gss; ++i) {
                if(seen[g[omap[i]]] == max[g[omap[i]]]+start) outjp[omap[i]] = fill;
                else {
                  outjp[omap[i]] = FUN(outjp[omap[i]], outjp[omap[i - np]]);
                  ++seen[g[omap[i]]];
                }
              }
            }
            out[pos++] = outjp;
            if(ds > 1) {
              NumericVector outtemp = Rf_shallow_duplicate(outjp);
              for(int q = 1; q != ds; ++q) {
                int dq = diff[q], L_dq = diff[q-1];
                if(dq <= L_dq) stop("differences must be passed in ascending order");
                for(int k = L_dq; k != dq; ++k) {
                  int start = np*(k+1);
                  std::vector<int> seen(ngp); // memset(seen, 0, memsize);
                  for(int i = 0; i != gss; ++i) {
                    if(seen[g[omap[i]]] == max[g[omap[i]]]+start) outtemp[omap[i]] = fill;
                    else {
                      outtemp[omap[i]] = FUN(outtemp[omap[i]], outtemp[omap[i - np]]);
                      ++seen[g[omap[i]]];
                    }
                  }
                }
                if(names) {
                  if(F1) nam[pos] = stub2 + diffc[q] + "." + na[j];
                  else nam[pos] = "F"+ nc[p] + stub + diffc[q] + "." + na[j];
                }
                out[pos++] = Rf_shallow_duplicate(outtemp);
              }
            }
          } else {
            if(names) nam[pos] = na[j];
            out[pos++] = column;
          }
        }
      }
    }
  }
  SHALLOW_DUPLICATE_ATTRIB(out, x);
  if(names) { // best way to code this ?
    Rf_namesgets(out, nam);
  } else if(ncol != l) {
    Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
  }
  return out;
}

// [[Rcpp::export]]
List fdiffgrowthlCpp(const List& x, const IntegerVector& n = 1, const IntegerVector& diff = 1,
                     double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0,
                     const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue,
                     int ret = 1, double rho = 1, bool names = true, double power = 1) {

  std::string stub;
  if(ret < 4) {
    double rho2;
    if(ret == 3) {
      rho2 = 1;
      if(power != 1) stop("High-powered log-difference growth rates are currently not supported");
      if(names) stub = "Dlog";
    } else {
      rho2 = rho;
      if(names) stub = (ret == 1 && rho == 1) ? "D" : (ret == 1) ? "QD" : (rho == 1) ? "Dlog" : "QDlog"; // QD for quasi-differences
    }
    return fdiffgrowthlCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho2](double y, double x) { return y-rho2*x; }); // return y-x; same efficiency as return y-rho*x; when rho = 1 -> smart compiler !, and reduced file size !!
  } else if (ret == 4) {
    if(names) stub = "G";
    if(power == 1) return fdiffgrowthlCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return (y-x)*(rho/x); }); // same speed as fixing 100 ! Faster using (y/x-1)*rho or (x*(1/x)-1)*rho ?
    return fdiffgrowthlCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho, power](double y, double x) { return (pow(y/x, power)-1)*rho; });
  } else stop("Unknown return option!");
}

// Old attempts without template ....
// #define FUN(y, x) (ret == 1 && rho1) ? ((y)-(x)) :
// (ret == 1) ? ((y)-rho*(x)) :
//   (ret == 2 && rho1) ? (log((y)*(1/(x)))) :
//   (ret == 2) ? (log((y)*(1/(rho*(x))))) :
//   (ret == 3) ? (((y)-(x))*(100/(x))) : (log((y)*(1/(x)))*100)

// #define rho1 (rho == 1)
// #define retm (ret)
//
// #if retm == 1 && rho1
// #define FUN(y, x) ((y)-(x))
// #elif retm == 1
// #define FUN(y, x) ((y)-rho*(x))
// #elif retm == 2 && rho1
// #define FUN(y, x) (log((y)*(1/(x))))
// #elif retm == 2
// #define FUN(y, x) (log((y)*(1/(rho*(x)))))
// #elif retm == 3
// #define FUN(y, x) (((y)-(x))*(100/(x)))
// #elif retm == 4
// #define FUN(y, x) (log((y)*(1/(x)))*100)
// #endif


// Previous: Internally computing log-differences--- compiled file was 648 kb, without debug info !!


// // [[Rcpp::export]]
// NumericVector fdiffgrowthCpp(const NumericVector& x, const IntegerVector& n = 1, const IntegerVector& diff = 1,
//                              double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0,
//                              const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue,
//                              int ret = 1, double rho = 1, bool names = true) {
//
//   std::string stub;
//   switch (ret)
//   {       // [rho] or [&rho] ?  // https://stackoverflow.com/questions/30217956/error-variable-cannot-be-implicitly-captured-because-no-default-capture-mode-h
//   case 1:
//     if(names) stub = (rho == 1) ? "D" : "QD"; // QD for quasi-differences !
//     return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return y-rho*x; }); // return y-x; same efficiency as return y-rho*x; when rho = 1 -> smart compiler !, and reduced file size !!
//   case 2:
//     if(rho == 1) goto fastld;
//     if(names) stub = "QDlog";
//     return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return R_log(y)-rho*R_log(x); }); // log(y*(1/(rho*x))) gives log(y) - log(rho*x), but we want log(y) - rho*log(x)
//   case 3:
//     if(names) stub = "G";
//     return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return (y-x)*(rho/x); }); // same speed as fixing 100 ! Faster using (y/x-1)*rho or (x*(1/x)-1)*rho ?
//   case 4:
//     fastld:
//     if(names) stub = "Dlog";
//     return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return rho*R_log(y*(1/x)); });
//   default: stop("Unknown return option!");
//   }
// }
