// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

// new setup: ret = 1L - differences, ret = 2L - log differences, ret = 3L - exact growth rates, ret = 4L - log-difference growth
// also: if rho != 1, quasi-differencing and log differencing with rho... i.e. for cochrane-orcutt regression

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
      LogicalVector ocheck(l, true);
      IntegerVector omap = no_init_vector(l); // int omap[l];
      for(int i = 0; i != l; ++i) {
        if(ord[i] > l || ord[i] < 1) stop("t needs to be a factor or integer vector of time-periods between 1 and length(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i;
        } else {
          stop("Repeated values in timevar");
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
          for(int i = np; i != l; ++i) outp[omap[i]] = FUN(x[omap[i]], x[omap[i - np]]);
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1)-1;
            for(int i = l-1; i != start; --i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
          }
          for(int i = end; i--; ) outp[omap[i]] = fill;
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
          for(int i = l+np; i--; ) outp[omap[i]] = FUN(x[omap[i]], x[omap[i - np]]);
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int final = l+np*(k+1);
            for(int i = 0; i != final; ++i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
          }
          for(int i = end; i != l; ++i) outp[omap[i]] = fill;
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
    IntegerVector gsv = Rf_isNull(gs) ? IntegerVector(ng) : as<IntegerVector>(gs); // no_init_vector(ng);
    if(Rf_isNull(t)) {
      if(maxdiff != 1) {
        if(Rf_isNull(gs)) {
          // gsv = IntegerVector(ng);
          // std::fill(gsv.begin(), gsv.end(), 0);
          for(int i = 0; i != l; ++i) ++gsv[g[i]-1];
        } else {
          // gsv = gs;
          if(ng != gsv.size()) stop("ng must match length(gs)");
        }
      }
      // int seen[ngp], memsize = sizeof(int)*(ngp);
      for(int p = 0; p != ns; ++p) {
        int np = n[p];
        if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
              if(seen[g[i]] == gsv[g[i]-1]-start) outp[i] = fill;
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
                  if(seen[g[i]] == gsv[g[i]-1]-start) outtemp[i] = fill;
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
              if(seen[g[i]] == gsv[g[i]-1]+start) outp[i] = fill;
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
                  if(seen[g[i]] == gsv[g[i]-1]+start) outtemp[i] = fill;
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
      if(l != ord.size()) stop("length(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX);
      IntegerVector ord2 = no_init_vector(l);
      if(Rf_isNull(gs)) {
        // gsv = IntegerVector(ng);
        // std::fill(gsv.begin(), gsv.end(), 0);
        for(int i = 0; i != l; ++i) {
          ++gsv[g[i]-1];
          if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
        }
      } else {
        // gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
      }
      IntegerVector omap(l), cgs = no_init_vector(ngp);
      // int cgs[ngp], seen[ngp], memsize = sizeof(int)*(ngp);
      cgs[1] = 0;
      for(int i = 1; i != ng; ++i) {
        cgs[i+1] = cgs[i] + gsv[i-1]; // or get "starts from forderv"
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?
      }
      if(min[ng] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?

      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]];
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        if(omap[cgs[g[i]]+ord2[i]] == 0) omap[cgs[g[i]]+ord2[i]] = i;
        else stop("Repeated values of timevar within one or more groups");
      }
      for(int p = 0; p != ns; ++p) {
        int np = n[p];
        if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
          for(int i = 0; i != l; ++i) {
            if(ord2[i] >= np) {
              outp[i] = FUN(x[i], x[omap[cgs[g[i]]+ord2[i]-np]]);
            } else {
              outp[i] = fill;
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1);
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            for(int i = l; i--; ) {
              if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outp[omap[i]] = fill;
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
                  if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outtemp[omap[i]] = fill;
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
          for(int i = 0; i != l; ++i) {
            if(ord2[i] < gsv[g[i]-1]+np) {
              outp[i] = FUN(x[i], x[omap[cgs[g[i]]+ord2[i]-np]]);
            } else {
              outp[i] = fill;
            }
          }
          if(d1 > 1) for(int k = 1; k != d1; ++k) {
            int start = np*(k+1);
            std::vector<int> seen(ngp); // memset(seen, 0, memsize);
            for(int i = 0; i != l; ++i) {
              if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outp[omap[i]] = fill;
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
                  if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outtemp[omap[i]] = fill;
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
  // if(ncol == 1) DUPLICATE_ATTRIB(out, x);
  // else if(names) out.attr("dimnames") = List::create(x.attr("names"), colnam);

  DUPLICATE_ATTRIB(out, x);
  if(ncol != 1) {
    Rf_setAttrib(out, R_NamesSymbol, R_NilValue); // if(x.hasAttribute("names")) out.attr("names") = R_NilValue;
    out.attr("dim") = Dimension(l, ncol);
    if(Rf_isObject(x)) {
      CharacterVector classes = out.attr("class");
      classes.push_back("matrix");
      out.attr("class") = classes;
    } else {
      out.attr("class") = "matrix";
    }
    if(names) out.attr("dimnames") = List::create(x.attr("names"), colnam);
  }

  return out;
}


// [[Rcpp::export]]
NumericVector fdiffgrowthCpp(const NumericVector& x, const IntegerVector& n = 1, const IntegerVector& diff = 1,
                         double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0,
                         const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue,
                         int ret = 1, double rho = 1, bool names = true) {

  std::string stub;
  switch (ret)
  {       // [rho] or [&rho] ?  // https://stackoverflow.com/questions/30217956/error-variable-cannot-be-implicitly-captured-because-no-default-capture-mode-h
  case 1:
    if(names) stub = (rho == 1) ? "D" : "QD"; // QD for quasi-differences !
    return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return y-rho*x; }); // return y-x; same efficiency as return y-rho*x; when rho = 1 -> smart compiler !, and reduced file size !!
  case 2:
    if(rho == 1) goto fastld;
    if(names) stub = "QDlog";
    return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return std::log(y)-rho*std::log(x); }); // log(y*(1/(rho*x))) gives log(y) - log(rho*x), but we want log(y) - rho*log(x)
  case 3:
    if(names) stub = "G";
    return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return (y-x)*(rho/x); }); // same speed as fixing 100 ! Faster using (y/x-1)*rho or (x*(1/x)-1)*rho ?
  case 4:
    fastld:
    if(names) stub = "Dlog";
    return fdiffgrowthCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return rho*std::log(y*(1/x)); });
  default: stop("Unknown return option!");
  }
}

inline SEXP coln_check(SEXP x) {
  if(Rf_isNull(x)) return NA_STRING;
  else return x;
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
      if(l != ord.size()) stop("nrow(x) must match length(t)");
      LogicalVector ocheck(l, true);
      IntegerVector omap = no_init_vector(l); // int omap[l];
      for(int i = 0; i != l; ++i) {
        if(ord[i] > l || ord[i] < 1) stop("t needs to be a factor or integer vector of time-periods between 1 and nrow(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i;
        } else {
          stop("Repeated values in timevar");
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
            for(int i = np; i != l; ++i) outp[omap[i]] = FUN(column[omap[i]], column[omap[i - np]]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1)-1;
              for(int i = l-1; i != start; --i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
            }
            for(int i = end; i--; ) outp[omap[i]] = fill;
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
            for(int i = l+np; i--; ) outp[omap[i]] = FUN(column[omap[i]], column[omap[i - np]]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int final = l+np*(k+1);
              for(int i = 0; i != final; ++i) outp[omap[i]] = FUN(outp[omap[i]], outp[omap[i - np]]);
            }
            for(int i = end; i != l; ++i) outp[omap[i]] = fill;
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
    IntegerVector gsv = Rf_isNull(gs) ? IntegerVector(ng) : as<IntegerVector>(gs); // no_init_vector(ng);
    if(Rf_isNull(t)) { // Ordered data
      if(maxdiff != 1) {
        if(Rf_isNull(gs)) {
          // gsv = IntegerVector(ng);
          // std::fill(gsv.begin(), gsv.end(), 0);
          for(int i = 0; i != l; ++i) ++gsv[g[i]-1];
        } else {
          // gsv = gs;
          if(ng != gsv.size()) stop("ng must match length(gs)");
        }
      }
      // int seen[ngp], memsize = sizeof(int)*(ngp);
      for(int j = 0; j != col; ++j) {
        NumericMatrix::ConstColumn column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
                if(seen[g[i]] == gsv[g[i]-1]-start) outp[i] = fill;
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
                    if(seen[g[i]] == gsv[g[i]-1]-start) outtemp[i] = fill;
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
                if(seen[g[i]] == gsv[g[i]-1]+start) outp[i] = fill;
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
                    if(seen[g[i]] == gsv[g[i]-1]+start) outtemp[i] = fill;
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
      if(l != ord.size()) stop("nrow(x) must match length(t)");
      IntegerVector min(ngp, INT_MAX);
      IntegerVector ord2 = no_init_vector(l);
      if(Rf_isNull(gs)) {
        // gsv = IntegerVector(ng);
        // std::fill(gsv.begin(), gsv.end(), 0);
        for(int i = 0; i != l; ++i) {
          ++gsv[g[i]-1];
          if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
        }
      } else {
        // gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != l; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
      }
      IntegerVector omap(l), cgs = no_init_vector(ngp), index = no_init_vector(l);
      // int cgs[ngp], seen[ngp], index[l], memsize = sizeof(int)*(ngp);
      cgs[1] = 0;
      for(int i = 1; i != ng; ++i) {
        cgs[i+1] = cgs[i] + gsv[i-1]; // or get "starts from forderv"
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?
      }
      if(min[ng] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?

      for(int i = 0; i != l; ++i) {
        ord2[i] = ord[i] - min[g[i]];
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        index[i] = cgs[g[i]]+ord2[i]; // index ??
        if(omap[index[i]] == 0) omap[index[i]] = i;
        else stop("Repeated values of timevar within one or more groups");
      }
      for(int j = 0; j != col; ++j) {
        NumericMatrix::ConstColumn column = x( _ , j);
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
            for(int i = 0; i != l; ++i) {
              if(ord2[i] >= np) {
                outp[i] = FUN(column[i], column[omap[index[i]-np]]);
              } else {
                outp[i] = fill;
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = l; i--; ) {
                if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outp[omap[i]] = fill;
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
                    if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outtemp[omap[i]] = fill;
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
            for(int i = 0; i != l; ++i) {
              if(ord2[i] < gsv[g[i]-1]+np) {
                outp[i] = FUN(column[i], column[omap[index[i]-np]]);
              } else {
                outp[i] = fill;
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = 0; i != l; ++i) {
                if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outp[omap[i]] = fill;
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
                    if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outtemp[omap[i]] = fill;
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
  //   if(ns*ds == 1) DUPLICATE_ATTRIB(out, x);
  //   // else rownames(out) = rownames(x); // redundant !!
  // }

  DUPLICATE_ATTRIB(out, x);
  if(ncol != col) out.attr("dim") = Dimension(l, ncol);
  if(names) {
    out.attr("dimnames") = List::create(rownames(x), colnam); // colnames(out) = colnam; also deletes rownames !
  } else if(ncol != col) {
    out.attr("dimnames") = R_NilValue;
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix fdiffgrowthmCpp(const NumericMatrix& x, const IntegerVector& n = 1, const IntegerVector& diff = 1,
                              double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0,
                              const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue,
                              int ret = 1, double rho = 1, bool names = true) {
  std::string stub;
  switch (ret)
  {       // [rho] or [&rho] ?  // https://stackoverflow.com/questions/30217956/error-variable-cannot-be-implicitly-captured-because-no-default-capture-mode-h
  case 1:
    if(names) stub = (rho == 1) ? "D" : "QD"; // QD for quasi-differences !
    return fdiffgrowthmCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return y-rho*x; }); // return y-x; same efficiency as return y-rho*x; when rho = 1 -> smart compiler !, and reduced file size !!
  case 2:
    if(rho == 1) goto fastld;
    if(names) stub = "QDlog";
    return fdiffgrowthmCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return std::log(y)-rho*std::log(x); }); // log(y*(1/(rho*x))) gives log(y) - log(rho*x), but we want log(y) - rho*log(x)
  case 3:
    if(names) stub = "G";
    return fdiffgrowthmCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return (y-x)*(rho/x); }); // same speed as fixing 100 !
  case 4:
    fastld:
    if(names) stub = "Dlog";
    return fdiffgrowthmCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return rho*std::log(y*(1/x)); });
  default: stop("Unknown return option!");
  }
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
  CharacterVector na = names ? coln_check(x.attr("names")) : NA_STRING;
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
      int os = ord.size(); // omap[os];
      IntegerVector omap = no_init_vector(os);
      LogicalVector ocheck(os, true);
      for(int i = 0; i != os; ++i) {
        if(ord[i] > os || ord[i] < 1) stop("t needs to be a factor or integer vector of time-periods between 1 and nrow(x)");
        if(ocheck[ord[i]-1]) {
          ocheck[ord[i]-1] = false;
          omap[ord[i]-1] = i;
        } else {
          stop("Repeated values in timevar");
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
            for(int i = np; i != os; ++i) outjp[omap[i]] = FUN(column[omap[i]], column[omap[i - np]]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1)-1;
              for(int i = os-1; i != start; --i) outjp[omap[i]] = FUN(outjp[omap[i]], outjp[omap[i - np]]);
            }
            for(int i = end; i--; ) outjp[omap[i]] = fill;
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
            for(int i = os+np; i--; ) outjp[omap[i]] = FUN(column[omap[i]], column[omap[i - np]]);
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int final = os+np*(k+1);
              for(int i = 0; i != final; ++i) outjp[omap[i]] = FUN(outjp[omap[i]], outjp[omap[i - np]]);
            }
            for(int i = end; i != os; ++i) outjp[omap[i]] = fill;
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
    IntegerVector gsv = Rf_isNull(gs) ? IntegerVector(ng) : as<IntegerVector>(gs); // no_init_vector(ng);
    if(Rf_isNull(t)) { // Ordered data
      if(maxdiff != 1) {
        if(Rf_isNull(gs)) {
          // gsv = IntegerVector(ng);
          // std::fill(gsv.begin(), gsv.end(), 0);
          for(int i = 0; i != gss; ++i) ++gsv[g[i]-1];
        } else {
          // gsv = gs;
          if(ng != gsv.size()) stop("ng must match length(gs)");
        }
      }
      // int seen[ngp], memsize = sizeof(int)*(ngp);
      for(int j = 0; j != l; ++j) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("nrow(x) must match length(g)");
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
                if(seen[g[i]] == gsv[g[i]-1]-start) outjp[i] = fill;
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
                    if(seen[g[i]] == gsv[g[i]-1]-start) outtemp[i] = fill;
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
                if(seen[g[i]] == gsv[g[i]-1]+start) outjp[i] = fill;
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
                    if(seen[g[i]] == gsv[g[i]-1]+start) outtemp[i] = fill;
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
      if(gss != ord.size()) stop("length(g) must match length(t)");
      IntegerVector min(ngp, INT_MAX);
      IntegerVector ord2 = no_init_vector(gss);
      if(Rf_isNull(gs)) {
        // gsv = IntegerVector(ng);
        // std::fill(gsv.begin(), gsv.end(), 0);
        for(int i = 0; i != gss; ++i) {
          ++gsv[g[i]-1];
          if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
        }
      } else {
        // gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != gss; ++i) if(ord[i] < min[g[i]]) min[g[i]] = ord[i];
      }
      IntegerVector omap(gss), cgs = no_init_vector(ngp), index = no_init_vector(gss);
      // int cgs[ngp], seen[ngp], index[gss], memsize = sizeof(int)*(ngp);
      cgs[1] = 0;
      for(int i = 1; i != ng; ++i) {
        cgs[i+1] = cgs[i] + gsv[i-1]; // or get "starts from forderv"
        if(min[i] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?
      }
      if(min[ng] == NA_INTEGER) stop("Timevar contains missing values"); // Fastest here ?

      for(int i = 0; i != gss; ++i) {
        ord2[i] = ord[i] - min[g[i]];
        if(ord2[i] >= gsv[g[i]-1]) stop("Gaps in timevar within one or more groups");
        index[i] = cgs[g[i]]+ord2[i]; // index ?
        if(omap[index[i]] == 0) omap[index[i]] = i;
        else stop("Repeated values of timevar within one or more groups");
      }
      for(int j = 0; j != l; ++j) {
        NumericVector column = x[j];
        if(gss != column.size()) stop("nrow(x) must match length(g)");
        for(int p = 0; p != ns; ++p) {
          int np = n[p];
          if(absn[p]*maxdiff > ags) warning("abs(n * diff) exceeds average group-size (%i). This could also be a result of unused factor levels. See #25.", ags);
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
            for(int i = 0; i != gss; ++i) {
              if(ord2[i] >= np) {
                outjp[i] = FUN(column[i], column[omap[index[i]-np]]);
              } else {
                outjp[i] = fill;
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = gss; i--; ) {
                if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outjp[omap[i]] = fill;
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
                    if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]-start) outtemp[omap[i]] = fill;
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
            for(int i = 0; i != gss; ++i) {
              if(ord2[i] < gsv[g[i]-1]+np) {
                outjp[i] = FUN(column[i], column[omap[index[i]-np]]);
              } else {
                outjp[i] = fill;
              }
            }
            if(d1 > 1) for(int k = 1; k != d1; ++k) {
              int start = np*(k+1);
              std::vector<int> seen(ngp); // memset(seen, 0, memsize);
              for(int i = 0; i != gss; ++i) {
                if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outjp[omap[i]] = fill;
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
                    if(seen[g[omap[i]]] == gsv[g[omap[i]]-1]+start) outtemp[omap[i]] = fill;
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
  DUPLICATE_ATTRIB(out, x);
  if(names) { // best way to code this ?
    out.attr("names") = nam;
  } else if(ncol != l) {
    out.attr("names") = R_NilValue;
  }
  return out;
}

// [[Rcpp::export]]
List fdiffgrowthlCpp(const List& x, const IntegerVector& n = 1, const IntegerVector& diff = 1,
                     double fill = NA_REAL, int ng = 0, const IntegerVector& g = 0,
                     const SEXP& gs = R_NilValue, const SEXP& t = R_NilValue,
                     int ret = 1, double rho = 1, bool names = true) {

  std::string stub;
  switch (ret)
  {       // [rho] or [&rho] ?  // https://stackoverflow.com/questions/30217956/error-variable-cannot-be-implicitly-captured-because-no-default-capture-mode-h
  case 1:
    if(names) stub = (rho == 1) ? "D" : "QD"; // QD for quasi-differences !
    return fdiffgrowthlCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return y-rho*x; }); // return y-x; same efficiency as return y-rho*x; when rho = 1 -> smart compiler !, and reduced file size !!
  case 2:
    if(rho == 1) goto fastld;
    if(names) stub = "QDlog";
    return fdiffgrowthlCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return std::log(y)-rho*std::log(x); }); // log(y*(1/(rho*x))) gives log(y) - log(rho*x), but we want log(y) - rho*log(x)
  case 3:
    if(names) stub = "G";
    return fdiffgrowthlCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return (y-x)*(rho/x); }); // same speed as fixing 100 !
  case 4:
    fastld:
    if(names) stub = "Dlog";
    return fdiffgrowthlCppImpl(x, n, diff, fill, ng, g, gs, t, stub, names, [rho](double y, double x) { return rho*std::log(y*(1/x)); });
  default: stop("Unknown return option!");
  }
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

