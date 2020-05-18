// [[Rcpp::plugins(cpp11)]]
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp ;

inline bool isnan2(double x) {
  return x != x;
}

// https://www.quantstart.com/articles/Passing-By-Reference-To-Const-in-C
// https://www.cs.fsu.edu/~myers/c++/notes/references.html
// https://stackoverflow.com/questions/24112893/within-c-functions-how-are-rcpp-objects-passed-to-other-functions-by-referen
// https://stackoverflow.com/questions/48363076/declare-a-variable-as-a-reference-in-rcpp?noredirect=1&lq=1
// See these articles, you cannot declare them as const references and then apply std::nth_element which modifies the data. you need to make a copy !!

//[[Rcpp::export]]
NumericVector fmedianCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0,
                         const SEXP& gs = R_NilValue, bool narm = true) {
  int l = x.size();

  if(ng == 0) {
    NumericVector med(1);
    if(narm) { // xd = xd[!is_na(xd)]; // super slow !
      NumericVector xd = no_init_vector(l);
      // NumericVector::iterator pend = std::remove_if(xd.begin(), xd.end(), isnan2);
      // auto pend = std::remove_if(xd.begin(), xd.end(), isnan2); // need remove_copy_if ?
      auto pend = std::remove_copy_if(x.begin(), x.end(), xd.begin(), isnan2); // faster !
      int sz = pend - xd.begin(), middle = sz/2-1; // std::distance(xd.begin(), pend)
      if(sz == 0) return NumericVector::create(x[0]);
      if(sz%2 == 0){
        std::nth_element(xd.begin(), xd.begin()+middle, pend);
        med = (xd[middle] + *(std::min_element(xd.begin()+middle+1, pend)))/2.0;
      } else {
        std::nth_element(xd.begin(), xd.begin()+middle+1, pend);
        med = xd[middle+1];
      }
    } else {
      for(int i = 0; i != l; ++i) if(std::isnan(x[i])) return NumericVector::create(x[i]);
      NumericVector xd = Rf_duplicate(x);
      int middle = l/2-1;
      if(l%2 == 0){
        std::nth_element(xd.begin(), xd.begin()+middle, xd.end());
        med = (xd[middle] + *(std::min_element(xd.begin()+middle+1, xd.end())))/2.0;
      } else{
        std::nth_element(xd.begin(), xd.begin()+middle+1, xd.end());
        med = xd[middle+1];
      }
    }
    // DUPLICATE_ATTRIB(med, x); // or keep double use wrap() ?
    return med;
  } else { // with groups
    if(l != g.size()) stop("length(g) must match length(x)");
    std::vector<std::vector<double> > gmap(ng);
    // IntegerVector gcount(ng); // array better ? 2.30 sec for 1e7 with 1e6 groups
    int ngp = ng+1; // , gcount[ngp]; // Yes! 2.25 sec !
    // memset(gcount, 0, sizeof(int)*ngp);
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
      NumericVector med(ng, NA_REAL);
      for(int i = 0; i != l; ++i) if(!std::isnan(x[i])) gmap[g[i]-1][gcount[g[i]]++] = x[i]; // good ?
      for(int i = 0; i != ng; ++i) {
        if(gcount[i+1] != 0) {
          int n = gcount[i+1]; // fastest !
          auto begin = gmap[i].begin(), mid = begin + n/2-1, end = begin + n; // gmap[i].end()-(gmap[i].size()-n); // (gs[i]-n) // good ?
          if(n%2 == 0){
            std::nth_element(begin, mid, end);
            med[i] = (*(mid) + *(std::min_element(mid+1, end)))/2.0; // gmap[i][middle] or *(begin+middle) ? -> second
          } else{
            std::nth_element(begin, mid+1, end);
            med[i] = *(mid+1); // gmap[i][middle+1] or *(begin+middle+1) ? -> second
          }
        }
      }
      DUPLICATE_ATTRIB(med, x);
      return med;
    } else {
      NumericVector med(ng); //  = no_init_vector // no init numerically unstable (calling isnan on a not-initialized vector is a bad idea)
     // for(int i = 0; i != l; ++i) gmap[g[i]-1][gcount[g[i]]++] = xd[i]; // good ?
      int ngs = 0;
      for(int i = 0; i != l; ++i) {
        // if(std::isnan(med[g[i]-1])) continue; // A lot faster for NWDI ! (reading into 2D structure takes time, unlike fsum)
        if(std::isnan(x[i])) {
          if(!std::isnan(med[g[i]-1])) { // however, if there are no missing values, this is a lot faster 2.21 vs. 3.09 seconds on 1e7 with 1e6 groups
            med[g[i]-1] = NA_REAL;
            ++ngs;
            if(ngs == ng) break;
          }
        } else {
          gmap[g[i]-1][gcount[g[i]]++] = x[i];
        }
      }
      for(int i = 0; i != ng; ++i) {
        if(std::isnan(med[i])) continue;
        int n = gcount[i+1]; //  gs[i] // fastest !
        auto begin = gmap[i].begin(), mid = begin + n/2-1, end = begin + n;
        if(n%2 == 0){
          std::nth_element(begin, mid, end);
          med[i] = (*(mid) + *(std::min_element(mid+1, end)))/2.0;
        } else{
          std::nth_element(begin, mid+1, end);
          med[i] = *(mid+1);
        }
      }
      DUPLICATE_ATTRIB(med, x);
      return med;
    }
  }
}




//[[Rcpp::export]]
SEXP fmedianmCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0,
                 const SEXP& gs = R_NilValue, bool narm = true, bool drop = true) {
  int l = x.nrow(), col = x.ncol();

  if(ng == 0) {
    NumericVector med = no_init_vector(col);
    if(narm) {
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn colum = x( _ , j);
        NumericVector column = no_init_vector(l);
        auto pend = std::remove_copy_if(colum.begin(), colum.end(), column.begin(), isnan2);
        int sz = pend - column.begin(), middle = sz/2-1; // std::distance(x.begin(), pend)
        if(sz%2 == 0){
          std::nth_element(column.begin(), column.begin()+middle, pend);
          med[j] = (column[middle] + *(std::min_element(column.begin()+middle+1, pend)))/2.0;
        } else{
          std::nth_element(column.begin(), column.begin()+middle+1, pend);
          med[j] = column[middle+1];
        }
      }
    } else {
      int middle = l/2-1;
      if(l%2 == 0) {
        for(int j = col; j--; ) {
          {
            NumericVector colum = x( _ , j);
            for(int i = 0; i != l; ++i) {
              if(std::isnan(colum[i])) {
                med[j] = colum[i];
                goto endloop;
              }
            }
            NumericVector column = Rf_duplicate(colum);
            std::nth_element(column.begin(), column.begin()+middle, column.end());
            med[j] = (column[middle] + *(std::min_element(column.begin()+middle+1, column.end())))/2.0;
          }
          endloop:;
        }
      } else {
        for(int j = col; j--; ) {
          {
            NumericVector colum = x( _ , j);
            for(int i = 0; i != l; ++i) {
              if(std::isnan(colum[i])) {
                med[j] = colum[i];
                goto endloop2;
              }
            }
            NumericVector column = Rf_duplicate(colum);
            std::nth_element(column.begin(), column.begin()+middle+1, column.end());
            med[j] = column[middle+1];
          }
          endloop2:;
        }
      }
    }
    if(drop) med.attr("names") = colnames(x);
    else {
      med.attr("dim") = Dimension(1, col);
      colnames(med) = colnames(x);
    }
    return med;
  } else { // with groups
    if(l != g.size()) stop("length(g) must match nrow(x)");
    std::vector<std::vector<double> > gmap(ng);
    int ngp = ng+1; // , gcount[ngp], memsize = sizeof(int)*ngp;
    std::vector<int> gcount(ngp);
    if(Rf_isNull(gs)) {
      // memset(gcount, 0, memsize);
      for(int i = 0; i != l; ++i) ++gcount[g[i]];
      for(int i = 0; i != ng; ++i) {
        if(gcount[i+1] == 0) stop("group size of 0 encountered");
        gmap[i] = std::vector<double> (gcount[i+1]);
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
      NumericMatrix med = no_init_matrix(ng, col);
      std::fill(med.begin(), med.end(), NA_REAL);
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column medj = med( _ , j);
        gcount.assign(ngp, 0); // memset(gcount, 0, memsize);
        for(int i = 0; i != l; ++i) if(!std::isnan(column[i])) gmap[g[i]-1][gcount[g[i]]++] = column[i]; // good ?
        for(int i = 0; i != ng; ++i) {
          if(gcount[i+1] != 0) {
            int n = gcount[i+1];
            auto begin = gmap[i].begin(), mid = begin + n/2-1, end = begin + n; // gmap[i].end()-(gs[i]-n);
            if(n%2 == 0){
              std::nth_element(begin, mid, end);
              medj[i] = (*(mid) + *(std::min_element(mid+1, end)))/2.0;
            } else{
              std::nth_element(begin, mid+1, end);
              medj[i] = *(mid+1);
            }
          }
        }
      }
      colnames(med) = colnames(x); // med.attr("dimnames") = List::create(R_NilValue, colnames(x)); // needed ?
      return med;
    } else {
      NumericMatrix med(ng, col); // no init numerically unstable
      for(int j = col; j--; ) {
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column medj = med( _ , j);
        gcount.assign(ngp, 0); // memset(gcount, 0, memsize);
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(column[i])) {
            if(!std::isnan(medj[g[i]-1])) {
              medj[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            gmap[g[i]-1][gcount[g[i]]++] = column[i];
          }
        }
        for(int i = 0; i != ng; ++i) {
          if(std::isnan(medj[i])) continue;
          int n = gcount[i+1];
          auto begin = gmap[i].begin(), mid = begin + n/2-1, end = begin + n;
          if(n%2 == 0){ // speed up by saving iterators ?? -> Yes
            std::nth_element(begin, mid, end);
            medj[i] = (*(mid) + *(std::min_element(mid+1, end)))/2.0;
          } else{
            std::nth_element(begin, mid+1, end);
            medj[i] = *(mid+1);
          }
        }
      }
      colnames(med) = colnames(x); // med.attr("dimnames") = List::create(R_NilValue, colnames(x)); // needed ?
      return med;
    }
  }
}




//[[Rcpp::export]]
SEXP fmedianlCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
                 const SEXP& gs = R_NilValue, bool narm = true, bool drop = true) {
  int l = x.size();

  if(ng == 0) {
    NumericVector med = no_init_vector(l);
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector colum = x[j];
        NumericVector column = no_init_vector(colum.size());
        auto pend = std::remove_copy_if(colum.begin(), colum.end(), column.begin(), isnan2);
        int sz = pend - column.begin(), middle = sz/2-1; // std::distance(x.begin(), pend)
        if(sz%2 == 0){
          std::nth_element(column.begin(), column.begin()+middle, pend);
          med[j] = (column[middle] + *(std::min_element(column.begin()+middle+1, pend)))/2.0;
        } else{
          std::nth_element(column.begin(), column.begin()+middle+1, pend);
          med[j] = column[middle+1];
        }
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = Rf_duplicate(x[j]);
        int row = column.size(), middle = row/2-1;
        for(int i = 0; i != row; ++i) {
          if(std::isnan(column[i])) {
            med[j] = column[i];
            goto endloop;
          }
        }
        if(row%2 == 0){
          std::nth_element(column.begin(), column.begin()+middle, column.end());
          med[j] = (column[middle] + *(std::min_element(column.begin()+middle+1, column.end())))/2.0;
        } else{
          std::nth_element(column.begin(), column.begin()+middle+1, column.end());
          med[j] = column[middle+1];
        }
        endloop:;
      }
    }
    // if(drop) return(Rf_setAttrib(Rf_coerceVector(med, REALSXP),CharacterVector::create("names"),x.attr("names")));
    // else {
    //   // List ax = ATTRIB(x);
    //   // if(ax["row.names"] != R_NilValue) ax["row.names"] = IntegerVector(1,1);
    //   SHALLOW_DUPLICATE_ATTRIB(med, x);
    //   if(x.hasAttribute("row.names")) med.attr("row.names") = 1;
    //   // SET_ATTRIB(med, Rf_coerceVector(ax, LISTSXP));
    // }
    if(drop) {
      med.attr("names") = x.attr("names");
      return med;
    } else {
      List out(l);
      for(int j = l; j--; ) {
        out[j] = med[j];
        SHALLOW_DUPLICATE_ATTRIB(out[j], x[j]);
      }
      DUPLICATE_ATTRIB(out, x);
      out.attr("row.names") = 1;
      return out;
    }

  } else { // with groups
    List med(l);
    std::vector<std::vector<double> > gmap(ng);
    int ngp = ng+1, gss = g.size(); //  , gcount[ngp], memsize = sizeof(int)*ngp;
    std::vector<int> gcount(ngp);
    if(Rf_isNull(gs)) {
      // memset(gcount, 0, memsize);
      for(int i = 0; i != gss; ++i) ++gcount[g[i]];
      for(int i = 0; i != ng; ++i) {
        if(gcount[i+1] == 0) stop("group size of 0 encountered");
        gmap[i] = std::vector<double> (gcount[i+1]);
      }
    } else {
      IntegerVector gsv = gs;
      if(ng != gsv.size()) stop("ng must match length(gs)");
      for(int i = 0; i != ng; ++i) {
        if(gsv[i] == 0) stop("group size of 0 encountered");
        gmap[i] = std::vector<double> (gsv[i]);
      }
    }
    if(narm) {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        NumericVector medj(ng, NA_REAL);
        if(gss != column.size()) stop("length(g) must match nrow(x)");
        gcount.assign(ngp, 0); // memset(gcount, 0, memsize);
        for(int i = 0; i != gss; ++i) if(!std::isnan(column[i])) gmap[g[i]-1][gcount[g[i]]++] = column[i]; // good ?
        for(int i = 0; i != ng; ++i) {
          if(gcount[i+1] != 0) {
            int n = gcount[i+1];
            auto begin = gmap[i].begin(), mid = begin + n/2-1, end = begin + n; // gmap[i].end()-(gs[i]-n);
            if(n%2 == 0){
              std::nth_element(begin, mid, end);
              medj[i] = (*(mid) + *(std::min_element(mid+1, end)))/2.0;
            } else{
              std::nth_element(begin, mid+1, end);
              medj[i] = *(mid+1);
            }
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(medj, column);
        med[j] = medj;
      }
    } else {
      for(int j = l; j--; ) {
        NumericVector column = x[j];
        NumericVector medj(ng); //  = no_init_vector // no init numerically instable !
        if(gss != column.size()) stop("length(g) must match nrow(x)");
        gcount.assign(ngp, 0); // memset(gcount, 0, memsize);
        int ngs = 0;
        for(int i = 0; i != gss; ++i) {
          if(std::isnan(column[i])) {
            if(!std::isnan(medj[g[i]-1])) {
              medj[g[i]-1] = NA_REAL;
              ++ngs;
              if(ngs == ng) break;
            }
          } else {
            gmap[g[i]-1][gcount[g[i]]++] = column[i];
          }
        }
        for(int i = 0; i != ng; ++i) {
          if(std::isnan(medj[i])) continue;
          int n = gcount[i+1];
          auto begin = gmap[i].begin(), mid = begin + n/2-1, end = begin + n;
          if(n%2 == 0){ // speed up by saving iterators ? -> Yes !
            std::nth_element(begin, mid, end);
            medj[i] = (*(mid) + *(std::min_element(mid+1, end)))/2.0;
          } else{
            std::nth_element(begin, mid+1, end);
            medj[i] = *(mid+1);
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(medj, column);
        med[j] = medj;
      }
    }
    DUPLICATE_ATTRIB(med, x);
    med.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng);
    return med;
  }
}

