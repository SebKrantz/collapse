// [[Rcpp::plugins(cpp11)]]
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp ;

inline bool isnan2(double x) {
  return x != x;
}

// Todo: Numeric Stability??
// https://www.quantstart.com/articles/Passing-By-Reference-To-Const-in-C
// https://www.cs.fsu.edu/~myers/c++/notes/references.html

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
    int ngp = ng+1, gcount[ngp], memsize = sizeof(int)*ngp; 
    if(Rf_isNull(gs)) {
      memset(gcount, 0, memsize);
      for(int i = 0; i != l; ++i) ++gcount[g[i]];
      for(int i = 0; i != ng; ++i) gmap[i] = std::vector<double> (gcount[i+1]); 
    } else {
      IntegerVector gsv = gs;
      if(ng != gsv.size()) stop("ng must match length(gs)");
      for(int i = 0; i != ng; ++i) gmap[i] = std::vector<double>(gsv[i]); 
    }
    if(narm) {
      NumericMatrix med = no_init_matrix(ng, col);  
      std::fill(med.begin(), med.end(), NA_REAL);
      for(int j = col; j--; ) { 
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column medj = med( _ , j);
        memset(gcount, 0, memsize);
        for(int i = 0; i != l; ++i) if(!std::isnan(column[i])) gmap[g[i]-1][gcount[g[i]]++] = column[i]; // good ??
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
      colnames(med) = colnames(x); // med.attr("dimnames") = List::create(R_NilValue, colnames(x)); // needed ?? 
      return med;
    } else {
      NumericMatrix med(ng, col); // no init numerically unstable !!
      for(int j = col; j--; ) { 
        NumericMatrix::ConstColumn column = x( _ , j);
        NumericMatrix::Column medj = med( _ , j);
        memset(gcount, 0, memsize);
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
          if(n%2 == 0){ // speed up by saving iterators ?? -> Yes !!
            std::nth_element(begin, mid, end);
            medj[i] = (*(mid) + *(std::min_element(mid+1, end)))/2.0;
          } else{
            std::nth_element(begin, mid+1, end);
            medj[i] = *(mid+1);
          }
        }
      }
      colnames(med) = colnames(x); // med.attr("dimnames") = List::create(R_NilValue, colnames(x)); // needed ?? 
      return med;
    }
  }
} 
