// [[Rcpp::plugins(cpp11)]]
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp ;

inline bool isnan2(double x) {
  return x != x;
}

// Todo: Numeric Stability??

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
    
  } else { // with groups !!!
    List med(l);
    std::vector<std::vector<double> > gmap(ng); 
    int ngp = ng+1, gcount[ngp], memsize = sizeof(int)*ngp, gss = g.size(); 
    if(Rf_isNull(gs)) {
      memset(gcount, 0, memsize);
      for(int i = 0; i != gss; ++i) ++gcount[g[i]];
      for(int i = 0; i != ng; ++i) gmap[i] = std::vector<double> (gcount[i+1]); 
    } else {
      IntegerVector gsv = gs;
      if(ng != gsv.size()) stop("ng must match length(gs)");
      for(int i = 0; i != ng; ++i) gmap[i] = std::vector<double> (gsv[i]); 
    }
    if(narm) {
      for(int j = l; j--; ) { 
        NumericVector column = x[j];
        NumericVector medj(ng, NA_REAL);
        if(gss != column.size()) stop("length(g) must match nrow(x)");
        memset(gcount, 0, memsize);
        for(int i = 0; i != gss; ++i) if(!std::isnan(column[i])) gmap[g[i]-1][gcount[g[i]]++] = column[i]; // good ??
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
        NumericVector medj(ng); //  = no_init_vector // no init numerically instable !!
        if(gss != column.size()) stop("length(g) must match nrow(x)");
        memset(gcount, 0, memsize);
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
          if(n%2 == 0){ // speed up by saving iterators ?? -> Yes !!
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
  med.attr("row.names") = NumericVector::create(NA_REAL, -ng);
  return med;
  }
} 

// A grouping option with ngp = ng+1 trick -> slightly faster !!
// if(ng != gs.size()) stop("ng must match length(gs)");
// int ngp = ng+1;
// std::vector<std::vector<double> > out(ngp); // NumericVector
// for(int i = 0; i != ng; ++i) out[i+1] = std::vector<double>(gs[i]); // no_init_vector(gs[i]);
// int gcount[ngp], memsize = sizeof(int)*ngp, row = 0, gss = g.size(); // good !! integervector gives error !!
// 
// if(narm) {
//   for(int j = l; j--; ) { 
//     NumericVector column = x[j];
//     memset(gcount, 0, memsize);
//     row = column.size();
//     if(row != gss) stop("length(g) must match nrow(x)");
//     NumericVector outj = no_init_vector(ng);
//     for(int i = 0; i != row; ++i) if(!std::isnan(column[i])) out[g[i]][gcount[g[i]]++] = column[i]; // good ??
//     for(int i = 1; i != ngp; ++i) {
//       int n = gcount[i], middle = n/2-1;
//       auto begin = out[i].begin(), end = out[i].end()-(gs[i-1]-n);
//       if(n%2==0){
//         std::nth_element(begin,begin+middle,end);
//         outj[i-1] = (out[i][middle]+*(std::min_element(begin+middle+1,end)))/2.0;
//       } else{
//         std::nth_element(begin,begin+middle+1,end);
//         outj[i-1] = out[i][middle+1];
//       }
//     }
//     SHALLOW_DUPLICATE_ATTRIB(outj, column);
//     med[j] = outj;  
//   }
// } else {
//   for(int j = l; j--; ) { 
//     NumericVector column = x[j];
//     memset(gcount, 0, memsize);
//     row = column.size();
//     if(row != gss) stop("length(g) must match nrow(x)");
//     NumericVector outj = no_init_vector(ng);
//     int ngs = 0;
//     for(int i = 0; i != row; ++i) {
//       if(std::isnan(column[i])) {
//         if(!std::isnan(outj[g[i]-1])) {
//           outj[g[i]-1] = NA_REAL;
//           ++ngs;
//           if(ngs == ng) break;
//         }
//       } else {
//         out[g[i]][gcount[g[i]]++] = column[i];
//       }
//     }
//     for(int i = 1; i != ngp; ++i) {
//       if(std::isnan(outj[i-1])) continue;
//       int n = gs[i-1], middle = n/2-1; 
//       if(n%2==0){
//         std::nth_element(out[i].begin(),out[i].begin()+middle,out[i].end());
//         outj[i-1] = (out[i][middle]+*(std::min_element(out[i].begin()+middle+1,out[i].end())))/2.0;
//       } else{
//         std::nth_element(out[i].begin(),out[i].begin()+middle+1,out[i].end());
//         outj[i-1] = out[i][middle+1];
//       }
//     }
//     SHALLOW_DUPLICATE_ATTRIB(outj, column);
//     med[j] = outj;  
//   }
// }
// DUPLICATE_ATTRIB(med, x); 
// if(x.hasAttribute("row.names")) med.attr("row.names") = seq_len(ng);
// }

// Another Grouping option -> saving indices and looping through groups: But slower !!
// if(ng != gs.size()) stop("ng must match length(gs)");
// std::vector<std::vector<double> > out(ng); // NumericVector
// for(int i = 0; i != ng; ++i) out[i] = std::vector<double>(gs[i]); // no_init_vector(gs[i]);
// int gcount[ng], memsize = sizeof(int)*ng, row = 0, gss = g.size(); // good !! integervector gives error !!
// memset(gcount, 0, memsize);
// for(int i = 0; i != gss; ++i) out[g[i]-1][gcount[g[i]-1]++] = i; // This is definitely slower !!
// 
// for(int j = l; j--; ) { 
//   NumericVector column = x[j];
//   memset(gcount, 0, memsize);
//   row = column.size();
//   if(row != gss) stop("length(g) must match nrow(x)");
//   NumericVector outj = no_init_vector(ng);
//   for(int i = 0; i != ng; ++i) {
//     int gsi = gs[i], count = 0;
//     NumericVector outi = no_init_vector(gsi);
//     for(int k = 0; k != gsi; ++k) if(!std::isnan(column[out[i][k]])) outi[count++] = column[out[i][k]];
//     int middle = count/2-1;
//     auto begin = outi.begin(), end = outi.end()-(gsi-count);
//     if(count%2==0){
//       std::nth_element(begin,begin+middle,end);
//       outj[i] = (outi[middle]+*(std::min_element(begin+middle+1,end)))/2.0;
//     } else{
//       std::nth_element(begin,begin+middle+1,end);
//       outj[i] = outi[middle+1];
//     }
//   }
//   SHALLOW_DUPLICATE_ATTRIB(outj, column);
//   med[j] = outj;  
// }