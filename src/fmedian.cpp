// [[Rcpp::plugins(cpp11)]]
#include <numeric>
#include <Rcpp.h>
using namespace Rcpp ;

inline bool isnan2(double x) {
  return x != x;
}

// Todo: Numeric Stability?? Note: THis changes X !! need to make a copy !!
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
    if(narm) { // xd = xd[!is_na(xd)]; // super slow !!
      NumericVector xd = no_init_vector(l);
      // NumericVector::iterator pend = std::remove_if(xd.begin(), xd.end(), isnan2);
      // auto pend = std::remove_if(xd.begin(), xd.end(), isnan2); // need remove_copy_if ??
      auto pend = std::remove_copy_if(x.begin(), x.end(), xd.begin(), isnan2); // faster !!
      int sz = pend - xd.begin(), middle = sz/2-1; // std::distance(xd.begin(), pend)
      if(sz%2 == 0){
        std::nth_element(xd.begin(), xd.begin()+middle, pend);
        med = (xd[middle] + *(std::min_element(xd.begin()+middle+1, pend)))/2.0;
      } else{
        std::nth_element(xd.begin(), xd.begin()+middle+1, pend);
        med = xd[middle+1];
      }
    } else {
      for(int i = 0; i != l; ++i) if(std::isnan(x[i])) return(NumericVector::create(x[i]));
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
    // DUPLICATE_ATTRIB(med, x); // or keep double use wrap() ??
    return med;
  } else { // with groups !!!
    if(l != g.size()) stop("length(g) must match length(x)");
    std::vector<std::vector<double> > gmap(ng); 
    //IntegerVector gcount(ng); // array better ?? 2.30 sec for 1e7 with 1e6 groups
    int ngp = ng+1, gcount[ngp]; // Yes! 2.25 sec !!
    memset(gcount, 0, sizeof(int)*ngp);
    if(Rf_isNull(gs)) {
      for(int i = 0; i != l; ++i) ++gcount[g[i]];
      for(int i = 0; i != ng; ++i) {
        gmap[i] = std::vector<double> (gcount[i+1]); 
        gcount[i+1] = 0;
      }
    } else {
      IntegerVector gsv = gs;
      if(ng != gsv.size()) stop("ng must match length(gs)");
      for(int i = 0; i != ng; ++i) gmap[i] = std::vector<double>(gsv[i]); 
    }
    
    if(narm) {
      NumericVector med(ng, NA_REAL); 
      for(int i = 0; i != l; ++i) if(!std::isnan(x[i])) gmap[g[i]-1][gcount[g[i]]++] = x[i]; // good ??
      for(int i = 0; i != ng; ++i) {
        if(gcount[i+1] != 0) {
          int n = gcount[i+1]; // fastest !!
          auto begin = gmap[i].begin(), mid = begin + n/2-1, end = begin + n; // gmap[i].end()-(gmap[i].size()-n); // (gs[i]-n) // good ??
          if(n%2 == 0){
            std::nth_element(begin, mid, end);
            med[i] = (*(mid) + *(std::min_element(mid+1, end)))/2.0; // gmap[i][middle] or *(begin+middle) ? -> second !!
          } else{
            std::nth_element(begin, mid+1, end);
            med[i] = *(mid+1); // gmap[i][middle+1] or *(begin+middle+1) ? -> second !!
          }
        }
      }
      DUPLICATE_ATTRIB(med, x);
      return med;
    } else {
      NumericVector med(ng); //  = no_init_vector // no init numerically unstable !!(calling isnan on a not-initialized vector is a bad idea)
     // for(int i = 0; i != l; ++i) gmap[g[i]-1][gcount[g[i]]++] = xd[i]; // good ??
      int ngs = 0;
      for(int i = 0; i != l; ++i) {
        // if(std::isnan(med[g[i]-1])) continue; // A lot faster for NWDI !! (reading into 2D structure takes time, unlike fsum !!)
        if(std::isnan(x[i])) {
          if(!std::isnan(med[g[i]-1])) { // however, if there are no missing values, this is a lot faster 2.21 vs. 3.09 seconds on 1e7 with 1e6 groups !!
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
        int n = gcount[i+1]; //  gs[i] // fastest !!
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


// //[[Rcpp::export]] // original: but changes underlying data !!
// long double med(NumericVector x){
//   
//   long double F;
//   int sz=x.size(),middle=sz/2-1;
//   if(sz%2==0){
//     std::nth_element(x.begin(),x.begin()+middle,x.end());
//     F=(x(middle)+*(std::min_element(x.begin()+middle+1,x.end())))/2.0;
//   }else{
//     std::nth_element(x.begin(),x.begin()+middle+1,x.end());
//     F=x(middle+1);
//   }
//   return F;
// }


// Previous Version:: using temp !!
// //[[Rcpp::export]]
// NumericVector fmedianCpp(NumericVector x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, bool narm = true) {
//   int l = x.size();
//   
//     if(ng == 0) {
//         NumericVector med(1);
//         if(narm) { // x = x[!is_na(x)]; // super slow !!
//           // NumericVector::iterator pend = std::remove_if(x.begin(), x.end(), isnan2);
//           auto pend = std::remove_if(x.begin(), x.end(), isnan2);
//           int sz = pend - x.begin(), middle = sz/2-1; // std::distance(x.begin(), pend)
//           if(sz%2==0){
//             std::nth_element(x.begin(),x.begin()+middle,pend);
//             med=(x(middle)+*(std::min_element(x.begin()+middle+1,pend)))/2.0;
//           } else{
//             std::nth_element(x.begin(),x.begin()+middle+1,pend);
//             med=x(middle+1);
//           }
//         } else {
//           int middle = l/2-1;
//           if(l%2==0){
//             std::nth_element(x.begin(),x.begin()+middle,x.end());
//             med=(x(middle)+*(std::min_element(x.begin()+middle+1,x.end())))/2.0;
//           } else{
//             std::nth_element(x.begin(),x.begin()+middle+1,x.end());
//             med=x(middle+1);
//           }
//         }
//       SHALLOW_DUPLICATE_ATTRIB(med, x); // or keep double use wrap() ?? 
//       return med;
//     } else { // with groups !!!
//       if(l != g.size()) stop("length(g) must match length(x)");
//       if(ng != gs.size()) stop("ng must match length(gs)");
//       std::vector<std::vector<double> > out(ng); // NumericVector
//       IntegerVector gcount(ng);
//       NumericVector med = no_init_vector(ng);
//       for(int i = 0; i != ng; ++i) out[i] = std::vector<double>(gs[i]); // no_init_vector(gs[i]); 
//       if(narm) {
//         for(int i = 0; i != l; ++i) if(!std::isnan(x[i])) out[g[i]-1][gcount[g[i]-1]++] = x[i]; // good ??
//         for(int i = 0; i != ng; ++i) {
//           std::vector<double> temp = out[i]; // faster using temp = out[i] ??
//           int n = gcount[i], middle = n/2-1;
//           auto end = temp.end()-(gs[i]-n);  
//           if(n%2==0){
//             std::nth_element(temp.begin(),temp.begin()+middle,end);
//             med[i] = (temp[middle]+*(std::min_element(temp.begin()+middle+1,end)))/2.0;
//           } else{
//             std::nth_element(temp.begin(),temp.begin()+middle+1,end);
//             med[i] = temp[middle+1];
//           }
//         }
//       } else {
//         for(int i = 0; i != l; ++i) out[g[i]-1][gcount[g[i]-1]++] = x[i]; // good ??
//         for(int i = 0; i != ng; ++i) {
//           std::vector<double> temp = out[i]; // faster using temp = out[i] ?? 
//           int n = gs[i], middle = n/2-1; 
//           if(n%2==0){
//             std::nth_element(temp.begin(),temp.begin()+middle,temp.end());
//             med[i] = (temp[middle]+*(std::min_element(temp.begin()+middle+1,temp.end())))/2.0;
//           } else{
//             std::nth_element(temp.begin(),temp.begin()+middle+1,temp.end());
//             med[i] = temp[middle+1];
//           }
//         }
//       }
//       SHALLOW_DUPLICATE_ATTRIB(med, x);
//       return med;
//     }
// } 



