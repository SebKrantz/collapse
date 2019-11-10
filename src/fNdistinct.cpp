// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp ;

// TODO: Perhaps redo everything with data pointers and 2d group indices (instead of filling the 2d structure every time !!): http://www.cplusplus.com/reference/vector/vector/data/
// https://stackoverflow.com/questions/1733143/converting-between-c-stdvector-and-c-array-without-copying?rq=1

// improve logical method ?? 

template <int RTYPE>
inline bool isnaNUM(typename Rcpp::traits::storage_type<RTYPE>::type x) { 
  return x != x;
}

template <int RTYPE>
inline bool isnaOTH(typename Rcpp::traits::storage_type<RTYPE>::type x) { 
  return x == Vector<RTYPE>::get_na();
}

template <int RTYPE>
IntegerVector fNdistinct2Impl(const Vector<RTYPE>& x, int ng, const IntegerVector& g, 
                              const SEXP& gs, bool narm) {
  int l = x.size();
  typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;
  
  if(ng == 0) {
    sugar::IndexHash<RTYPE> hash(x);
     if(narm) {
      for(int i = 0; i != l; ++i) {
        if(isnanT(x[i])) continue;
        unsigned int addr = hash.get_addr(x[i]);
        while(hash.data[addr] && hash.not_equal(x[hash.data[addr] - 1], x[i])) {
          ++addr;
          if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
        }
        if(!hash.data[addr]) {
          hash.data[addr] = i+1;
          ++hash.size_;
        }
      }
    } else {
     hash.fill();
    }
    return IntegerVector::create(hash.size_);
  } else {
    if(l != g.size()) stop("length(g) must match length(x)");
    std::vector<std::vector<storage_t> > gmap(ng);
    IntegerVector out(ng);
    if(Rf_isNull(gs)) {
      for(int i = 0; i != l; ++i) ++out[g[i]-1];
      for(int i = 0; i != ng; ++i) {
        gmap[i] = std::vector<storage_t> (out[i]); 
        out[i] = 0;
      }
      // memset(out, 0, sizeof(int)*ng); // Stable ??? -> Nope, gives error!!
    } else {
      IntegerVector gsv = gs;
      if(ng != gsv.size()) stop("ng must match length(gs)");
      for(int i = 0; i != ng; ++i) gmap[i] = std::vector<storage_t> (gsv[i]); 
    }
    for(int i = 0; i != l; ++i) gmap[g[i]-1][out[g[i]-1]++] = x[i];
    if(narm) {
      for(int gr = 0; gr != ng; ++gr) {
        const std::vector<storage_t>& temp = gmap[gr]; // good ?? // const Vector<RTYPE>& // wrap()
        sugar::IndexHash<RTYPE> hash(wrap(temp));
        for(int i = hash.n; i--; ) {
          if(isnanT(temp[i])) continue;
          unsigned int addr = hash.get_addr(temp[i]);
          while(hash.data[addr] && hash.not_equal(temp[hash.data[addr] - 1], temp[i])) {
            ++addr;
            if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
          }
          if(!hash.data[addr]) {
            hash.data[addr] = i+1;
            ++hash.size_;
          }
        }
        out[gr] = hash.size_;
      }
    } else {
      for(int gr = 0; gr != ng; ++gr) {
        sugar::IndexHash<RTYPE> hash(wrap(gmap[gr]));
        hash.fill();
        out[gr] = hash.size_;
      }
    }
    if(Rf_getAttrib(x, R_ClassSymbol) == R_NilValue) {
      SHALLOW_DUPLICATE_ATTRIB(out, x);
    } else {
      out.attr("label") = x.attr("label");
    }
    return out;
  }
}

template <> // No logical vector with sugar::IndexHash<RTYPE> !!!
IntegerVector fNdistinct2Impl(const Vector<LGLSXP>& x, int ng, const IntegerVector& g, 
                              const SEXP& gs, bool narm) {
  int l = x.size();
  
  if(ng == 0) { 
    int Nunique = 0;
    if(narm) {
      bool which = true;
      for(int i = 0; i != l; ++i) {
        if(x[i] == NA_LOGICAL) continue;
        if(x[i] == which) {
           Nunique = 1;
        } else {
          which = x[i];
          ++Nunique;
          if(Nunique == 2) break;
        }
      }
    } else {
      bool seen1 = true, seen2 = true, seen3 = true;
      for(int i = 0; i != l; ++i) { // better way?? 
        if(seen1 && x[i] == NA_LOGICAL) {
          ++Nunique;
          seen1 = false;
        } else if(seen2 && x[i] == true) {
          ++Nunique;
          seen2 = false;
        } else if(seen3 && x[i] == false) {
          ++Nunique;
          seen3 = false;
        }
        if(Nunique == 3) break;
      }
    }
    return IntegerVector::create(Nunique);
  } else {
    IntegerVector out(ng);
    if(narm) {
      LogicalVector which(ng);
      int ngs = 0;
      for(int i = 0; i != l; ++i) {
        if(x[i] == NA_LOGICAL) continue;
        if(x[i] == which[g[i]-1]) {
          out[g[i]-1] = 1;
        } else {
          which[g[i]-1] = x[i];
          ++out[g[i]-1];
          if(out[g[i]-1] == 2) {
            ++ngs;
            if(ngs == ng) break;
          }
        }
      }
    } else {
      LogicalVector seen1(ng, true), seen2(ng, true), seen3(ng, true);
      for(int i = 0; i != l; ++i) { // better way?? 
        if(seen1[g[i]-1] && x[i] == NA_LOGICAL) {
          ++out[g[i]-1];
          seen1[g[i]-1] = false;
        } else if(seen2[g[i]-1] && x[i] == true) {
          ++out[g[i]-1];
          seen2[g[i]-1] = false;
        } else if(seen3[g[i]-1] && x[i] == false) {
          ++out[g[i]-1];
          seen3[g[i]-1] = false;
        }
      }
    }
    SHALLOW_DUPLICATE_ATTRIB(out, x);
    return out;
  }
}

template <>
IntegerVector fNdistinct2Impl(const Vector<CPLXSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector fNdistinct2Impl(const Vector<VECSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector fNdistinct2Impl(const Vector<RAWSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector fNdistinct2Impl(const Vector<EXPRSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm) {
  stop("Not supported SEXP type!");
}


// [[Rcpp::export]]
SEXP fNdistinctCpp(SEXP x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, bool narm = true) {
  RCPP_RETURN_VECTOR(fNdistinct2Impl, x, ng, g, gs, narm);
}


// [[Rcpp::export]] // Better Solution ?? // What about string ??
SEXP fNdistinctlCpp(const List& x, int ng = 0, const IntegerVector& g = 0, 
                    const SEXP& gs = R_NilValue, bool narm = true, bool drop = true) {
  int l = x.size();
  List out(l);
  
  for(int j = l; j--; ) {
    switch(TYPEOF(x[j])) {
    case REALSXP:
      out[j] = fNdistinct2Impl<REALSXP>(x[j], ng, g, gs, narm);
      break;
    case INTSXP:
      out[j] = fNdistinct2Impl<INTSXP>(x[j], ng, g, gs, narm);
      break;
    case STRSXP:
      out[j] = fNdistinct2Impl<STRSXP>(x[j], ng, g, gs, narm);
      break;
    case LGLSXP:
      out[j] = fNdistinct2Impl<LGLSXP>(x[j], ng, g, gs, narm);
      break;
    default: stop("Not supported SEXP type !");
    }
  }
  if(drop && ng == 0) {
    IntegerVector res = no_init_vector(l);
    for(int i = l; i--; ) res[i] = out[i]; // Rf_coerceVector(out, INTSXP); // doesn't work !!
    res.attr("names") = x.attr("names");
    return res;
  } else {
    DUPLICATE_ATTRIB(out, x);
    if(ng == 0) out.attr("row.names") = 1;
    else out.attr("row.names") = NumericVector::create(NA_REAL, -ng);
    return out;
  }
}





// Previous std::sort based solution: Seems faster for Numeric Data !!! (about 100 milliseconds faster on NWDI grouped, but 300 milliseconds slower ungroupe) --------------
// // 
// inline bool isnan2(double x) {
//   return x != x;
// }
// 
// 
// template <int RTYPE>
// IntegerVector fNdistinctREAL(const Vector<RTYPE>& xx, int ng, IntegerVector g, SEXP gs, bool narm) {
//   int l = xx.size();
//   NumericVector x = Rf_duplicate(xx);
// 
//   if(ng == 0) {
//     int Nunique = 0;
//     if(narm) {
//       auto pend = std::remove_if(x.begin(), x.end(), isnan2);
//       std::sort(x.begin(), pend);
//       Nunique = std::unique(x.begin(), pend) - x.begin();
//     } else {
//       std::sort(x.begin(), x.end());
//       Nunique = std::unique(x.begin(), x.end()) - x.begin();
//     }
//     return IntegerVector::create(Nunique);
//   } else {
//     if(l != g.size()) stop("length(g) must match length(x)");
//     std::vector<std::vector<double> > gmap(ng);
//     IntegerVector out(ng);
//     if(Rf_isNull(gs)) {
//       for(int i = 0; i != l; ++i) ++out[g[i]-1];
//       for(int i = 0; i != ng; ++i) gmap[i] = std::vector<double> (out[i]); // better using vector of vectors (faster initialization) ??
//       memset(out, 0, sizeof(int)*ng);
//       // return out;
//     } else {
//       IntegerVector gsv = gs;
//       if(ng != gsv.size()) stop("ng must match length(gs)");
//       for(int i = 0; i != ng; ++i) gmap[i] = std::vector<double> (gsv[i]); // better using vector of vectors (faster initialization) ?? -> nope !!
//     }
//     for(int i = 0; i != l; ++i) gmap[g[i]-1][out[g[i]-1]++] = x[i];
//     if(narm) {
//       for(int gr = 0; gr != ng; ++gr) {
//         auto pend = std::remove_if(gmap[gr].begin(), gmap[gr].end(), isnan2);
//         std::sort(gmap[gr].begin(), pend);
//         out[gr] = std::unique(gmap[gr].begin(), pend) - gmap[gr].begin();
//       }
//     } else {
//       for(int gr = 0; gr != ng; ++gr) {
//         std::sort(gmap[gr].begin(), gmap[gr].end());
//         out[gr] = std::unique(gmap[gr].begin(), gmap[gr].end()) - gmap[gr].begin();
//       }
//     }
//     if(Rf_getAttrib(x, R_ClassSymbol) == R_NilValue) SHALLOW_DUPLICATE_ATTRIB(out, x);
//     return out;
//   }
// }
// 
// template <int RTYPE>
// IntegerVector fNdistinctImpl(const Vector<RTYPE>& xx, int ng, IntegerVector g, SEXP gs, bool narm) {
//   int l = xx.size();
//   Vector<RTYPE> x = Rf_duplicate(xx);
// 
//   typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
// 
//   if(ng == 0) {
//     int Nunique = 0;
//     if(narm) {
//       // auto pend = (TYPEOF(x) == REALSXP) ? std::remove_if(x.begin(), x.end(), isnan2)
//       //                                    : std::remove(x.begin(), x.end(), Vector<RTYPE>::get_na());
//       auto pend = std::remove(x.begin(), x.end(), Vector<RTYPE>::get_na());
//       std::sort(x.begin(), pend);
//       Nunique = std::unique(x.begin(), pend) - x.begin();
//     } else {
//       std::sort(x.begin(), x.end());
//       Nunique = std::unique(x.begin(), x.end()) - x.begin();
//     }
//     return IntegerVector::create(Nunique);
//   } else {
//     if(l != g.size()) stop("length(g) must match length(x)");
//     std::vector<std::vector<storage_t> > gmap(ng);
//     IntegerVector out(ng);
//     if(Rf_isNull(gs)) {
//       for(int i = 0; i != l; ++i) ++out[g[i]-1];
//       for(int i = 0; i != ng; ++i) gmap[i] = std::vector<storage_t> (out[i]); // better using vector of vectors (faster initialization) ??
//       memset(out, 0, sizeof(int)*ng);
//       // return out;
//     } else {
//       IntegerVector gsv = gs;
//       if(ng != gsv.size()) stop("ng must match length(gs)");
//       for(int i = 0; i != ng; ++i) gmap[i] = std::vector<storage_t> (gsv[i]); // better using vector of vectors (faster initialization) ?? -> nope !!
//     }
//     for(int i = 0; i != l; ++i) gmap[g[i]-1][out[g[i]-1]++] = x[i];
//     if(narm) {
//       // if(TYPEOF(x) == REALSXP) {
//       //   for(int gr = 0; gr != ng; ++gr) {
//       //     auto pend = std::remove_if(gmap[gr].begin(), gmap[gr].end(), isnan2);
//       //     std::sort(gmap[gr].begin(), pend);
//       //     out[gr] = std::unique(gmap[gr].begin(), pend) - gmap[gr].begin();
//       //   }
//       // } else {
//         for(int gr = 0; gr != ng; ++gr) {
//           auto pend = std::remove(gmap[gr].begin(), gmap[gr].end(), Vector<RTYPE>::get_na());
//           std::sort(gmap[gr].begin(), pend);
//           out[gr] = std::unique(gmap[gr].begin(), pend) - gmap[gr].begin();
//         }
//       // }
//     } else {
//       for(int gr = 0; gr != ng; ++gr) {
//         std::sort(gmap[gr].begin(), gmap[gr].end());
//         out[gr] = std::unique(gmap[gr].begin(), gmap[gr].end()) - gmap[gr].begin();
//       }
//     }
//     if(Rf_getAttrib(x, R_ClassSymbol) == R_NilValue) SHALLOW_DUPLICATE_ATTRIB(out, x);
//     return out;
//   }
// }
// 
// // [[Rcpp::export]]
// IntegerVector fNdistinctCpp(SEXP x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, bool narm = true) {
//  switch(TYPEOF(x)) {
//  case REALSXP: return fNdistinctREAL<REALSXP>(x, ng, g, gs, narm);
//  case INTSXP: return fNdistinctImpl<INTSXP>(x, ng, g, gs, narm);
//  case STRSXP: return fNdistinctImpl<STRSXP>(x, ng, g, gs, narm);
//  case LGLSXP: return fNdistinctImpl<LGLSXP>(x, ng, g, gs, narm);
//  default: stop("Not supported SEXP type !");
//  }
// }
// 
// // [[Rcpp::export]]
// SEXP fNdistinctmCpp(SEXP x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, bool narm = true, bool drop = true) {
// 
//     switch(TYPEOF(x)) {
//     case REALSXP: {
//       NumericMatrix xm = x;
//       int col = xm.ncol();
//       if(ng == 0) {
//         IntegerVector out = no_init_vector(col);
//         for(int j = col; j--; ) out[j] = fNdistinctREAL<REALSXP>(xm( _ ,j), ng, g, gs, narm)[0];
//         if(drop) out.attr("names") = colnames(xm);
//         else {
//          out.attr("dim") = Dimension(1, col);
//          colnames(out) = colnames(xm);
//         }
//         return out;
//       } else {
//         IntegerMatrix out = no_init_matrix(ng, col);
//         for(int j = col; j--; ) out(_ , j) = fNdistinctREAL<REALSXP>(xm( _ ,j), ng, g, gs, narm);
//         colnames(out) = colnames(xm);
//         return out;
//       }
//     }
//     case INTSXP: {
//       IntegerMatrix xm = x;
//       int col = xm.ncol();
//       if(ng == 0) {
//         IntegerVector out = no_init_vector(col);
//         for(int j = col; j--; ) out[j] = fNdistinctImpl<INTSXP>(xm(_ , j), ng, g, gs, narm)[0];
//         if(drop) out.attr("names") = colnames(xm);
//         else {
//           out.attr("dim") = Dimension(1, col);
//           colnames(out) = colnames(xm);
//         }
//         return out;
//       } else {
//         IntegerMatrix out = no_init_matrix(ng, col);
//         for(int j = col; j--; ) out(_ , j) = fNdistinctImpl<INTSXP>(xm(_ , j), ng, g, gs, narm);
//         colnames(out) = colnames(xm);
//         return out;
//       }
//     }
//     case STRSXP: {
//       CharacterMatrix xm = x;
//       int col = xm.ncol();
//       if(ng == 0) {
//         IntegerVector out = no_init_vector(col);
//         for(int j = col; j--; ) out[j] = fNdistinctImpl<STRSXP>(xm(_ , j), ng, g, gs, narm)[0];
//         if(drop) out.attr("names") = colnames(xm);
//         else {
//           out.attr("dim") = Dimension(1, col);
//           colnames(out) = colnames(xm);
//         }
//         return out;
//       } else {
//         IntegerMatrix out = no_init_matrix(ng, col);
//         for(int j = col; j--; ) out(_ , j) = fNdistinctImpl<STRSXP>(xm(_ , j), ng, g, gs, narm);
//         colnames(out) = colnames(xm);
//         return out;
//       }
//     }
//     case LGLSXP: {
//       LogicalMatrix xm = x;
//       int col = xm.ncol();
//       if(ng == 0) {
//         IntegerVector out = no_init_vector(col);
//         for(int j = col; j--; ) out[j] = fNdistinctImpl<LGLSXP>(xm(_ , j), ng, g, gs, narm)[0];
//         if(drop) out.attr("names") = colnames(xm);
//         else {
//           out.attr("dim") = Dimension(1, col);
//           colnames(out) = colnames(xm);
//         }
//       } else {
//         IntegerMatrix out = no_init_matrix(ng, col);
//         for(int j = col; j--; ) out(_ , j) = fNdistinctImpl<LGLSXP>(xm(_ , j), ng, g, gs, narm);
//         colnames(out) = colnames(xm);
//         return out;
//       }
//     }
//     default: stop("Not supported SEXP type !");
//     }
// }
// 
// // [[Rcpp::export]]
// List fNdistinctlCpp(List x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, bool narm = true) {
//   int l = x.size();
//   List out(l);
// 
//   for(int j = l; j--; ) {
//     switch(TYPEOF(x[j])) {
//     case REALSXP:
//       out[j] = fNdistinctREAL<REALSXP>(x[j], ng, g, gs, narm);
//       break;
//     case INTSXP:
//       out[j] = fNdistinctImpl<INTSXP>(x[j], ng, g, gs, narm);
//       break;
//     case STRSXP:
//       out[j] = fNdistinctImpl<STRSXP>(x[j], ng, g, gs, narm);
//       break;
//     case LGLSXP:
//       out[j] = fNdistinctImpl<LGLSXP>(x[j], ng, g, gs, narm);
//       break;
//     default: stop("Not supported SEXP type !");
//     }
//   }
//   DUPLICATE_ATTRIB(out, x);
//   if(ng == 0) out.attr("row.names") = 1;
//   else out.attr("row.names") = NumericVector::create(NA_REAL, -ng);
// 
//   return out;
// }
// 
// 
// 







// template <>
// IntegerVector fNdistinctImpl(Vector<CPLXSXP> x, int ng, IntegerVector g, SEXP gs, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// IntegerVector fNdistinctImpl(Vector<VECSXP> x, int ng, IntegerVector g, SEXP gs, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// IntegerVector fNdistinctImpl(Vector<RAWSXP> x, int ng, IntegerVector g, SEXP gs, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// template <>
// IntegerVector fNdistinctImpl(Vector<EXPRSXP> x, int ng, IntegerVector g, SEXP gs, bool narm) {
//   stop("Not supported SEXP type!");
// }
// 
// 
// // [[Rcpp::export]]
// SEXP fNdistinct(SEXP x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, bool narm = true) {
//   RCPP_RETURN_VECTOR(fNdistinctImpl, x, ng, g, gs, narm);
// }

// Version for Numeric vectors only:
// //[[Rcpp::export]]
// IntegerVector fNdistinct(NumericVector x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0, bool narm = true) {
//   int l = x.size();
//   if(ng == 0) {
//     int Nunique = 0;
//     if(narm) {
//       auto pend = std::remove_if(x.begin(), x.end(), isnan2);
//       std::sort(x.begin(), pend);
//       Nunique = std::unique(x.begin(), pend) - x.begin();      
//     } else {
//       std::sort(x.begin(), x.end());
//       Nunique = std::unique(x.begin(), x.end()) - x.begin();
//     }
//     return IntegerVector::create(Nunique);
//   } else {
//     if(l != g.size()) stop("length(g) must match length(x)");
//     if(ng != gs.size()) stop("ng must match length(gs)");
//     std::vector<std::vector<double> > gmap(ng);
//     IntegerVector out(ng);
//     for(int i = 0; i != ng; ++i) gmap[i] = std::vector<double> (gs[i]); // better using vector of vectors (faster initialization) ??
//     for(int i = 0; i != l; ++i) gmap[g[i]-1][out[g[i]-1]++] = x[i]; 
//     if(narm) {
//       for(int gr = 0; gr != ng; ++gr) {
//         auto pend = std::remove_if(gmap[gr].begin(), gmap[gr].end(), isnan2);
//         std::sort(gmap[gr].begin(), pend);
//         out[gr] = std::unique(gmap[gr].begin(), pend) - gmap[gr].begin();      
//       }
//     } else {
//       for(int gr = 0; gr != ng; ++gr) {
//         std::sort(gmap[gr].begin(), gmap[gr].end());
//         out[gr] = std::unique(gmap[gr].begin(), gmap[gr].end()) - gmap[gr].begin();      
//       }
//     }
//     if(Rf_getAttrib(x, R_ClassSymbol) == R_NilValue) SHALLOW_DUPLICATE_ATTRIB(out, x); 
//     return out;
//   }
// }



// Old idea: Trying to get an ordering vector from R
// //[[Rcpp::export]]
// IntegerVector fNdistinct(NumericVector x, int ng = 0, IntegerVector g = 0, IntegerVector r = 0, bool narm = true) {
//     int l = x.size();
//   if(ng == 0) {
//     int Nunique = 0;
//     if(narm) {
//       auto pend = std::remove_if(x.begin(), x.end(), isnan2);
//       std::sort(x.begin(), pend);
//       Nunique = std::unique(x.begin(), pend) - x.begin();      
//     } else {
//       std::sort(x.begin(), x.end());
//       Nunique = std::unique(x.begin(), x.end()) - x.begin();
//     }
//     return IntegerVector::create(Nunique);
//   } else {
//     if(r.size() != g.size()) stop("length(r) must match length(g)");
//     IntegerVector go = g[r];
//     NumericVector xo = x[r];
//     IntegerVector out(ng), last(ng);
//     if(narm) {
//       for(int i = 0; i != l; ++i) {
//         if(std::isnan(xo[i]) || xo[i] == last[g[i]-1]) continue;
//           last[g[i]-1] = xo[i];
//           ++out[g[i]-1];
//       }
//     } else {
//       for(int i = 0; i != l; ++i) {
//         if(xo[i] != last[g[i]-1]) {
//           last[g[i]-1] = xo[i];
//           ++out[g[i]-1];
//         }
//       }
//     }
//     return out;
//   }
// }
    