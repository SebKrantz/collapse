// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp ;

// TODO: Perhaps redo everything with data pointers and 2d group indices (instead of filling the 2d structure every time !): http://www.cplusplus.com/reference/vector/vector/data/
// https://stackoverflow.com/questions/1733143/converting-between-c-stdvector-and-c-array-without-copying?rq=1

// improve logical method ?

template <int RTYPE>
inline bool isnaNUM(typename Rcpp::traits::storage_type<RTYPE>::type x) {
  return x != x;
}

template <int RTYPE>
inline bool isnaOTH(typename Rcpp::traits::storage_type<RTYPE>::type x) {
  return x == Vector<RTYPE>::get_na();
}

template <int RTYPE>
IntegerVector fNdistinctImpl(const Vector<RTYPE>& x, int ng, const IntegerVector& g,
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
        if(out[i] == 0) stop("group size of 0 encountered");
        gmap[i] = std::vector<storage_t> (out[i]);
        out[i] = 0;
      }
      // memset(out, 0, sizeof(int)*ng); // Stable ? -> Nope, gives error
    } else {
      IntegerVector gsv = gs;
      if(ng != gsv.size()) stop("ng must match length(gs)");
      for(int i = 0; i != ng; ++i) {
        if(gsv[i] == 0) stop("group size of 0 encountered");
        gmap[i] = std::vector<storage_t> (gsv[i]);
      }
    }
    for(int i = 0; i != l; ++i) gmap[g[i]-1][out[g[i]-1]++] = x[i];
    if(narm) {
      for(int gr = 0; gr != ng; ++gr) {
        const std::vector<storage_t>& temp = gmap[gr]; // good ? // const Vector<RTYPE>& // wrap()
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

template <> // No logical vector with sugar::IndexHash<RTYPE> !
IntegerVector fNdistinctImpl(const Vector<LGLSXP>& x, int ng, const IntegerVector& g,
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
      for(int i = 0; i != l; ++i) { // better way?
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
      for(int i = 0; i != l; ++i) { // better way?
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
IntegerVector fNdistinctImpl(const Vector<CPLXSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector fNdistinctImpl(const Vector<VECSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector fNdistinctImpl(const Vector<RAWSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector fNdistinctImpl(const Vector<EXPRSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm) {
  stop("Not supported SEXP type!");
}


// [[Rcpp::export]]
SEXP fNdistinctCpp(SEXP x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, bool narm = true) {
  RCPP_RETURN_VECTOR(fNdistinctImpl, x, ng, g, gs, narm);
}


// [[Rcpp::export]] // Better Solution ? // What about string ?
SEXP fNdistinctlCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
                    const SEXP& gs = R_NilValue, bool narm = true, bool drop = true) {
  int l = x.size();
  List out(l);

  for(int j = l; j--; ) {
    switch(TYPEOF(x[j])) {
    case REALSXP:
      out[j] = fNdistinctImpl<REALSXP>(x[j], ng, g, gs, narm);
      break;
    case INTSXP:
      out[j] = fNdistinctImpl<INTSXP>(x[j], ng, g, gs, narm);
      break;
    case STRSXP:
      out[j] = fNdistinctImpl<STRSXP>(x[j], ng, g, gs, narm);
      break;
    case LGLSXP:
      out[j] = fNdistinctImpl<LGLSXP>(x[j], ng, g, gs, narm);
      break;
    default: stop("Not supported SEXP type !");
    }
  }
  if(drop && ng == 0) {
    IntegerVector res = no_init_vector(l);
    for(int i = l; i--; ) res[i] = out[i]; // Rf_coerceVector(out, INTSXP); // doesn't work
    res.attr("names") = x.attr("names");
    return res;
  } else {
    DUPLICATE_ATTRIB(out, x);
    if(ng == 0) out.attr("row.names") = 1;
    else out.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng);
    return out;
  }
}






template <int RTYPE>
SEXP fNdistinctmImpl(const Matrix<RTYPE>& x, int ng, const IntegerVector& g,
                      const SEXP& gs, bool narm, bool drop) {
  int l = x.nrow(), col = x.ncol();
  typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;

  if(ng == 0) {
    IntegerVector out = no_init_vector(col);
    if(narm) {
      for(int j = col; j--; ) {
        ConstMatrixColumn<RTYPE> column = x(_ , j);
        sugar::IndexHash<RTYPE> hash(wrap(column)); // why wrap needed ?
        for(int i = 0; i != l; ++i) {
          if(isnanT(column[i])) continue;
          unsigned int addr = hash.get_addr(column[i]);
          while(hash.data[addr] && hash.not_equal(column[hash.data[addr] - 1], column[i])) {
            ++addr;
            if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
          }
          if(!hash.data[addr]) {
            hash.data[addr] = i+1;
            ++hash.size_;
          }
        }
        out[j] = hash.size_;
      }
    } else {
      for(int j = col; j--; ) {
        ConstMatrixColumn<RTYPE> column = x(_ , j);
        sugar::IndexHash<RTYPE> hash(wrap(column)); // why wrap needed ?
        hash.fill();
        out[j] = hash.size_;
      }
    }
    if(drop) out.attr("names") = colnames(x);
    else {
      out.attr("dim") = Dimension(1, col);
      colnames(out) = colnames(x);
    }
    return out;
  } else {
    if(l != g.size()) stop("length(g) must match length(x)");
    std::vector<std::vector<storage_t> > gmap(ng);
    IntegerMatrix out = no_init_matrix(ng, col);
    std::vector<int> n(ng);
    if(Rf_isNull(gs)) {
      // memset(n, 0, sizeof(int)*ng);
      for(int i = 0; i != l; ++i) ++n[g[i]-1];
      for(int i = 0; i != ng; ++i) {
        if(n[i] == 0) stop("group size of 0 encountered");
        gmap[i] = std::vector<storage_t> (n[i]);
      }
    } else {
      IntegerVector gsv = gs;
      if(ng != gsv.size()) stop("ng must match length(gs)");
      for(int i = 0; i != ng; ++i) {
        if(gsv[i] == 0) stop("group size of 0 encountered");
        gmap[i] = std::vector<storage_t> (gsv[i]);
      }
    }
    if(narm) {
      for(int j = col; j--; ) {
        ConstMatrixColumn<RTYPE> column = x(_ , j);
        IntegerMatrix::Column outj = out(_, j);
        n.assign(ng, 0); // memset(n, 0, sizeof(int)*ng);
        for(int i = 0; i != l; ++i) gmap[g[i]-1][n[g[i]-1]++] = column[i]; // reading in all the values. Better way ?
        for(int gr = 0; gr != ng; ++gr) {
          const std::vector<storage_t>& temp = gmap[gr]; // good ? // const Vector<RTYPE>& // wrap()
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
          outj[gr] = hash.size_;
        }
      }
    } else {
      for(int j = col; j--; ) {
        ConstMatrixColumn<RTYPE> column = x(_ , j);
        IntegerMatrix::Column outj = out(_, j);
        n.assign(ng, 0); // memset(n, 0, sizeof(int)*ng);
        for(int i = 0; i != l; ++i) gmap[g[i]-1][n[g[i]-1]++] = column[i]; // reading in all the values. Better way ?
        for(int gr = 0; gr != ng; ++gr) {
          sugar::IndexHash<RTYPE> hash(wrap(gmap[gr]));
          hash.fill();
          outj[gr] = hash.size_;
        }
      }
    }
    colnames(out) = colnames(x);
    return out;
  }
}

template <> // No logical vector with sugar::IndexHash<RTYPE> !
SEXP fNdistinctmImpl(const Matrix<LGLSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm, bool drop) {
  int l = x.nrow(), col = x.ncol();

  if(ng == 0) {
    IntegerVector out(col);
    if(narm) {
      for(int j = col; j--; ) {
        LogicalMatrix::ConstColumn column = x(_ , j);
        bool which = true;
        for(int i = 0; i != l; ++i) {
          if(column[i] == NA_LOGICAL) continue;
          if(column[i] == which) {
            out[j] = 1;
          } else {
            which = column[i];
            ++out[j];
            if(out[j] == 2) break;
          }
        }
      }
    } else {
      for(int j = col; j--; ) {
        LogicalMatrix::ConstColumn column = x(_ , j);
        bool seen1 = true, seen2 = true, seen3 = true;
        for(int i = 0; i != l; ++i) { // better way?
          if(seen1 && column[i] == NA_LOGICAL) {
            ++out[j];
            seen1 = false;
          } else if(seen2 && column[i] == true) {
            ++out[j];
            seen2 = false;
          } else if(seen3 && column[i] == false) {
            ++out[j];
            seen3 = false;
          }
          if(out[j] == 3) break;
        }
      }
    }
    if(drop) out.attr("names") = colnames(x);
    else {
      out.attr("dim") = Dimension(1, col);
      colnames(out) = colnames(x);
    }
    return out;
  } else {
    IntegerMatrix out(ng, col); //  = no_init_matrix
    if(narm) {
      for(int j = col; j--; ) {
        LogicalMatrix::ConstColumn column = x(_ , j);
        IntegerMatrix::Column outj = out(_, j);
        LogicalVector which(ng);
        int ngs = 0;
        for(int i = 0; i != l; ++i) {
          if(column[i] == NA_LOGICAL) continue;
          if(column[i] == which[g[i]-1]) {
            outj[g[i]-1] = 1;
          } else {
            which[g[i]-1] = column[i];
            ++outj[g[i]-1];
            if(outj[g[i]-1] == 2) {
              ++ngs;
              if(ngs == ng) break;
            }
          }
        }
      }
    } else {
      for(int j = col; j--; ) {
        LogicalMatrix::ConstColumn column = x(_ , j);
        IntegerMatrix::Column outj = out(_, j);
        LogicalVector seen1(ng, true), seen2(ng, true), seen3(ng, true);
        for(int i = 0; i != l; ++i) { // better way?
          if(seen1[g[i]-1] && column[i] == NA_LOGICAL) {
            ++outj[g[i]-1];
            seen1[g[i]-1] = false;
          } else if(seen2[g[i]-1] && column[i] == true) {
            ++outj[g[i]-1];
            seen2[g[i]-1] = false;
          } else if(seen3[g[i]-1] && column[i] == false) {
            ++outj[g[i]-1];
            seen3[g[i]-1] = false;
          }
        }
      }
    }
    colnames(out) = colnames(x);
    return out;
  }
}

template <>
SEXP fNdistinctmImpl(const Matrix<CPLXSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fNdistinctmImpl(const Matrix<VECSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fNdistinctmImpl(const Matrix<RAWSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fNdistinctmImpl(const Matrix<EXPRSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}


// [[Rcpp::export]]
SEXP fNdistinctmCpp(SEXP x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, bool narm = true, bool drop = true) {
  RCPP_RETURN_MATRIX(fNdistinctmImpl, x, ng, g, gs, narm, drop);
}
