// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp ;

// General to do: Check if you can do it without unsigned int and n[l+1] but just with int and n[l]!!!!
// also:: perhaps redo everything with data pointers and 2d group indices (instead of filling the 2d structure every time !!): http://www.cplusplus.com/reference/vector/vector/data/
// https://stackoverflow.com/questions/1733143/converting-between-c-stdvector-and-c-array-without-copying?rq=1
// For named vectors, could add right name!!

template <int RTYPE>
inline bool isnaNUM(typename Rcpp::traits::storage_type<RTYPE>::type x) {
  return x != x;
}

template <int RTYPE>
inline bool isnaOTH(typename Rcpp::traits::storage_type<RTYPE>::type x) {
  return x == Vector<RTYPE>::get_na();
}

template <int RTYPE>
Vector<RTYPE> fmodeImpl(const Vector<RTYPE>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm) {
  int l = x.size();
  typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;

  if(Rf_isNull(w)) { // No Weights !!
    if(ng == 0) {
      sugar::IndexHash<RTYPE> hash(x);
      int max = 1, index = 0; //  n[l+1] is unstable !!
      IntegerVector n(l+1); //  = no_init_vector // better for valgrind !!
      storage_t mode = x[0];
      if(narm) {
        int i = 0, end = l-1;
        while(isnanT(mode) && i!=end) mode = x[++i];
        if(i!=end) for( ; i != l; ++i) {
          if(isnanT(x[i])) continue;
          unsigned int addr = hash.get_addr(x[i]);
          while(hash.data[addr] && hash.not_equal(x[hash.data[addr] - 1], x[i])) {
            ++addr;
            if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
          }
          if(!hash.data[addr]) {
            hash.data[addr] = i+1;
            ++hash.size_;
            n[i+1] = 1; // good?? stable ??
          } else {
            index = hash.data[addr];
            if(++n[index] > max) { // good, or create int index
              max = n[index];
              mode = x[i];
            }
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          unsigned int addr = hash.get_addr(x[i]);
          while(hash.data[addr] && hash.not_equal(x[hash.data[addr] - 1], x[i])) {
            ++addr;
            if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
          }
          if(!hash.data[addr]) {
            hash.data[addr] = i+1;
            ++hash.size_;
            n[i+1] = 1; // good?? stable ??
          } else {
            index = hash.data[addr];
            if(++n[index] > max) { // good, or create int index
              max = n[index];
              mode = x[i];
            }
          }
        }
      }
      Vector<RTYPE> out(1, mode);
      DUPLICATE_ATTRIB(out, x); // could add right name for names vectors !!
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
      std::vector<std::vector<storage_t> > gmap(ng); //
      Vector<RTYPE> out = no_init_vector(ng);
      int ngp = ng+1; // , n[ngp]; // good ?? stable ?? -> yes !!
      std::vector<int> n(ngp); // memset(n, 0, sizeof(int)*ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != l; ++i) ++n[g[i]];
        for(int i = 0; i != ng; ++i) {
          if(n[i+1] == 0) stop("group size of 0 encountered");
          gmap[i] = std::vector<storage_t> (n[i+1]); // Vector<RTYPE>
          n[i+1] = 0;
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("group size of 0 encountered");
          gmap[i] = std::vector<storage_t> (gsv[i]); // Vector<RTYPE>
        }
      }
      for(int i = 0; i != l; ++i) gmap[g[i]-1][n[g[i]]++] = x[i];
      if(narm) {
        for(int gr = 0; gr != ng; ++gr) {
          const std::vector<storage_t>& temp = gmap[gr]; // wrap() // good ?? // const Vector<RTYPE>& // better for character strings !!
          sugar::IndexHash<RTYPE> hash(wrap(temp));
          int i = 0, s = hash.n, end = s-1, max = 1; // n[s+1] // fastest ?? use n ??
          IntegerVector n(s+1); //  = no_init_vector // better for valgrind !!
          while(isnanT(temp[i]) && i!=end) ++i;
          out[gr] = temp[i]; // good !!
          if(i!=end) for( ; i != s; ++i) {
            if(isnanT(temp[i])) continue;
            unsigned int addr = hash.get_addr(temp[i]);
            while(hash.data[addr] && hash.not_equal(temp[hash.data[addr] - 1], temp[i])) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            }
            if(!hash.data[addr]) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = 1;
            } else {
              if(++n[hash.data[addr]] > max) { // good, or create int index
                max = n[hash.data[addr]];
                out[gr] = temp[i];
              }
            }
          }
        }
      } else {
        for(int gr = 0; gr != ng; ++gr) {
          const std::vector<storage_t>& temp = gmap[gr]; // good ?? // const Vector<RTYPE>& // wrap()
          sugar::IndexHash<RTYPE> hash(wrap(temp));
          out[gr] = temp[0];
          int s = hash.n, max = 1; // n[s+1] // fastest ?? use n ?? and reset partially ??
          IntegerVector n(s+1); //  = no_init_vector // better for valgrind !!
          for(int i = 0; i != s; ++i) {
            unsigned int addr = hash.get_addr(temp[i]);
            while(hash.data[addr] && hash.not_equal(temp[hash.data[addr] - 1], temp[i])) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            }
            if(!hash.data[addr]) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = 1;
            } else {
              if(++n[hash.data[addr]] > max) { // good, or create int index
                max = n[hash.data[addr]];
                out[gr] = temp[i];
              }
            }
          }
        }
      }
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    }
  } else { // With Weights !!
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match length(x)");

    if(ng == 0) {
      sugar::IndexHash<RTYPE> hash(x);
      double max = DBL_MIN;
      int index = 0;
      NumericVector n(l+1); //  = no_init_vector // better for valgrind !!
      storage_t mode = x[0];
      if(narm) {
        int i = 0, end = l-1;
        while((isnanT(mode) || std::isnan(wg[i])) && i!=end) mode = x[++i];
        if(i!=end) for( ; i != l; ++i) {
          if(isnanT(x[i]) || std::isnan(wg[i])) continue;
          unsigned int addr = hash.get_addr(x[i]);
          while(hash.data[addr] && hash.not_equal(x[hash.data[addr] - 1], x[i])) {
            ++addr;
            if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
          }
          if(!hash.data[addr]) {
            hash.data[addr] = i+1;
            ++hash.size_;
            n[i+1] = wg[i]; // good?? stable ??
            if(wg[i] > max) { // necessary, because second loop only entered for more than one occurrence of the same value !!
              max = wg[i];
              mode = x[i];
            }
          } else {
            index = hash.data[addr];
            n[index] += wg[i];
            if(n[index] > max) { // good, or create int index
              max = n[index];
              mode = x[i];
            }
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          unsigned int addr = hash.get_addr(x[i]);
          while(hash.data[addr] && hash.not_equal(x[hash.data[addr] - 1], x[i])) {
            ++addr;
            if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
          }
          if(!hash.data[addr]) {
            hash.data[addr] = i+1;
            ++hash.size_;
            n[i+1] = wg[i]; // good?? stable ??
            if(wg[i] > max) { // necessary, because second loop only entered for more than one occurrence of the same value !!
              max = wg[i];
              mode = x[i];
            }
          } else {
            index = hash.data[addr];
            n[index] += wg[i];
            if(n[index] > max) { // good, or create int index
              max = n[index];
              mode = x[i];
            }
          }
        }
      }
      Vector<RTYPE> out(1, mode);
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
      std::vector<std::vector<storage_t> > gmap(ng); //
      std::vector<std::vector<double> > wmap(ng);
      Vector<RTYPE> out = no_init_vector(ng);
      int ngp = ng+1; // n[ngp]; // good ?? stable ?? -> yes !!
      std::vector<int> n(ngp); // memset(n, 0, sizeof(int)*ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != l; ++i) ++n[g[i]];
        for(int i = 0; i != ng; ++i) {
          if(n[i+1] == 0) stop("group size of 0 encountered");
          gmap[i] = std::vector<storage_t> (n[i+1]);
          wmap[i] = std::vector<double> (n[i+1]);
          n[i+1] = 0;
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("group size of 0 encountered");
          gmap[i] = std::vector<storage_t> (gsv[i]);
          wmap[i] = std::vector<double> (gsv[i]);
        }
      }
      for(int i = 0; i != l; ++i) {
        gmap[g[i]-1][n[g[i]]] = x[i];
        wmap[g[i]-1][n[g[i]]++] = wg[i];
      }
      if(narm) {
        for(int gr = 0; gr != ng; ++gr) {
          std::vector<storage_t> temp = gmap[gr]; // good ?? // const Vector<RTYPE>& // wrap()
          const std::vector<double>& wtemp = wmap[gr];
          sugar::IndexHash<RTYPE> hash(wrap(temp));
          int i = 0, s = hash.n, end = s-1, index = 0;
          double max = DBL_MIN; // n[s+1]
          NumericVector n(s+1); //  = no_init_vector // better for valgrind !!
          while((isnanT(temp[i]) || std::isnan(wtemp[i])) && i!=end) ++i;
          out[gr] = temp[i]; // good !!
          if(i!=end) for( ; i != s; ++i) {
            if(isnanT(temp[i]) || std::isnan(wtemp[i])) continue;
            unsigned int addr = hash.get_addr(temp[i]);
            while(hash.data[addr] && hash.not_equal(temp[hash.data[addr] - 1], temp[i])) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            }
            if(!hash.data[addr]) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = wtemp[i];
              if(wtemp[i] > max) { // necessary, because second loop only entered for more than one occurrence of the same value !!
                max = wtemp[i];
                out[gr] = temp[i];
              }
            } else {
              index = hash.data[addr];
              n[index] += wtemp[i];
              if(n[index] > max) {
                max = n[index];
                out[gr] = temp[i];
              }
            }
          }
        }
      } else {
        for(int gr = 0; gr != ng; ++gr) {
          std::vector<storage_t> temp = gmap[gr]; // good ?? // const Vector<RTYPE>& // wrap()
          const std::vector<double>& wtemp = wmap[gr];
          sugar::IndexHash<RTYPE> hash(wrap(temp));
          out[gr] = temp[0];
          int s = hash.n, index = 0; // fastest ?? use n ?? and reset partially ??
          double max = DBL_MIN; // n[s+1];
          NumericVector n(s+1); //  = no_init_vector // better for valgrind !!
          for(int i = 0; i != s; ++i) {
            unsigned int addr = hash.get_addr(temp[i]);
            while(hash.data[addr] && hash.not_equal(temp[hash.data[addr] - 1], temp[i])) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            }
            if(!hash.data[addr]) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = wtemp[i];
              if(wtemp[i] > max) { // necessary, because second loop only entered for more than one occurrence of the same value !!
                max = wtemp[i];
                out[gr] = temp[i];
              }
            } else {
              index = hash.data[addr];
              n[index] += wtemp[i];
              if(n[index] > max) {
                max = n[index];
                out[gr] = temp[i];
              }
            }
          }
        }
      }
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    }
  }
}

template <> // No logical vector with sugar::IndexHash<RTYPE> !!!
Vector<LGLSXP> fmodeImpl(const Vector<LGLSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm) {
  int l = x.size();

 if(Rf_isNull(w)) {
  if(ng == 0) {
    int Ntrue = 0, Nfalse = 0;
    if(narm) {
      for(int i = 0; i != l; ++i) {
        if(x[i] == NA_LOGICAL) continue;
        if(x[i]) ++Ntrue;
        else ++Nfalse;
      }
      LogicalVector out(1, (Ntrue == 0 && Nfalse == 0) ? NA_LOGICAL : Ntrue >= Nfalse);
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    } else {
      int NNA = 0;
      for(int i = 0; i != l; ++i) {
        if(x[i] == NA_LOGICAL) ++NNA;
        else if(x[i]) ++Ntrue; // else if is crucial here !!
        else ++Nfalse;
      }
      LogicalVector out(1, (NNA > Ntrue && NNA > Nfalse) ? NA_LOGICAL : Ntrue >= Nfalse);
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    }
  } else {
    if(narm) {
      IntegerVector truefalse(ng);
      LogicalVector out(ng, NA_LOGICAL);
      for(int i = 0; i != l; ++i) {
        if(x[i] == NA_LOGICAL) continue;
        if(x[i]) {
          if(++truefalse[g[i]-1] >= 0) out[g[i]-1] = true;
        } else {
          if(--truefalse[g[i]-1] < 0) out[g[i]-1] = false;
        }
      }
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    } else {
      IntegerVector Ntrue(ng), Nfalse(ng), NNA(ng); // better way ??
      LogicalVector out = no_init_vector(ng);
      for(int i = 0; i != l; ++i) {
        if(x[i] == NA_LOGICAL) ++NNA[g[i]-1];
        else if(x[i]) ++Ntrue[g[i]-1];
        else ++Nfalse[g[i]-1];
      }
      for(int i = ng; i--; ) {
        if(NNA[i] > Ntrue[i] && NNA[i] > Nfalse[i]) out[i] = NA_LOGICAL;
        else out[i] = Ntrue[i] >= Nfalse[i];
      }
      DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    }
  }
 } else {
   NumericVector wg = w;
   if(l != wg.size()) stop("length(w) must match length(x)");

   if(ng == 0) {
     double sumwtrue = 0, sumwfalse = 0;
     if(narm) {
       for(int i = 0; i != l; ++i) {
         if(x[i] == NA_LOGICAL || std::isnan(wg[i])) continue;
         if(x[i]) sumwtrue += wg[i];
         else sumwfalse += wg[i];
       }
       LogicalVector out(1, (sumwtrue == 0 && sumwfalse == 0) ? NA_LOGICAL : sumwtrue >= sumwfalse);
       DUPLICATE_ATTRIB(out, x);
       if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
       return out;
     } else {
       double sumwNA = 0;
       for(int i = 0; i != l; ++i) {
         if(std::isnan(wg[i])) continue;
         if(x[i] == NA_LOGICAL) sumwNA += wg[i];
         else if(x[i]) sumwtrue += wg[i];
         else sumwfalse += wg[i];
       }
       LogicalVector out(1, (sumwNA > sumwtrue && sumwNA > sumwfalse) ? NA_LOGICAL : sumwtrue >= sumwfalse);
       DUPLICATE_ATTRIB(out, x);
       if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
       return out;
     }
   } else {
     if(narm) {
       NumericVector sumwtruefalse(ng);
       LogicalVector out(ng, NA_LOGICAL);
       for(int i = 0; i != l; ++i) {
         if(x[i] == NA_LOGICAL || std::isnan(wg[i])) continue;
         if(x[i]) {
           sumwtruefalse[g[i]-1] += wg[i];
           if(sumwtruefalse[g[i]-1] >= 0) out[g[i]-1] = true;
         } else {
           sumwtruefalse[g[i]-1] -= wg[i];
           if(sumwtruefalse[g[i]-1] < 0) out[g[i]-1] = false;
         }
       }
       DUPLICATE_ATTRIB(out, x);
       if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
       return out;
     } else {
       NumericVector sumwtrue(ng), sumwfalse(ng), sumwNA(ng); // better way ??
       LogicalVector out = no_init_vector(ng);
       for(int i = 0; i != l; ++i) {
         if(std::isnan(wg[i])) continue;
         if(x[i] == NA_LOGICAL) sumwNA[g[i]-1] += wg[i];
         else if(x[i]) sumwtrue[g[i]-1] += wg[i];
         else sumwfalse[g[i]-1] += wg[i];
       }
       for(int i = ng; i--; ) {
         if(sumwNA[i] > sumwtrue[i] && sumwNA[i] > sumwfalse[i]) out[i] = NA_LOGICAL;
         else out[i] = sumwtrue[i] >= sumwfalse[i];
       }
       DUPLICATE_ATTRIB(out, x);
       if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
       return out;
     }
   }
 }
}

template <>
Vector<CPLXSXP> fmodeImpl(const Vector<CPLXSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<VECSXP> fmodeImpl(const Vector<VECSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<RAWSXP> fmodeImpl(const Vector<RAWSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm) {
  stop("Not supported SEXP type!");
}

template <>
Vector<EXPRSXP> fmodeImpl(const Vector<EXPRSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
SEXP fmodeCpp(SEXP x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP w = R_NilValue, bool narm = true) {
  RCPP_RETURN_VECTOR(fmodeImpl, x, ng, g, gs, w, narm);
}


// Replicating weight 2d array all  the time is stupid !!
// [[Rcpp::export]] // Better Solution ?? // What about string ?? -> do like matrix, but keep vector LGLSXP method !!
SEXP fmodelCpp(const List& x, int ng = 0, const IntegerVector& g = 0, const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, bool narm = true) {
  int l = x.size();
  List out(l);

  for(int j = l; j--; ) {
    switch(TYPEOF(x[j])) {
    case REALSXP:
      out[j] = fmodeImpl<REALSXP>(x[j], ng, g, gs, w, narm);
      break;
    case INTSXP:
      out[j] = fmodeImpl<INTSXP>(x[j], ng, g, gs, w, narm);
      break;
    case STRSXP:
      out[j] = fmodeImpl<STRSXP>(x[j], ng, g, gs, w, narm);
      break;
    case LGLSXP:
      out[j] = fmodeImpl<LGLSXP>(x[j], ng, g, gs, w, narm);
      break;
    default: stop("Not supported SEXP type !");
    }
  }
  DUPLICATE_ATTRIB(out, x);
  if(ng == 0) out.attr("row.names") = 1;
  else out.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng); // NumericVector::create(NA_REAL, -ng);
  return out;
}





template <int RTYPE>
SEXP fmodemImpl(const Matrix<RTYPE>& x, int ng, const IntegerVector& g,
                const SEXP& gs, const SEXP& w, bool narm, bool drop) {
  int l = x.nrow(), col = x.ncol();
  typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
  auto isnanT = (RTYPE == REALSXP) ? isnaNUM<RTYPE> : isnaOTH<RTYPE>;

  if(Rf_isNull(w)) {
    if(ng == 0) {
      Vector<RTYPE> out = no_init_vector(col);
      if(narm) {
        for(int j = col; j--; ) {
          ConstMatrixColumn<RTYPE> column = x(_ , j);
          sugar::IndexHash<RTYPE> hash(wrap(column)); // why wrap needed ??
          int i = 0, end = l-1, max = 1, index = 0; // n[l+1]; // stable !!, but this not: IntegerVector n = no_init_vector(l);
          IntegerVector n(l+1); //  = no_init_vector // better for valgrind !!
          storage_t mode = column[0]; // best solution ??
          while(isnanT(mode) && i!=end) mode = column[++i];
          if(i!=end) for( ; i != l; ++i) {
            if(isnanT(column[i])) continue;
            unsigned int addr = hash.get_addr(column[i]);
            while(hash.data[addr] && hash.not_equal(column[hash.data[addr] - 1], column[i])) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            }
            if(!hash.data[addr]) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = 1; // good?? stable ??
            } else {
              index = hash.data[addr];
              if(++n[index] > max) { // good, or create int index
                max = n[index];
                mode = column[i];
              }
            }
          }
          out[j] = mode;
        }
      } else {
        for(int j = col; j--; ) {
          ConstMatrixColumn<RTYPE> column = x(_ , j);
          sugar::IndexHash<RTYPE> hash(wrap(column)); // why wrap needed ??
          int max = 1, index = 0; // n[l+1];
          IntegerVector n(l+1); //  = no_init_vector // better for valgrind !!
          storage_t mode = column[0]; // best solution ??
          for(int i = 0; i != l; ++i) {
            unsigned int addr = hash.get_addr(column[i]);
            while(hash.data[addr] && hash.not_equal(column[hash.data[addr] - 1], column[i])) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            }
            if(!hash.data[addr]) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = 1; // good?? stable ??
            } else {
              index = hash.data[addr];
              if(++n[index] > max) { // good, or create int index
                max = n[index];
                mode = column[i];
              }
            }
          }
          out[j] = mode;
        }
      }
      if(drop) out.attr("names") = colnames(x);
      else { // could do duplicate attrib, but problems, i.e. with ts matrices etc.. !! -> but what about factors or dates ??
        out.attr("dim") = Dimension(1, col);
        colnames(out) = colnames(x);
      }
      return out;
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
      std::vector<std::vector<storage_t> > gmap(ng);
      Matrix<RTYPE> out = no_init_matrix(ng, col);
      int ngp = ng+1; //  n[ngp];
      std::vector<int> n(ngp);
      if(Rf_isNull(gs)) {
        // memset(n, 0, sizeof(int)*ngp);
        for(int i = 0; i != l; ++i) ++n[g[i]];
        for(int i = 0; i != ng; ++i) {
          if(n[i+1] == 0) stop("group size of 0 encountered");
          gmap[i] = std::vector<storage_t> (n[i+1]);
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
          MatrixColumn<RTYPE> outj = out(_, j);
          n.assign(ngp, 0); // memset(n, 0, sizeof(int)*ngp);
          for(int i = 0; i != l; ++i) gmap[g[i]-1][n[g[i]]++] = column[i]; // reading in all the values. Better way ??
          for(int gr = 0; gr != ng; ++gr) {
            const std::vector<storage_t>& temp = gmap[gr]; // wrap() // good ?? // const Vector<RTYPE>&
            sugar::IndexHash<RTYPE> hash(wrap(temp));
            int i = 0, s = hash.n, end = s-1, max = 1; // ns[s+1] // fastest ?? use n ??
            IntegerVector ns(s+1); //  = no_init_vector // better for valgrind !!
            while(isnanT(temp[i]) && i!=end) ++i;
            outj[gr] = temp[i]; // good !!
            if(i!=end) for( ; i != s; ++i) {
              if(isnanT(temp[i])) continue;
              unsigned int addr = hash.get_addr(temp[i]);
              while(hash.data[addr] && hash.not_equal(temp[hash.data[addr] - 1], temp[i])) {
                ++addr;
                if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
              }
              if(!hash.data[addr]) {
                hash.data[addr] = i+1;
                ++hash.size_;
                ns[i+1] = 1;
              } else {
                if(++ns[hash.data[addr]] > max) {
                  max = ns[hash.data[addr]];
                  outj[gr] = temp[i];
                }
              }
            }
          }
        }
      } else {
        for(int j = col; j--; ) {
          ConstMatrixColumn<RTYPE> column = x(_ , j);
          MatrixColumn<RTYPE> outj = out(_, j);
          n.assign(ngp, 0); // memset(n, 0, sizeof(int)*ngp);
          for(int i = 0; i != l; ++i) gmap[g[i]-1][n[g[i]]++] = column[i]; // reading in all the values. Better way ??
          for(int gr = 0; gr != ng; ++gr) {
            const std::vector<storage_t>& temp = gmap[gr]; // wrap() // good ?? // const Vector<RTYPE>&
            sugar::IndexHash<RTYPE> hash(wrap(temp));
            outj[gr] = temp[0];
            int  s = hash.n, max = 1; // ns[s+1]
            IntegerVector ns(s+1); //  = no_init_vector // better for valgrind !!
            for(int i = 0; i != s; ++i) {
              unsigned int addr = hash.get_addr(temp[i]);
              while(hash.data[addr] && hash.not_equal(temp[hash.data[addr] - 1], temp[i])) {
                ++addr;
                if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
              }
              if(!hash.data[addr]) {
                hash.data[addr] = i+1;
                ++hash.size_;
                ns[i+1] = 1;
              } else {
                if(++ns[hash.data[addr]] > max) {
                  max = ns[hash.data[addr]];
                  outj[gr] = temp[i];
                }
              }
            }
          }
        }
      }
      colnames(out) = colnames(x);
      return out;
    }
  } else {
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match length(x)");

    if(ng == 0) {
      Vector<RTYPE> out = no_init_vector(col);
      if(narm) {
        for(int j = col; j--; ) {
          ConstMatrixColumn<RTYPE> column = x(_ , j);
          sugar::IndexHash<RTYPE> hash(wrap(column)); // why wrap needed ??
          int i = 0, end = l-1, index = 0;
          double max = DBL_MIN; // , n[l+1]
          NumericVector n(l+1); //  = no_init_vector // better for valgrind !!
          storage_t mode = column[0]; // best solution ??
          while((isnanT(mode) || std::isnan(wg[i])) && i!=end) mode = column[++i];
          if(i!=end) for( ; i != l; ++i) {
            if(isnanT(column[i]) || std::isnan(wg[i])) continue;
            unsigned int addr = hash.get_addr(column[i]);
            while(hash.data[addr] && hash.not_equal(column[hash.data[addr] - 1], column[i])) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            }
            if(!hash.data[addr]) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = wg[i];
              if(wg[i] > max) { // necessary, because second loop only entered for more than one occurrence of the same value !!
                max = wg[i];
                mode = column[i];
              }
            } else {
              index = hash.data[addr];
              n[index] += wg[i];
              if(n[index] > max) {
                max = n[index];
                mode = column[i];
              }
            }
          }
          out[j] = mode;
        }
      } else {
        for(int j = col; j--; ) {
          ConstMatrixColumn<RTYPE> column = x(_ , j);
          sugar::IndexHash<RTYPE> hash(wrap(column)); // why wrap needed ??
          int index = 0;
          double max = DBL_MIN; // , n[l+1]
          NumericVector n(l+1); //  = no_init_vector // better for valgrind !!
          storage_t mode = column[0]; // best solution ??
          for(int i = 0; i != l; ++i) {
            unsigned int addr = hash.get_addr(column[i]);
            while(hash.data[addr] && hash.not_equal(column[hash.data[addr] - 1], column[i])) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            }
            if(!hash.data[addr]) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = wg[i];
              if(wg[i] > max) { // necessary, because second loop only entered for more than one occurrence of the same value !!
                max = wg[i];
                mode = column[i];
              }
            } else {
              index = hash.data[addr];
              n[index] += wg[i];
              if(n[index] > max) {
                max = n[index];
                mode = column[i];
              }
            }
          }
          out[j] = mode;
        }
      }
      if(drop) out.attr("names") = colnames(x);
      else { // could do duplicate attrib, but problems, i.e. with ts matrices etc.. !! -> but what about factors or dates ??
        out.attr("dim") = Dimension(1, col);
        colnames(out) = colnames(x);
      }
      return out;
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
      std::vector<std::vector<storage_t> > gmap(ng);
      std::vector<std::vector<double> > wmap(ng);
      Matrix<RTYPE> out = no_init_matrix(ng, col);
      int ngp = ng+1; // , n[ngp];
      std::vector<int> n(ngp); // memset(n, 0, sizeof(int)*ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != l; ++i) ++n[g[i]];
        for(int i = 0; i != ng; ++i) {
          if(n[i+1] == 0) stop("group size of 0 encountered");
          gmap[i] = std::vector<storage_t> (n[i+1]);
          wmap[i] = std::vector<double> (n[i+1]);
          n[i+1] = 0;
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("group size of 0 encountered");
          gmap[i] = std::vector<storage_t> (gsv[i]);
          wmap[i] = std::vector<double> (gsv[i]);
        }
      }
      for(int i = 0; i != l; ++i) {
        gmap[g[i]-1][n[g[i]]] = x[i]; // good ?? not column 1 ??
        wmap[g[i]-1][n[g[i]]++] = wg[i];
      }
      if(narm) {
        for(int j = 0; j != col; ++j) {
          ConstMatrixColumn<RTYPE> column = x(_ , j);
          MatrixColumn<RTYPE> outj = out(_, j);
          if(j != 0) {
            n.assign(ngp, 0); // memset(n, 0, sizeof(int)*ngp);
            for(int i = 0; i != l; ++i) gmap[g[i]-1][n[g[i]]++] = column[i]; // reading in all the values. Better way ??
          }
          for(int gr = 0; gr != ng; ++gr) {
            const std::vector<storage_t>& temp = gmap[gr]; // wrap() // good ?? // const Vector<RTYPE>&
            const std::vector<double>& wtemp = wmap[gr];
            sugar::IndexHash<RTYPE> hash(wrap(temp));
            int i = 0, index = 0, s = hash.n, end = s-1; // fastest ?? use n ??
            double max = DBL_MIN; // , ns[s+1];
            NumericVector ns(s+1); //  = no_init_vector // better for valgrind !!
            while((isnanT(temp[i]) || std::isnan(wtemp[i])) && i!=end) ++i;
            outj[gr] = temp[i]; // good !!
            if(i!=end) for( ; i != s; ++i) {
              if(isnanT(temp[i]) || std::isnan(wtemp[i])) continue;
              unsigned int addr = hash.get_addr(temp[i]);
              while(hash.data[addr] && hash.not_equal(temp[hash.data[addr] - 1], temp[i])) {
                ++addr;
                if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
              }
              if(!hash.data[addr]) {
                hash.data[addr] = i+1;
                ++hash.size_;
                ns[i+1] = wtemp[i];
                if(wtemp[i] > max) { // necessary, because second loop only entered for more than one occurrence of the same value !!
                  max = wtemp[i];
                  outj[gr] = temp[i];
                }
              } else {
                index = hash.data[addr];
                ns[index] += wtemp[i];
                if(ns[index] > max) {
                  max = ns[index];
                  outj[gr] = temp[i];
                }
              }
            }
          }
        }
      } else {
        for(int j = 0; j != col; ++j) {
          ConstMatrixColumn<RTYPE> column = x(_ , j);
          MatrixColumn<RTYPE> outj = out(_, j);
          if(j != 0) {
            n.assign(ngp, 0); // memset(n, 0, sizeof(int)*ngp);
            for(int i = 0; i != l; ++i) gmap[g[i]-1][n[g[i]]++] = column[i]; // reading in all the values. Better way ??
          }
          for(int gr = 0; gr != ng; ++gr) {
            const std::vector<storage_t>& temp = gmap[gr]; // wrap() // good ?? // const Vector<RTYPE>&
            const std::vector<double>& wtemp = wmap[gr];
            sugar::IndexHash<RTYPE> hash(wrap(temp));
            int index = 0, s = hash.n; // fastest ?? use n ??
            double max = DBL_MIN; // ns[s+1];
            NumericVector ns(s+1); //  = no_init_vector // better for valgrind !!
            outj[gr] = temp[0];
            for(int i = 0; i != s; ++i) {
              unsigned int addr = hash.get_addr(temp[i]);
              while(hash.data[addr] && hash.not_equal(temp[hash.data[addr] - 1], temp[i])) {
                ++addr;
                if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
              }
              if(!hash.data[addr]) {
                hash.data[addr] = i+1;
                ++hash.size_;
                ns[i+1] = wtemp[i];
                if(wtemp[i] > max) { // necessary, because second loop only entered for more than one occurrence of the same value !!
                  max = wtemp[i];
                  outj[gr] = temp[i];
                }
              } else {
                index = hash.data[addr];
                ns[index] += wtemp[i];
                if(ns[index] > max) {
                  max = ns[index];
                  outj[gr] = temp[i];
                }
              }
            }
          }
        }
      }
      colnames(out) = colnames(x);
      return out;
    }

  }
}

template <> // No logical vector with sugar::IndexHash<RTYPE> !!!
SEXP fmodemImpl(const Matrix<LGLSXP>& x, int ng, const IntegerVector& g,
                const SEXP& gs, const SEXP& w, bool narm, bool drop) {
  int l = x.nrow(), col = x.ncol();

  if(Rf_isNull(w)) {
    if(ng == 0) {
      LogicalVector out = no_init_vector(col);
      if(narm) {
        for(int j = col; j--; ) {
          LogicalMatrix::ConstColumn column = x(_ , j);
          int Ntrue = 0, Nfalse = 0;
          for(int i = 0; i != l; ++i) {
            if(column[i] == NA_LOGICAL) continue;
            if(column[i]) ++Ntrue;
            else ++Nfalse;
          }
          if(Ntrue == 0 && Nfalse == 0) out[j] = NA_LOGICAL;
          else out[j] = Ntrue >= Nfalse;
        }
      } else {
        for(int j = col; j--; ) {
          LogicalMatrix::ConstColumn column = x(_ , j);
          int Ntrue = 0, Nfalse = 0, NNA = 0;
          for(int i = 0; i != l; ++i) {
            if(column[i] == NA_LOGICAL) ++NNA; // error before, should put contunie or else if aftert this !!
            else if(column[i]) ++Ntrue;
            else ++Nfalse;
          }
          if(NNA > Ntrue && NNA > Nfalse) out[j] = NA_LOGICAL;
          else out[j] = Ntrue >= Nfalse;
        }
      }
      if(drop) out.attr("names") = colnames(x);
      else {
        out.attr("dim") = Dimension(1, col);
        colnames(out) = colnames(x);
      }
      return out;
    } else {
      LogicalMatrix out = no_init_matrix(ng, col);
      if(narm) {
        std::fill(out.begin(), out.end(), NA_LOGICAL); // memset(out, NA_LOGICAL, sizeof(bool)*ng*col); // doesn't work !! -> not a matrix anymore afterwards !!
        for(int j = col; j--; ) {
          LogicalMatrix::ConstColumn column = x(_ , j);
          LogicalMatrix::Column outj = out(_, j);
          IntegerVector truefalse(ng);
          for(int i = 0; i != l; ++i) {
            if(column[i] == NA_LOGICAL) continue;
            if(column[i]) {
              if(++truefalse[g[i]-1] >= 0) outj[g[i]-1] = true;
            } else {
              if(--truefalse[g[i]-1] < 0) outj[g[i]-1] = false;
            }
          }
        }
      } else {
        for(int j = col; j--; ) {
          LogicalMatrix::ConstColumn column = x(_ , j);
          LogicalMatrix::Column outj = out(_, j);
          IntegerVector Ntrue(ng), Nfalse(ng), NNA(ng);
          for(int i = 0; i != l; ++i) {
            if(column[i] == NA_LOGICAL) ++NNA[g[i]-1];
            else if(column[i]) ++Ntrue[g[i]-1];
            else ++Nfalse[g[i]-1];
          }
          for(int i = ng; i--; ) {
            if(NNA[i] > Ntrue[i] && NNA[i] > Nfalse[i]) outj[i] = NA_LOGICAL;
            else outj[i] = Ntrue[i] >= Nfalse[i];
          }
        }
      }
      colnames(out) = colnames(x);
      return out;
    }
  } else {
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match length(x)");

    if(ng == 0) {
      LogicalVector out = no_init_vector(col);
      if(narm) {
        for(int j = col; j--; ) {
          LogicalMatrix::ConstColumn column = x(_ , j);
          double sumwtrue = 0, sumwfalse = 0;
          for(int i = 0; i != l; ++i) {
            if(column[i] == NA_LOGICAL || std::isnan(wg[i])) continue;
            if(column[i]) sumwtrue += wg[i];
            else sumwfalse += wg[i];
          }
          if(sumwtrue == 0 && sumwfalse == 0) out[j] = NA_LOGICAL;
          else out[j] = sumwtrue >= sumwfalse;
        }
      } else {
        for(int j = col; j--; ) {
          LogicalMatrix::ConstColumn column = x(_ , j);
          double sumwtrue = 0, sumwfalse = 0, sumwNA = 0;
          for(int i = 0; i != l; ++i) {
            if(std::isnan(wg[i])) continue;
            if(column[i] == NA_LOGICAL) sumwNA += wg[i]; // error before, should put contunie or else if aftert this !!
            else if(column[i]) sumwtrue += wg[i];
            else sumwfalse += wg[i];
          }
          if(sumwNA > sumwtrue && sumwNA > sumwfalse) out[j] = NA_LOGICAL;
          else out[j] = sumwtrue >= sumwfalse;
        }
      }
      if(drop) out.attr("names") = colnames(x);
      else {
        out.attr("dim") = Dimension(1, col);
        colnames(out) = colnames(x);
      }
      return out;
    } else {
      LogicalMatrix out = no_init_matrix(ng, col);
      if(narm) {
        std::fill(out.begin(), out.end(), NA_LOGICAL); // memset(out, NA_LOGICAL, sizeof(bool)*ng*col); // doesn't work !! -> not a matrix anymore afterwards !!
        for(int j = col; j--; ) {
          LogicalMatrix::ConstColumn column = x(_ , j);
          LogicalMatrix::Column outj = out(_, j);
          NumericVector sumwtruefalse(ng);
          for(int i = 0; i != l; ++i) {
            if(column[i] == NA_LOGICAL || std::isnan(wg[i])) continue;
            if(column[i]) {
              sumwtruefalse[g[i]-1] += wg[i];
              if(sumwtruefalse[g[i]-1] >= 0) outj[g[i]-1] = true;
            } else {
              sumwtruefalse[g[i]-1] -= wg[i];
              if(sumwtruefalse[g[i]-1] < 0) outj[g[i]-1] = false;
            }
          }
        }
      } else {
        for(int j = col; j--; ) {
          LogicalMatrix::ConstColumn column = x(_ , j);
          LogicalMatrix::Column outj = out(_, j);
          NumericVector sumwtrue(ng), sumwfalse(ng), sumwNA(ng);
          for(int i = 0; i != l; ++i) {
            if(std::isnan(wg[i])) continue;
            if(column[i] == NA_LOGICAL) sumwNA[g[i]-1] += wg[i];
            else if(column[i]) sumwtrue[g[i]-1] += wg[i];
            else sumwfalse[g[i]-1] += wg[i];
          }
          for(int i = ng; i--; ) {
            if(sumwNA[i] > sumwtrue[i] && sumwNA[i] > sumwfalse[i]) outj[i] = NA_LOGICAL;
            else outj[i] = sumwtrue[i] >= sumwfalse[i];
          }
        }
      }
      colnames(out) = colnames(x);
      return out;
    }
  }
}

template <>
SEXP fmodemImpl(const Matrix<CPLXSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fmodemImpl(const Matrix<VECSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fmodemImpl(const Matrix<RAWSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fmodemImpl(const Matrix<EXPRSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, bool drop) {
  stop("Not supported SEXP type!");
}


// [[Rcpp::export]]
SEXP fmodemCpp(SEXP x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP w = R_NilValue, bool narm = true, bool drop = true) {
  RCPP_RETURN_MATRIX(fmodemImpl, x, ng, g, gs, w, narm, drop);
}

