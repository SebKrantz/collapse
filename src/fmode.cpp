// [[Rcpp::plugins(cpp11)]]
#define STRICT_R_HEADERS
#include <cfloat>
#include <Rcpp.h>
using namespace Rcpp ;

// General to do: Check if you can do it without unsigned int and n[l+1] but just with int and n[l]
// also:: perhaps redo everything with data pointers and 2d group indices (instead of filling the 2d structure every time): http://www.cplusplus.com/reference/vector/vector/data/
// https://stackoverflow.com/questions/1733143/converting-between-c-stdvector-and-c-array-without-copying?rq=1
// For named vectors, could add right name!


template <int RTYPE>
Vector<RTYPE> fmodeImpl(const Vector<RTYPE>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, int ret) {
  int l = x.size();
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101
  bool minm = ret == 1, nfirstm = ret > 0, lastm = ret == 3;
  typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
  auto isnanT = (RTYPE == REALSXP) ? [](storage_t x) { return x != x; } :
    [](storage_t x) { return x == Vector<RTYPE>::get_na(); };

  unsigned int addr;

  if(Rf_isNull(w)) { // No Weights
    if(ng == 0) {
      sugar::IndexHash<RTYPE> hash(x);
      int max = 1, index = 0; //  n[l+1] is unstable
      std::vector<int> n(l+1); //  = no_init_vector // better for valgrind
      storage_t mode = hash.src[0];
      if(narm) {
        int i = 0, end = l-1;
        while(isnanT(mode) && i!=end) mode = hash.src[++i];
        if(i!=end) {
          for( ; i != l; ++i) {
          storage_t val = hash.src[i];
          if(isnanT(val)) continue;
          addr = hash.get_addr(val);
          index = hash.data[addr];
          while(index && hash.not_equal(hash.src[index - 1], val)) {
            ++addr;
            if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            index = hash.data[addr];
          }
          if(!index) {
            hash.data[addr] = i+1;
            ++hash.size_;
            n[i+1] = 1;
            if(nfirstm && max == 1) { // Could also do this at the end in a separate loop. What is faster ? -> This seems better !
              if(lastm) mode = val;
              else if(minm) {
                if(mode > val) mode = val;
              } else {
                if(mode < val) mode = val;
              }
            }
          } else {
            // if(++n[index] > max) { // good, or create int index
            //   max = n[index];
            //   mode = val;
            // }
            if(++n[index] >= max) {
              if(lastm || n[index] > max) {
                max = n[index];
                mode = val;
              } else if(nfirstm) {
                if(minm) {
                  if(mode > val) mode = val;
                } else {
                  if(mode < val) mode = val;
                }
              }
            }
          }
        }
          // if(nfirstm && max == 1) { // Above seems better !
          //   if(minm) {
          //     for(int i = 1; i != l; ++i) if(mode > x[i]) mode = x[i];
          //   } else {
          //     for(int i = 1; i != l; ++i) if(mode < x[i]) mode = x[i];
          //   }
          // }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          storage_t val = hash.src[i];
          addr = hash.get_addr(val);
          index = hash.data[addr];
          while(index && hash.not_equal(hash.src[index - 1], val)) {
            ++addr;
            if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            index = hash.data[addr];
          }
          if(!index) {
            hash.data[addr] = i+1;
            ++hash.size_;
            n[i+1] = 1;
            if(nfirstm && max == 1) { // Could also do this at the end in a separate loop. What is faster ? -> This seems better !
              if(lastm) mode = val;
              else if(minm) {
                if(mode > val) mode = val;
              } else {
                if(mode < val) mode = val;
              }
            }
          } else {
            // if(++n[index] > max) { // good, or create int index
            //   max = n[index];
            //   mode = val;
            // }
            if(++n[index] >= max) { // good, or create int index
              if(lastm || n[index] > max) {
                max = n[index];
                mode = val;
              } else if(nfirstm) {
                if(minm) {
                  if(mode > val) mode = val;
                } else {
                  if(mode < val) mode = val;
                }
              }
            }
          }
        }
        // if(nfirstm && max == 1) { // above seems better
        //   if(minm) {
        //     for(int i = 1; i != l; ++i) if(mode > x[i]) mode = x[i];
        //   } else {
        //     for(int i = 1; i != l; ++i) if(mode < x[i]) mode = x[i];
        //   }
        // }
      }
      Vector<RTYPE> out(1, mode);
      SHALLOW_DUPLICATE_ATTRIB(out, x); // could add right name for named vectors
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
      const int *pg = g.begin();
      int ngp = ng+1;
      std::vector<std::vector<storage_t> > gmap(ngp);
      Vector<RTYPE> out = no_init_vector(ng);
      std::vector<int> n(ngp); // memset(n, 0, sizeof(int)*ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != l; ++i) ++n[pg[i]];
        for(int i = 1; i != ngp; ++i) {
          if(n[i] == 0) stop("Group size of 0 encountered. This is probably due to unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i] = std::vector<storage_t> (n[i]); // Vector<RTYPE>
          n[i] = 0;
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("Group size of 0 encountered. This is probably due to unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i+1] = std::vector<storage_t> (gsv[i]); // Vector<RTYPE>
        }
      }
      for(int i = 0; i != l; ++i) gmap[pg[i]][n[pg[i]]++] = x[i];
      if(narm) {
        for(int gr = 0; gr != ng; ++gr) {
          // const std::vector<storage_t>& temp = gmap[gr]; // wrap() // good ? // const Vector<RTYPE>& // better for character strings
          sugar::IndexHash<RTYPE> hash(wrap(gmap[gr+1]));  // wrap(temp)
          int i = 0, s = hash.n, end = s-1, max = 1, index; // n[s+1] // fastest ? use n ?
          while(isnanT(hash.src[i]) && i!=end) ++i;
          out[gr] = hash.src[i]; // good
          if(i!=end) {
            std::vector<int> n(s+1); //  = no_init_vector // better for valgrind
            for( ; i != s; ++i) {
            storage_t val = hash.src[i];
            if(isnanT(val)) continue;
            addr = hash.get_addr(val);
            index = hash.data[addr];
            while(index && hash.not_equal(hash.src[index - 1], val)) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
              index = hash.data[addr];
            }
            if(!index) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = 1;
              if(nfirstm && max == 1) { // Could also do this at the end in a separate loop. What is faster ? -> This seems better !
                if(lastm) out[gr] = val;
                else if(minm) {
                  if(out[gr] > val) out[gr] = val;
                } else {
                  if(out[gr] < val) out[gr] = val;
                }
              }
            } else {
              // if(++n[hash.data[addr]] > max) { // good, or create int index
              //   max = n[hash.data[addr]];
              //   out[gr] = val;
              // }
              // index = hash.data[addr];
              if(++n[index] >= max) {
                if(lastm || n[index] > max) {
                  max = n[index];
                  out[gr] = val;
                } else if(nfirstm) {
                  if(minm) {
                    if(out[gr] > val) out[gr] = val;
                  } else {
                    if(out[gr] < val) out[gr] = val;
                  }
                }
              }
            }
          }
            // if(nfirstm && max == 1) { // Above seems better !
            //   if(minm) {
            //     for(int i = 1; i != s; ++i) if(out[gr] > hash.src[i]) out[gr] = hash.src[i];
            //   } else {
            //     for(int i = 1; i != s; ++i) if(out[gr] < hash.src[i]) out[gr] = hash.src[i];
            //   }
            // }
          }
        }
      } else {
        for(int gr = 0; gr != ng; ++gr) {
          // const std::vector<storage_t>& temp = gmap[gr]; // good ? // const Vector<RTYPE>& // wrap()
          sugar::IndexHash<RTYPE> hash(wrap(gmap[gr+1])); // wrap(temp)
          out[gr] = hash.src[0];
          int s = hash.n, max = 1, index; // n[s+1] // fastest ? use n ? and reset partially ?
          std::vector<int> n(s+1); //  = no_init_vector // better for valgrind
          for(int i = 0; i != s; ++i) {
            storage_t val = hash.src[i];
            addr = hash.get_addr(val);
            index = hash.data[addr];
            while(index && hash.not_equal(hash.src[index - 1], val)) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
              index = hash.data[addr];
            }
            if(!index) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = 1;
              if(nfirstm && max == 1) { // Could also do this at the end in a separate loop. What is faster ? -> This seems better !
                if(lastm) out[gr] = val;
                else if(minm) {
                  if(out[gr] > val) out[gr] = val;
                } else {
                  if(out[gr] < val) out[gr] = val;
                }
              }
            } else {
              // if(++n[hash.data[addr]] > max) { // good, or create int index
              //   max = n[hash.data[addr]];
              //   out[gr] = val;
              // }
              if(++n[index] >= max) {
                if(lastm || n[index] > max) {
                  max = n[index];
                  out[gr] = val;
                } else if(nfirstm) {
                  if(minm) {
                    if(out[gr] > val) out[gr] = val;
                  } else {
                    if(out[gr] < val) out[gr] = val;
                  }
                }
              }
            }
          }
          // if(nfirstm && max == 1) { // Above seems better !
          //   if(minm) {
          //     for(int i = 1; i != s; ++i) if(out[gr] > hash.src[i]) out[gr] = hash.src[i];
          //   } else {
          //     for(int i = 1; i != s; ++i) if(out[gr] < hash.src[i]) out[gr] = hash.src[i];
          //   }
          // }
        }
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    }
  } else { // With Weights
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match length(x)");
    double *pwg = wg.begin();

    if(ng == 0) {
      sugar::IndexHash<RTYPE> hash(x);
      double max = DBL_MIN;
      int index = 0;
      std::vector<double> n(l+1); //  = no_init_vector // better for valgrind
      storage_t mode = hash.src[0];
      if(narm) {
        int i = 0, end = l-1;
        while((isnanT(mode) || std::isnan(pwg[i])) && i!=end) mode = hash.src[++i];
        if(i!=end) for( ; i != l; ++i) {
          storage_t val = hash.src[i];
          if(isnanT(val) || std::isnan(pwg[i])) continue;
          addr = hash.get_addr(val);
          index = hash.data[addr];
          while(index && hash.not_equal(hash.src[index - 1], val)) {
            ++addr;
            if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            index = hash.data[addr];
          }
          if(!index) {
            hash.data[addr] = i+1;
            ++hash.size_;
            n[i+1] = pwg[i];
            if(pwg[i] >= max) { // necessary, because second loop only entered for more than one occurrence of the same value
              if(lastm || pwg[i] > max) {
                max = pwg[i];
                mode = val;
              } else if(nfirstm) { // Could also do this at the end in a separate loop. What is faster ??
                if(minm) {
                  if(mode > val) mode = val;
                } else {
                  if(mode < val) mode = val;
                }
              }
            }
          } else {
            n[index] += pwg[i];
            // if(n[index] > max) { // good, or create int index
            //   max = n[index];
            //   mode = val;
            // }
            if(n[index] >= max) {
              if(lastm || n[index] > max) {
                max = n[index];
                mode = val;
              } else if(nfirstm) {
                if(minm) {
                  if(mode > val) mode = val;
                } else {
                  if(mode < val) mode = val;
                }
              }
            }
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(std::isnan(pwg[i])) continue;
          storage_t val = hash.src[i];
          addr = hash.get_addr(val);
          index = hash.data[addr];
          while(index && hash.not_equal(hash.src[index - 1], val)) {
            ++addr;
            if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
            index = hash.data[addr];
          }
          if(!index) {
            hash.data[addr] = i+1;
            ++hash.size_;
            n[i+1] = pwg[i];
            if(pwg[i] >= max) { // necessary, because second loop only entered for more than one occurrence of the same value
              if(lastm || pwg[i] > max) {
                max = pwg[i];
                mode = val;
              } else if(nfirstm) { // Could also do this at the end in a separate loop. What is faster ??
                if(minm) {
                  if(mode > val) mode = val;
                } else {
                  if(mode < val) mode = val;
                }
              }
            }
          } else {
            n[index] += pwg[i];
            // if(n[index] > max) { // good, or create int index
            //   max = n[index];
            //   mode = val;
            // }
            if(n[index] >= max) { // good, or create int index
              if(lastm || n[index] > max) {
                max = n[index];
                mode = val;
              } else if(nfirstm) {
                if(minm) {
                  if(mode > val) mode = val;
                } else {
                  if(mode < val) mode = val;
                }
              }
            }
          }
        }
      }
      Vector<RTYPE> out(1, mode);
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
      const int *pg = g.begin();
      int ngp = ng+1;
      std::vector<std::vector<storage_t> > gmap(ngp);
      std::vector<std::vector<double> > wmap(ngp);
      Vector<RTYPE> out = no_init_vector(ng);
      std::vector<int> n(ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != l; ++i) ++n[pg[i]];
        for(int i = 1; i != ngp; ++i) {
          if(n[i] == 0) stop("Group size of 0 encountered. This is probably due to unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i] = std::vector<storage_t> (n[i]);
          wmap[i] = std::vector<double> (n[i]);
          n[i] = 0;
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("Group size of 0 encountered. This is probably due to unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i+1] = std::vector<storage_t> (gsv[i]);
          wmap[i+1] = std::vector<double> (gsv[i]);
        }
      }
      for(int i = 0; i != l; ++i) {
        int gi = pg[i];
        gmap[gi][n[gi]] = x[i];
        wmap[gi][n[gi]++] = pwg[i];
      }
      if(narm) {
        for(int gr = 0; gr != ng; ++gr) {
          // const std::vector<storage_t>& temp = gmap[gr]; // good ? // const Vector<RTYPE>& // wrap()
          const std::vector<double>& wtemp = wmap[gr+1];
          sugar::IndexHash<RTYPE> hash(wrap(gmap[gr+1])); // wrap(temp)
          int i = 0, s = hash.n, end = s-1, index;
          double max = DBL_MIN; // n[s+1]
          while((isnanT(hash.src[i]) || std::isnan(wtemp[i])) && i!=end) ++i;
          out[gr] = hash.src[i]; // good !
          if(i!=end) {
            std::vector<double> n(s+1); //  = no_init_vector // better for valgrind
            for( ; i != s; ++i) {
              storage_t val = hash.src[i];
              if(isnanT(val) || std::isnan(wtemp[i])) continue;
              addr = hash.get_addr(val);
              index = hash.data[addr];
              while(index && hash.not_equal(hash.src[index - 1], val)) {
                ++addr;
                if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
                index = hash.data[addr];
              }
              if(!index) {
                hash.data[addr] = i+1;
                ++hash.size_;
                n[i+1] = wtemp[i];
                if(wtemp[i] >= max) { // necessary, because second loop only entered for more than one occurrence of the same value
                  if(lastm || wtemp[i] > max) {
                    max = wtemp[i];
                    out[gr] = val;
                  } else if(nfirstm) { // Could also do this at the end in a separate loop. What is faster ??
                    if(minm) {
                      if(out[gr] > val) out[gr] = val;
                    } else {
                      if(out[gr] < val) out[gr] = val;
                    }
                  }
                }
              } else {
                n[index] += wtemp[i];
                // if(n[index] > max) {
                //   max = n[index];
                //   out[gr] = val;
                // }
                if(n[index] >= max) {
                  if(lastm || n[index] > max) {
                    max = n[index];
                    out[gr] = val;
                  } else if(nfirstm) {
                    if(minm) {
                      if(out[gr] > val) out[gr] = val;
                    } else {
                      if(out[gr] < val) out[gr] = val;
                    }
                  }
                }
              }
            }
          }
        }
      } else {
        for(int gr = 0; gr != ng; ++gr) {
          // const std::vector<storage_t>& temp = gmap[gr]; // good ? // const Vector<RTYPE>& // wrap()
          const std::vector<double>& wtemp = wmap[gr+1];
          sugar::IndexHash<RTYPE> hash(wrap(gmap[gr+1])); // wrap(temp)
          out[gr] = hash.src[0];
          int s = hash.n, index; // fastest ? use n ? and reset partially ?
          double max = DBL_MIN; // n[s+1];
          std::vector<double> n(s+1); //  = no_init_vector // better for valgrind
          for(int i = 0; i != s; ++i) {
            if(std::isnan(wtemp[i])) continue;
            storage_t val = hash.src[i];
            addr = hash.get_addr(val);
            index = hash.data[addr];
            while(index && hash.not_equal(hash.src[index - 1], val)) {
              ++addr;
              if(addr == static_cast<unsigned int>(hash.m)) addr = 0;
              index = hash.data[addr];
            }
            if(!index) {
              hash.data[addr] = i+1;
              ++hash.size_;
              n[i+1] = wtemp[i];
              if(wtemp[i] >= max) { // necessary, because second loop only entered for more than one occurrence of the same value
                if(lastm || wtemp[i] > max) {
                  max = wtemp[i];
                  out[gr] = val;
                } else if(nfirstm) { // Could also do this at the end in a separate loop. What is faster ??
                  if(minm) {
                    if(out[gr] > val) out[gr] = val;
                  } else {
                    if(out[gr] < val) out[gr] = val;
                  }
                }
              }
            } else {
              n[index] += wtemp[i];
              // if(n[index] > max) {
              //   max = n[index];
              //   out[gr] = val;
              // }
              if(n[index] >= max) {
                if(lastm || n[index] > max) {
                  max = n[index];
                  out[gr] = val;
                } else if(nfirstm) {
                  if(minm) {
                    if(out[gr] > val) out[gr] = val;
                  } else {
                    if(out[gr] < val) out[gr] = val;
                  }
                }
              }
            }
          }
        }
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    }
  }
}


IntegerVector fmodeFACT(const IntegerVector& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, int ret) {
  int l = x.size(), nlevp = Rf_nlevels(x)+1, val = 0;
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101
  bool minm = ret == 1, nfirstm = ret > 0, lastm = ret == 3;

  if(Rf_isNull(w)) { // No Weights
    if(ng == 0) { // No Groups
      std::vector<int> n(nlevp);
      int max = 1, mode = x[0];
      if(narm) {
        int i = 0, end = l-1;
        while(mode == NA_INTEGER && i!=end) mode = x[++i];
        if(i!=end) {
          for( ; i != l; ++i) {
            val = x[i];
            if(val == NA_INTEGER) continue;
            if(++n[val] >= max) {
              if(lastm || n[val] > max) {
                max = n[val];
                mode = val;
              } else if(nfirstm) {
                if(minm) {
                  if(mode > val) mode = val;
                } else {
                  if(mode < val) mode = val;
                }
              }
            }
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          val = x[i];
          if(val == NA_INTEGER) val = 0;
          if(++n[val] >= max) {
            if(lastm || n[val] > max) {
              max = n[val];
              mode = val;
            } else if(nfirstm) {
              if(minm) {
                if(mode > val) mode = val;
              } else {
                if(mode < val) mode = val;
              }
            }
          }
        }
        if(mode == 0) mode = NA_INTEGER;
      }
      IntegerVector out(1, mode);
      SHALLOW_DUPLICATE_ATTRIB(out, x); // could add right name for names vectors
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
      const int *pg = g.begin();
      int ngp = ng+1;
      std::vector<std::vector<int> > gmap(ngp);
      IntegerVector out = no_init_vector(ng);
      std::vector<int> n(ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != l; ++i) ++n[pg[i]];
        for(int i = 1; i != ngp; ++i) {
          if(n[i] == 0) stop("Group size of 0 encountered. This is probably due to unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i] = std::vector<int> (n[i]);
          n[i] = 0;
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("Group size of 0 encountered. This is probably due to unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i+1] = std::vector<int> (gsv[i]);
        }
      }
      for(int i = 0; i != l; ++i) gmap[pg[i]][n[pg[i]]++] = x[i];
      if(narm) {
        for(int gr = 0; gr != ng; ++gr) {
          const std::vector<int>& temp = gmap[gr+1];
          int i = 0, s = temp.size(), end = s-1, max = 1;
          while(temp[i] == NA_INTEGER && i!=end) ++i;
          out[gr] = temp[i];
          if(i!=end) {
            std::vector<int> n(nlevp);
            for( ; i != s; ++i) {
              val = temp[i];
              if(val == NA_INTEGER) continue;
              if(++n[val] >= max) {
                if(lastm || n[val] > max) {
                  max = n[val];
                  out[gr] = val;
                } else if(nfirstm) {
                  if(minm) {
                    if(out[gr] > val) out[gr] = val;
                  } else {
                    if(out[gr] < val) out[gr] = val;
                  }
                }
              }
            }
          }
        }
      } else {
        for(int gr = 0; gr != ng; ++gr) {
          const std::vector<int>& temp = gmap[gr+1];
          int tl = temp.size(), max = 1;
          std::vector<int> n(nlevp);
          out[gr] = temp[0];
          for(int i = 0; i != tl; ++i) {
            val = temp[i];
            if(val == NA_INTEGER) val = 0;
            if(++n[val] >= max) {
              if(lastm || n[val] > max) {
                max = n[val];
                out[gr] = val;
              } else if(nfirstm) {
                if(minm) {
                  if(out[gr] > val) out[gr] = val;
                } else {
                  if(out[gr] < val) out[gr] = val;
                }
              }
            }
          }
          if(out[gr] == 0) out[gr] = NA_INTEGER;
        }
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    }
  } else { // With Weights
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match length(x)");
    double *pwg = wg.begin();

    if(ng == 0) {
      double max = DBL_MIN;
      int mode = x[0];
      std::vector<double> n(nlevp);

      if(narm) {
        int i = 0, end = l-1;
        while((mode == NA_INTEGER || std::isnan(pwg[i])) && i!=end) mode = x[++i];
        if(i!=end) for( ; i != l; ++i) {
          val = x[i];
          if(val == NA_INTEGER || std::isnan(pwg[i])) continue;
          n[val] += pwg[i];
          if(n[val] >= max) {
            if(lastm || n[val] > max) {
              max = n[val];
              mode = val;
            } else if(nfirstm) {
              if(minm) {
                if(mode > val) mode = val;
              } else {
                if(mode < val) mode = val;
              }
            }
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(std::isnan(pwg[i])) continue;
          val = x[i];
          if(val == NA_INTEGER) val = 0;
          n[val] += pwg[i];
          if(n[val] >= max) {
            if(lastm || n[val] > max) {
              max = n[val];
              mode = val;
            } else if(nfirstm) {
              if(minm) {
                if(mode > val) mode = val;
              } else {
                if(mode < val) mode = val;
              }
            }
          }
        }
        if(mode == 0) mode = NA_INTEGER;
      }
      IntegerVector out(1, mode);
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
      const int *pg = g.begin();
      int ngp = ng+1;
      std::vector<std::vector<int> > gmap(ngp);
      std::vector<std::vector<double> > wmap(ngp);
      IntegerVector out = no_init_vector(ng);
      std::vector<int> n(ngp);
      if(Rf_isNull(gs)) {
        for(int i = 0; i != l; ++i) ++n[pg[i]];
        for(int i = 1; i != ngp; ++i) {
          if(n[i] == 0) stop("Group size of 0 encountered. This is probably due to unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i] = std::vector<int> (n[i]);
          wmap[i] = std::vector<double> (n[i]);
          n[i] = 0;
        }
      } else {
        IntegerVector gsv = gs;
        if(ng != gsv.size()) stop("ng must match length(gs)");
        for(int i = 0; i != ng; ++i) {
          if(gsv[i] == 0) stop("Group size of 0 encountered. This is probably due to unused factor levels. Use fdroplevels(f) to drop them.");
          gmap[i+1] = std::vector<int> (gsv[i]);
          wmap[i+1] = std::vector<double> (gsv[i]);
        }
      }
      for(int i = 0; i != l; ++i) {
        int gi = pg[i];
        gmap[gi][n[gi]] = x[i];
        wmap[gi][n[gi]++] = pwg[i];
      }
      if(narm) {
        for(int gr = 0; gr != ng; ++gr) {
          const std::vector<int>& temp = gmap[gr+1];
          const std::vector<double>& wtemp = wmap[gr+1];
          int i = 0, s = temp.size(), end = s-1;
          double max = DBL_MIN;
          while((temp[i] == NA_INTEGER || std::isnan(wtemp[i])) && i!=end) ++i;
          out[gr] = temp[i];
          if(i!=end) {
            std::vector<double> n(nlevp);
            for( ; i != s; ++i) {
              val = temp[i];
              if(val == NA_INTEGER || std::isnan(wtemp[i])) continue;
              n[val] += wtemp[i];
              if(n[val] >= max) {
                if(lastm || n[val] > max) {
                  max = n[val];
                  out[gr] = val;
                } else if(nfirstm) {
                  if(minm) {
                    if(out[gr] > val) out[gr] = val;
                  } else {
                    if(out[gr] < val) out[gr] = val;
                  }
                }
              }
            }
          }
        }
      } else {
        for(int gr = 0; gr != ng; ++gr) {
          const std::vector<int>& temp = gmap[gr+1];
          const std::vector<double>& wtemp = wmap[gr+1];
          int tl = temp.size();
          double max = DBL_MIN;
          std::vector<double> n(nlevp);
          out[gr] = temp[0];
          for(int i = 0; i != tl; ++i) {
            if(std::isnan(wtemp[i])) continue;
            val = temp[i];
            if(val == NA_INTEGER) val = 0;
            n[val] += wtemp[i];
            if(n[val] >= max) {
              if(lastm || n[val] > max) {
                max = n[val];
                out[gr] = val;
              } else if(nfirstm) {
                if(minm) {
                  if(out[gr] > val) out[gr] = val;
                } else {
                  if(out[gr] < val) out[gr] = val;
                }
              }
            }
          }
          if(out[gr] == 0) out[gr] = NA_INTEGER;
        }
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    }
  }
}



template <> // No logical vector with sugar::IndexHash<RTYPE> !
Vector<LGLSXP> fmodeImpl(const Vector<LGLSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, int ret) {
  int l = x.size();
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101
  bool maxm = ret != 1;

  if(Rf_isNull(w)) {
    if(ng == 0) {
      int Ntrue = 0, Nfalse = 0;
      LogicalVector out(1);
      if(narm) {
        for(int i = 0; i != l; ++i) {
          if(x[i] == NA_LOGICAL) continue;
          if(x[i]) ++Ntrue;
          else ++Nfalse;
        }
        out[0] = (Ntrue == 0 && Nfalse == 0) ? NA_LOGICAL : (maxm || Ntrue != Nfalse) ? Ntrue >= Nfalse : false;
      } else {
        int NNA = 0;
        for(int i = 0; i != l; ++i) {
          if(x[i] == NA_LOGICAL) ++NNA;
          else if(x[i]) ++Ntrue; // else if is crucial here !
          else ++Nfalse;
        }
        out[0] = (NNA > Ntrue && NNA > Nfalse) ? NA_LOGICAL : (maxm || Ntrue != Nfalse) ? Ntrue >= Nfalse : false;
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
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
        if(!maxm) for(int i = ng; i--; ) if(truefalse[i] == 0 && out[i] != NA_LOGICAL) out[i] = false;
        SHALLOW_DUPLICATE_ATTRIB(out, x);
        if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
        return out;
      } else {
        IntegerVector Ntrue(ng), Nfalse(ng), NNA(ng); // better way ?
        LogicalVector out = no_init_vector(ng);
        for(int i = 0; i != l; ++i) {
          if(x[i] == NA_LOGICAL) ++NNA[g[i]-1];
          else if(x[i]) ++Ntrue[g[i]-1];
          else ++Nfalse[g[i]-1];
        }
        if(maxm) {
          for(int i = ng; i--; ) {
            if(NNA[i] > Ntrue[i] && NNA[i] > Nfalse[i]) out[i] = NA_LOGICAL;
            else out[i] = Ntrue[i] >= Nfalse[i];
          }
        } else {
          for(int i = ng; i--; ) {
            if(NNA[i] > Ntrue[i] && NNA[i] > Nfalse[i]) out[i] = NA_LOGICAL;
            else out[i] = Ntrue[i] > Nfalse[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(out, x);
        if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
        return out;
      }
    }
  } else {
    NumericVector wg = w;
    if(l != wg.size()) stop("length(w) must match length(x)");

    if(ng == 0) {
      LogicalVector out(1);
      double sumwtrue = 0, sumwfalse = 0;
      if(narm) {
        for(int i = 0; i != l; ++i) {
          if(x[i] == NA_LOGICAL || std::isnan(wg[i])) continue;
          if(x[i]) sumwtrue += wg[i];
          else sumwfalse += wg[i];
        }
        out[0] = (sumwtrue == 0 && sumwfalse == 0) ? NA_LOGICAL : (maxm || sumwtrue != sumwfalse) ? sumwtrue >= sumwfalse : false;
      } else {
        double sumwNA = 0;
        for(int i = 0; i != l; ++i) {
          if(std::isnan(wg[i])) continue;
          if(x[i] == NA_LOGICAL) sumwNA += wg[i];
          else if(x[i]) sumwtrue += wg[i];
          else sumwfalse += wg[i];
        }         // important as w could be NA as well..
        out[0] = ((sumwNA > sumwtrue && sumwNA > sumwfalse) || (sumwtrue == 0 && sumwfalse == 0)) ? NA_LOGICAL : (maxm || sumwtrue != sumwfalse) ? sumwtrue >= sumwfalse : false;
      }
      SHALLOW_DUPLICATE_ATTRIB(out, x);
      if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
      return out;
    } else {
      if(l != g.size()) stop("length(g) must match length(x)");
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
        if(!maxm) for(int i = ng; i--; ) if(sumwtruefalse[i] == 0 && out[i] != NA_LOGICAL) out[i] = false;
        SHALLOW_DUPLICATE_ATTRIB(out, x);
        if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
        return out;
      } else {
        NumericVector sumwtrue(ng), sumwfalse(ng), sumwNA(ng); // better way ?
        LogicalVector out = no_init_vector(ng);
        for(int i = 0; i != l; ++i) {
          if(std::isnan(wg[i])) continue;
          if(x[i] == NA_LOGICAL) sumwNA[g[i]-1] += wg[i];
          else if(x[i]) sumwtrue[g[i]-1] += wg[i];
          else sumwfalse[g[i]-1] += wg[i];
        }
        if(maxm) {
          for(int i = ng; i--; ) {
            // important as w could be NA as well..
            if((sumwNA[i] > sumwtrue[i] && sumwNA[i] > sumwfalse[i]) || sumwtrue[i] + sumwfalse[i] == 0) out[i] = NA_LOGICAL;
            else out[i] = sumwtrue[i] >= sumwfalse[i];
          }
        } else {
          for(int i = ng; i--; ) {
            // important as w could be NA as well..
            if((sumwNA[i] > sumwtrue[i] && sumwNA[i] > sumwfalse[i]) || sumwtrue[i] + sumwfalse[i] == 0) out[i] = NA_LOGICAL;
            else out[i] = sumwtrue[i] > sumwfalse[i];
          }
        }
        SHALLOW_DUPLICATE_ATTRIB(out, x);
        if(Rf_getAttrib(x, R_NamesSymbol) != R_NilValue) Rf_setAttrib(out, R_NamesSymbol, R_NilValue);
        return out;
      }
    }
  }
}


// [[Rcpp::export]]
SEXP fmodeCpp(const SEXP& x, int ng = 0, const IntegerVector& g = 0, const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, bool narm = true, int ret = 0) {
  switch(TYPEOF(x)) {
  case REALSXP: return fmodeImpl<REALSXP>(x, ng, g, gs, w, narm, ret);
  case INTSXP:
    if(Rf_isFactor(x) && (ng == 0 || Rf_nlevels(x) < Rf_length(x) / ng * 3))
      return fmodeFACT(x, ng, g, gs, w, narm, ret);
    return fmodeImpl<INTSXP>(x, ng, g, gs, w, narm, ret);
  case STRSXP: return fmodeImpl<STRSXP>(x, ng, g, gs, w, narm, ret);
  case LGLSXP: return fmodeImpl<LGLSXP>(x, ng, g, gs, w, narm, ret);
  default: stop("Not supported SEXP type !");
  }
}


// Replicating weight 2d array all the time is stupid
// [[Rcpp::export]]   // Better Solution ? // What about string ? -> do like matrix, but keep vector LGLSXP method
SEXP fmodelCpp(const List& x, int ng = 0, const IntegerVector& g = 0, const SEXP& gs = R_NilValue, const SEXP& w = R_NilValue, bool narm = true, int ret = 0) {
  int l = x.size();
  List out(l);

  for(int j = l; j--; ) {
    switch(TYPEOF(x[j])) {
    case REALSXP:
      out[j] = fmodeImpl<REALSXP>(x[j], ng, g, gs, w, narm, ret);
      break;
    case INTSXP:
      if(Rf_isFactor(x[j]) && (ng == 0 || Rf_nlevels(x[j]) < Rf_length(x[j]) / ng * 3))
        out[j] = fmodeFACT(x[j], ng, g, gs, w, narm, ret);
      else out[j] = fmodeImpl<INTSXP>(x[j], ng, g, gs, w, narm, ret);
      break;
    case STRSXP:
      out[j] = fmodeImpl<STRSXP>(x[j], ng, g, gs, w, narm, ret);
      break;
    case LGLSXP:
      out[j] = fmodeImpl<LGLSXP>(x[j], ng, g, gs, w, narm, ret);
      break;
    default: stop("Not supported SEXP type !");
    }
  }
  SHALLOW_DUPLICATE_ATTRIB(out, x);
  if(ng == 0) Rf_setAttrib(out, R_RowNamesSymbol, Rf_ScalarInteger(1));
  else Rf_setAttrib(out, R_RowNamesSymbol, IntegerVector::create(NA_INTEGER, -ng));
  return out;
}



template <int RTYPE>
SEXP fmodemImpl(const Matrix<RTYPE>& x, int ng, const IntegerVector& g,
                const SEXP& gs, const SEXP& w, bool narm, bool drop, int ret) {
  int col = x.ncol();
  Matrix<RTYPE> out = (ng == 0) ? no_init_matrix(1, col) : no_init_matrix(ng, col);
  for(int j = col; j--; ) out(_, j) = fmodeImpl<RTYPE>(x(_, j), ng, g, gs, w, narm, ret);
  if(drop && ng == 0) {
    Rf_setAttrib(out, R_DimSymbol, R_NilValue); // Rf_dimgets(out, R_NilValue); -> Doesn't work !
    Rf_setAttrib(out, R_NamesSymbol, colnames(x));
  } else {
    colnames(out) = colnames(x);
    if(!Rf_isObject(x)) Rf_copyMostAttrib(x, out);
  }
  return out;
}


template <>
SEXP fmodemImpl(const Matrix<CPLXSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, bool drop, int ret) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fmodemImpl(const Matrix<VECSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, bool drop, int ret) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fmodemImpl(const Matrix<RAWSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, bool drop, int ret) {
  stop("Not supported SEXP type!");
}

template <>
SEXP fmodemImpl(const Matrix<EXPRSXP>& x, int ng, const IntegerVector& g, const SEXP& gs, const SEXP& w, bool narm, bool drop, int ret) {
  stop("Not supported SEXP type!");
}


// [[Rcpp::export]]
SEXP fmodemCpp(SEXP x, int ng = 0, IntegerVector g = 0, SEXP gs = R_NilValue, SEXP w = R_NilValue, bool narm = true, bool drop = true, int ret = 0) {
  RCPP_RETURN_MATRIX(fmodemImpl, x, ng, g, gs, w, narm, drop, ret);
}

