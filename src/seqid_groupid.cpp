#include <Rcpp.h>
using namespace Rcpp;

// TODO: Optimize !
// TODO: can do something about doubles using == ?
// TODO: Option na_fill ?

// Note: For x[i] == NA_INTEGER, which is equal to INT_MIN, cannot calculate x[i]-prev ! -> fixed in 1.2.1
// https://stackoverflow.com/questions/776624/whats-faster-iterating-an-stl-vector-with-vectoriterator-or-with-at


// [[Rcpp::export]]
IntegerVector seqid(const IntegerVector& x, const SEXP& o = R_NilValue, int del = 1, int start = 1,
                    bool na_skip = false, bool skip_seq = false, bool check_o = true) {
  int l = x.size(), id = start, prev;
  if(l < 1) return x; // Prevents seqfault for numeric(0) #101

  IntegerVector out = no_init_vector(l);
  if(Rf_isNull(o)) {
    if(na_skip) {
      int j = 0, end = l;
      while(x[j] == NA_INTEGER && j != end) out[j++] = NA_INTEGER;
      if(j != end) {
        prev = x[j];
        out[j] = id;
        for(int i = j+1; i != l; ++i) {
          if(x[i] != NA_INTEGER) {
            if(x[i] - prev != del) ++id; // x[i]-x[i-1]?
            prev = x[i];
            out[i] = id;
          } else { // Faster way ?
            out[i] = NA_INTEGER;
            if(skip_seq) prev += del;
          }
        }
      }
    } else {
      int nafill = INT_MAX - 1e7;
      prev = x[0];
      if(prev == NA_INTEGER) prev = nafill;
      out[0] = id;
      for(int i = 1; i != l; ++i) {
        if(x[i] == NA_INTEGER) {
          ++id;
          prev = nafill;
        } else {
          if(x[i] - prev != del) ++id;
          prev = x[i];
        }
        out[i] = id;
      }
    }
  } else {
    IntegerVector oo = o;
    if(oo.size() != l) stop("length(o) must match length(x)");
    int val(oo[0]-1);
    if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
    if(na_skip) {
      int j = 0, end = l-1;
      if(check_o) {
        while(x[val] == NA_INTEGER && j != end) {
          out[val] = NA_INTEGER;
          val = oo[++j]-1;
          if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
        }
        if(j != end) {
          prev = x[val];
          out[val] = id;
          for(int i = j+1; i != l; ++i) {
            val = oo[i]-1;
            if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
            if(x[val] != NA_INTEGER) {
              if(x[val] - prev != del) ++id; // x[i]-x[i-1]?
              prev = x[val];
              out[val] = id;
            } else {
              out[val] = NA_INTEGER;
              if(skip_seq) prev += del;
            }
          }
        }
      } else {
        while(x[val] == NA_INTEGER && j != end) {
          out[val] = NA_INTEGER;
          val = oo[++j]-1;
        }
        if(j != end) {
          prev = x[val];
          out[val] = id;
          for(int i = j+1; i != l; ++i) {
            val = oo[i]-1;
            if(x[val] != NA_INTEGER) {
              if(x[val] - prev != del) ++id; // x[i]-x[i-1]?
              prev = x[val];
              out[val] = id;
            } else {
              out[val] = NA_INTEGER;
              if(skip_seq) prev += del;
            }
          }
        }
      }
    } else {
      int nafill = INT_MAX - 1e7;
      prev = x[val];
      if(prev == NA_INTEGER) prev = nafill;
      out[val] = id; // faster than iterator ?
      if(check_o) {
        for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
          val = oo[i]-1;
          if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
          if(x[val] == NA_INTEGER) {
            ++id;
            prev = nafill;
          } else {
            if(x[val] - prev != del) ++id;
            prev = x[val];
          }
          out[val] = id;
        }
      } else {
        for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
          val = oo[i]-1;
          if(x[val] == NA_INTEGER) {
            ++id;
            prev = nafill;
          } else {
            if(x[val] - prev != del) ++id;
            prev = x[val];
          }
          out[val] = id;
        }
      }
    }
  }
  out.attr("N.groups") = id - start + 1;
  if(start == 1) Rf_classgets(out, na_skip ? CharacterVector::create("qG") : CharacterVector::create("qG", "na.included"));
  return out;
}

// TODO: Make unique argument and generalize to all vector input types !! Or starts ?? -> Nah, GRP already does that. need to think harder. First publish without..
// The problem with groups or starts is also that you either have to dynamically fill a vector or do a second iteration...
// Rather have it process starts attribute from radixorder...


template <int RTYPE>
IntegerVector groupidImpl(Vector<RTYPE> x, SEXP o, int start, bool na_skip, bool check_o) {
  int l = x.size(), id = start;
  if(l < 1) return IntegerVector(0); // Prevents seqfault for numeric(0) #101

  typedef typename Rcpp::traits::storage_type<RTYPE>::type storage_t;
  auto isnanT = (RTYPE == REALSXP) ? [](storage_t x) { return x != x; } :
    [](storage_t x) { return x == Vector<RTYPE>::get_na(); };

    storage_t prev;
    IntegerVector out = no_init_vector(l);
    if(Rf_isNull(o)) {
      if(na_skip) {
        int j = 0, end = l;
        while(isnanT(x[j]) && j != end) out[j++] = NA_INTEGER;
        if(j != end) {
          prev = x[j];
          out[j] = id;
          for(int i = j+1; i != l; ++i) {
            if(!isnanT(x[i])) {
              if(x[i] != prev) {
                ++id;
                prev = x[i];
              }
              out[i] = id;
            } else out[i] = NA_INTEGER;
          }
        }
      } else {
        prev = x[0];
        out[0] = id;
        if(RTYPE == REALSXP) {
          for(int i = 1; i != l; ++i) {
            if(x[i] != prev) {
              if(!(prev != prev && isnanT(x[i]))) ++id;
              prev = x[i];
            }
            out[i] = id;
          }
        } else {
          for(int i = 1; i != l; ++i) {
            if(x[i] != prev) {
              ++id;
              prev = x[i];
            }
            out[i] = id;
          }
        }
      }
    } else {
      IntegerVector oo = o;
      if(oo.size() != l) stop("length(o) must match length(x)");
      int val(oo[0]-1);
      if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
      if(na_skip) {
        int j = 0, end = l-1;
        if(check_o) {
          while(isnanT(x[val]) && j != end) {
            out[val] = NA_INTEGER;
            val = oo[++j]-1;
            if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
          }
          if(j != end) {
            prev = x[val];
            out[val] = id;
            for(int i = j+1; i != l; ++i) {
              val = oo[i]-1;
              if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
              if(!isnanT(x[val])) {
                if(x[val] != prev) {
                  ++id;
                  prev = x[val];
                }
                out[val] = id;
              } else out[val] = NA_INTEGER;
            }
          }
        } else {
          while(isnanT(x[val]) && j != end) {
            out[val] = NA_INTEGER;
            val = oo[++j]-1;
          }
          if(j != end) {
            prev = x[val];
            out[val] = id;
            for(int i = j+1; i != l; ++i) {
              val = oo[i]-1;
              if(!isnanT(x[val])) {
                if(x[val] != prev) {
                  ++id;
                  prev = x[val];
                }
                out[val] = id;
              } else out[val] = NA_INTEGER;
            }
          }
        }
      } else {
        prev = x[val];
        out[val] = id; // faster than iterator ?
        if(RTYPE == REALSXP) {
          if(check_o) {
            for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
              val = oo[i]-1;
              if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
              if(x[val] != prev) {
                if(!(prev != prev && isnanT(x[val]))) ++id;
                prev = x[val];
              }
              out[val] = id;
            }
          } else {
            for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
              val = oo[i]-1;
              if(x[val] != prev) {
                if(!(prev != prev && isnanT(x[val]))) ++id;
                prev = x[val];
              }
              out[val] = id;
            }
          }
        } else {
          if(check_o) {
            for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
              val = oo[i]-1;
              if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
              if(x[val] != prev) {
                ++id;
                prev = x[val];
              }
              out[val] = id;
            }
          } else {
            for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
              val = oo[i]-1;
              if(x[val] != prev) {
                ++id;
                prev = x[val];
              }
              out[val] = id;
            }
          }
        }
      }
    }
    out.attr("N.groups") = id - start + 1;
    if(start == 1) Rf_classgets(out, na_skip ? CharacterVector::create("qG") : CharacterVector::create("qG", "na.included"));
    return out;
}


template <>
IntegerVector groupidImpl(Vector<CPLXSXP> x, SEXP o, int start, bool na_skip, bool check_o) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector groupidImpl(Vector<VECSXP> x, SEXP o, int start, bool na_skip, bool check_o) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector groupidImpl(Vector<RAWSXP> x, SEXP o, int start, bool na_skip, bool check_o) {
  stop("Not supported SEXP type!");
}

template <>
IntegerVector groupidImpl(Vector<EXPRSXP> x, SEXP o, int start, bool na_skip, bool check_o) {
  stop("Not supported SEXP type!");
}

// [[Rcpp::export]]
IntegerVector groupid(const SEXP& x, const SEXP& o = R_NilValue, int start = 1,
                      bool na_skip = false, bool check_o = true) {
  RCPP_RETURN_VECTOR(groupidImpl, x, o, start, na_skip, check_o);
}



// Integer Version
// // [[Rcpp::export]]
// IntegerVector groupid(const IntegerVector& x, const SEXP& o = R_NilValue, int start = 1,
//                       bool na_skip = false, bool check_o = true) {
//   int l = x.size(), prev, id = start;
//   IntegerVector out = no_init_vector(l);
//   if(Rf_isNull(o)) {
//     if(na_skip) {
//       int j = 0, end = l-1;
//       while(x[j] == NA_INTEGER && j != end) out[j++] = NA_INTEGER;
//       if(j != end) {
//         prev = x[j];
//         out[j] = id;
//         for(int i = j+1; i != l; ++i) {
//           if(x[i] != NA_INTEGER) {
//             if(x[i] != prev) {
//               ++id;
//               prev = x[i];
//             }
//             out[i] = id;
//           } else out[i] = NA_INTEGER;
//         }
//       }
//     } else {
//       prev = x[0];
//       out[0] = id;
//       for(int i = 1; i != l; ++i) {
//         if(x[i] != prev) {
//           ++id;
//           prev = x[i];
//         }
//         out[i] = id;
//       }
//     }
//   } else {
//     IntegerVector oo = o;
//     int val(oo[0]-1);
//     if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
//     if(na_skip) {
//       int j = 0, end = l-1;
//       if(check_o) {
//         while(x[val] == NA_INTEGER && j != end) {
//           out[val] = NA_INTEGER;
//           val = oo[++j]-1;
//           if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
//         }
//         if(j != end) {
//           prev = x[val];
//           out[val] = id;
//           for(int i = j+1; i != l; ++i) {
//             val = oo[i]-1;
//             if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
//             if(x[val] != NA_INTEGER) {
//               if(x[val] != prev) {
//                 ++id;
//                 prev = x[val];
//               }
//               out[val] = id;
//             } else out[val] = NA_INTEGER;
//           }
//         }
//       } else {
//         while(x[val] == NA_INTEGER && j != end) {
//           out[val] = NA_INTEGER;
//           val = oo[++j]-1;
//         }
//         if(j != end) {
//           prev = x[val];
//           out[val] = id;
//           for(int i = j+1; i != l; ++i) {
//             val = oo[i]-1;
//             if(x[val] != NA_INTEGER) {
//               if(x[val] != prev) {
//                 ++id;
//                 prev = x[val];
//               }
//               out[val] = id;
//             } else out[val] = NA_INTEGER;
//           }
//         }
//       }
//     } else {
//       prev = x[val];
//       out[val] = id; // faster than iterator ??
//       if(check_o) {
//         for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
//           val = oo[i]-1;
//           if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
//           if(x[val] != prev) {
//             ++id;
//             prev = x[val];
//           }
//           out[val] = id;
//         }
//       } else {
//         for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
//           val = oo[i]-1;
//           if(x[val] != prev) {
//             ++id;
//             prev = x[val];
//           }
//           out[val] = id;
//         }
//       }
//     }
//   }
//   out.attr("N.groups") = id;
//   out.attr("class") = na_skip ? "qG" : CharacterVector::create("qG", "na.included");
//   return out;
// }
//
// Simple first versions
// // // [[Rcpp::export]]
// IntegerVector groupid(const IntegerVector& x, const SEXP& o = R_NilValue, bool check = true) {
//   int l = x.size(), prev, id = 1;
//   IntegerVector out = no_init_vector(l);
//   if(Rf_isNull(o)) {
//     prev = x[0];
//     out[0] = 1;
//     for(int i = 1; i != l; ++i) {
//       if(x[i] != prev) {
//         ++id;
//         prev = x[i];
//       }
//       out[i] = id;
//     }
//   } else {
//     IntegerVector oo = o;
//     int val(oo[0]-1);
//     prev = x[val]; // https://stackoverflow.com/questions/776624/whats-faster-iterating-an-stl-vector-with-vectoriterator-or-with-at
//     out[val] = 1; // faster than iterator ??
//     if(check) {
//       for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
//         val = oo[i]-1;
//         if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
//         if(x[val] != prev) {
//           ++id;
//           prev = x[val];
//         }
//         out[val] = id;
//       }
//     } else {
//       for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
//         val = oo[i]-1;
//         if(x[val] != prev) {
//           ++id;
//           prev = x[val];
//         }
//         out[val] = id;
//       }
//     }
//   }
//   out.attr("N.groups") = id;
//   out.attr("class") =  CharacterVector::create("qG","na.included");
//   return out;
// }

//
// // [[Rcpp::export]]
// IntegerVector groupid(const IntegerVector& x, const SEXP& o = R_NilValue, bool check = true) {
//   int l = x.size(), prev, id = 1;
//   IntegerVector out = no_init_vector(l);
//   if(Rf_isNull(o)) {
//     prev = x[0];
//     out[0] = 1;
//     for(int i = 1; i != l; ++i) {
//       if(x[i] != prev) {
//         ++id;
//         prev = x[i];
//       }
//       out[i] = id;
//     }
//   } else {
//     IntegerVector oo = o;
//     int val(oo[0]-1);
//     prev = x[val]; // https://stackoverflow.com/questions/776624/whats-faster-iterating-an-stl-vector-with-vectoriterator-or-with-at
//     out[val] = 1; // faster than iterator ??
//     if(check) {
//       for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
//         val = oo[i]-1;
//         if(val < 0 || val >= l) stop("o out of allowed range [1, length(x)]");
//         if(x[val] != prev) {
//           ++id;
//           prev = x[val];
//         }
//         out[val] = id;
//       }
//     } else {
//       for(int i = 1; i != l; ++i) { //   for(IntegerVector::iterator it = oo.begin()+1, end = oo.end(); it != end; ++it) { val = *it-1;
//         val = oo[i]-1;
//         if(x[val] != prev) {
//           ++id;
//           prev = x[val];
//         }
//         out[val] = id;
//       }
//     }
//   }
//   out.attr("N.groups") = id;
//   out.attr("class") =  CharacterVector::create("qG","na.included");
//   return out;
// }
