// #include <Rcpp.h>
// using namespace Rcpp;
//
//
// // [[Rcpp::export]]
// NumericVector fkurtosisCpp(const NumericVector& x, int ng = 0, const IntegerVector& g = 0,
//                         const SEXP& w = R_NilValue, bool narm = true) {
//   int l = x.size();
//   if(l == 1) return NumericVector::create(NA_REAL);
//
//   if(Rf_isNull(w)) { // No weights
//     if(ng == 0) {
//       if(narm) {
//         int j = l-1;
//         double n = 0;
//         long double mean = 0, d1 = 0, dn = 0, dn2 = 0, term1 = 0, M2 = 0, M3 = 0, M4 = 0;
//           while(std::isnan(x[j]) && j!=0) --j;
//           if(j != 0) {
//           for(int i = j+1; i--; ) {
//             if(std::isnan(x[i])) continue;
//             d1 = x[i]-mean;
//             dn = d1 * (1 / ++n);
//             mean += dn;
//             dn2 = dn * dn;
//             term1 = d1 * dn * (n-1);
//             M4 += term1*dn2*(n*n - 3*n + 3) + 6*dn2*M2 - 4*dn*M3;
//             M3 += term1*dn*(n - 2) - 3*dn*M2;
//             M2 += term1;
//           }
//           M4 = (n*M4)/(M2*M2); // kurtosis // Excess kurtosis: - 3;
//           if(std::isnan(M4)) M4 = NA_REAL;
//           return NumericVector::create((double)M4);
//         } else return NumericVector::create(NA_REAL);
//       } else {
//         double n = 0;
//         long double mean = 0, d1 = 0, dn = 0, dn2 = 0, term1 = 0, M2 = 0, M3 = 0, M4 = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(x[i])) {
//             return NumericVector::create(NA_REAL);
//           } else {
//             d1 = x[i]-mean;
//             dn = d1 * (1 / ++n);
//             mean += dn;
//             dn2 = dn * dn;
//             term1 = d1 * dn * (n-1);
//             M4 += term1*dn2*(n*n - 3*n + 3) + 6*dn2*M2 - 4*dn*M3;
//             M3 += term1*dn*(n - 2) - 3*dn*M2;
//             M2 += term1;
//           }
//         }
//         M4 = (n*M4)/(M2*M2); // kurtosis // Excess kurtosis: - 3;
//         if(std::isnan(M4)) M4 = NA_REAL;
//         return NumericVector::create((double)M4);
//       }
//     } else { // with groups
//       if(g.size() != l) stop("length(g) must match nrow(X)");
//       long double d1 = 0;
//       if(narm) {
//         double d1 = 0, dn = 0, dn2 = 0, term1 = 0;
//         int k = 0;
// // TODO:: Organize this code like fbstats (i.e. wrap weights inside !!!)
//         NumericVector M2(ng, NA_REAL);
//         NumericVector mean = no_init_vector(ng);
//         NumericVector M3 = no_init_vector(ng);
//         NumericVector M4 = no_init_vector(ng);
//         NumericVector n(ng, 1.0);
//         for(int i = l; i--; ) {
//           if(std::isnan(x[i])) continue;
//           if(std::isnan(M2[g[i]-1])) {
//             mean[g[i]-1] = x[i];
//             M2[g[i]-1] = 0;
//           } else {
//             d1 = x[i]-mean[g[i]-1];
//             mean[g[i]-1] += d1 * (1 / ++n[g[i]-1]);
//             M2[g[i]-1] += d1*(x[i]-mean[g[i]-1]);
//           }
//         }
//         if(sd) {
//           for(int i = ng; i--; ) {
//             if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             else {
//               M2[i] = sqrt(M2[i]/(n[i]-1));
//               if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             }
//           }
//         } else {
//           for(int i = ng; i--; ) {
//             if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             else {
//               M2[i] /= n[i]-1;
//               if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             }
//           }
//         }
//         DUPLICATE_ATTRIB(M2, x);
//         return M2;
//       } else {
//         NumericVector M2(ng), mean(ng), n(ng);
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(M2[g[i]-1])) continue;
//           if(std::isnan(x[i])) {
//             M2[g[i]-1] = NA_REAL;
//             ++ngs;
//             if(ngs == ng) break;
//           } else {
//             d1 = x[i]-mean[g[i]-1];
//             mean[g[i]-1] += d1 * (1 / ++n[g[i]-1]);
//             M2[g[i]-1] += d1*(x[i]-mean[g[i]-1]);
//           }
//         }
//         if(sd) {
//           for(int i = ng; i--; ) {
//             if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             else {
//               M2[i] = sqrt(M2[i]/(n[i]-1));
//               if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             }
//           }
//         } else {
//           for(int i = ng; i--; ) {
//             if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             else {
//               M2[i] /= n[i]-1;
//               if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             }
//           }
//         }
//         DUPLICATE_ATTRIB(M2, x);
//         return M2;
//       }
//     }
//   } else { // With weights
//     NumericVector wg = w;
//     if(l != wg.size()) stop("length(w) must match length(x)");
//     if(ng == 0) {
//       if(narm) {
//         int j = l-1;
//         long double sumw = 0, mean = 0, M2 = 0, d1 = 0;
//         while((std::isnan(x[j]) || std::isnan(wg[j])) && j!=0) --j;
//         if(j != 0) {
//           for(int i = j+1; i--; ) {
//             if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
//             sumw += wg[i];
//             d1 = x[i] - mean;
//             mean += d1 * (wg[i] / sumw);
//             M2 += wg[i] * d1 * (x[i] - mean);
//           }
//           M2 /= sumw-1;
//           if(sd) M2 = sqrt(M2);
//           if(std::isnan(M2)) M2 = NA_REAL;
//           return NumericVector::create((double)M2);
//         } else return NumericVector::create(NA_REAL);
//       } else {
//         long double sumw = 0, mean = 0, M2 = 0, d1 = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(x[i]) || std::isnan(wg[i])) {
//             return NumericVector::create(NA_REAL);
//           } else {
//             sumw += wg[i];
//             d1 = x[i] - mean;
//             mean += d1 * (wg[i] / sumw);
//             M2 += wg[i] * d1 * (x[i] - mean);
//           }
//         }
//         M2 /= sumw-1;
//         if(sd) M2 = sqrt(M2);
//         if(std::isnan(M2)) M2 = NA_REAL;
//         return NumericVector::create((double)M2);
//       }
//     } else { // with groups
//       if(g.size() != l) stop("length(g) must match nrow(X)");
//       long double d1 = 0;
//       if(narm) {
//         NumericVector M2(ng, NA_REAL);
//         NumericVector sumw = no_init_vector(ng);
//         NumericVector mean = no_init_vector(ng);
//         for(int i = l; i--; ) {
//           if(std::isnan(x[i]) || std::isnan(wg[i])) continue;
//           if(std::isnan(M2[g[i]-1])) {
//             sumw[g[i]-1] = wg[i];
//             mean[g[i]-1] = x[i];
//             M2[g[i]-1] = 0;
//           } else {
//             sumw[g[i]-1] += wg[i];
//             d1 = x[i] - mean[g[i]-1];
//             mean[g[i]-1] += d1 * (wg[i] / sumw[g[i]-1]);
//             M2[g[i]-1] += wg[i] * d1 * (x[i] - mean[g[i]-1]);
//           }
//         }
//         if(sd) {
//           for(int i = ng; i--; ) {
//             if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             else {
//               M2[i] = sqrt(M2[i]/(sumw[i]-1));
//               if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             }
//           }
//         } else {
//           for(int i = ng; i--; ) {
//             if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             else {
//               M2[i] /= sumw[i]-1;
//               if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             }
//           }
//         }
//         DUPLICATE_ATTRIB(M2, x);
//         return M2;
//       } else {
//         NumericVector M2(ng), sumw(ng), mean(ng);
//         int ngs = 0;
//         for(int i = 0; i != l; ++i) {
//           if(std::isnan(M2[g[i]-1])) continue;
//           if(std::isnan(x[i]) || std::isnan(wg[i])) {
//             M2[g[i]-1] = NA_REAL;
//             ++ngs;
//             if(ngs == ng) break;
//           } else {
//             sumw[g[i]-1] += wg[i];
//             d1 = x[i] - mean[g[i]-1];
//             mean[g[i]-1] += d1 * (wg[i] / sumw[g[i]-1]);
//             M2[g[i]-1] += wg[i] * d1 * (x[i] - mean[g[i]-1]);
//           }
//         }
//         if(sd) {
//           for(int i = ng; i--; ) {
//             if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             else {
//               M2[i] = sqrt(M2[i]/(sumw[i]-1));
//               if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             }
//           }
//         } else {
//           for(int i = ng; i--; ) {
//             if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             else {
//               M2[i] /= sumw[i]-1;
//               if(std::isnan(M2[i])) M2[i] = NA_REAL;
//             }
//           }
//         }
//         DUPLICATE_ATTRIB(M2, x);
//         return M2;
//       }
//     }
//   }
// }
//
//
//
//
//
// // [[Rcpp::export]]
// SEXP fkurtosismCpp(const NumericMatrix& x, int ng = 0, const IntegerVector& g = 0,
//                   const SEXP& w = R_NilValue, bool narm = true, bool drop = true) {
//   int l = x.nrow(), col = x.ncol();
//
//   if(Rf_isNull(w)) { // No weights
//     if(ng == 0) {
//       NumericVector out = no_init_vector(col);
//       if(narm) {
//         for(int j = col; j--; ) {
//           NumericMatrix::ConstColumn column = x( _ , j);
//           int k = l-1;
//           double ni = 0;
//           long double meani = 0, d1i = 0, M2i = 0;
//           while(std::isnan(column[k]) && k!=0) --k;
//           if(k != 0) {
//             for(int i = k+1; i--; ) {
//               if(std::isnan(column[i])) continue;
//               d1i = column[i]-meani;
//               meani += d1i * (1 / ++ni);
//               M2i += d1i*(column[i]-meani);
//             }
//             M2i /= ni-1;
//             if(sd) M2i = sqrt(M2i);
//             if(std::isnan(M2i)) M2i = NA_REAL;
//             out[j] = (double)M2i;
//           } else out[j] = NA_REAL;
//         }
//       } else {
//         for(int j = col; j--; ) {
//           NumericMatrix::ConstColumn column = x( _ , j);
//           double ni = 0;
//           long double meani = 0, d1i = 0, M2i = 0;
//           for(int i = 0; i != l; ++i) {
//             if(std::isnan(column[i])) {
//               M2i = NA_REAL;
//               break;
//             } else {
//               d1i = column[i]-meani;
//               meani += d1i * (1 / ++ni);
//               M2i += d1i*(column[i]-meani);
//             }
//           }
//           M2i /= l-1;
//           if(sd) M2i = sqrt(M2i);
//           if(std::isnan(M2i)) M2i = NA_REAL;
//           out[j] = (double)M2i;
//         }
//       }
//       if(drop) out.attr("names") = colnames(x);
//       else {
//         out.attr("dim") = Dimension(1, col);
//         colnames(out) = colnames(x);
//       }
//       return out;
//     } else { // with groups
//       if(g.size() != l) stop("length(g) must match nrow(X)");
//       if(narm) {
//         NumericMatrix M2 = no_init_matrix(ng, col);
//         std::fill(M2.begin(), M2.end(), NA_REAL);
//         for(int j = col; j--; ) {
//           NumericMatrix::ConstColumn column = x( _ , j);
//           NumericMatrix::Column M2j = M2( _ , j);
//           double meanj[ng], nj[ng], d1j = 0;
//           for(int i = l; i--; ) {
//             if(std::isnan(column[i])) continue;
//             if(std::isnan(M2j[g[i]-1])) {
//               meanj[g[i]-1] = column[i];
//               M2j[g[i]-1] = 0;
//               nj[g[i]-1] = 1;
//             } else {
//               d1j = column[i]-meanj[g[i]-1];
//               meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
//               M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
//             }
//           }
//           if(sd) {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] = sqrt(M2j[i]/(nj[i]-1));
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           } else {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] /= nj[i]-1;
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           }
//         }
//         colnames(M2) = colnames(x);
//         return M2;
//       } else {
//         NumericMatrix M2(ng, col);
//         for(int j = col; j--; ) {
//           NumericMatrix::ConstColumn column = x( _ , j);
//           NumericMatrix::Column M2j = M2( _ , j);
//           std::vector<double> meanj(ng), nj(ng);
//           double d1j = 0;
//           int ngs = 0;
//           for(int i = 0; i != l; ++i) {
//             if(std::isnan(M2j[g[i]-1])) continue;
//             if(std::isnan(column[i])) {
//               M2j[g[i]-1] = NA_REAL;
//               ++ngs;
//               if(ngs == ng) break;
//             } else {
//               d1j = column[i]-meanj[g[i]-1];
//               meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
//               M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
//             }
//           }
//           if(sd) {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] = sqrt(M2j[i]/(nj[i]-1));
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           } else {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] /= nj[i]-1;
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           }
//         }
//         colnames(M2) = colnames(x);
//         return M2;
//       }
//     }
//   } else { // With weights
//     NumericVector wg = w;
//     if(l != wg.size()) stop("length(w) must match nrow(X)");
//     if(ng == 0) {
//       NumericVector out = no_init_vector(col);
//       if(narm) {
//         for(int j = col; j--; ) {
//           NumericMatrix::ConstColumn column = x( _ , j);
//           int k = l-1;
//           long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
//           while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
//           if(k != 0) {
//             for(int i = k+1; i--; ) {
//               if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
//               sumwi += wg[i];
//               d1i = column[i] - meani;
//               meani += d1i * (wg[i] / sumwi);
//               M2i += wg[i] * d1i * (column[i] - meani);
//             }
//             M2i /= sumwi-1;
//             if(sd) M2i = sqrt(M2i);
//             if(std::isnan(M2i)) M2i = NA_REAL;
//             out[j] = (double)M2i;
//           } else out[j] = NA_REAL;
//         }
//       } else {
//         for(int j = col; j--; ) {
//           NumericMatrix::ConstColumn column = x( _ , j);
//           long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
//           for(int i = 0; i != l; ++i) {
//             if(std::isnan(column[i]) || std::isnan(wg[i])) {
//               M2i = NA_REAL;
//               break;
//             } else {
//               sumwi += wg[i];
//               d1i = column[i] - meani;
//               meani += d1i * (wg[i] / sumwi);
//               M2i += wg[i] * d1i * (column[i] - meani);
//             }
//           }
//           M2i /= sumwi-1;
//           if(sd) M2i = sqrt(M2i);
//           if(std::isnan(M2i)) M2i = NA_REAL;
//           out[j] = (double)M2i;
//         }
//       }
//       if(drop) out.attr("names") = colnames(x);
//       else {
//         out.attr("dim") = Dimension(1, col);
//         colnames(out) = colnames(x);
//       }
//       return out;
//     } else { // with groups and weights
//       if(g.size() != l) stop("length(g) must match nrow(X)");
//       if(narm) {
//         NumericMatrix M2 = no_init_matrix(ng, col);
//         std::fill(M2.begin(), M2.end(), NA_REAL);
//         for(int j = col; j--; ) {
//           NumericMatrix::ConstColumn column = x( _ , j);
//           NumericMatrix::Column M2j = M2( _ , j);
//           double meanj[ng], sumwj[ng], d1j = 0;
//           for(int i = l; i--; ) {
//             if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
//             if(std::isnan(M2j[g[i]-1])) {
//               sumwj[g[i]-1] = wg[i];
//               meanj[g[i]-1] = column[i];
//               M2j[g[i]-1] = 0;
//             } else {
//               sumwj[g[i]-1] += wg[i];
//               d1j = column[i] - meanj[g[i]-1];
//               meanj[g[i]-1] += d1j * (wg[i] / sumwj[g[i]-1]);
//               M2j[g[i]-1] += wg[i] * d1j * (column[i] - meanj[g[i]-1]);
//             }
//           }
//           if(sd) {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] = sqrt(M2j[i]/(sumwj[i]-1));
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           } else {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] /= sumwj[i]-1;
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           }
//         }
//         colnames(M2) = colnames(x);
//         return M2;
//       } else {
//         NumericMatrix M2(ng, col);
//         for(int j = col; j--; ) {
//           NumericMatrix::ConstColumn column = x( _ , j);
//           NumericMatrix::Column M2j = M2( _ , j);
//           std::vector<double> meanj(ng), sumwj(ng);
//           double d1j = 0;
//           int ngs = 0;
//           for(int i = 0; i != l; ++i) {
//             if(std::isnan(M2j[g[i]-1])) continue;
//             if(std::isnan(column[i]) || std::isnan(wg[i])) {
//               M2j[g[i]-1] = NA_REAL;
//               ++ngs;
//               if(ngs == ng) break;
//             } else {
//               sumwj[g[i]-1] += wg[i];
//               d1j = column[i] - meanj[g[i]-1];
//               meanj[g[i]-1] += d1j * (wg[i] / sumwj[g[i]-1]);
//               M2j[g[i]-1] += wg[i] * d1j * (column[i] - meanj[g[i]-1]);
//             }
//           }
//           if(sd) {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] = sqrt(M2j[i]/(sumwj[i]-1));
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           } else {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] /= sumwj[i]-1;
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           }
//         }
//         colnames(M2) = colnames(x);
//         return M2;
//       }
//     }
//   }
// }
//
//
//
//
// // [[Rcpp::export]]
// SEXP fkurtosislCpp(const List& x, int ng = 0, const IntegerVector& g = 0,
//                 const SEXP& w = R_NilValue, bool narm = true, bool drop = true) {
//   int l = x.size();
//
//   if(Rf_isNull(w)) { // No weights
//     if(ng == 0) {
//       NumericVector out(l);
//       if(narm) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int k = column.size()-1;
//           double ni = 0;
//           long double meani = 0, d1i = 0, M2i = 0;
//           while(std::isnan(column[k]) && k!=0) --k;
//           if(k != 0) {
//             for(int i = k+1; i--; ) {
//               if(std::isnan(column[i])) continue;
//               d1i = column[i]-meani;
//               meani += d1i * (1 / ++ni);
//               M2i += d1i*(column[i]-meani);
//             }
//             M2i /= ni-1;
//             if(sd) M2i = sqrt(M2i);
//             if(std::isnan(M2i)) M2i = NA_REAL;
//             out[j] = (double)M2i;
//           } else out[j] = NA_REAL;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           int row = column.size();
//           double ni = 0;
//           long double meani = 0, d1i = 0, M2i = 0;
//           for(int i = 0; i != row; ++i) {
//             if(std::isnan(column[i])) {
//               M2i = NA_REAL;
//               break;
//             } else {
//               d1i = column[i]-meani;
//               meani += d1i * (1 / ++ni);
//               M2i += d1i*(column[i]-meani);
//             }
//           }
//           M2i /= row-1;
//           if(sd) M2i = sqrt(M2i);
//           if(std::isnan(M2i)) M2i = NA_REAL;
//           out[j] = (double)M2i;
//         }
//       }
//       if(drop) {
//         out.attr("names") = x.attr("names");
//         return out;
//       } else {
//         List res(l);
//         for(int j = l; j--; ) {
//           res[j] = out[j];
//           SHALLOW_DUPLICATE_ATTRIB(res[j], x[j]);
//         }
//         DUPLICATE_ATTRIB(res, x);
//         res.attr("row.names") = 1;
//         return res;
//       }
//     } else { // With groups
//       List out(l);
//       int gss = g.size();
//       if(narm) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           if(gss != column.size()) stop("length(g) must match nrow(X)");
//           NumericVector M2j(ng, NA_REAL);
//           double meanj[ng], d1j = 0;
//           std::vector<double> nj(ng, 1.0);
//           for(int i = gss; i--; ) {
//             if(std::isnan(column[i])) continue;
//             if(std::isnan(M2j[g[i]-1])) {
//               meanj[g[i]-1] = column[i];
//               M2j[g[i]-1] = 0;
//             } else {
//               d1j = column[i]-meanj[g[i]-1];
//               meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
//               M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
//             }
//           }
//           if(sd) {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] = sqrt(M2j[i]/(nj[i]-1));
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           } else {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] /= nj[i]-1;
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           }
//           SHALLOW_DUPLICATE_ATTRIB(M2j, column);
//           out[j] = M2j;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           if(gss != column.size()) stop("length(g) must match nrow(X)");
//           NumericVector M2j(ng);
//           std::vector<double> meanj(ng), nj(ng);
//           double d1j = 0;
//           int ngs = 0;
//           for(int i = 0; i != gss; ++i) {
//             if(std::isnan(M2j[g[i]-1])) continue;
//             if(std::isnan(column[i])) {
//               M2j[g[i]-1] = NA_REAL;
//               ++ngs;
//               if(ngs == ng) break;
//             } else {
//               d1j = column[i]-meanj[g[i]-1];
//               meanj[g[i]-1] += d1j * (1 / ++nj[g[i]-1]);
//               M2j[g[i]-1] += d1j*(column[i]-meanj[g[i]-1]);
//             }
//           }
//           if(sd) {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] = sqrt(M2j[i]/(nj[i]-1));
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           } else {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] /= nj[i]-1;
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           }
//           SHALLOW_DUPLICATE_ATTRIB(M2j, column);
//           out[j] = M2j;
//         }
//       }
//       DUPLICATE_ATTRIB(out, x);
//       out.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng); // NumericVector::create(NA_REAL, -ng);
//       return out;
//     }
//   } else { // With weights
//     NumericVector wg = w;
//     int wgs = wg.size();
//     if(ng == 0) {
//       NumericVector out(l);
//       if(narm) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           if(column.size() != wgs) stop("length(w) must match nrow(X)");
//           int k = wgs-1;
//           long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
//           while((std::isnan(column[k]) || std::isnan(wg[k])) && k!=0) --k;
//           if(k != 0) {
//             for(int i = k+1; i--; ) {
//               if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
//               sumwi += wg[i];
//               d1i = column[i] - meani;
//               meani += d1i * (wg[i] / sumwi);
//               M2i += wg[i] * d1i * (column[i] - meani);
//             }
//             M2i /= sumwi-1;
//             if(sd) M2i = sqrt(M2i);
//             if(std::isnan(M2i)) M2i = NA_REAL;
//             out[j] = (double)M2i;
//           } else out[j] = NA_REAL;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           if(column.size() != wgs) stop("length(w) must match nrow(X)");
//           long double sumwi = 0, meani = 0, M2i = 0, d1i = 0;
//           for(int i = 0; i != wgs; ++i) {
//             if(std::isnan(column[i]) || std::isnan(wg[i])) {
//               M2i = NA_REAL;
//               break;
//             } else {
//               sumwi += wg[i];
//               d1i = column[i] - meani;
//               meani += d1i * (wg[i] / sumwi);
//               M2i += wg[i] * d1i * (column[i] - meani);
//             }
//           }
//           M2i /= sumwi-1;
//           if(sd) M2i = sqrt(M2i);
//           if(std::isnan(M2i)) M2i = NA_REAL;
//           out[j] = (double)M2i;
//         }
//       }
//       if(drop) {
//         out.attr("names") = x.attr("names");
//         return out;
//       } else {
//         List res(l);
//         for(int j = l; j--; ) {
//           res[j] = out[j];
//           SHALLOW_DUPLICATE_ATTRIB(res[j], x[j]);
//         }
//         DUPLICATE_ATTRIB(res, x);
//         res.attr("row.names") = 1;
//         return res;
//       }
//     } else {
//       List out(l);
//       int gss = g.size();
//       if(wgs != gss) stop("length(w) must match length(g)");
//       if(narm) {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           if(gss != column.size()) stop("length(g) must match nrow(X)");
//           NumericVector M2j(ng, NA_REAL);
//           double sumwj[ng], meanj[ng], d1j = 0;
//           for(int i = gss; i--; ) {
//             if(std::isnan(column[i]) || std::isnan(wg[i])) continue;
//             if(std::isnan(M2j[g[i]-1])) {
//               sumwj[g[i]-1] = wg[i];
//               meanj[g[i]-1] = column[i];
//               M2j[g[i]-1] = 0;
//             } else {
//               sumwj[g[i]-1] += wg[i];
//               d1j = column[i] - meanj[g[i]-1];
//               meanj[g[i]-1] += d1j * (wg[i] / sumwj[g[i]-1]);
//               M2j[g[i]-1] += wg[i] * d1j * (column[i] - meanj[g[i]-1]);
//             }
//           }
//           if(sd) {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] = sqrt(M2j[i]/(sumwj[i]-1));
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           } else {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] /= sumwj[i]-1;
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           }
//           SHALLOW_DUPLICATE_ATTRIB(M2j, column);
//           out[j] = M2j;
//         }
//       } else {
//         for(int j = l; j--; ) {
//           NumericVector column = x[j];
//           if(gss != column.size()) stop("length(g) must match nrow(X)");
//           NumericVector M2j(ng);
//           std::vector<double> sumwj(ng), meanj(ng);
//           double d1j = 0;
//           int ngs = 0;
//           for(int i = 0; i != gss; ++i) {
//             if(std::isnan(M2j[g[i]-1])) continue;
//             if(std::isnan(column[i]) || std::isnan(wg[i])) {
//               M2j[g[i]-1] = NA_REAL;
//               ++ngs;
//               if(ngs == ng) break;
//             } else {
//               sumwj[g[i]-1] += wg[i];
//               d1j = column[i] - meanj[g[i]-1];
//               meanj[g[i]-1] += d1j * (wg[i] / sumwj[g[i]-1]);
//               M2j[g[i]-1] += wg[i] * d1j * (column[i] - meanj[g[i]-1]);
//             }
//           }
//           if(sd) {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] = sqrt(M2j[i]/(sumwj[i]-1));
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           } else {
//             for(int i = ng; i--; ) {
//               if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               else {
//                 M2j[i] /= sumwj[i]-1;
//                 if(std::isnan(M2j[i])) M2j[i] = NA_REAL;
//               }
//             }
//           }
//           SHALLOW_DUPLICATE_ATTRIB(M2j, column);
//           out[j] = M2j;
//         }
//       }
//       DUPLICATE_ATTRIB(out, x);
//       out.attr("row.names") = IntegerVector::create(NA_INTEGER, -ng); // NumericVector::create(NA_REAL, -ng);
//       return out;
//     }
//   }
// }
