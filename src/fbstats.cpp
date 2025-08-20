#include <Rcpp.h>
using namespace Rcpp;


// TODO: Still check printing (naming and setting classes) options

// inline bool isnan2(double x) {
//   return x != x;
// }

CharacterVector get_stats_names(int n, bool panel = false) {
  String N = panel ? "N/T" : "N";
  switch(n) {
    case 5: return CharacterVector::create(N,"Mean","SD","Min","Max");
    case 6: return CharacterVector::create(N,"WeightSum","Mean","SD","Min","Max");
    case 7: return CharacterVector::create(N,"Mean","SD","Min","Max","Skew","Kurt");
    case 8: return CharacterVector::create(N,"WeightSum","Mean","SD","Min","Max","Skew","Kurt");
    default: stop("length of stats names needs to be between 5 and 8");
  }
}

// use constant references on the temp function also ?
NumericVector fbstatstemp(NumericVector x, bool ext = false, int ng = 0, IntegerVector g = 0, SEXP w = R_NilValue,
                          bool setn = true, bool stable_algo = true, SEXP gn = R_NilValue) {
  int l = x.size();
  bool weights = !Rf_isNull(w);

   if(!ext) {
    if(ng == 0) { // No groups
      if(l == 1) { // need this so that qsu(1) works properly
        NumericVector result = weights ? NumericVector::create(1,Rf_asReal(w),x[0],NA_REAL,x[0],x[0]) :
                                         NumericVector::create(1,x[0],NA_REAL,x[0],x[0]);
        if(setn) {
          Rf_namesgets(result, get_stats_names(5+weights));
          Rf_classgets(result, CharacterVector::create("qsu","table"));
        }
        return result;
      }
      int j = l-1;
      // double n = 0, min = R_PosInf, max = R_NegInf;
      // long double mean = 0, d1 = 0, M2 = 0;
      double n = 0, min = R_PosInf, max = R_NegInf, mean = 0, d1 = 0, M2 = 0;
      if(!weights) { // No weights
        while(std::isnan(x[j]) && j!=0) --j;
        if(j != 0) {  // if(j == 0) stop("Not enough non-mising obs.");
          if(stable_algo) {
            for(int i = j+1; i--; ) {
              if(std::isnan(x[i])) continue;
              d1 = x[i]-mean;
              mean += d1 * (1 / ++n);
              M2 += d1*(x[i]-mean);
              if(min > x[i]) min = x[i];
              if(max < x[i]) max = x[i];
            }
            M2 = sqrt(M2/(n-1));
          } else {
            int k = 0;
            long double sum = 0, sq_sum = 0;
            for(int i = j+1; i--; ) {
              d1 = x[i];
              if(std::isnan(d1)) continue;
              sum += d1;
              sq_sum += d1 * d1;
              if(min > d1) min = d1;
              if(max < d1) max = d1;
              ++k;
            }
            sum /= k;
            sq_sum -= sum*sum*k;
            M2 = (double)sqrt(sq_sum/(k-1));
            n = (double)k;
            mean = (double)sum;
          }
        } else mean = M2 = min = max = NA_REAL;
        if(std::isnan(M2)) M2 = NA_REAL;
        NumericVector result = NumericVector::create(n,mean,M2,min,max); // NumericVector::create(n,(double)mean,(double)M2,min,max);
        if(setn) {
          Rf_namesgets(result, CharacterVector::create("N","Mean","SD","Min","Max"));
          Rf_classgets(result, CharacterVector::create("qsu","table"));
        }
        return result;
      } else { // with weights
        NumericVector wg = w;
        if(l != wg.size()) stop("length(w) must match length(x)");
        // long double sumw = 0;
        double sumw = 0;
        while((std::isnan(x[j]) || std::isnan(wg[j]) || wg[j] == 0) && j!=0) --j;
         if(j != 0) { // if(j == 0) stop("Not enough non-mising obs.");
          for(int i = j+1; i--; ) {
            if(std::isnan(x[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
            sumw += wg[i];
            d1 = x[i] - mean;
            mean += d1 * (wg[i] / sumw);
            M2 += wg[i] * d1 * (x[i] - mean);
            ++n;
            if(min > x[i]) min = x[i];
            if(max < x[i]) max = x[i];
          }
          M2 = sqrt(M2/(sumw-1));
        } else mean = M2 = min = max = NA_REAL;
        if(std::isnan(M2)) M2 = NA_REAL;
        NumericVector result = NumericVector::create(n,sumw,mean,M2,min,max); // NumericVector::create(n,(double)mean,(double)M2,min,max);
        if(setn) {
          Rf_namesgets(result, CharacterVector::create("N","WeightSum","Mean","SD","Min","Max"));
          Rf_classgets(result, CharacterVector::create("qsu","table"));
        }
        return result;
      }
    } else { // with groups
      if(g.size() != l) stop("length(g) must match nrow(X)");
      // long double d1 = 0;
      double d1 = 0;
      int k = 0;
      NumericMatrix result(ng, 5+weights); //  = no_init_matrix initializing is better -> valgrind
      NumericMatrix::Column n = result( _ , 0);
      NumericMatrix::Column mean = result( _ , 1+weights);
      NumericMatrix::Column M2 = result( _ , 2+weights);
      NumericMatrix::Column min = result( _ , 3+weights);
      NumericMatrix::Column max = result( _ , 4+weights);
      std::fill(M2.begin(), M2.end(), NA_REAL);
      if(!weights) { // No weights
        if(stable_algo) {
          for(int i = l; i--; ) {
            if(std::isnan(x[i])) continue;
            k = g[i]-1;
            if(std::isnan(M2[k])) {
              mean[k] = min[k] = max[k] = x[i];
              M2[k] = 0.0;
              n[k] = 1.0;
            } else {
              d1 = x[i]-mean[k];
              mean[k] += d1 * (1 / ++n[k]);
              M2[k] += d1*(x[i]-mean[k]);
              if(min[k] > x[i]) min[k] = x[i];
              if(max[k] < x[i]) max[k] = x[i];
            }
          }
          for(int i = ng; i--; ) if(!std::isnan(M2[i])) M2[i] = sqrt(M2[i]/(n[i]-1));
        } else {
          for(int i = l; i--; ) {
            if(std::isnan(x[i])) continue;
            k = g[i]-1;
            if(std::isnan(M2[k])) {
              mean[k] = min[k] = max[k] = x[i];
              M2[k] = pow(x[i],2);
              n[k] = 1.0;
            } else {
              mean[k] += x[i];
              M2[k] += pow(x[i],2);
              if(min[k] > x[i]) min[k] = x[i];
              if(max[k] < x[i]) max[k] = x[i];
              ++n[k];
            }
          }
          for(int i = ng; i--; ) {
            if(std::isnan(M2[i])) continue;
            mean[i] /= n[i];
            M2[i] = sqrt((M2[i] - pow(mean[i],2)*n[i])/(n[i]-1));
            if(std::isnan(M2[i])) M2[i] = NA_REAL;
          }
        }
      } else { // with weights
        NumericVector wg = w;
        if(l != wg.size()) stop("length(w) must match length(x)");
        // NumericVector sumw(ng); // = no_init_vector(ng); // better for valgrind
        NumericMatrix::Column sumw = result( _ , 1);
        for(int i = l; i--; ) {
          if(std::isnan(x[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
          k = g[i]-1;
          if(std::isnan(M2[k])) {
            sumw[k] = wg[i];
            mean[k] = min[k] = max[k] = x[i];
            M2[k] = 0.0;
            n[k] = 1.0;
          } else {
            sumw[k] += wg[i];
            d1 = x[i] - mean[k];
            mean[k] += d1 * (wg[i] / sumw[k]);
            M2[k] += wg[i] * d1 * (x[i] - mean[k]);
            ++n[k];
            if(min[k] > x[i]) min[k] = x[i];
            if(max[k] < x[i]) max[k] = x[i];
          }
        }
        for(int i = ng; i--; ) {
          if(n[i] == 0) mean[i] = min[i] = max[i] = NA_REAL;
          else if(!std::isnan(M2[i])) M2[i] = sqrt(M2[i]/(sumw[i]-1));
        }
      }
      if(setn) {
        Rf_dimnamesgets(result, List::create(gn, get_stats_names(5+weights)));
        Rf_classgets(result, CharacterVector::create("qsu","matrix","table"));
      }
      return result;
    }
   } else {
     if(ng == 0) { // No groups
       int j = l-1;
       // double n = 0, min = R_PosInf, max = R_NegInf;
       // long double mean = 0, d1 = 0, dn = 0, dn2 = 0, term1 = 0, M2 = 0, M3 = 0, M4 = 0;
       double n = 0, min = R_PosInf, max = R_NegInf, mean = 0, d1 = 0, dn = 0, dn2 = 0, term1 = 0, M2 = 0, M3 = 0, M4 = 0;
       if(!weights) { // No weights
         while(std::isnan(x[j]) && j!=0) --j;
         if(j != 0) {  // if(j == 0) stop("Not enough non-mising obs.");
         for(int i = j+1; i--; ) {
           if(std::isnan(x[i])) continue;
           d1 = x[i]-mean;
           dn = d1 * (1 / ++n);
           mean += dn;
           dn2 = dn * dn;
           term1 = d1 * dn * (n-1);
           M4 += term1*dn2*(n*n - 3*n + 3) + 6*dn2*M2 - 4*dn*M3;
           M3 += term1*dn*(n - 2) - 3*dn*M2;
           M2 += term1;
           if(min > x[i]) min = x[i];
           if(max < x[i]) max = x[i];
         }
         M4 = (n*M4)/(M2*M2); // kurtosis // Excess kurtosis: - 3;
         M3 = (sqrt(n)*M3) / sqrt(pow(M2,3)); // Skewness
         M2 = sqrt(M2/(n-1)); // Standard Deviation
         } else mean = M2 = M3 = M4 = min = max = NA_REAL;
         NumericVector result = NumericVector::create(n,mean,M2,min,max,M3,M4); // NumericVector::create(n,(double)mean,(double)M2,min,max,(double)M3,(double)M4);
         if(setn) {
           Rf_namesgets(result, CharacterVector::create("N","Mean","SD","Min","Max","Skew","Kurt"));
           Rf_classgets(result, CharacterVector::create("qsu","table"));
         }
         return result;
       } else { // with weights
         NumericVector wg = w;
         if(l != wg.size()) stop("length(w) must match length(x)");
         // long double sumw = 0;
         double sumw = 0;
         while((std::isnan(x[j]) || std::isnan(wg[j]) || wg[j] == 0) && j!=0) --j;
         if(j != 0) {  // if(j == 0) stop("Not enough non-mising obs.");
         for(int i = j+1; i--; ) {
           if(std::isnan(x[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
           sumw += wg[i];
           mean += x[i] * wg[i]; ++n;
           if(min > x[i]) min = x[i];
           if(max < x[i]) max = x[i];
         }
         mean /= sumw;
         long double M2l = 0.0, M3l = 0.0, M4l = 0.0;
         for(int i = j+1; i--; ) {
           if(std::isnan(x[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
           d1 = x[i] - mean;
           dn = d1 * d1;
           dn2 = dn * dn;
           M2l += wg[i] * dn;
           M3l += wg[i] * dn * d1;
           M4l += wg[i] * dn2;
         }
         M4 = (sumw*M4l)/(M2l*M2l); // kurtosis // Excess kurtosis: - 3;
         M3 = (sqrt(sumw)*M3l) / sqrt(pow(M2l,3)); // Skewness
         M2 = sqrt(M2l/(sumw-1)); // Standard Deviation
         } else mean = M2 = M3 = M4 = min = max = NA_REAL;
         NumericVector result = NumericVector::create(n,sumw,mean,M2,min,max,M3,M4); // NumericVector::create(n,(double)mean,(double)M2,min,max,(double)M3,(double)M4);
         if(setn) {
           Rf_namesgets(result, CharacterVector::create("N","WeightSum","Mean","SD","Min","Max","Skew","Kurt"));
           Rf_classgets(result, CharacterVector::create("qsu","table"));
         }
         return result;
       }
     } else { // with groups
       if(g.size() != l) stop("length(g) must match nrow(X)");
       double d1 = 0, dn = 0, dn2 = 0, term1 = 0;
       int k = 0;
       NumericMatrix result(ng, 7+weights); //  = no_init_matrix // Initializing better -> valgrind
       NumericMatrix::Column n = result( _ , 0);
       NumericMatrix::Column mean = result( _ , 1+weights);
       NumericMatrix::Column M2 = result( _ , 2+weights);
       NumericMatrix::Column min = result( _ , 3+weights);
       NumericMatrix::Column max = result( _ , 4+weights);
       NumericMatrix::Column M3 = result( _ , 5+weights);
       NumericMatrix::Column M4 = result( _ , 6+weights);
       std::fill(M2.begin(), M2.end(), NA_REAL);
       if(!weights) { // No weights
         for(int i = l; i--; ) {
           if(std::isnan(x[i])) continue;
           k = g[i]-1;
           if(std::isnan(M2[k])) {
             mean[k] = min[k] = max[k] = x[i];
             M2[k] = M3[k] = M4[k] = 0.0;
             n[k] = 1.0;
           } else {
             d1 = x[i]-mean[k];
             dn = d1 * (1 / ++n[k]);
             mean[k] += dn;
             dn2 = dn * dn;
             term1 = d1 * dn * (n[k]-1);
             M4[k] += term1*dn2*(n[k]*n[k] - 3*n[k] + 3) + 6*dn2*M2[k] - 4*dn*M3[k];
             M3[k] += term1*dn*(n[k] - 2) - 3*dn*M2[k];
             M2[k] += term1;
             if(min[k] > x[i]) min[k] = x[i];
             if(max[k] < x[i]) max[k] = x[i];
           }
         }
         for(int i = ng; i--; ) {
           M4[i] = (n[i]*M4[i])/(M2[i]*M2[i]); // kurtosis // Excess kurtosis: - 3;
           M3[i] = (sqrt(n[i])*M3[i]) / sqrt(pow(M2[i],3)); // Skewness
           M2[i] = sqrt(M2[i]/(n[i]-1)); // Standard Deviation
         }
       } else { // with weights
         NumericVector wg = w;
         if(l != wg.size()) stop("length(w) must match length(x)");
         // NumericVector sumw(ng); // = no_init_vector(ng); // better for valgrind
         NumericMatrix::Column sumw = result( _ , 1);
         for(int i = l; i--; ) {
           if(std::isnan(x[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
           k = g[i]-1;
           if(std::isnan(M2[k])) {
             sumw[k] = wg[i];
             mean[k] = min[k] = max[k] = x[i];
             M2[k] = M3[k] = M4[k] = 0.0;
             n[k] = 1.0;
           } else {
             sumw[k] += wg[i];
             mean[k] += (x[i] - mean[k]) * (wg[i] / sumw[k]);
             ++n[k];
             if(min[k] > x[i]) min[k] = x[i];
             if(max[k] < x[i]) max[k] = x[i];
           }
         }
         for(int i = l; i--; ) {
           if(std::isnan(x[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
           k = g[i]-1;
           d1 = x[i] - mean[k];
           dn = d1 * d1;
           dn2 = dn * dn;
           M2[k] += wg[i] * dn;
           M3[k] += wg[i] * dn * d1;
           M4[k] += wg[i] * dn2;
         }
         for(int i = ng; i--; ) {
           if(n[i] == 0) mean[i] = min[i] = max[i] = M3[i] = M4[i] = NA_REAL;
           else {
             M4[i] = (sumw[i]*M4[i])/(M2[i]*M2[i]); // kurtosis // Excess kurtosis: - 3;
             M3[i] = (sqrt(sumw[i])*M3[i]) / sqrt(pow(M2[i],3)); // Skewness
             M2[i] = sqrt(M2[i]/(sumw[i]-1)); // Standard Deviation
           }
         }
       }
       if(setn) {
         Rf_dimnamesgets(result, List::create(gn, get_stats_names(7+weights)));
         Rf_classgets(result, CharacterVector::create("qsu","matrix","table"));
       }
       return result;
     }
   }
  // } else { // detailed summary: fully sorting. Note: This doesn't work grouped, groups must also be sorted -> need to sort within each group or compute ordering
  //   NumericVector y = no_init_vector(l);
  //   auto pend = std::remove_copy_if(x.begin(), x.end(), y.begin(), isnan2);
  //   l = pend - x.begin(); //  middle = sz/2-1;
  //   std::sort(y.begin(), pend); // good ??
  //
  //   if(dets == 1 && det[0] == 1) det = 5;
  // }
}

inline NumericVector replaceC12(NumericMatrix x, NumericVector y, bool div = false) {
  int nc = x.ncol();
  if(div) {
    NumericMatrix::Column C1 = x(_, 0); // best ?
    C1 = C1 / y;
    if(nc == 6 || nc == 8) { // WeightSum column
      NumericMatrix::Column C2 = x(_, 1);
      C2 = C2 / y;
    }
  } else {
    x(_, 0) = y; // best way ? use NumericMatrix::Column ?
    if(nc == 6 || nc == 8) { // WeightSum column
      x(_, 1) = y;
    }
  }
  return x;
}

// [[Rcpp::export]]
SEXP fbstatsCpp(const NumericVector& x, bool ext = false, int ng = 0, const IntegerVector& g = 0,
                int npg = 0, const IntegerVector& pg = 0, const SEXP& w = R_NilValue,
                bool stable_algo = true, bool array = true,
                bool setn = true, const SEXP& gn = R_NilValue) {

  if(npg == 0) { // No panel
    if(ng == 0) { // No groups
      return(fbstatstemp(x, ext, 0, 0, w, setn, stable_algo, gn));
    } else {
      return(fbstatstemp(x, ext, ng, g, w, setn, stable_algo, gn));
    }
  } else {
    int l = x.size();
    if(pg.size() != l) stop("length(pid) must match nrow(X)");
    bool weights = !Rf_isNull(w);
    int d = ((ext) ? 7 : 5) + weights;
    NumericVector sum(npg, NA_REAL);
    NumericVector sumw((weights) ? npg : 1); // no_init_vector(npg) : no_init_vector(1); // better for valgrind
    double osum = 0;

    if(!weights) {
      IntegerVector n(npg, 1);
      for(int i = l; i--; ) {
        if(!std::isnan(x[i])) {
          if(std::isnan(sum[pg[i]-1])) sum[pg[i]-1] = x[i];
          else {
            sum[pg[i]-1] += x[i];
            ++n[pg[i]-1];
          }
        }
      }
      int on = 0;
      for(int i = npg; i--; ) { // Problem: if one sum remained NA, osum becomes NA (also issue with B and W and TRA)
        if(std::isnan(sum[i])) continue; // solves the issue
        osum += sum[i];
        on += n[i];
        sum[i] /= n[i];
      }
      osum = osum/on;
    } else {
      NumericVector wg = w;
      if(l != wg.size()) stop("length(w) must match length(x)");
      // Note: Skipping zero weights is not really necessary here, but it might be numerically better and also faster if there are many.
      for(int i = l; i--; ) {
        if(std::isnan(x[i]) || std::isnan(wg[i]) || wg[i] == 0) continue;
        if(std::isnan(sum[pg[i]-1])) {
          sum[pg[i]-1] = x[i]*wg[i];
          sumw[pg[i]-1] = wg[i];
        } else {
          sum[pg[i]-1] += x[i]*wg[i];
          sumw[pg[i]-1] += wg[i];
        }
      }
      double osumw = 0;
      for(int i = npg; i--; ) {
        if(std::isnan(sum[i]) || sumw[i] == 0) continue; // solves the issue
        osum += sum[i];
        osumw += sumw[i];
        sum[i] /= sumw[i];
      }
      osum = osum/osumw;
      // for(int i = npg; i--; ) sumw[i] /= osumw;
    }

    NumericVector within = no_init_vector(l);

    if(ng == 0) { // No groups
      for(int i = 0; i != l; ++i) within[i] = x[i] - sum[pg[i]-1] + osum; // if-check for NA's is not faster
      NumericMatrix result = no_init_matrix(3, d);
      result(0, _) = fbstatstemp(x, ext, 0, 0, w, false, stable_algo);
      result(1, _) = (weights) ? fbstatstemp(sum, ext, 0, 0, sumw, false, stable_algo) : fbstatstemp(sum, ext, 0, 0, w, false, stable_algo);
      result(2, _) = fbstatstemp(within, ext, 0, 0, w, false, stable_algo);
      result[2] /= result[1];
      if(weights) {
        result[4] = result[1];
        result[5] /= result[1];
      }
      if(setn) {
        Rf_dimnamesgets(result, List::create(CharacterVector::create("Overall","Between","Within"),
                                             get_stats_names(d, true)));
        Rf_classgets(result, CharacterVector::create("qsu","matrix","table"));
      }
      return(result);
    } else {
      if(g.size() != l) stop("length(g) must match nrow(X)");
      NumericVector between = no_init_vector(l);
      // bool groupids[ng][npg]; // could do +1 trick, but that could be costly in term of memory if few g and many pg
      LogicalMatrix groupids = no_init_matrix(ng, npg);
      // memset(groupids, true, sizeof(bool)*ng*npg); // works ? necessary ?
      std::fill(groupids.begin(), groupids.end(), true);
      NumericVector gnpids(ng); // best ?
      for(int i = 0; i != l; ++i) {
        if(std::isnan(x[i])) { // important ? right ?
          between[i] = within[i] = NA_REAL; // x[i] ?
        } else {
          if(groupids(g[i]-1, pg[i]-1)) { // added this part
            ++gnpids[g[i]-1];
            groupids(g[i]-1, pg[i]-1) = false;
          }
          between[i] = sum[pg[i]-1];
          within[i] = x[i] - between[i] + osum;
        }
      }
      if(array) {
        NumericMatrix result = no_init_matrix(d*ng, 3);
        result(_,0) = fbstatstemp(x, ext, ng, g, w, false, stable_algo);
        result(_,1) = replaceC12(as<NumericMatrix>(fbstatstemp(between, ext, ng, g, w, false, stable_algo)), gnpids); // how to do this ? -> above best approach ?
        result(_,2) = replaceC12(as<NumericMatrix>(fbstatstemp(within, ext, ng, g, w, false, stable_algo)), gnpids, true);
        if(setn) {
          Rf_dimgets(result, Dimension(ng, d, 3));
          Rf_dimnamesgets(result, List::create(gn, get_stats_names(d, true), CharacterVector::create("Overall","Between","Within")));
          Rf_classgets(result, CharacterVector::create("qsu","array","table"));
        }
        return(result);
      } else {
        List result(3); // option array ?
        result[0] = fbstatstemp(x, ext, ng, g, w, true, stable_algo, gn);
        result[1] = replaceC12(as<NumericMatrix>(fbstatstemp(between, ext, ng, g, w, true, stable_algo, gn)), gnpids); // how to do this ? -> above best approach ?
        result[2] = replaceC12(as<NumericMatrix>(fbstatstemp(within, ext, ng, g, w, true, stable_algo, gn)), gnpids, true);
        Rf_namesgets(result, CharacterVector::create("Overall","Between","Within"));
        return(result);
      }
    }
  }
}


// [[Rcpp::export]]
SEXP fbstatsmCpp(const NumericMatrix& x, bool ext = false, int ng = 0, const IntegerVector& g = 0,
                 int npg = 0, const IntegerVector& pg = 0,
                 const SEXP& w = R_NilValue, bool stable_algo = true, bool array = true,
                 const SEXP& gn = R_NilValue) {
  bool weights = !Rf_isNull(w);
  int col = x.ncol(), d = ((ext) ? 7 : 5) + weights; // l = x.nrow(),

  if(npg == 0) { // No panel
    if(ng == 0) { // No groups
      NumericMatrix out = no_init_matrix(col, d);
      for(int j = col; j--; ) out(j, _) = fbstatstemp(x(_, j), ext, 0, 0, w, false, stable_algo);
      Rf_dimnamesgets(out, List::create(colnames(x), get_stats_names(d)));
      Rf_classgets(out, CharacterVector::create("qsu","matrix","table"));
      return out;
    } else {
      // if(g.size() != l) stop("length(g) must match nrow(X)"); // checked in fbstatstemp
      if(array) {
        NumericMatrix out = no_init_matrix(d*ng, col);
        for(int j = col; j--; ) out(_, j) = fbstatstemp(x(_, j), ext, ng, g, w, false, stable_algo);
        Rf_dimgets(out, Dimension(ng, d, col));
        Rf_dimnamesgets(out, List::create(gn, get_stats_names(d), colnames(x)));
        Rf_classgets(out, CharacterVector::create("qsu","array","table"));
        return out;
      } else {
        List out(col);
        for(int j = col; j--; ) out[j] = fbstatstemp(x(_, j), ext, ng, g, w, true, stable_algo, gn);
        Rf_setAttrib(out, R_NamesSymbol, colnames(x));
        return out;
      }
    }
  } else {
    if(ng == 0) {
      if(array) {
        NumericMatrix out = no_init_matrix(d*3, col);
        for(int j = col; j--; ) out(_, j) = as<NumericVector>(fbstatsCpp(x(_, j), ext, 0, 0, npg, pg, w, stable_algo, true, false)); // or Rf_coerce ?
        Rf_dimgets(out, Dimension(3, d, col));
        Rf_dimnamesgets(out, List::create(CharacterVector::create("Overall","Between","Within"),
                        get_stats_names(d, true), colnames(x)));
        Rf_classgets(out, CharacterVector::create("qsu","array","table"));
        return out;
      } else {
        List out(col);
        for(int j = col; j--; ) out[j] = fbstatsCpp(x(_, j), ext, 0, 0, npg, pg, w, stable_algo, false, true, gn);
        Rf_setAttrib(out, R_NamesSymbol, colnames(x));
        return out;
      }
    } else {
      if(array) {
        NumericMatrix out = no_init_matrix(d*3*ng, col);
        for(int j = col; j--; ) out(_, j) = as<NumericVector>(fbstatsCpp(x(_, j), ext, ng, g, npg, pg, w, stable_algo, true, false)); // or Rf_coerce ?
        Rf_dimgets(out, IntegerVector::create(ng, d, 3, col));
        Rf_dimnamesgets(out, List::create(gn, get_stats_names(d, true),
                        CharacterVector::create("Overall","Between","Within"), colnames(x)));
        Rf_classgets(out, CharacterVector::create("qsu","array","table"));
        return out;
      } else {
        List out(col);
        for(int j = col; j--; ) out[j] = fbstatsCpp(x(_, j), ext, ng, g, npg, pg, w, stable_algo, false, true, gn);
        Rf_setAttrib(out, R_NamesSymbol, colnames(x));
        return out;
      }
    }
  }
}


template <int RTYPE>
NumericVector fnobs5Impl(Vector<RTYPE> x, bool ext = false, int ng = 0, IntegerVector g = 0, SEXP w = R_NilValue, bool real = false, bool setn = false, SEXP gn = R_NilValue) {

  bool weights = !Rf_isNull(w);
  int l = x.size(), d = ((ext) ? 7 : 5) + weights;

  if(ng == 0) {
    int n = 0;
    double wsum = 0.0;
    NumericVector out(d, NA_REAL);
    if(weights) {
      NumericVector wg = w;
      if(real) {
        for(int i = 0; i != l; ++i) {
          if(x[i] == x[i] && wg[i] == wg[i] && wg[i] != 0) {
            wsum += wg[i]; ++n;
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(x[i] != Vector<RTYPE>::get_na() && wg[i] == wg[i] && wg[i] != 0) {
            wsum += wg[i]; ++n;
          }
        }
      }
      out[0] = (double)n;
      out[1] = wsum;
    } else {
      if(real) {
        for(int i = 0; i != l; ++i) if(x[i] == x[i]) ++n; // This loop is faster
      } else {
        for(int i = 0; i != l; ++i) if(x[i] != Vector<RTYPE>::get_na()) ++n;
      }
      out[0] = (double)n;
    }
    if(setn) {
      Rf_namesgets(out, get_stats_names(d));
      Rf_classgets(out, CharacterVector::create("qsu","table"));
    }
    return out;
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    NumericMatrix out = no_init_matrix(ng, d);
    std::fill_n(out.begin(), ng*(1+weights), 0.0); // works ?? -> yes
    std::fill(out.begin()+ng*(1+weights), out.end(), NA_REAL);
    NumericMatrix::Column n = out(_, 0);
    if(weights) {
      NumericVector wg = w;
      NumericMatrix::Column wsum = out(_, 1);
      if(real) {
        for(int i = 0; i != l; ++i) {
          if(x[i] == x[i] && wg[i] == wg[i] && wg[i] != 0) {
            wsum[g[i]-1] += wg[i]; ++n[g[i]-1];
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(x[i] != Vector<RTYPE>::get_na() && wg[i] == wg[i] && wg[i] != 0) {
            wsum[g[i]-1] += wg[i]; ++n[g[i]-1];
          }
        }
      }
    } else {
      if(real) {
        for(int i = 0; i != l; ++i) if(x[i] == x[i]) ++n[g[i]-1];
      } else {
        for(int i = 0; i != l; ++i) if(x[i] != Vector<RTYPE>::get_na()) ++n[g[i]-1];
      }
    }
    if(setn) {
      Rf_dimnamesgets(out, List::create(gn, get_stats_names(d)));
      Rf_classgets(out, CharacterVector::create("qsu","matrix","table"));
    }
    return out;
  }
}

template <int RTYPE>
NumericMatrix fnobs5pImpl(Vector<RTYPE> x, bool ext = false, int ng = 0, IntegerVector g = 0, int npg = 0, IntegerVector pg = 0, SEXP w = R_NilValue, bool real = false, bool array = true, SEXP gn = R_NilValue) {

  bool weights = !Rf_isNull(w);
  int l = x.size(), d = ((ext) ? 7 : 5) + weights;
  if(pg.size() != l) stop("length(pid) must match nrow(X)");

  if(ng == 0) {
    int n = 0, npgc = 0;
    // bool npgs[npg+1];
    // memset(npgs, true, sizeof(bool)*(npg+1));
    std::vector<bool> npgs(npg+1, true);
    double wsum = 0.0;
    if(weights) {
      NumericVector wg = w;
      if(real) {
        for(int i = 0; i != l; ++i) {
          if(x[i] == x[i] && wg[i] == wg[i] && wg[i] != 0) {
            wsum += wg[i]; ++n;
          }
          if(npgs[pg[i]-1]) {
            ++npgc;
            npgs[pg[i]-1] = false;
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(x[i] != Vector<RTYPE>::get_na() && wg[i] == wg[i] && wg[i] != 0) {
            wsum += wg[i]; ++n;
          }
          if(npgs[pg[i]-1]) {
            ++npgc;
            npgs[pg[i]-1] = false;
          }
        }
      }
    } else {
      if(real) {
        for(int i = 0; i != l; ++i) {
          if(x[i] == x[i]) ++n;
          if(npgs[pg[i]-1]) {
            ++npgc;
            npgs[pg[i]-1] = false;
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(x[i] != Vector<RTYPE>::get_na()) ++n;
          if(npgs[pg[i]-1]) {
            ++npgc;
            npgs[pg[i]-1] = false;
          }
        }
      }
    }
    NumericMatrix out = no_init_matrix(3, d);
    out[0] = (double)n;
    out[1] = (double)npgc;
    out[2] = out[0]/out[1];
    if(weights) {
      out[3] = (double)wsum;
      out[4] = (double)npgc;
      out[5] = out[3]/out[4];
    }
    std::fill(out.begin()+3*(1+weights), out.end(), NA_REAL);
    if(!array) {
      Rf_dimnamesgets(out, List::create(CharacterVector::create("Overall","Between","Within"), get_stats_names(d, true)));
      Rf_classgets(out, CharacterVector::create("qsu","matrix","table"));
    }
    return out;
  } else { // with groups
    if(g.size() != l) stop("length(g) must match nrow(X)");
    NumericMatrix out = no_init_matrix(ng*d, 3);
    std::fill_n(out.begin(), ng*(1+weights), 0.0); // works ? -> yes
    std::fill(out.begin()+ng*(1+weights), out.end(), NA_REAL);
    NumericMatrix::Column n = out(_, 0);
    NumericMatrix::Column gnpids = out(_, 1);
    std::fill_n(gnpids.begin(), ng, 0.0);
    // bool groupids[ng][npg]; // could do +1 trick, but that could be costly in term of memory, if few g and many pg
    // memset(groupids, true, sizeof(bool)*ng*npg);
    LogicalMatrix groupids = no_init_matrix(ng, npg);
    std::fill(groupids.begin(), groupids.end(), true);
    if(weights) {
      NumericVector wg = w;
      if(real) {
        for(int i = 0; i != l; ++i) {
          if(x[i] == x[i] && wg[i] == wg[i] && wg[i] != 0) {
            n[g[i]+ng-1] += wg[i]; ++n[g[i]-1];
            if(groupids(g[i]-1, pg[i]-1)) {
              ++gnpids[g[i]-1];
              groupids(g[i]-1, pg[i]-1) = false;
            }
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(x[i] != Vector<RTYPE>::get_na() && wg[i] == wg[i] && wg[i] != 0) {
            n[g[i]+ng-1] += wg[i]; ++n[g[i]-1];
            if(groupids(g[i]-1, pg[i]-1)) {
              ++gnpids[g[i]-1];
              groupids(g[i]-1, pg[i]-1) = false;
            }
          }
        }
      }
    } else {
      if(real) {
        for(int i = 0; i != l; ++i) {
          if(x[i] == x[i]) {
            ++n[g[i]-1];
            if(groupids(g[i]-1, pg[i]-1)) {
              ++gnpids[g[i]-1];
              groupids(g[i]-1, pg[i]-1) = false;
            }
          }
        }
      } else {
        for(int i = 0; i != l; ++i) {
          if(x[i] != Vector<RTYPE>::get_na()) {
            ++n[g[i]-1];
            if(groupids(g[i]-1, pg[i]-1)) {
              ++gnpids[g[i]-1];
              groupids(g[i]-1, pg[i]-1) = false;
            }
          }
        }
      }
    }
    NumericMatrix::Column nt = out(_, 2);
    if(weights) {
      for(int i = 0; i != ng; ++i) {
        gnpids[ng+i] = gnpids[i];
        nt[i] = n[i] / gnpids[i];
        nt[ng+i] = n[ng+i] / gnpids[i];
      }
    } else {
      for(int i = 0; i != ng; ++i) nt[i] = n[i] / gnpids[i];
    }
    if(!array) {
      Rf_dimgets(out, Dimension(ng, d, 3));
      Rf_dimnamesgets(out, List::create(gn, get_stats_names(d, true), CharacterVector::create("Overall","Between","Within")));
      Rf_classgets(out, CharacterVector::create("qsu","array","table"));
    }
    return out;
  }
}


// [[Rcpp::export]]
SEXP fbstatslCpp(const List& x, bool ext = false, int ng = 0, const IntegerVector& g = 0,
                 int npg = 0, const IntegerVector& pg = 0,
                 const SEXP& w = R_NilValue, bool stable_algo = true, bool array = true,
                 const SEXP& gn = R_NilValue) {
  bool weights = !Rf_isNull(w);
  int col = x.size(), d = ((ext) ? 7 : 5) + weights;

  if(npg == 0) { // No panel
    if(ng == 0) { // No groups
      NumericMatrix out = no_init_matrix(col, d);
      for(int j = col; j--; ) {
       switch(TYPEOF(x[j])) {
       case REALSXP:{
         NumericVector column = x[j];
         if(Rf_isObject(column)) out(j, _) = fnobs5Impl<REALSXP>(column, ext, 0, 0, w, true);
         else out(j, _) = fbstatstemp(column, ext, 0, 0, w, false, stable_algo);
         break;
       }
       case INTSXP: {
         IntegerVector column = x[j];
         if(Rf_isObject(column)) out(j, _) = fnobs5Impl<INTSXP>(column, ext, 0, 0, w);
         else out(j, _) = fbstatstemp(x[j], ext, 0, 0, w, false, stable_algo);
         break;
       }
       case STRSXP: out(j, _) = fnobs5Impl<STRSXP>(x[j], ext, 0, 0, w);
         break;
       case LGLSXP: out(j, _) = fnobs5Impl<LGLSXP>(x[j], ext, 0, 0, w);
         break;
       default: stop("Not supported SEXP type!");
       }
      }
      Rf_dimnamesgets(out, List::create(Rf_getAttrib(x, R_NamesSymbol), get_stats_names(d)));
      Rf_classgets(out, CharacterVector::create("qsu","matrix","table"));
      return out;
    } else {
      if(array) {
        NumericMatrix out = no_init_matrix(d*ng, col);
        for(int j = col; j--; ) {
          switch(TYPEOF(x[j])) {
          case REALSXP:{
            NumericVector column = x[j];
            if(Rf_isObject(column)) out(_, j) = fnobs5Impl<REALSXP>(column, ext, ng, g, w, true);
            else out(_, j) = fbstatstemp(column, ext, ng, g, w, false, stable_algo);
            break;
          }
          case INTSXP: {
            IntegerVector column = x[j];
            if(Rf_isObject(column)) out(_, j) = fnobs5Impl<INTSXP>(column, ext, ng, g, w);
            else out(_, j) = fbstatstemp(x[j], ext, ng, g, w, false, stable_algo);
            break;
          }
          case STRSXP: out(_, j) = fnobs5Impl<STRSXP>(x[j], ext, ng, g, w);
            break;
          case LGLSXP: out(_, j) = fnobs5Impl<LGLSXP>(x[j], ext, ng, g, w);
            break;
          default: stop("Not supported SEXP type!");
          }
        }
        Rf_dimgets(out, Dimension(ng, d, col));
        Rf_dimnamesgets(out, List::create(gn, get_stats_names(d), Rf_getAttrib(x, R_NamesSymbol)));
        Rf_classgets(out, CharacterVector::create("qsu","array","table"));
        return out;
      } else {
        List out(col);
        for(int j = col; j--; ) {
          switch(TYPEOF(x[j])) {
          case REALSXP:{
            NumericVector column = x[j];
            if(Rf_isObject(column)) out[j] = fnobs5Impl<REALSXP>(column, ext, ng, g, w, true, true, gn);
            else out[j] = fbstatstemp(column, ext, ng, g, w, true, stable_algo, gn);
            break;
          }
          case INTSXP: {
            IntegerVector column = x[j];
            if(Rf_isObject(column)) out[j] = fnobs5Impl<INTSXP>(column, ext, ng, g, w, false, true, gn);
            else out[j] = fbstatstemp(x[j], ext, ng, g, w, true, stable_algo, gn);
            break;
          }
          case STRSXP: out[j] = fnobs5Impl<STRSXP>(x[j], ext, ng, g, w, false, true, gn);
            break;
          case LGLSXP: out[j] = fnobs5Impl<LGLSXP>(x[j], ext, ng, g, w, false, true, gn);
            break;
          default: stop("Not supported SEXP type!");
          }
        }
        Rf_setAttrib(out, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
        return out;
      }
    }
  } else { // with panel
    if(ng == 0) {
      if(array) {
        NumericMatrix out = no_init_matrix(d*3, col);
        for(int j = col; j--; ) {
          switch(TYPEOF(x[j])) {
          case REALSXP:{
            NumericVector column = x[j];
            if(Rf_isObject(column)) out(_, j) = fnobs5pImpl<REALSXP>(column, ext, 0, 0, npg, pg, w, true);
            else out(_, j) = as<NumericVector>(fbstatsCpp(column, ext, 0, 0, npg, pg, w, stable_algo, true, false)); // or Rf_coerce ?
            break;
          }
          case INTSXP: {
            IntegerVector column = x[j];
            if(Rf_isObject(column)) out(_, j) = fnobs5pImpl<INTSXP>(column, ext, 0, 0, npg, pg, w);
            else out(_, j) = as<NumericVector>(fbstatsCpp(x[j], ext, 0, 0, npg, pg, w, stable_algo, true, false));
            break;
          }
          case STRSXP: out(_, j) = fnobs5pImpl<STRSXP>(x[j], ext, 0, 0, npg, pg, w);
            break;
          case LGLSXP: out(_, j) = fnobs5pImpl<LGLSXP>(x[j], ext, 0, 0, npg, pg, w);
            break;
          default: stop("Not supported SEXP type!");
          }
        }
        Rf_dimgets(out, Dimension(3, d, col));
        Rf_dimnamesgets(out, List::create(CharacterVector::create("Overall","Between","Within"),
                        get_stats_names(d, true), Rf_getAttrib(x, R_NamesSymbol)));
        Rf_classgets(out, CharacterVector::create("qsu","array","table"));
        return out;
      } else {
        List out(col);
        for(int j = col; j--; ) {
          switch(TYPEOF(x[j])) {
          case REALSXP:{
            NumericVector column = x[j];
            if(Rf_isObject(column)) out[j] = fnobs5pImpl<REALSXP>(column, ext, 0, 0, npg, pg, w, true, false, gn);
            else out[j] = fbstatsCpp(column, ext, 0, 0, npg, pg, w, stable_algo, false, true, gn);
            break;
          }
          case INTSXP: {
            IntegerVector column = x[j];
            if(Rf_isObject(column)) out[j] = fnobs5pImpl<INTSXP>(column, ext, 0, 0, npg, pg, w, false, false, gn);
            else out[j] = fbstatsCpp(x[j], ext, 0, 0, npg, pg, w, stable_algo, false, true, gn);
            break;
          }
          case STRSXP: out[j] = fnobs5pImpl<STRSXP>(x[j], ext, 0, 0, npg, pg, w, false, false, gn);
            break;
          case LGLSXP: out[j] = fnobs5pImpl<LGLSXP>(x[j], ext, 0, 0, npg, pg, w, false, false, gn);
            break;
          default: stop("Not supported SEXP type!");
          }
        }
        Rf_setAttrib(out, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
        return out;
      }
    } else {
      if(array) {
        NumericMatrix out = no_init_matrix(d*3*ng, col);
        for(int j = col; j--; ) {
          switch(TYPEOF(x[j])) {
          case REALSXP:{
            NumericVector column = x[j];
            if(Rf_isObject(column)) out(_, j) = fnobs5pImpl<REALSXP>(column, ext, ng, g, npg, pg, w, true);
            else out(_, j) = as<NumericVector>(fbstatsCpp(column, ext, ng, g, npg, pg, w, stable_algo, true, false)); // or Rf_coerce ?
            break;
          }
          case INTSXP: {
            IntegerVector column = x[j];
            if(Rf_isObject(column)) out(_, j) = fnobs5pImpl<INTSXP>(column, ext, ng, g, npg, pg, w);
            else out(_, j) = as<NumericVector>(fbstatsCpp(x[j], ext, ng, g, npg, pg, w, stable_algo, true, false));
            break;
          }
          case STRSXP: out(_, j) = fnobs5pImpl<STRSXP>(x[j], ext, ng, g, npg, pg, w);
            break;
          case LGLSXP: out(_, j) = fnobs5pImpl<LGLSXP>(x[j], ext, ng, g, npg, pg, w);
            break;
          default: stop("Not supported SEXP type!");
          }
        }
        Rf_dimgets(out, IntegerVector::create(ng, d, 3, col));
        Rf_dimnamesgets(out, List::create(gn, get_stats_names(d, true),
                 CharacterVector::create("Overall","Between","Within"), Rf_getAttrib(x, R_NamesSymbol)));
        Rf_classgets(out, CharacterVector::create("qsu","array","table"));
        return out;
      } else {
        List out(col);
        for(int j = col; j--; ) {
          switch(TYPEOF(x[j])) {
          case REALSXP:{
            NumericVector column = x[j];
            if(Rf_isObject(column)) out[j] = fnobs5pImpl<REALSXP>(column, ext, ng, g, npg, pg, w, true, false, gn);
            else out[j] = fbstatsCpp(column, ext, ng, g, npg, pg, w, stable_algo, false, true, gn);
            break;
          }
          case INTSXP: {
            IntegerVector column = x[j];
            if(Rf_isObject(column)) out[j] = fnobs5pImpl<INTSXP>(column, ext, ng, g, npg, pg, w, false, false, gn);
            else out[j] = fbstatsCpp(x[j], ext, ng, g, npg, pg, w, stable_algo, false, true, gn);
            break;
          }
          case STRSXP: out[j] = fnobs5pImpl<STRSXP>(x[j], ext, ng, g, npg, pg, w, false, false, gn);
            break;
          case LGLSXP: out[j] = fnobs5pImpl<LGLSXP>(x[j], ext, ng, g, npg, pg, w, false, false, gn);
            break;
          default: stop("Not supported SEXP type!");
          }
        }
        Rf_setAttrib(out, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
        return out;
      }
    }
  }
}


// Old / Experimental:
//
// template <>
// NumericVector fnobs5Impl(Vector<CPLXSXP> x, int ng, IntegerVector g, bool real, bool setn) {
//   stop("Not supported SEXP type!");
// }
//
// template <>
// NumericVector fnobs5Impl(Vector<VECSXP> x, int ng, IntegerVector g, bool real, bool setn) {
//   stop("Not supported SEXP type!");
// }
//
// template <>
// NumericVector fnobs5Impl(Vector<RAWSXP> x, int ng, IntegerVector g, bool real, bool setn) {
//   stop("Not supported SEXP type!");
// }
//
// template <>
// NumericVector fnobs5Impl(Vector<EXPRSXP> x, int ng, IntegerVector g, bool real, bool setn) {
//   stop("Not supported SEXP type!");
// }
//
// // [[Rcpp::export]]
// NumericVector fnobs5Cpp(SEXP x, int ng = 0, IntegerVector g = 0, bool real = false, bool setn = true){
//   RCPP_RETURN_VECTOR(fnobs5Impl, x, ng, g, real, setn);
// }




// // [[Rcpp::export]]
// SEXP fbstatsCpp(NumericVector x, int ng = 0, IntegerVector g = 0, IntegerVector gs = 0,
//                 int npg = 0, IntegerVector pg = 0, IntegerVector pgs = 0, // SEXP w,
//                 bool narm = true) {
//   int l = x.size();
//  if(ng == 0 && npg == 0) { // No groups, no panel !!
//     int n = 0;
//     double min = 0, max = 0, sum = 0, sq_sum = 0;
//     if(narm) {
//       int j = l-1;
//       while(std::isnan(x[j]) && j!=0) --j;
//       min = x[j]; max = x[j]; sum = x[j]; sq_sum = x[j];
//       if(j != 0) for(int i = j; i--; ) {
//           if(std::isnan(x[i])) continue;
//           sum += x[i];
//           sq_sum += x[i] * x[i];
//           if(min>x[i]) min = x[i];
//           if(max<x[i]) max = x[i];
//           ++n;
//         }
//       l = n;
//     } else {
//       for(int i = l; i--; ) {
//         sum += x[i];
//         sq_sum += x[i] * x[i];
//         if(min>x[i]) min = x[i];
//         if(max<x[i]) max = x[i];
//       }
//     }
//     NumericVector result(5);
//     result[0] = l;
//     result[1] = sum / l;
//     double ssr = (sq_sum - (result[1]*result[1])*l);
//     result[2] = sqrt(ssr/(l-1));
//     result[3] = min;
//     result[4] = max;
//     result.attr("names") = CharacterVector::create("N","Mean","SD","Min","Max");
//     return result;
//  } else if (npg == 0) {
//    int k = 0;
//    NumericMatrix result = no_init_matrix(ng, 5); // Only if combine.by !!!!
//    std::fill(result.begin(), result.end(), NA_REAL);
//    NumericMatrix::Column n = result( _ , 0);
//    NumericMatrix::Column sum = result( _ , 1);
//    NumericMatrix::Column sq_sum = result( _ , 2);
//    NumericMatrix::Column min = result( _ , 3);
//    NumericMatrix::Column max = result( _ , 4);
//    for(int i = l; i--; ) { // Best ??
//      k = g[i]-1;
//      if(!std::isnan(x[i])) { // faster way to code this ??? -> Not Bad at all
//        if(std::isnan(sum[k])) {
//          sum[k] = x[i];
//          sq_sum[k] = x[i]*x[i];
//          min[k] = x[i];
//          max[k] = x[i];
//          n[k] = 1;
//        } else { // integer for subsetting ??
//          sum[k] += x[i];
//          sq_sum[k] += x[i]*x[i];
//          if(min[k] > x[i]) min[k] = x[i];
//          if(max[k] < x[i]) max[k] = x[i];
//          ++n[k];
//        }
//      }
//    }
//    sum = sum / n;
//    sq_sum = sqrt((sq_sum - (sum*sum)*n)/(n-1));
//    return result;
//  } else if (ng == 0) {
//    // ....
//  }
//  return R_NilValue;
// }
//
// // [[Rcpp::export]]
// SEXP test(NumericVector x) {
//   int l = x.size();
//   int j = l-1;
//     while(std::isnan(x[j]) && j!=0) --j; // right -- before ??
//     return NumericVector::create(j);
// }
//


// #include <Rcpp.h>
// #include <unordered_set>
// using namespace Rcpp;
//
// // [[Rcpp::export]]
// NumericVector fbstats(NumericVector x, bool narm = false) { // possibly try quick conversion to factor??
//   int l = x.size();
//   //NumericVector un = unique(x); // fastest for now. see how constructed..
//   //std::sort(x.begin(), x.end());
//   //std::unordered_set<double> newvalue;
//   //std::unordered_map<double, bool> counts; // Also too slow!!
//   // https://stackoverflow.com/questions/23150905/effective-unique-on-unordered-elements
//   //std::vector<bool> set(1000000000); // simple: just put true if already occurred -> Needs to be positive integers!!
//   //int un = 0;
//   //NumericVector y = x * 100000;
//   double min = x[0]; // what about NA_RM of the first element in NA??
//   double max = x[0];
//   double sum = 0;
//   double sq_sum = 0;
//   //double c_sum = 0;
//   //double f_sum = 0;
//   if(narm) {
//     int n = 0;
//     for(int i = l; i--; ) {
//       if(ISNAN(x[i])) continue;
//       sum += x[i];
//       sq_sum += x[i] * x[i];
//       //c_sum += sq_sum * x[i];
//       //f_sum += c_sum * x[i];
//       if(min>x[i]) min = x[i];
//       if(max<x[i]) max = x[i];
//       //counts[x[i]]++;
//       //newvalue.insert(x[i]);
//       //if( ! set[y[i]]) {
//       //  set[y[i]] = true;
//       //  ++un;
//       //}
//       //if( ! counts[x[i]]) {
//       //  counts[x[i]] = true;
//       //  ++un;
//       //}
//       n++;
//     }
//     l = n;
//   } else {
//     for(int i = l; i--; ) {
//       sum += x[i];
//       sq_sum += x[i] * x[i];
//       //c_sum += sq_sum * x[i];
//       //f_sum += c_sum * x[i];
//       if(min>x[i]) min = x[i];
//       if(max<x[i]) max = x[i];
//       //counts[x[i]]++;
//       //newvalue.insert(x[i]);
//       //if( ! set[y[i]]) {
//       //  set[y[i]] = true;
//       //  ++un;
//       //}
//       //if( ! counts[x[i]]) {
//       //  counts[x[i]] = true;
//       //  ++un;
//       //}
//     }
//   }
//   //un = newvalue.size();
//   NumericVector result(5);
//   result[0] = l;
//   //result[1] = un;
//   result[1] = sum / l; // sqrt(sq_sum / l - result[1] * result[1]) The one below works!!
//   double ssr = (sq_sum - (result[1]*result[1])*l);
//   result[2] = sqrt(ssr/(l-1)); // https://sureshemre.wordpress.com/2012/04/21/how-to-compute-standard-deviation-in-one-pass/
//   result[3] = min;
//   result[4] = max;
//   //result[5] = ((c_sum - 3*result[1]*sq_sum + 3*pow(result[1],2)*sum - pow(result[1],3)*l)/l)/pow(ssr/(l-1),3/2); // These are not correct
//   //result[6] = l*(f_sum - (result[1]*result[1]*result[1]*result[1])*l)/(ssr*ssr);
//   return result;
// }
