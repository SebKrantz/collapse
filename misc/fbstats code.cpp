// for(int i = 0; i != l; ++i) within[i] = x[i] - sum[pg[i]-1] + osum; // if-check for NA's is not faster
NumericMatrix result = no_init_matrix(3, d);
// Added this block. It should be possible to simply do fbstatstemp(sum, ext, 0, 0, sumw, false) as below, but for some reason the result is not right for the standard deviation and higher moments..
if(weights) {
  NumericVector between = no_init_vector(l);
  for(int i = 0; i != l; ++i) {
    if(std::isnan(x[i])) {
      between[i] = within[i] = NA_REAL;
    } else {
      between[i] = sum[pg[i]-1];
      within[i] = x[i] - between[i] + osum;
    }
  }
  result(1, _) = fbstatstemp(between, ext, 0, 0, w, false);
  result[1] = npg;
  result[2] /= npg;
} else {
  for(int i = 0; i != l; ++i) within[i] = x[i] - sum[pg[i]-1] + osum;
  result(1, _) = fbstatstemp(sum, ext, 0, 0, w, false);
  result[2] /= result[1];
}
result(0, _) = fbstatstemp(x, ext, 0, 0, w, false);
// result(1, _) = (weights) ? fbstatstemp(sum, ext, 0, 0, sumw, false) : fbstatstemp(sum, ext, 0, 0, w, false);
result(2, _) = fbstatstemp(within, ext, 0, 0, w, false);
// result[2] /= result[1];



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
