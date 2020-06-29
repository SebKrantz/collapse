/*  R : A Computer Language for Statistical Data Analysis
*
  *  Copyright (C) 1999-2016	The R Core Team
*
  *  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
  *  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
  *  You should have received a copy of the GNU General Public License
*  along with this program; if not, a copy is available at
*  https://www.R-project.org/Licenses/.
*/

// #ifdef HAVE_CONFIG_H
// # include <config.h>
// #endif

// #include "data.table.h"
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

// #include <R.h>
// #include "ts.h"


/* cor is the autocorrelations starting from 0 lag*/
  static void uni_pacf(double *cor, double *p, int nlag)
{
  double a, b, c, *v, *w;

  v = (double*) R_alloc(nlag, sizeof(double));
  w = (double*) R_alloc(nlag, sizeof(double));
  w[0] = p[0] = cor[1];
  for(int ll = 1; ll < nlag; ll++) {
    a = cor[ll+1];
    b = 1.0;
    for(int i = 0; i < ll; i++) {
      a -= w[i] * cor[ll - i];
      b -= w[i] * cor[i + 1];
    }
    p[ll] = c = a/b;
    if(ll+1 == nlag) break;
    w[ll] = c;
    for(int i = 0; i < ll; i++)
      v[ll-i-1] = w[i];
    for(int i = 0; i < ll; i++)
      w[i] -= c*v[i];
  }
  }

SEXP pacf1(SEXP acf, SEXP lmax)
{
  int lagmax = asInteger(lmax);
  acf = PROTECT(coerceVector(acf, REALSXP));
  SEXP ans = PROTECT(allocVector(REALSXP, lagmax));
  uni_pacf(REAL(acf), REAL(ans), lagmax);
  SEXP d = PROTECT(allocVector(INTSXP, 3));
  INTEGER(d)[0] = lagmax;
  INTEGER(d)[1] = INTEGER(d)[2] = 1;
  setAttrib(ans, R_DimSymbol, d);
  UNPROTECT(3);
  return ans;
}
