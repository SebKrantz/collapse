/*
 This code is adapted from the kit package: https://github.com/2005m/kit
 and licensed under a GPL-3.0 license.
*/

#include <R.h>
#include <Rinternals.h>
#include <stdint.h> // needed for uintptr_t on linux

#define NOGE(x, l) ((x < 0 && x != NA_INTEGER) || (x >= l))
#define HASH(key, K)  (3141592653U * (unsigned int)(key) >> (32 - (K)))
#define HASHK(key, K)  (3141592653U * (unsigned int)(key) >> (K))
#define N_ISNAN(x, y) (!ISNAN(x) && !ISNAN(y))
#define B_IsNA(x, y)  (R_IsNA(x) && R_IsNA(y))
#define B_IsNaN(x, y) (R_IsNaN(x) && R_IsNaN(y))
#define B_ISNAN(x, y) (ISNAN(x) && ISNAN(y))
#define C_IsNA(x)     (R_IsNA(x.r) || R_IsNA(x.i))
#define C_IsNaN(x)    (R_IsNaN(x.r) || R_IsNaN(x.i))
#define C_ISNAN(x, y) (B_ISNAN(x, y) || (N_ISNAN(x, y) && x == y))
#define REQUAL(x, y)  (N_ISNAN(x, y) ? (x == y) : (B_IsNA(x, y) || B_IsNaN(x, y)))
#define CEQUAL(x, y) ((N_ISNAN(x.r, x.i) && N_ISNAN(y.r, y.i)) ? (x.r == y.r && x.i == y.i) : (C_IsNA(x) ? C_IsNA(y) : (C_IsNA(y) ? 0 : (C_ISNAN(x.r, y.r) && C_ISNAN(x.i, y.i)))))

union uno { double d; unsigned int u[2]; };

