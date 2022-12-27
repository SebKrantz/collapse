context("fquantile")

probs1 <- c(0, 0.25, 0.5, 0.75, 1)
probs2 <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

for(x in mtcars) {
  for(o in list(NULL, radixorder(x))) {
    for(Qprobs in list(probs1, probs2)) {
      for(t in 5:9) {
        expect_true(all_obj_equal(
                    fquantile(x, Qprobs, type = t, o = o),
                    fquantile(x, Qprobs, type = t, o = o, na.rm = FALSE),
                    quantile(x, Qprobs, type = t)))
        for(j in 1:3) {
          w = rep(j + rnorm(1, sd = 0.05), 32)
          expect_true(all_obj_equal(
                       fquantile(x, Qprobs, type = t),
                       fquantile(x, Qprobs, type = t, w = w, o = o, na.rm = FALSE),
                       fquantile(x, Qprobs, type = t, w = w, o = o)))
        }
      }
    }
  }
}

expect_equal(.quantile(1:2), c(1.00, 1.25, 1.50, 1.75, 2.00))
expect_equal(.quantile(1:3), c(1.0, 1.5, 2.0, 2.5, 3.0))
expect_equal(.quantile(1:2, na.rm = FALSE), c(1.00, 1.25, 1.50, 1.75, 2.00))
expect_equal(.quantile(1:3, na.rm = FALSE), c(1.0, 1.5, 2.0, 2.5, 3.0))

for(na_rm in c(TRUE, FALSE)) {
  for(t in 5:9) {
    expect_equal(.quantile(0, type = t, na.rm = na_rm), c(0,0,0,0,0))
    expect_equal(.quantile(c(0, 0), type = t, na.rm = na_rm), c(0,0,0,0,0))
    expect_equal(.quantile(c(0, 0, 0), type = t, na.rm = na_rm), c(0,0,0,0,0))
    expect_equal(.quantile(0L, type = t, na.rm = na_rm), rep.int(0L, 5))
    expect_equal(.quantile(c(0L, 0L), type = t, na.rm = na_rm), rep.int(0L, 5))
    expect_equal(.quantile(c(0L, 0L, 0L), type = t, na.rm = na_rm), rep.int(0L, 5))
    expect_equal(.quantile(numeric(0), type = t, na.rm = na_rm), rep(NA_real_, 5))
    expect_equal(.quantile(integer(0), type = t, na.rm = na_rm), rep(NA_real_, 5))
  }
}

for(x in na_insert(airquality, 0.05)) {
  for(o in list(NULL, radixorder(x))) {
    for(Qprobs in list(probs1, probs2)) {
      for(t in 5:9) {
        expect_equal(fquantile(x, Qprobs, type = t, o = o),
                      quantile(x, Qprobs, type = t, na.rm = TRUE))
        for(j in 1:3) {
          w = rep(j + rnorm(1, sd = 0.05), 153)
          expect_equal(fquantile(x, Qprobs, type = t),
                       fquantile(x, Qprobs, type = t, w = w, o = o))
        }
      }
    }
  }
}


