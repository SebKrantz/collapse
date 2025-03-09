context("fquantile, and quantiles with fnth")

test_zero_weights <- FALSE

probs1 <- c(0, 0.25, 0.5, 0.75, 1)
probs2 <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)

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
                       .quantile(x, Qprobs, type = t),
                       .quantile(x, Qprobs, type = t, w = w, o = o, na.rm = FALSE),
                       .quantile(x, Qprobs, type = t, w = w, o = o)))
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

    expect_equal(.quantile(0, w = 1, type = t, na.rm = na_rm), c(0,0,0,0,0))
    expect_equal(.quantile(c(0, 0), w = c(1, 1), type = t, na.rm = na_rm), c(0,0,0,0,0))
    expect_equal(.quantile(c(0, 0, 0), w = c(1, 1, 1), type = t, na.rm = na_rm), c(0,0,0,0,0))
    expect_equal(.quantile(0L, w = 1, type = t, na.rm = na_rm), rep.int(0L, 5))
    expect_equal(.quantile(c(0L, 0L), w = c(1, 1), type = t, na.rm = na_rm), rep.int(0L, 5))
    expect_equal(.quantile(c(0L, 0L, 0L), w = c(1, 1, 1), type = t, na.rm = na_rm), rep.int(0L, 5))

    expect_equal(.quantile(numeric(0), type = t, na.rm = na_rm), rep(NA_real_, 5))
    expect_equal(.quantile(integer(0), type = t, na.rm = na_rm), rep(NA_real_, 5))
    expect_equal(.quantile(NA_real_, type = t, na.rm = na_rm), rep(NA_real_, 5))
    expect_equal(.quantile(NA_integer_, type = t, na.rm = na_rm), rep(NA_real_, 5))
    expect_equal(.quantile(NA_real_, w = NA_real_, type = t, na.rm = na_rm), rep(NA_real_, 5))
    expect_equal(.quantile(NA_integer_, w = NA_real_, type = t, na.rm = na_rm), rep(NA_real_, 5))
    expect_equal(.quantile(1, w = 0, type = t, na.rm = na_rm), rep(NA_real_, 5))
    expect_equal(.quantile(1L, w = 0, type = t, na.rm = na_rm), rep(NA_real_, 5))
    expect_error(.quantile(1, w = NA_real_, type = t, na.rm = na_rm))
    expect_error(.quantile(1L, w = NA_real_, type = t, na.rm = na_rm))

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
          expect_equal(.quantile(x, Qprobs, type = t),
                       .quantile(x, Qprobs, type = t, w = w, o = o))
        }
      }
    }
  }
}

if(test_zero_weights) {

# Testing behavior with zero weights
for(x in mtcars) {
  for(o in list(NULL, radixorder(x))) {
    for(Qprobs in list(probs1, probs2)) {
      for(t in 5:9) {
        w = na_insert(abs(rnorm(32)), value = 0)
        wn0 = w[w > 0]
        xn0 = x[w > 0]
        on0 = if(length(o)) radixorder(xn0) else NULL
        expect_true(all_obj_equal(
          .quantile(x, Qprobs, type = t, w = w, o = o),
          .quantile(x, Qprobs, type = t, w = w, o = o, na.rm = FALSE),
          .quantile(xn0, Qprobs, type = t, w = wn0, o = on0),
          .quantile(xn0, Qprobs, type = t, w = wn0, o = on0, na.rm = FALSE)
        ))
      }
    }
  }
}

# Zero weights and NA's
for(x in na_insert(mtcars)) {
  for(o in list(NULL, radixorder(x))) {
    for(Qprobs in list(probs1, probs2)) {
      for(t in 5:9) {
        w = na_insert(abs(rnorm(32)), value = 0)
        wn0 = w[w > 0]
        xn0 = x[w > 0]
        on0 = if(length(o)) radixorder(xn0) else NULL
        expect_equal(.quantile(x, Qprobs, type = t, w = w, o = o),
                     .quantile(xn0, Qprobs, type = t, w = wn0, o = on0))
      }
    }
  }
}

}


# Testing with fnth:
.nthquantile <- function(x, probs = c(0.25, 0.5, 0.75), w = NULL, o = NULL, na.rm = TRUE,
                         type = 7L, check.o = is.null(attr(o, "sorted")), ...) {
    vapply(probs, fnth.default, 1.0, x = x, w = w, ties = type,
           o = o, na.rm = na.rm, check.o = check.o, USE.NAMES = FALSE,
           use.g.names = FALSE, ...)
}

probs <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

gmtc = GRP(rep(1L, 32))
gmtcus = gmtc
gmtcus$ordered %-=% 1L

for(g in list(NULL, gmtc, gmtcus)) {
  for(x in mtcars) {
    for(o in list(NULL, radixorder(x))) {
      for(t in 5:9) {
        expect_true(all_obj_equal(
          .quantile(x, probs, type = t, o = o),
          .nthquantile(x, probs, type = t, o = o, g = g),
          .nthquantile(x, probs, type = t, o = o, na.rm = FALSE, g = g)))
        for(j in 1:2) {
          w = rep(j + rnorm(1, sd = 0.05), 32)
          expect_true(all_obj_equal(
            .quantile(x, probs, type = t),
            .nthquantile(x, probs, type = t, w = w, o = o, na.rm = FALSE, g = g),
            .nthquantile(x, probs, type = t, w = w, o = o, g = g)))
        }
      }
    }
  }
}

for(g in list(NULL, rep(1L, 3L))) {
  expect_equal(.nthquantile(1:3, na.rm = FALSE), c(1.5, 2.0, 2.5), g = g)
  expect_equal(.nthquantile(1:3), c(1.5, 2.0, 2.5), g = g)
  for(na_rm in c(TRUE, FALSE)) {
    for(t in 5:9) {
      expect_equal(.nthquantile(c(0, 0, 0), type = t, na.rm = na_rm, g = g), c(0,0,0))
      expect_equal(.nthquantile(c(0L, 0L, 0L), type = t, na.rm = na_rm, g = g), rep.int(0L, 3))
      expect_equal(.nthquantile(c(0, 0, 0), w = c(1, 1, 1), type = t, na.rm = na_rm, g = g), c(0,0,0))
      expect_equal(.nthquantile(c(0L, 0L, 0L), w = c(1, 1, 1), type = t, na.rm = na_rm, g = g), rep.int(0L, 3))
    }
  }
}

for(g in list(NULL, rep(1L, 2L))) {
  expect_equal(.nthquantile(1:2), c(1.25, 1.50, 1.75), g = g)
  expect_equal(.nthquantile(1:2, na.rm = FALSE), c(1.25, 1.50, 1.75), g = g)
  for(na_rm in c(TRUE, FALSE)) {
    for(t in 5:9) {
      expect_equal(.nthquantile(c(0, 0), type = t, na.rm = na_rm, g = g), c(0,0,0))
      expect_equal(.nthquantile(c(0L, 0L), type = t, na.rm = na_rm, g = g), rep.int(0L, 3))
      expect_equal(.nthquantile(c(0, 0), w = c(1, 1), type = t, na.rm = na_rm, g = g), c(0,0,0))
      expect_equal(.nthquantile(c(0L, 0L), w = c(1, 1), type = t, na.rm = na_rm, g = g), rep.int(0L, 3))
    }
  }
}

for(g in list(NULL, 1L)) {
  for(na_rm in c(TRUE, FALSE)) {
    for(t in 5:9) {
      expect_equal(.nthquantile(0, type = t, na.rm = na_rm, g = g), c(0,0,0))
      expect_equal(.nthquantile(0L, type = t, na.rm = na_rm, g = g), rep.int(0L, 3))
      expect_equal(.nthquantile(0, w = 1, type = t, na.rm = na_rm, g = g), c(0,0,0))
      expect_equal(.nthquantile(0L, w = 1, type = t, na.rm = na_rm, g = g), rep.int(0L, 3))
      expect_equal(.nthquantile(NA_real_, type = t, na.rm = na_rm, g = g), rep(NA_real_, 3))
      expect_equal(.nthquantile(NA_integer_, type = t, na.rm = na_rm, g = g), rep(NA_real_, 3))
      expect_equal(.nthquantile(NA_real_, w = NA_real_, type = t, na.rm = na_rm, g = g), rep(NA_real_, 3))
      expect_equal(.nthquantile(NA_integer_, w = NA_real_, type = t, na.rm = na_rm, g = g), rep(NA_real_, 3))
      # expect_equal(.nthquantile(1, w = 0, type = t, na.rm = na_rm, g = g), rep(NA_real_, 3))
      # expect_equal(.nthquantile(1L, w = 0, type = t, na.rm = na_rm, g = g), rep(NA_real_, 3))
      # expect_error(.nthquantile(1, w = NA_real_, type = t, na.rm = na_rm, g = g))
      # expect_error(.nthquantile(1L, w = NA_real_, type = t, na.rm = na_rm, g = g))
    }
  }
}

gaq = GRP(rep(1L, fnrow(airquality)))
gaqus = gaq
gaqus$ordered %-=% 1L

for(g in list(NULL, gaq, gaqus)) {
  for(x in na_insert(airquality, 0.05)) {
    for(o in list(NULL, radixorder(x))) {
      for(t in 5:9) {
        expect_equal(.quantile(x, probs, type = t, o = o),
                     .nthquantile(x, probs, type = t, o = o, g = g))
        for(j in 1:3) {
          w = rep(j + rnorm(1, sd = 0.05), 153)
          expect_equal(.quantile(x, probs, type = t, o = o),
                       .nthquantile(x, probs, type = t, w = w, o = o, g = g))
        }
      }
    }
  }
}

if(test_zero_weights) {

# Testing behavior with zero weights
for(g in list(NULL, gmtc, gmtcus)) {
  for(x in mtcars) {
    for(o in list(NULL, radixorder(x))) {
      for(t in c(1:3, 5:9)) {
        w = fbetween(na_insert(abs(rnorm(32)), 0.15, value = 0), x) # averaging because R's quicksort is not stable
        wn0 = w[w > 0]
        xn0 = x[w > 0]
        on0 = if(length(o)) radixorder(xn0) else NULL
        if(t > 4L) {
          expect_true(all_obj_equal(
            .quantile(x, probs, type = t, w = w, o = o),
            .nthquantile(x, probs, type = t, w = w, o = o, g = g),
            .nthquantile(x, probs, type = t, w = w, o = o, na.rm = FALSE, g = g),
            .nthquantile(xn0, probs, type = t, w = wn0, o = on0),
            .nthquantile(xn0, probs, type = t, w = wn0, o = on0, na.rm = FALSE)
          ))
        } else {
          expect_true(all_obj_equal(
            .nthquantile(x, probs, type = t, w = w, o = o),
            .nthquantile(x, probs, type = t, w = w, o = o, g = g),
            .nthquantile(x, probs, type = t, w = w, o = o, na.rm = FALSE, g = g),
            .nthquantile(xn0, probs, type = t, w = wn0, o = on0),
            .nthquantile(xn0, probs, type = t, w = wn0, o = on0, na.rm = FALSE)
          ))
        }
      }
    }
  }
}

# Zero weights and NA's
for(g in list(NULL, gmtc, gmtcus)) {
  for(x in na_insert(mtcars)) {
    for(o in list(NULL, radixorder(x))) {
      for(t in c(1:3, 5:9)) {
        w = fbetween(na_insert(abs(rnorm(32)), 0.15, value = 0), x) # averaging because R's quicksort is not stable
        wn0 = w[w > 0]
        xn0 = x[w > 0]
        on0 = if(length(o)) radixorder(xn0) else NULL
        if(t > 4L) {
          expect_true(all_obj_equal(
                     .quantile(x, probs, type = t, w = w, o = o),
                     .nthquantile(x, probs, type = t, w = w, o = o, g = g),
                     .nthquantile(xn0, probs, type = t, w = wn0, o = on0)))
        } else {
          expect_true(all_obj_equal(
            .nthquantile(x, probs, type = t, w = w, o = o),
            .nthquantile(x, probs, type = t, w = w, o = o, g = g),
            .nthquantile(xn0, probs, type = t, w = wn0, o = on0)))
        }
      }
    }
  }
}

}
