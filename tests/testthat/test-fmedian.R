context("fmedian and fnth")



bmean <- base::mean
bsum <- base::sum
bmin <- base::min
bmax <- base::max
bmedian <- stats::median

# rm(list = ls())
set.seed(101)
x <- rnorm(100)
w <- as.integer(round(10*abs(rnorm(100))))  # -> Numeric precision issues in R
wdat <- as.integer(round(10*abs(rnorm(32))))
xNA <- x
wNA <- w
xNA[sample.int(100,20)] <- NA
wNA[is.na(xNA)] <- NA  # only missing weights if x also missing
f <- as.factor(sample.int(10, 100, TRUE))
g <- GRP(mtcars, ~ cyl + vs + am)
gf <- as_factor_GRP(g)
mtcNA <- na_insert(mtcars)
mtcNA[27, 1] <- NA # single group NA !!
m <- as.matrix(mtcars)
mNA <- as.matrix(mtcNA)
mNAc <- mNA
storage.mode(mNAc) <- "character"


nth <- function(x, n, na.rm = FALSE) {
  if(na.rm) {
    if(n > 1) n <- (n-1)/(length(x)-1L)
    x <- na_rm(x)
    if(!length(x)) return(NA_real_)
  } else {
    if(anyNA(x)) return(NA_real_)
  }
  if(n < 1) {
    n <- as.integer((length(x)-1L)*n)+1L
    if(n < 2L) return(bmin(x))
  }
  sort(x, partial = n)[n]
}

wnth <- function(x, n = 0.5, w, na.rm = FALSE, ties = "mean") {
  cc <- complete.cases(x, w)
  if(na.rm) {
    x <- x[cc]
    w <- w[cc]
    if(!length(x)) return(NA_real_)
  } else if(!all(cc)) return(NA_real_)
  sumwh <- bsum(w) * n
  if(sumwh == 0) return(NA_real_)
  if(length(x) < 2L) return(x)
  lp1 <- function(x) if(length(x)) x[length(x)] + 1L else 1L
  mean2 <- function(x) bsum(x) / length(x)
  o <- radixorder(x)
  csumw <- base::cumsum(w[o])
  if(csumw[1L] > sumwh) return(x[o[1L]])
  switch(ties,
         mean = mean2(x[o[lp1(which(csumw < sumwh)):lp1(which(csumw <= sumwh))]]),
         min = x[o[lp1(which(csumw < sumwh))]],
         max = x[o[lp1(which(csumw <= sumwh))]])
}

wmedian <- function(x, w, na.rm = FALSE) wnth(x, 0.5, w, na.rm, "mean")
  # matrixStats::weightedMedian(x, w, ties = ties) -> doesn't always properly average if ties = "mean"...


for (FUN in 1:2) {

  if(FUN == 2L) {
    if(Sys.getenv("OMP") == "TRUE") {
      fmedian <- function(x, ...) collapse::fmedian(x, ..., nthreads = 2L)
    } else break
  }

test_that("fmedian performs like base::median", {
  for(t in c(1L, 5:9)) { # All quantile methods should give the same median value estimate
  expect_equal(fmedian(NA, ties = t), as.double(bmedian(NA)))
  expect_equal(fmedian(NA, na.rm = FALSE, ties = t), as.double(bmedian(NA)))
  expect_equal(fmedian(1, ties = t), bmedian(1, na.rm = TRUE))
  expect_equal(fmedian(1:3, ties = t), bmedian(1:3, na.rm = TRUE))
  expect_equal(fmedian(-1:1, ties = t), bmedian(-1:1, na.rm = TRUE))
  expect_equal(fmedian(1, na.rm = FALSE, ties = t), bmedian(1))
  expect_equal(fmedian(1:3, na.rm = FALSE, ties = t), bmedian(1:3))
  expect_equal(fmedian(-1:1, na.rm = FALSE, ties = t), bmedian(-1:1))
  expect_equal(fmedian(x, ties = t), bmedian(x, na.rm = TRUE))
  expect_equal(fmedian(x, na.rm = FALSE, ties = t), bmedian(x))
  expect_equal(fmedian(xNA, na.rm = FALSE, ties = t), bmedian(xNA))
  expect_equal(fmedian(xNA, ties = t), bmedian(xNA, na.rm = TRUE))
  expect_equal(fmedian(mtcars, ties = t), fmedian(m))
  expect_equal(fmedian(m, ties = t), dapply(m, bmedian, na.rm = TRUE))
  expect_equal(fmedian(m, na.rm = FALSE, ties = t), dapply(m, bmedian))
  expect_equal(fmedian(mNA, na.rm = FALSE, ties = t), dapply(mNA, bmedian))
  expect_equal(fmedian(mNA, ties = t), dapply(mNA, bmedian, na.rm = TRUE))
  expect_equal(fmedian(mtcars, ties = t), dapply(mtcars, bmedian, na.rm = TRUE))
  expect_equal(fmedian(mtcars, na.rm = FALSE, ties = t), dapply(mtcars, bmedian))
  expect_equal(fmedian(mtcNA, na.rm = FALSE, ties = t), dapply(mtcNA, bmedian))
  expect_equal(fmedian(mtcNA, ties = t), dapply(mtcNA, bmedian, na.rm = TRUE))
  expect_equal(fmedian(x, f, ties = t), BY(x, f, bmedian, na.rm = TRUE))
  expect_equal(fmedian(x, f, na.rm = FALSE, ties = t), BY(x, f, bmedian))
  expect_equal(fmedian(xNA, f, na.rm = FALSE, ties = t), BY(xNA, f, bmedian))
  expect_equal(fmedian(xNA, f, ties = t), BY(xNA, f, bmedian, na.rm = TRUE))
  expect_equal(fmedian(m, g, ties = t), BY(m, g, bmedian, na.rm = TRUE))
  expect_equal(fmedian(m, g, na.rm = FALSE, ties = t), BY(m, g, bmedian))
  expect_equal(fmedian(mNA, g, na.rm = FALSE, ties = t), BY(mNA, g, bmedian))
  expect_equal(fmedian(mNA, g, ties = t), BY(mNA, g, bmedian, na.rm = TRUE))
  expect_equal(fmedian(mtcars, g, ties = t), BY(mtcars, g, bmedian, na.rm = TRUE))
  expect_equal(fmedian(mtcars, g, na.rm = FALSE, ties = t), BY(mtcars, g, bmedian))
  expect_equal(fmedian(mtcNA, g, na.rm = FALSE, ties = t), BY(mtcNA, g, bmedian))
  expect_equal(fmedian(mtcNA, g, ties = t), BY(mtcNA, g, bmedian, na.rm = TRUE))
  }
})

test_that("fmedian performs like fmedian with weights all equal", {
  expect_equal(fmedian(NA), fmedian(NA, w = 1))
  expect_equal(fmedian(NA, na.rm = FALSE), fmedian(NA, w = 1, na.rm = FALSE))
  expect_equal(fmedian(1), fmedian(1, w = 3))
  expect_equal(fmedian(1:3), fmedian(1:3, w = rep(1,3)))
  expect_equal(fmedian(-1:1), fmedian(-1:1, w = rep(4.2,3)))
  expect_equal(fmedian(1, na.rm = FALSE), fmedian(1, w = 5, na.rm = FALSE))
  expect_equal(fmedian(1:3, na.rm = FALSE), fmedian(1:3, w = rep(1, 3), na.rm = FALSE))
  expect_equal(fmedian(-1:1, na.rm = FALSE), fmedian(-1:1, w = rep(12, 3), na.rm = FALSE))
  expect_equal(fmedian(x), fmedian(x, w = rep(1,100)))
  expect_equal(fmedian(x, na.rm = FALSE), fmedian(x, w = rep(1, 100), na.rm = FALSE))
  expect_equal(fmedian(xNA, na.rm = FALSE), fmedian(xNA, w = rep(5, 100), na.rm = FALSE))
  expect_equal(fmedian(xNA), fmedian(xNA, w = rep(4, 100)))
  expect_equal(fmedian(m), fmedian(m, w = rep(6587, 32)))
  expect_equal(fmedian(m, na.rm = FALSE), fmedian(m, w = rep(6587, 32), na.rm = FALSE))
  expect_equal(fmedian(mNA, na.rm = FALSE), fmedian(mNA, w = rep(6587, 32), na.rm = FALSE))
  expect_equal(fmedian(mNA), fmedian(mNA, w = rep(6587, 32)))
  expect_equal(fmedian(mtcars), fmedian(mtcars, w = rep(6787, 32)))
  expect_equal(fmedian(mtcars, na.rm = FALSE), fmedian(mtcars, w = rep(6787, 32), na.rm = FALSE))
  expect_equal(fmedian(mtcNA, na.rm = FALSE), fmedian(mtcNA, w = rep(6787, 32), na.rm = FALSE))
  expect_equal(fmedian(mtcNA), fmedian(mtcNA, w = rep(6787, 32)))
  expect_equal(fmedian(x, f), fmedian(x, f, rep(547,100)))
  expect_equal(fmedian(x, f, na.rm = FALSE), fmedian(x, f, rep(6, 100), na.rm = FALSE))
  expect_equal(fmedian(xNA, f, na.rm = FALSE), fmedian(xNA, f, rep(52,100), na.rm = FALSE))
  expect_equal(fmedian(xNA, f), fmedian(xNA, f, rep(5997456,100)))
  expect_equal(fmedian(m, g), fmedian(m, g, rep(546,32)))
  expect_equal(fmedian(m, g, na.rm = FALSE), fmedian(m, g, rep(1,32), na.rm = FALSE))
  expect_equal(fmedian(mNA, g, na.rm = FALSE), fmedian(mNA, g, rep(5,32), na.rm = FALSE))
  expect_equal(fmedian(mNA, g), fmedian(mNA, g, rep(1,32)))
  expect_equal(fmedian(mtcars, g), fmedian(mtcars, g, rep(53,32)))
  expect_equal(fmedian(mtcars, g, na.rm = FALSE), fmedian(mtcars, g, rep(546,32), na.rm = FALSE))
  expect_equal(fmedian(mtcNA, g, na.rm = FALSE), fmedian(mtcNA, g, rep(1,32), na.rm = FALSE))
  expect_equal(fmedian(mtcNA, g), fmedian(mtcNA, g, rep(999,32)))
})

test_that("fmedian with weights performs like wmedian (defined above)", {
  # complete weights
  expect_equal(fmedian(NA, w = 1), wmedian(NA_real_, 1))
  expect_equal(fmedian(NA, w = 1, na.rm = FALSE), wmedian(NA_real_, 1))
  expect_equal(fmedian(1, w = 1), wmedian(1, w = 1))
  expect_equal(fmedian(1:3, w = 1:3), wmedian(1:3, 1:3))
  expect_equal(fmedian(-1:1, w = 1:3), wmedian(-1:1, 1:3))
  expect_equal(fmedian(1, w = 1, na.rm = FALSE), wmedian(1, 1))
  expect_equal(fmedian(1:3, w = c(0.99,3454,1.111), na.rm = FALSE), wmedian(1:3, c(0.99,3454,1.111)))
  expect_equal(fmedian(-1:1, w = 1:3, na.rm = FALSE), wmedian(-1:1, 1:3))
  expect_equal(fmedian(x, w = w), wmedian(x, w))
  expect_equal(fmedian(x, w = w, na.rm = FALSE), wmedian(x, w))
  expect_equal(fmedian(xNA, w = w, na.rm = FALSE), wmedian(xNA, w))
  expect_equal(fmedian(xNA, w = w), wmedian(xNA, w, na.rm = TRUE))
  expect_equal(fmedian(mtcars, w = wdat), fmedian(m, w = wdat))
  expect_equal(fmedian(m, w = wdat), dapply(m, wmedian, wdat, na.rm = TRUE))
  expect_equal(fmedian(m, w = wdat, na.rm = FALSE), dapply(m, wmedian, wdat))
  expect_equal(fmedian(mNA, w = wdat, na.rm = FALSE), dapply(mNA, wmedian, wdat))
  expect_equal(fmedian(mNA, w = wdat), dapply(mNA, wmedian, wdat, na.rm = TRUE))
  expect_equal(fmedian(mtcars, w = wdat), dapply(mtcars, wmedian, wdat, na.rm = TRUE))
  expect_equal(fmedian(mtcars, w = wdat, na.rm = FALSE), dapply(mtcars, wmedian, wdat))
  expect_equal(fmedian(mtcNA, w = wdat, na.rm = FALSE), dapply(mtcNA, wmedian, wdat))
  expect_equal(fmedian(mtcNA, w = wdat), dapply(mtcNA, wmedian, wdat, na.rm = TRUE))
  expect_equal(fmedian(x, f, w), BY(x, f, wmedian, w))
  expect_equal(fmedian(x, f, w, na.rm = FALSE), BY(x, f, wmedian, w))
  expect_equal(fmedian(xNA, f, w, na.rm = FALSE), BY(xNA, f, wmedian, w))
  expect_equal(fmedian(xNA, f, w), BY(xNA, f, wmedian, w, na.rm = TRUE))
  expect_equal(fmedian(m, g, wdat), BY(m, gf, wmedian, wdat))
  expect_equal(fmedian(m, g, wdat, na.rm = FALSE), BY(m, gf, wmedian, wdat))
  expect_equal(fmedian(mNA, g, wdat, na.rm = FALSE),  BY(mNA, gf, wmedian, wdat))
  expect_equal(fmedian(mNA, g, wdat), BY(mNA, gf, wmedian, wdat, na.rm = TRUE))
  expect_equal(fmedian(mtcars, g, wdat), BY(mtcars, gf, wmedian, wdat))
  expect_equal(fmedian(mtcars, g, wdat, na.rm = FALSE), BY(mtcars, gf, wmedian, wdat))
  expect_equal(fmedian(mtcNA, g, wdat, na.rm = FALSE), BY(mtcNA, gf, wmedian, wdat))
  expect_equal(fmedian(mtcNA, g, wdat), BY(mtcNA, gf, wmedian, wdat, na.rm = TRUE))
  # missing weights: Only supported if x is also missing...
  expect_equal(fmedian(NA, w = NA), wmedian(NA_real_, NA_real_))
  expect_equal(fmedian(NA, w = NA, na.rm = FALSE), wmedian(NA_real_, NA_real_))
  expect_equal(fmedian(xNA, w = wNA, na.rm = FALSE), wmedian(xNA, wNA))
  expect_equal(fmedian(xNA, w = wNA), wmedian(xNA, wNA, na.rm = TRUE))
  expect_equal(fmedian(xNA, f, wNA, na.rm = FALSE), BY(xNA, f, wmedian, wNA))
  expect_equal(fmedian(xNA, f, wNA), BY(xNA, f, wmedian, wNA, na.rm = TRUE))
})

test_that("fmedian performs numerically stable", {
  expect_true(all_obj_equal(replicate(50, fmedian(1), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(NA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(NA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(x), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(x, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(m), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(m, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcars), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcars, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(x, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(x, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(m, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(m, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mNA, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcars, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcars, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcNA, g), simplify = FALSE)))
})

test_that("fmedian with complete weights performs numerically stable", {
  expect_true(all_obj_equal(replicate(50, fmedian(1, w = 1), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(NA, w = 1), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(NA, w = 1, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(x, w = w), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(x, w = w, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, w = w, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, w = w), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(m, w = wdat), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(m, w = wdat, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mNA, w = wdat, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mNA, w = wdat), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcars, w = wdat), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcars, w = wdat, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcNA, w = wdat, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcNA, w = wdat), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(x, f, w), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(x, f, w, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, f, w, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, f, w), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(m, g, wdat), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(m, g, wdat, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mNA, g, wdat, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mNA, g, wdat), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcars, g, wdat), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcars, g, wdat, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcNA, g, wdat, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(mtcNA, g, wdat), simplify = FALSE)))
})

test_that("fmedian with missing weights performs numerically stable", {
  expect_true(all_obj_equal(replicate(50, fmedian(NA, w = NA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(NA, w = NA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, w = wNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, w = wNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, f, wNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmedian(xNA, f, wNA), simplify = FALSE)))
})

test_that("fmedian handles special values in the right way", {
  expect_equal(fmedian(NA), NA_real_)
  expect_equal(fmedian(NaN), NaN)
  expect_equal(fmedian(Inf), Inf)
  expect_equal(fmedian(-Inf), -Inf)
  expect_equal(fmedian(TRUE), 1)
  expect_equal(fmedian(FALSE), 0)
  expect_equal(fmedian(NA, na.rm = FALSE), NA_real_)
  expect_equal(fmedian(NaN, na.rm = FALSE), NaN)
  expect_equal(fmedian(Inf, na.rm = FALSE), Inf)
  expect_equal(fmedian(-Inf, na.rm = FALSE), -Inf)
  expect_equal(fmedian(TRUE, na.rm = FALSE), 1)
  expect_equal(fmedian(FALSE, na.rm = FALSE), 0)
  expect_equal(fmedian(c(1,NA)), 1)
  expect_equal(fmedian(c(1,NaN)), 1)
  expect_equal(fmedian(c(1,Inf)), Inf)
  expect_equal(fmedian(c(1,-Inf)), -Inf)
  expect_equal(fmedian(c(FALSE,TRUE)), 0.5)
  expect_equal(fmedian(c(FALSE,FALSE)), 0)
  expect_equal(fmedian(c(1,Inf), na.rm = FALSE), Inf)
  expect_equal(fmedian(c(1,-Inf), na.rm = FALSE), -Inf)
  expect_equal(fmedian(c(FALSE,TRUE), na.rm = FALSE), 0.5)
  expect_equal(fmedian(c(FALSE,FALSE), na.rm = FALSE), 0)
})

test_that("fmedian with weights handles special values in the right way", {
  expect_equal(fmedian(NA, w = 1), NA_real_)
  expect_equal(fmedian(NaN, w = 1), NaN)
  expect_equal(fmedian(Inf, w = 1), Inf)
  expect_equal(fmedian(-Inf, w = 1), -Inf)
  expect_equal(fmedian(TRUE, w = 1), 1)
  expect_equal(fmedian(FALSE, w = 1), 0)
  expect_equal(fmedian(NA, w = 1, na.rm = FALSE), NA_real_)
  expect_equal(fmedian(NaN, w = 1, na.rm = FALSE), NaN)
  expect_equal(fmedian(Inf, w = 1, na.rm = FALSE), Inf)
  expect_equal(fmedian(-Inf, w = 1, na.rm = FALSE), -Inf)
  expect_equal(fmedian(TRUE, w = 1, na.rm = FALSE), 1)
  expect_equal(fmedian(FALSE, w = 1, na.rm = FALSE), 0)
  expect_equal(fmedian(NA, w = NA), NA_real_)
  expect_equal(fmedian(NaN, w = NA), NA_real_)
  expect_equal(fmedian(Inf, w = NA), NA_real_)
  expect_equal(fmedian(-Inf, w = NA), NA_real_)
  expect_equal(fmedian(TRUE, w = NA), NA_real_)
  expect_equal(fmedian(FALSE, w = NA), NA_real_)
  expect_equal(fmedian(NA, w = NA, na.rm = FALSE), NA_real_)
  expect_equal(fmedian(NaN, w = NA, na.rm = FALSE), NA_real_)
  expect_equal(fmedian(Inf, w = NA, na.rm = FALSE), NA_real_)
  expect_equal(fmedian(-Inf, w = NA, na.rm = FALSE), NA_real_)
  expect_equal(fmedian(TRUE, w = NA, na.rm = FALSE), NA_real_)
  expect_equal(fmedian(FALSE, w = NA, na.rm = FALSE), NA_real_)
  # expect_equal(fmedian(1:3, w = c(1,Inf,3)), 2) # wmedian gives 2 !!!!!!
  # expect_equal(fmedian(1:3, w = c(1,-Inf,3)), 1) # wmedian gives 3 !!!!!!
  # expect_equal(fmedian(1:3, w = c(1,Inf,3), na.rm = FALSE), 2)
  # expect_equal(fmedian(1:3, w = c(1,-Inf,3), na.rm = FALSE), 3)
})

test_that("fmedian produces errors for wrong input", {
  expect_warning(fmedian("a"))
  expect_equal(fmedian(NA_character_), NA_real_)
  expect_error(fmedian(mNAc))
  expect_error(fmedian(mNAc, f))
  expect_error(fmedian(1:2,1:3))
  expect_error(fmedian(m,1:31))
  expect_error(fmedian(mtcars,1:31))
  expect_error(fmedian(mtcars, w = 1:31))
  expect_warning(fmedian("a", w = 1))
  expect_error(fmedian(1:2, w = 1:3))
  expect_equal(fmedian(NA_character_, w = 1), NA_real_)
  expect_error(fmedian(mNAc, w = wdat))
  expect_error(fmedian(mNAc, f, wdat))
  expect_error(fmedian(mNA, w = 1:33))
  expect_error(fmedian(1:2,1:2, 1:3))
  expect_error(fmedian(m,1:32,1:20))
  expect_error(fmedian(mtcars,1:32,1:10))
  expect_error(fmedian(1:2, w = c("a","b")))
  expect_error(fmedian(wlddev))
  expect_error(fmedian(wlddev, w = wlddev$year))
  expect_error(fmedian(wlddev, wlddev$iso3c))
  expect_error(fmedian(wlddev, wlddev$iso3c, wlddev$year))
})

}

# fnth

g <- GRP(mtcars, ~ cyl)
gf <- as_factor_GRP(g)

for (FUN in 1:2) {

  if(FUN == 2L) {
    if(Sys.getenv("OMP") == "TRUE") {
      fnth <- function(x, ...) collapse::fnth(x, ..., nthreads = 2L)
    } else break
  }


test_that("fnth gives a proper lower/upper/average weighted median on complete data", {

  expect_equal(fnth(1:3, w = c(3,1,1), ties = "mean"), 1)
  expect_true(all_identical(
    fnth(1:3, w = c(3,1,1), ties = "mean"),
    fnth(1:3, w = c(3,1,1), ties = "min"),
    fnth(1:3, w = c(3,1,1), ties = "max"),
    fnth(1:3, w = c(3,1,1), ties = "mean", na.rm = FALSE),
    fnth(1:3, w = c(3,1,1), ties = "min", na.rm = FALSE),
    fnth(1:3, w = c(3,1,1), ties = "max", na.rm = FALSE),
    fnth(1:3, g = rep(1,3), w = c(3,1,1), use.g.names = FALSE, ties = "mean"),
    fnth(1:3, g = rep(1,3), w = c(3,1,1), use.g.names = FALSE, ties = "min"),
    fnth(1:3, g = rep(1,3), w = c(3,1,1), use.g.names = FALSE, ties = "max"),
    fnth(1:3, g = rep(1,3), w = c(3,1,1), use.g.names = FALSE, ties = "mean", na.rm = FALSE),
    fnth(1:3, g = rep(1,3), w = c(3,1,1), use.g.names = FALSE, ties = "min", na.rm = FALSE),
    fnth(1:3, g = rep(1,3), w = c(3,1,1), use.g.names = FALSE, ties = "max", na.rm = FALSE)))

  expect_identical(fnth(1:3, w = c(1,1,3), ties = "mean"), 3)
  expect_true(all_identical(
    fnth(1:3, w = c(1,1,3), ties = "mean"),
    fnth(1:3, w = c(1,1,3), ties = "min"),
    fnth(1:3, w = c(1,1,3), ties = "max"),
    fnth(1:3, w = c(1,1,3), ties = "mean", na.rm = FALSE),
    fnth(1:3, w = c(1,1,3), ties = "min", na.rm = FALSE),
    fnth(1:3, w = c(1,1,3), ties = "max", na.rm = FALSE),
    fnth(1:3, g = rep(1,3), w = c(1,1,3), use.g.names = FALSE, ties = "mean"),
    fnth(1:3, g = rep(1,3), w = c(1,1,3), use.g.names = FALSE, ties = "min"),
    fnth(1:3, g = rep(1,3), w = c(1,1,3), use.g.names = FALSE, ties = "max"),
    fnth(1:3, g = rep(1,3), w = c(1,1,3), use.g.names = FALSE, ties = "mean", na.rm = FALSE),
    fnth(1:3, g = rep(1,3), w = c(1,1,3), use.g.names = FALSE, ties = "min", na.rm = FALSE),
    fnth(1:3, g = rep(1,3), w = c(1,1,3), use.g.names = FALSE, ties = "max", na.rm = FALSE)))

  w = c(0.15, 0.1, 0.2, 0.3, 0.25)
  y = seq_len(5) # [order(rnorm(5))]
  expect_identical(fnth(y, w = w, ties = "mean"), 4)
  expect_true(all_identical(4,
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "mean"),
              fnth(y, w = w, ties = "min"),
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min"),
              fnth(y, w = w, ties = "max"),
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max"),
              fnth(y, w = w, na.rm = FALSE, ties = "mean"),
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, na.rm = FALSE, ties = "mean"),
              fnth(y, w = w, ties = "min", na.rm = FALSE),
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min", na.rm = FALSE),
              fnth(y, w = w, ties = "max", na.rm = FALSE),
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max", na.rm = FALSE)))
  w = c(0.15, 0.2, 0.3, 0.25)
  y = seq_len(4) # [order(rnorm(4))]
  expect_identical(fnth(y, w = w, ties = "mean"), 3)
  expect_true(all_identical(3,
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "mean"),
              fnth(y, w = w, ties = "min"),
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min"),
              fnth(y, w = w, ties = "max"),
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max"),
              fnth(y, w = w, na.rm = FALSE, ties = "mean"),
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, na.rm = FALSE, ties = "mean"),
              fnth(y, w = w, ties = "min", na.rm = FALSE),
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min", na.rm = FALSE),
              fnth(y, w = w, ties = "max", na.rm = FALSE),
              fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max", na.rm = FALSE)))

  w = rep(0.25, 4)
  expect_identical(fnth(y, w = w, ties = "mean"), 2.5)
  expect_identical(2.5, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "mean"))
  expect_identical(fnth(y, w = w, ties = "min"), 2)
  expect_identical(2, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min"))
  expect_identical(fnth(y, w = w, ties = "max"), 3)
  expect_identical(3, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max"))
  expect_identical(fnth(y, w = w, na.rm = FALSE, ties = "mean"), 2.5)
  expect_identical(2.5, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, na.rm = FALSE, ties = "mean"))
  expect_identical(fnth(y, w = w, ties = "min", na.rm = FALSE), 2)
  expect_identical(2, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min", na.rm = FALSE))
  expect_identical(fnth(y, w = w, ties = "max", na.rm = FALSE), 3)
  expect_identical(3, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max", na.rm = FALSE))

  w = rep(0.25, 5)
  y = seq_len(5) #[order(rnorm(5))]
  expect_identical(fnth(y, w = w), 3)
  expect_true(all_identical(3,
                            fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "mean"),
                            fnth(y, w = w, ties = "min"),
                            fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min"),
                            fnth(y, w = w, ties = "max"),
                            fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max"),
                            fnth(y, w = w, na.rm = FALSE, ties = "mean"),
                            fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, na.rm = FALSE, ties = "mean"),
                            fnth(y, w = w, ties = "min", na.rm = FALSE),
                            fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min", na.rm = FALSE),
                            fnth(y, w = w, ties = "max", na.rm = FALSE),
                            fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max", na.rm = FALSE)))

  w = c(0.25, 0.25, 0, 0.25, 0.25)
  expect_identical(fnth(y, w = w, ties = "mean"), 3)
  expect_identical(3, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "mean"))
  expect_identical(fnth(y, w = w, ties = "min"), 2)
  expect_identical(2, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min"))
  expect_identical(fnth(y, w = w, ties = "max"), 4)
  expect_identical(4, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max"))
  expect_identical(fnth(y, w = w, na.rm = FALSE, ties = "mean"), 3)
  expect_identical(3, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, na.rm = FALSE, ties = "mean"))
  expect_identical(fnth(y, w = w, ties = "min", na.rm = FALSE), 2)
  expect_identical(2, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min", na.rm = FALSE))
  expect_identical(fnth(y, w = w, ties = "max", na.rm = FALSE), 4)
  expect_identical(4, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max", na.rm = FALSE))

  w = c(0.25, 0.25, 0, 0, 0.25, 0.25)
  y = seq_len(6) # [order(rnorm(6))]
  expect_identical(fnth(y, w = w, ties = "mean"), 3.5)
  expect_identical(3.5, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "mean"))
  expect_identical(fnth(y, w = w, ties = "min"), 2)
  expect_identical(2, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min"))
  expect_identical(fnth(y, w = w, ties = "max"), 5)
  expect_identical(5, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max"))
  expect_identical(fnth(y, w = w, na.rm = FALSE, ties = "mean"), 3.5)
  expect_identical(3.5, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, na.rm = FALSE, ties = "mean"))
  expect_identical(fnth(y, w = w, ties = "min", na.rm = FALSE), 2)
  expect_identical(2, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "min", na.rm = FALSE))
  expect_identical(fnth(y, w = w, ties = "max", na.rm = FALSE), 5)
  expect_identical(5, fnth(y, g = rep(1, length(y)), w = w, use.g.names = FALSE, ties = "max", na.rm = FALSE))


})

test_that("fnth performs like nth (defined above)", {
  n = 2
  expect_error(fnth(NA, n))
  expect_error(fnth(NA, n, na.rm = FALSE))
  expect_error(fnth(1, n))
  expect_equal(fnth(1:3, n), nth(1:3, n, na.rm = TRUE))
  expect_equal(fnth(-1:1, n), nth(-1:1, n, na.rm = TRUE))
  expect_equal(fnth(1:3, n, na.rm = FALSE), nth(1:3, n))
  expect_equal(fnth(-1:1, n, na.rm = FALSE), nth(-1:1, n))
  expect_equal(fnth(x, n), nth(x, n, na.rm = TRUE))
  expect_equal(fnth(x, n, na.rm = FALSE), nth(x, n))
  expect_equal(fnth(xNA, n, na.rm = FALSE), nth(xNA, n))
  expect_equal(fnth(xNA, n), nth(xNA, n, na.rm = TRUE))
  expect_equal(fnth(mtcars, n), fnth(m, n))
  expect_equal(fnth(m, n), dapply(m, nth, n, na.rm = TRUE)) # failed on oldrel-windows-ix86+x86_64
  expect_equal(fnth(m, n, na.rm = FALSE), dapply(m, nth, n)) # failed on oldrel-windows-ix86+x86_64
  expect_equal(fnth(mNA, n, na.rm = FALSE), dapply(mNA, nth, n))
  expect_equal(fnth(mNA, n), dapply(mNA, nth, n, na.rm = TRUE))
  expect_equal(fnth(mtcars, n), dapply(mtcars, nth, n, na.rm = TRUE)) # failed on oldrel-windows-ix86+x86_64
  expect_equal(fnth(mtcars, n, na.rm = FALSE), dapply(mtcars, nth, n)) # failed on oldrel-windows-ix86+x86_64
  expect_equal(fnth(mtcNA, n, na.rm = FALSE), dapply(mtcNA, nth, n))
  expect_equal(fnth(mtcNA, n), dapply(mtcNA, nth, n, na.rm = TRUE))
  f2 <- as.factor(rep(1:10, each = 10)[order(rnorm(100))])
  expect_equal(fnth(x, n, f2), BY(x, f2, nth, n, na.rm = TRUE)) # failed on oldrel-windows-ix86+x86_64
  expect_equal(fnth(x, n, f2, na.rm = FALSE), BY(x, f2, nth, n)) # failed on oldrel-windows-ix86+x86_64
  g2 <- GRP(rep(1:2, each = 16)[order(rnorm(32))])
  expect_equal(fnth(m, n, g2), BY(m, g2, nth, n, na.rm = TRUE)) # failed on oldrel-windows-ix86+x86_64
  expect_equal(fnth(m, n, g2, na.rm = FALSE), BY(m, g2, nth, n)) # failed on oldrel-windows-ix86+x86_64
  expect_equal(fnth(mtcars, n, g2), BY(mtcars, g2, nth, n, na.rm = TRUE)) # failed on oldrel-windows-ix86+x86_64
  expect_equal(fnth(mtcars, n, g2, na.rm = FALSE), BY(mtcars, g2, nth, n)) # failed on oldrel-windows-ix86+x86_64
  for(i in 1:5) {
  n = runif(1, min = 1, max = 999) / 1000 # Probability needed for nth to work with groups
  expect_equal(fnth(1:3, n, ties = "min"), nth(1:3, n, na.rm = TRUE))
  expect_equal(fnth(-1:1, n, ties = "min"), nth(-1:1, n, na.rm = TRUE))
  expect_equal(fnth(1:3, n, na.rm = FALSE, ties = "min"), nth(1:3, n))
  expect_equal(fnth(-1:1, n, na.rm = FALSE, ties = "min"), nth(-1:1, n))
  expect_equal(fnth(x, n, ties = "min"), nth(x, n, na.rm = TRUE))
  expect_equal(fnth(x, n, na.rm = FALSE, ties = "min"), nth(x, n))
  expect_equal(fnth(xNA, n, na.rm = FALSE, ties = "min"), nth(xNA, n))
  expect_equal(fnth(xNA, n, ties = "min"), nth(xNA, n, na.rm = TRUE))
  expect_equal(fnth(mtcars, n, ties = "min"), fnth(m, n, ties = "min"))
  expect_equal(fnth(m, n, ties = "min"), dapply(m, nth, n, na.rm = TRUE))
  expect_equal(fnth(m, n, na.rm = FALSE, ties = "min"), dapply(m, nth, n))
  expect_equal(fnth(mNA, n, na.rm = FALSE, ties = "min"), dapply(mNA, nth, n))
  expect_equal(fnth(mNA, n, ties = "min"), dapply(mNA, nth, n, na.rm = TRUE))
  expect_equal(fnth(mtcars, n, ties = "min"), dapply(mtcars, nth, n, na.rm = TRUE))
  expect_equal(fnth(mtcars, n, na.rm = FALSE, ties = "min"), dapply(mtcars, nth, n))
  expect_equal(fnth(mtcNA, n, na.rm = FALSE, ties = "min"), dapply(mtcNA, nth, n))
  expect_equal(fnth(mtcNA, n, ties = "min"), dapply(mtcNA, nth, n, na.rm = TRUE))
  expect_equal(fnth(xNA, n, f2, na.rm = FALSE, ties = "min"), BY(xNA, f2, nth, n))
  expect_equal(fnth(xNA, n, f2, ties = "min"), BY(xNA, f2, nth, n, na.rm = TRUE))
  expect_equal(fnth(m, n, g, ties = "min"), BY(m, g, nth, n, na.rm = TRUE))
  expect_equal(fnth(m, n, g, na.rm = FALSE, ties = "min"), BY(m, g, nth, n))
  expect_equal(fnth(mNA, n, g, na.rm = FALSE, ties = "min"), BY(mNA, g, nth, n))
  expect_equal(fnth(mNA, n, g, ties = "min"), BY(mNA, g, nth, n, na.rm = TRUE))
  expect_equal(fnth(mtcars, n, g, ties = "min"), BY(mtcars, g, nth, n, na.rm = TRUE))
  expect_equal(fnth(mtcars, n, g, na.rm = FALSE, ties = "min"), BY(mtcars, g, nth, n))
  expect_equal(fnth(mtcNA, n, g, na.rm = FALSE, ties = "min"), BY(mtcNA, g, nth, n))
  expect_equal(fnth(mtcNA, n, g, ties = "min"), BY(mtcNA, g, nth, n, na.rm = TRUE))
  }
})

test_that("fnth matrix and data.frame method work alike", {
  for(i in 1:3) {
  n = runif(1, min = 1, max = 999) / 1000
  expect_equal(fnth(mtcars, n, ties = "min"), fnth(m, n, ties = "min"))
  expect_equal(fnth(mtcars, n), fnth(m, n))
  expect_equal(fnth(mtcars, n, ties = "max"), fnth(m, n, ties = "max"))
  expect_equal(fnth(mtcNA, n, ties = "min"), fnth(mNA, n, ties = "min"))
  expect_equal(fnth(mtcNA, n), fnth(mNA, n))
  expect_equal(fnth(mtcNA, n, ties = "max"), fnth(mNA, n, ties = "max"))
  expect_equal(qM(fnth(mtcars, n, g, ties = "min")), fnth(m, n, g, ties = "min"))
  expect_equal(qM(fnth(mtcars, n, g)), fnth(m, n, g))
  expect_equal(qM(fnth(mtcars, n, g, ties = "max")), fnth(m, n, g, ties = "max"))
  expect_equal(qM(fnth(mtcNA, n, g, ties = "min")), fnth(mNA, n, g, ties = "min"))
  expect_equal(qM(fnth(mtcNA, n, g)), fnth(mNA, n, g))
  expect_equal(qM(fnth(mtcNA, n, g, ties = "max")), fnth(mNA, n, g, ties = "max"))

  expect_equal(fnth(mtcars, n, w = wdat, ties = "min"), fnth(m, n, w = wdat, ties = "min"))
  expect_equal(fnth(mtcars, n, w = wdat), fnth(m, n, w = wdat))
  expect_equal(fnth(mtcars, n, w = wdat, ties = "max"), fnth(m, n, w = wdat, ties = "max"))
  expect_equal(fnth(mtcNA, n, w = wdat, ties = "min"), fnth(mNA, n, w = wdat, ties = "min"))
  expect_equal(fnth(mtcNA, n, w = wdat), fnth(mNA, n, w = wdat))
  expect_equal(fnth(mtcNA, n, w = wdat, ties = "max"), fnth(mNA, n, w = wdat, ties = "max"))
  expect_equal(qM(fnth(mtcars, n, g, wdat, ties = "min")), fnth(m, n, g, wdat, ties = "min"))
  expect_equal(qM(fnth(mtcars, n, g, wdat)), fnth(m, n, g, wdat))
  expect_equal(qM(fnth(mtcars, n, g, wdat, ties = "max")), fnth(m, n, g, wdat, ties = "max"))
  expect_equal(qM(fnth(mtcNA, n, g, wdat, ties = "min")), fnth(mNA, n, g, wdat, ties = "min"))
  expect_equal(qM(fnth(mtcNA, n, g, wdat)), fnth(mNA, n, g, wdat))
  expect_equal(qM(fnth(mtcNA, n, g, wdat, ties = "max")), fnth(mNA, n, g, wdat, ties = "max"))
 }
})

test_that("fnth performs like fnth with weights all equal", {
  for(t in c("min","max")) { # "mean", # already tested above..
    # for(i in 1:3) {
      n = 0.5 # round(runif(1, min = 1, max = 999) / 1000, 3) # other numbers than 0.5 do not work and cannot work..
      expect_equal(fnth(NA, n, ties = t), fnth(NA, n, w = 1, ties = t))
      expect_equal(fnth(NA, n, na.rm = FALSE, ties = t), fnth(NA, n, w = 1, na.rm = FALSE, ties = t))
      expect_equal(fnth(1, n, ties = t), fnth(1, n, w = 3, ties = t))
      expect_equal(fnth(1:3, n, ties = t), fnth(1:3, n, w = rep(1,3), ties = t))
      expect_equal(fnth(-1:1, n, ties = t), fnth(-1:1, n, w = rep(4.2,3), ties = t))
      expect_equal(fnth(1, n, na.rm = FALSE, ties = t), fnth(1, n, w = 5, na.rm = FALSE, ties = t))
      expect_equal(fnth(1:3, n, na.rm = FALSE, ties = t), fnth(1:3, n, w = rep(1, 3), na.rm = FALSE, ties = t))
      expect_equal(fnth(-1:1, n, na.rm = FALSE, ties = t), fnth(-1:1, n, w = rep(12, 3), na.rm = FALSE, ties = t))
      expect_equal(fnth(x, n, ties = t), fnth(x, n, w = rep(1,100), ties = t))
      expect_equal(fnth(x, n, na.rm = FALSE, ties = t), fnth(x, n, w = rep(1, 100), na.rm = FALSE, ties = t))
      expect_equal(fnth(xNA, n, na.rm = FALSE, ties = t), fnth(xNA, n, w = rep(5, 100), na.rm = FALSE, ties = t))
      expect_equal(fnth(xNA, n, ties = t), fnth(xNA, n, w = rep(4, 100), ties = t))
      expect_equal(fnth(m, n, ties = t), fnth(m, n, w = rep(6587, 32), ties = t))
      expect_equal(fnth(m, n, na.rm = FALSE, ties = t), fnth(m, n, w = rep(6587, 32), na.rm = FALSE, ties = t))
      expect_equal(fnth(mNA, n, na.rm = FALSE, ties = t), fnth(mNA, n, w = rep(6587, 32), na.rm = FALSE, ties = t))
      expect_equal(fnth(mNA, n, ties = t), fnth(mNA, n, w = rep(6587, 32), ties = t))
      expect_equal(fnth(mtcars, n, ties = t), fnth(mtcars, n, w = rep(6787, 32), ties = t))
      expect_equal(fnth(mtcars, n, na.rm = FALSE, ties = t), fnth(mtcars, n, w = rep(6787, 32), na.rm = FALSE, ties = t))
      expect_equal(fnth(mtcNA, n, na.rm = FALSE, ties = t), fnth(mtcNA, n, w = rep(6787, 32), na.rm = FALSE, ties = t))
      expect_equal(fnth(mtcNA, n, ties = t), fnth(mtcNA, n, w = rep(6787, 32), ties = t))
      expect_equal(fnth(x, n, f, ties = t), fnth(x, n, f, rep(547,100), ties = t))
      expect_equal(fnth(x, n, f, na.rm = FALSE, ties = t), fnth(x, n, f, rep(6, 100), na.rm = FALSE, ties = t))
      expect_equal(fnth(xNA, n, f, na.rm = FALSE, ties = t), fnth(xNA, n, f, rep(52,100), na.rm = FALSE, ties = t))
      expect_equal(fnth(xNA, n, f, ties = t), fnth(xNA, n, f, rep(5997456,100), ties = t))
      expect_equal(fnth(m, n, g, ties = t), fnth(m, n, g, rep(546,32), ties = t))
      expect_equal(fnth(m, n, g, na.rm = FALSE, ties = t), fnth(m, n, g, rep(1,32), na.rm = FALSE, ties = t))
      expect_equal(fnth(mNA, n, g, na.rm = FALSE, ties = t), fnth(mNA, n, g, rep(5,32), na.rm = FALSE, ties = t))
      expect_equal(fnth(mNA, n, g, ties = t), fnth(mNA, n, g, rep(1,32), ties = t))
      expect_equal(fnth(mtcars, n, g, ties = t), fnth(mtcars, n, g, rep(53,32), ties = t))
      expect_equal(fnth(mtcars, n, g, na.rm = FALSE, ties = t), fnth(mtcars, n, g, rep(546,32), na.rm = FALSE, ties = t))
      expect_equal(fnth(mtcNA, n, g, na.rm = FALSE, ties = t), fnth(mtcNA, n, g, rep(1,32), na.rm = FALSE, ties = t))
      expect_equal(fnth(mtcNA, n, g, ties = t), fnth(mtcNA, n, g, rep(999,32), ties = t))
    #}
  }
})

test_that("fnth with weights performs like wnth (defined above)", {
  for(t in c("mean","min","max")) {
    # print(t)
    for(i in 1:3) {
    n = round(runif(1, min = 1, max = 999) / 1000, 3)
    # complete weights
    expect_equal(fnth(NA, n, w = 1, ties = t), wnth(NA_real_, n, 1, ties = t))
    expect_equal(fnth(NA, n, w = 1, na.rm = FALSE, ties = t), wnth(NA_real_, n, 1, ties = t))
    expect_equal(fnth(1, n, w = 1, ties = t), wnth(1, n, w = 1, ties = t))
    expect_equal(fnth(1:3, n, w = 1:3, ties = t), wnth(1:3, n, 1:3, ties = t))
    expect_equal(fnth(-1:1, n, w = 1:3, ties = t), wnth(-1:1, n, 1:3, ties = t))
    expect_equal(fnth(1, n, w = 1, na.rm = FALSE, ties = t), wnth(1, n, 1, ties = t))
    expect_equal(fnth(1:3, n, w = c(0.99,3454,1.111), na.rm = FALSE, ties = t), wnth(1:3, n, c(0.99,3454,1.111), ties = t))
    expect_equal(fnth(-1:1, n, w = 1:3, na.rm = FALSE, ties = t), wnth(-1:1, n, 1:3, ties = t))
    expect_equal(fnth(x, n, w = w, ties = t), wnth(x, n, w, ties = t))
    expect_equal(fnth(x, n, w = w, na.rm = FALSE, ties = t), wnth(x, n, w, ties = t))
    expect_equal(fnth(xNA, n, w = w, na.rm = FALSE, ties = t), wnth(xNA, n, w, ties = t))
    expect_equal(fnth(xNA, n, w = w, ties = t), wnth(xNA, n, w, na.rm = TRUE, ties = t))
    expect_equal(fnth(mtcars, n, w = wdat, ties = t), fnth(m, n, w = wdat, ties = t))
    expect_equal(fnth(m, n, w = wdat, ties = t), dapply(m, wnth, n, wdat, na.rm = TRUE, ties = t))
    expect_equal(fnth(m, n, w = wdat, na.rm = FALSE, ties = t), dapply(m, wnth, n, wdat, ties = t))
    expect_equal(fnth(mNA, n, w = wdat, na.rm = FALSE, ties = t), dapply(mNA, wnth, n, wdat, ties = t))
    expect_equal(fnth(mNA, n, w = wdat, ties = t), dapply(mNA, wnth, n, wdat, na.rm = TRUE, ties = t))
    expect_equal(fnth(mtcars, n, w = wdat, ties = t), dapply(mtcars, wnth, n, wdat, na.rm = TRUE, ties = t))
    expect_equal(fnth(mtcars, n, w = wdat, na.rm = FALSE, ties = t), dapply(mtcars, wnth, n, wdat, ties = t))
    expect_equal(fnth(mtcNA, n, w = wdat, na.rm = FALSE, ties = t), dapply(mtcNA, wnth, n, wdat, ties = t))
    expect_equal(fnth(mtcNA, n, w = wdat, ties = t), dapply(mtcNA, wnth, n, wdat, na.rm = TRUE, ties = t))
    expect_equal(fnth(x, n, f, w, ties = t), BY(x, f, wnth, n = n, w = w, ties = t))
    expect_equal(fnth(x, n, f, w, na.rm = FALSE, ties = t), BY(x, f, wnth, n = n, w = w, ties = t))
    expect_equal(fnth(xNA, n, f, w, na.rm = FALSE, ties = t), BY(xNA, f, wnth, n = n, w = w, ties = t))
    expect_equal(fnth(xNA, n, f, w, ties = t), BY(xNA, f, wnth, n = n, w = w, na.rm = TRUE, ties = t))
    expect_equal(fnth(m, n, g, wdat, ties = t), BY(m, gf, wnth, n = n, w = wdat, ties = t))
    expect_equal(fnth(m, n, g, wdat, na.rm = FALSE, ties = t), BY(m, gf, wnth, n = n, w = wdat, ties = t))
    expect_equal(fnth(mNA, n, g, wdat, na.rm = FALSE, ties = t),  BY(mNA, gf, wnth, n = n, w = wdat, ties = t))
    expect_equal(fnth(mNA, n, g, wdat, ties = t), BY(mNA, gf, wnth, n = n, w = wdat, na.rm = TRUE, ties = t))
    expect_equal(fnth(mtcars, n, g, wdat, ties = t), BY(mtcars, gf, wnth, n = n, w = wdat, ties = t))
    expect_equal(fnth(mtcars, n, g, wdat, na.rm = FALSE, ties = t), BY(mtcars, gf, wnth, n = n, w = wdat, ties = t))
    expect_equal(fnth(mtcNA, n, g, wdat, na.rm = FALSE, ties = t), BY(mtcNA, gf, wnth, n = n, w = wdat, ties = t))
    expect_equal(fnth(mtcNA, n, g, wdat, ties = t), BY(mtcNA, gf, wnth, w = wdat, n = n, na.rm = TRUE, ties = t))
    # missing weights: Only supported if x is also missing...
    expect_equal(fnth(NA, n, w = NA, ties = t), wnth(NA_real_, n, NA_real_, ties = t))
    expect_equal(fnth(NA, n, w = NA, na.rm = FALSE, ties = t), wnth(NA_real_, n, NA_real_, ties = t))
    expect_equal(fnth(xNA, n, w = wNA, na.rm = FALSE, ties = t), wnth(xNA, n, wNA, ties = t))
    expect_equal(fnth(xNA, n, w = wNA, ties = t), wnth(xNA, n, wNA, na.rm = TRUE, ties = t))
    expect_equal(fnth(xNA, n, f, wNA, na.rm = FALSE, ties = t), BY(xNA, f, wnth, n = n, w = w, ties = t))
    expect_equal(fnth(xNA, n, f, wNA, ties = t), BY(xNA, f, wnth, n = n, w = w, na.rm = TRUE, ties = t))
    }
  }
})

test_that("fnth properly deals with missing data", {
  expect_equal(fnth(NA), NA_real_)
  expect_equal(fnth(NA, na.rm = FALSE), NA_real_)
  expect_equal(fnth(rep(NA, 2), w = 1:2), NA_real_)
  expect_equal(fnth(rep(NA, 2), w = 1:2), NA_real_)
  expect_equal(fnth(NA, w = 1), NA_real_)
  expect_equal(fnth(NA, w = 1, na.rm = FALSE), NA_real_)
  expect_equal(fnth(1), 1)
  expect_equal(fnth(1, na.rm = FALSE), 1)
  expect_error(fnth(1:2, w = rep(NA, 2)))
  expect_error(fnth(1:2, w = c(1, NA)))
  expect_error(fnth(1:2, w = c(NA, 1)))
})

}

