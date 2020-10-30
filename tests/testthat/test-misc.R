context("Misc")

# rm(list = ls())

m <- na_insert(qM(mtcars))

test_that("descr, pwcor, pwcov, pwNobs", {
  expect_visible(descr(wlddev))
  expect_visible(as.data.frame(descr(wlddev)))
  expect_output(print(descr(wlddev)))
  expect_visible(descr(GGDC10S))
  expect_output(print(pwcor(nv(wlddev))))
  expect_output(print(pwcor(nv(wlddev), N = TRUE)))
  expect_output(print(pwcor(nv(wlddev), P = TRUE)))
  expect_output(print(pwcor(nv(wlddev), N = TRUE, P = TRUE)))
  expect_output(print(pwcor(nv(wlddev), N = TRUE, P = TRUE, use = "complete.obs")))
  expect_visible(pwcor(nv(GGDC10S)))
  expect_visible(pwcov(nv(wlddev)))
  expect_output(print(pwcov(nv(wlddev))))
  expect_output(print(pwcov(nv(wlddev), N = TRUE)))
  expect_output(print(pwcov(nv(wlddev), P = TRUE)))
  expect_output(print(pwcov(nv(wlddev), N = TRUE, P = TRUE)))
  expect_output(print(pwcov(nv(wlddev), N = TRUE, P = TRUE, use = "complete.obs")))

  expect_visible(pwNobs(wlddev))
  expect_visible(pwNobs(GGDC10S))

  expect_visible(descr(m))
  expect_visible(pwcor(m))
  expect_visible(pwcov(m))
  expect_visible(pwNobs(m))

})

test_that("deep matrix dispatch works well", {
  tsm <- EuStockMarkets
  class(tsm) <- setdiff(class(tsm), "matrix")
  f <- qF(sample.int(5, nrow(tsm), TRUE))
  NCOL2 <- function(x) if(length(d <- dim(x)) > 1L) d[2L] else length(x)

  for(i in setdiff(c(.FAST_FUN, .OPERATOR_FUN), c("fnth","flag","L","F", "fdiff","D","Dlog", "fgrowth","G")))
      expect_equal(NCOL2(match.fun(i)(tsm, f)), 4L)

  expect_equal(NCOL2(fnth(tsm, 0.5, f)), 4L)
  expect_equal(NCOL2(BY(tsm, f, sum)), 4L)
  expect_equal(nrow(qsu(tsm)), 4L)

  for(i in c("flag", "L", "fdiff", "D", "Dlog", "fgrowth", "G"))
      expect_true(all(is.na(match.fun(i)(tsm)[1L, ])))

})
