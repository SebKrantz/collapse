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

m <- qM(mtcars)
v <- mtcars$mpg
f <- qF(mtcars$cyl)
fcc <- qF(mtcars$cyl, na.exclude = FALSE)
g <- GRP(mtcars, ~ cyl)
gl <- mtcars["cyl"]
gmtc <- fgroup_by(mtcars, cyl)

test_that("fast functions give same result using different grouping mechanisms", {

 for(i in .FAST_STAT_FUN) {

   FUN <- match.fun(i)
   expect_true(all_obj_equal(FUN(v, g = mtcars$cyl), FUN(v, g = f), FUN(v, g = fcc), FUN(v, g = g), FUN(v, g = gl)))
   expect_true(all_obj_equal(FUN(m, g = mtcars$cyl), FUN(m, g = f), FUN(m, g = fcc), FUN(m, g = g), FUN(m, g = gl)))
   expect_true(all_obj_equal(FUN(mtcars, g = mtcars$cyl), FUN(mtcars, g = f), FUN(mtcars, g = fcc), FUN(mtcars, g = g), FUN(mtcars, g = gl)))
   expect_true(all_obj_equal(gv(FUN(mtcars, g = mtcars$cyl, use.g.names = FALSE), -2), gv(FUN(gmtc), -1), FUN(gmtc, keep.group_vars = FALSE)))

   expect_true(all_obj_equal(FUN(v, g = mtcars$cyl, TRA = 1L), TRA(v, FUN(v, g = mtcars$cyl), 1L, mtcars$cyl),
                             FUN(v, g = f, TRA = 1L), TRA(v, FUN(v, g = f), 1L, f),
                             FUN(v, g = fcc, TRA = 1L), TRA(v, FUN(v, g = fcc), 1L, fcc),
                             FUN(v, g = g, TRA = 1L), TRA(v, FUN(v, g = g), 1L, g),
                             FUN(v, g = gl, TRA = 1L), TRA(v, FUN(v, g = gl), 1L, gl)))

   expect_true(all_obj_equal(FUN(m, g = mtcars$cyl, TRA = 1L), TRA(m, FUN(m, g = mtcars$cyl), 1L, mtcars$cyl),
                             FUN(m, g = f, TRA = 1L), TRA(m, FUN(m, g = f), 1L, f),
                             FUN(m, g = fcc, TRA = 1L), TRA(m, FUN(m, g = fcc), 1L, fcc),
                             FUN(m, g = g, TRA = 1L), TRA(m, FUN(m, g = g), 1L, g),
                             FUN(m, g = gl, TRA = 1L), TRA(m, FUN(m, g = gl), 1L, gl)))

   expect_true(all_obj_equal(FUN(mtcars, g = mtcars$cyl, TRA = 1L), TRA(mtcars, FUN(mtcars, g = mtcars$cyl), 1L, mtcars$cyl),
                             FUN(mtcars, g = f, TRA = 1L), TRA(mtcars, FUN(mtcars, g = f), 1L, f),
                             FUN(mtcars, g = fcc, TRA = 1L), TRA(mtcars, FUN(mtcars, g = fcc), 1L, fcc),
                             FUN(mtcars, g = g, TRA = 1L), TRA(mtcars, FUN(mtcars, g = g), 1L, g),
                             FUN(mtcars, g = gl, TRA = 1L), TRA(mtcars, FUN(mtcars, g = gl), 1L, gl)))

   expect_equal(colorder(FUN(gmtc, TRA = 1L), mpg, cyl), TRA(gmtc, FUN(gmtc), 1L))
   expect_equal(FUN(fselect(gmtc, -cyl), TRA = 1L), TRA(fselect(gmtc, -cyl), FUN(gmtc, keep.group_vars = FALSE), 1L))
 }

  for(i in setdiff(.FAST_FUN, c(.FAST_STAT_FUN, "fHDbetween", "fHDwithin"))) {

    FUN <- match.fun(i)
    expect_true(all_obj_equal(FUN(v, g = mtcars$cyl), FUN(v, g = f), FUN(v, g = fcc), FUN(v, g = g), FUN(v, g = gl)))
    expect_true(all_obj_equal(FUN(m, g = mtcars$cyl), FUN(m, g = f), FUN(m, g = fcc), FUN(m, g = g), FUN(m, g = gl)))
    expect_true(all_obj_equal(FUN(mtcars, g = mtcars$cyl), FUN(mtcars, g = f), FUN(mtcars, g = fcc), FUN(mtcars, g = g), FUN(mtcars, g = gl)))

  }

  for(i in c("STD", "B", "W", "L", "D", "Dlog", "G")) {

    FUN <- match.fun(i)
    expect_true(all_obj_equal(FUN(v, g = mtcars$cyl), FUN(v, g = f), FUN(v, g = fcc), FUN(v, g = g), FUN(v, g = gl)))
    expect_true(all_obj_equal(FUN(m, g = mtcars$cyl), FUN(m, g = f), FUN(m, g = fcc), FUN(m, g = g), FUN(m, g = gl)))
    expect_true(all_obj_equal(FUN(mtcars, by = mtcars$cyl), FUN(mtcars, by = f), FUN(mtcars, by = fcc), FUN(mtcars, by = g), FUN(mtcars, by = gl)))

  }

})
