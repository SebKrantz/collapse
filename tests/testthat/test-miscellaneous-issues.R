context("miscellaneous issues")



# rm(list = ls())

options(warn = -1)

F <- getNamespace("collapse")$F

if(identical(Sys.getenv("NCRAN"), "TRUE")) {

test_that("Using a factor with unused levels does not pose a problem to flag, fdiff or fgrowth (#25)", {
  wlddev2 <- subset(wlddev, iso3c %in% c("ALB", "AFG", "DZA"))
  wlddev3 <- droplevels(wlddev2)
  expect_identical(L(wlddev3, 1, LIFEEX~iso3c, ~year), L(wlddev3, 1, ~iso3c, ~year, cols="LIFEEX"))
  expect_identical(L(wlddev3, -1:1, LIFEEX~iso3c, ~year), L(wlddev3, -1:1, ~iso3c, ~year, cols="LIFEEX"))
  expect_identical(droplevels(L(wlddev2, 1, ~iso3c, ~year, cols="LIFEEX")), L(wlddev3, 1, ~iso3c, ~year, cols="LIFEEX"))
  expect_identical(droplevels(L(wlddev2, -1:1, ~iso3c, ~year, cols="LIFEEX")), L(wlddev3, -1:1, ~iso3c, ~year, cols="LIFEEX"))
  expect_identical(droplevels(D(wlddev2, 1, 1, ~iso3c, ~year, cols="LIFEEX")), D(wlddev3, 1, 1, ~iso3c, ~year, cols="LIFEEX"))
  expect_identical(droplevels(D(wlddev2, -1:1, 1:2, ~iso3c, ~year, cols="LIFEEX")), D(wlddev3, -1:1, 1:2, ~iso3c, ~year, cols="LIFEEX"))
  expect_identical(droplevels(Dlog(wlddev2, 1, 1, ~iso3c, ~year, cols="LIFEEX")), Dlog(wlddev3, 1, 1, ~iso3c, ~year, cols="LIFEEX"))
  expect_identical(droplevels(Dlog(wlddev2, -1:1, 1:2, ~iso3c, ~year, cols="LIFEEX")), Dlog(wlddev3, -1:1, 1:2, ~iso3c, ~year, cols="LIFEEX"))
  expect_identical(droplevels(D(wlddev2, 1, 1, ~iso3c, ~year, cols="LIFEEX", rho = 0.95)), D(wlddev3, 1, 1, ~iso3c, ~year, cols="LIFEEX", rho = 0.95))
  expect_identical(droplevels(D(wlddev2, -1:1, 1:2, ~iso3c, ~year, cols="LIFEEX", rho = 0.95)), D(wlddev3, -1:1, 1:2, ~iso3c, ~year, cols="LIFEEX", rho = 0.95))
  expect_identical(droplevels(Dlog(wlddev2, 1, 1, ~iso3c, ~year, cols="LIFEEX", rho = 0.95)), Dlog(wlddev3, 1, 1, ~iso3c, ~year, cols="LIFEEX", rho = 0.95))
  expect_identical(droplevels(Dlog(wlddev2, -1:1, 1:2, ~iso3c, ~year, cols="LIFEEX", rho = 0.95)), Dlog(wlddev3, -1:1, 1:2, ~iso3c, ~year, cols="LIFEEX", rho = 0.95))
  expect_identical(droplevels(G(wlddev2, 1, 1, ~iso3c, ~year, cols="LIFEEX")), G(wlddev3, 1, 1, ~iso3c, ~year, cols="LIFEEX"))
  expect_identical(droplevels(G(wlddev2, -1:1, 1:2, ~iso3c, ~year, cols="LIFEEX")), G(wlddev3, -1:1, 1:2, ~iso3c, ~year, cols="LIFEEX"))

  expect_identical(L(wlddev3, 1, LIFEEX~iso3c), L(wlddev3, 1, ~iso3c, cols="LIFEEX"))
  expect_identical(L(wlddev3, -1:1, LIFEEX~iso3c), L(wlddev3, -1:1, ~iso3c, cols="LIFEEX"))
  expect_identical(droplevels(L(wlddev2, 1, ~iso3c, cols="LIFEEX")), L(wlddev3, 1, ~iso3c, cols="LIFEEX"))
  expect_identical(droplevels(L(wlddev2, -1:1, ~iso3c, cols="LIFEEX")), L(wlddev3, -1:1, ~iso3c, cols="LIFEEX"))
  expect_identical(droplevels(D(wlddev2, 1, 1, ~iso3c, cols="LIFEEX")), D(wlddev3, 1, 1, ~iso3c, cols="LIFEEX"))
  expect_identical(droplevels(D(wlddev2, -1:1, 1:2, ~iso3c, cols="LIFEEX")), D(wlddev3, -1:1, 1:2, ~iso3c, cols="LIFEEX"))
  expect_identical(droplevels(Dlog(wlddev2, 1, 1, ~iso3c, cols="LIFEEX")), Dlog(wlddev3, 1, 1, ~iso3c, cols="LIFEEX"))
  expect_identical(droplevels(Dlog(wlddev2, -1:1, 1:2, ~iso3c, cols="LIFEEX")), Dlog(wlddev3, -1:1, 1:2, ~iso3c, cols="LIFEEX"))
  expect_identical(droplevels(D(wlddev2, 1, 1, ~iso3c, cols="LIFEEX", rho = 0.95)), D(wlddev3, 1, 1, ~iso3c, cols="LIFEEX", rho = 0.95))
  expect_identical(droplevels(D(wlddev2, -1:1, 1:2, ~iso3c, cols="LIFEEX", rho = 0.95)), D(wlddev3, -1:1, 1:2, ~iso3c, cols="LIFEEX", rho = 0.95))
  expect_identical(droplevels(Dlog(wlddev2, 1, 1, ~iso3c, cols="LIFEEX", rho = 0.95)), Dlog(wlddev3, 1, 1, ~iso3c, cols="LIFEEX", rho = 0.95))
  expect_identical(droplevels(Dlog(wlddev2, -1:1, 1:2, ~iso3c, cols="LIFEEX", rho = 0.95)), Dlog(wlddev3, -1:1, 1:2, ~iso3c, cols="LIFEEX", rho = 0.95))
  expect_identical(droplevels(G(wlddev2, 1, 1, ~iso3c, cols="LIFEEX")), G(wlddev3, 1, 1, ~iso3c, cols="LIFEEX"))
  expect_identical(droplevels(G(wlddev2, -1:1, 1:2, ~iso3c, cols="LIFEEX")), G(wlddev3, -1:1, 1:2, ~iso3c, cols="LIFEEX"))

})

test_that("Using a factor with unused levels does not pose a problem to statistical functions", {
    wlddev2 <- fsubset(wlddev, iso3c %in% c("ALB", "AFG", "DZA"))
    d <- nv(wlddev2)
    m <- qM(d)
    v <- d$PCGDP
    w <- rep(1, length(v))
    f <- wlddev2$iso3c
    lev <- levels(f)
    fd <- fdroplevels(f)
    levd <- levels(fd)

    # Testing BY:
    expect_equal(attr(BY(d, f, sum), "row.names"), lev)
    expect_equal(dimnames(BY(m, f, sum))[[1L]], lev)
    expect_equal(names(BY(v, f, sum)), lev)

    # Fast Statistical Functions
    for(i in .FAST_STAT_FUN) {
      # print(i)
      FUN <- match.fun(i)
      expect_equal(attr(FUN(d, g = f), "row.names"), lev)
      expect_equal(dimnames(FUN(m, g = f))[[1L]], lev)
      expect_equal(names(FUN(v, g = f)), lev)
      expect_equal(attr(FUN(d, g = fd), "row.names"), levd)
      expect_equal(dimnames(FUN(m, g = fd))[[1L]], levd)
      expect_equal(names(FUN(v, g = fd)), levd)
      if(i != "fnobs") {
        expect_equal(attr(FUN(d, g = f, na.rm = FALSE), "row.names"), lev)
        expect_equal(dimnames(FUN(m, g = f, na.rm = FALSE))[[1L]], lev)
        expect_equal(names(FUN(v, g = f, na.rm = FALSE)), lev)
      }
      if(i %in% c("fsum", "fprod", "fmean", "fmedian", "fnth", "fmode", "fvar", "fsd")) {
        expect_equal(attr(FUN(d, g = f, w = w), "row.names"), lev)
        expect_equal(dimnames(FUN(m, g = f, w = w))[[1L]], lev)
        expect_equal(names(FUN(v, g = f, w = w)), lev)
        expect_equal(attr(FUN(d, g = f, w = w, na.rm = FALSE), "row.names"), lev)
        expect_equal(dimnames(FUN(m, g = f, w = w, na.rm = FALSE))[[1L]], lev)
        expect_equal(names(FUN(v, g = f, w = w, na.rm = FALSE)), lev)
        expect_equal(FUN(d, g = f, w = w), FUN(d, g = f))
        expect_equal(FUN(m, g = f, w = w), FUN(m, g = f))
        expect_equal(FUN(v, g = f, w = w), FUN(v, g = f))
      }
    }

    # Other Statistical Functions
    for(i in setdiff(c(.FAST_FUN, .OPERATOR_FUN), .FAST_STAT_FUN)) {
      # print(i)
      FUN <- match.fun(i)
      if(grepl("hd", i, ignore.case = TRUE)) {
        expect_equal(FUN(d, fl = f), FUN(d, fl = fd))
        expect_equal(FUN(m, fl = f), FUN(m, fl = fd))
        expect_equal(FUN(v, fl = f), FUN(v, fl = fd))
        expect_equal(FUN(d, fl = f, na.rm = FALSE), FUN(d, fl = fd, na.rm = FALSE))
        expect_equal(FUN(m, fl = f, na.rm = FALSE), FUN(m, fl = fd, na.rm = FALSE))
        expect_equal(FUN(v, fl = f, na.rm = FALSE), FUN(v, fl = fd, na.rm = FALSE))
        expect_equal(FUN(d, fl = f, w = w), FUN(d, fl = fd))
        expect_equal(FUN(m, fl = f, w = w), FUN(m, fl = fd))
        expect_equal(FUN(v, fl = f, w = w), FUN(v, fl = fd))
        expect_equal(FUN(d, fl = f, w = w, na.rm = FALSE), FUN(d, fl = fd, na.rm = FALSE))
        expect_equal(FUN(m, fl = f, w = w, na.rm = FALSE), FUN(m, fl = fd, na.rm = FALSE))
        expect_equal(FUN(v, fl = f, w = w, na.rm = FALSE), FUN(v, fl = fd, na.rm = FALSE))
      } else {
        expect_equal(FUN(d, g = f), FUN(d, g = fd))
        expect_equal(FUN(m, g = f), FUN(m, g = fd))
        expect_equal(FUN(v, g = f), FUN(v, g = fd))
        expect_equal(FUN(d, g = f, na.rm = FALSE), FUN(d, g = fd, na.rm = FALSE))
        expect_equal(FUN(m, g = f, na.rm = FALSE), FUN(m, g = fd, na.rm = FALSE))
        expect_equal(FUN(v, g = f, na.rm = FALSE), FUN(v, g = fd, na.rm = FALSE))
        if(i %in% c("fscale", "STD", "fbetween", "B", "fwithin", "W")) {
          expect_equal(FUN(d, g = f, w = w), FUN(d, g = fd))
          expect_equal(FUN(m, g = f, w = w), FUN(m, g = fd))
          expect_equal(FUN(v, g = f, w = w), FUN(v, g = fd))
          expect_equal(FUN(d, g = f, w = w, na.rm = FALSE), FUN(d, g = fd, na.rm = FALSE))
          expect_equal(FUN(m, g = f, w = w, na.rm = FALSE), FUN(m, g = fd, na.rm = FALSE))
          expect_equal(FUN(v, g = f, w = w, na.rm = FALSE), FUN(v, g = fd, na.rm = FALSE))
        }
      }
    }
})

test_that("Testing grouped_df methods", {
  skip_if_not_installed("magrittr")
  library(magrittr)
  for(sortg in c(TRUE, FALSE)) {
    for(retgrp in c(TRUE, FALSE)) {
      gdf <- wlddev %>% fsubset(year > 1990, region, income, PCGDP:ODA) %>% fgroup_by(region, income, return.groups = retgrp, sort = sortg)
      gdf[["wgt"]] <- round(abs(10*rnorm(fnrow(gdf))), 1)
      expect_visible(gdf %>% fmean)
      expect_visible(gdf %>% fmean(wgt))
      expect_equal(gdf %>% fmean(wgt) %>% slt(-sum.wgt), gdf %>% fmean(wgt, keep.w = FALSE))
      expect_visible(gdf %>% fmedian)
      expect_visible(gdf %>% fmedian(wgt))
      expect_equal(gdf %>% fmedian(wgt) %>% slt(-sum.wgt), gdf %>% fmedian(wgt, keep.w = FALSE))
      expect_visible(gdf %>% fnth)
      expect_visible(gdf %>% fnth(0.75))
      expect_visible(gdf %>% fnth(0.75, wgt))
      expect_equal(gdf %>% fnth(0.75, wgt) %>% slt(-sum.wgt), gdf %>% fnth(0.75, wgt, keep.w = FALSE))
      expect_visible(gdf %>% fmode)
      expect_visible(gdf %>% fmode(wgt))
      expect_equal(gdf %>% fmode(wgt) %>% slt(-sum.wgt), gdf %>% fmode(wgt, keep.w = FALSE))
      expect_visible(gdf %>% fsum)
      expect_visible(gdf %>% fsum(wgt))
      expect_equal(gdf %>% fsum(wgt) %>% slt(-sum.wgt), gdf %>% fsum(wgt, keep.w = FALSE))
      expect_visible(gdf %>% fprod)
      expect_visible(gdf %>% fprod(wgt))
      expect_equal(gdf %>% fprod(wgt) %>% slt(-prod.wgt), gdf %>% fprod(wgt, keep.w = FALSE))
      expect_visible(gdf %>% fsd)
      expect_visible(gdf %>% fsd(wgt))
      expect_equal(gdf %>% fsd(wgt) %>% slt(-sum.wgt), gdf %>% fsd(wgt, keep.w = FALSE))
      expect_visible(gdf %>% fvar)
      expect_visible(gdf %>% fvar(wgt))
      expect_equal(gdf %>% fvar(wgt) %>% slt(-sum.wgt), gdf %>% fvar(wgt, keep.w = FALSE))
      expect_visible(gdf %>% fmin)
      expect_visible(gdf %>% fmax)
      expect_visible(gdf %>% ffirst)
      expect_visible(gdf %>% flast)
      expect_visible(gdf %>% fnobs)
      expect_visible(gdf %>% fndistinct)
      expect_visible(gdf %>% collapg)
      expect_visible(gdf %>% varying)
      expect_visible(gdf %>% varying(any_group = FALSE))
      expect_visible(gdf %>% fmean(w = wgt)) # good?
      expect_equal(gdf %>% collapg(w = wgt) %>% slt(-wgt), gdf %>% collapg(w = wgt, keep.w = FALSE))
      expect_visible(gdf %>% fscale)
      expect_visible(gdf %>% fscale(wgt))
      expect_equal(gdf %>% fscale(wgt) %>% slt(-wgt), gdf %>% fscale(wgt, keep.w = FALSE))
      expect_visible(gdf %>% STD)
      expect_visible(gdf %>% STD(wgt))
      expect_equal(gdf %>% STD(wgt) %>% slt(-wgt), gdf %>% STD(wgt, keep.w = FALSE))
      expect_equal(gdf %>% fscale, gdf %>% STD(stub = FALSE))
      expect_visible(gdf %>% fbetween)
      expect_visible(gdf %>% fbetween(wgt))
      expect_equal(gdf %>% fbetween(wgt) %>% slt(-wgt), gdf %>% fbetween(wgt, keep.w = FALSE))
      expect_visible(gdf %>% B)
      expect_visible(gdf %>% B(wgt))
      expect_equal(gdf %>% B(wgt) %>% slt(-wgt), gdf %>% B(wgt, keep.w = FALSE))
      expect_equal(gdf %>% fbetween, gdf %>% B(stub = FALSE))
      expect_visible(gdf %>% fwithin)
      expect_visible(gdf %>% fwithin(wgt))
      expect_equal(gdf %>% fwithin(wgt) %>% slt(-wgt), gdf %>% fwithin(wgt, keep.w = FALSE))
      expect_visible(gdf %>% W)
      expect_visible(gdf %>% W(wgt))
      expect_equal(gdf %>% W(wgt) %>% slt(-wgt), gdf %>% W(wgt, keep.w = FALSE))
      expect_equal(gdf %>% fwithin, gdf %>% W(stub = FALSE))
      expect_visible(gdf %>% fcumsum)
      expect_visible(gdf %>% flag)
      expect_visible(gdf %>% L)
      expect_visible(gdf %>% F)
      expect_true(all_obj_equal(gdf %>% flag, gdf %>% L(stubs = FALSE), gdf %>% F(-1, stubs = FALSE)))
      expect_true(all_obj_equal(gdf %>% flag(-3:3), gdf %>% L(-3:3), gdf %>% F(3:-3)))
      expect_visible(gdf %>% fdiff)
      expect_visible(gdf %>% D)
      expect_true(all_obj_equal(gdf %>% fdiff, gdf %>% D(stubs = FALSE)))
      expect_equal(gdf %>% fdiff(-2:2, 1:2), gdf %>% D(-2:2, 1:2))
      expect_visible(gdf %>% fdiff(rho = 0.95))
      expect_visible(gdf %>% fdiff(-2:2, 1:2, rho = 0.95))
      expect_visible(gdf %>% fdiff(log = TRUE))
      expect_visible(gdf %>% fdiff(-2:2, 1:2, log = TRUE))
      expect_visible(gdf %>% fdiff(log = TRUE, rho = 0.95))
      expect_visible(gdf %>% fdiff(-2:2, 1:2, log = TRUE, rho = 0.95))
      expect_visible(gdf %>% fgrowth)
      expect_visible(gdf %>% G)
      expect_true(all_obj_equal(gdf %>% fgrowth, gdf %>% G(stubs = FALSE)))
      expect_equal(gdf %>% fgrowth(-2:2, 1:2), gdf %>% G(-2:2, 1:2))
      expect_visible(gdf %>% fgrowth(scale = 1))
      expect_visible(gdf %>% fgrowth(-2:2, 1:2, scale = 1))
      expect_visible(gdf %>% fgrowth(logdiff = TRUE))
      expect_visible(gdf %>% fgrowth(-2:2, 1:2, logdiff = TRUE))
      expect_visible(gdf %>% fgrowth(logdiff = TRUE, scale = 1))
      expect_visible(gdf %>% fgrowth(-2:2, 1:2, logdiff = TRUE, scale = 1))
      expect_equal(BY(gby(iris,Species), sum), BY(nv(gby(iris,Species)), sum))
    }
  }
})

# Also better not run on CRAN...
test_that("0-length vectors give expected output", {
  funs <- .c(fsum, fprod, fmean, fmedian, fmin, fmax, fnth, fcumsum, fbetween, fwithin, fscale)
  for(i in funs) {
    FUN <- match.fun(i)
    if(i %!in% .c(fsum, fmin, fmax, fcumsum, fprod, fmean, fmedian, fnth)) {
      expect_true(all_identical(FUN(numeric(0)), FUN(integer(0)), numeric(0)))
    } else {
      expect_identical(FUN(numeric(0)), numeric(0))
      if(i %in% .c(fmean, fprod, fnth, fmedian)) expect_identical(FUN(integer(0)), NA_real_)
      else expect_identical(FUN(integer(0)), integer(0))
    }
  }
  funs <- .c(fmode, ffirst, flast)
  for(i in funs) {
    FUN <- match.fun(i)
    expect_identical(FUN(numeric(0)), numeric(0))
    expect_identical(FUN(integer(0)), integer(0))
    expect_identical(FUN(character(0)), character(0))
    expect_identical(FUN(logical(0)), logical(0))
    expect_identical(FUN(factor(0)), factor(0))
  }
  funs <- .c(fvar, fsd)
  for(i in funs) {
    FUN <- match.fun(i)
    expect_identical(FUN(numeric(0)), NA_real_)
    expect_identical(FUN(integer(0)), NA_real_)
  }
  funs <- .c(fnobs, fndistinct)
  for(i in funs) {
    FUN <- match.fun(i)
    expect_identical(FUN(numeric(0)), 0L)
    expect_identical(FUN(integer(0)), 0L)
  }
  funs <- .c(flag, fdiff, fgrowth)
  for(i in funs) {
    FUN <- match.fun(i)
    expect_error(FUN(numeric(0)))
    expect_error(FUN(integer(0)))
  }
  funs <- .c(groupid, seqid)
  for(i in funs) {
    FUN <- match.fun(i)
    expect_identical(FUN(numeric(0)), integer(0))
    expect_identical(FUN(integer(0)), integer(0))
  }
  expect_identical(varying(numeric(0)), FALSE)
  expect_identical(TRA(numeric(0), 1), numeric(0))
})

}

X <- matrix(rnorm(1000), ncol = 10)
g <- qG(sample.int(10, 100, TRUE))
gf <- as_factor_qG(g)
funs <- grep("hd|log", c(.FAST_FUN, .OPERATOR_FUN), ignore.case = TRUE, invert = TRUE, value = TRUE)

test_that("functions work on plain matrices", {
  F <- getNamespace("collapse")$F
  for(i in funs) {
    expect_visible(match.fun(i)(X))
    expect_visible(match.fun(i)(X, g = g))
    expect_visible(match.fun(i)(X, g = gf))
    expect_visible(match.fun(i)(X, g = g, use.g.names = FALSE))
    expect_visible(match.fun(i)(X, g = gf, use.g.names = FALSE))
  }
})

Xl <- mctl(X)

test_that("functions work on plain lists", {
  F <- getNamespace("collapse")$F
  for(i in funs) {
    expect_visible(match.fun(i)(Xl))
    expect_visible(match.fun(i)(Xl, g = g, by = g))
    expect_visible(match.fun(i)(Xl, g = gf, by = gf))
    expect_visible(match.fun(i)(X, g = g, by = g, use.g.names = FALSE))
    expect_visible(match.fun(i)(X, g = gf, by = gf, use.g.names = FALSE))
  }
})

test_that("time series functions work inside lm", {
  expect_equal(unname(coef(lm(mpg ~ L(cyl, 0:2), mtcars))), unname(coef(lm(mpg ~ cyl + L(cyl, 1) + L(cyl, 2), mtcars))))
  expect_equal(unname(coef(lm(mpg ~ F(cyl, 0:2), mtcars))), unname(coef(lm(mpg ~ cyl + F(cyl, 1) + F(cyl, 2), mtcars))))
  expect_equal(unname(coef(lm(mpg ~ D(cyl, 0:2), mtcars))), unname(coef(lm(mpg ~ cyl + D(cyl, 1) + D(cyl, 2), mtcars))))
  expect_equal(unname(coef(lm(mpg ~ G(cyl, 0:2), mtcars))), unname(coef(lm(mpg ~ cyl + G(cyl, 1) + G(cyl, 2), mtcars))))

  expect_equal(unname(coef(lm(mpg ~ L(L(cyl, 0:2)), mtcars))), unname(coef(lm(mpg ~ L(cyl) + L(cyl, 2) + L(cyl, 3), mtcars))))
  expect_equal(unname(coef(lm(mpg ~ L(F(cyl, 0:2)), mtcars))), unname(coef(lm(mpg ~ L(cyl) + cyl + F(cyl, 1), mtcars))))
  expect_equal(unname(coef(lm(mpg ~ L(D(cyl, 0:2)), mtcars))), unname(coef(lm(mpg ~ L(cyl) + L(D(cyl)) + L(D(cyl, 2)), mtcars))))
  expect_equal(unname(coef(lm(mpg ~ L(G(cyl, 0:2)), mtcars))), unname(coef(lm(mpg ~ L(cyl) + L(G(cyl)) + L(G(cyl, 2)), mtcars))))

})

test_that("functions using welfords method properly deal with zero weights", {
  for(g in list(NULL, rep(1L, 3))) {
    expect_equal(unattrib(fvar(x = c(2, 1, 0), g = g, w = c(1, 1, 0), na.rm = TRUE)), 0.5)
    expect_equal(unattrib(fvar(x = c(2, 1, 3), g = g, w = c(0, 1, 1), na.rm = FALSE)), 2)
    expect_equal(unattrib(fsd(x = c(2, 1, 0), g = g, w = c(1, 1, 0), na.rm = TRUE)), sqrt(0.5))
    expect_equal(unattrib(fsd(x = c(2, 1, 3), g = g, w = c(0, 1, 1), na.rm = FALSE)), sqrt(2))
    expect_equal(unattrib(fscale(x = c(2, 1, 0), g = g, w = c(1, 1, 0), na.rm = TRUE)), (c(2, 1, 0)-1.5)/sqrt(0.5))
    expect_equal(unattrib(fscale(x = c(2, 1, 3), g = g, w = c(0, 1, 1), na.rm = FALSE)), (c(2, 1, 3)-2)/sqrt(2))
    expect_equal(unattrib(qsu(x = c(2, 1, 0), g = g, w = c(1, 1, 0))), c(2, 1.5, sqrt(0.5), 1, 2))
    expect_equal(unattrib(qsu(x = c(2, 1, 3), g = g, w = c(0, 1, 1))), c(2, 2, sqrt(2), 1, 3))
    expect_equal(unattrib(qsu(x = c(2, 1, 0), g = g, w = c(1, 1, 0), higher = TRUE))[1:5], c(2, 1.5, sqrt(0.5), 1, 2))
    expect_equal(unattrib(qsu(x = c(2, 1, 3), g = g, w = c(0, 1, 1), higher = TRUE))[1:5], c(2, 2, sqrt(2), 1, 3))

  }
})


test_that("singleton groups are handled properly by all statistical functions", {
  w <- rep(1, fnrow(wlddev))
  # Ordered
  g <- GRP(seq_row(wlddev), return.groups = FALSE)
  expect_equal(fmode(wlddev, g), wlddev)
  expect_equal(fmode(wlddev, g, w), wlddev)
  expect_equal(ffirst(wlddev, g), wlddev)
  expect_equal(flast(wlddev, g), wlddev)
  expect_equal(dapply(fndistinct(wlddev, g), unattrib), dapply(wlddev, function(x) as.integer(!is.na(x))))
  expect_equal(fmode(wlddev, g, na.rm = FALSE), wlddev)
  expect_equal(fmode(wlddev, g, w, na.rm = FALSE), wlddev)
  expect_equal(ffirst(wlddev, g, na.rm = FALSE), wlddev)
  expect_equal(flast(wlddev, g, na.rm = FALSE), wlddev)
  expect_equal(dapply(fndistinct(wlddev, g, na.rm = FALSE), unattrib), dapply(wlddev, function(x) rep(1L, length(x))))
  for(FUN in list(fmean, fmedian, fnth, fsum, fprod, fmin, fmax, fbetween, fcumsum)) {
    # print(FUN)
    expect_equal(FUN(nv(wlddev), g = g), nv(wlddev))
    expect_equal(FUN(nv(wlddev), g = g, na.rm = FALSE), nv(wlddev))
    expect_equal(FUN(nv(wlddev), g = g, w = w), nv(wlddev))
    expect_equal(FUN(nv(wlddev), g = g, w = w, na.rm = FALSE), nv(wlddev))
  }
  for(FUN in list(fvar, fsd, fscale, flag, fdiff, fgrowth)) {
    expect_true(all(dapply(FUN(nv(wlddev), g = g), allNA)))
    expect_true(all(dapply(FUN(nv(wlddev), g = g, na.rm = FALSE), allNA)))
    expect_true(all(dapply(FUN(nv(wlddev), g = g, w = w, n = -1), allNA)))
    expect_true(all(dapply(FUN(nv(wlddev), g = g, w = w, n = -1, na.rm = FALSE), allNA)))
  }
  # Unordered
  o <- radixorder(rnorm(fnrow(wlddev)))
  g <- GRP(o, return.groups = FALSE)
  wlduo <- setRownames(ss(wlddev, radixorder(o)))
  expect_equal(fmode(wlddev, g), wlduo)
  expect_equal(fmode(wlddev, g, w), wlduo)
  expect_equal(ffirst(wlddev, g), wlduo)
  expect_equal(flast(wlddev, g), wlduo)
  expect_equal(dapply(fndistinct(wlddev, g), unattrib), dapply(wlduo, function(x) as.integer(!is.na(x))))
  expect_equal(fmode(wlddev, g, na.rm = FALSE), wlduo)
  expect_equal(fmode(wlddev, g, w, na.rm = FALSE), wlduo)
  expect_equal(ffirst(wlddev, g, na.rm = FALSE), wlduo)
  expect_equal(flast(wlddev, g, na.rm = FALSE), wlduo)
  expect_equal(dapply(fndistinct(wlddev, g, na.rm = FALSE), unattrib), dapply(wlduo, function(x) rep(1L, length(x))))
  for(FUN in list(fmean, fmedian, fnth, fsum, fprod, fmin, fmax)) {
    # print(FUN)
    expect_equal(FUN(nv(wlddev), g = g), nv(wlduo))
    expect_equal(FUN(nv(wlddev), g = g, na.rm = FALSE), nv(wlduo))
    expect_equal(FUN(nv(wlddev), g = g, w = w), nv(wlduo))
    expect_equal(FUN(nv(wlddev), g = g, w = w, na.rm = FALSE), nv(wlduo))
  }
  for(FUN in list(fbetween, fcumsum)) {
    expect_equal(FUN(nv(wlddev), g), nv(wlddev))
    expect_equal(FUN(nv(wlddev), g, na.rm = FALSE), nv(wlddev))
    expect_equal(FUN(nv(wlddev), g, w), nv(wlddev))
    expect_equal(FUN(nv(wlddev), g, w, na.rm = FALSE), nv(wlddev))
  }
  for(FUN in list(fvar, fsd, fscale, flag, fdiff, fgrowth)) {
    expect_true(all(dapply(FUN(nv(wlddev), g = g), allNA)))
    expect_true(all(dapply(FUN(nv(wlddev), g = g, na.rm = FALSE), allNA)))
    expect_true(all(dapply(FUN(nv(wlddev), g = g, w = w, n = -1), allNA)))
    expect_true(all(dapply(FUN(nv(wlddev), g = g, w = w, n = -1, na.rm = FALSE), allNA)))
  }
})


test_that("functions work for data frames with zero rows", {
  mtc0 <- qDF(mtcars)[NULL, ]
  expect_equal(mtc0, funique(mtc0))
  expect_equal(mtc0, funique(mtc0, sort = TRUE))
  expect_equal(mtc0, roworderv(mtc0))
  expect_visible(colorder(mtc0, mpg, hp))
  expect_visible(GRP(mtc0))
  expect_visible(fgroup_by(mtc0, cyl, vs, am))
  expect_visible(GRP(mtc0, sort = FALSE))
  expect_visible(fgroup_by(mtc0, cyl, vs, am, sort = FALSE))
  expect_visible(fduplicated(mtc0))
  expect_false(any_duplicated(mtc0))
  expect_visible(fselect(mtc0, hp, carb))
  expect_visible(get_vars(mtc0, 9:8))
})

test_that("issue with integer followed by NA #432", {
    for (f in setdiff(.FAST_STAT_FUN, c("fvar", "fsd", "fnobs", "fndistinct"))) {
      # if(!isTRUE(all.equal(match.fun(f)(c(10L, NA)), 10L))) print(f)
      expect_equal(match.fun(f)(c(10L, NA)), 10L)
      expect_equal(match.fun(f)(c(NA, 10L)), 10L)
      expect_equal(match.fun(f)(c(10, NA)), 10)
      expect_equal(match.fun(f)(c(NA, 10)), 10)
      expect_equal(match.fun(f)(c(10L, NA), g = rep(1L, 2), use.g.names = FALSE), 10L)
      expect_equal(match.fun(f)(c(NA, 10L), g = rep(1L, 2), use.g.names = FALSE), 10L)
      expect_equal(match.fun(f)(c(10, NA), g = rep(1L, 2), use.g.names = FALSE), 10)
      expect_equal(match.fun(f)(c(NA, 10), g = rep(1L, 2), use.g.names = FALSE), 10)
      # na.rm = FALSE
      if(f %!in% c("fmode", "ffirst")) expect_equal(match.fun(f)(c(10L, NA), na.rm = FALSE), NA_integer_)
      if(f != "flast") expect_equal(match.fun(f)(c(NA, 10L), na.rm = FALSE), NA_integer_)
      if(f %!in% c("fmode", "ffirst")) expect_equal(match.fun(f)(c(10, NA), na.rm = FALSE), NA_real_)
      if(f != "flast") expect_equal(match.fun(f)(c(NA, 10), na.rm = FALSE), NA_real_)
      # Some functions are optimized and don't check here
      # expect_equal(match.fun(f)(c(10L, NA), g = rep(1L, 2), na.rm = FALSE, use.g.names = FALSE), NA_integer_)
      # expect_equal(match.fun(f)(c(NA, 10L), g = rep(1L, 2), na.rm = FALSE, use.g.names = FALSE), NA_integer_)
      if(f %!in% c("fmode", "ffirst")) expect_equal(match.fun(f)(c(10, NA), g = rep(1L, 2), na.rm = FALSE, use.g.names = FALSE), NA_real_)
      if(f != "flast") expect_equal(match.fun(f)(c(NA, 10), g = rep(1L, 2), na.rm = FALSE, use.g.names = FALSE), NA_real_)
    }
  skip_if_not(Sys.getenv("OMP") == "TRUE")
  for (f in c("fsum", "fmean", "fmode", "fnth", "fmedian")) {
    expect_equal(match.fun(f)(c(10L, rep(NA_integer_, 1e5)), nthreads = 2L), 10L)
    expect_equal(match.fun(f)(c(rep(NA_integer_, 1e5), 10L), nthreads = 2L), 10L)
    expect_equal(match.fun(f)(c(10, rep(NA_real_, 1e5)), nthreads = 2L), 10)
    expect_equal(match.fun(f)(c(rep(NA_real_, 1e5), 10), nthreads = 2L), 10)
    expect_equal(match.fun(f)(c(10L, rep(NA_integer_, 1e5)), g = rep(1L, 1e5+1), nthreads = 2L, use.g.names = FALSE), 10L)
    expect_equal(match.fun(f)(c(rep(NA_integer_, 1e5), 10L), g = rep(1L, 1e5+1), nthreads = 2L, use.g.names = FALSE), 10L)
    expect_equal(match.fun(f)(c(10, rep(NA_real_, 1e5)), g = rep(1L, 1e5+1), nthreads = 2L, use.g.names = FALSE), 10)
    expect_equal(match.fun(f)(c(rep(NA_real_, 1e5), 10), g = rep(1L, 1e5+1), nthreads = 2L, use.g.names = FALSE), 10)
    # na.rm = FALSE
    expect_equal(match.fun(f)(c(10L, rep(NA_integer_, 1e5)), na.rm = FALSE, nthreads = 2L), NA_integer_)
    expect_equal(match.fun(f)(c(rep(NA_integer_, 1e5), 10L), na.rm = FALSE, nthreads = 2L), NA_integer_)
    expect_equal(match.fun(f)(c(10, rep(NA_real_, 1e5)), na.rm = FALSE, nthreads = 2L), NA_real_)
    expect_equal(match.fun(f)(c(rep(NA_real_, 1e5), 10), na.rm = FALSE, nthreads = 2L), NA_real_)
    # Some functions are optimized and don't check here
    # expect_equal(match.fun(f)(c(10L, rep(NA_integer_, 1e5)), g = rep(1L, 1e5+1), na.rm = FALSE, nthreads = 2L, use.g.names = FALSE), NA_integer_)
    # expect_equal(match.fun(f)(c(rep(NA_integer_, 1e5), 10L), g = rep(1L, 1e5+1), na.rm = FALSE, nthreads = 2L, use.g.names = FALSE), NA_integer_)
    expect_equal(match.fun(f)(c(10, rep(NA_real_, 1e5)), g = rep(1L, 1e5+1), na.rm = FALSE, nthreads = 2L, use.g.names = FALSE), NA_real_)
    expect_equal(match.fun(f)(c(rep(NA_real_, 1e5), 10), g = rep(1L, 1e5+1), na.rm = FALSE, nthreads = 2L, use.g.names = FALSE), NA_real_)
  }
})

test_that("fmedian ties handled properly with weights", {
  x <- c(1, 2, 3, 4)
  w <- c(2.5, 2.4, 3.8, 1.1)
  expect_equal(c(fmedian(x, w = w, ties = "mean"), fmedian(x, w = w, ties = "min"), fmedian(x, w = w, ties = "max")),
               c(2.5, 2, 3))
  w <- c(2.5, 2.4, 3.7, 1.2)
  expect_equal(c(fmedian(x, w = w, ties = "mean"), fmedian(x, w = w, ties = "min"), fmedian(x, w = w, ties = "max")),
               c(2.5, 2, 3))
})

options(warn = 1)
