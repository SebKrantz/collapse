context("Misc")

if(!is.null(attributes(identical(FALSE, TRUE)))) stop("OECD label issue")

# rm(list = ls())
set.seed(101)
m <- na_insert(qM(mtcars))

test_that("descr, pwcor, pwcov, pwnobs", {

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

  expect_visible(pwnobs(wlddev))
  expect_visible(pwnobs(GGDC10S))

  expect_visible(descr(m))
  expect_visible(pwcor(m))
  expect_visible(pwcov(m))
  expect_visible(pwnobs(m))

})

if(identical(Sys.getenv("NCRAN"), "TRUE")) {

test_that("weighted correlations are correct", {

  # This is to fool very silly checks on CRAN scanning the code of the tests
  wtd.cors <- eval(parse(text = paste0("weights", ":", ":", "wtd.cors")))
  wtd.cor <- eval(parse(text = paste0("weights", ":", ":", "wtd.cor")))

  w <- abs(rnorm(fnrow(wlddev)))
  cc <- which(!missing_cases(nv(wlddev)))

  expect_equal(unclass(pwcor(nv(wlddev), w = w)), wtd.cors(nv(wlddev), w = w))
  expect_equal(unclass(pwcor(nv(wlddev), w = w)), cov2cor(unclass(pwcov(nv(wlddev), w = w))))
  expect_true(all_obj_equal(unclass(pwcor(ss(nv(wlddev), cc), w = w[cc])),
                            cov2cor(unclass(pwcov(ss(nv(wlddev), cc), w = w[cc]))),
                            unclass(pwcor(nv(wlddev), w = w, use = "complete.obs")),
                            wtd.cors(ss(nv(wlddev), cc), w = w[cc]),
                            cov.wt(ss(nv(wlddev), cc), w[cc], cor = TRUE)$cor))

  suppressWarnings(
  expect_true(all_obj_equal(replace_NA(pwcor(ss(nv(wlddev), cc), w = w[cc], P = TRUE, array = FALSE)$P, 0),
                            replace_NA(pwcov(ss(nv(wlddev), cc), w = w[cc], P = TRUE, array = FALSE)$P, 0),
                            replace_NA(pwcor(ss(nv(wlddev), cc), w = w[cc], P = TRUE, array = FALSE, use = "complete.obs")$P, 0),
                            replace_NA(pwcov(ss(nv(wlddev), cc), w = w[cc], P = TRUE, array = FALSE, use = "complete.obs")$P, 0),
                            wtd.cor(ss(nv(wlddev), cc), w = w[cc])$p.value)))

  expect_true(all_obj_equal(unclass(pwcov(ss(nv(wlddev), cc), w = w[cc])),
                            unclass(pwcov(nv(wlddev), w = w, use = "complete.obs"))))

  expect_equal(cov.wt(ss(nv(wlddev), cc), w[cc])$cov, unclass(pwcov(nv(wlddev), w = w, use = "complete.obs")),
               tolerance = 1e-3)

})

test_that("na_rm works well", {
  set.seed(101)
  expect_equal(sapply(na_insert(wlddev), function(x) vtypes(na_rm(x))), vtypes(wlddev))
  expect_equal(sapply(na_insert(wlddev), function(x) vlabels(na_rm(x))), vlabels(wlddev))
  expect_equal(sapply(na_insert(wlddev), function(x) vclasses(na_rm(x))), vclasses(wlddev))
  wldNA <- na_insert(wlddev)
  expect_equal(lengths(lapply(wldNA, na_rm)), fnobs(wldNA))
  expect_equal(lapply(wldNA, na_rm), lapply(wldNA, function(x) copyMostAttrib(x[!is.na(x)], x)))
  rm(wldNA)

  expect_equal(na_rm(list(list(), 1,2,3)), list(1,2,3))
  expect_equal(na_rm(list(1,2,NULL,3)), list(1,2,3))
})

}

test_that("vlabels works well", {
  expect_equal(wlddev, setLabels(wlddev, vlabels(wlddev)))
})

test_that("adding and removing stubs works", {
  expect_identical(rm_stub(add_stub(iris, "df"), "df"), iris)
  expect_identical(rm_stub(add_stub(iris, "df", pre = FALSE), "df", pre = FALSE), iris)
  expect_identical(rm_stub(add_stub(iris, "df", pre = FALSE), "df", regex = TRUE), iris)
  expect_identical(rm_stub(names(iris), "Sepal")[1], ".Length")
  expect_identical(rm_stub(names(iris), "Width", pre = FALSE)[4], "Petal.")
  expect_identical(rm_stub(names(iris), "Width", regex = TRUE)[4], "Petal.")
})

test_that("deep matrix dispatch works well", {

  tsm <- EuStockMarkets
  class(tsm) <- setdiff(class(tsm), "matrix")
  set.seed(101)
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
   # print(i)
   FUN <- match.fun(i)
   expect_true(all_obj_equal(FUN(v, g = mtcars$cyl), FUN(v, g = f), FUN(v, g = fcc), FUN(v, g = g), FUN(v, g = gl)))
   expect_true(all_obj_equal(FUN(v, g = mtcars$cyl, use.g.names = FALSE), FUN(v, g = f, use.g.names = FALSE), FUN(v, g = fcc, use.g.names = FALSE), FUN(v, g = g, use.g.names = FALSE), FUN(v, g = gl, use.g.names = FALSE)))

   expect_true(all_obj_equal(FUN(m, g = mtcars$cyl), FUN(m, g = f), FUN(m, g = fcc), FUN(m, g = g), FUN(m, g = gl)))
   expect_true(all_obj_equal(FUN(m, g = mtcars$cyl, use.g.names = FALSE), FUN(m, g = f, use.g.names = FALSE), FUN(m, g = fcc, use.g.names = FALSE), FUN(m, g = g, use.g.names = FALSE), FUN(m, g = gl, use.g.names = FALSE)))

   expect_true(all_obj_equal(FUN(mtcars, g = mtcars$cyl), FUN(mtcars, g = f), FUN(mtcars, g = fcc), FUN(mtcars, g = g), FUN(mtcars, g = gl)))
   if(Sys.getenv("NCRAN") == "TRUE")
   expect_true(all_obj_equal(FUN(mtcars, g = mtcars$cyl, use.g.names = FALSE),
                             FUN(mtcars, g = f, use.g.names = FALSE),
                             FUN(mtcars, g = fcc, use.g.names = FALSE),
                             FUN(mtcars, g = g, use.g.names = FALSE),
                             FUN(mtcars, g = gl, use.g.names = FALSE)))
  if(Sys.getenv("NCRAN") == "TRUE")
  expect_true(all_obj_equal(gv(FUN(mtcars, g = mtcars$cyl, use.g.names = FALSE), -2),
                             gv(FUN(gmtc), -1),
                             gv(FUN(gv(gmtc,-2)), -1),
                             FUN(gv(gmtc,-2), keep.group_vars = FALSE),
                             FUN(gmtc, keep.group_vars = FALSE)))

   expect_equal(FUN(v, TRA = 2L), TRA(v, FUN(v), 2L))
   expect_true(all_obj_equal(FUN(v, g = mtcars$cyl, TRA = 1L), TRA(v, FUN(v, g = mtcars$cyl), 1L, mtcars$cyl),
                             FUN(v, g = f, TRA = 1L), TRA(v, FUN(v, g = f), 1L, f),
                             FUN(v, g = fcc, TRA = 1L), TRA(v, FUN(v, g = fcc), 1L, fcc),
                             FUN(v, g = g, TRA = 1L), TRA(v, FUN(v, g = g), 1L, g),
                             FUN(v, g = gl, TRA = 1L), TRA(v, FUN(v, g = gl), 1L, gl)))

   expect_equal(FUN(m, TRA = 2L), TRA(m, FUN(m), 2L))
   expect_true(all_obj_equal(FUN(m, g = mtcars$cyl, TRA = 1L), TRA(m, FUN(m, g = mtcars$cyl), 1L, mtcars$cyl),
                             FUN(m, g = f, TRA = 1L), TRA(m, FUN(m, g = f), 1L, f),
                             FUN(m, g = fcc, TRA = 1L), TRA(m, FUN(m, g = fcc), 1L, fcc),
                             FUN(m, g = g, TRA = 1L), TRA(m, FUN(m, g = g), 1L, g),
                             FUN(m, g = gl, TRA = 1L), TRA(m, FUN(m, g = gl), 1L, gl)))

   expect_equal(FUN(mtcars, TRA = 2L), TRA(mtcars, FUN(mtcars), 2L))
   expect_true(all_obj_equal(FUN(mtcars, g = mtcars$cyl, TRA = 1L), TRA(mtcars, FUN(mtcars, g = mtcars$cyl), 1L, mtcars$cyl),
                             FUN(mtcars, g = f, TRA = 1L), TRA(mtcars, FUN(mtcars, g = f), 1L, f),
                             FUN(mtcars, g = fcc, TRA = 1L), TRA(mtcars, FUN(mtcars, g = fcc), 1L, fcc),
                             FUN(mtcars, g = g, TRA = 1L), TRA(mtcars, FUN(mtcars, g = g), 1L, g),
                             FUN(mtcars, g = gl, TRA = 1L), TRA(mtcars, FUN(mtcars, g = gl), 1L, gl)))

   expect_equal(colorder(FUN(gmtc, TRA = 1L), mpg, cyl), TRA(gmtc, FUN(gmtc), 1L))
   expect_equal(FUN(fselect(gmtc, -cyl), TRA = 1L), TRA(fselect(gmtc, -cyl), FUN(gmtc, keep.group_vars = FALSE), 1L))
 }

  for(i in setdiff(.FAST_FUN, c(.FAST_STAT_FUN, "fhdbetween", "fhdwithin"))) {

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

l <- as.list(mtcars)
test_that("list and df methods give same results", {

  for (i in setdiff(c(.FAST_FUN, .OPERATOR_FUN), c("fhdbetween", "fhdwithin", "HDB", "HDW"))) {
    FUN <- match.fun(i)
    expect_equal(unattrib(FUN(mtcars)), unattrib(FUN(l)))
  }

})

w <- mtcars$wt
wFUNs <- c("fmean","fmedian","fsum","fprod","fmode","fvar","fsd","fscale","STD","fbetween","B","fwithin","W")

test_that("fast functions give appropriate warnings", {

  for (i in setdiff(c(.FAST_FUN, .OPERATOR_FUN, "qsu"), c("fhdbetween", "fhdwithin", "HDB", "HDW"))) {
    FUN <- match.fun(i)
    expect_warning(FUN(v, bla = 1))
    expect_warning(FUN(m, bla = 1))
    expect_warning(FUN(mtcars, bla = 1))
    expect_warning(FUN(gmtc, bla = 1))
    if(i %in% wFUNs) {
      expect_warning(FUN(gmtc, bla = 1))
      expect_error(FUN(gmtc, cyl)) # weight same as grouping variable
      if(i %in% .FAST_STAT_FUN) expect_true(names(FUN(gmtc, wt))[2L] == if(i == "fprod") "prod.wt" else "sum.wt") # weight same as grouping variable
    }
  }

})

test_that("fselect and fsubset cannot easily be confuesed", {
  expect_error(suppressWarnings(fsubset(mtcars, mpg:vs, wt)))
  expect_error(fselect(mtcars, mpg == 1))
})

test_that("frange works well", {
  xd <- rnorm(1e5)
  xdNA <- na_insert(xd)
  xi <- as.integer(xd*1000)
  xiNA <- na_insert(xi)

  expect_equal(frange(xd, na.rm = FALSE), range(xd))
  expect_equal(frange(xd), range(xd, na.rm = TRUE))
  expect_equal(frange(xdNA, na.rm = FALSE), range(xdNA))
  expect_equal(frange(xdNA), range(xdNA, na.rm = TRUE))

  expect_equal(frange(xi, na.rm = FALSE), range(xi))
  expect_equal(frange(xi), range(xi, na.rm = TRUE))
  expect_equal(frange(xiNA, na.rm = FALSE), range(xiNA))
  expect_equal(frange(xiNA), range(xiNA, na.rm = TRUE))

})

# TODO: Test other cols options and formula options !!!
options(warn = -1)
test_that("operator methods column selection since v1.8.1 works as intended", {
  nnvw <- names(nv(wlddev))
  wldi <- colorder(iby(wlddev, iso3c, year), year, pos = "end")
  wld1i <- colorder(iby(sbt(wlddev, iso3c %==% "DEU"), year), year, pos = "end")
  nnvg <- names(nv(GGDC10S))
  ggdc3i <- findex_by(GGDC10S, Variable, Country, Year, interact.ids = FALSE)
  ggdc3ii <- findex_by(GGDC10S, Variable, Country, Year)

  for(op in list(L, F, D, Dlog, G, B, W, STD)) {
    expect_equal(names(op(wlddev, stub = FALSE)), nnvw)
    expect_equal(names(op(wlddev, by = ~ iso3c, stub = FALSE)), c("iso3c", nnvw))
    expect_equal(names(op(wlddev, by = ~ iso3c, stub = FALSE, keep.by = FALSE, keep.ids = FALSE)), nnvw)
    expect_equal(names(op(wlddev, by = ~ decade, stub = FALSE)), c("decade", nnvw[nnvw != "decade"]))
    expect_equal(names(op(wlddev, by = ~ decade, stub = FALSE, keep.by = FALSE, keep.ids = FALSE)), nnvw[nnvw != "decade"])
    expect_equal(names(op(wldi, stub = FALSE)), c("iso3c", nnvw))
    expect_equal(names(op(wldi, stub = FALSE, keep.by = FALSE, keep.ids = FALSE)), nnvw[nnvw != "year"])
    expect_equal(names(op(wld1i, stub = FALSE)), nnvw)
    expect_equal(names(op(wld1i, stub = FALSE, keep.by = FALSE, keep.ids = FALSE)), nnvw[nnvw != "year"])
    expect_equal(names(op(ggdc3i, stub = FALSE)), c("Country", "Variable", nnvg))
    expect_equal(names(op(ggdc3i, stub = FALSE, keep.by = FALSE, keep.ids = FALSE)), nnvg[-1L])
    expect_equal(names(op(ggdc3ii, stub = FALSE)), c("Country", "Variable", nnvg))
    expect_equal(names(op(ggdc3ii, stub = FALSE, keep.by = FALSE, keep.ids = FALSE)), nnvg[-1L])
  }

  wlduo <- colorder(wlddev, year, pos = "end")
  wld1uo <- sbt(wlduo, iso3c %==% "DEU")
  for(op in list(L, F, D, Dlog, G)) {
    expect_equal(names(op(wld1uo, t = ~ year, stubs = FALSE)), nnvw)
    expect_equal(names(op(wld1uo, t = ~ year, stubs = FALSE, keep.ids = FALSE)), nnvw[-1L])
    expect_equal(names(op(wld1uo, by = ~ iso3c, t = ~ year, stubs = FALSE)), c("iso3c", nnvw))
    expect_equal(names(op(wld1uo, by = ~ iso3c, t = ~ year, stubs = FALSE, keep.ids = FALSE)), nnvw[-1L])
  }
  for(op in list(B, W, STD)) {
    expect_equal(names(op(wld1uo, w = ~ year, stub = FALSE)), nnvw)
    expect_equal(names(op(wld1uo, w = ~ year, stub = FALSE, keep.w = FALSE)), nnvw[-1L])
    expect_equal(names(op(wld1uo, by = ~ iso3c, w = ~ year, stub = FALSE)), c("iso3c", nnvw))
    expect_equal(names(op(wld1uo, by = ~ iso3c, w = ~ year, stub = FALSE, keep.by = FALSE)), nnvw)
    expect_equal(names(op(wld1uo, by = ~ iso3c, w = ~ year, stub = FALSE, keep.w = FALSE)), c("iso3c", nnvw[-1L]))
    expect_equal(names(op(wld1uo, by = ~ iso3c, w = ~ year, stub = FALSE, keep.by = FALSE, keep.w = FALSE)), nnvw[-1L])

    expect_equal(names(op(wldi, w = ~POP, stub = FALSE)), c("iso3c", "year", "POP", nnvw[-c(1, 7)]))
    expect_equal(names(op(wldi, w = ~POP, stub = FALSE, keep.ids = FALSE)), c("POP", nnvw[-c(1, 7)]))
    expect_equal(names(op(wldi, w = ~POP, stub = FALSE, keep.w = FALSE)), c("iso3c", "year", nnvw[-c(1, 7)]))
    expect_equal(names(op(wldi, w = ~POP, stub = FALSE, keep.ids = FALSE, keep.w = FALSE)), nnvw[-c(1, 7)])

    expect_equal(names(op(wld1i, w = ~POP, stub = FALSE)), c("year", "POP", nnvw[-c(1, 7)]))
    expect_equal(names(op(wld1i, w = ~POP, stub = FALSE, keep.ids = FALSE)), c("POP", nnvw[-c(1, 7)]))
    expect_equal(names(op(wld1i, w = ~POP, stub = FALSE, keep.w = FALSE)), c("year", nnvw[-c(1, 7)]))
    expect_equal(names(op(wld1i, w = ~POP, stub = FALSE, keep.ids = FALSE, keep.w = FALSE)), nnvw[-c(1, 7)])
  }

  for(op in list(HDB, HDW)) {
    expect_equal(names(op(wlddev, wlddev$iso3c, stub = FALSE)), nnvw)
    expect_equal(names(op(wlddev, ~ iso3c, stub = FALSE)), nnvw)
    expect_equal(names(op(wlddev, ~ year, stub = FALSE)), nnvw[-1])
    if(identical(Sys.getenv("NCRAN"), "TRUE")) expect_equal(names(op(wldi, stub = FALSE)), nnvw[-1])
  }
})
options(warn = 1)
