context("Attribute Handling")



v <- wlddev$PCGDP
date <- wlddev$date
fac <- wlddev$region
g1 <- GRP(wlddev$country)
m <- qM(mtcars)
gmtc <- fgroup_by(mtcars, cyl, vs, am)
gm <- qM(gmtc, keep.attr = TRUE)
g2 <- GRP(mtcars, ~ cyl + vs + am)

# gDTmtc <- fgroup_by(qDT(mtcars), cyl, vs, am)
set.seed(101)
f1 <- sample.int(5, length(AirPassengers), replace = TRUE)
f2 <- sample.int(5, nrow(EuStockMarkets), replace = TRUE)

# numFUN <- setdiff(.FAST_STAT_FUN, c("fnth", "fmode", "ffirst", "flast", "fmin", "fmax"))
countFUN <- c("fnobs", "fndistinct")

F <- getNamespace("collapse")$F

test_that("statistical functions handle attributes properly", {

  for(i in setdiff(.FAST_STAT_FUN, "fnth")) {
    # print(i)
    FUN <- match.fun(i)
    if(i %!in% c("fvar", "fsd", countFUN)) {
      expect_identical(attributes(FUN(v)), attributes(v))
      expect_identical(attributes(FUN(date)), attributes(date))
    }
    if(i %!in% c("fsum", "fvar", "fsd", countFUN)) expect_identical(attributes(FUN(fac)), attributes(fac))
    if(i != "fmode") expect_true(is.null(attributes(FUN(AirPassengers))))
    expect_identical(attributes(FUN(EuStockMarkets)), list(names = colnames(EuStockMarkets)))
    expect_identical(attributes(FUN(EuStockMarkets, drop = FALSE)), list(dim = c(1L, 4L), dimnames = list(NULL, colnames(EuStockMarkets))))
    expect_identical(attributes(FUN(m)), list(names = colnames(m)))
    expect_identical(attributes(FUN(m, drop = FALSE)), list(dim = c(1L, 11L), dimnames = list(NULL, colnames(m))))
    expect_identical(attributes(FUN(gm)), list(names = colnames(m)))
    expect_identical(attributes(FUN(gm, drop = FALSE)), list(dim = c(1L, 11L), dimnames = list(NULL, colnames(m)), groups = attr(gm, "groups")))
    expect_identical(attributes(FUN(mtcars)), list(names = names(mtcars)))
    expect_identical(attributes(FUN(mtcars, drop = FALSE)), `[[<-`(attributes(mtcars), "row.names", 1L))
    expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"))), list(names = names(mtcars)))
    expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), drop = FALSE)), `[[<-`(`[[<-`(attributes(gmtc), "row.names", 1L), "class", "data.frame"))

    # Grouped
    expect_identical(attributes(FUN(v, g1, use.g.names = FALSE)), attributes(v))
    expect_identical(attributes(FUN(v, g1)), c(attributes(v), list(names = unattrib(GRPnames(g1)))))
    expect_identical(attributes(FUN(date, g1, use.g.names = FALSE)), if(i %!in% countFUN) attributes(date) else list(label = vlabels(date)))
    expect_identical(attributes(FUN(date, g1)), if(i %!in% countFUN) c(attributes(date), list(names = unattrib(GRPnames(g1)))) else list(label = vlabels(date), names = unattrib(GRPnames(g1))))
    expect_identical(attributes(FUN(fac, g1, use.g.names = FALSE)), if(i %!in% countFUN) attributes(fac) else list(label = vlabels(fac)))
    expect_identical(attributes(FUN(fac, g1)), if(i %!in% countFUN) c(attributes(fac), list(names = unattrib(GRPnames(g1)))) else list(label = vlabels(fac), names = unattrib(GRPnames(g1))))
    if(i != "fmode") expect_identical(attributes(FUN(AirPassengers, f1, use.g.names = FALSE)), NULL)
    if(i != "fmode") expect_identical(attributes(FUN(AirPassengers, f1)), list(names = as.character(1:5)))
    expect_identical(attributes(FUN(EuStockMarkets, f2, use.g.names = FALSE)), list(dim = c(5L, 4L), dimnames = dimnames(EuStockMarkets)))
    expect_identical(attributes(FUN(EuStockMarkets, f2)), list(dim = c(5L, 4L), dimnames = list(as.character(1:5), colnames(EuStockMarkets))))
    expect_identical(attributes(FUN(m, g2, use.g.names = FALSE)), list(dim = c(7L, 11L), dimnames = list(NULL, colnames(m))))
    expect_identical(attributes(FUN(m, g2)), list(dim = c(7L, 11L), dimnames = list(GRPnames(g2), colnames(m))))
    expect_identical(attributes(FUN(gm, attr(gm, "groups"), use.g.names = FALSE)), list(dim = c(7L, 11L), dimnames = list(NULL, colnames(m)), groups = attr(gm, "groups")))
    expect_identical(attributes(FUN(gm, attr(gm, "groups"))), list(dim = c(7L, 11L), dimnames = list(GRPnames(g2), colnames(m)), groups = attr(gm, "groups")))
    if(Sys.getenv("NCRAN") == "TRUE") {
    expect_identical(attributes(FUN(mtcars, g2, use.g.names = FALSE)), `[[<-`(attributes(mtcars), "row.names", value = seq_len(g2[[1L]])))
    expect_identical(attributes(FUN(mtcars, g2)), `[[<-`(attributes(mtcars), "row.names", value = GRPnames(g2)))
    }
    expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), g2, use.g.names = FALSE)),  `[[<-`(`[[<-`(attributes(gmtc), "row.names", value = seq_len(g2[[1L]])), "class", "data.frame"))
    expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), g2)), `[[<-`(`[[<-`(attributes(gmtc), "row.names", GRPnames(g2)), "class", "data.frame"))
    expect_identical(attributes(FUN(gmtc)), list(names = c(fgroup_vars(gmtc, "names"), names(mtcars)[-c(2,8:9)]), row.names = seq_len(g2[[1L]]), class = "data.frame"))
    expect_identical(attributes(FUN(gmtc, use.g.names = TRUE)), list(names = c(fgroup_vars(gmtc, "names"), names(mtcars)[-c(2,8:9)]), row.names = GRPnames(g2), class = "data.frame"))
    expect_identical(attributes(FUN(gmtc, keep.group_vars = FALSE)), list(names = names(mtcars)[-c(2,8:9)], row.names = seq_len(g2[[1L]]), class = "data.frame"))
    expect_identical(attributes(FUN(gmtc, keep.group_vars = FALSE, use.g.names = TRUE)), list(names = names(mtcars)[-c(2,8:9)], row.names = GRPnames(g2), class = "data.frame"))
  }

  for(i in c("fmode", "ffirst", "flast")) {
    # print(i)
    FUN <- match.fun(i)
    for(k in names(wlddev)) expect_identical(attributes(FUN(wlddev[[k]])), attributes(wlddev[[k]]))
    if(i != "fmode") expect_identical(attributes(FUN(AirPassengers)), NULL)
    expect_identical(attributes(FUN(EuStockMarkets)), list(names = colnames(EuStockMarkets)))
    expect_identical(attributes(FUN(EuStockMarkets, drop = FALSE)), list(dim = c(1L, 4L), dimnames = list(NULL, colnames(EuStockMarkets))))
    expect_identical(attributes(FUN(m)), list(names = colnames(m)))
    expect_identical(attributes(FUN(m, drop = FALSE)), list(dim = c(1L, 11L), dimnames = list(NULL, colnames(m))))
    expect_identical(attributes(FUN(gm)), list(names = colnames(m)))
    expect_identical(attributes(FUN(gm, drop = FALSE)), list(dim = c(1L, 11L), dimnames = list(NULL, colnames(m)), groups = attr(gm, "groups")))
    expect_identical(attributes(FUN(wlddev)), list(names = names(wlddev)))
    expect_identical(attributes(FUN(wlddev, drop = FALSE)), `[[<-`(attributes(wlddev), "row.names", 1L))
    expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"))), list(names = names(mtcars)))
    expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), drop = FALSE)), `[[<-`(`[[<-`(attributes(gmtc), "row.names", 1L), "class", "data.frame"))

    # Grouped
    for(k in names(wlddev)) expect_identical(attributes(FUN(wlddev[[k]], g1, use.g.names = FALSE)), attributes(wlddev[[k]]))
    for(k in names(wlddev)) expect_identical(attributes(FUN(wlddev[[k]], g1)), c(attributes(wlddev[[k]]), list(names = unattrib(GRPnames(g1)))))
    if(i != "fmode") expect_identical(attributes(FUN(AirPassengers, f1, use.g.names = FALSE)), NULL)
    if(i != "fmode") expect_identical(attributes(FUN(AirPassengers, f1)), list(names = as.character(1:5)))
    expect_identical(attributes(FUN(EuStockMarkets, f2, use.g.names = FALSE)), list(dim = c(5L, 4L), dimnames = dimnames(EuStockMarkets)))
    expect_identical(attributes(FUN(EuStockMarkets, f2)), list(dim = c(5L, 4L), dimnames = list(as.character(1:5), colnames(EuStockMarkets))))
    expect_identical(attributes(FUN(m, g2, use.g.names = FALSE)), list(dim = c(7L, 11L), dimnames = list(NULL, colnames(m))))
    expect_identical(attributes(FUN(m, g2)), list(dim = c(7L, 11L), dimnames = list(GRPnames(g2), colnames(m))))
    expect_identical(attributes(FUN(gm, attr(gm, "groups"), use.g.names = FALSE)), list(dim = c(7L, 11L), dimnames = list(NULL, colnames(m)), groups = attr(gm, "groups")))
    expect_identical(attributes(FUN(gm, attr(gm, "groups"))), list(dim = c(7L, 11L), dimnames = list(GRPnames(g2), colnames(m)), groups = attr(gm, "groups")))
    expect_identical(attributes(FUN(wlddev, g1, use.g.names = FALSE)), `[[<-`(attributes(wlddev), "row.names", value = seq_len(g1[[1L]])))
    expect_identical(attributes(FUN(wlddev, g1)), `[[<-`(attributes(wlddev), "row.names", value = GRPnames(g1)))
    if(Sys.getenv("NCRAN") == "TRUE") {
    expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), g2, use.g.names = FALSE)),  `[[<-`(`[[<-`(attributes(gmtc), "row.names", value = seq_len(g2[[1L]])), "class", "data.frame"))
    expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), g2)), `[[<-`(`[[<-`(attributes(gmtc), "row.names", GRPnames(g2)), "class", "data.frame"))
    }
    expect_identical(attributes(FUN(gmtc)), list(names = c(fgroup_vars(gmtc, "names"), names(mtcars)[-c(2,8:9)]), row.names = seq_len(g2[[1L]]), class = "data.frame"))
    expect_identical(attributes(FUN(gmtc, use.g.names = TRUE)), list(names = c(fgroup_vars(gmtc, "names"), names(mtcars)[-c(2,8:9)]), row.names = GRPnames(g2), class = "data.frame"))
    expect_identical(attributes(FUN(gmtc, keep.group_vars = FALSE)), list(names = names(mtcars)[-c(2,8:9)], row.names = seq_len(g2[[1L]]), class = "data.frame"))
    expect_identical(attributes(FUN(gmtc, keep.group_vars = FALSE, use.g.names = TRUE)), list(names = names(mtcars)[-c(2,8:9)], row.names = GRPnames(g2), class = "data.frame"))
  }

  for(i in c("fmin", "fmax")) {
    # print(i)
    FUN <- match.fun(i)
    for(k in num_vars(wlddev, "names")) expect_identical(attributes(FUN(wlddev[[k]])), attributes(wlddev[[k]]))
    expect_identical(attributes(FUN(AirPassengers)), NULL)
    expect_identical(attributes(FUN(EuStockMarkets)), list(names = colnames(EuStockMarkets)))
    expect_identical(attributes(FUN(EuStockMarkets, drop = FALSE)), list(dim = c(1L, 4L), dimnames = list(NULL, colnames(EuStockMarkets))))
    expect_identical(attributes(FUN(m)), list(names = colnames(m)))
    expect_identical(attributes(FUN(m, drop = FALSE)), list(dim = c(1L, 11L), dimnames = list(NULL, colnames(m))))
    expect_identical(attributes(FUN(gm)), list(names = colnames(m)))
    expect_identical(attributes(FUN(gm, drop = FALSE)), list(dim = c(1L, 11L), dimnames = list(NULL, colnames(m)), groups = attr(gm, "groups")))
    expect_identical(attributes(FUN(nv(wlddev))), list(names = nv(wlddev, "names")))
    expect_identical(attributes(FUN(nv(wlddev), drop = FALSE)), `[[<-`(attributes(nv(wlddev)), "row.names", 1L))
    expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"))), list(names = names(mtcars)))
    expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), drop = FALSE)), `[[<-`(`[[<-`(attributes(gmtc), "row.names", 1L), "class", "data.frame"))

    # Grouped
    for(k in num_vars(wlddev, "names")) expect_identical(attributes(FUN(wlddev[[k]], g1, use.g.names = FALSE)), attributes(wlddev[[k]]))
    for(k in num_vars(wlddev, "names")) expect_identical(attributes(FUN(wlddev[[k]], g1)), c(attributes(wlddev[[k]]), list(names = unattrib(GRPnames(g1)))))
    expect_identical(attributes(FUN(AirPassengers, f1, use.g.names = FALSE)), NULL)
    expect_identical(attributes(FUN(AirPassengers, f1)), list(names = as.character(1:5)))
    expect_identical(attributes(FUN(EuStockMarkets, f2, use.g.names = FALSE)), list(dim = c(5L, 4L), dimnames = dimnames(EuStockMarkets)))
    expect_identical(attributes(FUN(EuStockMarkets, f2)), list(dim = c(5L, 4L), dimnames = list(as.character(1:5), colnames(EuStockMarkets))))
    expect_identical(attributes(FUN(m, g2, use.g.names = FALSE)), list(dim = c(7L, 11L), dimnames = list(NULL, colnames(m))))
    expect_identical(attributes(FUN(m, g2)), list(dim = c(7L, 11L), dimnames = list(GRPnames(g2), colnames(m))))
    expect_identical(attributes(FUN(gm, attr(gm, "groups"), use.g.names = FALSE)), list(dim = c(7L, 11L), dimnames = list(NULL, colnames(m)), groups = attr(gm, "groups")))
    expect_identical(attributes(FUN(gm, attr(gm, "groups"))), list(dim = c(7L, 11L), dimnames = list(GRPnames(g2), colnames(m)), groups = attr(gm, "groups")))
    expect_identical(attributes(FUN(nv(wlddev), g1, use.g.names = FALSE)), `[[<-`(attributes(nv(wlddev)), "row.names", value = seq_len(g1[[1L]])))
    expect_identical(attributes(FUN(nv(wlddev), g1)), `[[<-`(attributes(nv(wlddev)), "row.names", value = GRPnames(g1)))
    if(Sys.getenv("NCRAN") == "TRUE") {
      expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), g2, use.g.names = FALSE)),  `[[<-`(`[[<-`(attributes(gmtc), "row.names", value = seq_len(g2[[1L]])), "class", "data.frame"))
      expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), g2)), `[[<-`(`[[<-`(attributes(gmtc), "row.names", GRPnames(g2)), "class", "data.frame"))
    }
    expect_identical(attributes(FUN(gmtc)), list(names = c(fgroup_vars(gmtc, "names"), names(mtcars)[-c(2,8:9)]), row.names = seq_len(g2[[1L]]), class = "data.frame"))
    expect_identical(attributes(FUN(gmtc, use.g.names = TRUE)), list(names = c(fgroup_vars(gmtc, "names"), names(mtcars)[-c(2,8:9)]), row.names = GRPnames(g2), class = "data.frame"))
    expect_identical(attributes(FUN(gmtc, keep.group_vars = FALSE)), list(names = names(mtcars)[-c(2,8:9)], row.names = seq_len(g2[[1L]]), class = "data.frame"))
    expect_identical(attributes(FUN(gmtc, keep.group_vars = FALSE, use.g.names = TRUE)), list(names = names(mtcars)[-c(2,8:9)], row.names = GRPnames(g2), class = "data.frame"))
  }

})

transFUN <- setdiff(c(.FAST_FUN, .OPERATOR_FUN), c(.FAST_STAT_FUN, "fhdbetween", "fhdwithin", "HDB", "HDW"))

options(collapse_unused_arg_action = "none", warn = -1)

test_that("preservation of difftime (and related classes)", {
  v <- diff(wlddev$date)
  av <- attributes(v)
  v <- c(NA, v)
  attributes(v) <- av
  vd <- qDT(v)
  g <- group(wlddev$iso3c)
  w <- abs(rnorm(length(v))) + 5

  for(i in setdiff(.FAST_STAT_FUN, c("fnobs", "fndistinct"))) {
    # print(i)
    FUN <- match.fun(i)
    for(t in list(NULL, "replace_fill")) {
      expect_identical(attributes(FUN(v, TRA = t)), av)
      expect_identical(attributes(FUN(v, g = g, use.g.names = FALSE, TRA = t)), av)
      expect_identical(attributes(FUN(v, na.rm = FALSE, TRA = t)), av)
      expect_identical(attributes(FUN(v, g = g, na.rm = FALSE, use.g.names = FALSE, TRA = t)), av)
      expect_identical(attributes(FUN(vd, drop = FALSE, TRA = t)[[1L]]), av)
      expect_identical(attributes(FUN(vd, g = g, use.g.names = FALSE, TRA = t)[[1L]]), av)
      expect_identical(attributes(FUN(vd, drop = FALSE, na.rm = FALSE, TRA = t)[[1L]]), av)
      expect_identical(attributes(FUN(vd, g = g, na.rm = FALSE, use.g.names = FALSE, TRA = t)[[1L]]), av)
      if(i %in% c("fsum", "fprod", "fmean", "fmedian", "fnth", "fmode", "fvar", "fsd")) {
        expect_identical(attributes(FUN(v, w = w, TRA = t)), av)
        expect_identical(attributes(FUN(v, g = g, w = w, use.g.names = FALSE, TRA = t)), av)
        expect_identical(attributes(FUN(v, w = w, na.rm = FALSE, TRA = t)), av)
        expect_identical(attributes(FUN(v, g = g, w = w, na.rm = FALSE, use.g.names = FALSE, TRA = t)), av)
        expect_identical(attributes(FUN(vd, drop = FALSE, w = w, TRA = t)[[1L]]), av)
        expect_identical(attributes(FUN(vd, g = g, w = w, use.g.names = FALSE, TRA = t)[[1L]]), av)
        expect_identical(attributes(FUN(vd, drop = FALSE, w = w, na.rm = FALSE, TRA = t)[[1L]]), av)
        expect_identical(attributes(FUN(vd, g = g, w = w, na.rm = FALSE, use.g.names = FALSE, TRA = t)[[1L]]), av)
      }
    }
  }

  for(i in c("fnobs", "fndistinct")) {
    FUN <- match.fun(i)
    for(t in list(NULL, "replace_fill")) {
      expect_false(identical(attributes(FUN(v, TRA = t)), av))
      expect_false(identical(attributes(FUN(v, g, use.g.names = FALSE, TRA = t)), av))
      expect_false(identical(attributes(FUN(vd, TRA = t)[[1L]]), av))
      expect_false(identical(attributes(FUN(vd, g, use.g.names = FALSE, TRA = t)[[1L]]), av))
    }
  }

  for(i in setdiff(c(.FAST_FUN, .OPERATOR_FUN), c(.FAST_STAT_FUN, "fhdbetween", "HDB", "fhdwithin", "HDW", "Dlog"))) {
    # print(i)
    FUN <- match.fun(i)
    expect_identical(attributes(FUN(v)), av)
    expect_identical(attributes(FUN(v, g = g, by = g)), av)
    expect_identical(attributes(FUN(vd, cols = NULL)[[1L]]), av)
    expect_identical(attributes(FUN(vd, g = g, by = g, cols = NULL)[[1L]]), av)
    if(i %in% c("fscale", "STD", "fbetween", "B", "fwithin", "W")) {
      expect_identical(attributes(FUN(v, w = w)), av)
      expect_identical(attributes(FUN(v, g = g, by = g, w = w)), av)
      expect_identical(attributes(FUN(vd, w = w, cols = NULL)[[1L]]), av)
      expect_identical(attributes(FUN(vd, g = g, by = g, w = w, cols = NULL)[[1L]]), av)
    }
  }

})

test_that("transformation functions preserve all attributes", {

 for(i in transFUN) {
  # print(i)
  FUN <- match.fun(i)
  for(k in if(i %in% c("flag","L","F")) names(wlddev) else num_vars(wlddev, "names")) expect_identical(attributes(FUN(wlddev[[k]])), attributes(wlddev[[k]]))
  expect_identical(attributes(FUN(AirPassengers)), attributes(AirPassengers))
  expect_identical(attributes(FUN(EuStockMarkets, stubs = FALSE, stub = FALSE)), attributes(EuStockMarkets))
  expect_identical(attributes(FUN(m, stubs = FALSE, stub = FALSE)), attributes(m))
  expect_identical(attributes(FUN(gm, stubs = FALSE, stub = FALSE)), attributes(gm))
  expect_identical(attributes(FUN(if(i == "flag") wlddev else num_vars(wlddev), stubs = FALSE, stub = FALSE)), attributes(if(i == "flag") wlddev else num_vars(wlddev)))
  expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), stubs = FALSE, stub = FALSE)), attributes(`oldClass<-`(gmtc, "data.frame")))

  # Grouped
  for(k in if(i == "flag") names(wlddev) else num_vars(wlddev, "names")) expect_identical(attributes(FUN(wlddev[[k]], g = g1)), attributes(wlddev[[k]]))
  expect_identical(attributes(FUN(AirPassengers, g = f1)), attributes(AirPassengers))
  expect_identical(attributes(FUN(EuStockMarkets, g = f2, stubs = FALSE, stub = FALSE)), attributes(EuStockMarkets))
  expect_identical(attributes(FUN(m, g = g2, stubs = FALSE, stub = FALSE)), attributes(m))
  expect_identical(attributes(FUN(gm, g = attr(gm, "groups"), stubs = FALSE, stub = FALSE)), attributes(gm))
  expect_identical(attributes(FUN(if(i == "flag") wlddev else num_vars(wlddev), g = g1, by = g1, stubs = FALSE, stub = FALSE)), attributes(if(i == "flag") wlddev else num_vars(wlddev)))
  expect_identical(attributes(FUN(`oldClass<-`(gmtc, "data.frame"), g = g2, by = g2, stubs = FALSE, stub = FALSE)), `[[<-`(attributes(gmtc), "class", "data.frame"))
  expect_identical(attributes(if(i %in% c("B","W", "STD")) FUN(gmtc, stub = FALSE) else FUN(gmtc, stubs = FALSE)), `[[<-`(attributes(gmtc), "names", c(fgroup_vars(gmtc, "names"), names(mtcars)[-c(2L ,8:9)])))
  if(i %in% c("fcumsum", "flag", "L", "F", "fdiff", "D", "Dlog", "fgrowth", "G"))
    expect_identical(attributes(FUN(gmtc, keep.ids = FALSE, stubs = FALSE)), `[[<-`(attributes(gmtc), "names", names(mtcars)[-c(2L ,8:9)])) else
    expect_identical(attributes(FUN(gmtc, keep.group_vars = FALSE, stub = FALSE)), `[[<-`(attributes(gmtc), "names", names(mtcars)[-c(2L ,8:9)]))
 }

})

options(collapse_unused_arg_action = "warning", warn = 1)

test_that("TRA attribute preservation works well", {
  # Default Vector Method
  expect_equal(attributes(TRA(AirPassengers, 1, "replace_NA")), attributes(AirPassengers))      # Double
  expect_equal(attributes(TRA(AirPassengers, 1L, "replace_NA"))[[1]], tsp(AirPassengers))       # Integer -> Change of type !!
  expect_equal(attributes(TRA(AirPassengers, 1, "replace")), attributes(AirPassengers))      # Double
  expect_equal(attributes(TRA(AirPassengers, 1L, "replace"))[[1]], tsp(AirPassengers))       # Integer -> Change of type !!
  expect_equal(attributes(TRA(AirPassengers, 1, "replace_fill")), attributes(AirPassengers)) # Double
  expect_equal(attributes(TRA(AirPassengers, 1L, "replace_fill"))[[1]], tsp(AirPassengers))  # Integer -> Change of type !!
  expect_equal(attributes(TRA(AirPassengers, 1, "-")), attributes(AirPassengers))            # Double
  expect_equal(attributes(TRA(AirPassengers, 1L, "-")), attributes(AirPassengers))           # Integer -> Coerced to double in numeric operation
  set.seed(101)
  f <- qF(sample.int(5L, length(AirPassengers), TRUE), na.exclude = FALSE)
  num <- unclass(fmean(AirPassengers, f)); int <- fnobs(AirPassengers, f)
  expect_equal(attributes(TRA(AirPassengers, num, "replace_NA", f)), attributes(AirPassengers))      # Double
  expect_equal(attributes(TRA(AirPassengers, int, "replace_NA", f))[[1]], tsp(AirPassengers))        # Integer -> Change of type !!
  expect_equal(attributes(TRA(AirPassengers, num, "replace", f)), attributes(AirPassengers))      # Double
  expect_equal(attributes(TRA(AirPassengers, int, "replace", f))[[1]], tsp(AirPassengers))        # Integer -> Change of type !!
  expect_equal(attributes(TRA(AirPassengers, num, "replace_fill", f)), attributes(AirPassengers)) # Double
  expect_equal(attributes(TRA(AirPassengers, int, "replace_fill", f))[[1]], tsp(AirPassengers))   # Integer -> Change of type !!
  expect_equal(attributes(TRA(AirPassengers, num, "-", f)), attributes(AirPassengers))            # Double
  expect_equal(attributes(TRA(AirPassengers, int, "-", f)), attributes(AirPassengers))            # Integer -> Coerced to double in numeric operation

  # Matrix Method
  expect_equal(attributes(TRA(EuStockMarkets, rep(1, 4L), "replace_NA")), attributes(EuStockMarkets))          # Double
  expect_equal(attributes(TRA(EuStockMarkets, rep(1L, 4L), "replace_NA"))[["tsp"]], tsp(EuStockMarkets))       # Integer -> Change of type !!
  expect_equal(attributes(TRA(EuStockMarkets, rep(1, 4L), "replace")), attributes(EuStockMarkets))          # Double
  expect_equal(attributes(TRA(EuStockMarkets, rep(1L, 4L), "replace"))[["tsp"]], tsp(EuStockMarkets))       # Integer -> Change of type !!
  expect_equal(attributes(TRA(EuStockMarkets, rep(1, 4L), "replace_fill")), attributes(EuStockMarkets))     # Double
  expect_equal(attributes(TRA(EuStockMarkets, rep(1L, 4L), "replace_fill"))[["tsp"]], tsp(EuStockMarkets))  # Integer -> Change of type !!
  expect_equal(attributes(TRA(EuStockMarkets, rep(1, 4L), "-")), attributes(EuStockMarkets))                # Double
  expect_equal(attributes(TRA(EuStockMarkets, rep(1L, 4L), "-")), attributes(EuStockMarkets))               # Integer -> Coerced to double in numeric operation
  set.seed(101)
  f <- qF(sample.int(5L, nrow(EuStockMarkets), TRUE), na.exclude = FALSE)
  num <- unclass(fmean(EuStockMarkets, f)); int <- fnobs(EuStockMarkets, f)
  expect_equal(attributes(TRA(EuStockMarkets, num, "replace_NA", f)), attributes(EuStockMarkets))         # Double
  expect_equal(attributes(TRA(EuStockMarkets, int, "replace_NA", f))[["tsp"]], tsp(EuStockMarkets))       # Integer -> Change of type !!
  expect_equal(attributes(TRA(EuStockMarkets, num, "replace", f)), attributes(EuStockMarkets))         # Double
  expect_equal(attributes(TRA(EuStockMarkets, int, "replace", f))[["tsp"]], tsp(EuStockMarkets))       # Integer -> Change of type !!
  expect_equal(attributes(TRA(EuStockMarkets, num, "replace_fill", f)), attributes(EuStockMarkets))    # Double
  expect_equal(attributes(TRA(EuStockMarkets, int, "replace_fill", f))[["tsp"]], tsp(EuStockMarkets))  # Integer -> Change of type !!
  expect_equal(attributes(TRA(EuStockMarkets, num, "-", f)), attributes(EuStockMarkets))               # Double
  expect_equal(attributes(TRA(EuStockMarkets, int, "-", f)), attributes(EuStockMarkets))               # Integer -> Coerced to double in numeric operation

  # Data Frame Method
  # CATEGORICAL
  # Simple
  expect_equal(vclasses(unattrib(fndistinct(wlddev, TRA = "replace_fill"))), rep("integer", 13L))
  expect_equal(vclasses(unattrib(fndistinct(wlddev, TRA = "replace"))), rep("integer", 13L))
  expect_equal(vclasses(unattrib(fnobs(wlddev, TRA = "replace_fill"))), rep("integer", 13L))
  expect_equal(vclasses(unattrib(fnobs(wlddev, TRA = "replace"))), rep("integer", 13L))
  expect_equal(lapply(ffirst(wlddev, TRA = "replace_fill"), attributes), lapply(wlddev, attributes))
  expect_equal(lapply(ffirst(wlddev, TRA = "replace"), attributes), lapply(wlddev, attributes))
  expect_equal(lapply(flast(wlddev, TRA = "replace_fill"), attributes), lapply(wlddev, attributes))
  expect_equal(lapply(flast(wlddev, TRA = "replace"), attributes), lapply(wlddev, attributes))
  expect_equal(lapply(fmode(wlddev, TRA = "replace_fill"), attributes), lapply(wlddev, attributes))
  expect_equal(lapply(fmode(wlddev, TRA = "replace"), attributes), lapply(wlddev, attributes))
  # Grouped
  expect_equal(vclasses(unattrib(fndistinct(wlddev, wlddev$iso3c, TRA = "replace_fill"))), rep("integer", 13L))
  expect_equal(vclasses(unattrib(fndistinct(wlddev, wlddev$iso3c, TRA = "replace"))), rep("integer", 13L))
  expect_equal(vclasses(unattrib(fnobs(wlddev, wlddev$iso3c, TRA = "replace_fill"))), rep("integer", 13L))
  expect_equal(vclasses(unattrib(fnobs(wlddev, wlddev$iso3c, TRA = "replace"))), rep("integer", 13L))
  expect_equal(lapply(ffirst(wlddev, wlddev$iso3c, TRA = "replace_fill"), attributes), lapply(wlddev, attributes))
  expect_equal(lapply(ffirst(wlddev, wlddev$iso3c, TRA = "replace"), attributes), lapply(wlddev, attributes))
  expect_equal(lapply(flast(wlddev, wlddev$iso3c, TRA = "replace_fill"), attributes), lapply(wlddev, attributes))
  expect_equal(lapply(flast(wlddev, wlddev$iso3c, TRA = "replace"), attributes), lapply(wlddev, attributes))
  expect_equal(lapply(fmode(wlddev, wlddev$iso3c, TRA = "replace_fill"), attributes), lapply(wlddev, attributes))
  expect_equal(lapply(fmode(wlddev, wlddev$iso3c, TRA = "replace"), attributes), lapply(wlddev, attributes))

  # Numeric
  nwld <- num_vars(wlddev)
  # Simple
  expect_equal(vclasses(fndistinct(nwld, TRA = "replace_NA")), vclasses(nwld))
  expect_equal(vclasses(unattrib(fndistinct(nwld, TRA = "replace_fill"))), rep("integer", fncol(nwld)))
  expect_equal(vclasses(unattrib(fndistinct(nwld, TRA = "replace"))), rep("integer", fncol(nwld)))
  expect_equal(vclasses(unattrib(fndistinct(nwld, TRA = "-"))), rep("numeric", fncol(nwld)))
  expect_equal(vclasses(fnobs(nwld, TRA = "replace_NA")), vclasses(nwld))
  expect_equal(vclasses(unattrib(fnobs(nwld, TRA = "replace_fill"))), rep("integer", fncol(nwld)))
  expect_equal(vclasses(unattrib(fnobs(nwld, TRA = "replace"))), rep("integer", fncol(nwld)))
  expect_equal(vclasses(unattrib(fnobs(nwld, TRA = "%%"))), rep("numeric", fncol(nwld)))
  expect_equal(lapply(fmean(nwld, TRA = "replace_NA"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmean(nwld, TRA = "replace_fill"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmean(nwld, TRA = "replace"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmean(nwld, TRA = "+"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fsd(nwld, TRA = "replace_NA"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fsd(nwld, TRA = "replace_fill"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fsd(nwld, TRA = "replace"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fsd(nwld, TRA = "/"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmedian(nwld, TRA = "replace_NA"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmedian(nwld, TRA = "replace_fill"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmedian(nwld, TRA = "replace"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmedian(nwld, TRA = "-"), attributes), lapply(nwld, attributes))
  # Grouped
  expect_equal(vclasses(fndistinct(nwld, wlddev$iso3c, TRA = "replace_NA")), vclasses(nwld))
  expect_equal(vclasses(unattrib(fndistinct(nwld, wlddev$iso3c, TRA = "replace_fill"))), rep("integer", fncol(nwld)))
  expect_equal(vclasses(unattrib(fndistinct(nwld, wlddev$iso3c, TRA = "replace"))), rep("integer", fncol(nwld)))
  expect_equal(vclasses(unattrib(fndistinct(nwld, wlddev$iso3c, TRA = "-%%"))), rep("numeric", fncol(nwld)))
  expect_equal(vclasses(fnobs(nwld, wlddev$iso3c, TRA = "replace_NA")), vclasses(nwld))
  expect_equal(vclasses(unattrib(fnobs(nwld, wlddev$iso3c, TRA = "replace_fill"))), rep("integer", fncol(nwld)))
  expect_equal(vclasses(unattrib(fnobs(nwld, wlddev$iso3c, TRA = "replace"))), rep("integer", fncol(nwld)))
  expect_equal(vclasses(unattrib(fnobs(nwld, wlddev$iso3c, TRA = "*"))), rep("numeric", fncol(nwld)))
  expect_equal(lapply(fmean(nwld, wlddev$iso3c, TRA = "replace_NA"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmean(nwld, wlddev$iso3c, TRA = "replace_fill"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmean(nwld, wlddev$iso3c, TRA = "replace"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmean(nwld, wlddev$iso3c, TRA = "-+"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fsd(nwld, wlddev$iso3c, TRA = "replace_NA"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fsd(nwld, wlddev$iso3c, TRA = "replace_fill"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fsd(nwld, wlddev$iso3c, TRA = "replace"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fsd(nwld, wlddev$iso3c, TRA = "/"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmedian(nwld, wlddev$iso3c, TRA = "replace_NA"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmedian(nwld, wlddev$iso3c, TRA = "replace_fill"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmedian(nwld, wlddev$iso3c, TRA = "replace"), attributes), lapply(nwld, attributes))
  expect_equal(lapply(fmedian(nwld, wlddev$iso3c, TRA = "+"), attributes), lapply(nwld, attributes))

})

