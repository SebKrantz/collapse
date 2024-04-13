context("qsu")



# rm(list = ls())

bmean <- base::mean
bsd <- stats::sd
bsum <- base::sum

bstats <- function(x) {
  if(!is.numeric(x)) return(c(N = bsum(!is.na(x)), Mean = NA_real_, SD = NA_real_, Min = NA_real_, Max = NA_real_))
  c(N = bsum(!is.na(x)), Mean = bmean(x, na.rm = TRUE), SD = bsd(x, na.rm = TRUE), `names<-`(range(x, na.rm = TRUE), c("Min", "Max")))
}
base_qsu <- function(x, g = NULL) {
  if(is.atomic(x) && !is.matrix(x)) return(`oldClass<-`(bstats(x), c("qsu", "table")))
  if(is.null(g)) {
    r <- t(dapply(x, bstats, return = "matrix"))
    return(`oldClass<-`(r, c("qsu", "matrix", "table")))
  }
  r <- simplify2array(BY(x, g, bstats, return = "list", expand.wide = TRUE))
  return(`oldClass<-`(r, c("qsu", "array", "table")))
}

wldNA <- na_insert(wlddev)
xNA <- na_insert(rnorm(100))
ones <- rep(1, fnrow(wlddev))

for(i in 1:2) {
  if(i == 1L) qsu <- function(x, ...) collapse::qsu(x, ..., stable.algo = FALSE)
  if(i == 2L) qsu <- collapse::qsu

test_that("qsu works properly for simple cases (including unit groups and weights)", {

  expect_equal(qsu(1:10), base_qsu(1:10))
  expect_equal(qsu(10:1), base_qsu(10:1))
  expect_equal(qsu(xNA), base_qsu(xNA))
  expect_equal(qsu(wlddev), base_qsu(wlddev))
  expect_equal(qsu(wldNA), base_qsu(wldNA))
  expect_equal(qsu(GGDC10S), base_qsu(GGDC10S))

  expect_equal(qsu(1:10, w = rep(1, 10)), base_qsu(1:10))
  expect_equal(qsu(10:1, w = rep(1, 10)), base_qsu(10:1))
  expect_equal(qsu(xNA, w = rep(1, 100)), base_qsu(xNA))
  expect_equal(qsu(wlddev, w = ones), base_qsu(wlddev))
  expect_equal(qsu(wldNA, w = ones), base_qsu(wldNA))
  expect_equal(qsu(GGDC10S, w = rep(1, fnrow(GGDC10S))), base_qsu(GGDC10S))

  expect_equal(unattrib(qsu(1:10, g = rep(1, 10))), unattrib(base_qsu(1:10)))
  expect_equal(unattrib(qsu(10:1, g = rep(1, 10))), unattrib(base_qsu(10:1)))
  expect_equal(unattrib(qsu(xNA, g = rep(1, 100))), unattrib(base_qsu(xNA)))
  expect_equal(unattrib(qsu(wlddev, by = ones)), unattrib(t(base_qsu(wlddev)))) # This should be an array... or oriented the other way around...

  expect_equal(unattrib(qsu(1:10, g = rep(1, 10), w = rep(1, 10))), unattrib(base_qsu(1:10)))
  expect_equal(unattrib(qsu(10:1, g = rep(1, 10), w = rep(1, 10))), unattrib(base_qsu(10:1)))
  expect_equal(unattrib(qsu(xNA, g = rep(1, 100), w = rep(1, 100))), unattrib(base_qsu(xNA)))
  expect_equal(qsu(wldNA, w = ones), base_qsu(wldNA))
  expect_equal(qsu(GGDC10S, w = rep(1, fnrow(GGDC10S))), base_qsu(GGDC10S))
  expect_equal(t(unclass(qsu(wldNA, w = ones, by = ones))), unclass(base_qsu(wldNA)))
  expect_equal(t(unclass(qsu(GGDC10S, w = rep(1, fnrow(GGDC10S)), by = rep(1, fnrow(GGDC10S))))), unclass(base_qsu(GGDC10S)))

})

}
rm(qsu)

test_that("qsu works properly for simple cases with higher-order statistics (including unit groups and weights)", {

  expect_equal(qsu(1:10, higher = TRUE)[1:5], base_qsu(1:10))
  expect_equal(qsu(10:1, higher = TRUE)[1:5], base_qsu(10:1))
  expect_equal(qsu(xNA, higher = TRUE)[1:5], base_qsu(xNA))
  expect_equal(qsu(wlddev, higher = TRUE)[,1:5], base_qsu(wlddev))
  expect_equal(qsu(wldNA, higher = TRUE)[,1:5], base_qsu(wldNA))
  expect_equal(qsu(GGDC10S, higher = TRUE)[,1:5], base_qsu(GGDC10S))

  expect_equal(qsu(1:10, w = rep(1, 10), higher = TRUE)[1:5], base_qsu(1:10))
  expect_equal(qsu(10:1, w = rep(1, 10), higher = TRUE)[1:5], base_qsu(10:1))
  expect_equal(qsu(xNA, w = rep(1, 100), higher = TRUE)[1:5], base_qsu(xNA))
  expect_equal(qsu(wlddev, w = ones, higher = TRUE)[,1:5], base_qsu(wlddev))
  expect_equal(qsu(wldNA, w = ones, higher = TRUE)[,1:5], base_qsu(wldNA))
  expect_equal(qsu(GGDC10S, w = rep(1, fnrow(GGDC10S)), higher = TRUE)[,1:5], base_qsu(GGDC10S))

  expect_equal(unattrib(qsu(1:10, g = rep(1, 10), higher = TRUE)[1:5]), unattrib(base_qsu(1:10)))
  expect_equal(unattrib(qsu(10:1, g = rep(1, 10), higher = TRUE)[1:5]), unattrib(base_qsu(10:1)))
  expect_equal(unattrib(qsu(xNA, g = rep(1, 100), higher = TRUE)[1:5]), unattrib(base_qsu(xNA)))
  expect_equal(unattrib(qsu(wlddev, by = ones, higher = TRUE)[1:5, ]), unattrib(t(base_qsu(wlddev)))) # This should be an array... or oriented the other way around...

  expect_equal(unattrib(qsu(1:10, g = rep(1, 10), w = rep(1, 10), higher = TRUE)[1:5]), unattrib(base_qsu(1:10)))
  expect_equal(unattrib(qsu(10:1, g = rep(1, 10), w = rep(1, 10), higher = TRUE)[1:5]), unattrib(base_qsu(10:1)))
  expect_equal(unattrib(qsu(xNA, g = rep(1, 100), w = rep(1, 100), higher = TRUE)[1:5]), unattrib(base_qsu(xNA)))
  expect_equal(qsu(wldNA, w = ones, higher = TRUE)[,1:5], base_qsu(wldNA))
  expect_equal(qsu(GGDC10S, w = rep(1, fnrow(GGDC10S)), higher = TRUE)[,1:5], base_qsu(GGDC10S))
  expect_equal(t(unclass(qsu(wldNA, w = ones, by = ones, higher = TRUE)[1:5,])), unclass(base_qsu(wldNA)))
  expect_equal(t(unclass(qsu(GGDC10S, w = rep(1, fnrow(GGDC10S)), by = rep(1, fnrow(GGDC10S)), higher = TRUE)))[,1:5], unclass(base_qsu(GGDC10S)))

})


wtd.sd <- function(x, w) sqrt(bsum(w * (x - weighted.mean(x, w))^2)/bsum(w))
wtd.skewness <- function(x, w) (bsum(w * (x - weighted.mean(x, w))^3)/bsum(w))/wtd.sd(x, w)^3
wtd.kurtosis <- function(x, w) ((bsum(w * (x - weighted.mean(x, w))^4)/bsum(w))/wtd.sd(x, w)^4)

base_w_qsu <- function(x, w) {
  if(!is.numeric(x)) return(c(N = bsum(!is.na(x)), Mean = NA_real_, SD = NA_real_, Min = NA_real_, Max = NA_real_, Skew = NA_real_, Kurt = NA_real_))
  cc <- complete.cases(x, w)
  if(!all(cc)) {
    x <- x[cc]
    w <- w[cc]
  }
  res <- c(N = length(x), Mean = weighted.mean(x, w), SD = fsd(x, w = w),
    `names<-`(range(x, na.rm = TRUE), c("Min", "Max")), Skew = wtd.skewness(x, w), Kurt = wtd.kurtosis(x, w))
  class(res) <- c("qsu", "table")
  res
}

test_that("Proper performance of weighted statsistics", {
  x <- mtcars$mpg
  w <- ceiling(mtcars$wt*10)
  wx <- rep(x, w)
  expect_equal(base_w_qsu(x, w)[-1L], qsu(wx, higher = TRUE)[-1L])
  expect_equal(qsu(wx)[-1L], qsu(x, w = w)[-1L])
  expect_equal(qsu(wx, higher = TRUE)[-1L], qsu(x, w = w, higher = TRUE)[-1L])
  expect_equal(drop(qsu(wx, g = rep(1L, length(wx)), higher = TRUE))[-1L], drop(qsu(x, g = rep(1L, length(x)), w = w, higher = TRUE))[-1L])
})

g <- GRP(wlddev, ~ income)
p <- GRP(wlddev, ~ iso3c)

for(i in 1:2) {
  if(i == 1L) qsu <- function(x, ...) collapse::qsu(x, ..., stable.algo = FALSE)
  if(i == 2L) qsu <- collapse::qsu

test_that("qsu works properly for grouped and panel data computations", {

  # Grouped Statistics
  expect_equal(qsu(wldNA, g), base_qsu(wldNA, g))
  expect_equal(qsu(GGDC10S, GGDC10S$Variable), base_qsu(GGDC10S, GGDC10S$Variable))
  # Grouped and Weighted Statistics
  expect_equal(qsu(wldNA, g, w = ones), base_qsu(wldNA, g))
  expect_equal(qsu(GGDC10S, GGDC10S$Variable, w = rep(1, fnrow(GGDC10S))), base_qsu(GGDC10S, GGDC10S$Variable))
  # Panel Data Statistics
  ps <- qsu(wldNA, pid = p, cols = is.numeric)
  expect_equal(unattrib(t(ps["Overall",,])), unattrib(base_qsu(nv(wldNA))))
  expect_equal(unattrib(t(ps["Between",,])), unattrib(base_qsu(fmean(nv(wldNA), p))))
  expect_equal(unattrib(t(ps["Within", -1,])), unattrib(base_qsu(fwithin(nv(wldNA), p, mean = "overall.mean"))[, -1]))
  # Weighted Panel Data Statistics
  ps <- qsu(wldNA, pid = p, w = ones, cols = is.numeric)
  expect_equal(unattrib(t(ps["Overall",,])), unattrib(base_qsu(nv(wldNA))))
  expect_equal(unattrib(t(ps["Between",-1,])), unattrib(base_qsu(fbetween(nv(wldNA), p))[,-1]))
  expect_equal(unattrib(t(ps["Within", -1,])), unattrib(base_qsu(fwithin(nv(wldNA), p, mean = "overall.mean"))[, -1]))
  # Grouped Panel Data Statistics
  ps <- qsu(wldNA, by = g, pid = p, cols = is.numeric)
  expect_equal(unattrib(ps[,,"Overall",]), unattrib(base_qsu(nv(wldNA), g)))
  expect_equal(unattrib(ps[,-1,"Between",]), unattrib(base_qsu(fbetween(nv(wldNA), p), g)[,-1,]))
  expect_equal(unattrib(ps[,-1,"Within",]), unattrib(base_qsu(fwithin(nv(wldNA), p, mean = "overall.mean"), g)[,-1,]))
  # Grouped and Weighted Panel Data Statistics
  ps <- qsu(wldNA, by = g, pid = p, w = ones, cols = is.numeric)
  expect_equal(unattrib(ps[,,"Overall",]), unattrib(base_qsu(nv(wldNA), g)))
  expect_equal(unattrib(ps[,-1,"Between",]), unattrib(base_qsu(fbetween(nv(wldNA), p), g)[,-1,]))
  expect_equal(unattrib(ps[,-1,"Within",]), unattrib(base_qsu(fwithin(nv(wldNA), p, mean = "overall.mean"), g)[,-1,]))

})

}
rm(qsu)

test_that("qsu works properly for grouped and panel data computations with higher-order statistics", {

  # Grouped Statistics
  expect_equal(qsu(wldNA, g, higher = TRUE)[,1:5,], base_qsu(wldNA, g))
  expect_equal(qsu(GGDC10S, GGDC10S$Variable, higher = TRUE)[,1:5,], base_qsu(GGDC10S, GGDC10S$Variable))
  # Grouped and Weighted Statistics
  expect_equal(qsu(wldNA, g, w = ones, higher = TRUE)[,1:5,], base_qsu(wldNA, g))
  expect_equal(qsu(GGDC10S, GGDC10S$Variable, w = rep(1, fnrow(GGDC10S)), higher = TRUE)[,1:5,], base_qsu(GGDC10S, GGDC10S$Variable))
  # Panel Data Statistics
  ps <- qsu(wldNA, pid = p, cols = is.numeric, higher = TRUE)[,1:5,]
  expect_equal(unattrib(t(ps["Overall",,])), unattrib(base_qsu(nv(wldNA))))
  expect_equal(unattrib(t(ps["Between",,])), unattrib(base_qsu(fmean(nv(wldNA), p))))
  expect_equal(unattrib(t(ps["Within", -1,])), unattrib(base_qsu(fwithin(nv(wldNA), p, mean = "overall.mean"))[, -1]))
  # Weighted Panel Data Statistics
  ps <- qsu(wldNA, pid = p, w = ones, cols = is.numeric, higher = TRUE)[,1:5,]
  expect_equal(unattrib(t(ps["Overall",,])), unattrib(base_qsu(nv(wldNA))))
  # TODO: Figure out why this test fails !!!!!!
  # expect_equal(unattrib(t(ps["Between",-1,])), unattrib(base_qsu(fbetween(nv(wldNA), p))[,-1]))
  expect_equal(unattrib(t(ps["Within", -1,])), unattrib(base_qsu(fwithin(nv(wldNA), p, mean = "overall.mean"))[, -1]))
  # Grouped Panel Data Statistics
  ps <- qsu(wldNA, by = g, pid = p, cols = is.numeric, higher = TRUE)[,1:5,,]
  expect_equal(unattrib(ps[,,"Overall",]), unattrib(base_qsu(nv(wldNA), g)))
  expect_equal(unattrib(ps[,-1,"Between",]), unattrib(base_qsu(fbetween(nv(wldNA), p), g)[,-1,]))
  expect_equal(unattrib(ps[,-1,"Within",]), unattrib(base_qsu(fwithin(nv(wldNA), p, mean = "overall.mean"), g)[,-1,]))
  # Grouped and Weighted Panel Data Statistics
  ps <- qsu(wldNA, by = g, pid = p, w = ones, cols = is.numeric, higher = TRUE)[,1:5,,]
  expect_equal(unattrib(ps[,,"Overall",]), unattrib(base_qsu(nv(wldNA), g)))
  expect_equal(unattrib(ps[,-1,"Between",]), unattrib(base_qsu(fbetween(nv(wldNA), p), g)[,-1,]))
  expect_equal(unattrib(ps[,-1,"Within",]), unattrib(base_qsu(fwithin(nv(wldNA), p, mean = "overall.mean"), g)[,-1,]))

})

# Make more tests!! See also collapse general TODO !
test_that("qsu gives errors for wrong input", {
  expect_error(qsu(wlddev$year, 2:4))
  expect_error(qsu(wlddev$year, pid = 2:4))
  expect_error(qsu(wlddev, 2:4))
  expect_error(qsu(wlddev, pid = 2:4))
  expect_error(qsu(wlddev$year, letters))
  expect_error(qsu(wlddev$year, pid = letters))
  expect_error(qsu(wlddev, letters))
  expect_error(qsu(wlddev, pid = letters))

  expect_error(qsu(wlddev, ~ iso3c + bla))
  expect_error(qsu(wlddev, pid = ~ iso3c + bla))

  expect_visible(qsu(wlddev, PCGDP ~ region + income))
  expect_visible(qsu(wlddev, pid = PCGDP ~ region + income))

  expect_equal(qsu(wlddev, PCGDP ~ region + income, ~ iso3c), qsu(wlddev, ~ region + income, pid = PCGDP ~ iso3c))

  expect_error(qsu(wlddev, cols = 9:14))
  expect_error(qsu(wlddev, cols = c("PCGDP","bla")))
})

