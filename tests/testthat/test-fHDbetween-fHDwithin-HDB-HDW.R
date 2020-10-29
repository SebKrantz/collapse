context("fHDbetween / HDB and fHDwithin / HDW")

# rm(list = ls())

options(warn = -1)

x <- rnorm(100)
w <- abs(100*rnorm(100))
wdat <- abs(100*rnorm(32))
xNA <- x
wNA <- w
wdatNA <- wdat
xNA[sample.int(100,20)] <- NA
wNA[sample.int(100,20)] <- NA
wdatNA[sample.int(32, 5)] <- NA
f <- as.factor(rep(1:10, each = 10))
g <- as.factor(rep(c(1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,7,7,7,7,7,7,7,10,10,10,10,10,10,10,10,10,10)))
mtcNA <- na_insert(mtcars)
mtcNA[1,1] <- NA # single group NA !!
m <- as.matrix(mtcars)
mNA <- as.matrix(mtcNA)
mNAc <- mNA
storage.mode(mNAc) <- "character"

baseresid <- function(y, X, na.rm = FALSE) {
  y <- qM(y)
  if(is.list(X)) X <- do.call(cbind, X)
  X <- cbind(Intercept = 1L, X)
  if(na.rm) {
    cc <- complete.cases(y, X)
    y <- y[cc, , drop = FALSE]
    X <- X[cc, , drop = FALSE]
  }
  drop(qr.resid(qr.default(X), y))
}

basefitted <- function(y, X, na.rm = FALSE) {
  y <- qM(y)
  if(is.list(X)) X <- do.call(cbind, X)
  X <- cbind(Intercept = 1L, X)
  if(na.rm) {
    cc <- complete.cases(y, X)
    y <- y[cc, , drop = FALSE]
    X <- X[cc, , drop = FALSE]
  }
  drop(qr.fitted(qr.default(X), y))
}

# fHDbetween and fHDwithin

test_that("fHDbetween with one factor performs like fbetween", {
  expect_equal(fHDbetween(x, f), fbetween(x, f))
  expect_equal(fHDbetween(x, f, na.rm = FALSE), fbetween(x, f, na.rm = FALSE))
  expect_equal(fHDbetween(xNA, f, na.rm = FALSE), fbetween(xNA, f, na.rm = FALSE))
  expect_equal(`attributes<-`(fHDbetween(xNA, f, fill = TRUE), NULL), fbetween(xNA, f))
  expect_equal(fHDbetween(m, g), fbetween(m, g))
  expect_equal(fHDbetween(m, g, na.rm = FALSE), fbetween(m, g, na.rm = FALSE))
  expect_equal(fHDbetween(mNA, g, na.rm = FALSE), fbetween(mNA, g, na.rm = FALSE))
  # expect_equal(fHDbetween(mNA, g, fill = TRUE), fbetween(mNA, g)) # not matching, fHDbetween matrix is not variable.wise
  expect_equal(fHDbetween(mtcars, g), fbetween(mtcars, g))
  expect_equal(fHDbetween(mtcars, g, na.rm = FALSE), fbetween(mtcars, g, na.rm = FALSE))
  expect_equal(fHDbetween(mtcNA, g, na.rm = FALSE), fbetween(mtcNA, g, na.rm = FALSE))
  expect_equal(fHDbetween(mtcNA, g, variable.wise = TRUE), fbetween(mtcNA, g))
})

test_that("fHDwithin with one factor performs like fwithin", {
  expect_equal(fHDwithin(x, f), fwithin(x, f))
  expect_equal(fHDwithin(x, f, na.rm = FALSE), fwithin(x, f, na.rm = FALSE))
  expect_equal(fHDwithin(xNA, f, na.rm = FALSE), fwithin(xNA, f, na.rm = FALSE))
  expect_equal(`attributes<-`(fHDwithin(xNA, f, fill = TRUE), NULL), fwithin(xNA, f))
  expect_equal(fHDwithin(m, g), fwithin(m, g))
  expect_equal(fHDwithin(m, g, na.rm = FALSE), fwithin(m, g, na.rm = FALSE))
  expect_equal(fHDwithin(mNA, g, na.rm = FALSE), fwithin(mNA, g, na.rm = FALSE))
  # expect_equal(fHDwithin(mNA, g, fill = TRUE), fwithin(mNA, g)) # not matching, fHDwithin matrix is not variable.wise
  expect_equal(fHDwithin(mtcars, g), fwithin(mtcars, g))
  expect_equal(fHDwithin(mtcars, g, na.rm = FALSE), fwithin(mtcars, g, na.rm = FALSE))
  expect_equal(fHDwithin(mtcNA, g, na.rm = FALSE), fwithin(mtcNA, g, na.rm = FALSE))
  expect_equal(fHDwithin(mtcNA, g, variable.wise = TRUE), fwithin(mtcNA, g))
})

f2 <- qF(sample.int(10, 100, TRUE))
fl <- list(f, f2)

g2 <- qF(sample.int(5, 32, TRUE))
gl <- list(g, g2)

test_that("fHDbetween with two factors performs like lfe::demeanlist", {
  expect_equal(fHDbetween(x, fl), lfe::demeanlist(x, fl, means = TRUE), tolerance = 1e-5)
  expect_equal(fHDbetween(xNA, fl), lfe::demeanlist(xNA, fl, means = TRUE, na.rm = TRUE), tolerance = 1e-5)
  expect_visible(fHDbetween(xNA, fl, fill = TRUE))
  expect_equal(fHDbetween(m, gl), lfe::demeanlist(m, gl, means = TRUE), tolerance = 1e-5)
  expect_equal(fHDbetween(mNA, gl, na.rm = FALSE), lfe::demeanlist(mNA, gl, means = TRUE), tolerance = 1e-5)
  expect_equal(fHDbetween(mNA, gl), lfe::demeanlist(mNA, gl, means = TRUE, na.rm = TRUE), tolerance = 1e-5)
  expect_visible(fHDbetween(mNA, gl, fill = TRUE))
  expect_equal(fHDbetween(mtcars, gl), lfe::demeanlist(mtcars, gl, means = TRUE), tolerance = 1e-5)
  expect_equal(fHDbetween(mtcNA, gl, na.rm = FALSE), lfe::demeanlist(mtcNA, gl, means = TRUE), tolerance = 1e-5)
  expect_equal(setRownames(fHDbetween(mtcNA, gl)), lfe::demeanlist(mtcNA, gl, means = TRUE, na.rm = TRUE), tolerance = 1e-5)
  expect_visible(fHDbetween(mtcNA, gl, fill = TRUE))
  expect_visible(fHDbetween(mtcNA, gl, variable.wise = TRUE))
})

test_that("fHDwithin with two factors performs like lfe::demeanlist", {
  expect_equal(fHDwithin(x, fl), lfe::demeanlist(x, fl), tolerance = 1e-5)
  expect_equal(fHDwithin(xNA, fl), lfe::demeanlist(xNA, fl, na.rm = TRUE), tolerance = 1e-5)
  expect_visible(fHDwithin(xNA, fl, fill = TRUE))
  expect_equal(fHDwithin(m, gl), lfe::demeanlist(m, gl), tolerance = 1e-5)
  expect_equal(fHDwithin(mNA, gl, na.rm = FALSE), lfe::demeanlist(mNA, gl), tolerance = 1e-5)
  expect_equal(fHDwithin(mNA, gl), lfe::demeanlist(mNA, gl, na.rm = TRUE), tolerance = 1e-5)
  expect_visible(fHDwithin(mNA, gl, fill = TRUE))
  expect_equal(fHDwithin(mtcars, gl), lfe::demeanlist(mtcars, gl), tolerance = 1e-5)
  expect_equal(fHDwithin(mtcNA, gl, na.rm = FALSE), lfe::demeanlist(mtcNA, gl), tolerance = 1e-5)
  expect_equal(setRownames(fHDwithin(mtcNA, gl)), lfe::demeanlist(mtcNA, gl, na.rm = TRUE), tolerance = 1e-5)
  expect_visible(fHDwithin(mtcNA, gl, fill = TRUE))
  expect_visible(fHDwithin(mtcNA, gl, variable.wise = TRUE))
})

x2 <- 3 * x + rnorm(100)

test_that("fHDbetween with only continuous variables performs like basefitted (defined above)", {
  expect_equal(fHDbetween(x, x2), basefitted(x, x2), tolerance = 1e-5)
  expect_equal(`attr<-`(fHDbetween(xNA, x2), "na.rm", NULL), basefitted(xNA, x2, na.rm = TRUE), tolerance = 1e-5)
  expect_visible(fHDbetween(xNA, x2, fill = TRUE))
  expect_equal(fHDbetween(m, m), fHDbetween(m, mtcars), tolerance = 1e-5)
  expect_equal(fHDbetween(m, m), basefitted(m, m), tolerance = 1e-5)
  expect_equal(`attr<-`(fHDbetween(mNA, m), "na.rm", NULL), basefitted(mNA, m, na.rm = TRUE), tolerance = 1e-5)
  expect_equal(fHDbetween(mNA, m, fill = TRUE), fHDbetween(mNA, mtcars, fill = TRUE), tolerance = 1e-5)
  expect_equal(fHDbetween(mtcars, mtcars), fHDbetween(mtcars, m), tolerance = 1e-5)
  expect_equal(fHDbetween(mtcars, mtcars), qDF(basefitted(mtcars, mtcars)), tolerance = 1e-5)
  expect_equal(`attr<-`(fHDbetween(mtcNA, mtcars), "na.rm", NULL), qDF(basefitted(mtcNA, mtcars, na.rm = TRUE)), tolerance = 1e-5)
  expect_equal(fHDbetween(mtcNA, mtcars, fill = TRUE), fHDbetween(mtcNA, m, fill = TRUE), tolerance = 1e-5)
  expect_equal(fHDbetween(mtcNA, mtcars, variable.wise = TRUE), fHDbetween(mtcNA, m, variable.wise = TRUE), tolerance = 1e-5)
})

test_that("fHDwithin with only continuous variables performs like baseresid (defined above)", {
  expect_equal(fHDwithin(x, x2), baseresid(x, x2), tolerance = 1e-5)
  expect_equal(`attr<-`(fHDwithin(xNA, x2), "na.rm", NULL), baseresid(xNA, x2, na.rm = TRUE), tolerance = 1e-5)
  expect_visible(fHDwithin(xNA, x2, fill = TRUE))
  expect_equal(fHDwithin(m, m), fHDwithin(m, mtcars), tolerance = 1e-5)
  expect_equal(fHDwithin(m, m), baseresid(m, m), tolerance = 1e-5)
  expect_equal(`attr<-`(fHDwithin(mNA, m), "na.rm", NULL), baseresid(mNA, m, na.rm = TRUE), tolerance = 1e-5)
  expect_equal(fHDwithin(mNA, m, fill = TRUE), fHDwithin(mNA, mtcars, fill = TRUE), tolerance = 1e-5)
  expect_equal(fHDwithin(mtcars, mtcars), fHDwithin(mtcars, m), tolerance = 1e-5)
  expect_equal(fHDwithin(mtcars, mtcars), qDF(baseresid(mtcars, mtcars)), tolerance = 1e-5)
  expect_equal(`attr<-`(fHDwithin(mtcNA, mtcars), "na.rm", NULL), qDF(baseresid(mtcNA, mtcars, na.rm = TRUE)), tolerance = 1e-5)
  expect_equal(fHDwithin(mtcNA, mtcars, fill = TRUE), fHDwithin(mtcNA, m, fill = TRUE), tolerance = 1e-5)
  expect_equal(fHDwithin(mtcNA, mtcars, variable.wise = TRUE), fHDwithin(mtcNA, m, variable.wise = TRUE), tolerance = 1e-5)
})

data <- wlddev
data$year <- qF(data$year)
data <- get_vars(data, c("iso3c","year","region","income","PCGDP","LIFEEX","ODA"))

test_that("fHDbetween with multiple variables performs like lm", {
  expect_equal(fHDbetween(iris$Sepal.Length, iris[-1]), `names<-`(fitted(lm(Sepal.Length ~., iris)), NULL), tolerance = 1e-5)
  expect_equal(fHDbetween(iris[1], iris[-1])[[1]], `names<-`(fitted(lm(Sepal.Length ~., iris)), NULL), tolerance = 1e-5)
  expect_equal(setRownames(qM(fHDbetween(iris[1:2], iris[-(1:2)]))), fitted(lm(cbind(Sepal.Length, Sepal.Width) ~., iris)), tolerance = 1e-5)

  expect_equal(`attributes<-`(fHDbetween(data$PCGDP, data[-5]), NULL), `attributes<-`(fitted(lm(PCGDP ~., data)), NULL), tolerance = 1e-5)
  expect_visible(fHDbetween(data$PCGDP, data[-5], fill = TRUE))
  expect_equal(`attributes<-`(fHDbetween(data[5], data[-5])[[1]], NULL), `attributes<-`(fitted(lm(PCGDP ~., data)), NULL), tolerance = 1e-5)
  expect_visible(fHDbetween(data[5], data[-5], fill = TRUE))

  expect_equal(setRownames(qM(fHDbetween(data[5:6], data[-(5:6)]))), setRownames(fitted(lm(cbind(PCGDP, LIFEEX)  ~., data))), tolerance = 1e-5)
  expect_visible(fHDbetween(data[5:6], data[-(5:6)], fill = TRUE))
  expect_visible(fHDbetween(data[5:6], data[-(5:6)], variable.wise = TRUE))

  expect_equal(setRownames(qM(fHDbetween(data[5:7], data[-(5:7)]))), setRownames(fitted(lm(cbind(PCGDP, LIFEEX, ODA)  ~., data))), tolerance = 1e-5)
  expect_visible(fHDbetween(data[5:7], data[-(5:7)], fill = TRUE))
  expect_visible(fHDbetween(data[5:7], data[-(5:7)], variable.wise = TRUE))

  expect_equal(setRownames(qM(fHDbetween(data[5:6], data$ODA))), setRownames(fitted(lm(cbind(PCGDP, LIFEEX)  ~., data[5:7]))), tolerance = 1e-5)
  expect_equal(fHDbetween(data[5:6], data[7], fill = TRUE), fHDbetween(data[5:6], data$ODA, fill = TRUE), tolerance = 1e-5)
  expect_equal(fHDbetween(data[5:6], data[7], variable.wise = TRUE), fHDbetween(data[5:6], data$ODA, variable.wise = TRUE), tolerance = 1e-5)
})

test_that("fHDwithin with multiple variables performs like lm", {
  expect_equal(fHDwithin(iris$Sepal.Length, iris[-1]), `names<-`(resid(lm(Sepal.Length ~., iris)), NULL), tolerance = 1e-5)
  expect_equal(fHDwithin(iris[1], iris[-1])[[1]], `names<-`(resid(lm(Sepal.Length ~., iris)), NULL), tolerance = 1e-5)
  expect_equal(setRownames(qM(fHDwithin(iris[1:2], iris[-(1:2)]))), resid(lm(cbind(Sepal.Length, Sepal.Width) ~., iris)), tolerance = 1e-5)

  expect_equal(`attributes<-`(fHDwithin(data$PCGDP, data[-5]), NULL), `attributes<-`(resid(lm(PCGDP ~., data)), NULL), tolerance = 1e-5)
  expect_visible(fHDwithin(data$PCGDP, data[-5], fill = TRUE))
  expect_equal(`attributes<-`(fHDwithin(data[5], data[-5])[[1]], NULL), `attributes<-`(resid(lm(PCGDP ~., data)), NULL), tolerance = 1e-5)
  expect_visible(fHDwithin(data[5], data[-5], fill = TRUE))

  expect_equal(setRownames(qM(fHDwithin(data[5:6], data[-(5:6)]))), setRownames(resid(lm(cbind(PCGDP, LIFEEX)  ~., data))), tolerance = 1e-5)
  expect_visible(fHDwithin(data[5:6], data[-(5:6)], fill = TRUE))
  expect_visible(fHDwithin(data[5:6], data[-(5:6)], variable.wise = TRUE))

  expect_equal(setRownames(qM(fHDwithin(data[5:7], data[-(5:7)]))), setRownames(resid(lm(cbind(PCGDP, LIFEEX, ODA)  ~., data))), tolerance = 1e-5)
  expect_visible(fHDwithin(data[5:7], data[-(5:7)], fill = TRUE))
  expect_visible(fHDwithin(data[5:7], data[-(5:7)], variable.wise = TRUE))

  expect_equal(setRownames(qM(fHDwithin(data[5:6], data$ODA))), setRownames(resid(lm(cbind(PCGDP, LIFEEX)  ~., data[5:7]))), tolerance = 1e-5)
  expect_equal(fHDwithin(data[5:6], data[7], fill = TRUE), fHDwithin(data[5:6], data$ODA, fill = TRUE), tolerance = 1e-5)
  expect_equal(fHDwithin(data[5:6], data[7], variable.wise = TRUE), fHDwithin(data[5:6], data$ODA, variable.wise = TRUE), tolerance = 1e-5)
})

test_that("fHDbetween produces errors for wrong input", {
  expect_visible(fHDbetween(1:2,1:2))
  expect_error(fHDbetween("a", 1))
  expect_error(fHDbetween(mNAc, f))
  expect_error(fHDbetween(1:2,1:3))
  expect_error(fHDbetween(m,1:31))
  expect_error(fHDbetween(mNA,1:31))
  expect_error(fHDbetween(mtcars,1:31))
  expect_warning(fHDbetween(1:2, 1:2, bla = 1))
  expect_error(fHDbetween(wlddev, list(wlddev$iso3c, wlddev$income[1:10000])))
  expect_visible(fHDbetween(1:2,1:2, na.rm = FALSE))
  expect_error(fHDbetween("a", 1, na.rm = FALSE))
  expect_error(fHDbetween(mNAc, f, na.rm = FALSE))
  expect_error(fHDbetween(1:2,1:3, na.rm = FALSE))
  expect_error(fHDbetween(m,1:31, na.rm = FALSE))
  expect_error(fHDbetween(mNA,1:31, na.rm = FALSE))
  expect_error(fHDbetween(mtcars,1:31, na.rm = FALSE))
  expect_warning(fHDbetween(1:2, 1:2, bla = 1, na.rm = FALSE))
  expect_error(fHDbetween(wlddev, list(wlddev$iso3c, wlddev$income[1:10000]), na.rm = FALSE))
})

test_that("fHDwithin produces errors for wrong input", {
  expect_visible(fHDwithin(1:2,1:2))
  expect_error(fHDwithin("a", 1))
  expect_error(fHDwithin(mNAc, f))
  expect_error(fHDwithin(1:2,1:3))
  expect_error(fHDwithin(m,1:31))
  expect_error(fHDwithin(mNA,1:31))
  expect_error(fHDwithin(mtcars,1:31))
  expect_warning(fHDwithin(1:2, 1:2, bla = 1))
  expect_error(fHDwithin(wlddev, list(wlddev$iso3c, wlddev$income[1:10000])))
  expect_visible(fHDwithin(1:2,1:2, na.rm = FALSE))
  expect_error(fHDwithin("a", 1, na.rm = FALSE))
  expect_error(fHDwithin(mNAc, f, na.rm = FALSE))
  expect_error(fHDwithin(1:2,1:3, na.rm = FALSE))
  expect_error(fHDwithin(m,1:31, na.rm = FALSE))
  expect_error(fHDwithin(mNA,1:31, na.rm = FALSE))
  expect_error(fHDwithin(mtcars,1:31, na.rm = FALSE))
  expect_warning(fHDwithin(1:2, 1:2, bla = 1, na.rm = FALSE))
  expect_error(fHDwithin(wlddev, list(wlddev$iso3c, wlddev$income[1:10000]), na.rm = FALSE))
})

# HDB and HDW
test_that("HDW data.frame method (formula input) performs properly", {
  # simple lm, continuous vars
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb + gear + wt, mtcars))[2:3], tolerance = 1e-3)
  # continuous interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ carb*gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear + wt, mtcars))[2:3], tolerance = 1e-3)
  # continuous 3-way interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ carb*gear*wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear*wt, mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am), mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects + factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am), mtcars))[2:3], tolerance = 1e-3)
  # 3 way factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):factor(am), mtcars))[2:3], tolerance = 1e-3)
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp , W(mtcars, ~ cyl + vs + am, stub = FALSE)))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variable
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variables and factor-continuous interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):gear + factor(am):carb + wt, mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variables and full interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs)*gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs)*gear + factor(am):carb + wt, mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variables and factor-continuous interactions + factor interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs) + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am):carb + wt, mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variables and polynomaial interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, mtcars))[2:3], tolerance = 1e-3)
  # 3-way interaction continuous-factor: error
  # expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):vs:gear + factor(am):carb + wt, stub = FALSE)))[2:3],
  #             coef(lm(mpg ~ hp + disp + factor(cyl):vs:gear + factor(am):carb + wt, mtcars))[2:3])
  # 3-way interaction factor-continuous: error
  # expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
  #              coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):gear + factor(am):carb + wt, mtcars))[2:3])

})

test_that("HDW data.frame method (formula input) with 2-sided formula performs properly", {
  # simple lm, continuous vars
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb + gear + wt, mtcars))[2:3], tolerance = 1e-3)
  # continuous interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ carb*gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear + wt, mtcars))[2:3], tolerance = 1e-3)
  # continuous 3-way interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ carb*gear*wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear*wt, mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am), mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects + factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am), mtcars))[2:3], tolerance = 1e-3)
  # 3 way factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs):factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):factor(am), mtcars))[2:3], tolerance = 1e-3)
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs):factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp , W(mtcars, ~ cyl + vs + am, stub = FALSE)))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variable
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variables and factor-continuous interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):gear + factor(am):carb + wt, mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variables and full interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs)*gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs)*gear + factor(am):carb + wt, mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variables and factor-continuous interactions + factor interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am):carb + wt, mtcars))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variables and polynomaial interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, mtcars))[2:3], tolerance = 1e-3)
  # 3-way interaction continuous-factor: error
  # expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):vs:gear + factor(am):carb + wt, stub = FALSE)))[2:3],
  #             coef(lm(mpg ~ hp + disp + factor(cyl):vs:gear + factor(am):carb + wt, mtcars))[2:3])
  # 3-way interaction factor-continuous: error
  # expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
  #              coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):gear + factor(am):carb + wt, mtcars))[2:3])

})

test_that("HDW data.frame method (formula input) with 2-sided formula and missing values performs properly", {
  # simple lm, continuous vars
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb + gear + wt, mtcNA))[2:3], tolerance = 1e-3)
  # continuous interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ carb*gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear + wt, mtcNA))[2:3], tolerance = 1e-3)
  # continuous 3-way interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ carb*gear*wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear*wt, mtcNA))[2:3], tolerance = 1e-3)
  # HD fixed effects
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am), mtcNA))[2:3], tolerance = 1e-3)
  # HD fixed effects + factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am), mtcNA))[2:3], tolerance = 1e-3)
  # 3 way factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs):factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):factor(am), mtcNA))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variable
  # expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, stub = FALSE)))[2:3],
  #              coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, mtcNA))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variables and factor-continuous interactions : Somestimes test fails, I don't know why (maybe demeanlist numeric problem)
  # expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
  #              coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):gear + factor(am):carb + wt, mtcNA))[2:3], tolerance = 1e-3)
  # HD fixed effects and continuous variables and full interactions
  #expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs)*gear + factor(am):carb + wt, stub = FALSE)))[2:3],
   #            coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs)*gear + factor(am):carb + wt, mtcNA))[2:3], tolerance = 1) # faile R CMD Arch i386 (32 Bit)
  # HD fixed effects and continuous variables and factor-continuous interactions + factor interactions
  # expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am):carb + wt, stub = FALSE)))[2:3],
  #             coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am):carb + wt, mtcNA))[2:3], tolerance = 1e-2)
  # 3-way interaction continuous-factor: error
  # expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, ~ factor(cyl):vs:gear + factor(am):carb + wt, stub = FALSE)))[2:3],
  #             coef(lm(mpg ~ hp + disp + factor(cyl):vs:gear + factor(am):carb + wt, mtcNA))[2:3])
  # 3-way interaction factor-continuous: error
  # expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, ~ factor(cyl):factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
  #              coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):gear + factor(am):carb + wt, mtcNA))[2:3])

})

test_that("HDB data.frame method (formula input) throw errors", {
  expect_error(HDB(mtcars, ~ cyl + vs1))
  expect_error(HDB(mtcars, mpg1 + hp ~ cyl + vs))
  expect_error(HDB(mtcars, ~ cyl + vs, cols = 13))
  expect_error(HDB(mtcars, ~ cyl + vs, cols = "mpg2"))
})

test_that("HDW data.frame method (formula input) throw errors", {
  expect_error(HDW(mtcars, ~ cyl + vs1))
  expect_error(HDW(mtcars, mpg1 + hp ~ cyl + vs))
  expect_error(HDW(mtcars, ~ cyl + vs, cols = 13))
  expect_error(HDW(mtcars, ~ cyl + vs, cols = "mpg2"))
})

options(warn = 1)
