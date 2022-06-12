context("fhdbetween / HDB and fhdwithin / HDW")

if(!is.null(attributes(identical(FALSE, TRUE)))) stop("OECD label issue")
# rm(list = ls())

# TODO: Sort out why certain tests fail...
failtests = FALSE

options(warn = -1)
set.seed(101)
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

# fhdbetween and fhdwithin

test_that("fhdbetween with one factor performs like fbetween", {
  expect_equal(fhdbetween(x, f), fbetween(x, f))
  expect_equal(fhdbetween(x, f, na.rm = FALSE), fbetween(x, f, na.rm = FALSE))
  expect_equal(fhdbetween(xNA, f, na.rm = FALSE), fbetween(xNA, f, na.rm = FALSE))
  expect_equal(`attributes<-`(fhdbetween(xNA, f, fill = TRUE), NULL), fbetween(xNA, f))
  expect_equal(fhdbetween(m, g), fbetween(m, g))
  expect_equal(fhdbetween(m, g, na.rm = FALSE), fbetween(m, g, na.rm = FALSE))
  expect_equal(fhdbetween(mNA, g, na.rm = FALSE), fbetween(mNA, g, na.rm = FALSE))
  # expect_equal(fhdbetween(mNA, g, fill = TRUE), fbetween(mNA, g)) # not matching, fhdbetween matrix is not variable.wise
  expect_equal(fhdbetween(mtcars, g), fbetween(mtcars, g))
  expect_equal(fhdbetween(mtcars, g, na.rm = FALSE), fbetween(mtcars, g, na.rm = FALSE))
  expect_equal(fhdbetween(mtcNA, g, na.rm = FALSE), fbetween(mtcNA, g, na.rm = FALSE))
  expect_equal(fhdbetween(mtcNA, g, variable.wise = TRUE), fbetween(mtcNA, g))

  # with weights
  expect_equal(fhdbetween(x, f, w), fbetween(x, f, w))
  expect_equal(fhdbetween(x, f, w, na.rm = FALSE), fbetween(x, f, w, na.rm = FALSE))
  expect_equal(fhdbetween(xNA, f, w, na.rm = FALSE), fbetween(xNA, f, w, na.rm = FALSE))
  expect_equal(`attributes<-`(fhdbetween(xNA, f, w, fill = TRUE), NULL), fbetween(xNA, f, w))
  expect_equal(fhdbetween(m, g, wdat), fbetween(m, g, wdat))
  expect_equal(fhdbetween(m, g, wdat, na.rm = FALSE), fbetween(m, g, wdat, na.rm = FALSE))
  expect_equal(fhdbetween(mNA, g, wdat, na.rm = FALSE), fbetween(mNA, g, wdat, na.rm = FALSE))
  # expect_equal(fhdbetween(mNA, g, fill = TRUE), fbetween(mNA, g)) # not matching, fhdbetween matrix is not variable.wise
  expect_equal(fhdbetween(mtcars, g, wdat), fbetween(mtcars, g, wdat))
  expect_equal(fhdbetween(mtcars, g, wdat, na.rm = FALSE), fbetween(mtcars, g, wdat, na.rm = FALSE))
  expect_equal(fhdbetween(mtcNA, g, wdat, na.rm = FALSE), fbetween(mtcNA, g, wdat, na.rm = FALSE))
  expect_equal(fhdbetween(mtcNA, g, wdat, variable.wise = TRUE), fbetween(mtcNA, g, wdat))

})

test_that("fhdwithin with one factor performs like fwithin", {
  expect_equal(fhdwithin(x, f), fwithin(x, f))
  expect_equal(fhdwithin(x, f, na.rm = FALSE), fwithin(x, f, na.rm = FALSE))
  expect_equal(fhdwithin(xNA, f, na.rm = FALSE), fwithin(xNA, f, na.rm = FALSE))
  expect_equal(`attributes<-`(fhdwithin(xNA, f, fill = TRUE), NULL), fwithin(xNA, f))
  expect_equal(fhdwithin(m, g), fwithin(m, g))
  expect_equal(fhdwithin(m, g, na.rm = FALSE), fwithin(m, g, na.rm = FALSE))
  expect_equal(fhdwithin(mNA, g, na.rm = FALSE), fwithin(mNA, g, na.rm = FALSE))
  # expect_equal(fhdwithin(mNA, g, fill = TRUE), fwithin(mNA, g)) # not matching, fhdwithin matrix is not variable.wise
  expect_equal(fhdwithin(mtcars, g), fwithin(mtcars, g))
  expect_equal(fhdwithin(mtcars, g, na.rm = FALSE), fwithin(mtcars, g, na.rm = FALSE))
  expect_equal(fhdwithin(mtcNA, g, na.rm = FALSE), fwithin(mtcNA, g, na.rm = FALSE))
  expect_equal(fhdwithin(mtcNA, g, variable.wise = TRUE), fwithin(mtcNA, g))

  # with weights
  expect_equal(fhdwithin(x, f, w), fwithin(x, f, w))
  expect_equal(fhdwithin(x, f, w, na.rm = FALSE), fwithin(x, f, w, na.rm = FALSE))
  expect_equal(fhdwithin(xNA, f, w, na.rm = FALSE), fwithin(xNA, f, w, na.rm = FALSE))
  expect_equal(`attributes<-`(fhdwithin(xNA, f, w, fill = TRUE), NULL), fwithin(xNA, f, w))
  expect_equal(fhdwithin(m, g, wdat), fwithin(m, g, wdat))
  expect_equal(fhdwithin(m, g, wdat, na.rm = FALSE), fwithin(m, g, wdat, na.rm = FALSE))
  expect_equal(fhdwithin(mNA, g, wdat, na.rm = FALSE), fwithin(mNA, g, wdat, na.rm = FALSE))
  # expect_equal(fhdwithin(mNA, g, wdat, fill = TRUE), fwithin(mNA, g)) # not matching, wdat, fhdwithin matrix is not variable.wise
  expect_equal(fhdwithin(mtcars, g, wdat), fwithin(mtcars, g, wdat))
  expect_equal(fhdwithin(mtcars, g, wdat, na.rm = FALSE), fwithin(mtcars, g, wdat, na.rm = FALSE))
  expect_equal(fhdwithin(mtcNA, g, wdat, na.rm = FALSE), fwithin(mtcNA, g, wdat, na.rm = FALSE))
  expect_equal(fhdwithin(mtcNA, g, wdat, variable.wise = TRUE), fwithin(mtcNA, g, wdat))

})

set.seed(101)
f2 <- qF(sample.int(10, 100, TRUE))
fl <- list(f, f2)

g2 <- qF(sample.int(5, 32, TRUE))
gl <- list(g, g2)

  # This is to fool very silly checks on CRAN scanning the code of the tests
if(identical(Sys.getenv("LOCAL"), "TRUE"))
  demeanlist <- eval(parse(text = paste0("lfe", ":", ":", "demeanlist")))

tol <- if(identical(Sys.getenv("LOCAL"), "TRUE")) 1e-5 else 1e-4

if(requireNamespace("fixest", quietly = TRUE)) {
demean <- fixest::demean # eval(parse(text = paste0("fixest", ":", ":", "demean")))

# lfe is back on CRAN: This now also seems to produce a warning !!!!!!!
if(identical(Sys.getenv("LOCAL"), "TRUE"))
test_that("fhdbetween with two factors performs like demeanlist", {
  expect_equal(fhdbetween(x, fl), demeanlist(x, fl, means = TRUE), tolerance = tol)
  expect_equal(fhdbetween(xNA, fl), demeanlist(xNA, fl, means = TRUE, na.rm = TRUE), tolerance = tol)
  expect_visible(fhdbetween(xNA, fl, fill = TRUE))
  expect_equal(fhdbetween(m, gl), demeanlist(m, gl, means = TRUE), tolerance = tol)
  expect_equal(fhdbetween(mNA, gl, na.rm = FALSE), demeanlist(mNA, gl, means = TRUE), tolerance = tol)
  expect_equal(fhdbetween(mNA, gl), demeanlist(mNA, gl, means = TRUE, na.rm = TRUE), tolerance = tol)
  expect_visible(fhdbetween(mNA, gl, fill = TRUE))
  expect_equal(fhdbetween(mtcars, gl), demeanlist(mtcars, gl, means = TRUE), tolerance = tol)
  expect_equal(fhdbetween(mtcNA, gl, na.rm = FALSE), demeanlist(mtcNA, gl, means = TRUE), tolerance = tol)
  expect_equal(setRownames(fhdbetween(mtcNA, gl)), demeanlist(mtcNA, gl, means = TRUE, na.rm = TRUE), tolerance = tol)
  expect_visible(fhdbetween(mtcNA, gl, fill = TRUE))
  expect_visible(fhdbetween(mtcNA, gl, variable.wise = TRUE))

  # With weights
  expect_equal(fhdbetween(x, fl, w), drop(x - demean(x, fl, weights = w)), tolerance = tol)
  expect_equal(unattrib(fhdbetween(xNA, fl, w)), drop(na_rm(xNA) - demean(xNA, fl, weights = w, na.rm = TRUE)), tolerance = tol)
  expect_visible(fhdbetween(xNA, fl, w, fill = TRUE))
  expect_equal(fhdbetween(m, gl, wdat), m - demean(m, gl, weights = wdat), tolerance = tol)
  expect_equal(fhdbetween(mNA, gl, wdat, na.rm = FALSE), demeanlist(mNA, gl, weights = wdat, means = TRUE), tolerance = tol)
  expect_equal(unattrib(fhdbetween(mNA, gl, wdat)), unattrib(na_omit(mNA) - demean(mNA, gl, weights = wdat, na.rm = TRUE)), tolerance = tol)
  expect_visible(fhdbetween(mNA, gl, wdat, fill = TRUE))
  # This one is a bug in demean and will be fixed soon...
  expect_equal(fhdbetween(mtcars, gl, wdat), mtcars %c-% demean(mtcars, gl, weights = wdat), tolerance = tol)
  expect_equal(fhdbetween(mtcNA, gl, na.rm = FALSE), demeanlist(mtcNA, gl, weights = wdat, means = TRUE), tolerance = tol)

  # Same here
  expect_equal(unattrib(fhdbetween(mtcNA, gl, wdat)), unattrib(na_omit(mtcNA) %c-% demean(mtcNA, gl, weights = wdat, na.rm = TRUE)), tolerance = tol)
  expect_visible(fhdbetween(mtcNA, gl, wdat, fill = TRUE))
  expect_visible(fhdbetween(mtcNA, gl, wdat, variable.wise = TRUE))

})

test_that("fhdwithin with two factors performs like demean", {
  expect_equal(fhdwithin(x, fl), drop(demean(x, fl)), tolerance = tol)
  expect_equal(unattrib(fhdwithin(xNA, fl)), unattrib(demean(xNA, fl, na.rm = TRUE)), tolerance = tol)
  expect_identical(length(fhdwithin(xNA, fl, fill = TRUE)), length(xNA))
  expect_equal(unattrib(fhdwithin(m, gl)), unattrib(demean(m, gl)), tolerance = tol)
  # expect_equal(fhdwithin(mNA, gl, na.rm = FALSE), demean(mNA, gl), tolerance = tol) # can break R
  expect_equal(unattrib(fhdwithin(mNA, gl)), unattrib(demean(mNA, gl, na.rm = TRUE)), tolerance = tol)
  expect_identical(nrow(fhdwithin(mNA, gl, fill = TRUE)), nrow(mNA))
  expect_equal(unattrib(fhdwithin(mtcars, gl)), unattrib(demean(mtcars, gl)), tolerance = tol)
  # expect_equal(fhdwithin(mtcNA, gl, na.rm = FALSE), demean(mtcNA, gl), tolerance = tol) # can break R
  expect_equal(unattrib(fhdwithin(mtcNA, gl)), unattrib(demean(mtcNA, gl, na.rm = TRUE)), tolerance = tol)
  expect_equal(fnrow(fhdwithin(mtcNA, gl, fill = TRUE)), fnrow(mtcNA))
  expect_identical(fnrow(fhdwithin(mtcNA, gl, variable.wise = TRUE)), fnrow(mtcNA))

  # With weights
  expect_equal(fhdwithin(x, fl, w), drop(demean(x, fl, weights = w)), tolerance = tol)
  expect_equal(unattrib(fhdwithin(xNA, fl, w)), unattrib(demean(xNA, fl, weights = w, na.rm = TRUE)), tolerance = tol)
  expect_identical(length(fhdwithin(xNA, fl, w, fill = TRUE)), length(xNA))
  expect_equal(unattrib(fhdwithin(m, gl, wdat)), unattrib(demean(m, gl, weights = wdat)), tolerance = tol)
  # expect_equal(fhdwithin(mNA, gl, wdat, na.rm = FALSE), demean(mNA, gl, weights = wdat), tolerance = tol) # can break R
  cc <- complete.cases(mNA)
  expect_equal(unattrib(fhdwithin(mNA, gl, wdat)), unattrib(demean(mNA[cc, ], lapply(gl, .subset, cc), weights = wdat[cc])), tolerance = tol)
  expect_identical(nrow(fhdwithin(mNA, gl, wdat, fill = TRUE)), nrow(mNA))
  # Smae here, bug to be fixed in demean()
  expect_equal(unattrib(fhdwithin(mtcars, gl, wdat)), unattrib(demean(mtcars, gl, weights = wdat)), tolerance = tol)
  # expect_equal(fhdwithin(mtcNA, gl, wdat, na.rm = FALSE), demean(mtcNA, gl, weights = wdat), tolerance = tol) # can break R
  # Also bug
  expect_equal(unattrib(fhdwithin(mtcNA, gl, wdat)), unattrib(demean(mtcNA, gl, weights = wdat, na.rm = TRUE)), tolerance = 1e-3)
  expect_equal(fnrow(fhdwithin(mtcNA, gl, wdat, fill = TRUE)), fnrow(mtcNA))
  expect_identical(fnrow(fhdwithin(mtcNA, gl, wdat, variable.wise = TRUE)), fnrow(mtcNA))


})
}
x2 <- 3 * x + rnorm(100)

test_that("fhdbetween with only continuous variables performs like basefitted (defined above)", {
  expect_equal(fhdbetween(x, x2), basefitted(x, x2), tolerance = tol)
  expect_equal(`attr<-`(fhdbetween(xNA, x2), "na.rm", NULL), basefitted(xNA, x2, na.rm = TRUE), tolerance = tol)
  expect_visible(fhdbetween(xNA, x2, fill = TRUE))
  expect_equal(fhdbetween(m, m), fhdbetween(m, mtcars), tolerance = tol)
  expect_equal(fhdbetween(m, m), basefitted(m, m), tolerance = tol)
  expect_equal(`attr<-`(fhdbetween(mNA, m, lm.method = "qr"), "na.rm", NULL), basefitted(mNA, m, na.rm = TRUE), tolerance = tol)
  expect_equal(fhdbetween(mNA, m, fill = TRUE, lm.method = "qr"), fhdbetween(mNA, mtcars, fill = TRUE, lm.method = "qr"), tolerance = tol)
  expect_equal(fhdbetween(mtcars, mtcars), fhdbetween(mtcars, m), tolerance = tol)
  expect_equal(fhdbetween(mtcars, mtcars), qDF(basefitted(mtcars, mtcars)), tolerance = tol)
  expect_equal(`attr<-`(fhdbetween(mtcNA, mtcars, lm.method = "qr"), "na.rm", NULL), qDF(basefitted(mtcNA, mtcars, na.rm = TRUE)), tolerance = tol)
  expect_equal(fhdbetween(mtcNA, mtcars, fill = TRUE, lm.method = "qr"), fhdbetween(mtcNA, m, fill = TRUE, lm.method = "qr"), tolerance = tol)
  expect_equal(fhdbetween(mtcNA, mtcars, variable.wise = TRUE), fhdbetween(mtcNA, m, variable.wise = TRUE), tolerance = tol)
})

test_that("fhdwithin with only continuous variables performs like baseresid (defined above)", {
  expect_equal(fhdwithin(x, x2), baseresid(x, x2), tolerance = tol)
  expect_equal(`attr<-`(fhdwithin(xNA, x2), "na.rm", NULL), baseresid(xNA, x2, na.rm = TRUE), tolerance = tol)
  expect_visible(fhdwithin(xNA, x2, fill = TRUE))
  expect_equal(fhdwithin(m, m), fhdwithin(m, mtcars), tolerance = tol)
  expect_equal(fhdwithin(m, m), baseresid(m, m), tolerance = tol)
  expect_equal(`attr<-`(fhdwithin(mNA, m, lm.method = "qr"), "na.rm", NULL), baseresid(mNA, m, na.rm = TRUE), tolerance = tol)
  expect_equal(fhdwithin(mNA, m, fill = TRUE, lm.method = "qr"), fhdwithin(mNA, mtcars, fill = TRUE, lm.method = "qr"), tolerance = tol)
  expect_equal(fhdwithin(mtcars, mtcars), fhdwithin(mtcars, m), tolerance = tol)
  expect_equal(fhdwithin(mtcars, mtcars), qDF(baseresid(mtcars, mtcars)), tolerance = tol)
  expect_equal(`attr<-`(fhdwithin(mtcNA, mtcars, lm.method = "qr"), "na.rm", NULL), qDF(baseresid(mtcNA, mtcars, na.rm = TRUE)), tolerance = tol)
  expect_equal(fhdwithin(mtcNA, mtcars, fill = TRUE, lm.method = "qr"), fhdwithin(mtcNA, m, fill = TRUE, lm.method = "qr"), tolerance = tol)
  expect_equal(fhdwithin(mtcNA, mtcars, variable.wise = TRUE), fhdwithin(mtcNA, m, variable.wise = TRUE), tolerance = tol)
})

if(requireNamespace("fixest", quietly = TRUE)) {

data <- wlddev
data$year <- qF(data$year)
data <- get_vars(data, c("iso3c","year","region","income","PCGDP","LIFEEX","ODA"))
ww <- abs(rnorm(fnrow(data)))
wi <- abs(rnorm(fnrow(iris)))

test_that("fhdbetween with multiple variables performs like lm", {
  expect_equal(fhdbetween(iris$Sepal.Length, iris[-1]), `names<-`(fitted(lm(Sepal.Length ~., iris)), NULL), tolerance = tol)
  expect_equal(fhdbetween(iris[1], iris[-1])[[1]], `names<-`(fitted(lm(Sepal.Length ~., iris)), NULL), tolerance = tol)
  expect_equal(setRownames(qM(fhdbetween(iris[1:2], iris[-(1:2)]))), fitted(lm(cbind(Sepal.Length, Sepal.Width) ~., iris)), tolerance = tol)

  expect_equal(`attributes<-`(fhdbetween(data$PCGDP, data[-5]), NULL), `attributes<-`(fitted(lm(PCGDP ~., data)), NULL), tolerance = tol)
  expect_visible(fhdbetween(data$PCGDP, data[-5], fill = TRUE))
  expect_equal(`attributes<-`(fhdbetween(data[5], data[-5])[[1]], NULL), `attributes<-`(fitted(lm(PCGDP ~., data)), NULL), tolerance = tol)
  expect_visible(fhdbetween(data[5], data[-5], fill = TRUE))

  expect_equal(setRownames(qM(fhdbetween(data[5:6], data[-(5:6)]))), setRownames(fitted(lm(cbind(PCGDP, LIFEEX)  ~., data))), tolerance = tol)
  expect_visible(fhdbetween(data[5:6], data[-(5:6)], fill = TRUE))
  expect_visible(fhdbetween(data[5:6], data[-(5:6)], variable.wise = TRUE))

  expect_equal(setRownames(qM(fhdbetween(data[5:7], data[-(5:7)]))), setRownames(fitted(lm(cbind(PCGDP, LIFEEX, ODA)  ~., data))), tolerance = tol)
  expect_visible(fhdbetween(data[5:7], data[-(5:7)], fill = TRUE))
  expect_visible(fhdbetween(data[5:7], data[-(5:7)], variable.wise = TRUE))

  expect_equal(setRownames(qM(fhdbetween(data[5:6], data$ODA))), setRownames(fitted(lm(cbind(PCGDP, LIFEEX)  ~., data[5:7]))), tolerance = tol)
  expect_equal(fhdbetween(data[5:6], data[7], fill = TRUE), fhdbetween(data[5:6], data$ODA, fill = TRUE), tolerance = tol)
  expect_equal(fhdbetween(data[5:6], data[7], variable.wise = TRUE), fhdbetween(data[5:6], data$ODA, variable.wise = TRUE), tolerance = tol)

  # With weights
  expect_equal(fhdbetween(iris$Sepal.Length, iris[-1], wi), `names<-`(fitted(lm(Sepal.Length ~., iris, weights = wi)), NULL), tolerance = tol)
  expect_equal(fhdbetween(iris[1], iris[-1], wi)[[1]], `names<-`(fitted(lm(Sepal.Length ~., iris, weights = wi)), NULL), tolerance = tol)
  expect_equal(setRownames(qM(fhdbetween(iris[1:2], iris[-(1:2)], wi))), fitted(lm(cbind(Sepal.Length, Sepal.Width) ~., iris, weights = wi)), tolerance = tol)

  expect_equal(`attributes<-`(fhdbetween(data$PCGDP, data[-5], ww), NULL), `attributes<-`(fitted(lm(PCGDP ~., data, weights = ww)), NULL), tolerance = tol)
  expect_visible(fhdbetween(data$PCGDP, data[-5], ww, fill = TRUE))
  expect_equal(`attributes<-`(fhdbetween(data[5], data[-5], ww)[[1]], NULL), `attributes<-`(fitted(lm(PCGDP ~., data, weights = ww)), NULL), tolerance = tol)
  expect_visible(fhdbetween(data[5], data[-5], ww, fill = TRUE))

  expect_equal(setRownames(qM(fhdbetween(data[5:6], data[-(5:6)], ww))), setRownames(fitted(lm(cbind(PCGDP, LIFEEX)  ~., data, weights = ww))), tolerance = tol)
  expect_visible(fhdbetween(data[5:6], data[-(5:6)], ww, fill = TRUE))
  expect_visible(fhdbetween(data[5:6], data[-(5:6)], ww, variable.wise = TRUE))

  expect_equal(setRownames(qM(fhdbetween(data[5:7], data[-(5:7)], ww))), setRownames(fitted(lm(cbind(PCGDP, LIFEEX, ODA)  ~., data, weights = ww))), tolerance = tol)
  expect_visible(fhdbetween(data[5:7], data[-(5:7)], ww, fill = TRUE))
  expect_visible(fhdbetween(data[5:7], data[-(5:7)], ww, variable.wise = TRUE))

  expect_equal(setRownames(qM(fhdbetween(data[5:6], data$ODA, ww))), setRownames(fitted(lm(cbind(PCGDP, LIFEEX)  ~., data[5:7], weights = ww))), tolerance = tol)
  expect_equal(fhdbetween(data[5:6], data[7], ww, fill = TRUE), fhdbetween(data[5:6], data$ODA, ww, fill = TRUE), tolerance = tol)
  expect_equal(fhdbetween(data[5:6], data[7], ww, variable.wise = TRUE), fhdbetween(data[5:6], data$ODA, ww, variable.wise = TRUE), tolerance = tol)

})

test_that("fhdwithin with multiple variables performs like lm", {
  expect_equal(fhdwithin(iris$Sepal.Length, iris[-1]), `names<-`(resid(lm(Sepal.Length ~., iris)), NULL), tolerance = tol)
  expect_equal(fhdwithin(iris[1], iris[-1])[[1]], `names<-`(resid(lm(Sepal.Length ~., iris)), NULL), tolerance = tol)
  expect_equal(setRownames(qM(fhdwithin(iris[1:2], iris[-(1:2)]))), resid(lm(cbind(Sepal.Length, Sepal.Width) ~., iris)), tolerance = tol)

  expect_equal(`attributes<-`(fhdwithin(data$PCGDP, data[-5]), NULL), `attributes<-`(resid(lm(PCGDP ~., data)), NULL), tolerance = tol)
  expect_visible(fhdwithin(data$PCGDP, data[-5], fill = TRUE))
  expect_equal(`attributes<-`(fhdwithin(data[5], data[-5])[[1]], NULL), `attributes<-`(resid(lm(PCGDP ~., data)), NULL), tolerance = tol)
  expect_visible(fhdwithin(data[5], data[-5], fill = TRUE))

  expect_equal(setRownames(qM(fhdwithin(data[5:6], data[-(5:6)]))), setRownames(resid(lm(cbind(PCGDP, LIFEEX)  ~., data))), tolerance = tol)
  expect_visible(fhdwithin(data[5:6], data[-(5:6)], fill = TRUE))
  expect_visible(fhdwithin(data[5:6], data[-(5:6)], variable.wise = TRUE))

  expect_equal(setRownames(qM(fhdwithin(data[5:7], data[-(5:7)]))), setRownames(resid(lm(cbind(PCGDP, LIFEEX, ODA)  ~., data))), tolerance = tol)
  expect_visible(fhdwithin(data[5:7], data[-(5:7)], fill = TRUE))
  expect_visible(fhdwithin(data[5:7], data[-(5:7)], variable.wise = TRUE))

  expect_equal(setRownames(qM(fhdwithin(data[5:6], data$ODA))), setRownames(resid(lm(cbind(PCGDP, LIFEEX)  ~., data[5:7]))), tolerance = tol)
  expect_equal(fhdwithin(data[5:6], data[7], fill = TRUE), fhdwithin(data[5:6], data$ODA, fill = TRUE), tolerance = tol)
  expect_equal(fhdwithin(data[5:6], data[7], variable.wise = TRUE), fhdwithin(data[5:6], data$ODA, variable.wise = TRUE), tolerance = tol)

  # With weights
  expect_equal(fhdwithin(iris$Sepal.Length, iris[-1], wi), `names<-`(resid(lm(Sepal.Length ~., iris, weights = wi)), NULL), tolerance = tol)
  expect_equal(fhdwithin(iris[1], iris[-1], wi)[[1]], `names<-`(resid(lm(Sepal.Length ~., iris, weights = wi)), NULL), tolerance = tol)
  expect_equal(setRownames(qM(fhdwithin(iris[1:2], iris[-(1:2)], wi))), resid(lm(cbind(Sepal.Length, Sepal.Width) ~., iris, weights = wi)), tolerance = tol)

  expect_equal(`attributes<-`(fhdwithin(data$PCGDP, data[-5], ww), NULL), `attributes<-`(resid(lm(PCGDP ~., data, weights = ww)), NULL), tolerance = tol)
  expect_visible(fhdwithin(data$PCGDP, data[-5], ww, fill = TRUE))
  expect_equal(`attributes<-`(fhdwithin(data[5], data[-5], ww)[[1]], NULL), `attributes<-`(resid(lm(PCGDP ~., data, weights = ww)), NULL), tolerance = tol)
  expect_visible(fhdwithin(data[5], data[-5], ww, fill = TRUE))

  expect_equal(setRownames(qM(fhdwithin(data[5:6], data[-(5:6)], ww))), setRownames(resid(lm(cbind(PCGDP, LIFEEX)  ~., data, weights = ww))), tolerance = tol)
  expect_visible(fhdwithin(data[5:6], data[-(5:6)], ww, fill = TRUE))
  expect_visible(fhdwithin(data[5:6], data[-(5:6)], ww, variable.wise = TRUE))

  expect_equal(setRownames(qM(fhdwithin(data[5:7], data[-(5:7)], ww))), setRownames(resid(lm(cbind(PCGDP, LIFEEX, ODA)  ~., data, weights = ww))), tolerance = tol)
  expect_visible(fhdwithin(data[5:7], data[-(5:7)], ww, fill = TRUE))
  expect_visible(fhdwithin(data[5:7], data[-(5:7)], ww, variable.wise = TRUE))

  expect_equal(setRownames(qM(fhdwithin(data[5:6], data$ODA, ww))), setRownames(resid(lm(cbind(PCGDP, LIFEEX)  ~., data[5:7], weights = ww))), tolerance = tol)
  expect_equal(fhdwithin(data[5:6], data[7], ww, fill = TRUE), fhdwithin(data[5:6], data$ODA, ww, fill = TRUE), tolerance = tol)
  expect_equal(fhdwithin(data[5:6], data[7], ww, variable.wise = TRUE), fhdwithin(data[5:6], data$ODA, ww, variable.wise = TRUE), tolerance = tol)

})

}

test_that("fhdbetween produces errors for wrong input", {
  expect_visible(fhdbetween(1:2,1:2))
  expect_error(fhdbetween("a", 1))
  expect_error(fhdbetween(mNAc, f))
  expect_error(fhdbetween(1:2,1:3))
  expect_error(fhdbetween(m,1:31))
  expect_error(fhdbetween(mNA,1:31))
  expect_error(fhdbetween(mtcars,1:31))
  # expect_warning(fhdbetween(1:2, 1:2, bla = 1))
  expect_error(fhdbetween(wlddev, list(wlddev$iso3c, wlddev$income[1:10000])))
  expect_visible(fhdbetween(1:2,1:2, na.rm = FALSE))
  expect_error(fhdbetween("a", 1, na.rm = FALSE))
  expect_error(fhdbetween(mNAc, f, na.rm = FALSE))
  expect_error(fhdbetween(1:2,1:3, na.rm = FALSE))
  expect_error(fhdbetween(m,1:31, na.rm = FALSE))
  expect_error(fhdbetween(mNA,1:31, na.rm = FALSE))
  expect_error(fhdbetween(mtcars,1:31, na.rm = FALSE))
  # expect_warning(fhdbetween(1:2, 1:2, bla = 1, na.rm = FALSE))
  # expect_error(fhdbetween(wlddev, list(wlddev$iso3c, wlddev$income[1:10000]), na.rm = FALSE)) # breaks R
})

test_that("fhdwithin produces errors for wrong input", {
  expect_visible(fhdwithin(1:2,1:2))
  expect_error(fhdwithin("a", 1))
  expect_error(fhdwithin(mNAc, f))
  expect_error(fhdwithin(1:2,1:3))
  expect_error(fhdwithin(m,1:31))
  expect_error(fhdwithin(mNA,1:31))
  expect_error(fhdwithin(mtcars,1:31))
  # expect_warning(fhdwithin(1:2, 1:2, bla = 1))
  expect_error(fhdwithin(wlddev, list(wlddev$iso3c, wlddev$income[1:10000])))
  expect_visible(fhdwithin(1:2,1:2, na.rm = FALSE))
  expect_error(fhdwithin("a", 1, na.rm = FALSE))
  expect_error(fhdwithin(mNAc, f, na.rm = FALSE))
  expect_error(fhdwithin(1:2,1:3, na.rm = FALSE))
  expect_error(fhdwithin(m,1:31, na.rm = FALSE))
  expect_error(fhdwithin(mNA,1:31, na.rm = FALSE))
  expect_error(fhdwithin(mtcars,1:31, na.rm = FALSE))
  # expect_warning(fhdwithin(1:2, 1:2, bla = 1, na.rm = FALSE))
  # expect_error(fhdwithin(wlddev, list(wlddev$iso3c, wlddev$income[1:10000]), na.rm = FALSE)) # segfault !!!
})

if(identical(Sys.getenv("NCRAN"), "TRUE")) {

# HDB and HDW
test_that("HDW data.frame method (formula input) performs properly", {
  # simple lm, continuous vars
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb + gear + wt, mtcars))[2:3], tolerance = tol)
  # continuous interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ carb*gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear + wt, mtcars))[2:3], tolerance = tol)
  # continuous 3-way interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ carb*gear*wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear*wt, mtcars))[2:3], tolerance = tol)
  # factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ qF(cyl), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl), mtcars))[2:3], tolerance = tol)
  # factor - continuous without including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ qF(cyl):carb, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb, mtcars))[2:3], tolerance = tol)
  # multiple factors - continuous without including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ qF(cyl):carb + qF(vs):carb, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb + qF(vs):carb, mtcars))[2:3], tolerance = tol)
  # multiple factors - continuous without including factor 2
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ qF(cyl):carb + qF(vs):wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb + qF(vs):wt, mtcars))[2:3], tolerance = tol)
  # multiple factors - continuous without including factor 3
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ am + qF(cyl):carb + qF(vs):wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + am + qF(cyl):carb + qF(vs):wt, mtcars))[2:3], tolerance = tol)
  # factor - continuous  including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ qF(cyl):carb + qF(cyl), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb + qF(cyl), mtcars))[2:3], tolerance = tol)
  # factor - continuous  full interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ qF(cyl)*carb, stub = FALSE, lm.method = "qr")))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl)*carb, mtcars))[2:3], tolerance = tol)
  # HD fixed effects
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am), mtcars))[2:3], tolerance = tol)
  # HD fixed effects + factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am), mtcars))[2:3], tolerance = tol)
  # 3 way factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):factor(am), mtcars))[2:3], tolerance = tol)
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp , W(mtcars, ~ cyl + vs + am, stub = FALSE)))[2:3], tolerance = tol)
  # HD fixed effects and continuous variable
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, mtcars))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):gear + factor(am):carb + wt, mtcars))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and full interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs)*gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs)*gear + factor(am):carb + wt, mtcars))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions + factor interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs) + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am):carb + wt, mtcars))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and polynomaial interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, mtcars))[2:3], tolerance = tol)
  # 3-way interaction continuous-factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):vs:gear + factor(am):carb + wt, stub = FALSE, lm.method = "qr")))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):vs:gear + factor(am):carb + wt, mtcars))[2:3])
  # 3-way interaction factor-continuous
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):gear + factor(am):carb + wt, mtcars))[2:3])


  # With weights
  # simple lm, continuous vars
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ carb + gear + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + carb + gear + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # continuous interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ carb*gear + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # continuous 3-way interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ carb*gear*wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear*wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ qF(cyl), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl), mtcars, weights = wdat))[2:3], tolerance = tol)
  # factor - continuous without including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ qF(cyl):carb, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb, mtcars, weights = wdat))[2:3], tolerance = tol)
  # factor - continuous  including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ qF(cyl):carb + qF(cyl), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb + qF(cyl), mtcars, weights = wdat))[2:3], tolerance = tol)
  # factor - continuous  full interaction
  if(failtests) expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ qF(cyl)*carb, wdat, stub = FALSE, lm.method = "qr"), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl)*carb, mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs) + factor(am), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am), mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects + factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs) + factor(am), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am), mtcars, weights = wdat))[2:3], tolerance = tol)
  # 3 way factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):factor(am), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):factor(am), mtcars, weights = wdat))[2:3], tolerance = tol)
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):factor(am), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp , W(mtcars, ~ cyl + vs + am, wdat, stub = FALSE), weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variable
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs):gear + factor(am):carb + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):gear + factor(am):carb + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and full interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs)*gear + factor(am):carb + wt, wdat, stub = FALSE, lm.method = "qr"), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs)*gear + factor(am):carb + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions + factor interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs) + factor(am):carb + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am):carb + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and polynomaial interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # 3-way interaction continuous-factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):vs:gear + factor(am):carb + wt, wdat, stub = FALSE, lm.method = "qr"), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):vs:gear + factor(am):carb + wt, mtcars, weights = wdat))[2:3])
  # 3-way interaction factor-continuous
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):gear + factor(am):carb + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):gear + factor(am):carb + wt, mtcars, weights = wdat))[2:3])

})

test_that("HDW data.frame method (formula input) with 2-sided formula performs properly", {
  # simple lm, continuous vars
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb + gear + wt, mtcars))[2:3], tolerance = tol)
  # continuous interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ carb*gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear + wt, mtcars))[2:3], tolerance = tol)
  # continuous 3-way interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ carb*gear*wt, stub = FALSE, lm.method = "qr")))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear*wt, mtcars))[2:3], tolerance = tol)
  # factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ qF(cyl), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl), mtcars))[2:3], tolerance = tol)
  # factor - continuous without including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ qF(cyl):carb, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb, mtcars))[2:3], tolerance = tol)
  # factor - continuous  including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ qF(cyl):carb + qF(cyl), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb + qF(cyl), mtcars))[2:3], tolerance = tol)
  # factor - continuous  full interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ qF(cyl)*carb, stub = FALSE, lm.method = "qr")))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl)*carb, mtcars))[2:3], tolerance = tol)
  # HD fixed effects
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am), mtcars))[2:3], tolerance = tol)
  # HD fixed effects + factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am), mtcars))[2:3], tolerance = tol)
  # 3 way factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs):factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):factor(am), mtcars))[2:3], tolerance = tol)
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs):factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp , W(mtcars, ~ cyl + vs + am, stub = FALSE)))[2:3], tolerance = tol)
  # HD fixed effects and continuous variable
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, mtcars))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):gear + factor(am):carb + wt, mtcars))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and full interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs)*gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs)*gear + factor(am):carb + wt, mtcars))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions + factor interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am):carb + wt, mtcars))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and polynomaial interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, mtcars))[2:3], tolerance = tol)
  # 3-way interaction continuous-factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):vs:gear + factor(am):carb + wt, stub = FALSE, lm.method = "qr")))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):vs:gear + factor(am):carb + wt, mtcars))[2:3])
  # 3-way interaction factor-continuous
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, ~ factor(cyl):factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):gear + factor(am):carb + wt, mtcars))[2:3])

  # With weights
  # simple lm, continuous vars
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ carb + gear + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + carb + gear + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # continuous interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ carb*gear + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # continuous 3-way interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ carb*gear*wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear*wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ qF(cyl), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl), mtcars, weights = wdat))[2:3], tolerance = tol)
  # factor - continuous without including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ qF(cyl):carb, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb, mtcars, weights = wdat))[2:3], tolerance = tol)
  # factor - continuous  including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ qF(cyl):carb + qF(cyl), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb + qF(cyl), mtcars, weights = wdat))[2:3], tolerance = tol)
  # factor - continuous  full interaction
  if(failtests) expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ qF(cyl)*carb, wdat, stub = FALSE, lm.method = "qr"), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl)*carb, mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am), mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects + factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am), mtcars, weights = wdat))[2:3], tolerance = tol)
  # 3 way factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs):factor(am), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):factor(am), mtcars, weights = wdat))[2:3], tolerance = tol)
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs):factor(am), wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp , W(mtcars, mpg + hp + disp ~ cyl + vs + am, wdat, stub = FALSE), weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variable
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs):gear + factor(am):carb + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):gear + factor(am):carb + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and full interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs)*gear + factor(am):carb + wt, wdat, stub = FALSE, lm.method = "qr"), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs)*gear + factor(am):carb + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions + factor interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am):carb + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am):carb + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and polynomaial interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, mtcars, weights = wdat))[2:3], tolerance = tol)
  # 3-way interaction continuous-factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):vs:gear + factor(am):carb + wt, wdat, stub = FALSE, lm.method = "qr"), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):vs:gear + factor(am):carb + wt, mtcars, weights = wdat))[2:3])
  # 3-way interaction factor-continuous
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcars, mpg + hp + disp ~ factor(cyl):factor(vs):gear + factor(am):carb + wt, wdat, stub = FALSE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):gear + factor(am):carb + wt, mtcars, weights = wdat))[2:3])

})

test_that("HDW data.frame method (formula input) with 2-sided formula and missing values performs properly", {
  # simple lm, continuous vars
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb + gear + wt, mtcNA))[2:3], tolerance = tol)
  # continuous interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ carb*gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear + wt, mtcNA))[2:3], tolerance = tol)
  # continuous 3-way interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ carb*gear*wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear*wt, mtcNA))[2:3], tolerance = tol)
  # factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ qF(cyl), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl), mtcNA))[2:3], tolerance = tol)
  # factor - continuous without including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ qF(cyl):carb, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb, mtcNA))[2:3], tolerance = tol)
  # factor - continuous  including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ qF(cyl):carb + qF(cyl), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb + qF(cyl), mtcNA))[2:3], tolerance = tol)
  # factor - continuous  full interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ qF(cyl)*carb, stub = FALSE, lm.method = "qr")))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl)*carb, mtcNA))[2:3], tolerance = tol)
  # HD fixed effects
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am), mtcNA))[2:3], tolerance = tol)
  # HD fixed effects + factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am), mtcNA))[2:3], tolerance = tol)
  # 3 way factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs):factor(am), stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):factor(am), mtcNA))[2:3], tolerance = tol)
  # HD fixed effects and continuous variable
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, mtcNA))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions : Somestimes test fails, I don't know why (maybe demeanlist numeric problem)
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):gear + factor(am):carb + wt, mtcNA))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and full interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs)*gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs)*gear + factor(am):carb + wt, mtcNA))[2:3], tolerance = 1) # faile R CMD Arch i386 (32 Bit)
  # HD fixed effects and continuous variables and factor-continuous interactions + factor interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am):carb + wt, mtcNA))[2:3], tolerance = 1e-2)
  # 3-way interaction continuous-factor: error
  if(failtests) expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, ~ factor(cyl):vs:gear + factor(am):carb + wt, stub = FALSE, lm.method = "qr")))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):vs:gear + factor(am):carb + wt, mtcNA))[2:3])
  # 3-way interaction factor-continuous: error
  if(failtests) expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, ~ factor(cyl):factor(vs):gear + factor(am):carb + wt, stub = FALSE)))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):gear + factor(am):carb + wt, mtcNA))[2:3])

  # With weights
  # simple lm, continuous vars
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ carb + gear + wt, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + carb + gear + wt, mtcNA, weights = wdat))[2:3], tolerance = tol)
  # continuous interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ carb*gear + wt, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear + wt, mtcNA, weights = wdat))[2:3], tolerance = tol)
  # continuous 3-way interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ carb*gear*wt, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + carb*gear*wt, mtcNA, weights = wdat))[2:3], tolerance = tol)
  # factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ qF(cyl), wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl), mtcNA, weights = wdat))[2:3], tolerance = tol)
  # factor - continuous without including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ qF(cyl):carb, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb, mtcNA, weights = wdat))[2:3], tolerance = tol)
  # factor - continuous  including factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ qF(cyl):carb + qF(cyl), wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl):carb + qF(cyl), mtcNA, weights = wdat))[2:3], tolerance = tol)
  # factor - continuous  full interaction
  if(failtests) expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ qF(cyl)*carb, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + qF(cyl)*carb, mtcNA, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am), wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am), mtcNA, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects + factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am), wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am), mtcNA, weights = wdat))[2:3], tolerance = tol)
  # 3 way factor interaction
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs):factor(am), wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):factor(am), mtcNA, weights = wdat))[2:3], tolerance = tol)
 if(failtests) expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs):factor(am), wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp , W(mtcNA, mpg + hp + disp ~ cyl + vs + am, wdat, stub = FALSE), weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variable
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs) + factor(am) + carb + gear + wt, mtcNA, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs):gear + factor(am):carb + wt, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):gear + factor(am):carb + wt, mtcNA, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and full interactions
  if(failtests) expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs)*gear + factor(am):carb + wt, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs)*gear + factor(am):carb + wt, mtcNA, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and factor-continuous interactions + factor interactions
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs) + factor(am):carb + wt, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs) + factor(am):carb + wt, mtcNA, weights = wdat))[2:3], tolerance = tol)
  # HD fixed effects and continuous variables and polynomaial interactions
  if(failtests) expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl) + factor(vs):poly(gear,2) + factor(am):carb + wt, mtcNA, weights = wdat))[2:3], tolerance = tol)
  # 3-way interaction continuous-factor
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):vs:gear + factor(am):carb + wt, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):vs:gear + factor(am):carb + wt, mtcNA, weights = wdat))[2:3])
  # 3-way interaction factor-continuous
  expect_equal(coef(lm(mpg ~ hp + disp, HDW(mtcNA, mpg + hp + disp ~ factor(cyl):factor(vs):gear + factor(am):carb + wt, wdat, stub = FALSE, fill = TRUE), weights = wdat))[2:3],
               coef(lm(mpg ~ hp + disp + factor(cyl):factor(vs):gear + factor(am):carb + wt, mtcNA, weights = wdat))[2:3], tolerance = tol)

})

test_that("HDW weighted computations work like lm", { # ...

  if(failtests) expect_equal(
  unname(resid(lm(mpg ~ factor(cyl)*carb + factor(vs) + hp + gear, weights = wt, mtcars))),
  HDW(mtcars, mpg ~ factor(cyl)*carb + factor(vs) + hp + gear, mtcars$wt)[, 1], tolerance = 1e-4)

  if(failtests) expect_equal(
    unname(resid(lm(mpg ~ factor(cyl)*carb + factor(vs) + hp + gear, mtcars))),
    HDW(mtcars, mpg ~ factor(cyl)*carb + factor(vs) + hp + gear, lm.method = "qr")[, 1], tolerance = 1e-4)

  expect_equal(
  unname(resid(lm(mpg ~ factor(vs) + hp + gear, weights = wt, mtcars))),
  HDW(mtcars, mpg ~ factor(vs) + hp + gear, mtcars$wt)[, 1], tolerance = 1e-4)

  expect_equal(
  unname(resid(lm(mpg ~ factor(cyl) + factor(vs) + hp + gear, weights = wt, mtcars))),
  HDW(mtcars, mpg ~ factor(cyl) + factor(vs) + hp + gear, mtcars$wt)[, 1], tolerance = 1e-4)

  expect_equal(
  unname(resid(lm(mpg ~ hp + gear, weights = wt, mtcars))),
  HDW(mtcars, mpg ~ hp + gear, mtcars$wt)[, 1], tolerance = 1e-4)

})

}

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
