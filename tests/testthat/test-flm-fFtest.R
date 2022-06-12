context("flm and fFtest")

if(!is.null(attributes(identical(FALSE, TRUE)))) stop("OECD label issue")

y <- mtcars$mpg
x <- qM(mtcars[c("cyl","vs","am","carb","hp")])
w <- mtcars$wt

lmr <- lm(mpg ~ cyl + vs + am + carb + hp, mtcars)
lmw <- lm(mpg ~ cyl + vs + am + carb + hp, weights = wt, mtcars)

NCRAN <- identical(Sys.getenv("NCRAN"), "TRUE")

test_that("flm works as intended", {

  if(NCRAN) for(i in 1:6) expect_equal(drop(flm(y, x, add.icpt = TRUE, method = i)), coef(lmr))
  if(NCRAN) for(i in 1:6) expect_equal(drop(flm(y, x, w, add.icpt = TRUE, method = i)), coef(lmw))
  expect_equal(flm(y, x, method = 1L, return.raw = TRUE), .lm.fit(x, y))
  expect_equal(flm(y, x, method = 2L, return.raw = TRUE), solve(crossprod(x), crossprod(x, y)))
  expect_equal(flm(y, x, method = 3L, return.raw = TRUE), qr.coef(qr(x), y))
  expect_equal(flm(y, x, method = 5L, return.raw = TRUE), cinv(crossprod(x)) %*% crossprod(x, y))
  if(NCRAN) {
    # This is to fool very silly checks on CRAN scanning the code of the tests
    afmlp <- eval(parse(text = paste0("RcppArmadillo", ":", ":", "fastLmPure")))
    efmlp <- eval(parse(text = paste0("RcppEigen", ":", ":", "fastLmPure")))

    expect_equal(flm(y, x, method = 4L, return.raw = TRUE), afmlp(x, y))
    expect_equal(flm(y, x, method = 6L, return.raw = TRUE), efmlp(x, y, 3L))
  }
  if(NCRAN) for(i in 1:6) expect_visible(flm(y, x, w, method = i, return.raw = TRUE))
  ym <- cbind(y, y)
  for(i in c(1:3, 5L)) expect_visible(flm(ym, x, w, method = i))

  expect_error(flm(y[-1L], x, w))
  expect_error(flm(y, x, w[-1L]))
  expect_error(flm(y, x[-1L, ], w))

})


test_that("fFtest works as intended", {

  r <- fFtest(iris$Sepal.Length, gv(iris, -1L))
  rlm <- summary(lm(Sepal.Length ~., iris))
  expect_equal(unattrib(r)[1:4], unattrib(c(rlm$r.squared, rlm$fstatistic[c(2:3, 1L)])))
  # Same with weights:
  w <- abs(rnorm(fnrow(iris)))
  r <- fFtest(iris$Sepal.Length, gv(iris, -1L), w = w)
  rlm <- summary(lm(Sepal.Length ~., weights = w, iris))
  expect_equal(unattrib(r)[1:4], unattrib(c(rlm$r.squared, rlm$fstatistic[c(2:3, 1L)])))

  # Repeat with missing values
  set.seed(101)
  iris <- na_insert(iris)
  r <- fFtest(iris$Sepal.Length, gv(iris, -1L))
  rlm <- summary(lm(Sepal.Length ~., iris))
  expect_equal(unattrib(r)[1:4], unattrib(c(rlm$r.squared, rlm$fstatistic[c(2:3, 1L)])))
  # Same with weights:
  set.seed(101)
  w <- na_insert(w)
  r <- fFtest(iris$Sepal.Length, gv(iris, -1L), w = w)
  rlm <- summary(lm(Sepal.Length ~., weights = w, iris))
  expect_equal(unattrib(r)[1:4], unattrib(c(rlm$r.squared, rlm$fstatistic[c(2:3, 1L)])))
  rm(iris)

  if(NCRAN) {
  r <- fFtest(wlddev$PCGDP, qF(wlddev$year), wlddev[c("iso3c","LIFEEX")])
  # Same test done using lm:
  data <- na_omit(get_vars(wlddev, c("iso3c","year","PCGDP","LIFEEX")), na.attr = TRUE)
  full <- lm(PCGDP ~ LIFEEX + iso3c + qF(year), data)
  rest <- lm(PCGDP ~ LIFEEX + iso3c, data)
  ranv <- anova(rest, full)

  expect_equal(unattrib(r[1L, 1:4]), unlist(summary(full)[c("r.squared", "fstatistic")],
                                            use.names = FALSE)[c(1L, 3:4, 2L)])
  expect_equal(unattrib(r[2L, 1:4]), unlist(summary(rest)[c("r.squared", "fstatistic")],
                                            use.names = FALSE)[c(1L, 3:4, 2L)])
  expect_equal(rev(unattrib(r[1:2, 3L])), ranv$Res.Df)
  expect_equal(r[3L, 2L], na_rm(ranv$Df))
  expect_equal(r[3L, 4L], na_rm(ranv$F))
  expect_equal(r[3L, 5L], na_rm(ranv$`Pr(>F)`))

  # Same with weights:
  w <- abs(rnorm(fnrow(wlddev)))
  r <- fFtest(wlddev$PCGDP, qF(wlddev$year), wlddev[c("iso3c","LIFEEX")], w)
  full <- lm(PCGDP ~ LIFEEX + iso3c + qF(year), weights = w[-attr(data, "na.action")], data)
  rest <- lm(PCGDP ~ LIFEEX + iso3c, weights = w[-attr(data, "na.action")], data)
  ranv <- anova(rest, full)

  expect_equal(unattrib(r[1L, 1:4]), unlist(summary(full)[c("r.squared", "fstatistic")],
                                            use.names = FALSE)[c(1L, 3:4, 2L)])
  expect_equal(unattrib(r[2L, 1:4]), unlist(summary(rest)[c("r.squared", "fstatistic")],
                                            use.names = FALSE)[c(1L, 3:4, 2L)])
  expect_equal(rev(unattrib(r[1:2, 3L])), ranv$Res.Df)
  expect_equal(r[3L, 2L], na_rm(ranv$Df))
  expect_equal(r[3L, 4L], na_rm(ranv$F))
  expect_equal(r[3L, 5L], na_rm(ranv$`Pr(>F)`))

  }
})
