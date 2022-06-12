context("dapply")

# rm(list = ls())
if(!is.null(attributes(identical(FALSE, TRUE)))) stop("OECD label issue")

test_that("All common uses of dapply can be performed, as per examples", {
  # data.frame
  expect_equal(dapply(mtcars, force), mtcars)
  expect_equal(dapply(`attr<-`(mtcars, "bla", 1), force), `attr<-`(mtcars, "bla", 1))
  expect_equal(dapply(`attr<-`(mtcars, "bla", 1), force, MARGIN = 1), `attr<-`(mtcars, "bla", 1))
  expect_visible(dapply(mtcars, log))
  expect_true(is.matrix(dapply(mtcars, log, return = "matrix")))

  # matrix
  m <- as.matrix(mtcars)
  expect_equal(dapply(m, force), m)
  expect_equal(dapply(EuStockMarkets, force), EuStockMarkets)
  expect_equal(dapply(EuStockMarkets, force, MARGIN = 1), EuStockMarkets)
  expect_visible(dapply(m, log))
  expect_true(is.data.frame(dapply(m, log, return = "data.frame")))

  # matrix <> data.frame conversions
  expect_equal(dapply(mtcars, log, return = "matrix"), dapply(m, log))
  expect_equal(dapply(mtcars, log, return = "matrix", MARGIN = 1), dapply(m, log, MARGIN = 1))
  expect_equal(dapply(m, log, return = "data.frame"), dapply(mtcars, log))
  expect_equal(dapply(m, log, return = "data.frame", MARGIN = 1), dapply(mtcars, log, MARGIN = 1))
  expect_equal(dapply(mtcars, quantile, return = "matrix"), dapply(m, quantile))
  expect_equal(dapply(mtcars, quantile, return = "matrix", MARGIN = 1), dapply(m, quantile, MARGIN = 1))
  expect_equal(dapply(m, quantile, return = "data.frame"), dapply(mtcars, quantile))
  expect_equal(dapply(m, quantile, return = "data.frame", MARGIN = 1), dapply(mtcars, quantile, MARGIN = 1))

  # scalar function gives atomic vector
  expect_true(is.atomic(dapply(mtcars, sum)))
  expect_equal(dapply(m, sum), dapply(mtcars, sum))
  expect_true(is.atomic(dapply(mtcars, sum, MARGIN = 1)))
  expect_equal(dapply(m, sum, MARGIN = 1), dapply(mtcars, sum, MARGIN = 1))

  # drop = FALSE retains object structure
  expect_true(is.data.frame(dapply(mtcars, sum, drop = FALSE)))
  expect_true(is.data.frame(dapply(mtcars, sum, MARGIN = 1, drop = FALSE)))
  expect_true(is.matrix(dapply(m, sum, drop = FALSE)))
  expect_true(is.matrix(dapply(m, sum, MARGIN = 1, drop = FALSE)))

  # matrix <> data.frame conversions without drop dimensions
  expect_equal(dapply(m, sum, drop = FALSE), dapply(mtcars, sum, return = "matrix", drop = FALSE))
  expect_equal(dapply(mtcars, sum, drop = FALSE), dapply(m, sum, return = "data.frame", drop = FALSE))

  # ... but if function is vector value, drop = FALSE does nothing
  expect_true(is.data.frame(dapply(mtcars, log, drop = FALSE)))
  expect_true(is.data.frame(dapply(mtcars, log, MARGIN = 1, drop = FALSE)))
  expect_true(is.data.frame(dapply(mtcars, quantile, drop = FALSE)))
  expect_true(is.data.frame(dapply(mtcars, quantile, MARGIN = 1, drop = FALSE)))
  expect_true(is.matrix(dapply(m, log, drop = FALSE)))
  expect_true(is.matrix(dapply(m, log, MARGIN = 1, drop = FALSE)))
  expect_true(is.matrix(dapply(m, quantile, drop = FALSE)))
  expect_true(is.matrix(dapply(m, quantile, MARGIN = 1, drop = FALSE)))

  # passing additional arguments works:
  dapply(mtcars, weighted.mean, mtcars$hp, na.rm = TRUE)
  dapply(m, weighted.mean, mtcars$hp, na.rm = TRUE)
})


test_that("dapply produces errors for wrong input", {
  expect_error(dapply("a", sum))
  expect_error(dapply(~ y, sum))
  expect_error(dapply(iris3, sum))
  expect_error(dapply(mtcars, sum2))
  expect_error(dapply(mtcars, sum, MARGIN = 3))
  expect_error(dapply(mtcars, sum, MARGIN = 1:2))
  expect_error(dapply(mtcars, sum, MARGIN = "a"))
  expect_error(dapply(mtcars, sum, return = "bla", drop = FALSE))
})
