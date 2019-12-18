context("BY")

x <- rnorm(100)
xNA <- x
xNA[sample.int(100,20)] <- NA
f <- as.factor(sample.int(10, 100, TRUE))
g <- GRP(mtcars, ~ cyl + vs + am)
m <- as.matrix(mtcars)

options(warn = -1)

test_that("All common uses of BY can be performed, as per examples", {

  BY(x, f, sum)
  BY(x, f, sum, return = "list")
  BY(x, f, scale)
  BY(x, f, scale, use.g.names = FALSE)
  BY(x, f, quantile)
  BY(x, f, quantile, expand.wide = TRUE)
  BY(x, f, quantile, expand.wide = TRUE, return = "list")

  BY(xNA, f, sum, na.rm = TRUE)
  BY(xNA, f, quantile, na.rm = TRUE)
  BY(xNA, f, quantile, na.rm = TRUE, expand.wide = TRUE)

  ## matrix method
  BY(m, g, scale)
  BY(m, g, quantile)
  BY(m, g, quantile, expand.wide = TRUE)
  BY(m, g, quantile, expand.wide = TRUE,
     return = "list")

  BY(mtcars, g, scale)
  BY(mtcars, g, scale, use.g.names = FALSE)
  BY(mtcars, g, quantile)
  BY(mtcars, g, quantile, expand.wide = TRUE)
  BY(mtcars, g, quantile, expand.wide = TRUE, return = "list")

  expect_true(is.atomic(BY(x, f, sum)))
  expect_true(is.atomic(BY(xNA, f, sum, na.rm = TRUE)))
  expect_true(is.matrix(BY(mtcars, g, sum, return = "matrix")))
  expect_true(is.data.frame(BY(m, g, sum, return = "data.frame")))

  expect_equal(BY(mtcars, g, quantile, return = "list", expand.wide = TRUE), BY(m, g, quantile, return = "list", expand.wide = TRUE))

})

test_that("BY matrix <> data.frame conversions run seamlessly", {
  expect_equal(BY(mtcars, g, sum, return = "matrix"), BY(m, g, sum))
  expect_equal(BY(mtcars, g, sum, return = "matrix", use.g.names = FALSE), BY(m, g, sum, use.g.names = FALSE))
  expect_equal(BY(m, g, sum, return = "data.frame"), BY(mtcars, g, sum))
  expect_equal(BY(m, g, sum, return = "data.frame", use.g.names = FALSE), BY(mtcars, g, sum, use.g.names = FALSE))

  expect_equal(BY(mtcars, g, log, return = "matrix"), BY(m, g, log))
  expect_equal(BY(mtcars, g, log, return = "matrix", use.g.names = FALSE), BY(m, g, log, use.g.names = FALSE))
  expect_equal(BY(m, g, log, return = "data.frame"), BY(mtcars, g, log))
  expect_equal(BY(m, g, log, return = "data.frame", use.g.names = FALSE), BY(mtcars, g, log, use.g.names = FALSE))

  expect_equal(BY(mtcars, g, quantile, return = "matrix"), BY(m, g, quantile))
  expect_equal(BY(mtcars, g, quantile, return = "matrix", use.g.names = FALSE), BY(m, g, quantile, use.g.names = FALSE))
  expect_equal(BY(m, g, quantile, return = "data.frame"), BY(mtcars, g, quantile))
  expect_equal(BY(m, g, quantile, return = "data.frame", use.g.names = FALSE), BY(mtcars, g, quantile, use.g.names = FALSE))
})

test_that("BY produces errors for wrong input", {
  expect_error(BY(letters, g, sum))
  expect_warning(BY(1:31, g, sum))
  expect_error(BY(mtcars, g, sum2))
  expect_error(BY(mtcars, g, log, bla = 1))
  expect_error(BY(mtcars, g, sum, return = "bla"))
})

options(warn = 1)
