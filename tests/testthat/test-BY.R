context("BY")

# rm(list = ls())
set.seed(101)
x <- rnorm(100)
xNA <- x
xNA[sample.int(100,20)] <- NA
f <- as.factor(sort(sample.int(10, 100, TRUE)))
g <- GRP(mtcars, ~ cyl + vs + am)
f2 <- as.factor(sort(sample.int(6, 32, TRUE)))
m <- as.matrix(mtcars)
mNA <- na_insert(m)
mtcNA <- na_insert(mtcars)

na20 <- function(x) {
  x[is.na(x)] <- 0
  x
}

myscale <- function(x, na.rm = FALSE) (x - mean.default(x, na.rm = na.rm)) / sd(x, na.rm = na.rm)
mysumf <- function(x, na.rm = FALSE) c(N = sum(!is.na(x)), Mean = mean(x, na.rm = na.rm),
                                       SD = sd(x, na.rm = na.rm), Min = min(x, na.rm = na.rm),
                                       Max = max(x, na.rm = na.rm))

options(warn = -1)

test_that("BY.default works as intended", {

  # No missing values
  expect_equal(BY(x, f, sum), fsum(x, f))
  expect_equal(BY(x, f, sum, return = "list"), as.list(fsum(x, f)))
  expect_equal(BY(x, f, mean), fmean(x, f))
  expect_equal(BY(x, f, mean, return = "list"), as.list(fmean(x, f)))
  # BY(x, f, scale)
  expect_equal(BY(x, f, scale, use.g.names = FALSE),   fscale(x, f))
  expect_equal(BY(x, f, log, use.g.names = FALSE),   log(x))
  expect_equal(BY(x, f, quantile), unlist(lapply(split(x, f), quantile)))
  expect_equal(BY(x, f, quantile, expand.wide = TRUE),   t(sapply(split(x, f), quantile)))
  expect_equal(BY(x, f, quantile, return = "list"), lapply(split(x, f), quantile))
  expect_equal(BY(x, f, quantile, return = "list", expand.wide = TRUE), lapply(split(x, f), quantile)) # This should have no effect !!

  # Missing values removed
  expect_equal(BY(xNA, f, sum, na.rm = TRUE), na20(fsum(xNA, f)))
  expect_equal(BY(xNA, f, sum, return = "list", na.rm = TRUE), as.list(na20(fsum(xNA, f))))
  expect_equal(BY(xNA, f, mean, na.rm = TRUE), fmean(xNA, f))
  expect_equal(BY(xNA, f, mean, return = "list", na.rm = TRUE), as.list(fmean(xNA, f)))
  expect_equal(BY(xNA, f, scale, use.g.names = FALSE),   fscale(xNA, f))
  expect_equal(BY(xNA, f, quantile, na.rm = TRUE), unlist(lapply(split(xNA, f), quantile, na.rm = TRUE)))
  expect_equal(BY(xNA, f, quantile, expand.wide = TRUE, na.rm = TRUE),   t(sapply(split(xNA, f), quantile, na.rm = TRUE)))
  expect_equal(BY(xNA, f, quantile, return = "list", na.rm = TRUE), lapply(split(xNA, f), quantile, na.rm = TRUE))
  expect_equal(BY(xNA, f, quantile, return = "list", expand.wide = TRUE, na.rm = TRUE), lapply(split(xNA, f), quantile, na.rm = TRUE)) # This should have no effect !!

  # Missing values kept
  expect_equal(BY(xNA, f, sum), fsum(xNA, f, na.rm = FALSE))
  expect_equal(BY(xNA, f, sum, return = "list"), as.list(fsum(xNA, f, na.rm = FALSE)))
  expect_equal(BY(xNA, f, mean), fmean(xNA, f, na.rm = FALSE))
  expect_equal(BY(xNA, f, mean, return = "list"), as.list(fmean(xNA, f, na.rm = FALSE)))
  expect_equal(BY(xNA, f, myscale, use.g.names = FALSE),   fscale(xNA, f, na.rm = FALSE))
  expect_equal(BY(xNA, f, mysumf), unlist(lapply(split(xNA, f), mysumf)))
  expect_equal(BY(xNA, f, mysumf, expand.wide = TRUE),   t(sapply(split(xNA, f), mysumf)))
  expect_equal(BY(xNA, f, mysumf, return = "list"), lapply(split(xNA, f), mysumf))
  expect_equal(BY(xNA, f, mysumf, return = "list", expand.wide = TRUE), lapply(split(xNA, f), mysumf)) # This should have no effect !!

})

test_that("BY.matrix works as intended", {

  # No missing values
  expect_equal(BY(m, g, sum), fsum(m, g))
  expect_equal(BY(m, g, sum, return = "data.frame"), qDF(fsum(m, g)))
  expect_equal(BY(m, g, mean), fmean(m, g))
  expect_equal(BY(m, g, mean, return = "data.frame"), qDF(fmean(m, g)))
  # BY(m, g, scale)
  expect_equal(BY(m, f2, scale, use.g.names = FALSE),   fscale(m, f2))
  expect_equal(BY(m, f2, log, use.g.names = FALSE),   log(m))
  expect_equal(BY(m, f2, quantile), qM(lapply(mctl(m, names = TRUE), function(x) unlist(lapply(split(x, f2), quantile)))))
  expect_equal(setDimnames(BY(m, f2, quantile, expand.wide = TRUE), NULL), setDimnames(do.call(cbind, lapply(mctl(m, names = TRUE), function(x) do.call(rbind, lapply(split(x, f2), quantile)))), NULL))
  expect_equal(BY(m, f2, quantile, return = "data.frame"), qDF(qM(lapply(mctl(m, names = TRUE), function(x) unlist(lapply(split(x, f2), quantile))))))
  expect_equal(unname(BY(m, f2, quantile, return = "data.frame", expand.wide = TRUE)), unname(qDF(do.call(cbind, lapply(mctl(m, names = TRUE), function(x) do.call(rbind, lapply(split(x, f2), quantile)))))))

  # Missing values removed
  expect_equal(BY(mNA, g, sum, na.rm = TRUE), na20(fsum(mNA, g)))
  expect_equal(BY(mNA, g, sum, return = "data.frame", na.rm = TRUE), qDF(na20(fsum(mNA, g))))
  expect_equal(BY(mNA, g, mean, na.rm = TRUE), fmean(mNA, g))
  expect_equal(BY(mNA, g, mean, return = "data.frame", na.rm = TRUE), qDF(fmean(mNA, g)))
  expect_equal(BY(mNA, f2, scale, use.g.names = FALSE),   fscale(mNA, f2))
  expect_equal(BY(mNA, f2, log, use.g.names = FALSE),   log(mNA))
  expect_equal(BY(mNA, f2, quantile, na.rm = TRUE), qM(lapply(mctl(mNA, names = TRUE), function(x) unlist(lapply(split(x, f2), quantile, na.rm = TRUE)))))
  expect_equal(setDimnames(BY(mNA, f2, quantile, expand.wide = TRUE, na.rm = TRUE), NULL), setDimnames(do.call(cbind, lapply(mctl(mNA, names = TRUE), function(x) do.call(rbind, lapply(split(x, f2), quantile, na.rm = TRUE)))), NULL))
  expect_equal(BY(mNA, f2, quantile, return = "data.frame", na.rm = TRUE), qDF(qM(lapply(mctl(mNA, names = TRUE), function(x) unlist(lapply(split(x, f2), quantile, na.rm = TRUE))))))
  expect_equal(unname(BY(mNA, f2, quantile, return = "data.frame", expand.wide = TRUE, na.rm = TRUE)), unname(qDF(do.call(cbind, lapply(mctl(mNA, names = TRUE), function(x) do.call(rbind, lapply(split(x, f2), quantile, na.rm = TRUE)))))))

  # Missing values kept
  expect_equal(BY(mNA, g, sum), fsum(mNA, g, na.rm = FALSE))
  expect_equal(BY(mNA, g, sum, return = "data.frame"), qDF(fsum(mNA, g, na.rm = FALSE)))
  expect_equal(BY(mNA, g, mean), fmean(mNA, g, na.rm = FALSE))
  expect_equal(BY(mNA, g, mean, return = "data.frame"), qDF(fmean(mNA, g, na.rm = FALSE)))
  expect_equal(BY(mNA, f2, mysumf), qM(lapply(mctl(mNA, names = TRUE), function(x) unlist(lapply(split(x, f2), mysumf)))))
  expect_equal(setDimnames(BY(mNA, f2, mysumf, expand.wide = TRUE), NULL), setDimnames(do.call(cbind, lapply(mctl(mNA, names = TRUE), function(x) do.call(rbind, lapply(split(x, f2), mysumf)))), NULL))
  expect_equal(BY(mNA, f2, mysumf, return = "data.frame"), qDF(qM(lapply(mctl(mNA, names = TRUE), function(x) unlist(lapply(split(x, f2), mysumf))))))
  expect_equal(unname(BY(mNA, f2, mysumf, return = "data.frame", expand.wide = TRUE)), unname(qDF(do.call(cbind, lapply(mctl(mNA, names = TRUE), function(x) do.call(rbind, lapply(split(x, f2), mysumf)))))))

})

test_that("BY.data.frame works as intended", {

  # No missing values
  expect_equal(BY(mtcars, g, sum), fsum(mtcars, g))
  expect_equal(BY(mtcars, g, sum, return = "matrix"), qM(fsum(mtcars, g)))
  expect_equal(BY(mtcars, g, mean), fmean(mtcars, g))
  expect_equal(BY(mtcars, g, mean, return = "matrix"), qM(fmean(mtcars, g)))
  # BY(mtcars, g, scale)
  expect_equal(BY(mtcars, f2, scale, use.g.names = FALSE),   fscale(mtcars, f2))
  expect_equal(BY(mtcars, f2, log, use.g.names = FALSE),   log(mtcars))
  expect_equal(BY(mtcars, f2, quantile), qDF(qM(lapply(mtcars, function(x) unlist(lapply(split(x, f2), quantile))))))
  expect_equal(unname(BY(mtcars, f2, quantile, expand.wide = TRUE)), unname(qDF(do.call(cbind, lapply(mtcars, function(x) do.call(rbind, lapply(split(x, f2), quantile)))))))
  expect_equal(BY(mtcars, f2, quantile, return = "matrix"), qM(lapply(mtcars, function(x) unlist(lapply(split(x, f2), quantile)))))
  expect_equal(setDimnames(BY(mtcars, f2, quantile, return = "matrix", expand.wide = TRUE), NULL), setDimnames(do.call(cbind, lapply(mtcars, function(x) do.call(rbind, lapply(split(x, f2), quantile)))), NULL))

  # Missing values removed
  expect_equal(BY(mtcNA, g, sum, na.rm = TRUE), na20(fsum(mtcNA, g)))
  expect_equal(BY(mtcNA, g, sum, return = "matrix", na.rm = TRUE), na20(qM(fsum(mtcNA, g))))
  expect_equal(BY(mtcNA, g, mean, na.rm = TRUE), fmean(mtcNA, g))
  expect_equal(BY(mtcNA, g, mean, return = "matrix", na.rm = TRUE), qM(fmean(mtcNA, g)))
  expect_equal(BY(mtcNA, f2, scale, use.g.names = FALSE),   fscale(mtcNA, f2))
  expect_equal(BY(mtcNA, f2, log, use.g.names = FALSE),   log(mtcNA))
  expect_equal(BY(mtcNA, f2, quantile, na.rm = TRUE), qDF(qM(lapply(mtcNA, function(x) unlist(lapply(split(x, f2), quantile, na.rm = TRUE))))))
  expect_equal(unname(BY(mtcNA, f2, quantile, expand.wide = TRUE, na.rm = TRUE)), unname(qDF(do.call(cbind, lapply(mtcNA, function(x) do.call(rbind, lapply(split(x, f2), quantile, na.rm = TRUE)))))))
  expect_equal(BY(mtcNA, f2, quantile, return = "matrix", na.rm = TRUE), qM(lapply(mtcNA, function(x) unlist(lapply(split(x, f2), quantile, na.rm = TRUE)))))
  expect_equal(setDimnames(BY(mtcNA, f2, quantile, return = "matrix", expand.wide = TRUE, na.rm = TRUE), NULL), setDimnames(do.call(cbind, lapply(mtcNA, function(x) do.call(rbind, lapply(split(x, f2), quantile, na.rm = TRUE)))), NULL))

  # Missing values kept
  expect_equal(BY(mtcNA, g, sum), fsum(mtcNA, g, na.rm = FALSE))
  expect_equal(BY(mtcNA, g, sum, return = "matrix"), qM(fsum(mtcNA, g, na.rm = FALSE)))
  expect_equal(BY(mtcNA, g, mean), fmean(mtcNA, g, na.rm = FALSE))
  expect_equal(BY(mtcNA, g, mean, return = "matrix"), qM(fmean(mtcNA, g, na.rm = FALSE)))
  expect_equal(BY(mtcNA, f2, mysumf), qDF(qM(lapply(mtcNA, function(x) unlist(lapply(split(x, f2), mysumf))))))
  expect_equal(unname(BY(mtcNA, f2, mysumf, expand.wide = TRUE)), unname(qDF(do.call(cbind, lapply(mtcNA, function(x) do.call(rbind, lapply(split(x, f2), mysumf)))))))
  expect_equal(BY(mtcNA, f2, mysumf, return = "matrix"), qM(lapply(mtcNA, function(x) unlist(lapply(split(x, f2), mysumf)))))
  expect_equal(setDimnames(BY(mtcNA, f2, mysumf, return = "matrix", expand.wide = TRUE), NULL), setDimnames(do.call(cbind, lapply(mtcNA, function(x) do.call(rbind, lapply(split(x, f2), mysumf)))), NULL))

})

test_that("Output type is as expected", {

  expect_true(is.atomic(BY(x, f, sum)))
  expect_true(is.atomic(BY(xNA, f, sum, na.rm = TRUE)))
  expect_true(is.matrix(BY(mtcars, g, sum, return = "matrix")))
  expect_true(is.data.frame(BY(m, g, sum, return = "data.frame")))
  # BY(mtcars, g, quantile, expand.wide = TRUE, return = "list")
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

  expect_error(BY(~bla, g, sum)) # Not supported type
  expect_warning(BY(1, g, sum)) # This only gives a wrning in split.default: g is too long
  expect_warning(BY(x, g, sum)) # This only gives a wrning in split.default: g is too short
  expect_error(BY(letters, sample.int(5, length(letters), TRUE), sum)) # wrong type
  expect_error(BY(x, f, sum2)) # unknown object
  expect_error(BY(x, f, "sum2")) # unknown object
  expect_error(BY(x, f, log, bla = 1)) # unknown function argument
  expect_error(BY(x, f, sum, return = "bla")) # unknown return option
  expect_error(BY(m, g, sum2)) # unknown object
  expect_error(BY(m, g, "sum2")) # unknown object
  expect_error(BY(m, g, log, bla = 1)) # unknown function argument
  expect_error(BY(m, g, sum, return = "bla")) # unknown return option
  expect_error(BY(mtcars, g, sum2)) # unknown object
  expect_error(BY(mtcars, g, "sum2")) # unknown object
  expect_error(BY(mtcars, g, log, bla = 1)) # unknown function argument
  expect_error(BY(mtcars, g, sum, return = "bla")) # unknown return option
  expect_error(BY(mtcars, ~g, sum)) # Not supported type
  expect_error(BY(m, ~g, sum)) # Not supported type
  expect_error(BY(x, ~g, sum)) # Not supported type

})

options(warn = 1)
