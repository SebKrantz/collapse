context("BY")



bmean <- base::mean
bsd <- stats::sd
bsum <- base::sum
bmin <- base::min
bmax <- base::max
bscale <- base::scale

# rm(list = ls())
set.seed(101)
x <- rnorm(100)
xNA <- x
xNA[sample.int(100,20)] <- NA
fuo <- sample.int(10, 100, TRUE)
fo <- as.factor(sort(fuo))
fuo <- as.factor(fuo)
g <- GRP(mtcars, ~ cyl + vs + am)
f2uo <- sample.int(6, 32, TRUE)
f2o <- as.factor(sort(f2uo))
f2uo <- as.factor(f2uo)
m <- as.matrix(mtcars)
mNA <- na_insert(m)
mtcNA <- na_insert(mtcars)

na20 <- function(x) {
  x[is.na(x)] <- 0
  x
}

myscale <- function(x, na.rm = FALSE) (x - mean.default(x, na.rm = na.rm)) / bsd(x, na.rm = na.rm)
mysumf <- function(x, na.rm = FALSE) c(N = bsum(!is.na(x)), Mean = bmean(x, na.rm = na.rm),
                                       SD = bsd(x, na.rm = na.rm), Min = bmin(x, na.rm = na.rm),
                                       Max = bmax(x, na.rm = na.rm))

options(warn = -1)

test_that("BY.default works as intended", {

  for (f in list(fuo, fo)) {

  # No missing values
  expect_equal(BY(x, f, bsum), fsum(x, f))
  expect_equal(BY(x, f, bsum, return = "list"), as.list(fsum(x, f)))
  expect_equal(BY(x, f, bmean), fmean(x, f))
  expect_equal(BY(x, f, bmean, return = "list"), as.list(fmean(x, f)))
  # BY(x, f, bscale)
  expect_equal(BY(x, f, bscale, use.g.names = FALSE), fscale(x, f))
  expect_equal(BY(x, f, log, use.g.names = FALSE), log(x))
  expect_equal(BY(x, f, quantile), unlist(lapply(split(x, f), quantile)))
  expect_equal(BY(x, f, quantile, expand.wide = TRUE),   t(sapply(split(x, f), quantile)))
  expect_equal(BY(x, f, quantile, return = "list"), lapply(split(x, f), quantile))
  expect_equal(BY(x, f, quantile, return = "list", expand.wide = TRUE), lapply(split(x, f), quantile)) # This should have no effect !!

  # Missing values removed
  expect_equal(BY(xNA, f, bsum, na.rm = TRUE), na20(fsum(xNA, f)))
  expect_equal(BY(xNA, f, bsum, return = "list", na.rm = TRUE), as.list(na20(fsum(xNA, f))))
  expect_equal(BY(xNA, f, bmean, na.rm = TRUE), fmean(xNA, f))
  expect_equal(BY(xNA, f, bmean, return = "list", na.rm = TRUE), as.list(fmean(xNA, f)))
  expect_equal(BY(xNA, f, bscale, use.g.names = FALSE),   fscale(xNA, f))
  expect_equal(BY(xNA, f, quantile, na.rm = TRUE), unlist(lapply(split(xNA, f), quantile, na.rm = TRUE)))
  expect_equal(BY(xNA, f, quantile, expand.wide = TRUE, na.rm = TRUE),   t(sapply(split(xNA, f), quantile, na.rm = TRUE)))
  expect_equal(BY(xNA, f, quantile, return = "list", na.rm = TRUE), lapply(split(xNA, f), quantile, na.rm = TRUE))
  expect_equal(BY(xNA, f, quantile, return = "list", expand.wide = TRUE, na.rm = TRUE), lapply(split(xNA, f), quantile, na.rm = TRUE)) # This should have no effect !!

  # Missing values kept
  expect_equal(BY(xNA, f, bsum), fsum(xNA, f, na.rm = FALSE))
  expect_equal(BY(xNA, f, bsum, return = "list"), as.list(fsum(xNA, f, na.rm = FALSE)))
  expect_equal(BY(xNA, f, bmean), fmean(xNA, f, na.rm = FALSE))
  expect_equal(BY(xNA, f, bmean, return = "list"), as.list(fmean(xNA, f, na.rm = FALSE)))
  expect_equal(BY(xNA, f, myscale, use.g.names = FALSE),   fscale(xNA, f, na.rm = FALSE))
  expect_equal(BY(xNA, f, mysumf), unlist(lapply(split(xNA, f), mysumf)))
  expect_equal(BY(xNA, f, mysumf, expand.wide = TRUE),   t(sapply(split(xNA, f), mysumf)))
  expect_equal(BY(xNA, f, mysumf, return = "list"), lapply(split(xNA, f), mysumf))
  expect_equal(BY(xNA, f, mysumf, return = "list", expand.wide = TRUE), lapply(split(xNA, f), mysumf)) # This should have no effect !!

  }

})

test_that("BY.matrix works as intended", {

  for (f in list(g, f2uo, f2o)) {
  # No missing values
  expect_equal(BY(m, f, bsum), fsum(m, f))
  expect_equal(BY(m, f, bsum, return = "data.frame"), qDF(fsum(m, f)))
  expect_equal(BY(m, f, bmean), fmean(m, f))
  expect_equal(BY(m, f, bmean, return = "data.frame"), qDF(fmean(m, f)))

  expect_true(all_obj_equal(BY(m, f, bscale), BY(m, f, bscale, use.g.names = FALSE), fscale(m, f)))
  expect_true(all_obj_equal(BY(m, f, log), BY(m, f, log, use.g.names = FALSE), log(m)))

  # Missing values kept
  expect_equal(BY(mNA, f, bsum), fsum(mNA, f, na.rm = FALSE))
  expect_equal(BY(mNA, f, bsum, return = "data.frame"), qDF(fsum(mNA, f, na.rm = FALSE)))
  expect_equal(BY(mNA, f, bmean), fmean(mNA, f, na.rm = FALSE))
  expect_equal(BY(mNA, f, bmean, return = "data.frame"), qDF(fmean(mNA, f, na.rm = FALSE)))

  }

  for (f in list(f2uo, f2o)) {

  expect_equal(BY(m, f, quantile), qM(lapply(mctl(m, names = TRUE), function(x) unlist(lapply(split(x, f), quantile)))))
  expect_equal(setDimnames(BY(m, f, quantile, expand.wide = TRUE), NULL), setDimnames(do.call(cbind, lapply(mctl(m, names = TRUE), function(x) do.call(rbind, lapply(split(x, f), quantile)))), NULL))
  expect_equal(BY(m, f, quantile, return = "data.frame"), qDF(qM(lapply(mctl(m, names = TRUE), function(x) unlist(lapply(split(x, f), quantile))))))
  expect_equal(unname(BY(m, f, quantile, return = "data.frame", expand.wide = TRUE)), unname(qDF(do.call(cbind, lapply(mctl(m, names = TRUE), function(x) do.call(rbind, lapply(split(x, f), quantile)))))))

  # Missing values removed
  expect_equal(BY(mNA, f, bsum, na.rm = TRUE), na20(fsum(mNA, f)))
  expect_equal(BY(mNA, f, bsum, return = "data.frame", na.rm = TRUE), qDF(na20(fsum(mNA, f))))
  expect_equal(BY(mNA, f, bmean, na.rm = TRUE), fmean(mNA, f))
  expect_equal(BY(mNA, f, bmean, return = "data.frame", na.rm = TRUE), qDF(fmean(mNA, f)))
  expect_true(all_obj_equal(BY(mNA, f, bscale), BY(mNA, f, bscale, use.g.names = FALSE), fscale(mNA, f)))
  expect_true(all_obj_equal(BY(mNA, f, log), BY(mNA, f, log, use.g.names = FALSE), log(mNA)))

  expect_equal(BY(mNA, f, quantile, na.rm = TRUE), qM(lapply(mctl(mNA, names = TRUE), function(x) unlist(lapply(split(x, f), quantile, na.rm = TRUE)))))
  expect_equal(setDimnames(BY(mNA, f, quantile, expand.wide = TRUE, na.rm = TRUE), NULL), setDimnames(do.call(cbind, lapply(mctl(mNA, names = TRUE), function(x) do.call(rbind, lapply(split(x, f), quantile, na.rm = TRUE)))), NULL))
  expect_equal(BY(mNA, f, quantile, return = "data.frame", na.rm = TRUE), qDF(qM(lapply(mctl(mNA, names = TRUE), function(x) unlist(lapply(split(x, f), quantile, na.rm = TRUE))))))
  expect_equal(unname(BY(mNA, f, quantile, return = "data.frame", expand.wide = TRUE, na.rm = TRUE)), unname(qDF(do.call(cbind, lapply(mctl(mNA, names = TRUE), function(x) do.call(rbind, lapply(split(x, f), quantile, na.rm = TRUE)))))))

  # Missing values kept
  expect_equal(BY(mNA, f, mysumf), qM(lapply(mctl(mNA, names = TRUE), function(x) unlist(lapply(split(x, f), mysumf)))))
  expect_equal(setDimnames(BY(mNA, f, mysumf, expand.wide = TRUE), NULL), setDimnames(do.call(cbind, lapply(mctl(mNA, names = TRUE), function(x) do.call(rbind, lapply(split(x, f), mysumf)))), NULL))
  expect_equal(BY(mNA, f, mysumf, return = "data.frame"), qDF(qM(lapply(mctl(mNA, names = TRUE), function(x) unlist(lapply(split(x, f), mysumf))))))
  expect_equal(unname(BY(mNA, f, mysumf, return = "data.frame", expand.wide = TRUE)), unname(qDF(do.call(cbind, lapply(mctl(mNA, names = TRUE), function(x) do.call(rbind, lapply(split(x, f), mysumf)))))))

  }

})

test_that("BY.data.frame works as intended", {

  condsetrn <- function(x) if(is.null(rownames(x))) setRownames(x) else x

  for (f in list(g, f2uo, f2o)) {

  # No missing values
  expect_equal(BY(mtcars, f, bsum), fsum(mtcars, f))
  expect_equal(BY(mtcars, f, bsum, return = "matrix"), condsetrn(qM(fsum(mtcars, f))))
  expect_equal(BY(mtcars, f, bmean), fmean(mtcars, f))
  expect_equal(BY(mtcars, f, bmean, return = "matrix"), condsetrn(qM(fmean(mtcars, f))))
  # BY(mtcars, f, bscale)
  expect_equal(BY(mtcars, f, bscale, use.g.names = FALSE),  fscale(mtcars, f))
  expect_equal(BY(mtcars, f, log, use.g.names = FALSE),  log(mtcars))

  # Missing values removed
  expect_equal(BY(mtcNA, f, bsum, na.rm = TRUE), na20(fsum(mtcNA, f)))
  expect_equal(BY(mtcNA, f, bsum, return = "matrix", na.rm = TRUE), condsetrn(na20(qM(fsum(mtcNA, f)))))
  expect_equal(BY(mtcNA, f, bmean, na.rm = TRUE), fmean(mtcNA, f))
  expect_equal(BY(mtcNA, f, bmean, return = "matrix", na.rm = TRUE), condsetrn(qM(fmean(mtcNA, f))))
  expect_equal(BY(mtcNA, f, bscale, use.g.names = FALSE), fscale(mtcNA, f))
  expect_equal(BY(mtcNA, f, log, use.g.names = FALSE), log(mtcNA))

  # Missing values kept
  expect_equal(BY(mtcNA, f, bsum), fsum(mtcNA, f, na.rm = FALSE))
  expect_equal(BY(mtcNA, f, bsum, return = "matrix"), condsetrn(qM(fsum(mtcNA, f, na.rm = FALSE))))
  expect_equal(BY(mtcNA, f, bmean), fmean(mtcNA, f, na.rm = FALSE))
  expect_equal(BY(mtcNA, f, bmean, return = "matrix"), condsetrn(qM(fmean(mtcNA, f, na.rm = FALSE))))

  }

  for (f in list(f2uo, f2o)) {

  # No missing values
  expect_equal(BY(mtcars, f, quantile), qDF(qM(lapply(mtcars, function(x) unlist(lapply(split(x, f), quantile))))))
  expect_equal(unname(BY(mtcars, f, quantile, expand.wide = TRUE)), unname(qDF(do.call(cbind, lapply(mtcars, function(x) do.call(rbind, lapply(split(x, f), quantile)))))))
  expect_equal(BY(mtcars, f, quantile, return = "matrix"), qM(lapply(mtcars, function(x) unlist(lapply(split(x, f), quantile)))))
  expect_equal(setDimnames(BY(mtcars, f, quantile, return = "matrix", expand.wide = TRUE), NULL), setDimnames(do.call(cbind, lapply(mtcars, function(x) do.call(rbind, lapply(split(x, f), quantile)))), NULL))

  # Missing values removed
  expect_equal(BY(mtcNA, f, quantile, na.rm = TRUE), qDF(qM(lapply(mtcNA, function(x) unlist(lapply(split(x, f), quantile, na.rm = TRUE))))))
  expect_equal(unname(BY(mtcNA, f, quantile, expand.wide = TRUE, na.rm = TRUE)), unname(qDF(do.call(cbind, lapply(mtcNA, function(x) do.call(rbind, lapply(split(x, f), quantile, na.rm = TRUE)))))))
  expect_equal(BY(mtcNA, f, quantile, return = "matrix", na.rm = TRUE), qM(lapply(mtcNA, function(x) unlist(lapply(split(x, f), quantile, na.rm = TRUE)))))
  expect_equal(setDimnames(BY(mtcNA, f, quantile, return = "matrix", expand.wide = TRUE, na.rm = TRUE), NULL), setDimnames(do.call(cbind, lapply(mtcNA, function(x) do.call(rbind, lapply(split(x, f), quantile, na.rm = TRUE)))), NULL))

  # Missing values kept
  expect_equal(BY(mtcNA, f, mysumf), qDF(qM(lapply(mtcNA, function(x) unlist(lapply(split(x, f), mysumf))))))
  expect_equal(unname(BY(mtcNA, f, mysumf, expand.wide = TRUE)), unname(qDF(do.call(cbind, lapply(mtcNA, function(x) do.call(rbind, lapply(split(x, f), mysumf)))))))
  expect_equal(BY(mtcNA, f, mysumf, return = "matrix"), qM(lapply(mtcNA, function(x) unlist(lapply(split(x, f), mysumf)))))
  expect_equal(setDimnames(BY(mtcNA, f, mysumf, return = "matrix", expand.wide = TRUE), NULL), setDimnames(do.call(cbind, lapply(mtcNA, function(x) do.call(rbind, lapply(split(x, f), mysumf)))), NULL))

  }

})

test_that("Output type is as expected", {

  expect_true(is.atomic(BY(x, fuo, bsum)))
  expect_true(is.atomic(BY(xNA, fuo, bsum, na.rm = TRUE)))
  expect_true(is.matrix(BY(mtcars, g, bsum, return = "matrix")))
  expect_true(is.data.frame(BY(m, g, bsum, return = "data.frame")))
  # BY(mtcars, g, quantile, expand.wide = TRUE, return = "list")
  expect_equal(BY(mtcars, g, quantile, return = "list", expand.wide = TRUE), BY(m, g, quantile, return = "list", expand.wide = TRUE))

})

test_that("BY matrix <> data.frame conversions run seamlessly", {
  expect_equal(BY(mtcars, g, bsum, return = "matrix"), BY(m, g, bsum))
  expect_equal(BY(mtcars, g, bsum, return = "matrix", use.g.names = FALSE), BY(m, g, bsum, use.g.names = FALSE))
  expect_equal(BY(m, g, bsum, return = "data.frame"), BY(mtcars, g, bsum))
  expect_equal(BY(m, g, bsum, return = "data.frame", use.g.names = FALSE), BY(mtcars, g, bsum, use.g.names = FALSE))

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

  expect_error(BY(~bla, g, bsum)) # Not supported type
  expect_error(BY(1, g, bsum)) # This only gives a warning in gsplit: g is too long
  expect_error(BY(x, g, bsum)) # This only gives a warning in gsplit: g is too short
  expect_error(BY(letters, sample.int(5, length(letters), TRUE), bsum)) # wrong type
  expect_error(BY(x, f, sum2)) # unknown object
  expect_error(BY(x, f, "sum2")) # unknown object
  expect_error(BY(x, f, log, bla = 1)) # unknown function argument
  expect_error(BY(x, f, bsum, return = "bla")) # unknown return option
  expect_error(BY(m, g, sum2)) # unknown object
  expect_error(BY(m, g, "sum2")) # unknown object
  expect_error(BY(m, g, log, bla = 1)) # unknown function argument
  expect_error(BY(m, g, bsum, return = "bla")) # unknown return option
  expect_error(BY(mtcars, g, sum2)) # unknown object
  expect_error(BY(mtcars, g, "sum2")) # unknown object
  expect_error(BY(mtcars, g, log, bla = 1)) # unknown function argument
  expect_error(BY(mtcars, g, bsum, return = "bla")) # unknown return option
  expect_error(BY(mtcars, ~g, bsum)) # Not supported type
  expect_error(BY(m, ~g, bsum)) # Not supported type
  expect_error(BY(x, ~g, bsum)) # Not supported type

})

test_that("no row-names are generated for data.table's (only)", {

  mtcDT <- qDT(mtcars)

  for(FUN in list(bsum, quantile, identity)) {

    expect_false(is.character(attr(BY(mtcDT, g, FUN), "row.names")))
    if(!identical(FUN, identity)) {
      expect_true(is.character(attr(BY(mtcDT, g, FUN, return = "data.frame"), "row.names")))
      expect_true(is.character(dimnames(BY(mtcDT, g, FUN, return = "matrix"))[[1L]]))
    }
    expect_false(is.character(attr(BY(mtcDT, g, FUN, use.g.names = FALSE), "row.names")))
    expect_false(is.character(attr(BY(mtcDT, g, FUN, use.g.names = FALSE, return = "data.frame"), "row.names")))
    expect_false(is.character(dimnames(BY(mtcDT, g, FUN, use.g.names = FALSE, return = "matrix"))[[1L]]))

  }

})

options(warn = 1)
