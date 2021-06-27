context("fmin and fmax")

# rm(list = ls())
set.seed(101)
x <- rnorm(100) * 10000
xNA <- x
xNA[sample.int(100,20)] <- NA
f <- as.factor(sample.int(10, 100, TRUE))
g <- GRP(mtcars, ~ cyl + vs + am)
mtcNA <- na_insert(mtcars)
mtcNA[27,1] <- NA # single group NA !!
m <- as.matrix(mtcars)
mNA <- as.matrix(mtcNA)
mNAc <- mNA
storage.mode(mNAc) <- "character"

inf2NA <- function(x) {
  if(is.atomic(x)) {
   x[is.infinite(x)] <- NA
  } else {
   x[do.call(cbind, lapply(x, is.infinite))] <- NA
  }
  x
}

options(warn = -1)

# fmin double

test_that("fmin performs like base::min", {
  expect_equal(fmin(NA), min(NA))
  expect_equal(fmin(NA, na.rm = FALSE), min(NA))
  expect_equal(fmin(1), min(1, na.rm = TRUE))
  expect_equal(fmin(1:3), min(1:3, na.rm = TRUE))
  expect_equal(fmin(-1:1), min(-1:1, na.rm = TRUE))
  expect_equal(fmin(1, na.rm = FALSE), min(1))
  expect_equal(fmin(1:3, na.rm = FALSE), min(1:3))
  expect_equal(fmin(-1:1, na.rm = FALSE), min(-1:1))
  expect_equal(fmin(x), min(x, na.rm = TRUE))
  expect_equal(fmin(x, na.rm = FALSE), min(x))
  expect_equal(fmin(xNA, na.rm = FALSE), min(xNA))
  expect_equal(fmin(xNA), min(xNA, na.rm = TRUE))
  expect_equal(fmin(mtcars), fmin(m))
  expect_equal(fmin(m), dapply(m, min, na.rm = TRUE))
  expect_equal(fmin(m, na.rm = FALSE), dapply(m, min))
  expect_equal(fmin(mNA, na.rm = FALSE), dapply(mNA, min))
  expect_equal(fmin(mNA), dapply(mNA, min, na.rm = TRUE))
  expect_equal(fmin(mtcars), dapply(mtcars, min, na.rm = TRUE))
  expect_equal(fmin(mtcars, na.rm = FALSE), dapply(mtcars, min))
  expect_equal(fmin(mtcNA, na.rm = FALSE), dapply(mtcNA, min))
  expect_equal(fmin(mtcNA), dapply(mtcNA, min, na.rm = TRUE))
  expect_equal(fmin(x, f), BY(x, f, min, na.rm = TRUE))
  expect_equal(fmin(x, f, na.rm = FALSE), BY(x, f, min))
  expect_equal(fmin(xNA, f, na.rm = FALSE), BY(xNA, f, min))
  expect_equal(fmin(xNA, f), inf2NA(BY(xNA, f, min, na.rm = TRUE)))
  expect_equal(fmin(m, g), BY(m, g, min, na.rm = TRUE))
  expect_equal(fmin(m, g, na.rm = FALSE), BY(m, g, min))
  expect_equal(fmin(mNA, g, na.rm = FALSE), BY(mNA, g, min))
  expect_equal(fmin(mNA, g), inf2NA(BY(mNA, g, min, na.rm = TRUE))) # min(NA, na.rm = TRUE) gives Inf
  expect_equal(fmin(mtcars, g), BY(mtcars, g, min, na.rm = TRUE))
  expect_equal(fmin(mtcars, g, na.rm = FALSE), BY(mtcars, g, min))
  expect_equal(fmin(mtcNA, g, na.rm = FALSE), BY(mtcNA, g, min))
  expect_equal(fmin(mtcNA, g), inf2NA(BY(mtcNA, g, min, na.rm = TRUE))) # min(NA, na.rm = TRUE) gives Inf
})

test_that("fmin performs numerically stable", {
  expect_true(all_obj_equal(replicate(50, fmin(1), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(NA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(NA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(x), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(x, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(xNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(xNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(m), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(m, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mtcars), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mtcars, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mtcNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mtcNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(x, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(x, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(xNA, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(xNA, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(m, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(m, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mNA, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mtcars, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mtcars, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mtcNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmin(mtcNA, g), simplify = FALSE)))
})

test_that("fmin handles special values in the right way", {
  expect_equal(fmin(NA), NA_real_)
  expect_equal(fmin(NaN), NaN)
  expect_equal(fmin(Inf), Inf)
  expect_equal(fmin(-Inf), -Inf)
  expect_equal(fmin(TRUE), 1)
  expect_equal(fmin(FALSE), 0)
  expect_equal(fmin(NA, na.rm = FALSE), NA_real_)
  expect_equal(fmin(NaN, na.rm = FALSE), NaN)
  expect_equal(fmin(Inf, na.rm = FALSE), Inf)
  expect_equal(fmin(-Inf, na.rm = FALSE), -Inf)
  expect_equal(fmin(TRUE, na.rm = FALSE), 1)
  expect_equal(fmin(FALSE, na.rm = FALSE), 0)
})

test_that("fmin produces errors for wrong input", {
  expect_error(fmin("a"))
  expect_error(fmin(NA_character_))
  expect_error(fmin(mNAc))
  expect_error(fmin(mNAc, f))
  expect_error(fmin(1:2,1:3))
  expect_error(fmin(m,1:31))
  expect_error(fmin(mtcars,1:31))
  expect_error(fmin(wlddev))
  expect_error(fmin(wlddev, wlddev$iso3c))

})


# fmax double

test_that("fmax performs like base::max", {
  expect_equal(fmax(NA), max(NA))
  expect_equal(fmax(NA, na.rm = FALSE), max(NA))
  expect_equal(fmax(1), max(1, na.rm = TRUE))
  expect_equal(fmax(1:3), max(1:3, na.rm = TRUE))
  expect_equal(fmax(-1:1), max(-1:1, na.rm = TRUE))
  expect_equal(fmax(1, na.rm = FALSE), max(1))
  expect_equal(fmax(1:3, na.rm = FALSE), max(1:3))
  expect_equal(fmax(-1:1, na.rm = FALSE), max(-1:1))
  expect_equal(fmax(x), max(x, na.rm = TRUE))
  expect_equal(fmax(x, na.rm = FALSE), max(x))
  expect_equal(fmax(xNA, na.rm = FALSE), max(xNA))
  expect_equal(fmax(xNA), max(xNA, na.rm = TRUE))
  expect_equal(fmax(mtcars), fmax(m))
  expect_equal(fmax(m), dapply(m, max, na.rm = TRUE))
  expect_equal(fmax(m, na.rm = FALSE), dapply(m, max))
  expect_equal(fmax(mNA, na.rm = FALSE), dapply(mNA, max))
  expect_equal(fmax(mNA), dapply(mNA, max, na.rm = TRUE))
  expect_equal(fmax(mtcars), dapply(mtcars, max, na.rm = TRUE))
  expect_equal(fmax(mtcars, na.rm = FALSE), dapply(mtcars, max))
  expect_equal(fmax(mtcNA, na.rm = FALSE), dapply(mtcNA, max))
  expect_equal(fmax(mtcNA), dapply(mtcNA, max, na.rm = TRUE))
  expect_equal(fmax(x, f), BY(x, f, max, na.rm = TRUE))
  expect_equal(fmax(x, f, na.rm = FALSE), BY(x, f, max))
  expect_equal(fmax(xNA, f, na.rm = FALSE), BY(xNA, f, max))
  expect_equal(fmax(xNA, f), inf2NA(BY(xNA, f, max, na.rm = TRUE)))
  expect_equal(fmax(m, g), BY(m, g, max, na.rm = TRUE))
  expect_equal(fmax(m, g, na.rm = FALSE), BY(m, g, max))
  expect_equal(fmax(mNA, g, na.rm = FALSE), BY(mNA, g, max))
  expect_equal(fmax(mNA, g), inf2NA(BY(mNA, g, max, na.rm = TRUE))) # max(NA, na.rm = TRUE) gives -Inf
  expect_equal(fmax(mtcars, g), BY(mtcars, g, max, na.rm = TRUE))
  expect_equal(fmax(mtcars, g, na.rm = FALSE), BY(mtcars, g, max))
  expect_equal(fmax(mtcNA, g, na.rm = FALSE), BY(mtcNA, g, max))
  expect_equal(fmax(mtcNA, g), inf2NA(BY(mtcNA, g, max, na.rm = TRUE))) # max(NA, na.rm = TRUE) gives -Inf
})

test_that("fmax performs numerically stable", {
  expect_true(all_obj_equal(replicate(50, fmax(1), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(NA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(NA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(x), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(x, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(xNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(xNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(m), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(m, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mtcars), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mtcars, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mtcNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mtcNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(x, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(x, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(xNA, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(xNA, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(m, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(m, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mNA, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mtcars, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mtcars, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mtcNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fmax(mtcNA, g), simplify = FALSE)))
})

test_that("fmax handles special values in the right way", {
  expect_equal(fmax(NA), NA_real_)
  expect_equal(fmax(NaN), NaN)
  expect_equal(fmax(Inf), Inf)
  expect_equal(fmax(-Inf), -Inf)
  expect_equal(fmax(TRUE), 1)
  expect_equal(fmax(FALSE), 0)
  expect_equal(fmax(NA, na.rm = FALSE), NA_real_)
  expect_equal(fmax(NaN, na.rm = FALSE), NaN)
  expect_equal(fmax(Inf, na.rm = FALSE), Inf)
  expect_equal(fmax(-Inf, na.rm = FALSE), -Inf)
  expect_equal(fmax(TRUE, na.rm = FALSE), 1)
  expect_equal(fmax(FALSE, na.rm = FALSE), 0)
})

test_that("fmax produces errors for wrong input", {
  expect_error(fmax("a"))
  expect_error(fmax(NA_character_))
  expect_error(fmax(mNAc))
  expect_error(fmax(mNAc, f))
  expect_error(fmax(1:2,1:3))
  expect_error(fmax(m,1:31))
  expect_error(fmax(mtcars,1:31))
  expect_error(fmax(wlddev))
  expect_error(fmax(wlddev, wlddev$iso3c))

})



# fmin int

x <- as.integer(x)
xNA <- as.integer(xNA)
mtcNA <- dapply(mtcNA, as.integer)
mtcars <- dapply(mtcars, as.integer)
storage.mode(m) <- "integer"
storage.mode(mNA) <- "integer"

toint <- function(x) {
  storage.mode(x) <- "integer"
  x
}

test_that("fmin with integers performs like base::min", {
  expect_identical(fmin(x), min(x, na.rm = TRUE))
  expect_identical(fmin(x, na.rm = FALSE), min(x))
  expect_identical(fmin(xNA, na.rm = FALSE), min(xNA))
  expect_identical(fmin(xNA), min(xNA, na.rm = TRUE))
  expect_identical(toint(fmin(mtcars)), fmin(m))
  expect_identical(fmin(m), dapply(m, min, na.rm = TRUE))
  expect_identical(fmin(m, na.rm = FALSE), dapply(m, min))
  expect_identical(fmin(mNA, na.rm = FALSE), dapply(mNA, min))
  expect_identical(fmin(mNA), dapply(mNA, min, na.rm = TRUE))
  expect_identical(toint(fmin(mtcars)), dapply(mtcars, min, na.rm = TRUE))
  expect_identical(toint(fmin(mtcars, na.rm = FALSE)), dapply(mtcars, min))
  expect_identical(toint(fmin(mtcNA, na.rm = FALSE)), dapply(mtcNA, min))
  expect_identical(toint(fmin(mtcNA)), dapply(mtcNA, min, na.rm = TRUE))
  expect_identical(fmin(x, f), BY(x, f, min, na.rm = TRUE))
  expect_identical(fmin(x, f, na.rm = FALSE), BY(x, f, min))
  expect_identical(fmin(xNA, f, na.rm = FALSE), BY(xNA, f, min))
  expect_identical(fmin(xNA, f), inf2NA(BY(xNA, f, min, na.rm = TRUE)))
  expect_identical(fmin(m, g), BY(m, g, min, na.rm = TRUE))
  expect_identical(fmin(m, g, na.rm = FALSE), BY(m, g, min))
  expect_identical(fmin(mNA, g, na.rm = FALSE), BY(mNA, g, min))
  expect_identical(fmin(mNA, g), toint(inf2NA(BY(mNA, g, min, na.rm = TRUE)))) # min(NA, na.rm = TRUE) gives Inf
  expect_identical(fmin(mtcars, g), BY(mtcars, g, min, na.rm = TRUE))
  expect_identical(fmin(mtcars, g, na.rm = FALSE), BY(mtcars, g, min))
  expect_identical(fmin(mtcNA, g, na.rm = FALSE), BY(mtcNA, g, min))
  expect_identical(fmin(mtcNA, g), dapply(inf2NA(BY(mtcNA, g, min, na.rm = TRUE)), toint)) # min(NA, na.rm = TRUE) gives Inf
})

test_that("fmin with integers performs numerically stable", {
  expect_true(all_identical(replicate(50, fmin(x), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(x, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(xNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(xNA), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(m), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(m, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mNA), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mtcars), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mtcars, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mtcNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mtcNA), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(x, f), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(x, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(xNA, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(xNA, f), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(m, g), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(m, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mNA, g), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mtcars, g), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mtcars, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mtcNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmin(mtcNA, g), simplify = FALSE)))
})

test_that("fmin with integers produces errors for wrong input", {
  expect_error(fmin(m,1:31))
  expect_error(fmin(mtcars,1:31))
})


# fmax int

test_that("fmax with integers performs like base::max", {
  expect_identical(fmax(x), max(x, na.rm = TRUE))
  expect_identical(fmax(x, na.rm = FALSE), max(x))
  expect_identical(fmax(xNA, na.rm = FALSE), max(xNA))
  expect_identical(fmax(xNA), max(xNA, na.rm = TRUE))
  expect_identical(toint(fmax(mtcars)), fmax(m))
  expect_identical(fmax(m), dapply(m, max, na.rm = TRUE))
  expect_identical(fmax(m, na.rm = FALSE), dapply(m, max))
  expect_identical(fmax(mNA, na.rm = FALSE), dapply(mNA, max))
  expect_identical(fmax(mNA), dapply(mNA, max, na.rm = TRUE))
  expect_identical(toint(fmax(mtcars)), dapply(mtcars, max, na.rm = TRUE))
  expect_identical(toint(fmax(mtcars, na.rm = FALSE)), dapply(mtcars, max))
  expect_identical(toint(fmax(mtcNA, na.rm = FALSE)), dapply(mtcNA, max))
  expect_identical(toint(fmax(mtcNA)), dapply(mtcNA, max, na.rm = TRUE))
  expect_identical(fmax(x, f), BY(x, f, max, na.rm = TRUE))
  expect_identical(fmax(x, f, na.rm = FALSE), BY(x, f, max))
  expect_identical(fmax(xNA, f, na.rm = FALSE), BY(xNA, f, max))
  expect_identical(fmax(xNA, f), inf2NA(BY(xNA, f, max, na.rm = TRUE)))
  expect_identical(fmax(m, g), BY(m, g, max, na.rm = TRUE))
  expect_identical(fmax(m, g, na.rm = FALSE), BY(m, g, max))
  expect_identical(fmax(mNA, g, na.rm = FALSE), BY(mNA, g, max))
  expect_identical(fmax(mNA, g), toint(inf2NA(BY(mNA, g, max, na.rm = TRUE)))) # max(NA, na.rm = TRUE) gives -Inf
  expect_identical(fmax(mtcars, g), BY(mtcars, g, max, na.rm = TRUE))
  expect_identical(fmax(mtcars, g, na.rm = FALSE), BY(mtcars, g, max))
  expect_identical(fmax(mtcNA, g, na.rm = FALSE), BY(mtcNA, g, max))
  expect_identical(fmax(mtcNA, g), dapply(inf2NA(BY(mtcNA, g, max, na.rm = TRUE)), toint)) # max(NA, na.rm = TRUE) gives -Inf
})

test_that("fmax with integers performs numerically stable", {
  expect_true(all_identical(replicate(50, fmax(x), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(x, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(xNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(xNA), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(m), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(m, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mNA), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mtcars), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mtcars, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mtcNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mtcNA), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(x, f), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(x, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(xNA, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(xNA, f), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(m, g), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(m, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mNA, g), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mtcars, g), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mtcars, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mtcNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_identical(replicate(50, fmax(mtcNA, g), simplify = FALSE)))
})


test_that("fmax with integers produces errors for wrong input", {
  expect_error(fmax(m,1:31))
  expect_error(fmax(mtcars,1:31))

})

options(warn = 1)
