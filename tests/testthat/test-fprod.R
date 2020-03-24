context("fprod")

x <- rnorm(100)
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

na21 <- function(x) {
  x[is.na(x)] <- 1
  x
}

test_that("fprod performs like base::prod", {
  expect_equal(fprod(NA), as.double(prod(NA)))
  expect_equal(fprod(NA, na.rm = FALSE), as.double(prod(NA)))
  expect_equal(fprod(1), prod(1, na.rm = TRUE))
  expect_equal(fprod(1:3), prod(1:3, na.rm = TRUE))
  expect_equal(fprod(-1:1), prod(-1:1, na.rm = TRUE))
  expect_equal(fprod(1, na.rm = FALSE), prod(1))
  expect_equal(fprod(1:3, na.rm = FALSE), prod(1:3))
  expect_equal(fprod(-1:1, na.rm = FALSE), prod(-1:1))
  expect_equal(fprod(x), prod(x, na.rm = TRUE))
  expect_equal(fprod(x, na.rm = FALSE), prod(x))
  expect_equal(fprod(xNA, na.rm = FALSE), prod(xNA))
  expect_equal(fprod(xNA), prod(xNA, na.rm = TRUE))
  expect_equal(fprod(mtcars), fprod(m))
  expect_equal(fprod(m), dapply(m, prod, na.rm = TRUE))
  expect_equal(fprod(m, na.rm = FALSE), dapply(m, prod))
  expect_equal(fprod(mNA, na.rm = FALSE), dapply(mNA, prod))
  expect_equal(fprod(mNA), dapply(mNA, prod, na.rm = TRUE))
  expect_equal(fprod(x, f), BY(x, f, prod, na.rm = TRUE))
  expect_equal(fprod(x, f, na.rm = FALSE), BY(x, f, prod))
  expect_equal(fprod(xNA, f, na.rm = FALSE), BY(xNA, f, prod))
  expect_equal(na21(fprod(xNA, f)), BY(xNA, f, prod, na.rm = TRUE))
  expect_equal(fprod(m, g), BY(m, g, prod, na.rm = TRUE))
  expect_equal(fprod(m, g, na.rm = FALSE), BY(m, g, prod))
  expect_equal(fprod(mNA, g, na.rm = FALSE), BY(mNA, g, prod))
  expect_equal(na21(fprod(mNA, g)), BY(mNA, g, prod, na.rm = TRUE)) # prod(NA, na.rm = TRUE) gives 1
  expect_equal(fprod(mtcars, g), BY(mtcars, g, prod, na.rm = TRUE))
  expect_equal(fprod(mtcars, g, na.rm = FALSE), BY(mtcars, g, prod))
  expect_equal(fprod(mtcNA, g, na.rm = FALSE), BY(mtcNA, g, prod))
  expect_equal(na21(fprod(mtcNA, g)), BY(mtcNA, g, prod, na.rm = TRUE)) # prod(NA, na.rm = TRUE) gives 1
})

test_that("fprod performs numerically stable", {
  expect_true(all_obj_equal(replicate(50, fprod(1), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(NA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(NA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(x), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(x, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(xNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(xNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(m), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(m, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mtcars), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mtcars, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mtcNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mtcNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(x, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(x, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(xNA, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(xNA, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(m, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(m, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mNA, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mtcars, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mtcars, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mtcNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fprod(mtcNA, g), simplify = FALSE)))
})

test_that("fprod handles special values in the right way", {
  expect_equal(fprod(NA), NA_real_)
  expect_equal(fprod(NaN), NaN)
  expect_equal(fprod(Inf), Inf)
  expect_equal(fprod(-Inf), -Inf)
  expect_equal(fprod(TRUE), 1)
  expect_equal(fprod(FALSE), 0)
  expect_equal(fprod(NA, na.rm = FALSE), NA_real_)
  expect_equal(fprod(NaN, na.rm = FALSE), NaN)
  expect_equal(fprod(Inf, na.rm = FALSE), Inf)
  expect_equal(fprod(-Inf, na.rm = FALSE), -Inf)
  expect_equal(fprod(TRUE, na.rm = FALSE), 1)
  expect_equal(fprod(FALSE, na.rm = FALSE), 0)
})

test_that("fprod produces errors for wrong input", {
  expect_error(fprod("a"))
  expect_error(fprod(NA_character_))
  expect_error(fprod(mNAc))
  expect_error(fprod(mNAc, f))
  expect_error(fprod(1:2,1:3))
  expect_error(fprod(m,1:31))
  expect_error(fprod(mtcars,1:31))
  expect_error(fprod(wlddev))
  expect_error(fprod(wlddev, wlddev$iso3c))

})
