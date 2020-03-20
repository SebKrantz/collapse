context("fsum")

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

na20 <- function(x) {
  x[is.na(x)] <- 0
  x
}


test_that("fsum performs like base::sum and base::colSums", {
  expect_equal(fsum(NA), sum(NA))
  expect_equal(fsum(NA, na.rm = FALSE), sum(NA))
  expect_equal(fsum(1), sum(1, na.rm = TRUE))
  expect_equal(fsum(1:3), sum(1:3, na.rm = TRUE))
  expect_equal(fsum(-1:1), sum(-1:1, na.rm = TRUE))
  expect_equal(fsum(1, na.rm = FALSE), sum(1))
  expect_equal(fsum(1:3, na.rm = FALSE), sum(1:3))
  expect_equal(fsum(-1:1, na.rm = FALSE), sum(-1:1))
  expect_equal(fsum(x), sum(x, na.rm = TRUE))
  expect_equal(fsum(x, na.rm = FALSE), sum(x))
  expect_equal(fsum(xNA, na.rm = FALSE), sum(xNA))
  expect_equal(fsum(xNA), sum(xNA, na.rm = TRUE))
  expect_equal(fsum(mtcars), fsum(m))
  expect_equal(fsum(m), colSums(m, na.rm = TRUE))
  expect_equal(fsum(m, na.rm = FALSE), colSums(m))
  expect_equal(fsum(mNA, na.rm = FALSE), colSums(mNA))
  expect_equal(fsum(mNA), colSums(mNA, na.rm = TRUE))
  expect_equal(fsum(x, f), BY(x, f, sum, na.rm = TRUE))
  expect_equal(fsum(x, f, na.rm = FALSE), BY(x, f, sum))
  expect_equal(fsum(xNA, f, na.rm = FALSE), BY(xNA, f, sum))
  expect_equal(fsum(xNA, f), BY(xNA, f, sum, na.rm = TRUE))
  expect_equal(fsum(m, g), BY(m, g, sum, na.rm = TRUE))
  expect_equal(fsum(m, g, na.rm = FALSE), BY(m, g, sum))
  expect_equal(fsum(mNA, g, na.rm = FALSE), BY(mNA, g, sum))
  expect_equal(na20(fsum(mNA, g)), BY(mNA, g, sum, na.rm = TRUE)) # error, sum(NA) give 0
  expect_equal(fsum(mtcars, g), BY(mtcars, g, sum, na.rm = TRUE))
  expect_equal(fsum(mtcars, g, na.rm = FALSE), BY(mtcars, g, sum))
  expect_equal(fsum(mtcNA, g, na.rm = FALSE), BY(mtcNA, g, sum))
  expect_equal(na20(fsum(mtcNA, g)), BY(mtcNA, g, sum, na.rm = TRUE)) # error, sum(NA) give 0
})

test_that("fsum performs numerically stable", {
  expect_true(all_obj_equal(replicate(50, fsum(1), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(NA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(NA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(x), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(x, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(xNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(xNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(m), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(m, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mtcars), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mtcars, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mtcNA, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mtcNA), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(x, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(x, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(xNA, f, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(xNA, f), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(m, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(m, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mNA, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mtcars, g), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mtcars, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mtcNA, g, na.rm = FALSE), simplify = FALSE)))
  expect_true(all_obj_equal(replicate(50, fsum(mtcNA, g), simplify = FALSE)))
})

test_that("fsum handles special values in the right way", {
  expect_equal(fsum(NA), NA_real_)
  expect_equal(fsum(NaN), NaN)
  expect_equal(fsum(Inf), Inf)
  expect_equal(fsum(-Inf), -Inf)
  expect_equal(fsum(TRUE), 1)
  expect_equal(fsum(FALSE), 0)
  expect_equal(fsum(NA, na.rm = FALSE), NA_real_)
  expect_equal(fsum(NaN, na.rm = FALSE), NaN)
  expect_equal(fsum(Inf, na.rm = FALSE), Inf)
  expect_equal(fsum(-Inf, na.rm = FALSE), -Inf)
  expect_equal(fsum(TRUE, na.rm = FALSE), 1)
  expect_equal(fsum(FALSE, na.rm = FALSE), 0)
})

test_that("fsum produces errors for wrong input", {
  expect_error(fsum("a"))
  expect_error(fsum(NA_character_))
  expect_error(fsum(mNAc))
  expect_error(fsum(mNAc, f))
  expect_error(fsum(1:2,1:3))
  expect_error(fsum(m,1:31))
  expect_error(fsum(mtcars,1:31))
  expect_error(fsum(wlddev))
  expect_error(fsum(wlddev, wlddev$iso3c))
})

