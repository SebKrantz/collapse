context("quick-conversion")

x <- rnorm(10)
xNA <- x
xNA[c(3,10)] <- NA
f <- sample.int(3, 10, TRUE)
fNA <- f
fNA[c(3,10)] <- NA
l1 <- replicate(10, rnorm(10), simplify = FALSE)
l2 <- as.list(mtcars)
m <- as.matrix(mtcars)
m2 <- replicate(10, rnorm(10))

setdfdt <- function(x) {
  attr(x, "row.names") <- .set_row_names(length(x[[1L]]))
  class(x) <- c("data.table","data.frame")
  x
}

# TODO: Possibly add keep.attributes argument to qM, qDF and qDT

test_that("conversions to factor run smoothly", {
  expect_identical(ordered(as.factor(x)), qF(x, ordered = TRUE))
  expect_identical(ordered(as.factor(f)), qF(f, ordered = TRUE))
  expect_identical(as.integer(as.factor(xNA)), as.integer(qF(xNA, ordered = TRUE)))
  # expect_identical(as.integer(as.factor(fNA)), as.integer(qF(fNA, ordered = TRUE))) # Problem here: For integers na.exclude = TRUE doesn't work !!
  expect_identical(as.integer(as.factor(x)), as.integer(qG(x, ordered = TRUE)))
  expect_identical(as.integer(as.factor(f)), as.integer(qF(f, ordered = TRUE)))
  expect_identical(as.integer(as.factor(xNA)), as.integer(qG(xNA, ordered = TRUE)))
  expect_identical(as.integer(qF(fNA, ordered = TRUE)), as.integer(qG(fNA, ordered = TRUE)))
})

test_that("conversions to matrix run smoothly", {
  expect_identical(do.call(cbind, l1), qM(l1))
  expect_identical(do.call(cbind, l2), qM(l2))
  expect_identical(as.matrix(mtcars), qM(mtcars))
  expect_identical(`dimnames<-`(as.matrix(x), list(NULL, "x")), qM(x))
  expect_identical(qM(m), m)
  expect_identical(qM(m2), m2)
})

test_that("conversions to data.frame / data.table run smoothly", {
  expect_identical(setNames(as.data.frame(l1), paste0("V",1:10)), qDF(l1))
  expect_identical(as.data.frame(l2), qDF(l2))
  expect_identical(as.data.frame(m), qDF(m))
  expect_identical(as.data.frame(m2), qDF(m2))
  expect_identical(as.data.frame(x), qDF(x))
  expect_identical(qDF(mtcars), mtcars)

  expect_identical(setdfdt(setNames(as.data.frame(l1), paste0("V",1:10))), qDT(l1))
  expect_identical(setdfdt(as.data.frame(l2)), qDT(l2))
  expect_identical(setdfdt(as.data.frame(m)), qDT(m))
  expect_identical(setdfdt(as.data.frame(m2)), qDT(m2))
  expect_identical(setdfdt(as.data.frame(x)), qDT(x))
  expect_identical(qDT(mtcars), setdfdt(mtcars))
})

test_that("double-conversions are ok", {
  expect_identical(qDF(qDT(mtcars)), setRownames(mtcars))
  expect_identical(qM(qDT(m)), setRownames(m, NULL))
  expect_identical(qM(qDF(m)), m)
})
