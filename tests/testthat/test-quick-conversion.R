context("quick-conversion")

# rm(list = ls())
set.seed(101)
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

# Test this (plain matrix)
# X = sweep(d, 2L, colMeans(qM(d), na.rm = TRUE), "replace_fill")

setdfdt <- function(x) {
  attr(x, "row.names") <- .set_row_names(length(x[[1L]]))
  class(x) <- c("data.table","data.frame")
  alc(x)
}


test_that("conversions to factor run smoothly", {
  expect_identical(ordered(as.factor(x)), qF(x, ordered = TRUE))
  expect_identical(ordered(as.factor(f)), qF(f, ordered = TRUE))
  expect_identical(as.integer(as.factor(xNA)), as.integer(qF(xNA, ordered = TRUE)))
  expect_identical(as.integer(as.factor(fNA)), as.integer(qF(fNA, ordered = TRUE)))
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

test_that("mrtl and mctl work well", {
  expect_equal(mctl(m), lapply(seq_col(m), function(i) unattrib(m[, i])))
  expect_equal(mctl(m, TRUE), setNames(lapply(seq_col(m), function(i) unattrib(m[, i])), colnames(m)))
  expect_equal(mctl(m, TRUE, "data.frame"), mtcars)
  expect_equal(mctl(m, TRUE, "data.table"), qDT(mtcars))
  expect_equal(mctl(m, FALSE, "data.frame"), setRownames(setNames(mtcars, paste0("V", seq_col(mtcars)))))
  expect_equal(mctl(m, FALSE, "data.table"), qDT(setNames(mtcars, paste0("V", seq_col(mtcars)))))

  expect_equal(mrtl(m), lapply(seq_row(m), function(i) unattrib(m[i, ])))
  expect_equal(mrtl(m, TRUE), setNames(lapply(seq_row(m), function(i) unattrib(m[i, ])), rownames(m)))
  expect_equal(mrtl(m, TRUE, "data.frame"), as.data.frame(t(m)))
  expect_equal(mrtl(m, TRUE, "data.table"), qDT(as.data.frame(t(m))))
  expect_equal(mrtl(m, FALSE, "data.frame"), setRownames(setNames(as.data.frame(t(m)), paste0("V", seq_row(mtcars)))))
  expect_equal(mrtl(m, FALSE, "data.table"), qDT(setNames(as.data.frame(t(m)), paste0("V", seq_row(mtcars)))))
})

test_that("qM keep.attr and class options work as intended", {
  expect_identical(qM(m), m)
  expect_identical(qM(m, keep.attr = TRUE), m)
  expect_identical(qM(m, keep.attr = TRUE, class = "matrix"), `oldClass<-`(m, "matrix"))
  expect_identical(qM(m, class = "matrix"), `oldClass<-`(m, "matrix"))

  expect_identical(qM(mtcars), m)
  expect_identical(qM(mtcars, keep.attr = TRUE), m)
  expect_identical(qM(mtcars, keep.attr = TRUE, class = "matrix"), `oldClass<-`(m, "matrix"))
  expect_identical(qM(mtcars, class = "matrix"), `oldClass<-`(m, "matrix"))

  gmtcars <- `attr<-`(fgroup_by(mtcars, cyl, vs, am), "was.tibble", NULL)
  expect_identical(qM(gmtcars), m)
  expect_identical(qM(gmtcars, keep.attr = TRUE), `attr<-`(m, "groups", attr(gmtcars, "groups")))
  expect_identical(qM(gmtcars, keep.attr = TRUE, class = "matrix"), `oldClass<-`(`attr<-`(m, "groups", attr(gmtcars, "groups")), "matrix"))
  expect_identical(qM(gmtcars, class = "matrix"), `oldClass<-`(m, "matrix"))

  expect_identical(qM(EuStockMarkets, keep.attr = TRUE), EuStockMarkets)
  expect_identical(qM(EuStockMarkets), unclass(`attr<-`(EuStockMarkets, "tsp", NULL)))
  expect_false(identical(qM(EuStockMarkets), EuStockMarkets))
  expect_false(identical(qM(EuStockMarkets, keep.attr = TRUE, class = "matrix"), EuStockMarkets))

  tsl <- list(a = AirPassengers, b = AirPassengers)
  expect_identical(qM(tsl, keep.attr = TRUE), do.call(cbind, tsl))
  expect_identical(qM(tsl), unclass(`attr<-`(do.call(cbind, tsl), "tsp", NULL)))
  expect_false(identical(qM(tsl), do.call(cbind, tsl)))
  expect_false(identical(qM(tsl, keep.attr = TRUE, class = "matrix"), do.call(cbind, tsl)))

})
