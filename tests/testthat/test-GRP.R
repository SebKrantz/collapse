context("GRP")

mtcNA <- na_insert(mtcars)

test_that("GRP works as intended", {
 expect_visible(GRP(unname(as.list(mtcars))))
 expect_visible(GRP(unname(as.list(mtcars)), 8:9))
 expect_visible(GRP(mtcars$mpg))
 expect_visible(GRP(1:10))
 expect_visible(GRP(1:10, order = -1))
 expect_visible(GRP(mtcNA$mpg))
 expect_visible(GRP(mtcNA$mpg, return.groups = FALSE))
 expect_visible(GRP(mtcNA$mpg, return.groups = FALSE, return.order = TRUE))
 expect_visible(GRP(mtcNA$mpg, na.last = FALSE))
 expect_visible(GRP(mtcNA$mpg, na.last = FALSE, order = -1))
 expect_visible(GRP(mtcNA$mpg, na.last = FALSE, order = -1, return.order = TRUE))
 expect_visible(GRP(list(a = 1:3, b = 1:3)))
 expect_visible(GRP(mtcars))
 expect_visible(GRP(mtcNA))
 expect_visible(GRP(mtcNA, return.groups = FALSE))
 expect_visible(GRP(mtcNA, return.groups = FALSE, return.order = TRUE))
 expect_visible(GRP(mtcNA, na.last = FALSE))
 expect_visible(GRP(mtcNA, na.last = FALSE, order = -1))
 expect_visible(GRP(mtcNA, na.last = FALSE, order = -1, return.order = TRUE))
 expect_visible(GRP(wlddev))
 expect_visible(GRP(wlddev, return.groups = FALSE))
 expect_true(all_obj_equal(GRP(mtcars, ~ cyl + vs + am)[1:7],
                           GRP(mtcars, c("cyl","vs","am"))[1:7],
                           GRP(mtcars, c(2,8:9))[1:7]))
})

test_that("GRP gives errors for wrong input", {
  expect_error(GRP(~ bla))
  expect_error(GRP(1:10, 1))
  expect_error(GRP(1:10, ~ cyl))
  expect_error(GRP(1:10, "cyl"))
  expect_error(GRP(1:10, order = 0))
  expect_error(GRP(1:10, order = 2))
  expect_error(GRP(mtcars, TRUE))
  expect_error(GRP(mtcars, ~ cyl + bla))
  expect_error(GRP(mtcars, c("bal","cyl")))
  expect_error(GRP(mtcars, 11:12))
  expect_error(GRP(list(a = 1:3, b = 1:4)))
})

test_that("GRP <> factor conversions run seamlessly", {
  expect_identical(unclass(iris$Species), unclass(as.factor.GRP(GRP(iris$Species)))) # as.factor.GRP always adds class "na.included"
  expect_identical(unclass(`vlabels<-`(wlddev$iso3c, "label", NULL)), unclass(as.factor.GRP(GRP(wlddev$iso3c))))
  int <- sample.int(10,100,TRUE)
  expect_identical(unclass(qF(int)), unclass(as.factor.GRP(GRP(int))))
  expect_identical(unclass(qF(int)), unclass(as.factor.GRP(GRP(qF(int)))))
  intNA <- int
  intNA[sample(100,20)] <- NA
  expect_identical(unclass(qF(intNA)), unclass(as.factor.GRP(GRP(intNA)))) # NA_INTEGER is a very small integer -> last in both cases
  expect_identical(unclass(qF(intNA)), unclass(as.factor.GRP(GRP(qF(intNA)))))
  dblNA <- as.double(intNA)
  expect_false(identical(unclass(qF(dblNA)), unclass(as.factor.GRP(GRP(dblNA))))) # qF with na.exclude = TRUE retains double NA's...
  expect_false(identical(unclass(qF(dblNA)), unclass(as.factor.GRP(GRP(qF(dblNA))))))
  expect_identical(qF(dblNA, na.exclude = FALSE), as.factor.GRP(GRP(dblNA)))
  expect_identical(qF(dblNA, na.exclude = FALSE), as.factor.GRP(GRP(qF(dblNA))))
  expect_identical(qF(dblNA, na.exclude = FALSE), as.factor.GRP(GRP(qF(dblNA, na.exclude = FALSE))))
})

# could also do qG to GRP, but qG is same as factor.. and is a programmers function anyway..

